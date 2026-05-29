/**
 * QMAT Command Line Interface
 *
 * A standalone CLI to run QMAT's core medial axis computation and simplification
 * without Qt GUI dependencies.
 *
 * Supported file formats:
 *   .off  - Object File Format
 *   .obj  - Wavefront OBJ format
 *
 * Usage:
 *   qmat_cli <input.off|input.obj> [options]
 *
 * Options:
 *   --simplify <N>     Simplify to N vertices (default: no simplification)
 *   --k <value>        K factor for slab initialization (default: 0.00001)
 *   --output <prefix>  Output file prefix (default: input filename without extension)
 *   --help             Show this help message
 *
 * Examples:
 *   qmat_cli model.off
 *   qmat_cli model.obj
 *   qmat_cli model.obj --simplify 500 --k 0.0001 --output result
 */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstring>
#include <ctime>
#include <chrono>
#include <thread>
#include <filesystem>
#include <unordered_map>
#include "ThreeDimensionalShape.h"
#include "ObjLoader.h"
#ifdef QMAT_WITH_VCGLIB
#  include "VcgDirectSimplifier.h"
#endif

// Export the current (live, non-compact) slab mesh state as an OFF file.
// Vertices are written in their original indexed order; deleted vertices are
// skipped and a compact remapping is built on the fly.
// Non-static so QmatVisualizer can call it from its ImGui panel.
void ExportMatAsOff(const SlabMesh& sm, const std::string& path)
{
    std::unordered_map<unsigned, unsigned> vid_to_idx;
    vid_to_idx.reserve(sm.numVertices);
    std::vector<unsigned> active_vids;
    active_vids.reserve(sm.numVertices);
    for (unsigned i = 0; i < sm.vertices.size(); ++i)
    {
        if (!sm.vertices[i].first) continue;
        vid_to_idx[i] = (unsigned)active_vids.size();
        active_vids.push_back(i);
    }

    unsigned n_faces = 0;
    for (unsigned i = 0; i < sm.faces.size(); ++i)
        if (sm.faces[i].first) ++n_faces;

    std::ofstream f(path);
    if (!f) { std::cerr << "[ExportMatAsOff] Cannot open " << path << "\n"; return; }

    f << "OFF\n";
    f << active_vids.size() << " " << n_faces << " 0\n";

    const double scale = sm.pmesh ? sm.pmesh->bb_diagonal_length : 1.0;
    f << std::fixed << std::setprecision(10);
    for (unsigned vid : active_vids)
    {
        const auto& c = sm.vertices[vid].second->sphere.center;
        f << c[0]*scale << " " << c[1]*scale << " " << c[2]*scale << "\n";
    }

    for (unsigned i = 0; i < sm.faces.size(); ++i)
    {
        if (!sm.faces[i].first) continue;
        const auto& vset = sm.faces[i].second->vertices_;
        f << vset.size();
        for (unsigned v : vset)
            f << " " << vid_to_idx.at(v);
        f << "\n";
    }
    f.close();
    std::cout << "[Export] OFF written to: " << path << "\n";
}

// Export everything needed to re-create the simplified-MAT visualization in an
// external project: per-primitive ids (struct_ids, cluster type, topo type,
// rejection reason) plus the id->name/color legends and the struct_id color
// formula. One self-describing JSON file, usable from any language.
//
// Must be called BEFORE SlabMesh::Export() — that routine calls AdjustStorage()
// which compacts storage and invalidates edge_last_rejection.
static void ExportSimpVisualizeInfo(const SlabMesh& sm, const std::string& path)
{
    // Local copy of the cluster-type palette so this function is usable in
    // CLI-only builds (the main_cli.cpp palette lives behind
    // #ifdef QMAT_WITH_POLYSCOPE).  Keep in sync with kClusterTypeColors
    // elsewhere in this file.
    static constexpr std::array<const char*, 15> ct_names = {{
        "T0", "T1_spike", "T2", "T3", "T4", "T5", "T1_non_spike",
        "MS_Unknown", "MS_Sheet", "MS_Seam", "MS_Boundary", "MS_Junction",
        "MS_Sheet_Boundary", "MS_Seam_Boundary", "MS_Junction_Boundary",
    }};
    static constexpr std::array<std::array<int,3>, 15> ct_rgb = {{
        {230,   0, 230}, { 26,  26,  26}, {  0, 217, 255}, {255, 128,   0},
        {255,  26,  26}, {255, 255, 255}, {140, 140, 140}, { 89,  89,  89},
        {  0, 255,  77}, {255, 230,   0}, {  0, 128, 255}, {255,   0, 128},
        {  0, 230, 255}, {255,  89,   0}, {153,   0, 255},
    }};
    // Edge topo_type names + colors. Source of truth: SlabEdge::TopoType (SlabMesh.h)
    // and kEdgeTopoTypeColors (above, polyscope-only). Kept as a local copy here so
    // the export works in CLI-only builds. Indices match SlabEdge::TopoType values.
    static constexpr std::array<const char*, 6> et_names = {{
        "Unknown", "Sheet", "Seam", "Boundary", "Seam_Boundary", "Orphan",
    }};
    static constexpr std::array<std::array<int,3>, 6> et_rgb = {{
        {140, 140, 140}, {  0, 255,  77}, {255, 230,   0},
        {  0, 128, 255}, {255,  89,   0}, {230,   0, 230},
    }};

    std::ofstream f(path);
    if (!f) {
        std::cerr << "[ExportSimpVisualizeInfo] cannot open: " << path << "\n";
        return;
    }

    // Compact old->new index maps for active vertices/edges/faces.
    std::vector<unsigned> newv(sm.vertices.size(), UINT_MAX);
    unsigned cv = 0, ce = 0, cf = 0;
    for (unsigned i = 0; i < (unsigned)sm.vertices.size(); ++i)
        if (sm.vertices[i].first) newv[i] = cv++;
    for (unsigned i = 0; i < (unsigned)sm.edges.size(); ++i)
        if (sm.edges[i].first) ++ce;
    for (unsigned i = 0; i < (unsigned)sm.faces.size(); ++i)
        if (sm.faces[i].first) ++cf;

    const double scale = sm.pmesh ? sm.pmesh->bb_diagonal_length : 1.0;

    auto write_rgb_i = [&](const std::array<int,3>& c) {
        f << "[" << c[0] << "," << c[1] << "," << c[2] << "]";
    };
    auto write_rgb_u8 = [&](const std::array<uint8_t,3>& c) {
        f << "[" << (int)c[0] << "," << (int)c[1] << "," << (int)c[2] << "]";
    };
    auto write_set = [&](const std::set<int>& s) {
        f << "[";
        bool first = true;
        for (int id : s) { if (!first) f << ","; f << id; first = false; }
        f << "]";
    };

    f << std::fixed << std::setprecision(10);
    f << "{\n";

    // ── legends ─────────────────────────────────────────────────────────────
    f << "  \"legends\": {\n";

    f << "    \"cluster_types\": [\n";
    for (size_t i = 0; i < ct_names.size(); ++i) {
        f << "      {\"id\": " << i
          << ", \"name\": \"" << ct_names[i] << "\""
          << ", \"rgb\": ";
        write_rgb_i(ct_rgb[i]);
        f << "}";
        if (i + 1 < ct_names.size()) f << ",";
        f << "\n";
    }
    f << "    ],\n";

    f << "    \"edge_topo_types\": [\n";
    for (size_t i = 0; i < et_names.size(); ++i) {
        f << "      {\"id\": " << i
          << ", \"name\": \"" << et_names[i] << "\""
          << ", \"rgb\": ";
        write_rgb_i(et_rgb[i]);
        f << "}";
        if (i + 1 < et_names.size()) f << ",";
        f << "\n";
    }
    f << "    ],\n";

    f << "    \"rejection_reasons\": [\n";
    const size_t num_reasons = static_cast<size_t>(SlabMesh::RejectionReason::Count);
    for (size_t i = 0; i < num_reasons; ++i) {
        auto rr  = static_cast<SlabMesh::RejectionReason>(i);
        auto rgb = SlabMesh::RejectionReasonColorU8(rr);
        f << "      {\"id\": " << i
          << ", \"name\": \"" << SlabMesh::RejectionReasonName(rr) << "\""
          << ", \"rgb\": ";
        write_rgb_u8(rgb);
        f << "}";
        if (i + 1 < num_reasons) f << ",";
        f << "\n";
    }
    f << "    ],\n";

    f << "    \"struct_id_color\": {\n";
    f << "      \"formula\": \"golden_ratio_hsv\",\n";
    f << "      \"note\": \"For struct_id >= 0: hue = fmod(struct_id * 0.618033988749895, 1.0); rgb = HSV(hue, saturation=0.85, value=0.95). For struct_id < 0 (no struct): rgb = [128,128,128] grey.\"\n";
    f << "    }\n";

    f << "  },\n";

    // ── vertices ────────────────────────────────────────────────────────────
    f << "  \"vertices\": [\n";
    {
        bool first = true;
        for (unsigned i = 0; i < (unsigned)sm.vertices.size(); ++i) {
            if (!sm.vertices[i].first) continue;
            if (!first) f << ",\n";
            first = false;
            const SlabVertex* sv = sm.vertices[i].second;
            const auto& c = sv->sphere.center;
            f << "    {\"pos\": ["
              << c.X() * scale << "," << c.Y() * scale << "," << c.Z() * scale
              << "], \"struct_ids\": ";
            write_set(sv->struct_ids);
            f << ", \"cluster_type\": "
              << (int)static_cast<uint8_t>(sv->nmn_cluster_type);
            // Original (pre-simplification) MAT vertex ids that collapsed into
            // this surviving vertex.  Indexes into the top-level
            // "original_positions" array.
            f << ", \"original_ancestors\": [";
            {
                bool fa = true;
                for (unsigned id : sv->original_ancestors) {
                    if (!fa) f << ",";
                    f << id;
                    fa = false;
                }
            }
            f << "]}";
        }
    }
    f << "\n  ],\n";

    // ── edges ───────────────────────────────────────────────────────────────
    f << "  \"edges\": [\n";
    {
        bool first = true;
        for (unsigned i = 0; i < (unsigned)sm.edges.size(); ++i) {
            if (!sm.edges[i].first) continue;
            if (!first) f << ",\n";
            first = false;
            const SlabEdge* se = sm.edges[i].second;
            unsigned a = newv[se->vertices_.first];
            unsigned b = newv[se->vertices_.second];
            f << "    {\"v\": [" << a << "," << b << "]"
              << ", \"struct_ids\": ";
            write_set(se->struct_ids);
            f << ", \"topo_type\": "
              << (int)static_cast<uint8_t>(se->topo_type);
            auto rit = sm.edge_last_rejection.find(i);
            if (rit != sm.edge_last_rejection.end())
                f << ", \"rejection_reason\": "
                  << (int)static_cast<uint8_t>(rit->second);
            else
                f << ", \"rejection_reason\": null";
            f << "}";
        }
    }
    f << "\n  ],\n";

    // ── faces ───────────────────────────────────────────────────────────────
    f << "  \"faces\": [\n";
    {
        bool first = true;
        for (unsigned i = 0; i < (unsigned)sm.faces.size(); ++i) {
            if (!sm.faces[i].first) continue;
            if (!first) f << ",\n";
            first = false;
            const SlabFace* sf = sm.faces[i].second;
            auto it = sf->vertices_.begin();
            unsigned a = newv[*it++];
            unsigned b = newv[*it++];
            unsigned c = newv[*it];
            f << "    {\"v\": [" << a << "," << b << "," << c << "]"
              << ", \"struct_id\": " << sf->struct_id << "}";
        }
    }
    f << "\n  ],\n";

    // ── original positions ──────────────────────────────────────────────────
    // Snapshot of pre-simplification MAT vertex positions, indexed by id.
    // Each surviving vertex's "original_ancestors" list contains ids into this
    // array.  Positions follow the same convention as "vertices[].pos":
    // origin-centered, world-scale (multiplied by bb_diagonal_length).
    f << "  \"original_positions\": [\n";
    {
        bool first = true;
        for (size_t i = 0; i < sm.original_positions.size(); ++i) {
            if (!first) f << ",\n";
            first = false;
            const auto& p = sm.original_positions[i];
            f << "    [" << p[0] * scale << ","
                         << p[1] * scale << ","
                         << p[2] * scale << "]";
        }
    }
    f << "\n  ]\n";

    f << "}\n";

    std::cout << "[ExportSimpVisualizeInfo] wrote " << cv << " verts, "
              << ce << " edges, " << cf << " faces to " << path << "\n";
}

// Export the initial (pre-simplification) MAT geometry plus per-primitive
// struct_ids — vertex positions, edge and face connectivity, and the
// struct_ids set on each.  Stripped-down sibling of ExportSimpVisualizeInfo:
// no cluster types, rejection reasons, or ancestry, since none of those
// have meaning before simplification runs.  The struct_id colour legend is
// kept so visualizers can render the same colormap as the simplified export.
static void ExportInitialVisualizeInfo(const SlabMesh& sm, const std::string& path)
{
    std::ofstream f(path);
    if (!f) {
        std::cerr << "[ExportInitialVisualizeInfo] cannot open: " << path << "\n";
        return;
    }

    // Compact old->new index maps for active vertices/edges/faces.
    std::vector<unsigned> newv(sm.vertices.size(), UINT_MAX);
    unsigned cv = 0, ce = 0, cf = 0;
    for (unsigned i = 0; i < (unsigned)sm.vertices.size(); ++i)
        if (sm.vertices[i].first) newv[i] = cv++;
    for (unsigned i = 0; i < (unsigned)sm.edges.size(); ++i)
        if (sm.edges[i].first) ++ce;
    for (unsigned i = 0; i < (unsigned)sm.faces.size(); ++i)
        if (sm.faces[i].first) ++cf;

    const double scale = sm.pmesh ? sm.pmesh->bb_diagonal_length : 1.0;

    auto write_set = [&](const std::set<int>& s) {
        f << "[";
        bool first = true;
        for (int id : s) { if (!first) f << ","; f << id; first = false; }
        f << "]";
    };

    f << std::fixed << std::setprecision(10);
    f << "{\n";

    // ── legend (struct_id colour formula only) ─────────────────────────────
    f << "  \"legends\": {\n";
    f << "    \"struct_id_color\": {\n";
    f << "      \"formula\": \"golden_ratio_hsv\",\n";
    f << "      \"note\": \"For struct_id >= 0: hue = fmod(struct_id * 0.618033988749895, 1.0); rgb = HSV(hue, saturation=0.85, value=0.95). For struct_id < 0 (no struct): rgb = [128,128,128] grey.\"\n";
    f << "    }\n";
    f << "  },\n";

    // ── vertices ────────────────────────────────────────────────────────────
    f << "  \"vertices\": [\n";
    {
        bool first = true;
        for (unsigned i = 0; i < (unsigned)sm.vertices.size(); ++i) {
            if (!sm.vertices[i].first) continue;
            if (!first) f << ",\n";
            first = false;
            const SlabVertex* sv = sm.vertices[i].second;
            const auto& c = sv->sphere.center;
            f << "    {\"pos\": ["
              << c.X() * scale << "," << c.Y() * scale << "," << c.Z() * scale
              << "], \"struct_ids\": ";
            write_set(sv->struct_ids);
            f << "}";
        }
    }
    f << "\n  ],\n";

    // ── edges ───────────────────────────────────────────────────────────────
    f << "  \"edges\": [\n";
    {
        bool first = true;
        for (unsigned i = 0; i < (unsigned)sm.edges.size(); ++i) {
            if (!sm.edges[i].first) continue;
            if (!first) f << ",\n";
            first = false;
            const SlabEdge* se = sm.edges[i].second;
            unsigned a = newv[se->vertices_.first];
            unsigned b = newv[se->vertices_.second];
            f << "    {\"v\": [" << a << "," << b << "]"
              << ", \"struct_ids\": ";
            write_set(se->struct_ids);
            f << "}";
        }
    }
    f << "\n  ],\n";

    // ── faces ───────────────────────────────────────────────────────────────
    f << "  \"faces\": [\n";
    {
        bool first = true;
        for (unsigned i = 0; i < (unsigned)sm.faces.size(); ++i) {
            if (!sm.faces[i].first) continue;
            if (!first) f << ",\n";
            first = false;
            const SlabFace* sf = sm.faces[i].second;
            auto it = sf->vertices_.begin();
            unsigned a = newv[*it++];
            unsigned b = newv[*it++];
            unsigned c = newv[*it];
            f << "    {\"v\": [" << a << "," << b << "," << c << "]"
              << ", \"struct_id\": " << sf->struct_id << "}";
        }
    }
    f << "\n  ]\n";

    f << "}\n";

    std::cout << "[ExportInitialVisualizeInfo] wrote " << cv << " verts, "
              << ce << " edges, " << cf << " faces to " << path << "\n";
}

#ifdef QMAT_WITH_POLYSCOPE
#  include "QmatVisualizer.h"
#  ifdef QMAT_WITH_VCGLIB
#    include "VdeVisualizer.h"
#  endif
#endif

// Simple command line argument parsing
struct CLIOptions {
    std::string inputFile;
    std::string outputPrefix;
    int simplifyTarget = -1;  // --simplify <N>  : vertex target (-1 = no simplification)
    int faceTarget = -1;      // --nf <N>        : face target (-1 = unset)

    // --simplifier picks which decimator runs.  All three read the same target
    // inputs above.  qmat/vcg-port both go through SlabMesh::Simplify and which
    // one runs is fixed at compile time by ONLY_USE_QEM_CONDITION_CHECKS — a
    // mismatched flag at runtime warns.  vcg-direct uses vcglib directly and
    // requires -DQMAT_WITH_VCGLIB=ON.
    enum class SimplifierMode { Qmat, VcgPort, VcgDirect };
#if defined(ONLY_USE_QEM_CONDITION_CHECKS)
    SimplifierMode simplifier = SimplifierMode::VcgPort;
#else
    SimplifierMode simplifier = SimplifierMode::Qmat;
#endif
    double k = 0.00001;
    double featureAngleDeg = 10.0; // turning-angle threshold for sharp feature protection
    // VCG-faithful QEM thresholds (ONLY_USE_QEM_CONDITION_CHECKS path).
    // Defaults mirror VcgQuadricParameter's struct defaults.
    double qemQualityThr = 0.3;     // --qem-quality-thr  → VcgQuadricParameter::QualityThr
    double qemBoundaryWeight = 1.0; // --qem-boundary-weight = MeshLab UI "Boundary
                                    // Preserving Weight" slider; effective weight is
                                    // 0.5 * this (meshfilter.cpp:14).  Default 1.0 → 0.5.
    // VCG-faithful QEM toggles (MeshLab QEM checkboxes), default to MeshLab's GUI
    // defaults; override with --qem-* <0|1> so no recompile is needed.
    bool qemOptimalPlacement = true; // --qem-optimal-placement (MeshLab "Optimal position")
    bool qemPreserveTopology = true; // --qem-preserve-topology (MeshLab "Preserve Topology")
    bool qemNormalCheck      = true; // --qem-normal-check       (MeshLab "Preserve Normal")
    bool visualize = false;
    bool showHelp = false;
    bool valid = true;
    std::string errorMessage;
#ifdef USE_MATSTRUCT_INITIALIZATION
    std::string matstructFile;  // path to typed .ma file from external tool
#endif
};

void printUsage(const char* programName) {
    std::cout << "QMAT Command Line Interface\n"
              << "Compute medial axis and optionally simplify.\n\n"
              << "Usage:\n"
              << "  " << programName << " <input.off|input.obj> [options]\n\n"
              << "Supported formats:\n"
              << "  .off               Object File Format\n"
              << "  .obj               Wavefront OBJ format\n\n"
              << "Options:\n"
              << "  --simplify <N>     Simplify to N vertices (default: no simplification)\n"
              << "  --nf <N>           Stop QEM when MAT face count reaches N (requires --simplify; default: off)\n"
              << "  --qem-quality-thr <v>     QEM QualityThr threshold (default: 0.3)\n"
              << "  --qem-boundary-weight <v> MeshLab UI 'Boundary Preserving Weight' slider;\n"
              << "                            effective weight = 0.5*v (default: 1.0 -> 0.5)\n"
              << "  --qem-optimal-placement <0|1>  Optimal vertex position (default: 1)\n"
              << "  --qem-preserve-topology <0|1>  Preserve topology / link condition (default: 1)\n"
              << "  --qem-normal-check <0|1>       Preserve normal (default: 1)\n"
              << "  --simplifier <mode>  qmat | vcg-port | vcg-direct\n"
              << "                       qmat/vcg-port: SlabMesh::Simplify (compile-time pick).\n"
              << "                       vcg-direct: vcglib TriEdgeCollapseQuadric (needs QMAT_WITH_VCGLIB).\n"
              << "  --k <value>        K factor for slab initialization (default: 0.00001)\n"
              << "  --feature-angle <deg>  Turning-angle threshold for sharp feature protection (default: 30)\n"
              << "  --output <prefix>  Output file prefix (default: input filename)\n"
#ifdef QMAT_WITH_POLYSCOPE
              << "  --visualize        Open live Polyscope viewer during simplification\n"
#endif
#ifdef USE_MATSTRUCT_INITIALIZATION
              << "  --matstruct <file> Path to typed .ma file from external MAT tool\n"
#endif
              << "  --help             Show this help message\n\n"
              << "Examples:\n"
              << "  " << programName << " model.off\n"
              << "  " << programName << " model.obj\n"
              << "  " << programName << " model.obj --simplify 500 --k 0.0001 --output result\n";
}

CLIOptions parseArguments(int argc, char* argv[]) {
    CLIOptions options;

    if (argc < 2) {
        options.valid = false;
        options.errorMessage = "No input file specified.";
        return options;
    }

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];

        if (arg == "--help" || arg == "-h") {
            options.showHelp = true;
            return options;
        }
        else if (arg == "--simplify") {
            if (i + 1 >= argc) {
                options.valid = false;
                options.errorMessage = "--simplify requires a value.";
                return options;
            }
            try {
                options.simplifyTarget = std::stoi(argv[++i]);
                if (options.simplifyTarget <= 0) {
                    options.valid = false;
                    options.errorMessage = "--simplify value must be positive.";
                    return options;
                }
            } catch (...) {
                options.valid = false;
                options.errorMessage = "Invalid value for --simplify.";
                return options;
            }
        }
        else if (arg == "--nf") {
            if (i + 1 >= argc) {
                options.valid = false;
                options.errorMessage = "--nf requires a value.";
                return options;
            }
            try {
                options.faceTarget = std::stoi(argv[++i]);
                if (options.faceTarget <= 0) {
                    options.valid = false;
                    options.errorMessage = "--nf value must be positive.";
                    return options;
                }
            } catch (...) {
                options.valid = false;
                options.errorMessage = "Invalid value for --nf.";
                return options;
            }
        }
        else if (arg == "--qem-quality-thr") {
            if (i + 1 >= argc) {
                options.valid = false;
                options.errorMessage = "--qem-quality-thr requires a value.";
                return options;
            }
            try {
                options.qemQualityThr = std::stod(argv[++i]);
                if (options.qemQualityThr <= 0) {
                    options.valid = false;
                    options.errorMessage = "--qem-quality-thr value must be positive.";
                    return options;
                }
            } catch (...) {
                options.valid = false;
                options.errorMessage = "Invalid value for --qem-quality-thr.";
                return options;
            }
        }
        else if (arg == "--qem-boundary-weight") {
            if (i + 1 >= argc) {
                options.valid = false;
                options.errorMessage = "--qem-boundary-weight requires a value.";
                return options;
            }
            try {
                options.qemBoundaryWeight = std::stod(argv[++i]);
                if (options.qemBoundaryWeight < 0) {
                    options.valid = false;
                    options.errorMessage = "--qem-boundary-weight value must be non-negative.";
                    return options;
                }
            } catch (...) {
                options.valid = false;
                options.errorMessage = "Invalid value for --qem-boundary-weight.";
                return options;
            }
        }
        else if (arg == "--simplifier") {
            if (i + 1 >= argc) {
                options.valid = false;
                options.errorMessage = "--simplifier requires a value (qmat|vcg-port|vcg-direct).";
                return options;
            }
            std::string v = argv[++i];
            if      (v == "qmat")       options.simplifier = CLIOptions::SimplifierMode::Qmat;
            else if (v == "vcg-port")   options.simplifier = CLIOptions::SimplifierMode::VcgPort;
            else if (v == "vcg-direct") options.simplifier = CLIOptions::SimplifierMode::VcgDirect;
            else {
                options.valid = false;
                options.errorMessage = "--simplifier must be qmat, vcg-port, or vcg-direct.";
                return options;
            }
        }
        else if (arg == "--qem-optimal-placement" || arg == "--qem-preserve-topology"
                 || arg == "--qem-normal-check") {
            if (i + 1 >= argc) {
                options.valid = false;
                options.errorMessage = arg + " requires a value (0|1).";
                return options;
            }
            std::string v = argv[++i];
            bool on;
            if (v == "1" || v == "true" || v == "on")        on = true;
            else if (v == "0" || v == "false" || v == "off") on = false;
            else {
                options.valid = false;
                options.errorMessage = arg + " must be 0 or 1.";
                return options;
            }
            if      (arg == "--qem-optimal-placement") options.qemOptimalPlacement = on;
            else if (arg == "--qem-preserve-topology") options.qemPreserveTopology = on;
            else                                       options.qemNormalCheck      = on;
        }
        else if (arg == "--k") {
            if (i + 1 >= argc) {
                options.valid = false;
                options.errorMessage = "--k requires a value.";
                return options;
            }
            try {
                options.k = std::stod(argv[++i]);
                if (options.k <= 0) {
                    options.valid = false;
                    options.errorMessage = "--k value must be positive.";
                    return options;
                }
            } catch (...) {
                options.valid = false;
                options.errorMessage = "Invalid value for --k.";
                return options;
            }
        }
        else if (arg == "--output") {
            if (i + 1 >= argc) {
                options.valid = false;
                options.errorMessage = "--output requires a value.";
                return options;
            }
            options.outputPrefix = argv[++i];
        }
#ifdef USE_MATSTRUCT_INITIALIZATION
        else if (arg == "--matstruct") {
            if (i + 1 >= argc) {
                options.valid = false;
                options.errorMessage = "--matstruct requires a file path.";
                return options;
            }
            options.matstructFile = argv[++i];
        }
#endif
        else if (arg == "--feature-angle") {
            if (i + 1 >= argc) {
                options.valid = false;
                options.errorMessage = "--feature-angle requires a value in degrees.";
                return options;
            }
            try { options.featureAngleDeg = std::stod(argv[++i]); std::cout << "Feature angle threshold set to " << options.featureAngleDeg << " degrees.\n"; }
            catch (...) {
                options.valid = false;
                options.errorMessage = "Invalid value for --feature-angle.";
                return options;
            }
        }
        else if (arg == "--visualize") {
#ifdef QMAT_WITH_POLYSCOPE
            options.visualize = true;
#else
            std::cerr << "Warning: --visualize ignored (built without QMAT_WITH_POLYSCOPE)\n";
#endif
        }
        else if (arg[0] == '-') {
            options.valid = false;
            options.errorMessage = "Unknown option: " + arg;
            return options;
        }
        else {
            // Positional argument - input file
            if (options.inputFile.empty()) {
                options.inputFile = arg;
            } else {
                options.valid = false;
                options.errorMessage = "Multiple input files specified.";
                return options;
            }
        }
    }

    if (options.inputFile.empty()) {
        options.valid = false;
        options.errorMessage = "No input file specified.";
        return options;
    }

    // --nf only refines an existing --simplify request; ignore it otherwise.
    if (options.faceTarget > 0 && options.simplifyTarget <= 0) {
        std::cerr << "Warning: --nf ignored (requires --simplify to also be set).\n";
        options.faceTarget = -1;
    }

    // Set default output prefix from input filename
    if (options.outputPrefix.empty()) {
        options.outputPrefix = options.inputFile;
        // Remove .off or .obj extension if present
        size_t dotPos = options.outputPrefix.rfind('.');
        if (dotPos != std::string::npos) {
            std::string ext = options.outputPrefix.substr(dotPos);
            if (ext == ".off" || ext == ".OFF" || ext == ".obj" || ext == ".OBJ") {
                options.outputPrefix = options.outputPrefix.substr(0, dotPos);
            }
        }
    }

    // ── redirect all output into a subfolder named after the mesh stem ────────
    // e.g.  path/to/cube_subdiv_500  →  create  path/to/cube_subdiv_500/
    //                                   prefix  path/to/cube_subdiv_500/cube_subdiv_500
    {
        namespace fs = std::filesystem;
        fs::path p(options.outputPrefix);
        std::string stem = p.filename().string();   // final component only
        fs::path outDir  = p.parent_path() / stem;  // sibling folder named after stem
        std::error_code ec;
        fs::create_directories(outDir, ec);
        if (ec)
            std::cerr << "Warning: could not create output directory "
                      << outDir << ": " << ec.message() << "\n";
        options.outputPrefix = (outDir / stem).string();
    }

    return options;
}

// RenderVcgDirectSnapshot moved to VdeVisualizer.cpp (Phase 2).

int main(int argc, char* argv[]) {

    std::cout << "argc = " << argc << "\n";
    // for (int i = 0; i < argc; ++i)
    //     std::cout << "argv[" << i << "] = [" << argv[i] << "]\n";
    // Parse command line arguments
    CLIOptions options = parseArguments(argc, argv);

    if (options.showHelp) {
        printUsage(argv[0]);
        return 0;
    }

    if (!options.valid) {
        std::cerr << "Error: " << options.errorMessage << std::endl;
        std::cerr << "Use --help for usage information." << std::endl;
        return 1;
    }

    std::cout << "QMAT CLI - Medial Axis Computation" << std::endl;
    std::cout << "===================================" << std::endl;
    std::cout << "Input file: " << options.inputFile << std::endl;
    std::cout << "Output prefix: " << options.outputPrefix << std::endl;
    std::cout << "K value: " << options.k << std::endl;
    if (options.simplifyTarget > 0) {
        std::cout << "Simplify target: " << options.simplifyTarget << " vertices" << std::endl;
        if (options.faceTarget > 0)
            std::cout << "Face target (--nf): stop QEM at " << options.faceTarget << " faces" << std::endl;
        std::cout << "QEM QualityThr: " << options.qemQualityThr
                  << ", QEM BoundaryWeight: " << options.qemBoundaryWeight << std::endl;
    }
    std::cout << std::endl;

    // Create ThreeDimensionalShape object
    ThreeDimensionalShape shape;

    // Step 1: Load the mesh file (OFF or OBJ)
    std::cout << "Loading mesh from " << options.inputFile << "..." << std::endl;
    long startTime = clock();

    bool loadSuccess = false;
    if (IsObjFile(options.inputFile)) {
        // Load OBJ file using tinyobjloader
        std::string objError;
        loadSuccess = LoadObjFile(options.inputFile, shape.input, objError);
        if (!loadSuccess) {
            std::cerr << "Error loading OBJ file: " << objError << std::endl;
            return 1;
        }
    } else if (IsOffFile(options.inputFile)) {
        // Load OFF file using CGAL
        std::ifstream stream(options.inputFile.c_str());
        if (!stream) {
            std::cerr << "Error: Could not open file " << options.inputFile << std::endl;
            return 1;
        }
        stream >> shape.input;
        stream.close();
        loadSuccess = true;
    } else {
        std::cerr << "Error: Unsupported file format. Use .off or .obj files." << std::endl;
        return 1;
    }

    // Compute mesh properties
    shape.input.computebb();
    shape.input.GenerateList();
    shape.input.GenerateRandomColor();
    shape.input.compute_normals();

    long loadTime = clock() - startTime;
    std::cout << "  Loaded mesh with " << shape.input.size_of_vertices() << " vertices, "
              << shape.input.size_of_facets() << " faces" << std::endl;
    std::cout << "  Load time: " << loadTime << " ms" << std::endl;

    // Step 2: Create CGAL mesh domain for inside/outside queries
    std::cout << "Creating mesh domain..." << std::endl;
    Polyhedron pol;
    if (IsObjFile(options.inputFile)) {
        // Load OBJ file for mesh domain
        std::string objError;
        if (!LoadObjFile(options.inputFile, pol, objError)) {
            std::cerr << "Error loading OBJ file for mesh domain: " << objError << std::endl;
            return 1;
        }
    } else {
        // Load OFF file for mesh domain
        std::ifstream streamPol(options.inputFile.c_str());
        streamPol >> pol;
        streamPol.close();
    }

    Mesh_domain* domain = new Mesh_domain(pol);
    shape.input.domain = domain;
    shape.input_nmm.domain = domain;
    shape.input_nmm.pmesh = &shape.input;
    shape.input_nmm.meshname = options.outputPrefix;

#ifndef USE_MATSTRUCT_INITIALIZATION
    // Step 3: Compute Delaunay Triangulation and Medial Axis
    // If the .ma file already exists on disk, skip the expensive DT + MAT
    // computation and reuse it directly.
    const std::string maFile = options.outputPrefix + ".ma";
    if (std::filesystem::exists(maFile)) {
        std::cout << "Found existing MA file: " << maFile << std::endl;
        std::cout << "  Skipping DT and MAT computation." << std::endl;
    } else {
        std::cout << "Computing Delaunay Triangulation..." << std::endl;
        startTime = clock();
        shape.input.computedt();
        long dtTime = clock() - startTime;
        std::cout << "  DT computation time: " << dtTime << " ms" << std::endl;

        // Cluster boundary sample points (input mesh vertices) by position + normal.
        // Output: <outputPrefix>_boundary_clusters.txt




        std::cout << "Computing Medial Axis..." << std::endl;
        startTime = clock();
        shape.ComputeInputNMM();
        long maTime = clock() - startTime;
        std::cout << "  MA computation time: " << maTime << " ms" << std::endl;
        std::cout << "  Raw MA exported to: " << maFile << std::endl;
    }
#endif // USE_MATSTRUCT_INITIALIZATION

#ifdef USE_MATSTRUCT_INITIALIZATION
    // Step 3b: Also compute the Voronoi-based initial MAT for comparison.
    // This is the raw MAT derived from the Delaunay triangulation of the input
    // surface samples — the same one produced without USE_MATSTRUCT.
    // All output goes into a dedicated init_vor_mat/ subfolder:
    //   init_vor_mat/<stem>.ma                        — raw .ma format
    //   init_vor_mat/<stem>.ply                       — ASCII PLY (vertex/face/edge)
    //   init_vor_mat/<stem>_sampledpoints.txt
    //   init_vor_mat/<stem>_surf_2_vor_mat_map.txt
    //   init_vor_mat/<stem>_mat_vertices_with_fields.txt
    {
        namespace fs = std::filesystem;
        fs::path pp(options.outputPrefix);

        // Create the dedicated subfolder.
        fs::path vorDir = pp.parent_path() / "init_vor_mat";
        fs::create_directories(vorDir);
        const std::string vorPrefix = (vorDir / pp.filename()).string();

        std::cout << std::endl << "Computing Voronoi MAT for comparison..." << std::endl;
        shape.input_nmm.meshname = vorPrefix;   // Export() and sidecar writers use this

        startTime = clock();
        shape.input.computedt();
        shape.ComputeInputNMM();                // writes .ma and all sidecar files
        long vorTime = clock() - startTime;
        std::cout << "  Voronoi MAT computation time: " << vorTime << " ms" << std::endl;
        std::cout << "  Voronoi MAT written to: " << vorPrefix << ".ma" << std::endl;

        // ── Export as ASCII PLY for MeshLab visualization ─────────────────────
        // PLY natively supports vertex, face, and edge as first-class elements.
        //
        // Bug fix: faces in the MAT are stored as std::set<unsigned> of vertex
        // indices. In degenerate cases two of the three circumcenter tags can be
        // equal, leaving the set with fewer than 3 entries. Writing the actual
        // set size (not a hardcoded 3) keeps the face list correctly aligned so
        // MeshLab does not lose sync and drop subsequent faces/edges.
        const std::string vorPly = vorPrefix + ".ply";
        {
            std::ofstream ply(vorPly);
            if (ply.is_open()) {
                const auto& nmverts = shape.input_nmm.vertices;
                const auto& nmfaces = shape.input_nmm.faces;
                const auto& nmedges = shape.input_nmm.edges;

                // ── Header ────────────────────────────────────────────────────
                ply << "ply\n"
                    << "format ascii 1.0\n"
                    << "element vertex " << nmverts.size() << "\n"
                    << "property double x\n"
                    << "property double y\n"
                    << "property double z\n";

                if (!nmfaces.empty())
                    ply << "element face " << nmfaces.size() << "\n"
                        << "property list uchar int vertex_indices\n";

                if (!nmedges.empty())
                    ply << "element edge " << nmedges.size() << "\n"
                        << "property int vertex1\n"
                        << "property int vertex2\n";

                ply << "end_header\n";

                // ── Vertices ──────────────────────────────────────────────────
                ply << std::fixed << std::setprecision(15);
                for (const auto& bvp : nmverts) {
                    const auto& c = bvp.second->sphere.center;
                    ply << c[0] << " " << c[1] << " " << c[2] << "\n";
                }

                // ── Faces ─────────────────────────────────────────────────────
                if (!nmfaces.empty()) {
                    for (const auto& bfp : nmfaces) {
                        ply << bfp.second->vertices_.size();
                        for (unsigned vi : bfp.second->vertices_)
                            ply << " " << vi;
                        ply << "\n";
                    }
                }

                // ── Edges ─────────────────────────────────────────────────────
                if (!nmedges.empty()) {
                    for (const auto& bep : nmedges)
                        ply << bep.second->vertices_.first
                            << " " << bep.second->vertices_.second << "\n";
                }

                std::cout << "  Voronoi MAT PLY written to: " << vorPly
                          << "  (v=" << nmverts.size()
                          << " f=" << nmfaces.size()
                          << " e=" << nmedges.size() << ")" << std::endl;
            } else {
                std::cerr << "  Warning: could not open " << vorPly << " for writing\n";
            }
        }

        // Restore meshname so the rest of the pipeline uses the normal prefix.
        shape.input_nmm.meshname = options.outputPrefix;
    }
#endif // USE_MATSTRUCT_INITIALIZATION


    // Step 4: If simplification requested, load into slab mesh and simplify
    if (options.simplifyTarget > 0) {
        std::cout << std::endl << "Loading MA for simplification..." << std::endl;

        // Setup slab mesh
        shape.slab_mesh.pmesh = &shape.input;
        shape.slab_mesh.type = 1;
        shape.slab_mesh.k = options.k;
        shape.slab_mesh.bound_weight = 1.0;
        shape.slab_mesh.export_prefix = options.outputPrefix;

        // VCG-faithful QEM thresholds + optional face-count stop (--nf).
        // Consumed only by the ONLY_USE_QEM_CONDITION_CHECKS path in Simplify().
        shape.slab_mesh.qem_quality_thr     = options.qemQualityThr;
        shape.slab_mesh.qem_boundary_weight = options.qemBoundaryWeight;
        shape.slab_mesh.qem_face_target     = options.faceTarget;
        shape.slab_mesh.qem_optimal_placement = options.qemOptimalPlacement;
        shape.slab_mesh.qem_preserve_topology = options.qemPreserveTopology;
        shape.slab_mesh.qem_normal_check      = options.qemNormalCheck;

        // Initialize slab mesh settings (same as GUI initialize())
        shape.slab_mesh.preserve_boundary_method = 0;
        shape.slab_mesh.hyperbolic_weight_type = 3;
        shape.slab_mesh.compute_hausdorff = false;
        shape.slab_mesh.boundary_compute_scale = 0;
        shape.slab_mesh.prevent_inversion = true;

#ifdef USE_MATSTRUCT_INITIALIZATION
        // Load the typed .ma file from the external MAT tool.
        shape.LoadMatstructMA(options.matstructFile);
        // struct_ids sets are populated during LoadMatstructMA; collapsibility is
        // checked inline in CanMerge — no separate pre-computation needed.
        shape.ComputeFeatureEdges();          // detect sharp/concave edges on input mesh
        shape.ExportSurfacemeshFeatureEdges(options.outputPrefix);
        // ClusterNMNBplist() is skipped — T-types are already loaded from file.
#else
        // Load the MA file into the slab mesh (computed above or loaded from cache)
        shape.LoadInputNMM(maFile);
        shape.ComputeFeatureEdges();          // detect sharp/concave edges on input mesh
        shape.ExportSurfacemeshFeatureEdges(options.outputPrefix);
        shape.slab_mesh.ClusterNMNBplist();   // must run before DetermineTopology (T1 filtering)
#endif
        // shape.slab_mesh.DetermineTopology();  // uses T-types to correctly classify topology
        shape.slab_mesh.feature_angle_threshold = options.featureAngleDeg;
        shape.slab_mesh.MarkSharpFeatureVertices(options.featureAngleDeg);

        std::cout << "  Loaded slab mesh with " << shape.slab_mesh.numVertices << " vertices" << std::endl;

        // Initialize slab mesh for simplification
        std::cout << "Initializing slab mesh..." << std::endl;
        startTime = clock();
        long initTime = shape.LoadSlabMesh();
        std::cout << "  Initialization time: " << initTime << " ms" << std::endl;

        // Simplify
        int currentVertices = shape.slab_mesh.numVertices;

        // ── vcg-direct branch: vcg's own TriEdgeCollapseQuadric via vcglib ──
        // Bypasses SlabMesh::Simplify entirely.  Takes --simplify (vertex) and/or
        // --nf (face) as targets; --nf wins if both are set.  Output is shown in
        // polyscope as "VCG Direct Simplified MAT"; QMAT post-export steps are
        // skipped because the slab mesh itself is left untouched.
        if (options.simplifier == CLIOptions::SimplifierMode::VcgDirect) {
#ifdef QMAT_WITH_VCGLIB
            if (options.simplifyTarget <= 0 && options.faceTarget <= 0) {
                std::cout << "Warning: --simplifier vcg-direct requires --simplify or --nf.\n";
            } else {
                startTime = clock();
                shape.slab_mesh.CleanIsolatedVertices();
                ExportMatAsOff(shape.slab_mesh, options.outputPrefix + "_mat_initial.off");
                ExportInitialVisualizeInfo(shape.slab_mesh,
                                           options.outputPrefix + "_initial_visualize_info.json");

                VcgDirectParams p;
                p.TargetFaceNum         = options.faceTarget;
                p.TargetVertexNum       = options.simplifyTarget;
                p.QualityThr            = options.qemQualityThr;
                p.BoundaryQuadricWeight = options.qemBoundaryWeight;
                p.OptimalPlacement      = options.qemOptimalPlacement;
                p.PreserveTopology      = options.qemPreserveTopology;
                p.NormalCheck           = options.qemNormalCheck;

                LiveUpdateCallback live_cb;
#ifdef QMAT_WITH_POLYSCOPE
                VdeVisualizer vde;
                ViewerState& vs = vde.State();
                if (options.visualize) {
                    vs.outputPrefix = options.outputPrefix;
                    vde.Setup(shape.slab_mesh);
                    live_cb = vde.MakeLiveCallback();
                }
#endif

                VcgDirectResult res;
                bool ok = RunVcgDirectSimplify(shape.slab_mesh, p, res, live_cb);
                std::cout << "  vcg-direct time: " << (clock() - startTime) << " ms\n";

#ifdef QMAT_WITH_POLYSCOPE
                if (ok && options.visualize) {
                    vde.Render(res.snapshot);
                    std::cout << "vcg-direct done. Close the viewer to exit.\n";
                    vde.Show();
                }
#endif
            }
#else
            std::cerr << "--simplifier vcg-direct requires -DQMAT_WITH_VCGLIB=ON.\n";
#endif
        }
        else if (options.simplifyTarget >= currentVertices) {
            std::cout << "Warning: Target vertex count (" << options.simplifyTarget
                      << ") >= current count (" << currentVertices << "). Skipping simplification." << std::endl;
        } else {
            int reductionCount = currentVertices - options.simplifyTarget;
            std::cout << "Simplifying from " << currentVertices << " to " << options.simplifyTarget
                      << " vertices (removing " << reductionCount << ")..." << std::endl;

            // Label every vertex with its topological class before simplification.
            // Labels (mirrors check_non_manifold in read_qmat_output.py):
            //   REGULAR    – interior vertex (not on any boundary or non-manifold edge)
            //   BOUNDARY   – on a boundary edge (edge shared by exactly 1 face)
            //   NM_EDGE    – on a non-manifold edge (>2 faces), but not a corner
            //   NM_CORNER  – on both a boundary and a non-manifold edge
            // Simplify will only collapse an edge whose two endpoints share the same label.

            startTime = clock();
            shape.slab_mesh.CleanIsolatedVertices();

            // Export MAT before simplification starts
            ExportMatAsOff(shape.slab_mesh, options.outputPrefix + "_mat_initial.off");

            // JSON snapshot of the initial MAT: vertex positions, edge/face
            // connectivity, and per-primitive struct_ids.  Companion to the
            // post-simplification _simp_visualize_info.json so external tools
            // can compare structure assignments before vs. after simplifying.
            ExportInitialVisualizeInfo(shape.slab_mesh,
                                       options.outputPrefix + "_initial_visualize_info.json");

#ifdef QMAT_WITH_POLYSCOPE
            QmatVisualizer qmat;
            ViewerState& vs = qmat.State();
            vs.outputPrefix = options.outputPrefix;
            if (options.visualize)
                qmat.Setup(shape.slab_mesh);
#endif
            // qmat/vcg-port: which path runs is fixed at compile time by
            // ONLY_USE_QEM_CONDITION_CHECKS — warn if --simplifier disagrees.
#if defined(ONLY_USE_QEM_CONDITION_CHECKS)
            if (options.simplifier == CLIOptions::SimplifierMode::Qmat)
                std::cerr << "[warn] --simplifier qmat ignored: binary built with ONLY_USE_QEM_CONDITION_CHECKS → vcg-port runs.\n";
#else
            if (options.simplifier == CLIOptions::SimplifierMode::VcgPort)
                std::cerr << "[warn] --simplifier vcg-port ignored: binary built without ONLY_USE_QEM_CONDITION_CHECKS → qmat runs.\n";
#endif
            shape.slab_mesh.allow_steep_collapse = true;
            shape.slab_mesh.Simplify(reductionCount);
            long simplifyTime = clock() - startTime;

            std::cout << "  Simplification time: " << simplifyTime << " ms" << std::endl;
            std::cout << "  Final vertex count: " << shape.slab_mesh.numVertices << std::endl;

#ifdef QMAT_WITH_POLYSCOPE
            if (options.visualize) {
                qmat.RenderFinal(shape.slab_mesh);
                std::cout << "Simplification done. Close the viewer window to exit.\n";
                qmat.Show();
            }
#endif


            // Must run before Export(): Export() calls AdjustStorage() which
            // remaps all edge IDs, invalidating edge_last_rejection.
            // radius is a fraction of the geometry's bbox diagonal (computed inside
            // ExportSkeletonPLY); the old 1.0f was an absolute value that, against the
            // normalised ≈[-1,1] coordinates, made every cylinder model-sized.
            shape.slab_mesh.ExportSkeletonPLY(options.outputPrefix + "_rejection_skeleton.ply");
            std::cout << "  Skeleton PLY exported: " << options.outputPrefix << "_rejection_skeleton.ply\n";

            // JSON with per-primitive struct_ids / cluster-type / topo-type /
            // rejection-reason plus id->color legends, for external
            // re-visualization tools.  Must run before Export() (same
            // AdjustStorage dependency as the skeleton PLY above).
            ExportSimpVisualizeInfo(shape.slab_mesh,
                                    options.outputPrefix + "_simp_visualize_info.json");

            // .mat_typed: same as .ma but each vertex line has its ClusterType
            // name appended, e.g. "v x y z r MS_Seam_Boundary".
            shape.slab_mesh.ExportTypedMA(options.outputPrefix + ".mat_typed");
            std::cout << "  Typed MA exported: " << options.outputPrefix << ".mat_typed\n";

            // PLY with per-vertex colours matching the cluster-type palette.
            shape.slab_mesh.ExportClusterPLY(options.outputPrefix + "_cluster.ply");
            std::cout << "  Cluster PLY exported: " << options.outputPrefix << "_cluster.ply\n";

            // Compute final mesh properties
            shape.slab_mesh.ComputeFacesNormal();
            shape.slab_mesh.ComputeVerticesNormal();
            shape.slab_mesh.ComputeEdgesCone();
            shape.slab_mesh.ComputeFacesSimpleTriangles();

            // Export simplified mesh
            std::cout << "Exporting simplified MA..." << std::endl;
            shape.slab_mesh.Export(options.outputPrefix);
            std::cout << "  Simplified MA exported with prefix: " << options.outputPrefix << std::endl;
            ExportMatAsOff(shape.slab_mesh, options.outputPrefix + "_mat_simplified.off");


        }
    }

    std::cout << std::endl << "Done!" << std::endl;

    return 0;
}
