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

// Export the current (live, non-compact) slab mesh state as an OFF file.
// Vertices are written in their original indexed order; deleted vertices are
// skipped and a compact remapping is built on the fly.
static void ExportMatAsOff(const SlabMesh& sm, const std::string& path)
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
#  include "polyscope/polyscope.h"
#  include "polyscope/surface_mesh.h"
#  include "polyscope/curve_network.h"
#  include "polyscope/point_cloud.h"
#  include "polyscope/pick.h"
#  include "imgui.h"
#  include "QemRejectionViz.h"   // VCG-path QEM rejection overlay (flag-guarded)

// ── Polyscope live-simplification viewer ─────────────────────────────────────

// Builds contiguous arrays from the currently active MAT in `sm`.
// Also fills idx_to_vid: contiguous point-cloud index → original vertex id.
struct MatArrays {
    std::vector<std::array<double,3>>  verts;
    std::vector<std::array<size_t,2>>  edges;
    std::vector<std::array<size_t,2>>  boundary_edges; // edges where both endpoints are boundary
    std::vector<std::array<size_t,2>>  spike_edges;    // edges with at least one T1_spike endpoint
    std::vector<std::array<size_t,3>>  faces;
    std::vector<std::array<size_t,3>>  spike_faces;    // faces with at least one T1_spike vertex
    // Spike-free temp mesh: spike faces removed first, then remaining spike edges removed.
    std::vector<std::array<double,3>>  ns_verts;       // compact vertex list for no-spike mesh
    std::vector<std::array<size_t,2>>  ns_edges;       // no-spike edges (remapped to ns_verts)
    std::vector<std::array<size_t,3>>  ns_faces;       // no-spike faces (remapped to ns_verts)
    std::vector<unsigned>              idx_to_vid;       // index in verts → slab vertex id
    std::vector<unsigned>              idx_to_eid;       // MAT Edges edge slot → slab edge id
    std::vector<unsigned>              idx_to_fid;       // MAT Faces face slot → slab face id
    std::vector<std::array<float,3>>   vert_colors;           // per-vertex color by cluster type
    std::vector<std::array<double,3>>  unknown_ttype_verts;   // positions of T0/T5 vertices only
    // Per-cluster filter clouds: [0]=MS_Sheet, [1]=MS_Seam, [2]=MS_Boundary, [3]=MS_Junction,
    //                            [4]=MS_Sheet_Boundary, [5]=MS_Seam_Boundary, [6]=MS_Junction_Boundary
    std::array<std::vector<std::array<double,3>>, 7> cluster_filter_verts;
    // Vertices marked sharpNotContractable by MarkSharpFeatureVertices().
    std::vector<std::array<double,3>> sharp_verts;
    // Per-edge color for struct membership on "MAT Edges":
    //   green  = struct edge, struct_ids match both endpoints (collapsible)
    //   red    = struct edge, struct_ids mismatch (not collapsible)
    //   grey   = not a struct edge
    std::vector<std::array<float,3>> edge_structure_collapsible_colors;
    // Per-edge color by SlabEdge::topo_type (palette: kEdgeTopoTypeColors).
    // Sourced from the imported MAT structure data, not from face-count
    // geometry — survives collapses via CombineTopoType union.
    std::vector<std::array<float,3>> edge_topo_type_colors;
    // Per-face color by struct_id (grey = no struct).
    std::vector<std::array<float,3>> face_struct_id_colors;
};

// When true, the ImGui edge thickness slider overrides all curve network radii.
// When false, each network keeps its hardcoded radius.
static constexpr bool kModifyGlobalEdgeThickness = true;

// Colors for each ClusterType (index = uint8_t value of the enum).
// T0=invalid: magenta, T1: yellow-green, T2: cyan, T3: orange, T4: red, T5=invalid: white
// MS_* values (7-14) are for the MatStruct external-file path.
static constexpr std::array<std::array<float,3>, 15> kClusterTypeColors = {{
    {0.9f, 0.0f, 0.9f},   // 0  T0                   — invalid (magenta)
    {0.1f, 0.1f, 0.1f},   // 1  T1_spike             — true spike (near-black)
    {0.0f, 0.85f, 1.0f},  // 2  T2                   — 2 clusters (bright cyan)
    {1.0f, 0.5f, 0.0f},   // 3  T3                   — 3 clusters (vivid orange)
    {1.0f, 0.1f, 0.1f},   // 4  T4                   — 4 clusters (vivid red)
    {1.0f, 1.0f, 1.0f},   // 5  T5                   — invalid (white)
    {0.55f, 0.55f, 0.55f},// 6  T1_non_spike         — 1 cluster (grey)
    {0.35f, 0.35f, 0.35f},// 7  MS_Unknown           — unknown (dark grey)
    {0.0f, 1.0f, 0.3f},   // 8  MS_Sheet             — interior sheet (vivid green)
    {1.0f, 0.9f, 0.0f},   // 9  MS_Seam              — internal feature (vivid yellow)
    {0.0f, 0.5f, 1.0f},   // 10 MS_Boundary          — boundary only (vivid blue)
    {1.0f, 0.0f, 0.5f},   // 11 MS_Junction          — junction (vivid pink)
    {0.0f, 0.9f, 1.0f},   // 12 MS_Sheet_Boundary    — sheet+boundary (bright cyan)
    {1.0f, 0.35f, 0.0f},  // 13 MS_Seam_Boundary     — seam+boundary (vivid orange-red)
    {0.6f, 0.0f, 1.0f},   // 14 MS_Junction_Boundary — junction+boundary (vivid purple)
}};

// Colors for each SlabEdge::TopoType (index = uint8_t value of the enum).
// Mirrors the vertex MS_* palette so an edge's topo color visually matches
// the vertex cluster color of the same kind.
//   Unknown:       grey
//   Sheet:         vivid green   (matches MS_Sheet)
//   Seam:          vivid yellow  (matches MS_Seam)
//   Boundary:      vivid blue    (matches MS_Boundary)
//   Seam_Boundary: orange-red    (matches MS_Seam_Boundary)
//   Orphan:        magenta       (visual outlier — should be rare)
static constexpr std::array<std::array<float,3>, 6> kEdgeTopoTypeColors = {{
    {0.55f, 0.55f, 0.55f},  // 0 Unknown
    {0.0f,  1.0f,  0.3f },  // 1 Sheet
    {1.0f,  0.9f,  0.0f },  // 2 Seam
    {0.0f,  0.5f,  1.0f },  // 3 Boundary
    {1.0f,  0.35f, 0.0f },  // 4 Seam_Boundary
    {0.9f,  0.0f,  0.9f },  // 5 Orphan
}};

// MS_* names only (indices 0..7).  To look up a name from a ClusterType enum
// value (where MS_* occupies 7..14 and T-types 0..6), subtract 7: the matstruct
// flow only produces MS_* values so this is the intended indexing.
static constexpr std::array<const char*, 8> kClusterTypeNames = {{
    "MS_Unknown",           // ct_idx 7
    "MS_Sheet",             // ct_idx 8
    "MS_Seam",              // ct_idx 9
    "MS_Boundary",          // ct_idx 10
    "MS_Junction",          // ct_idx 11
    "MS_Sheet_Boundary",    // ct_idx 12
    "MS_Seam_Boundary",     // ct_idx 13
    "MS_Junction_Boundary", // ct_idx 14
}};

// Convert HSV (all in [0,1]) to RGB float3.
static std::array<float,3> HsvToRgb(float h, float s, float v)
{
    float r = 0, g = 0, b = 0;
    int   i = (int)(h * 6.0f);
    float f = h * 6.0f - i;
    float p = v * (1.0f - s);
    float q = v * (1.0f - f * s);
    float t = v * (1.0f - (1.0f - f) * s);
    switch (i % 6) {
        case 0: r=v; g=t; b=p; break;
        case 1: r=q; g=v; b=p; break;
        case 2: r=p; g=v; b=t; break;
        case 3: r=p; g=q; b=v; break;
        case 4: r=t; g=p; b=v; break;
        case 5: r=v; g=p; b=q; break;
    }
    return {r, g, b};
}

// Map a struct_id to a visually distinct colour using the golden-ratio hue step.
static std::array<float,3> StructIdColor(int struct_id)
{
    if (struct_id < 0) return {0.5f, 0.5f, 0.5f};
    const float golden = 0.618033988749895f;
    float hue = std::fmod(struct_id * golden, 1.0f);
    return HsvToRgb(hue, 0.85f, 0.95f);
}

static MatArrays BuildMatArrays(const SlabMesh& sm)
{
    MatArrays out;
    std::map<unsigned,size_t> vid_map;
    for (unsigned i = 0; i < sm.vertices.size(); ++i) {
        if (!sm.vertices[i].first) continue;
        vid_map[i] = out.verts.size();
        const auto& c = sm.vertices[i].second->sphere.center;
        out.verts.push_back({c.X(), c.Y(), c.Z()});
        out.idx_to_vid.push_back(i);
        const auto ct = sm.vertices[i].second->nmn_cluster_type;
        const auto idx = static_cast<uint8_t>(ct);
        out.vert_colors.push_back(kClusterTypeColors[idx < 15 ? idx : 5]);

        // Collect positions of unknown-T-type (T0/T5) vertices separately.
        using CT2 = SlabVertex::ClusterType;
        if (ct == CT2::T0 || ct == CT2::T5)
            out.unknown_ttype_verts.push_back({c.X(), c.Y(), c.Z()});

        // Collect into per-cluster filter arrays for the isolated-type view.
        switch (ct) {
            case CT2::MS_Sheet:            out.cluster_filter_verts[0].push_back({c.X(), c.Y(), c.Z()}); break;
            case CT2::MS_Seam:             out.cluster_filter_verts[1].push_back({c.X(), c.Y(), c.Z()}); break;
            case CT2::MS_Boundary:         out.cluster_filter_verts[2].push_back({c.X(), c.Y(), c.Z()}); break;
            case CT2::MS_Junction:         out.cluster_filter_verts[3].push_back({c.X(), c.Y(), c.Z()}); break;
            case CT2::MS_Sheet_Boundary:   out.cluster_filter_verts[4].push_back({c.X(), c.Y(), c.Z()}); break;
            case CT2::MS_Seam_Boundary:    out.cluster_filter_verts[5].push_back({c.X(), c.Y(), c.Z()}); break;
            case CT2::MS_Junction_Boundary:out.cluster_filter_verts[6].push_back({c.X(), c.Y(), c.Z()}); break;
            default: break;
        }

        // Sharp feature vertices (sharpNotContractable).
        if (sm.vertices[i].second->sharpNotContractable)
            out.sharp_verts.push_back({c.X(), c.Y(), c.Z()});
    }
    using CT = SlabVertex::ClusterType;
    for (unsigned i = 0; i < sm.edges.size(); ++i) {
        if (!sm.edges[i].first) continue;
        size_t a = vid_map.at(sm.edges[i].second->vertices_.first);
        size_t b = vid_map.at(sm.edges[i].second->vertices_.second);
        out.edges.push_back({a, b});
        out.idx_to_eid.push_back(i);
        const SlabVertex* va = sm.vertices[sm.edges[i].second->vertices_.first].second;
        const SlabVertex* vb = sm.vertices[sm.edges[i].second->vertices_.second].second;
        // Structure collapsibility color: grey=not a struct edge, green=collapsible, red=not
        {
            const SlabEdge* se = sm.edges[i].second;
            if (se->struct_ids.empty()) {
                out.edge_structure_collapsible_colors.push_back({0.55f, 0.55f, 0.55f}); // grey: not a struct edge
            } else {
                const SlabVertex* va = sm.vertices[se->vertices_.first].second;
                const SlabVertex* vb = sm.vertices[se->vertices_.second].second;
                bool match = (se->struct_ids == va->struct_ids) && (se->struct_ids == vb->struct_ids);
                out.edge_structure_collapsible_colors.push_back(match
                    ? std::array<float,3>{0.1f, 0.9f, 0.1f}   // green: collapsible
                    : std::array<float,3>{0.9f, 0.1f, 0.1f});  // red: not collapsible
            }
        }
        // Edge topo type color (from imported structure data; see SlabEdge::TopoType).
        {
            const auto tt_idx = static_cast<uint8_t>(sm.edges[i].second->topo_type);
            out.edge_topo_type_colors.push_back(
                kEdgeTopoTypeColors[tt_idx < kEdgeTopoTypeColors.size() ? tt_idx : 0]);
        }
        if (va->topo_is_boundary && vb->topo_is_boundary)
            out.boundary_edges.push_back({a, b});
        if (va->nmn_cluster_type == CT::T1_spike || vb->nmn_cluster_type == CT::T1_spike)
            out.spike_edges.push_back({a, b});
    }

    // Build compact vertex list for no-spike mesh.
    std::unordered_map<size_t,size_t> ns_remap;
    auto ns_vid = [&](size_t vid) -> size_t {
        auto it = ns_remap.find(vid);
        if (it != ns_remap.end()) return it->second;
        size_t idx = out.ns_verts.size();
        ns_remap[vid] = idx;
        out.ns_verts.push_back(out.verts[vid]);
        return idx;
    };

    for (unsigned i = 0; i < sm.faces.size(); ++i) {
        if (!sm.faces[i].first) continue;
        auto it = sm.faces[i].second->vertices_.begin();
        size_t a = vid_map.at(*it++);
        size_t b = vid_map.at(*it++);
        size_t c = vid_map.at(*it);
        out.faces.push_back({a, b, c});
        out.idx_to_fid.push_back(i);
        out.face_struct_id_colors.push_back(StructIdColor(sm.faces[i].second->struct_id));
        bool is_spike_face = false;
        for (unsigned fvid : sm.faces[i].second->vertices_) {
            if (sm.vertices[fvid].second->nmn_cluster_type == CT::T1_spike) {
                is_spike_face = true;
                break;
            }
        }
        if (is_spike_face)
            out.spike_faces.push_back({a, b, c});
        else
            // Step 1: exclude spike faces.
            out.ns_faces.push_back({ns_vid(a), ns_vid(b), ns_vid(c)});
    }

    // Step 2: exclude spike edges (edges where either endpoint is T1_spike).
    for (unsigned i = 0; i < sm.edges.size(); ++i) {
        if (!sm.edges[i].first) continue;
        const SlabVertex* va = sm.vertices[sm.edges[i].second->vertices_.first].second;
        const SlabVertex* vb = sm.vertices[sm.edges[i].second->vertices_.second].second;
        if (va->nmn_cluster_type == CT::T1_spike || vb->nmn_cluster_type == CT::T1_spike) continue;
        size_t a = vid_map.at(sm.edges[i].second->vertices_.first);
        size_t b = vid_map.at(sm.edges[i].second->vertices_.second);
        // Only add if both endpoints already exist in ns_verts (i.e., they appear in ns_faces).
        auto ia = ns_remap.find(a), ib = ns_remap.find(b);
        if (ia == ns_remap.end() || ib == ns_remap.end()) continue;
        out.ns_edges.push_back({ia->second, ib->second});
    }

    return out;
}

// Returns the bounding-box center of the input mesh (used to center all
// input-mesh-derived positions before dividing by bb_diagonal_length).
static std::array<double,3> InputMeshCenter(const MPMesh* pmesh)
{
    return {
        (pmesh->m_min[0] + pmesh->m_max[0]) * 0.5,
        (pmesh->m_min[1] + pmesh->m_max[1]) * 0.5,
        (pmesh->m_min[2] + pmesh->m_max[2]) * 0.5
    };
}

// Distinct colours cycled over clusters.
static const std::array<std::array<float,3>, 8> kClusterPalette = {{
    {0.0f,  1.0f,  0.85f},  // cyan
    {1.0f,  0.45f, 0.1f},   // orange
    {0.4f,  1.0f,  0.3f},   // green
    {1.0f,  0.2f,  0.75f},  // pink
    {1.0f,  1.0f,  0.2f},   // yellow
    {0.55f, 0.3f,  1.0f},   // purple
    {1.0f,  0.2f,  0.2f},   // red
    {0.2f,  0.65f, 1.0f},   // blue
}};

// Names of point clouds the click/clear handlers manage themselves; they are
// excluded from the snapshot so a stale "on" state doesn't get restored.
static constexpr const char* kClickManagedNames[] = {
    "unsimp_mat_crspnd_points",
    "MAT Vert Selected",
};

// Walk every registered Polyscope structure and record (type, name) for those
// currently enabled.  Used to capture the layer-panel state at the instant a
// MAT vertex is clicked so it can be restored when the selection is cleared.
static void SnapshotEnabledPolyscopeStructures(
    std::vector<std::pair<std::string,std::string>>& snapshot)
{
    snapshot.clear();
    for (const auto& [type, named] : polyscope::state::structures) {
        for (const auto& [name, sptr] : named) {
            if (!sptr || !sptr->isEnabled()) continue;
            bool managed = false;
            for (const char* m : kClickManagedNames) {
                if (name == m) { managed = true; break; }
            }
            if (managed) continue;
            snapshot.emplace_back(type, name);
        }
    }
}

// Disable every registered Polyscope structure (any type).
static void DisableAllPolyscopeStructures()
{
    for (const auto& [type, named] : polyscope::state::structures)
        for (const auto& [name, sptr] : named)
            if (sptr) sptr->setEnabled(false);
}

// Re-enable every (type, name) recorded by SnapshotEnabledPolyscopeStructures.
// Silently skips entries whose structure no longer exists (e.g. re-registered
// under a new name during simplification).
static void RestoreEnabledPolyscopeStructures(
    const std::vector<std::pair<std::string,std::string>>& snapshot)
{
    for (const auto& [type, name] : snapshot)
        if (polyscope::hasStructure(type, name))
            polyscope::getStructure(type, name)->setEnabled(true);
}

// Register "unsimp_mat_crspnd_points" — the positions of the original
// (pre-simplification) MAT vertices that collapsed into the selected vertex.
// Positions come from the load-time snapshot sm.original_positions, keyed by
// the original ids stored in SlabVertex::original_ancestors.  Also spawns the
// single-point "MAT Vert Selected" cloud to highlight the chosen vertex.
//
// Side-effect: snapshots the currently-enabled Polyscope structures into
// `enabled_snapshot`, then disables every structure so the click view is
// unambiguous.  The clear-selection / auto-clear handlers call
// RestoreEnabledPolyscopeStructures(enabled_snapshot) to undo this.
static void ShowUnsimpMatCrspndPoints(
    const SlabMesh& sm, unsigned vid,
    std::vector<std::pair<std::string,std::string>>& enabled_snapshot)
{
    if (!sm.vertices[vid].first) return;
    const SlabVertex& sv = *sm.vertices[vid].second;

    std::vector<std::array<double,3>> pts;
    pts.reserve(sv.original_ancestors.size());
    for (unsigned aid : sv.original_ancestors) {
        if (aid < sm.original_positions.size())
            pts.push_back(sm.original_positions[aid]);
    }

    if (pts.empty()) {
        if (polyscope::hasPointCloud("unsimp_mat_crspnd_points"))
            polyscope::getPointCloud("unsimp_mat_crspnd_points")->setEnabled(false);
        return;
    }

    // Capture pre-click scene state, then clear the scene so only the
    // click-related primitives below remain visible.  Only re-snapshot when
    // there's no active selection — otherwise an A-then-B click sequence
    // would overwrite the true pre-A snapshot with the post-A scene state.
    if (enabled_snapshot.empty())
        SnapshotEnabledPolyscopeStructures(enabled_snapshot);
    DisableAllPolyscopeStructures();

    auto* uc = polyscope::registerPointCloud("unsimp_mat_crspnd_points", pts);
    uc->setPointRadius(0.0020, true);
    uc->setPointColor(glm::vec3(0.0f, 1.0f, 0.85f)); // cyan
    uc->setEnabled(true);

    // Spawn a single-point cloud for the selected MAT vertex itself,
    // coloured by its T-type and drawn with a larger radius.
    {
        const auto& c = sv.sphere.center;
        std::vector<std::array<double,3>> mpt = {{ {c.X(), c.Y(), c.Z()} }};
        const auto ct_idx = static_cast<uint8_t>(sv.nmn_cluster_type);
        const auto& col = kClusterTypeColors[ct_idx < 12 ? ct_idx : 5];
        auto* mpc = polyscope::registerPointCloud("MAT Vert Selected", mpt);
        mpc->setPointColor(glm::vec3(col[0], col[1], col[2]));
        mpc->setPointRadius(0.0040, true);
        mpc->setEnabled(true);
    }

    // Show the static initial-MAT face mesh as a spatial reference.
    if (polyscope::hasSurfaceMesh("Initial MAT Faces"))
    {
        auto ps_mesh = polyscope::getSurfaceMesh("Initial MAT Faces");
        ps_mesh->setEnabled(true);
        ps_mesh->setEdgeWidth(1.24f);
    }
}

// Registers the input surface mesh (pmesh) into Polyscope.
// Uses the CGAL Polyhedron face/vertex iterators directly so it works
// regardless of whether GenerateFaceList() was called.
static void RegisterInputMesh(const SlabMesh& sm)
{
    namespace ps = polyscope;
    if (!sm.pmesh) return;

    // Center at origin and normalise: (p - bb_center) / bb_diagonal_length.
    const double inv_diag = 1.0 / sm.pmesh->bb_diagonal_length;
    const auto cm = InputMeshCenter(sm.pmesh);
    std::vector<std::array<double,3>> verts;
    verts.reserve(sm.pmesh->pVertexList.size());
    for (unsigned i = 0; i < sm.pmesh->pVertexList.size(); ++i) {
        const auto& p = sm.pmesh->pVertexList[i]->point();
        verts.push_back({(p[0]-cm[0])*inv_diag, (p[1]-cm[1])*inv_diag, (p[2]-cm[2])*inv_diag});
    }

    // Faces — iterate directly via CGAL halfedge structure.
    // Each facet_begin() halfedge circulates the face vertices in order.
    std::vector<std::array<size_t,3>> faces;
    for (auto fi = sm.pmesh->facets_begin();
         fi != sm.pmesh->facets_end(); ++fi)
    {
        auto h = fi->facet_begin();
        size_t v0 = h->vertex()->id;
        size_t v1 = h->next()->vertex()->id;
        size_t v2 = h->next()->next()->vertex()->id;
        faces.push_back({v0, v1, v2});
    }

    bool en = ps::hasSurfaceMesh("surface mesh")
              ? ps::getSurfaceMesh("surface mesh")->isEnabled() : false;
    auto* mm = ps::registerSurfaceMesh("surface mesh", verts, faces);
    mm->setSurfaceColor(glm::vec3(0.55f, 0.70f, 0.85f));
    mm->setTransparency(0.55f);
    mm->setEnabled(en);

    // ── Feature edges / corners ───────────────────────────────────────────────
    // Sharp edges  → orange curve network
    // Concave edges → purple curve network
    // Corners       → magenta point cloud
    // Each uses preserved enabled state (default off — overlay on request).

    auto registerFeatureEdges = [&](
        const std::set<std::array<int,2>>& edge_set,
        const char* name,
        glm::vec3 color,
        bool default_enabled)
    {
        std::vector<std::array<double,3>> nodes;
        std::vector<std::array<size_t,2>> segs;
        // Collect unique node positions with a compact remap.
        std::unordered_map<int, size_t> id_to_idx;
        for (const auto& e : edge_set)
        {
            for (int vid : e)
            {
                if (id_to_idx.find(vid) == id_to_idx.end())
                {
                    id_to_idx[vid] = nodes.size();
                    const auto& p = sm.pmesh->pVertexList[vid]->point();
                    nodes.push_back({(p[0]-cm[0])*inv_diag, (p[1]-cm[1])*inv_diag, (p[2]-cm[2])*inv_diag});
                }
            }
            segs.push_back({id_to_idx[e[0]], id_to_idx[e[1]]});
        }
        if (nodes.empty()) return;
        bool fen = ps::hasCurveNetwork(name)
                   ? ps::getCurveNetwork(name)->isEnabled() : default_enabled;
        auto* cn = ps::registerCurveNetwork(name, nodes, segs);
        cn->setColor(color);
        cn->setRadius(0.003f);
        cn->setEnabled(fen);
    };

    registerFeatureEdges(sm.sharp_edges,   "Sharp Edges",   glm::vec3(1.0f, 0.6f, 0.0f), false);
    registerFeatureEdges(sm.concave_edges, "Concave Edges", glm::vec3(0.6f, 0.1f, 0.9f), false);

    // Corners → magenta point cloud
    if (!sm.feature_corners.empty())
    {
        std::vector<std::array<double,3>> cpts;
        for (int vid : sm.feature_corners)
        {
            const auto& p = sm.pmesh->pVertexList[vid]->point();
            cpts.push_back({(p[0]-cm[0])*inv_diag, (p[1]-cm[1])*inv_diag, (p[2]-cm[2])*inv_diag});
        }
        bool cen = ps::hasPointCloud("Input Corners")
                   ? ps::getPointCloud("Input Corners")->isEnabled() : false;
        auto* pc = ps::registerPointCloud("Input Corners", cpts);
        pc->setPointColor(glm::vec3(1.0f, 0.0f, 0.8f));
        pc->setPointRadius(0.005f);
        pc->setEnabled(cen);
    }
}

struct ViewerState {
    int  collapse_count = 0;
    bool paused         = true;
    bool step_once      = false;
    int  update_every   = 1;
    // After this many collapses each subsequent collapse adds a delay.
    // Set to -1 to disable (no delay ever).
    int  collapse_delay_after = 3000;
    int  collapse_delay_ms    = 1;
    std::chrono::steady_clock::time_point last_frame =
        std::chrono::steady_clock::now();

    // Maps the contiguous "MAT Verts" point cloud index → original slab vertex id.
    // Rebuilt every time BuildMatArrays is called.
    std::vector<unsigned> idx_to_vid;

    // Maps curve-network edge slot → original slab edge ID.
    // Polyscope's pick handler returns a network-local edge index (0, 1, 2 …) with no
    // knowledge of slab IDs.  We fill this vector in the same loop that appends segments
    // to the "MAT Rejection Edges" network so we can translate the pick back instantly:
    //   slab_eid = rejection_eid_order[pick_edge_slot]
    // Rebuilt every call to UpdateRejectionEdgeColors.
    std::vector<unsigned> rejection_eid_order;

    // Number of nodes registered in the "MAT Rejection Edges" curve network.
    // Polyscope's pick index space for a curve network is flat:
    //   slots [0, rejection_node_count)           → node picks
    //   slots [rejection_node_count, ...)          → edge picks
    // We subtract this offset so that a raw local_idx from getSelection() maps
    // directly into rejection_eid_order[local_idx - rejection_node_count].
    // Rebuilt every call to UpdateRejectionEdgeColors.
    size_t rejection_node_count = 0;

    // Slab edge ID of the currently highlighted rejection edge (-1 = none selected).
    int selected_rejection_eid = -1;

    // Whether to visualize endpoint/target spheres when a rejection edge is selected.
    bool show_rejection_spheres = false;

#if defined(ONLY_USE_QEM_CONDITION_CHECKS)
    // ── VCG-path QEM rejection overlay state (rendered into the SAME Polyscope
    // window by the QemRejectionViz module — not a separate viewer). ──────────
    qemviz::State qem_viz;
#endif

    // Global edge thickness for all curve networks (used when kModifyGlobalEdgeThickness=true).
    float edge_thickness = 0.0008f;

    // Whether to show struct-ID color overlay (seam/boundary/junction structs).
    bool show_struct_colors = false;
    // Whether to show the initial (pre-collapse) MAT struct viz instead of the current one.
    bool show_initial_struct = false;

    // Tracks the enabled state of the "Structure Collapsible" edge color quantity
    // on "MAT Edges" so it survives across rebuilds of the curve network.
    bool struct_collapsible_quantity_enabled = true;

    // Currently selected MAT vertex (-1 = none).
    int selected_vid = -1;

    // Currently selected MAT edge (-1 = none).
    int selected_eid = -1;

    // Maps curve-network edge slot → slab edge ID for "MAT Edges".
    // Rebuilt every time UpdateMatStructures is called.
    std::vector<unsigned> idx_to_eid;
    // Number of nodes in the "MAT Edges" network (used as pick index offset).
    size_t mat_edge_node_count = 0;

    // Currently selected face on "Initial MAT Faces" (-1 = none).
    // Selection is intentionally only supported on the initial (pre-simplification)
    // MAT: live face ids are meaningless after collapses recreate faces.
    // The selected id is also the imported .ma file's face index since
    // LoadMatstructMA places file face elem_id into slab_mesh.faces[elem_id].
    int selected_init_fid = -1;

    // Maps surface-mesh face slot on "Initial MAT Faces" → slab face ID.
    // Polyscope SurfaceMesh pick order is [verts, faces, edges, halfedges]
    // (see surface_mesh.h:377-380), so face_slot = local_idx - init_mat_face_vert_count.
    // Populated once in SetupSimplificationViewer; never overwritten thereafter.
    std::vector<unsigned> init_idx_to_fid;
    // # of verts in "Initial MAT Faces" (pick offset). Set once at setup.
    size_t init_mat_face_vert_count = 0;

    // Output prefix used for naming exported files (set from CLIOptions).
    std::string outputPrefix;

    // ── Vertex color mode ─────────────────────────────────────────────────────
    enum class ColorMode { ClusterType, UnknownTType } color_mode = ColorMode::ClusterType;

    // ── Cluster type filter ───────────────────────────────────────────────────
    // -1 = show all (MAT Verts cloud), 8/9/10/11 = show only MS_Sheet/Seam/Boundary/Junction
    // as an isolated, slightly-larger point cloud.
    int cluster_filter = -1;

    // ── Double-click detection ────────────────────────────────────────────────
    // A double-click is two picks of the same vertex within kDoubleClickMs ms.
    int  last_picked_vid  = -1;
    std::chrono::steady_clock::time_point last_pick_time = {};
    static constexpr int kDoubleClickMs = 300;

    // ── Polyscope scene snapshot ─────────────────────────────────────────────
    // (type, name) of every Polyscope structure that was enabled at the moment
    // a MAT vertex was clicked.  Captured by ShowUnsimpMatCrspndPoints so that
    // when the selection is cleared we can restore the layer-panel state the
    // user had before clicking.  Empty when no vertex selection is active.
    std::vector<std::pair<std::string,std::string>> enabled_snapshot;
};

// Re-registers (or updates) the live MAT structures in Polyscope.
// Also updates vs.idx_to_vid from arr.
// The enabled state of each structure is preserved across updates: if the
static void ApplyGlobalEdgeThickness(float r); // forward declaration

// structure already exists its current enabled flag is read back and restored
// after re-registration, so user toggles in the Polyscope UI survive.
static void UpdateMatStructures(const MatArrays& arr, ViewerState& vs)
{
    namespace ps = polyscope;

    vs.idx_to_vid  = arr.idx_to_vid;
    vs.idx_to_eid  = arr.idx_to_eid;
    vs.mat_edge_node_count = arr.verts.size();

    if (!arr.faces.empty()) {
        // Default: on at first registration; preserved thereafter.
        bool en = ps::hasSurfaceMesh("MAT Faces")
                  ? ps::getSurfaceMesh("MAT Faces")->isEnabled() : true;
        auto* mm = ps::registerSurfaceMesh("MAT Faces", arr.verts, arr.faces);
        mm->setSurfaceColor(glm::vec3(0.9f, 0.6f, 0.2f));
        mm->setTransparency(1.0f);
        if (!arr.face_struct_id_colors.empty()) {
            auto* q = mm->addFaceColorQuantity("Struct ID", arr.face_struct_id_colors);
            q->setEnabled(vs.show_struct_colors);
        }
        mm->setEnabled(en);
    } else if (ps::hasSurfaceMesh("MAT Faces")) {
        ps::removeStructure("MAT Faces");
    }

    if (!arr.edges.empty()) {
        // Default: off at first registration; preserved thereafter.
        bool en = ps::hasCurveNetwork("MAT Edges")
                  ? ps::getCurveNetwork("MAT Edges")->isEnabled() : true;
        auto* cn = ps::registerCurveNetwork("MAT Edges", arr.verts, arr.edges);
        cn->setColor(glm::vec3(1.0f, 0.80f, 0.30f));
        cn->setRadius(0.0008f, true);
        if (!arr.edge_structure_collapsible_colors.empty())
            cn->addEdgeColorQuantity("Structure Collapsible", arr.edge_structure_collapsible_colors)
              ->setEnabled(vs.struct_collapsible_quantity_enabled);
        if (!arr.edge_topo_type_colors.empty())
            cn->addEdgeColorQuantity("Edge Topo Type", arr.edge_topo_type_colors);
        cn->setEnabled(en);
    }

    if (!arr.boundary_edges.empty()) {
        // Build a compact node list with only boundary edge endpoints.
        std::vector<std::array<double,3>> bnodes;
        std::vector<std::array<size_t,2>> bsegs;
        std::unordered_map<size_t,size_t> remap;
        for (const auto& e : arr.boundary_edges) {
            for (size_t vid : e) {
                if (remap.find(vid) == remap.end()) {
                    remap[vid] = bnodes.size();
                    bnodes.push_back(arr.verts[vid]);
                }
            }
            bsegs.push_back({remap[e[0]], remap[e[1]]});
        }
        bool en = ps::hasCurveNetwork("MAT Boundary Edges")
                  ? ps::getCurveNetwork("MAT Boundary Edges")->isEnabled() : false;
        auto* cn = ps::registerCurveNetwork("MAT Boundary Edges", bnodes, bsegs);
        cn->setColor(glm::vec3(1.0f, 0.15f, 0.15f));
        cn->setRadius(0.0015f, true);
        cn->setEnabled(en);
    } else if (ps::hasCurveNetwork("MAT Boundary Edges")) {
        ps::removeStructure("MAT Boundary Edges");
    }

    if (!arr.spike_edges.empty()) {
        std::vector<std::array<double,3>> snodes;
        std::vector<std::array<size_t,2>> ssegs;
        std::unordered_map<size_t,size_t> remap;
        for (const auto& e : arr.spike_edges) {
            for (size_t vid : e) {
                if (remap.find(vid) == remap.end()) {
                    remap[vid] = snodes.size();
                    snodes.push_back(arr.verts[vid]);
                }
            }
            ssegs.push_back({remap[e[0]], remap[e[1]]});
        }
        bool en = ps::hasCurveNetwork("Spike Edges")
                  ? ps::getCurveNetwork("Spike Edges")->isEnabled() : false;
        auto* cn = ps::registerCurveNetwork("Spike Edges", snodes, ssegs);
        cn->setColor(glm::vec3(0.0f, 0.0f, 0.0f));   // black — matches T1 vertex color
        cn->setRadius(0.0012f, true);
        cn->setEnabled(en);
    } else if (ps::hasCurveNetwork("Spike Edges")) {
        ps::removeStructure("Spike Edges");
    }

    if (!arr.spike_faces.empty()) {
        std::vector<std::array<double,3>> sfnodes;
        std::vector<std::array<size_t,3>> sftris;
        std::unordered_map<size_t,size_t> sfremap;
        for (const auto& f : arr.spike_faces) {
            std::array<size_t,3> tri;
            for (int k = 0; k < 3; ++k) {
                size_t vid = f[k];
                if (sfremap.find(vid) == sfremap.end()) {
                    sfremap[vid] = sfnodes.size();
                    sfnodes.push_back(arr.verts[vid]);
                }
                tri[k] = sfremap[vid];
            }
            sftris.push_back(tri);
        }
        bool en = ps::hasSurfaceMesh("Spike Faces")
                  ? ps::getSurfaceMesh("Spike Faces")->isEnabled() : false;
        auto* mm = ps::registerSurfaceMesh("Spike Faces", sfnodes, sftris);
        mm->setSurfaceColor(glm::vec3(0.0f, 0.0f, 0.0f));  // black
        mm->setTransparency(0.45f);
        mm->setEnabled(en);
    } else if (ps::hasSurfaceMesh("Spike Faces")) {
        ps::removeStructure("Spike Faces");
    }

    // Spike-free temp MAT: spike faces removed first, then remaining spike edges removed.
    if (!arr.ns_faces.empty()) {
        bool en = ps::hasSurfaceMesh("MAT Faces (No Spikes)")
                  ? ps::getSurfaceMesh("MAT Faces (No Spikes)")->isEnabled() : false;
        auto* mm = ps::registerSurfaceMesh("MAT Faces (No Spikes)", arr.ns_verts, arr.ns_faces);
        mm->setSurfaceColor(glm::vec3(0.3f, 0.8f, 1.0f));  // light blue
        mm->setTransparency(0.25f);
        mm->setEnabled(en);
    } else if (ps::hasSurfaceMesh("MAT Faces (No Spikes)")) {
        ps::removeStructure("MAT Faces (No Spikes)");
    }
    // if (!arr.ns_edges.empty()) {
    //     bool en = ps::hasCurveNetwork("MAT Edges (No Spikes)")
    //               ? ps::getCurveNetwork("MAT Edges (No Spikes)")->isEnabled() : t;
    //     auto* cn = ps::registerCurveNetwork("MAT Edges (No Spikes)", arr.ns_verts, arr.ns_edges);
    //     cn->setColor(glm::vec3(0.2f, 0.6f, 1.0f));
    //     cn->setRadius(0.0008f, true);
    //     cn->setEnabled(en);
    // } 
    
    else if (ps::hasCurveNetwork("MAT Edges (No Spikes)")) {
        ps::removeStructure("MAT Edges (No Spikes)");
    }

    if (!arr.verts.empty()) {
        // Default: off at first registration; preserved thereafter.
        bool en = ps::hasPointCloud("MAT Verts")
                  ? ps::getPointCloud("MAT Verts")->isEnabled() : true;
        auto* pc = ps::registerPointCloud("MAT Verts", arr.verts);
        pc->setPointRadius(0.00297, true);
        // Hide main cloud when an unknown-only mode is active or a filter is selected.
        using CM = ViewerState::ColorMode;
        bool main_verts_visible = (vs.color_mode == CM::ClusterType);
        bool filter_active = (vs.cluster_filter != -1);
        pc->setEnabled(en && main_verts_visible && !filter_active);
        if (!arr.vert_colors.empty())
            pc->addColorQuantity("Cluster Type", arr.vert_colors)
              ->setEnabled(vs.color_mode == CM::ClusterType);
        if (arr.vert_colors.empty())
            pc->setPointColor(glm::vec3(1.0f, 1.0f, 0.4f));
    }

    // Unknown T-Type cloud — only contains T0/T5 vertices.
    if (!arr.unknown_ttype_verts.empty()) {
        bool en = ps::hasPointCloud("MAT Verts (Unknown TType)")
                  ? ps::getPointCloud("MAT Verts (Unknown TType)")->isEnabled()
                  : (vs.color_mode == ViewerState::ColorMode::UnknownTType);
        auto* upc = ps::registerPointCloud("MAT Verts (Unknown TType)", arr.unknown_ttype_verts);
        upc->setPointColor(glm::vec3(1.0f, 0.0f, 1.0f));  // magenta
        upc->setPointRadius(0.0015, true);
        upc->setEnabled(en && vs.color_mode == ViewerState::ColorMode::UnknownTType);
    } else if (ps::hasPointCloud("MAT Verts (Unknown TType)")) {
        ps::getPointCloud("MAT Verts (Unknown TType)")->setEnabled(false);
    }

    // Sharp feature vertices — always visible, drawn on top in red at a larger radius.
    if (!arr.sharp_verts.empty()) {
        auto* spc = ps::registerPointCloud("Sharp Feature Verts", arr.sharp_verts);
        spc->setPointColor(glm::vec3(1.0f, 0.15f, 0.15f));  // vivid red
        spc->setPointRadius(0.004, true);
        spc->setEnabled(false);
    } else if (ps::hasPointCloud("Sharp Feature Verts")) {
        ps::getPointCloud("Sharp Feature Verts")->setEnabled(false);
    }

    // Per-cluster-type filter clouds: one cloud per MS_* type, slightly larger points.
    // Shown only when vs.cluster_filter matches that type; hidden otherwise.
    static const char* kCFCloudNames[7] = {
        "MAT Verts [MS_Sheet]", "MAT Verts [MS_Seam]",
        "MAT Verts [MS_Boundary]", "MAT Verts [MS_Junction]",
        "MAT Verts [MS_Sheet_Boundary]", "MAT Verts [MS_Seam_Boundary]",
        "MAT Verts [MS_Junction_Boundary]"
    };
    static const int kCFCtIdx[7] = { 8, 9, 10, 11, 12, 13, 14 }; // index into kClusterTypeColors
    for (int fi = 0; fi < 7; ++fi) {
        const auto& fv = arr.cluster_filter_verts[fi];
        const char* name = kCFCloudNames[fi];
        bool want_visible = (vs.cluster_filter == kCFCtIdx[fi]);
        if (!fv.empty()) {
            auto* fpc = ps::registerPointCloud(name, fv);
            const auto& col = kClusterTypeColors[kCFCtIdx[fi]];
            fpc->setPointColor(glm::vec3(col[0], col[1], col[2]));
            fpc->setPointRadius(0.003, true); // slightly bigger than the normal 0.0015
            fpc->setEnabled(want_visible);
        } else if (ps::hasPointCloud(name)) {
            ps::getPointCloud(name)->setEnabled(false);
        }
    }

    if (kModifyGlobalEdgeThickness)
        ApplyGlobalEdgeThickness(vs.edge_thickness);
}

// Register / refresh "MAT Struct Edges" (seam+boundary) and "MAT Struct Verts"
// Applies radius r to every active curve network in the scene.
static void ApplyGlobalEdgeThickness(float r)
{
    namespace ps = polyscope;
    static const char* kNetworks[] = {
        "MAT Edges", "MAT Boundary Edges", "Spike Edges",
        "MAT Struct Edges", "Initial MAT Struct Edges", "Initial MAT Edges",
        "MAT Rejection Edges", "Rejection Edges",
        "Sharp Edges", "Concave Edges",
#if defined(ONLY_USE_QEM_CONDITION_CHECKS)
        "MAT QEM Rejection Edges", "QEM Rejection Edges",
#endif
    };
    for (const char* name : kNetworks) {
        if (!ps::hasCurveNetwork(name)) continue;
        auto* cn = ps::getCurveNetwork(name);
        if (cn->isEnabled())
            cn->setRadius(r, true);
    }
}

// (junction) coloured by struct_id.  Pass enabled=false to hide them.
static void UpdateStructColorVisualization(const SlabMesh& sm, bool enabled)
{
    namespace ps = polyscope;

    if (ps::hasSurfaceMesh("MAT Faces"))
        if (auto* q = ps::getSurfaceMesh("MAT Faces")->getQuantity("Struct ID"))
            q->setEnabled(enabled);

    if (!enabled) {
        if (ps::hasCurveNetwork("MAT Struct Edges"))
            ps::getCurveNetwork("MAT Struct Edges")->setEnabled(false);
        if (ps::hasPointCloud("MAT Struct Verts"))
            ps::getPointCloud("MAT Struct Verts")->setEnabled(false);
        return;
    }

    // ── Seam + boundary edges coloured by struct_id ───────────────────────────
    {
        std::vector<std::array<double,3>> nodes;
        std::vector<std::array<size_t,2>> segs;
        std::vector<std::array<float,3>>  edge_colors;
        std::unordered_map<unsigned,size_t> remap;

        auto addV = [&](unsigned vid) -> size_t {
            auto it = remap.find(vid);
            if (it != remap.end()) return it->second;
            size_t idx = nodes.size();
            remap[vid] = idx;
            const auto& c = sm.vertices[vid].second->sphere.center;
            nodes.push_back({c.X(), c.Y(), c.Z()});
            return idx;
        };

        for (unsigned i = 0; i < (unsigned)sm.edges.size(); ++i) {
            if (!sm.edges[i].first) continue;
            if (sm.edges[i].second->struct_ids.empty()) continue;  // not a named struct edge
            const unsigned v0 = sm.edges[i].second->vertices_.first;
            const unsigned v1 = sm.edges[i].second->vertices_.second;
            if (!sm.vertices[v0].first || !sm.vertices[v1].first) continue;
            segs.push_back({addV(v0), addV(v1)});
            auto col = StructIdColor(*sm.edges[i].second->struct_ids.begin());
            edge_colors.push_back(col);
        }

        if (!segs.empty()) {
            auto* cn = ps::registerCurveNetwork("MAT Struct Edges", nodes, segs);
            cn->addEdgeColorQuantity("Struct ID", edge_colors)->setEnabled(true);
            cn->setRadius(0.003f, true);
            cn->setEnabled(true);
        }
    }

    // ── Junction vertices coloured by struct_id ───────────────────────────────
    {
        std::vector<std::array<double,3>> pts;
        std::vector<std::array<float,3>>  pt_colors;

        for (unsigned i = 0; i < (unsigned)sm.vertices.size(); ++i) {
            if (!sm.vertices[i].first) continue;
            if (sm.vertices[i].second->struct_ids.empty()) continue;
            const auto& c = sm.vertices[i].second->sphere.center;
            pts.push_back({c.X(), c.Y(), c.Z()});
            pt_colors.push_back(StructIdColor(*sm.vertices[i].second->struct_ids.begin()));
        }

        if (!pts.empty()) {
            auto* pc = ps::registerPointCloud("MAT Struct Verts", pts);
            pc->addColorQuantity("Struct ID", pt_colors)->setEnabled(true);
            pc->setPointRadius(0.005, true);
            pc->setEnabled(true);
        }
    }

}

// Captures the initial MAT struct colorization (faces, edges, verts) and
// registers it as persistent Polyscope structures (initially disabled).
// Call once after the initial MAT is registered in SetupSimplificationViewer.
static void RegisterInitialStructViz(const SlabMesh& sm, const MatArrays& init_arrays)
{
    namespace ps = polyscope;

    // Add "Struct ID" face color to "Initial MAT Faces" (already registered).
    if (ps::hasSurfaceMesh("Initial MAT Faces") && !init_arrays.face_struct_id_colors.empty()) {
        auto* mm = ps::getSurfaceMesh("Initial MAT Faces");
        mm->addFaceColorQuantity("Struct ID", init_arrays.face_struct_id_colors)->setEnabled(false);
    }

    // Build struct edge curve network.
    {
        std::vector<std::array<double,3>> nodes;
        std::vector<std::array<size_t,2>> segs;
        std::vector<std::array<float,3>>  edge_colors;
        std::unordered_map<unsigned,size_t> remap;

        auto addV = [&](unsigned vid) -> size_t {
            auto it = remap.find(vid);
            if (it != remap.end()) return it->second;
            size_t idx = nodes.size();
            remap[vid] = idx;
            const auto& c = sm.vertices[vid].second->sphere.center;
            nodes.push_back({c.X(), c.Y(), c.Z()});
            return idx;
        };

        for (unsigned i = 0; i < (unsigned)sm.edges.size(); ++i) {
            if (!sm.edges[i].first) continue;
            if (sm.edges[i].second->struct_ids.empty()) continue;
            const unsigned v0 = sm.edges[i].second->vertices_.first;
            const unsigned v1 = sm.edges[i].second->vertices_.second;
            if (!sm.vertices[v0].first || !sm.vertices[v1].first) continue;
            segs.push_back({addV(v0), addV(v1)});
            edge_colors.push_back(StructIdColor(*sm.edges[i].second->struct_ids.begin()));
        }

        if (!segs.empty()) {
            auto* cn = ps::registerCurveNetwork("Initial MAT Struct Edges", nodes, segs);
            cn->addEdgeColorQuantity("Struct ID", edge_colors)->setEnabled(true);
            cn->setRadius(0.003f, true);
            cn->setEnabled(false);
        }
    }

    // Build struct vert point cloud.
    {
        std::vector<std::array<double,3>> pts;
        std::vector<std::array<float,3>>  pt_colors;

        for (unsigned i = 0; i < (unsigned)sm.vertices.size(); ++i) {
            if (!sm.vertices[i].first) continue;
            if (sm.vertices[i].second->struct_ids.empty()) continue;
            const auto& c = sm.vertices[i].second->sphere.center;
            pts.push_back({c.X(), c.Y(), c.Z()});
            pt_colors.push_back(StructIdColor(*sm.vertices[i].second->struct_ids.begin()));
        }

        if (!pts.empty()) {
            auto* pc = ps::registerPointCloud("Initial MAT Struct Verts", pts);
            pc->addColorQuantity("Struct ID", pt_colors)->setEnabled(true);
            pc->setPointRadius(0.005, true);
            pc->setEnabled(false);
        }
    }
}

// Toggles between initial MAT struct viz and current MAT struct viz.
// show_initial=true  → enable initial structures, disable current struct overlay
// show_initial=false → disable initial structures, enable current struct overlay
static void ApplyInitialStructToggle(bool show_initial, bool show_struct_colors)
{
    namespace ps = polyscope;

    // Initial structures.
    if (ps::hasSurfaceMesh("Initial MAT Faces")) {
        auto* mm = ps::getSurfaceMesh("Initial MAT Faces");
        mm->setEnabled(show_initial);
        if (auto* q = mm->getQuantity("Struct ID"))
            q->setEnabled(show_initial);
    }
    if (ps::hasCurveNetwork("Initial MAT Struct Edges"))
        ps::getCurveNetwork("Initial MAT Struct Edges")->setEnabled(show_initial);
    if (ps::hasPointCloud("Initial MAT Struct Verts"))
        ps::getPointCloud("Initial MAT Struct Verts")->setEnabled(show_initial);

    // Current MAT Faces + Edges: hidden while viewing initial, restored when back.
    if (ps::hasSurfaceMesh("MAT Faces")) {
        auto* mm = ps::getSurfaceMesh("MAT Faces");
        mm->setEnabled(!show_initial);
        if (auto* q = mm->getQuantity("Struct ID"))
            q->setEnabled(!show_initial && show_struct_colors);
    }
    if (ps::hasCurveNetwork("MAT Edges"))
        ps::getCurveNetwork("MAT Edges")->setEnabled(!show_initial);
    // Current struct edge/vert overlays: only when not showing initial.
    bool cur = !show_initial && show_struct_colors;
    if (ps::hasCurveNetwork("MAT Struct Edges"))
        ps::getCurveNetwork("MAT Struct Edges")->setEnabled(cur);
    if (ps::hasPointCloud("MAT Struct Verts"))
        ps::getPointCloud("MAT Struct Verts")->setEnabled(cur);
}

// ─────────────────────────────────────────────────────────────────────────────
// Thin float [0,1] wrapper around the shared uint8 colour table in SlabMesh.h.
static std::array<float,3> RejectionReasonColor(SlabMesh::RejectionReason rr)
{
    auto c = SlabMesh::RejectionReasonColorU8(rr);
    return { c[0] / 255.0f, c[1] / 255.0f, c[2] / 255.0f };
}

// Register / refresh "MAT Rejection Edges" curve network.
// Each active edge is coloured by its last collapse-rejection reason
// (white = never attempted).  Enabled-state is preserved across refreshes.
// Also rebuilds vs.rejection_eid_order and vs.rejection_node_count for pick mapping.
static void UpdateRejectionEdgeColors(const SlabMesh& sm, ViewerState& vs)
{
    namespace ps = polyscope;

    std::map<unsigned, size_t> vid_map;
    std::vector<std::array<double,3>> nodes;
    std::vector<std::array<size_t,2>> segs;
    std::vector<std::array<float,3>>  edge_colors;
    vs.rejection_eid_order.clear();

    for (unsigned i = 0; i < (unsigned)sm.vertices.size(); ++i) {
        if (!sm.vertices[i].first) continue;
        vid_map[i] = nodes.size();
        const auto& c = sm.vertices[i].second->sphere.center;
        nodes.push_back({c.X(), c.Y(), c.Z()});
    }

    for (unsigned i = 0; i < (unsigned)sm.edges.size(); ++i) {
        if (!sm.edges[i].first) continue;
        size_t a = vid_map.at(sm.edges[i].second->vertices_.first);
        size_t b = vid_map.at(sm.edges[i].second->vertices_.second);
        segs.push_back({a, b});
        vs.rejection_eid_order.push_back(i);  // slot j → slab edge ID i

        auto it = sm.edge_last_rejection.find(i);
        edge_colors.push_back(it != sm.edge_last_rejection.end()
            ? RejectionReasonColor(it->second)
            : std::array<float,3>{1.0f, 1.0f, 1.0f});  // white = never attempted
    }

    if (segs.empty()) return;

    vs.rejection_node_count = nodes.size();  // edge picks start at this offset

    bool en = ps::hasCurveNetwork("MAT Rejection Edges")
              ? ps::getCurveNetwork("MAT Rejection Edges")->isEnabled() : false;
    auto* cn = ps::registerCurveNetwork("MAT Rejection Edges", nodes, segs);
    cn->setRadius(0.0010f, true);
    cn->setEnabled(en);
    cn->addEdgeColorQuantity("Rejection Reason", edge_colors)->setEnabled(true);
}

// ─────────────────────────────────────────────────────────────────────────────
// Hide all "Rejection *" overlay structures (called when selection changes).
static void ClearRejectionPrimitives()
{
    namespace ps = polyscope;
    if (ps::hasPointCloud("Rejection Verts"))
        ps::getPointCloud("Rejection Verts")->setEnabled(false);
    if (ps::hasCurveNetwork("Rejection Edges"))
        ps::getCurveNetwork("Rejection Edges")->setEnabled(false);
    if (ps::hasSurfaceMesh("Rejection Faces"))
        ps::getSurfaceMesh("Rejection Faces")->setEnabled(false);
    if (ps::hasPointCloud("Rejection Target"))
        ps::getPointCloud("Rejection Target")->setEnabled(false);
    if (ps::hasSurfaceMesh("Flipped Face Before"))
        ps::getSurfaceMesh("Flipped Face Before")->setEnabled(false);
    if (ps::hasSurfaceMesh("Flipped Face After"))
        ps::getSurfaceMesh("Flipped Face After")->setEnabled(false);
    if (ps::hasPointCloud("Rejection Spheres"))
        ps::getPointCloud("Rejection Spheres")->setEnabled(false);
}

// Visualise the ReasonPrimitives stored for slab edge `eid`.
// Vertices → yellow point cloud   ("Rejection Verts")
// Edges    → yellow curve network  ("Rejection Edges")
// Faces    → yellow translucent mesh ("Rejection Faces")
static void ShowRejectionPrimitives(const SlabMesh& sm, unsigned eid, bool show_spheres = false)
{
    namespace ps = polyscope;

    auto it = sm.edge_reason_primitives.find(eid);
    if (it == sm.edge_reason_primitives.end()) { ClearRejectionPrimitives(); return; }
    const auto& prims = it->second;

    // ── helper: get world position for a slab vertex ─────────────────────────
    auto vertPos = [&](unsigned vid) -> std::array<double,3> {
        if (vid < sm.vertices.size() && sm.vertices[vid].first) {
            const auto& c = sm.vertices[vid].second->sphere.center;
            return {c.X(), c.Y(), c.Z()};
        }
        return {0.0, 0.0, 0.0};
    };

    // ── Vertices ─────────────────────────────────────────────────────────────
    if (!prims.vertices.empty()) {
        std::vector<std::array<double,3>> pts;
        pts.reserve(prims.vertices.size());
        for (unsigned vid : prims.vertices)
            pts.push_back(vertPos(vid));
        auto* pc = ps::registerPointCloud("Rejection Verts", pts);
        pc->setPointColor(glm::vec3(1.0f, 1.0f, 0.0f));  // bright yellow
        pc->setPointRadius(0.006, true);
        pc->setEnabled(true);
    } else if (ps::hasPointCloud("Rejection Verts")) {
        ps::getPointCloud("Rejection Verts")->setEnabled(false);
    }

    // ── Edges ─────────────────────────────────────────────────────────────────
    if (!prims.edges.empty()) {
        std::unordered_map<unsigned,size_t> remap;
        std::vector<std::array<double,3>> nodes;
        std::vector<std::array<size_t,2>> segs;
        auto addV = [&](unsigned vid) -> size_t {
            auto jt = remap.find(vid);
            if (jt != remap.end()) return jt->second;
            size_t idx = nodes.size();
            remap[vid] = idx;
            nodes.push_back(vertPos(vid));
            return idx;
        };
        for (const auto& e : prims.edges)
            segs.push_back({addV(e[0]), addV(e[1])});
        auto* cn = ps::registerCurveNetwork("Rejection Edges", nodes, segs);
        cn->setColor(glm::vec3(1.0f, 1.0f, 0.0f));  // bright yellow
        cn->setRadius(0.004f, true);
        cn->setEnabled(true);
    } else if (ps::hasCurveNetwork("Rejection Edges")) {
        ps::getCurveNetwork("Rejection Edges")->setEnabled(false);
    }

    // ── Faces ─────────────────────────────────────────────────────────────────
    if (!prims.faces.empty()) {
        std::unordered_map<unsigned,size_t> remap;
        std::vector<std::array<double,3>> nodes;
        std::vector<std::array<size_t,3>> tris;
        auto addV = [&](unsigned vid) -> size_t {
            auto jt = remap.find(vid);
            if (jt != remap.end()) return jt->second;
            size_t idx = nodes.size();
            remap[vid] = idx;
            nodes.push_back(vertPos(vid));
            return idx;
        };
        for (const auto& f : prims.faces)
            tris.push_back({addV(f[0]), addV(f[1]), addV(f[2])});
        auto* mm = ps::registerSurfaceMesh("Rejection Faces", nodes, tris);
        mm->setSurfaceColor(glm::vec3(1.0f, 1.0f, 0.0f));  // bright yellow
        mm->setTransparency(0.35f);
        mm->setEnabled(true);
    } else if (ps::hasSurfaceMesh("Rejection Faces")) {
        ps::getSurfaceMesh("Rejection Faces")->setEnabled(false);
    }

    // ── Target vertex ─────────────────────────────────────────────────────────
    if (prims.targ_ver.has_value()) {
        std::vector<std::array<double,3>> pts = { *prims.targ_ver };
        auto* pc = ps::registerPointCloud("Rejection Target", pts);
        pc->setPointColor(glm::vec3(0.0f, 1.0f, 0.0f));  // bright green
        pc->setPointRadius(0.008, true);
        pc->setEnabled(true);
    } else if (ps::hasPointCloud("Rejection Target")) {
        ps::getPointCloud("Rejection Target")->setEnabled(false);
    }

    // ── Flipped face (InversionWouldOccur) ───────────────────────────────────
    if (prims.flipped_face.has_value()) {
        const auto& f = *prims.flipped_face;
        // Before collapse — red
        {
            std::vector<std::array<double,3>> nodes = { f[0][0], f[0][1], f[0][2] };
            std::vector<std::array<size_t,3>> tris  = { {0, 1, 2} };
            auto* mm = ps::registerSurfaceMesh("Flipped Face Before", nodes, tris);
            mm->setSurfaceColor(glm::vec3(1.0f, 0.0f, 0.0f));  // red
            mm->setTransparency(1.0f);
            mm->setEnabled(true);
        }
        // After collapse — orange
        {
            std::vector<std::array<double,3>> nodes = { f[1][0], f[1][1], f[1][2] };
            std::vector<std::array<size_t,3>> tris  = { {0, 1, 2} };
            auto* mm = ps::registerSurfaceMesh("Flipped Face After", nodes, tris);
            mm->setSurfaceColor(glm::vec3(1.0f, 0.55f, 0.0f));  // orange
            mm->setTransparency(1.0f);
            mm->setEnabled(true);
        }
    } else {
        if (ps::hasSurfaceMesh("Flipped Face Before"))
            ps::getSurfaceMesh("Flipped Face Before")->setEnabled(false);
        if (ps::hasSurfaceMesh("Flipped Face After"))
            ps::getSurfaceMesh("Flipped Face After")->setEnabled(false);
    }

    // ── Endpoint & target spheres ─────────────────────────────────────────────
    if (show_spheres && eid < sm.edges.size() && sm.edges[eid].first) {
        const unsigned v1 = sm.edges[eid].second->vertices_.first;
        const unsigned v2 = sm.edges[eid].second->vertices_.second;

        std::vector<std::array<double,3>> pts;
        std::vector<double> radii;

        auto addSphere = [&](unsigned vid) {
            if (vid < sm.vertices.size() && sm.vertices[vid].first) {
                const auto& c = sm.vertices[vid].second->sphere.center;
                pts.push_back({c.X(), c.Y(), c.Z()});
                float radius_ptr = sm.vertices[vid].second->sphere.radius;
                if (radius_ptr < 1e-6f) radius_ptr = 1e-4f; // avoid zero-radius spheres
                radii.push_back(radius_ptr);
            }
        };
        addSphere(v1);
        addSphere(v2);

        // target position (from targ_ver if available)
        if (prims.targ_ver.has_value()) {
            pts.push_back(*prims.targ_ver);
            // target radius comes from the edge's collapsed sphere
            radii.push_back(sm.edges[eid].second->sphere.radius);
        }

        if (!pts.empty()) {
            auto* pc = ps::registerPointCloud("Rejection Spheres", pts);
            pc->setPointColor(glm::vec3(0.4f, 0.8f, 1.0f));  // light blue
            pc->setTransparency(0.5f);
            // use per-point radius quantity so each sphere has its real radius
            pc->addScalarQuantity("radius", radii)->setEnabled(false);
            pc->setPointRadiusQuantity("radius", false);
            pc->setEnabled(true);
        }
    } else if (ps::hasPointCloud("Rejection Spheres")) {
        ps::getPointCloud("Rejection Spheres")->setEnabled(false);
    }
}

// ─────────────────────────────────────────────────────────────────────────────

// Call once before Simplify(). Initialises Polyscope, registers the initial
// MAT, installs the ImGui callback, and wires up sm.on_collapse_cb.
static void SetupSimplificationViewer(SlabMesh& sm, ViewerState& vs)
{
    namespace ps = polyscope;

    ps::init();
    ps::options::programName = "QMAT Simplification Viewer";
    ps::view::bgColor = {0.10f, 0.10f, 0.14f, 1.0f};
    ps::options::groundPlaneMode = ps::GroundPlaneMode::None;

    // Register the input surface mesh so bplist points can be verified
    // visually against it.
    RegisterInputMesh(sm);

    // Register initial MAT (faded, hidden by default — shown via layer panel)
    MatArrays init = BuildMatArrays(sm);
    // Snapshot the face-pick mapping for "Initial MAT Faces". This is the only
    // surface mesh whose face ids carry imported-.ma meaning, so we capture it
    // once here and never refresh it (BuildMatArrays calls during simplification
    // describe the live MAT, where face ids drift away from the imported indices).
    vs.init_idx_to_fid = init.idx_to_fid;
    vs.init_mat_face_vert_count = init.verts.size();
    if (!init.faces.empty()) {
        bool en = ps::hasSurfaceMesh("Initial MAT Faces")
                  ? ps::getSurfaceMesh("Initial MAT Faces")->isEnabled() : false;
        auto* mm = ps::registerSurfaceMesh("Initial MAT Faces", init.verts, init.faces);
        mm->setSurfaceColor(glm::vec3(0.55f, 0.55f, 0.55f));
        mm->setTransparency(1.0f);
        mm->setEnabled(en);
    }
    if (!init.edges.empty()) {
        bool en = ps::hasCurveNetwork("Initial MAT Edges")
                  ? ps::getCurveNetwork("Initial MAT Edges")->isEnabled() : false;
        auto* cn = ps::registerCurveNetwork("Initial MAT Edges", init.verts, init.edges);
        cn->setColor(glm::vec3(0.6f, 0.6f, 0.6f));
        cn->setRadius(0.0005f, true);
        cn->setEnabled(en);
    }

    // Capture initial struct colorization (persists for the lifetime of the viewer).
    RegisterInitialStructViz(sm, init);

    // Live MAT — updated every N collapses
    UpdateMatStructures(init, vs);

#if defined(ONLY_USE_QEM_CONDITION_CHECKS)
    // Seed the QEM rejection overlay so it appears in the structure list from the
    // start (refreshed per collapse and after Simplify by the QemRejectionViz module).
    qemviz::UpdateEdgeColors(sm, vs.qem_viz);
#endif

    // Placeholder: highlighted edge (v1 → v2), shown/updated per collapse
    {
        std::vector<std::array<double,3>>  p = {{0,0,0},{0,0,0}};
        std::vector<std::array<size_t,2>>  e = {{0,1}};
        auto* ce = ps::registerCurveNetwork("Collapsed Edge", p, e);
        ce->setColor(glm::vec3(1.0f, 0.15f, 0.15f));
        ce->setRadius(0.0030f, true);
        ce->setEnabled(false);
    }
    // Placeholder point clouds: v1 (red), v2 (orange), result (green)
    {
        std::vector<std::array<double,3>> p = {{0,0,0}};

        auto* pv1 = ps::registerPointCloud("v1", p);
        pv1->setPointColor(glm::vec3(1.0f, 0.20f, 0.20f));
        pv1->setPointRadius(0.0040, true);
        pv1->setEnabled(false);

        auto* pv2 = ps::registerPointCloud("v2", p);
        pv2->setPointColor(glm::vec3(1.0f, 0.55f, 0.10f));
        pv2->setPointRadius(0.0040, true);
        pv2->setEnabled(false);

        auto* pr = ps::registerPointCloud("result", p);
        pr->setPointColor(glm::vec3(0.20f, 1.0f, 0.30f));
        pr->setPointRadius(0.0050, true);
        pr->setEnabled(false);
    }
    // Placeholder click-selection clouds (populated by ShowUnsimpMatCrspndPoints)
    {
        std::vector<std::array<double,3>> p = {{0,0,0}};

        // For manual pick selection
        auto* bpSel = ps::registerPointCloud("unsimp_mat_crspnd_points", p);
        bpSel->setPointColor(glm::vec3(0.0f, 1.0f, 0.85f));
        bpSel->setPointRadius(0.0020, true);
        bpSel->setEnabled(false);

        // Placeholder for the single selected MAT vertex highlight
        auto* mSel = ps::registerPointCloud("MAT Vert Selected", p);
        mSel->setPointRadius(0.0040, true);
        mSel->setEnabled(false);
    }

    // ImGui panel — assigned to polyscope::state::userCallback
    polyscope::state::userCallback = [&vs, &sm]() {
        // Sync "Structure Collapsible" quantity enabled state from Polyscope UI into vs,
        // so it survives the next curve network rebuild.
        if (polyscope::hasCurveNetwork("MAT Edges")) {
            auto* cn = polyscope::getCurveNetwork("MAT Edges");
            auto* q = cn->getQuantity("Structure Collapsible");
            if (q) vs.struct_collapsible_quantity_enabled = q->isEnabled();
        }

        ImGui::PushItemWidth(230);
        ImGui::Text("QMAT Simplification Viewer");
        ImGui::Separator();
        ImGui::Text("Collapses: %d", vs.collapse_count);
        ImGui::Separator();
        ImGui::Checkbox("Pause", &vs.paused);
        if (vs.paused) {
            ImGui::SameLine();
            if (ImGui::Button("Step")) vs.step_once = true;
        }
        ImGui::Separator();
        ImGui::SliderInt("Update every N", &vs.update_every, 1, 3000);
        ImGui::Separator();
        ImGui::InputInt("Delay after N collapses (-1=off)", &vs.collapse_delay_after);
        if (vs.collapse_delay_after >= 0) {
            ImGui::InputInt("Delay (ms)", &vs.collapse_delay_ms);
            if (vs.collapse_delay_ms < 0) vs.collapse_delay_ms = 0;
        }
        ImGui::Separator();

        // ── pick handling ─────────────────────────────────────────────────────
        if (polyscope::pick::haveSelection()) {
            auto [struct_ptr, local_idx] = polyscope::pick::getSelection();

#if defined(ONLY_USE_QEM_CONDITION_CHECKS)
            // ── QEM rejection-edge pick (handled by the QemRejectionViz module) ──
            if (qemviz::HandlePick(sm, vs.qem_viz, struct_ptr, local_idx)) {
                polyscope::pick::resetSelection();
            } else
#endif
            // ── Rejection-edge pick: click "MAT Rejection Edges" to highlight cause ──
            if (ps::hasCurveNetwork("MAT Rejection Edges") &&
                struct_ptr == ps::getCurveNetwork("MAT Rejection Edges"))
            {
                // Polyscope flat pick index: [0, nNodes) = node, [nNodes, ...) = edge.
                if (local_idx >= vs.rejection_node_count) {
                    size_t edge_slot = local_idx - vs.rejection_node_count;
                    if (edge_slot < vs.rejection_eid_order.size()) {
                        unsigned eid = vs.rejection_eid_order[edge_slot];
                        if ((int)eid != vs.selected_rejection_eid) {
                            vs.selected_rejection_eid = (int)eid;
                            ClearRejectionPrimitives();
                            ShowRejectionPrimitives(sm, eid, vs.show_rejection_spheres);
                        }
                    }
                }
                polyscope::pick::resetSelection();
            }
            // ── MAT vertex pick: click to see bplist ─────────────────────────
            else if (struct_ptr == polyscope::getPointCloud("MAT Verts") &&
                local_idx < vs.idx_to_vid.size())
            {
                unsigned vid = vs.idx_to_vid[local_idx];

                // ── double-click detection ────────────────────────────────────
                using clock = std::chrono::steady_clock;
                auto now = clock::now();
                auto ms  = std::chrono::duration_cast<std::chrono::milliseconds>(
                               now - vs.last_pick_time).count();
                bool is_double_click = ((int)vid == vs.last_picked_vid) &&
                                       (ms < ViewerState::kDoubleClickMs);
                vs.last_picked_vid  = (int)vid;
                vs.last_pick_time   = now;

                if (is_double_click) {
                    // Center and orbit around this vertex, preserving camera distance.
                    const auto& c = sm.vertices[vid].second->sphere.center;
                    glm::vec3 target = {(float)c.X(), (float)c.Y(), (float)c.Z()};
                    glm::vec3 camPos = polyscope::view::getCameraWorldPosition();
                    // Keep the camera at its current position, just repoint at vertex.
                    // This matches MeshLab's behaviour: pivot moves, camera stays.
                    polyscope::view::lookAt(camPos, target, /*flyTo=*/true);
                    polyscope::pick::resetSelection(); // consume the pick
                } else if ((int)vid != vs.selected_vid) {
                    vs.selected_vid           = (int)vid;
                    vs.selected_eid           = -1;   // deselect any edge
                    vs.selected_init_fid      = -1;   // deselect any initial face
                    vs.selected_rejection_eid = -1;
                    ClearRejectionPrimitives();
                    ShowUnsimpMatCrspndPoints(sm, vid, vs.enabled_snapshot);
                }
            }
            // ── MAT edge pick: click to see struct_ids ────────────────────────
            else if (ps::hasCurveNetwork("MAT Edges") &&
                     struct_ptr == ps::getCurveNetwork("MAT Edges"))
            {
                if (local_idx >= vs.mat_edge_node_count) {
                    size_t edge_slot = local_idx - vs.mat_edge_node_count;
                    if (edge_slot < vs.idx_to_eid.size()) {
                        vs.selected_eid = (int)vs.idx_to_eid[edge_slot];
                        vs.selected_vid = -1;       // deselect any vertex
                        vs.selected_init_fid = -1;  // deselect any initial face
                    }
                }
                polyscope::pick::resetSelection();
            }
            // ── Initial-MAT face pick: click to see struct_id / imported face idx ─
            // SurfaceMesh pick layout (surface_mesh.h:377-380):
            //   [0, nVerts) = vertex pick, [nVerts, nVerts+nFaces) = face pick.
            // Selection is restricted to "Initial MAT Faces" because only those
            // face ids correspond to indices in the imported .ma struct file;
            // live "MAT Faces" face ids drift as collapses recreate faces.
            else if (ps::hasSurfaceMesh("Initial MAT Faces") &&
                     struct_ptr == ps::getSurfaceMesh("Initial MAT Faces"))
            {
                if (local_idx >= vs.init_mat_face_vert_count &&
                    local_idx <  vs.init_mat_face_vert_count + vs.init_idx_to_fid.size())
                {
                    size_t face_slot = local_idx - vs.init_mat_face_vert_count;
                    vs.selected_init_fid = (int)vs.init_idx_to_fid[face_slot];
                    vs.selected_vid = -1;   // deselect any vertex
                    vs.selected_eid = -1;   // deselect any edge
                }
                polyscope::pick::resetSelection();
            }
        }

        // Auto-clear selection if the selected vertex was collapsed away.
        if (vs.selected_vid >= 0 &&
            ((unsigned)vs.selected_vid >= sm.vertices.size() ||
             !sm.vertices[vs.selected_vid].first))
        {
            vs.selected_vid = -1;
            polyscope::pick::resetSelection();
            if (polyscope::hasPointCloud("unsimp_mat_crspnd_points"))
                polyscope::getPointCloud("unsimp_mat_crspnd_points")->setEnabled(false);
            if (polyscope::hasPointCloud("MAT Vert Selected"))
                polyscope::getPointCloud("MAT Vert Selected")->setEnabled(false);
            RestoreEnabledPolyscopeStructures(vs.enabled_snapshot);
            vs.enabled_snapshot.clear();
            // Initial MAT Faces is forced off on clear regardless of pre-click state.
            if (polyscope::hasSurfaceMesh("Initial MAT Faces"))
                polyscope::getSurfaceMesh("Initial MAT Faces")->setEnabled(false);
        }

        // Show info about selected vertex
        if (vs.selected_vid >= 0 &&
            (unsigned)vs.selected_vid < sm.vertices.size() &&
            sm.vertices[vs.selected_vid].first)
        {
            const auto& sv = *sm.vertices[vs.selected_vid].second;
            const auto ct_idx = static_cast<uint8_t>(sv.nmn_cluster_type);
            ImGui::Text("Selected vertex: %d", vs.selected_vid);
            ImGui::Text("  T-type: %s",
                (ct_idx >= 7 && ct_idx <= 14) ? kClusterTypeNames[ct_idx - 7] : "MS_Unknown");
            ImGui::Text("  nmn_bplist size: %d", (int)sv.nmn_bplist.size());
            ImGui::Text("  clusters: %d", (int)sv.nmn_bplist_clusters.size());
            ImGui::Text("  (unsimp_mat_crspnd_points: ancestors from initial MAT)");
            if (sv.struct_ids.empty()) {
                ImGui::Text("  struct_ids: (none)");
            } else {
                std::string s;
                for (int id : sv.struct_ids) { if (!s.empty()) s += ", "; s += std::to_string(id); }
                ImGui::Text("  struct_ids: {%s}", s.c_str());
            }
            // Original (pre-simplification) MAT vertex ids that collapsed into this vertex.
            ImGui::Text("  original_ancestors: %d", (int)sv.original_ancestors.size());
            if (!sv.original_ancestors.empty()) {
                constexpr size_t kMaxShown = 32;
                std::string s;
                size_t shown = 0;
                for (unsigned id : sv.original_ancestors) {
                    if (shown >= kMaxShown) { s += ", ..."; break; }
                    if (!s.empty()) s += ", ";
                    s += std::to_string(id);
                    ++shown;
                }
                ImGui::TextWrapped("    {%s}", s.c_str());
            }
            if (ImGui::Button("Clear selection")) {
                vs.selected_vid = -1;
                polyscope::pick::resetSelection();   // prevent next-frame re-trigger
                if (polyscope::hasPointCloud("unsimp_mat_crspnd_points"))
                    polyscope::getPointCloud("unsimp_mat_crspnd_points")->setEnabled(false);
                if (polyscope::hasPointCloud("MAT Vert Selected"))
                    polyscope::getPointCloud("MAT Vert Selected")->setEnabled(false);
                RestoreEnabledPolyscopeStructures(vs.enabled_snapshot);
                vs.enabled_snapshot.clear();
                // Initial MAT Faces is forced off on clear regardless of pre-click state.
                if (polyscope::hasSurfaceMesh("Initial MAT Faces"))
                    polyscope::getSurfaceMesh("Initial MAT Faces")->setEnabled(false);
            }
        } else {
            ImGui::TextDisabled("Click a MAT vertex to see its bplist");
        }

        // ── Selected MAT edge info ────────────────────────────────────────────
        ImGui::Separator();
        if (vs.selected_eid >= 0 &&
            (unsigned)vs.selected_eid < sm.edges.size() &&
            sm.edges[vs.selected_eid].first)
        {
            const SlabEdge& se = *sm.edges[vs.selected_eid].second;
            ImGui::Text("Selected edge: %d", vs.selected_eid);
            ImGui::Text("  v0=%u  v1=%u", se.vertices_.first, se.vertices_.second);
            ImGui::Text("  collapse cost: %.6f", se.collapse_cost);
            ImGui::Text("  topo_type: %s", SlabEdge::TopoTypeName(se.topo_type));
            if (se.struct_ids.empty()) {
                ImGui::Text("  struct_ids: (none)");
            } else {
                std::string s;
                for (int id : se.struct_ids) { if (!s.empty()) s += ", "; s += std::to_string(id); }
                ImGui::Text("  struct_ids: {%s}", s.c_str());
            }
            if (ImGui::Button("Clear edge selection"))
                vs.selected_eid = -1;
        } else {
            ImGui::TextDisabled("Click a MAT edge to see its struct_ids / topo_type");
            if (vs.selected_eid >= 0) vs.selected_eid = -1; // stale: edge was deleted
        }

        // ── Selected initial-MAT face info ────────────────────────────────────
        // Only available on "Initial MAT Faces"; live MAT face ids are not stable
        // across simplification, so the imported .ma face index has no meaning there.
        ImGui::Separator();
        if (vs.selected_init_fid >= 0 &&
            (unsigned)vs.selected_init_fid < sm.faces.size() &&
            sm.faces[vs.selected_init_fid].first)
        {
            const SlabFace& sf = *sm.faces[vs.selected_init_fid].second;
            ImGui::Text("Selected initial face: %d", vs.selected_init_fid);
            ImGui::Text("  init_mat_face_idx: %d", vs.selected_init_fid);
            ImGui::Text("  struct_id: %d", sf.struct_id);
            if (ImGui::Button("Clear face selection"))
                vs.selected_init_fid = -1;
        } else {
            ImGui::TextDisabled("Click an Initial MAT face to see its struct_id");
            if (vs.selected_init_fid >= 0) vs.selected_init_fid = -1; // stale
        }

        // ── Selected rejection edge info ──────────────────────────────────────
        ImGui::Separator();
        if (vs.selected_rejection_eid >= 0) {
            unsigned eid = (unsigned)vs.selected_rejection_eid;
            ImGui::Text("Rejection edge: %u", eid);
            if (eid < sm.edges.size() && sm.edges[eid].first)
                ImGui::Text("  Collapse cost: %.6f", sm.edges[eid].second->collapse_cost);
            auto rit = sm.edge_last_rejection.find(eid);
            if (rit != sm.edge_last_rejection.end()) {
                ImGui::Text("  Reason: %s", SlabMesh::RejectionReasonName(rit->second));
                if (rit->second == SlabMesh::RejectionReason::struct_ids_sets_different &&
                    eid < sm.edges.size() && sm.edges[eid].first)
                {
                    auto setStr = [](const std::set<int>& s) -> std::string {
                        if (s.empty()) return "{}";
                        std::string out = "{";
                        for (int id : s) { if (out.size() > 1) out += ","; out += std::to_string(id); }
                        return out + "}";
                    };
                    const SlabEdge* e = sm.edges[eid].second;
                    ImGui::Text("    edge  struct_ids: %s", setStr(e->struct_ids).c_str());
                    unsigned v0 = e->vertices_.first, v1 = e->vertices_.second;
                    if (v0 < sm.vertices.size() && sm.vertices[v0].first)
                        ImGui::Text("    v%u struct_ids: %s", v0, setStr(sm.vertices[v0].second->struct_ids).c_str());
                    if (v1 < sm.vertices.size() && sm.vertices[v1].first)
                        ImGui::Text("    v%u struct_ids: %s", v1, setStr(sm.vertices[v1].second->struct_ids).c_str());
                }
            }
            auto pit = sm.edge_reason_primitives.find(eid);
            if (pit != sm.edge_reason_primitives.end()) {
                const auto& p = pit->second;
                ImGui::Text("  Cause: %d vert(s)  %d edge(s)  %d face(s)",
                    (int)p.vertices.size(), (int)p.edges.size(), (int)p.faces.size());
                for (const auto& m : p.metrics)
                    ImGui::Text("    %s = %.4f", m.first.c_str(), m.second);
            }
            if (ImGui::Checkbox("Show endpoint/target spheres", &vs.show_rejection_spheres)) {
                // checkbox toggled — refresh the visualization immediately
                ShowRejectionPrimitives(sm, eid, vs.show_rejection_spheres);
            }
            if (ImGui::Button("Clear rejection selection")) {
                vs.selected_rejection_eid = -1;
                ClearRejectionPrimitives();
            }
        } else {
            ImGui::TextDisabled("Click a 'MAT Rejection Edges' edge to see cause");
        }

#if defined(ONLY_USE_QEM_CONDITION_CHECKS)
        // ── QEM rejection panel (VCG path) — rendered by the QemRejectionViz module
        qemviz::DrawPanel(sm, vs.qem_viz);
#endif

        // ── Edge thickness ────────────────────────────────────────────────────
        if (kModifyGlobalEdgeThickness) {
            ImGui::Separator();
            if (ImGui::SliderFloat("Edge Thickness", &vs.edge_thickness, 0.0001f, 0.005f, "%.4f"))
                ApplyGlobalEdgeThickness(vs.edge_thickness);
        }

        // ── Struct color overlay ──────────────────────────────────────────────
        ImGui::Separator();
        if (ImGui::Checkbox("Show struct colors (seam/boundary/junction)", &vs.show_struct_colors)) {
            vs.show_initial_struct = false;
            UpdateStructColorVisualization(sm, vs.show_struct_colors);
            ApplyInitialStructToggle(false, vs.show_struct_colors);
        }
        if (vs.show_struct_colors || vs.show_initial_struct) {
            ImGui::SameLine();
            const char* lbl = vs.show_initial_struct ? "Current MAT##structtog" : "Initial MAT##structtog";
            if (ImGui::Button(lbl)) {
                vs.show_initial_struct = !vs.show_initial_struct;
                ApplyInitialStructToggle(vs.show_initial_struct, vs.show_struct_colors);
                UpdateStructColorVisualization(sm, !vs.show_initial_struct && vs.show_struct_colors);
            }
        }

        // ── Vertex color mode toggle ──────────────────────────────────────────
        ImGui::Separator();
        ImGui::Text("MAT Vertex Coloring:");
        using CM = ViewerState::ColorMode;
        bool use_cluster      = vs.color_mode == CM::ClusterType;
        bool use_unk_ttype    = vs.color_mode == CM::UnknownTType;

        // Helper: switch main MAT Verts cloud on/off and update active quantity.
        auto setColorMode = [&](CM new_mode) {
            vs.color_mode = new_mode;
            // Keep MAT Verts hidden when a cluster filter is active.
            bool show_main = (new_mode == CM::ClusterType)
                             && (vs.cluster_filter == -1);
            if (polyscope::hasPointCloud("MAT Verts")) {
                auto* pc = polyscope::getPointCloud("MAT Verts");
                pc->setEnabled(show_main);
                if (show_main) {
                    auto* q_ct = pc->getQuantity("Cluster Type");
                    if (q_ct) q_ct->setEnabled(new_mode == CM::ClusterType);
                }
            }
            if (polyscope::hasPointCloud("MAT Verts (Unknown TType)"))
                polyscope::getPointCloud("MAT Verts (Unknown TType)")->setEnabled(new_mode == CM::UnknownTType);
        };

        ImGui::SameLine();
        if (ImGui::RadioButton("Cluster Type##cm", use_cluster))  setColorMode(CM::ClusterType);
        ImGui::SameLine();
        if (ImGui::RadioButton("Unknown T-Type##cm", use_unk_ttype)) setColorMode(CM::UnknownTType);

        // Legend for cluster type mode
        if (use_cluster) {
            for (int k = 0; k < 8; ++k) {
                const auto& col = kClusterTypeColors[7+k];
                ImGui::TextColored(ImVec4(col[0], col[1], col[2], 1.0f), "  %s", kClusterTypeNames[k]);
            }
        }

        // ── Cluster type filter ───────────────────────────────────────────────
        // Shows only vertices of the selected MS_* type as a separate, slightly-
        // larger point cloud (same colour).  "All" restores the normal MAT Verts.
        ImGui::Separator();
        ImGui::Text("Cluster Type Filter:");
        static const char* kCFLabels[7] = {
            "MS_Sheet", "MS_Seam", "MS_Boundary", "MS_Junction",
            "MS_Sheet_Boundary", "MS_Seam_Boundary", "MS_Junction_Boundary"
        };
        static const char* kCFNames[7]  = {
            "MAT Verts [MS_Sheet]", "MAT Verts [MS_Seam]",
            "MAT Verts [MS_Boundary]", "MAT Verts [MS_Junction]",
            "MAT Verts [MS_Sheet_Boundary]", "MAT Verts [MS_Seam_Boundary]",
            "MAT Verts [MS_Junction_Boundary]"
        };
        static const int kCFVals[7] = { 8, 9, 10, 11, 12, 13, 14 };

        auto applyFilter = [&](int new_filter) {
            vs.cluster_filter = new_filter;
            bool showing_all = (new_filter == -1);
            // Main MAT Verts: visible only when no filter and colour mode supports it.
            if (polyscope::hasPointCloud("MAT Verts"))
                polyscope::getPointCloud("MAT Verts")->setEnabled(
                    showing_all && (vs.color_mode == CM::ClusterType));
            // Filter clouds: exactly one visible (or none when showing all).
            for (int j = 0; j < 7; ++j)
                if (polyscope::hasPointCloud(kCFNames[j]))
                    polyscope::getPointCloud(kCFNames[j])->setEnabled(
                        !showing_all && (kCFVals[j] == new_filter));
        };

        if (ImGui::RadioButton("All##cf", vs.cluster_filter == -1))
            applyFilter(-1);
        for (int fi = 0; fi < 7; ++fi) {
            ImGui::SameLine();
            char lbl[32];
            std::snprintf(lbl, sizeof(lbl), "%s##cf", kCFLabels[fi]);
            if (ImGui::RadioButton(lbl, vs.cluster_filter == kCFVals[fi]))
                applyFilter(kCFVals[fi]);
        }

        ImGui::Separator();
        if (ImGui::Button("Export MAT as OFF")) {
            std::string path = vs.outputPrefix
                + "_snapshot_" + std::to_string(vs.collapse_count) + ".off";
            ExportMatAsOff(sm, path);
        }

        ImGui::Separator();
        ImGui::TextDisabled("Left panel layers:");
        ImGui::TextDisabled("  MAT Faces/Edges – orange/yellow");
        ImGui::TextDisabled("  Collapsed Edge  – red");
        ImGui::TextDisabled("  v1/v2 bplist    – blue/orange dots");
        ImGui::TextDisabled("  unsimp_mat_crspnd_points – initial-MAT ancestors");
        ImGui::PopItemWidth();
    };

    // Warm up the GLFW/WGL context: render a few frames so the window is fully
    // initialized before we start calling frameTick() from the collapse loop.
    for (int i = 0; i < 5; ++i)
        polyscope::frameTick();

    // Per-collapse callback — fired just before each MergeVertices call
    sm.on_collapse_cb = [&vs, &sm](
        unsigned raw_v1, const Wm4::Vector3d& v1p, double /*v1r*/,
        unsigned raw_v2, const Wm4::Vector3d& v2p, double /*v2r*/,
        const Sphere& result)
    {
        namespace ps = polyscope;
        using clock  = std::chrono::steady_clock;

        ++vs.collapse_count;

        // ── stop if the user closed the window ───────────────────────────────
        if (ps::windowRequestsClose()) {
            sm.on_collapse_cb = nullptr;   // disable further callbacks
            return;
        }

        // ── pause / step ─────────────────────────────────────────────────────
        // Track whether the previous collapse was a manual step so we can
        // force a rebuild on this callback (which now fires post-prev-merge,
        // meaning vid_tgt from the previous step is already in the scene).
        static bool prev_was_manual_step = false;
        bool force_update = prev_was_manual_step;
        prev_was_manual_step = false;

        if (vs.paused) {
            while (vs.paused && !vs.step_once && !ps::windowRequestsClose()) {
                std::this_thread::sleep_for(std::chrono::milliseconds(16));
                ps::frameTick();
            }
            if (vs.step_once) {
                vs.step_once = false;
                prev_was_manual_step = true; // rebuild next callback to show result
            }
        }

        // ── per-collapse delay after N collapses ─────────────────────────────
        // Once collapse_count reaches the threshold, every subsequent collapse
        // forces a visual rebuild and sleeps for collapse_delay_ms milliseconds.
        // Set collapse_delay_after = -1 to disable.
        if (vs.collapse_delay_after >= 0 && vs.collapse_count == vs.collapse_delay_after) {
            vs.paused = true; // auto-pause when threshold is reached, so user can inspect
        }
        if (vs.collapse_delay_after >= 0 && vs.collapse_count >= vs.collapse_delay_after) {
            force_update = true;
            std::this_thread::sleep_for(std::chrono::milliseconds(vs.collapse_delay_ms));
            ps::frameTick();
        }

        // ── rate-limit: skip render if too soon and no structural update due ─
        static const auto kMinInterval = std::chrono::milliseconds(33); // ~30 fps
        bool do_update = (vs.collapse_count % vs.update_every == 0) || force_update;
        bool time_ok   = (clock::now() - vs.last_frame) >= kMinInterval;

        if (!do_update && !time_ok)
            return;

        vs.last_frame = clock::now();

        if (do_update) {
            // Rebuild live MAT (v1 + v2 still active at this point — pre-merge)
            UpdateRejectionEdgeColors(sm, vs);
#if defined(ONLY_USE_QEM_CONDITION_CHECKS)
            qemviz::UpdateEdgeColors(sm, vs.qem_viz);
#endif
            UpdateMatStructures(BuildMatArrays(sm), vs);
            UpdateStructColorVisualization(sm, !vs.show_initial_struct && vs.show_struct_colors);

            // Highlighted edge: v1 → v2
            std::vector<std::array<double,3>>  ep = {
                {v1p.X(), v1p.Y(), v1p.Z()},
                {v2p.X(), v2p.Y(), v2p.Z()}
            };
            std::vector<std::array<size_t,2>> ec = {{0, 1}};
            auto* ce = ps::registerCurveNetwork("Collapsed Edge", ep, ec);
            auto collapsed_edge_color = glm::vec3(67, 26, 250) / 255.0f;
            ce->setColor(collapsed_edge_color);
            ce->setRadius(0.0030f, true);
            ce->setEnabled(true);

            // v1, v2, result sphere-centre markers
            // std::vector<std::array<double,3>> p1 = {{v1p.X(), v1p.Y(), v1p.Z()}};
            // std::vector<std::array<double,3>> p2 = {{v2p.X(), v2p.Y(), v2p.Z()}};
            // std::vector<std::array<double,3>> pr = {{result.center.X(),
            //                                          result.center.Y(),
            //                                          result.center.Z()}};
            // ps::getPointCloud("v1")->updatePointPositions(p1);
            // ps::getPointCloud("v1")->setEnabled(true);
            // ps::getPointCloud("v2")->updatePointPositions(p2);
            // ps::getPointCloud("v2")->setEnabled(true);
            // ps::getPointCloud("result")->updatePointPositions(pr);
            // ps::getPointCloud("result")->setEnabled(true);

        }

        ps::frameTick();
    };
}
#endif  // QMAT_WITH_POLYSCOPE

// Simple command line argument parsing
struct CLIOptions {
    std::string inputFile;
    std::string outputPrefix;
    int simplifyTarget = -1;  // -1 means no simplification
    int faceTarget = -1;      // --nf: QEM stops at this MAT face count; -1 = unset (refines --simplify)
    double k = 0.00001;
    double featureAngleDeg = 10.0; // turning-angle threshold for sharp feature protection
    // VCG-faithful QEM thresholds (ONLY_USE_QEM_CONDITION_CHECKS path).
    // Defaults mirror VcgQuadricParameter's struct defaults.
    double qemQualityThr = 0.3;     // --qem-quality-thr  → VcgQuadricParameter::QualityThr
    double qemBoundaryWeight = 1.0; // --qem-boundary-weight = MeshLab UI "Boundary
                                    // Preserving Weight" slider; effective weight is
                                    // 0.5 * this (meshfilter.cpp:14).  Default 1.0 → 0.5.
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
        if (options.simplifyTarget >= currentVertices) {
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
            ViewerState vs;
            vs.outputPrefix = options.outputPrefix;
            if (options.visualize)
                SetupSimplificationViewer(shape.slab_mesh, vs);
#endif
            shape.slab_mesh.allow_steep_collapse = true;
            shape.slab_mesh.Simplify(reductionCount);
            long simplifyTime = clock() - startTime;

            std::cout << "  Simplification time: " << simplifyTime << " ms" << std::endl;
            std::cout << "  Final vertex count: " << shape.slab_mesh.numVertices << std::endl;

#ifdef QMAT_WITH_POLYSCOPE
            if (options.visualize) {
                // Clear the per-collapse callback so frameTick is no longer
                // driven by a running Simplify loop.
                shape.slab_mesh.on_collapse_cb = nullptr;
                // Register the final simplified MAT and hand control to the
                // Polyscope window for interactive inspection.
                UpdateRejectionEdgeColors(shape.slab_mesh, vs);
#if defined(ONLY_USE_QEM_CONDITION_CHECKS)
                qemviz::UpdateEdgeColors(shape.slab_mesh, vs.qem_viz);
#endif
                UpdateMatStructures(BuildMatArrays(shape.slab_mesh), vs);
                UpdateStructColorVisualization(shape.slab_mesh, !vs.show_initial_struct && vs.show_struct_colors);
                if (polyscope::hasCurveNetwork("Collapsed Edge"))
                    polyscope::getCurveNetwork("Collapsed Edge")->setEnabled(false);
                // Keep the original interactive callback intact so the user can
                // still click vertices, export OFF snapshots, etc.
                std::cout << "Simplification done. Close the viewer window to exit.\n";
                polyscope::show();
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
