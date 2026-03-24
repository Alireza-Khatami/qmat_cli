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

#ifdef QMAT_WITH_POLYSCOPE
#  include "polyscope/polyscope.h"
#  include "polyscope/surface_mesh.h"
#  include "polyscope/curve_network.h"
#  include "polyscope/point_cloud.h"
#  include "polyscope/pick.h"
#  include "imgui.h"

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
    std::vector<std::array<float,3>>   vert_colors;      // per-vertex color by cluster type
    std::vector<std::array<float,3>>   topo_vert_colors; // per-vertex color by topo type
};

// Colors for each ClusterType (index = uint8_t value of the enum).
// T0=invalid: magenta, T1: yellow-green, T2: cyan, T3: orange, T4: red, T5=invalid: white
static constexpr std::array<std::array<float,3>, 7> kClusterTypeColors = {{
    {1.0f, 0.0f, 1.0f},   // T0           — invalid (magenta)
    {0.0f, 0.0f, 0.0f},   // T1_spike     — true spike, 4 bpoints (black)
    {0.2f, 0.8f, 1.0f},   // T2           — 2 clusters / sheet (cyan-blue)
    {1.0f, 0.55f, 0.1f},  // T3           — 3 clusters / seam (orange)
    {1.0f, 0.15f, 0.15f}, // T4           — 4 clusters / junction (red)
    {1.0f, 1.0f, 1.0f},   // T5           — invalid (white)
    {0.4f, 0.4f, 0.4f},   // T1_non_spike — 1 cluster, >4 bpoints (grey)
}};

static constexpr std::array<const char*, 7> kClusterTypeNames = {{
    "T0 (invalid)",
    "T1_spike (true spike)",
    "T2 (sheet)",
    "T3 (seam)",
    "T4 (junction)",
    "T5 (invalid)",
    "T1_non_spike (boundary)",
}};

// Colors for each TopoType (index = uint8_t value of the enum, 0–4).
static constexpr std::array<std::array<float,3>, 5> kTopoTypeColors = {{
    {0.5f, 0.5f, 0.5f},   // 0 Unknown  — grey
    {0.2f, 0.8f, 1.0f},   // 1 Sheet    — cyan
    {0.0f, 0.35f, 1.0f},  // 2 Boundary — dark blue
    {1.0f, 0.75f, 0.1f},  // 3 Seam     — yellow-orange
    {1.0f, 0.15f, 0.15f}, // 4 Junction — red
}};

static constexpr std::array<const char*, 5> kTopoTypeNames = {{
    "Unknown",
    "Sheet",
    "Boundary",
    "Seam",
    "Junction",
}};

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
        out.vert_colors.push_back(kClusterTypeColors[idx < 7 ? idx : 5]);
        const auto tt = sm.vertices[i].second->topo_type;
        const auto tidx = static_cast<uint8_t>(tt);
        out.topo_vert_colors.push_back(kTopoTypeColors[tidx < 5 ? tidx : 0]);
    }
    using CT = SlabVertex::ClusterType;
    for (unsigned i = 0; i < sm.edges.size(); ++i) {
        if (!sm.edges[i].first) continue;
        size_t a = vid_map.at(sm.edges[i].second->vertices_.first);
        size_t b = vid_map.at(sm.edges[i].second->vertices_.second);
        out.edges.push_back({a, b});
        const SlabVertex* va = sm.vertices[sm.edges[i].second->vertices_.first].second;
        const SlabVertex* vb = sm.vertices[sm.edges[i].second->vertices_.second].second;
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

// Get the 3-D positions of all nmn_bplist entries for a slab vertex.
// Points come from the input mesh (sm.pmesh->pVertexList).
static std::vector<std::array<double,3>> BplistPositions(
    const SlabMesh& sm, unsigned vid)
{
    std::vector<std::array<double,3>> pts;
    if (!sm.vertices[vid].first) return pts;
    const auto& bps = sm.vertices[vid].second->nmn_bplist;
    pts.reserve(bps.size());
    // MAT vertices are stored normalized by bb_diagonal_length (see LoadInputNMM).
    // Divide bp positions by the same factor so they land in the same space.
    const double inv_diag = 1.0 / sm.pmesh->bb_diagonal_length;
    for (unsigned bp : bps) {
        const auto& p = sm.pmesh->pVertexList[bp]->point();
        pts.push_back({p[0] * inv_diag, p[1] * inv_diag, p[2] * inv_diag});
    }
    return pts;
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

// Register "BPList selected" with per-cluster colours.
// Falls back to solid cyan when no cluster data is available.
static void ShowBplistClusters(const SlabMesh& sm, unsigned vid)
{
    if (!sm.vertices[vid].first) return;
    const SlabVertex& sv = *sm.vertices[vid].second;
    const double inv_diag = 1.0 / sm.pmesh->bb_diagonal_length;

    std::vector<std::array<double,3>> pts;
    std::vector<std::array<float,3>>  colors;

    const auto& clusters = sv.nmn_bplist_clusters;
    if (!clusters.empty()) {
        for (unsigned ci = 0; ci < (unsigned)clusters.size(); ++ci) {
            const auto& col = kClusterPalette[ci % kClusterPalette.size()];
            for (unsigned bp : clusters[ci]) {
                const auto& p = sm.pmesh->pVertexList[bp]->point();
                pts.push_back({p[0]*inv_diag, p[1]*inv_diag, p[2]*inv_diag});
                colors.push_back(col);
            }
        }
    } else {
        // No cluster info — fall back to solid cyan
        for (unsigned bp : sv.nmn_bplist) {
            const auto& p = sm.pmesh->pVertexList[bp]->point();
            pts.push_back({p[0]*inv_diag, p[1]*inv_diag, p[2]*inv_diag});
            colors.push_back({0.0f, 1.0f, 0.85f});
        }
    }

    if (pts.empty()) {
        polyscope::getPointCloud("BPList selected")->setEnabled(false);
        return;
    }

    auto* bpSel = polyscope::registerPointCloud("BPList selected", pts);
    bpSel->setPointRadius(0.0020, true);
    bpSel->addColorQuantity("cluster", colors)->setEnabled(true);
    bpSel->setEnabled(true);

    // Spawn a single-point cloud for the selected MAT vertex itself,
    // coloured by its T-type and drawn with a larger radius.
    {
        const auto& c = sv.sphere.center;
        std::vector<std::array<double,3>> mpt = {{ {c.X(), c.Y(), c.Z()} }};
        const auto ct_idx = static_cast<uint8_t>(sv.nmn_cluster_type);
        const auto& col = kClusterTypeColors[ct_idx < 7 ? ct_idx : 5];
        auto* mpc = polyscope::registerPointCloud("MAT Vert Selected", mpt);
        mpc->setPointColor(glm::vec3(col[0], col[1], col[2]));
        mpc->setPointRadius(0.0040, true);
        mpc->setEnabled(true);
    }

    // Hide the full MAT Verts cloud so the selected point stands out.
    if (polyscope::hasPointCloud("MAT Verts"))
        polyscope::getPointCloud("MAT Verts")->setEnabled(false);
    if (polyscope::hasSurfaceMesh("Input Mesh"))
        {
            auto ps_mesh = polyscope::getSurfaceMesh("Input Mesh");
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

    // Vertices — normalize by bb_diagonal_length to match the MAT coordinate space.
    const double inv_diag = 1.0 / sm.pmesh->bb_diagonal_length;
    std::vector<std::array<double,3>> verts;
    verts.reserve(sm.pmesh->pVertexList.size());
    for (unsigned i = 0; i < sm.pmesh->pVertexList.size(); ++i) {
        const auto& p = sm.pmesh->pVertexList[i]->point();
        verts.push_back({p[0] * inv_diag, p[1] * inv_diag, p[2] * inv_diag});
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

    bool en = ps::hasSurfaceMesh("Input Mesh")
              ? ps::getSurfaceMesh("Input Mesh")->isEnabled() : false;
    auto* mm = ps::registerSurfaceMesh("Input Mesh", verts, faces);
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
                    nodes.push_back({p[0]*inv_diag, p[1]*inv_diag, p[2]*inv_diag});
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
            cpts.push_back({p[0]*inv_diag, p[1]*inv_diag, p[2]*inv_diag});
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
    bool paused         = false;
    bool step_once      = false;
    int  update_every   = 1;
    std::chrono::steady_clock::time_point last_frame =
        std::chrono::steady_clock::now();

    // Maps the contiguous "MAT Verts" point cloud index → original slab vertex id.
    // Rebuilt every time BuildMatArrays is called.
    std::vector<unsigned> idx_to_vid;

    // Currently selected MAT vertex (-1 = none).
    int selected_vid = -1;

    // Output prefix used for naming exported files (set from CLIOptions).
    std::string outputPrefix;

    // ── History viewer state ──────────────────────────────────────────────────
    // When >= 0 the UI is showing a historical keyframe instead of the live mesh.
    int  history_step        = -1;   // -1 = live mode
    int  lineage_cursor      = -1;   // index into lineage records for selected_vid
    bool show_ancestry       = false; // whether to visualise the ancestry cloud

    // ── Vertex color mode ─────────────────────────────────────────────────────
    enum class ColorMode { ClusterType, TopoType } color_mode = ColorMode::ClusterType;

    // ── Double-click detection ────────────────────────────────────────────────
    // A double-click is two picks of the same vertex within kDoubleClickMs ms.
    int  last_picked_vid  = -1;
    std::chrono::steady_clock::time_point last_pick_time = {};
    static constexpr int kDoubleClickMs = 300;
};

// Re-registers (or updates) the live MAT structures in Polyscope.
// Also updates vs.idx_to_vid from arr.
// The enabled state of each structure is preserved across updates: if the
// structure already exists its current enabled flag is read back and restored
// after re-registration, so user toggles in the Polyscope UI survive.
static void UpdateMatStructures(const MatArrays& arr, ViewerState& vs)
{
    namespace ps = polyscope;

    vs.idx_to_vid = arr.idx_to_vid;

    if (!arr.faces.empty()) {
        // Default: on at first registration; preserved thereafter.
        bool en = ps::hasSurfaceMesh("MAT Faces")
                  ? ps::getSurfaceMesh("MAT Faces")->isEnabled() : true;
        auto* mm = ps::registerSurfaceMesh("MAT Faces", arr.verts, arr.faces);
        mm->setSurfaceColor(glm::vec3(0.9f, 0.6f, 0.2f));
        mm->setTransparency(0.35f);
        mm->setEnabled(en);
    } else if (ps::hasSurfaceMesh("MAT Faces")) {
        ps::removeStructure("MAT Faces");
    }

    if (!arr.edges.empty()) {
        // Default: off at first registration; preserved thereafter.
        bool en = ps::hasCurveNetwork("MAT Edges")
                  ? ps::getCurveNetwork("MAT Edges")->isEnabled() : false;
        auto* cn = ps::registerCurveNetwork("MAT Edges", arr.verts, arr.edges);
        cn->setColor(glm::vec3(1.0f, 0.80f, 0.30f));
        cn->setRadius(0.0008f, true);
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
                  ? ps::getSurfaceMesh("MAT Faces (No Spikes)")->isEnabled() : true;
        auto* mm = ps::registerSurfaceMesh("MAT Faces (No Spikes)", arr.ns_verts, arr.ns_faces);
        mm->setSurfaceColor(glm::vec3(0.3f, 0.8f, 1.0f));  // light blue
        mm->setTransparency(0.25f);
        mm->setEnabled(en);
    } else if (ps::hasSurfaceMesh("MAT Faces (No Spikes)")) {
        ps::removeStructure("MAT Faces (No Spikes)");
    }
    if (!arr.ns_edges.empty()) {
        bool en = ps::hasCurveNetwork("MAT Edges (No Spikes)")
                  ? ps::getCurveNetwork("MAT Edges (No Spikes)")->isEnabled() : false;
        auto* cn = ps::registerCurveNetwork("MAT Edges (No Spikes)", arr.ns_verts, arr.ns_edges);
        cn->setColor(glm::vec3(0.2f, 0.6f, 1.0f));
        cn->setRadius(0.0008f, true);
        cn->setEnabled(en);
    } else if (ps::hasCurveNetwork("MAT Edges (No Spikes)")) {
        ps::removeStructure("MAT Edges (No Spikes)");
    }

    if (!arr.verts.empty()) {
        // Default: off at first registration; preserved thereafter.
        bool en = ps::hasPointCloud("MAT Verts")
                  ? ps::getPointCloud("MAT Verts")->isEnabled() : false;
        auto* pc = ps::registerPointCloud("MAT Verts", arr.verts);
        pc->setPointRadius(0.0015, true);
        pc->setEnabled(en);
        // Register both color quantities; enable the one matching color_mode.
        if (!arr.vert_colors.empty())
            pc->addColorQuantity("Cluster Type", arr.vert_colors)
              ->setEnabled(vs.color_mode == ViewerState::ColorMode::ClusterType);
        if (!arr.topo_vert_colors.empty())
            pc->addColorQuantity("Topo Type", arr.topo_vert_colors)
              ->setEnabled(vs.color_mode == ViewerState::ColorMode::TopoType);
        if (arr.vert_colors.empty() && arr.topo_vert_colors.empty())
            pc->setPointColor(glm::vec3(1.0f, 1.0f, 0.4f));
    }
}

// ── History helpers ───────────────────────────────────────────────────────────

// Visualise the ancestry of the selected vertex as a point cloud.
// Shows original positions of all vertices that were ever merged into vid.
static void ShowAncestry(const SlabMesh& sm, unsigned vid)
{
    const auto& hist = sm.history;
    const std::vector<unsigned> ancestors = hist.GetAncestors(vid);

    // Build a direct map: vid → pre-merge position from the collapse log.
    // Each record stores the exact sphere centres of both sources before the merge.
    std::unordered_map<unsigned, std::array<double,3>> pos_map;
    for (const auto& rec : hist.Records()) {
        pos_map.emplace(rec.vid_src1, rec.pos_src1);
        pos_map.emplace(rec.vid_src2, rec.pos_src2);
    }

    std::vector<std::array<double,3>> pts;
    pts.reserve(ancestors.size());
    for (unsigned av : ancestors) {
        if (av < sm.vertices.size() && sm.vertices[av].first) {
            // Still alive — use current position.
            const auto& c = sm.vertices[av].second->sphere.center;
            pts.push_back({c.X(), c.Y(), c.Z()});
        } else {
            // Removed — use the exact position captured before it was merged.
            auto it = pos_map.find(av);
            if (it != pos_map.end())
                pts.push_back(it->second);
        }
    }

    if (pts.empty()) return;
    auto* pc = polyscope::registerPointCloud("Ancestry", pts);
    pc->setPointColor(glm::vec3(0.2f, 1.0f, 0.5f));  // green
    pc->setPointRadius(0.0018, true);
    pc->setEnabled(true);
}

// Visualise one step from the lineage of the selected vertex.
// Shows the bplists of vid_src1 and vid_src2 before the merge, and the
// combined bplist of vid_tgt after, all on the input mesh surface.
static void ShowLineageStep(const SlabMesh& sm, const CollapseRecord& rec)
{
    const double inv_diag = 1.0 / sm.pmesh->bb_diagonal_length;

    auto bpToPositions = [&](const std::vector<unsigned>& bps)
        -> std::vector<std::array<double,3>>
    {
        std::vector<std::array<double,3>> out;
        out.reserve(bps.size());
        for (unsigned bp : bps) {
            if (bp >= sm.pmesh->pVertexList.size()) continue;
            const auto& p = sm.pmesh->pVertexList[bp]->point();
            out.push_back({p[0]*inv_diag, p[1]*inv_diag, p[2]*inv_diag});
        }
        return out;
    };

    namespace ps = polyscope;

    // src1 bplist — blue
    auto pts1 = bpToPositions(rec.bplist_src1);
    if (!pts1.empty()) {
        auto* pc = ps::registerPointCloud("Lineage BPList src1", pts1);
        pc->setPointColor(glm::vec3(0.3f, 0.5f, 1.0f));
        pc->setPointRadius(0.0020, true);
        pc->setEnabled(true);
    }
    // src2 bplist — orange
    auto pts2 = bpToPositions(rec.bplist_src2);
    if (!pts2.empty()) {
        auto* pc = ps::registerPointCloud("Lineage BPList src2", pts2);
        pc->setPointColor(glm::vec3(1.0f, 0.55f, 0.15f));
        pc->setPointRadius(0.0020, true);
        pc->setEnabled(true);
    }
    // merged bplist — cyan
    auto ptsT = bpToPositions(rec.bplist_after);
    if (!ptsT.empty()) {
        auto* pc = ps::registerPointCloud("Lineage BPList merged", ptsT);
        pc->setPointColor(glm::vec3(0.0f, 1.0f, 0.85f));
        pc->setPointRadius(0.0022, true);
        pc->setEnabled(false);
    }
    // disable BPList selected 
    if (ps::hasPointCloud("BPList selected"))
        ps::getPointCloud("BPList selected")->setEnabled(false);
}

static void HideLineageStructures()
{
    namespace ps = polyscope;
    for (const char* name : {"Lineage BPList src1", "Lineage BPList src2", "Lineage BPList merged"})
        if (ps::hasPointCloud(name)) ps::getPointCloud(name)->setEnabled(false);
}

// Register a snapshot as the displayed MAT (separate from the live MAT).
static void ShowHistorySnapshot(const MeshSnapshot& snap)
{
    std::vector<std::array<double,3>> verts;
    verts.reserve(snap.verts.size());
    for (const auto& v : snap.verts)
        verts.push_back(v.pos);

    if (verts.empty()) return;

    bool en = polyscope::hasPointCloud("History Snapshot")
              ? polyscope::getPointCloud("History Snapshot")->isEnabled() : true;
    auto* pc = polyscope::registerPointCloud("History Snapshot", verts);
    pc->setPointColor(glm::vec3(0.85f, 0.85f, 0.85f));
    pc->setPointRadius(0.0012, true);
    pc->setEnabled(en);
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
    if (!init.faces.empty()) {
        bool en = ps::hasSurfaceMesh("Initial MAT Faces")
                  ? ps::getSurfaceMesh("Initial MAT Faces")->isEnabled() : false;
        auto* mm = ps::registerSurfaceMesh("Initial MAT Faces", init.verts, init.faces);
        mm->setSurfaceColor(glm::vec3(0.55f, 0.55f, 0.55f));
        mm->setTransparency(0.70f);
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

    // Live MAT — updated every N collapses
    UpdateMatStructures(init, vs);

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
    // Placeholder bplist clouds (per-collapse: v1 bps blue, v2 bps orange)
    {
        std::vector<std::array<double,3>> p = {{0,0,0}};

        auto* bp1 = ps::registerPointCloud("BPList v1", p);
        bp1->setPointColor(glm::vec3(0.30f, 0.50f, 1.0f));
        bp1->setPointRadius(0.0020, true);
        bp1->setEnabled(false);

        auto* bp2 = ps::registerPointCloud("BPList v2", p);
        bp2->setPointColor(glm::vec3(1.0f, 0.55f, 0.10f));
        bp2->setPointRadius(0.0020, true);
        bp2->setEnabled(false);

        // For manual pick selection
        auto* bpSel = ps::registerPointCloud("BPList selected", p);
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
        ImGui::SliderInt("Update every N", &vs.update_every, 1, 500);
        ImGui::Separator();

        // ── pick handling: click a MAT vertex to see its nmn_bplist ──────────
        if (polyscope::pick::haveSelection()) {
            auto [struct_ptr, local_idx] = polyscope::pick::getSelection();
            // Only react when the user clicked the live MAT vertex cloud
            if (struct_ptr == polyscope::getPointCloud("MAT Verts") &&
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
                    vs.selected_vid    = (int)vid;
                    vs.lineage_cursor  = 0;
                    vs.show_ancestry   = false;
                    if (polyscope::hasPointCloud("Ancestry"))
                        polyscope::getPointCloud("Ancestry")->setEnabled(false);
                    HideLineageStructures();
                    ShowBplistClusters(sm, vid);
                }
            }
        }

        // Auto-clear selection if the selected vertex was collapsed away.
        if (vs.selected_vid >= 0 &&
            ((unsigned)vs.selected_vid >= sm.vertices.size() ||
             !sm.vertices[vs.selected_vid].first))
        {
            vs.selected_vid = -1;
            polyscope::pick::resetSelection();
            if (polyscope::hasPointCloud("BPList selected"))
                polyscope::getPointCloud("BPList selected")->setEnabled(false);
            if (polyscope::hasPointCloud("MAT Vert Selected"))
                polyscope::getPointCloud("MAT Vert Selected")->setEnabled(false);
            if (polyscope::hasPointCloud("MAT Verts"))
                polyscope::getPointCloud("MAT Verts")->setEnabled(true);
        }

        // Show info about selected vertex
        if (vs.selected_vid >= 0 &&
            (unsigned)vs.selected_vid < sm.vertices.size() &&
            sm.vertices[vs.selected_vid].first)
        {
            const auto& sv = *sm.vertices[vs.selected_vid].second;
            const auto ct_idx = static_cast<uint8_t>(sv.nmn_cluster_type);
            const auto tt_idx = static_cast<uint8_t>(sv.topo_type);
            ImGui::Text("Selected vertex: %d", vs.selected_vid);
            ImGui::Text("  T-type: %s", kClusterTypeNames[ct_idx < 7 ? ct_idx : 5]);
            ImGui::Text("  Topo type: %s  (nf=%u)", kTopoTypeNames[tt_idx < 5 ? tt_idx : 0], sv.nf);
            ImGui::Text("  nmn_bplist size: %d", (int)sv.nmn_bplist.size());
            ImGui::Text("  clusters: %d", (int)sv.nmn_bplist_clusters.size());
            ImGui::Text("  (bplist coloured by cluster)");
            if (ImGui::Button("Clear selection")) {
                vs.selected_vid = -1;
                polyscope::pick::resetSelection();   // prevent next-frame re-trigger
                if (polyscope::hasPointCloud("BPList selected"))
                    polyscope::getPointCloud("BPList selected")->setEnabled(false);
                if (polyscope::hasPointCloud("MAT Vert Selected"))
                    polyscope::getPointCloud("MAT Vert Selected")->setEnabled(false);
                if (polyscope::hasPointCloud("MAT Verts"))
                    polyscope::getPointCloud("MAT Verts")->setEnabled(true);
                if (polyscope::hasSurfaceMesh("Input Mesh"))
                    polyscope::getSurfaceMesh("Input Mesh")->setEnabled(false);
            }
        } else {
            ImGui::TextDisabled("Click a MAT vertex to see its bplist");
        }

        // ── Vertex color mode toggle ──────────────────────────────────────────
        ImGui::Separator();
        ImGui::Text("MAT Vertex Coloring:");
        ImGui::SameLine();
        bool use_cluster = vs.color_mode == ViewerState::ColorMode::ClusterType;
        bool use_topo    = vs.color_mode == ViewerState::ColorMode::TopoType;
        if (ImGui::RadioButton("Cluster Type##cm", use_cluster)) {
            vs.color_mode = ViewerState::ColorMode::ClusterType;
            if (polyscope::hasPointCloud("MAT Verts")) {
                auto* pc = polyscope::getPointCloud("MAT Verts");
                auto* q_ct = pc->getQuantity("Cluster Type");
                auto* q_tt = pc->getQuantity("Topo Type");
                if (q_ct) q_ct->setEnabled(true);
                if (q_tt) q_tt->setEnabled(false);
            }
        }
        ImGui::SameLine();
        if (ImGui::RadioButton("Topo Type##cm", use_topo)) {
            vs.color_mode = ViewerState::ColorMode::TopoType;
            if (polyscope::hasPointCloud("MAT Verts")) {
                auto* pc = polyscope::getPointCloud("MAT Verts");
                auto* q_ct = pc->getQuantity("Cluster Type");
                auto* q_tt = pc->getQuantity("Topo Type");
                if (q_ct) q_ct->setEnabled(false);
                if (q_tt) q_tt->setEnabled(true);
            }
        }

        // Legend for current mode
        if (use_topo) {
            for (int k = 0; k < 5; ++k) {
                const auto& col = kTopoTypeColors[k];
                ImGui::TextColored(ImVec4(col[0], col[1], col[2], 1.0f), "  %s", kTopoTypeNames[k]);
            }
        }

        ImGui::Separator();
        if (ImGui::Button("Export MAT as OFF")) {
            std::string path = vs.outputPrefix
                + "_snapshot_" + std::to_string(vs.collapse_count) + ".off";
            ExportMatAsOff(sm, path);
        }

        // ── History panel ─────────────────────────────────────────────────────
        ImGui::Separator();
        ImGui::Text("History  (%u collapses, %u keyframes)",
                    sm.history.TotalSteps(),
                    (unsigned)sm.history.Keyframes().size());

        // Keyframe scrubber
        {
            int total = (int)sm.history.TotalSteps();
            if (total > 0) {
                int cur = (vs.history_step < 0) ? total : vs.history_step;
                if (ImGui::SliderInt("Step##hist", &cur, 0, total)) {
                    vs.history_step = cur;
                    const MeshSnapshot* kf = sm.history.GetKeyframeBefore((unsigned)cur);
                    if (kf) ShowHistorySnapshot(*kf);
                }
                ImGui::SameLine();
                if (ImGui::Button("Live##hist")) {
                    vs.history_step = -1;
                    if (polyscope::hasPointCloud("History Snapshot"))
                        polyscope::getPointCloud("History Snapshot")->setEnabled(false);
                }
            } else {
                ImGui::TextDisabled("(no history yet)");
            }
        }

        // Per-vertex lineage (only shown when a vertex is selected)
        if (vs.selected_vid >= 0) {
            ImGui::Separator();
            ImGui::Text("Lineage for vertex %d", vs.selected_vid);

            // Ancestry toggle
            if (ImGui::Checkbox("Show ancestry cloud", &vs.show_ancestry)) {
                if (vs.show_ancestry)
                    ShowAncestry(sm, (unsigned)vs.selected_vid);
                else if (polyscope::hasPointCloud("Ancestry"))
                    polyscope::getPointCloud("Ancestry")->setEnabled(false);
            }

            // Lineage step browser
            std::vector<CollapseRecord> lineage =
                sm.history.GetLineage((unsigned)vs.selected_vid);
            int n_lin = (int)lineage.size();
            if (n_lin == 0) {
                ImGui::TextDisabled("  (no collapses yet for this vertex)");
            } else {
                ImGui::Text("  %d merge steps", n_lin);
                if (vs.lineage_cursor < 0) vs.lineage_cursor = 0;
                if (vs.lineage_cursor >= n_lin) vs.lineage_cursor = n_lin - 1;

                ImGui::SliderInt("Merge step##lin", &vs.lineage_cursor, 0, n_lin - 1);

                ImGui::SameLine();
                if (ImGui::Button("<##lin") && vs.lineage_cursor > 0) {
                    --vs.lineage_cursor;
                }
                ImGui::SameLine();
                if (ImGui::Button(">##lin") && vs.lineage_cursor < n_lin - 1) {
                    ++vs.lineage_cursor;
                }

                const CollapseRecord& rec = lineage[(unsigned)vs.lineage_cursor];
                ImGui::Text("  step %u: v%u + v%u → v%u",
                            rec.step, rec.vid_src1, rec.vid_src2, rec.vid_tgt);
                ImGui::Text("  bplist: %d + %d → %d",
                            (int)rec.bplist_src1.size(),
                            (int)rec.bplist_src2.size(),
                            (int)rec.bplist_after.size());

                if (ImGui::Button("Show bplists at this step"))
                    ShowLineageStep(sm, rec);
                ImGui::SameLine();
                if (ImGui::Button("Hide bplists"))
                    HideLineageStructures();
            }
        }

        ImGui::Separator();
        ImGui::TextDisabled("Left panel layers:");
        ImGui::TextDisabled("  MAT Faces/Edges – orange/yellow");
        ImGui::TextDisabled("  Collapsed Edge  – red");
        ImGui::TextDisabled("  v1/v2 bplist    – blue/orange dots");
        ImGui::TextDisabled("  BPList selected – cluster colours");
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
        if (vs.paused) {
            while (vs.paused && !vs.step_once && !ps::windowRequestsClose()) {
                std::this_thread::sleep_for(std::chrono::milliseconds(16));
                ps::frameTick();
            }
            vs.step_once = false;
        }

        // ── rate-limit: skip render if too soon and no structural update due ─
        static const auto kMinInterval = std::chrono::milliseconds(33); // ~30 fps
        bool do_update = (vs.collapse_count % vs.update_every == 0);
        bool time_ok   = (clock::now() - vs.last_frame) >= kMinInterval;

        if (!do_update && !time_ok)
            return;

        vs.last_frame = clock::now();

        if (do_update) {
            // Rebuild live MAT (v1 + v2 still active at this point — pre-merge)
            UpdateMatStructures(BuildMatArrays(sm), vs);

            // Highlighted edge: v1 → v2
            std::vector<std::array<double,3>>  ep = {
                {v1p.X(), v1p.Y(), v1p.Z()},
                {v2p.X(), v2p.Y(), v2p.Z()}
            };
            std::vector<std::array<size_t,2>> ec = {{0, 1}};
            auto* ce = ps::registerCurveNetwork("Collapsed Edge", ep, ec);
            ce->setColor(glm::vec3(1.0f, 0.15f, 0.15f));
            ce->setRadius(0.0030f, true);
            ce->setEnabled(true);

            // v1, v2, result sphere-centre markers
            std::vector<std::array<double,3>> p1 = {{v1p.X(), v1p.Y(), v1p.Z()}};
            std::vector<std::array<double,3>> p2 = {{v2p.X(), v2p.Y(), v2p.Z()}};
            std::vector<std::array<double,3>> pr = {{result.center.X(),
                                                     result.center.Y(),
                                                     result.center.Z()}};
            ps::getPointCloud("v1")->updatePointPositions(p1);
            ps::getPointCloud("v1")->setEnabled(true);
            ps::getPointCloud("v2")->updatePointPositions(p2);
            ps::getPointCloud("v2")->setEnabled(true);
            ps::getPointCloud("result")->updatePointPositions(pr);
            ps::getPointCloud("result")->setEnabled(true);

            // ── bplist of v1 and v2 on the input mesh surface ────────────────
            auto bp1_pts = BplistPositions(sm, raw_v1);
            auto bp2_pts = BplistPositions(sm, raw_v2);

            if (!bp1_pts.empty()) {
                auto* bp1 = ps::registerPointCloud("BPList v1", bp1_pts);
                bp1->setPointColor(glm::vec3(0.30f, 0.50f, 1.0f));
                bp1->setPointRadius(0.0020, true);
                bp1->setEnabled(true);
            }
            if (!bp2_pts.empty()) {
                auto* bp2 = ps::registerPointCloud("BPList v2", bp2_pts);
                bp2->setPointColor(glm::vec3(1.0f, 0.55f, 0.10f));
                bp2->setPointRadius(0.0020, true);
                bp2->setEnabled(true);
            }
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
    double k = 0.00001;
    bool visualize = false;
    bool showHelp = false;
    bool valid = true;
    std::string errorMessage;
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
              << "  --k <value>        K factor for slab initialization (default: 0.00001)\n"
              << "  --output <prefix>  Output file prefix (default: input filename)\n"
#ifdef QMAT_WITH_POLYSCOPE
              << "  --visualize        Open live Polyscope viewer during simplification\n"
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


    // Step 4: If simplification requested, load into slab mesh and simplify
    if (options.simplifyTarget > 0) {
        std::cout << std::endl << "Loading MA for simplification..." << std::endl;

        // Setup slab mesh
        shape.slab_mesh.pmesh = &shape.input;
        shape.slab_mesh.type = 1;
        shape.slab_mesh.k = options.k;
        shape.slab_mesh.bound_weight = 1.0;
        shape.slab_mesh.export_prefix = options.outputPrefix;

        // Initialize slab mesh settings (same as GUI initialize())
        shape.slab_mesh.preserve_boundary_method = 0;
        shape.slab_mesh.hyperbolic_weight_type = 3;
        shape.slab_mesh.compute_hausdorff = false;
        shape.slab_mesh.boundary_compute_scale = 0;
        shape.slab_mesh.prevent_inversion = false;

        // Load the MA file into the slab mesh (computed above or loaded from cache)
        shape.LoadInputNMM(maFile);
        shape.ComputeFeatureEdges();          // detect sharp/concave edges on input mesh
        shape.slab_mesh.ClusterNMNBplist();   // must run before DetermineTopology (T1 filtering)
        shape.slab_mesh.DetermineTopology();  // uses T-types to correctly classify topology

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

#ifdef QMAT_WITH_POLYSCOPE
            if (options.visualize) {
                // Clear the per-collapse callback so frameTick is no longer
                // driven by a running Simplify loop.
                shape.slab_mesh.on_collapse_cb = nullptr;
                // Register the final simplified MAT and hand control to the
                // Polyscope window for interactive inspection.
                UpdateMatStructures(BuildMatArrays(shape.slab_mesh), vs);
                if (polyscope::hasCurveNetwork("Collapsed Edge"))
                    polyscope::getCurveNetwork("Collapsed Edge")->setEnabled(false);
                // Keep the original interactive callback intact so the user can
                // still click vertices, export OFF snapshots, etc.
                std::cout << "Simplification done. Close the viewer window to exit.\n";
                polyscope::show();
            }
#endif
        }
    }

    std::cout << std::endl << "Done!" << std::endl;

    return 0;
}
