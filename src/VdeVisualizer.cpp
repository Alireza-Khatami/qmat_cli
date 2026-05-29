#include "VdeVisualizer.h"

#if defined(QMAT_WITH_POLYSCOPE) && defined(QMAT_WITH_VCGLIB)

#include <array>
#include <chrono>
#include <cstdint>
#include <thread>
#include <unordered_map>
#include <vector>

#include "SlabMesh.h"
#include "QmatVisualizer.h"   // for InstallSharedScene

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"

namespace {

// In-place overwrite of "MAT Faces"/"MAT Edges"/"MAT Verts" plus the
// snapshot-derived struct/boundary overlays.  Also drops slab-only structures
// that have no snapshot counterpart so they don't draw stale data.
void RenderVcgDirectSnapshot(const VcgDirectSnapshot& snap,
                             bool show_struct_colors)
{
    namespace ps = polyscope;
    using F3 = std::array<size_t, 3>;
    using E2 = std::array<size_t, 2>;
    using C3 = std::array<float,  3>;

    // MAT Faces ────────────────────────────────────────────────────────────
    bool mm_enabled = ps::hasSurfaceMesh("MAT Faces")
                      ? ps::getSurfaceMesh("MAT Faces")->isEnabled() : true;
    std::vector<F3> faces; faces.reserve(snap.faces.size());
    for (const auto& f : snap.faces)
        faces.push_back({ (size_t)f[0], (size_t)f[1], (size_t)f[2] });
    auto* mm = ps::registerSurfaceMesh("MAT Faces", snap.vertices, faces);
    mm->setSurfaceColor(glm::vec3(0.9f, 0.6f, 0.2f));
    mm->setTransparency(1.0f);
    if (!snap.face_struct_id.empty()) {
        std::vector<C3> face_colors; face_colors.reserve(snap.face_struct_id.size());
        for (int sid : snap.face_struct_id) face_colors.push_back(StructIdColor(sid));
        mm->addFaceColorQuantity("Struct ID", face_colors)->setEnabled(show_struct_colors);
    }
    mm->setEnabled(mm_enabled);

    // MAT Edges ────────────────────────────────────────────────────────────
    if (!snap.edges.empty()) {
        bool cn_enabled = ps::hasCurveNetwork("MAT Edges")
                          ? ps::getCurveNetwork("MAT Edges")->isEnabled() : true;
        std::vector<E2> edges; edges.reserve(snap.edges.size());
        for (const auto& e : snap.edges)
            edges.push_back({ (size_t)e[0], (size_t)e[1] });
        auto* cn = ps::registerCurveNetwork("MAT Edges", snap.vertices, edges);
        cn->setColor(glm::vec3(1.0f, 0.80f, 0.30f));
        cn->setRadius(0.0008f, true);

        if (!snap.edge_topo_type.empty()) {
            std::vector<C3> topo_colors; topo_colors.reserve(snap.edge_topo_type.size());
            for (uint8_t tt : snap.edge_topo_type)
                topo_colors.push_back(kEdgeTopoTypeColors[tt < kEdgeTopoTypeColors.size() ? tt : 0]);
            cn->addEdgeColorQuantity("Edge Topo Type", topo_colors)->setEnabled(true);
        }
        if (!snap.edge_struct_match.empty() && !snap.edge_first_struct_id.empty()) {
            std::vector<C3> col; col.reserve(snap.edge_struct_match.size());
            for (size_t i = 0; i < snap.edge_struct_match.size(); ++i) {
                const int  sid   = snap.edge_first_struct_id[i];
                const bool match = snap.edge_struct_match[i] != 0;
                if (sid < 0)        col.push_back({0.55f, 0.55f, 0.55f}); // grey: not a struct edge
                else if (match)     col.push_back({0.1f,  0.9f,  0.1f});  // green: collapsible
                else                col.push_back({0.9f,  0.1f,  0.1f});  // red: mismatch
            }
            cn->addEdgeColorQuantity("Structure Collapsible", col);
        }
        cn->setEnabled(cn_enabled);
    } else if (ps::hasCurveNetwork("MAT Edges")) {
        ps::removeStructure("MAT Edges");
    }

    // MAT Verts ────────────────────────────────────────────────────────────
    if (!snap.vertices.empty()) {
        bool pc_enabled = ps::hasPointCloud("MAT Verts")
                          ? ps::getPointCloud("MAT Verts")->isEnabled() : true;
        auto* pc = ps::registerPointCloud("MAT Verts", snap.vertices);
        pc->setPointRadius(0.00297, true);
        if (!snap.vertex_cluster_type.empty()) {
            std::vector<C3> vc; vc.reserve(snap.vertex_cluster_type.size());
            for (uint8_t ct : snap.vertex_cluster_type)
                vc.push_back(kClusterTypeColors[ct < kClusterTypeColors.size() ? ct : 5]);
            pc->addColorQuantity("Cluster Type", vc)->setEnabled(true);
        }
        pc->setEnabled(pc_enabled);
    }

    // MAT Boundary Edges — both endpoints with topo_flags bit 3 (boundary). ─
    {
        std::vector<std::array<double,3>> bnodes;
        std::vector<E2>                   bsegs;
        std::unordered_map<int,size_t>    remap;
        auto addV = [&](int vid) -> size_t {
            auto it = remap.find(vid);
            if (it != remap.end()) return it->second;
            size_t idx = bnodes.size();
            remap[vid] = idx;
            bnodes.push_back(snap.vertices[vid]);
            return idx;
        };
        for (const auto& e : snap.edges) {
            const int a = e[0], b = e[1];
            if (a < 0 || b < 0) continue;
            const bool ba = ((size_t)a < snap.vertex_topo_flags.size())
                            && (snap.vertex_topo_flags[a] & 0x8);
            const bool bb = ((size_t)b < snap.vertex_topo_flags.size())
                            && (snap.vertex_topo_flags[b] & 0x8);
            if (ba && bb) bsegs.push_back({addV(a), addV(b)});
        }
        if (!bsegs.empty()) {
            bool en = ps::hasCurveNetwork("MAT Boundary Edges")
                      ? ps::getCurveNetwork("MAT Boundary Edges")->isEnabled() : false;
            auto* cn = ps::registerCurveNetwork("MAT Boundary Edges", bnodes, bsegs);
            cn->setColor(glm::vec3(1.0f, 0.15f, 0.15f));
            cn->setRadius(0.0015f, true);
            cn->setEnabled(en);
        } else if (ps::hasCurveNetwork("MAT Boundary Edges")) {
            ps::removeStructure("MAT Boundary Edges");
        }
    }

    // MAT Struct Edges — edges with non-empty struct_ids. ──────────────────
    {
        std::vector<std::array<double,3>> nodes;
        std::vector<E2>                   segs;
        std::vector<C3>                   edge_colors;
        std::unordered_map<int,size_t>    remap;
        auto addV = [&](int vid) -> size_t {
            auto it = remap.find(vid);
            if (it != remap.end()) return it->second;
            size_t idx = nodes.size();
            remap[vid] = idx;
            nodes.push_back(snap.vertices[vid]);
            return idx;
        };
        for (size_t i = 0; i < snap.edges.size(); ++i) {
            if (i >= snap.edge_first_struct_id.size()) break;
            const int sid = snap.edge_first_struct_id[i];
            if (sid < 0) continue;
            const int a = snap.edges[i][0], b = snap.edges[i][1];
            if (a < 0 || b < 0) continue;
            segs.push_back({addV(a), addV(b)});
            edge_colors.push_back(StructIdColor(sid));
        }
        if (!segs.empty()) {
            auto* cn = ps::registerCurveNetwork("MAT Struct Edges", nodes, segs);
            cn->addEdgeColorQuantity("Struct ID", edge_colors)->setEnabled(true);
            cn->setRadius(0.003f, true);
            cn->setEnabled(show_struct_colors);
        } else if (ps::hasCurveNetwork("MAT Struct Edges")) {
            ps::removeStructure("MAT Struct Edges");
        }
    }

    // MAT Struct Verts — vertices with non-empty struct_ids. ───────────────
    {
        std::vector<std::array<double,3>> pts;
        std::vector<C3>                   pt_colors;
        for (size_t i = 0; i < snap.vertices.size(); ++i) {
            if (i >= snap.vertex_first_struct_id.size()) break;
            const int sid = snap.vertex_first_struct_id[i];
            if (sid < 0) continue;
            pts.push_back(snap.vertices[i]);
            pt_colors.push_back(StructIdColor(sid));
        }
        if (!pts.empty()) {
            auto* pc = ps::registerPointCloud("MAT Struct Verts", pts);
            pc->addColorQuantity("Struct ID", pt_colors)->setEnabled(true);
            pc->setPointRadius(0.005, true);
            pc->setEnabled(show_struct_colors);
        } else if (ps::hasPointCloud("MAT Struct Verts")) {
            ps::removeStructure("MAT Struct Verts");
        }
    }

    // Drop stale slab-only structures (no snapshot counterpart).
    static const char* kStaleSurfaceMeshes[] = {
        "Spike Faces", "MAT Faces (No Spikes)", "Initial MAT Faces",
    };
    static const char* kStalePointClouds[] = {
        "MAT Verts (Unknown TType)", "Sharp Feature Verts",
        "MAT Verts [MS_Sheet]",          "MAT Verts [MS_Seam]",
        "MAT Verts [MS_Boundary]",       "MAT Verts [MS_Junction]",
        "MAT Verts [MS_Sheet_Boundary]", "MAT Verts [MS_Seam_Boundary]",
        "MAT Verts [MS_Junction_Boundary]",
        "Initial MAT Struct Verts",
    };
    static const char* kStaleCurveNetworks[] = {
        "Spike Edges", "MAT Edges (No Spikes)",
        "Initial MAT Edges", "Initial MAT Struct Edges",
    };
    for (const char* n : kStaleSurfaceMeshes)
        if (ps::hasSurfaceMesh(n)) ps::removeStructure(n);
    for (const char* n : kStalePointClouds)
        if (ps::hasPointCloud(n))  ps::removeStructure(n);
    for (const char* n : kStaleCurveNetworks)
        if (ps::hasCurveNetwork(n)) ps::removeStructure(n);

    ps::frameTick();
}

} // namespace

void VdeVisualizer::Setup(SlabMesh& sm)
{
    InstallSharedScene(sm, vs_);
    sm.on_collapse_cb = nullptr;
    vs_.vcg_direct_active = true;
}

LiveUpdateCallback VdeVisualizer::MakeLiveCallback()
{
    // Captures vs_ by reference (visualizer must outlive the simplify call).
    return [this](const VcgDirectSnapshot& snap) {
        RenderVcgDirectSnapshot(snap, vs_.show_struct_colors);
        vs_.collapse_count++;
        // Pause/step spin — same pattern as QMAT's on_collapse_cb.
        if (vs_.paused) {
            while (vs_.paused && !vs_.step_once && !polyscope::windowRequestsClose()) {
                std::this_thread::sleep_for(std::chrono::milliseconds(16));
                polyscope::frameTick();
            }
            if (vs_.step_once) vs_.step_once = false;
        }
    };
}

void VdeVisualizer::Render(const VcgDirectSnapshot& snap)
{
    RenderVcgDirectSnapshot(snap, vs_.show_struct_colors);
}

#endif  // QMAT_WITH_POLYSCOPE && QMAT_WITH_VCGLIB
