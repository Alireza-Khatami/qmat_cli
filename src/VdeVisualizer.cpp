#include "VdeVisualizer.h"

#if defined(QMAT_WITH_POLYSCOPE) && defined(QMAT_WITH_VCGLIB)

#include <array>
#include <chrono>
#include <cstdint>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

#include "SlabMesh.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"
#include "polyscope/pick.h"
#include "imgui.h"

// ExportMatAsOff is defined non-static in main_cli.cpp.
void ExportMatAsOff(const SlabMesh& sm, const std::string& path);

namespace {

// In-place overwrite of "MAT Faces"/"MAT Edges"/"MAT Verts" plus the
// snapshot-derived struct/boundary overlays.
void RenderVcgDirectSnapshot(const VcgDirectSnapshot& snap,
                             bool show_struct_colors)
{
    namespace ps = polyscope;
    using F3 = std::array<size_t, 3>;
    using E2 = std::array<size_t, 2>;
    using C3 = std::array<float,  3>;

    // MAT Faces.
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

    // MAT Edges.
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

    // MAT Verts.
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

    // MAT Boundary Edges — both endpoints with vertex_topo_flags bit 3 set.
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

    // MAT Struct Edges — edges with non-empty struct_ids.
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

    // MAT Struct Verts — vertices with non-empty struct_ids.
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

    ps::frameTick();
}

// Show the initial-MAT ancestors of snapshot vertex `vid` as a point cloud.
// Mirrors QmatVisualizer::ShowUnsimpMatCrspndPoints but reads from the
// snapshot's vertex_original_ancestors + original_positions (compact indices,
// vcg-direct never exposes slab ids past BuildFromSlab).
void ShowVdeAncestorPoints(
    const VcgDirectSnapshot& snap, unsigned vid,
    std::vector<std::pair<std::string,std::string>>& enabled_snapshot)
{
    if (vid >= snap.vertex_original_ancestors.size()) return;
    const auto& ancestors = snap.vertex_original_ancestors[vid];

    std::vector<std::array<double,3>> pts;
    pts.reserve(ancestors.size());
    for (unsigned aid : ancestors) {
        if (aid < snap.original_positions.size())
            pts.push_back(snap.original_positions[aid]);
    }

    if (pts.empty()) {
        if (polyscope::hasPointCloud("unsimp_mat_crspnd_points"))
            polyscope::getPointCloud("unsimp_mat_crspnd_points")->setEnabled(false);
        return;
    }

    // Snapshot once per active selection so A->B click chains don't lose the
    // true pre-A scene state.  Mirrors the QMAT-side helper.
    if (enabled_snapshot.empty())
        SnapshotEnabledPolyscopeStructures(enabled_snapshot);
    DisableAllPolyscopeStructures();

    auto* uc = polyscope::registerPointCloud("unsimp_mat_crspnd_points", pts);
    uc->setPointRadius(0.0020, true);
    uc->setPointColor(glm::vec3(0.0f, 1.0f, 0.85f));
    uc->setEnabled(true);

    // Marker on the currently selected (live) MAT vertex.
    {
        const auto& v = snap.vertices[vid];
        std::vector<std::array<double,3>> mpt = {{ {v[0], v[1], v[2]} }};
        const uint8_t ct = (vid < snap.vertex_cluster_type.size())
                           ? snap.vertex_cluster_type[vid] : 5;
        const auto& col = kClusterTypeColors[ct < kClusterTypeColors.size() ? ct : 5];
        auto* mpc = polyscope::registerPointCloud("MAT Vert Selected", mpt);
        mpc->setPointColor(glm::vec3(col[0], col[1], col[2]));
        mpc->setPointRadius(0.0040, true);
        mpc->setEnabled(true);
    }
}

// VDE ImGui panel — pause/step, struct-color toggle, edge thickness, export,
// plus a "MAT Verts" pick handler that shows each vertex's initial-MAT
// ancestors (mirrors QmatVisualizer's ShowUnsimpMatCrspndPoints flow).
void InstallVdePanel(SlabMesh& sm, ViewerState& vs, VdeVisualizer& self)
{
    polyscope::state::userCallback = [&vs, &sm, &self]() {
        ImGui::PushItemWidth(230);
        ImGui::Text("QMAT vcg-direct Simplification Viewer");
        ImGui::Separator();
        ImGui::Text("Collapses: %d", vs.collapse_count);
        ImGui::Separator();
        ImGui::Checkbox("Pause", &vs.paused);
        if (vs.paused) {
            ImGui::SameLine();
            if (ImGui::Button("Step")) vs.step_once = true;
        }
        ImGui::Separator();

        if (kModifyGlobalEdgeThickness) {
            if (ImGui::SliderFloat("Edge Thickness", &vs.edge_thickness, 0.0001f, 0.005f, "%.4f"))
                ApplyGlobalEdgeThickness(vs.edge_thickness);
            ImGui::Separator();
        }

        // Overlays are already registered against the snapshot by
        // RenderVcgDirectSnapshot; the toggle just flips their visibility.
        if (ImGui::Checkbox("Show struct colors (seam/boundary/junction)", &vs.show_struct_colors)) {
            namespace ps = polyscope;
            if (ps::hasSurfaceMesh("MAT Faces"))
                if (auto* q = ps::getSurfaceMesh("MAT Faces")->getQuantity("Struct ID"))
                    q->setEnabled(vs.show_struct_colors);
            if (ps::hasCurveNetwork("MAT Struct Edges"))
                ps::getCurveNetwork("MAT Struct Edges")->setEnabled(vs.show_struct_colors);
            if (ps::hasPointCloud("MAT Struct Verts"))
                ps::getPointCloud("MAT Struct Verts")->setEnabled(vs.show_struct_colors);
        }

        ImGui::Separator();

        // Pick handling — "MAT Verts" click shows initial-MAT ancestors.
        // local_idx maps 1:1 to snap vid because RenderVcgDirectSnapshot
        // registers the cloud with snap.vertices in compact order.
        const VcgDirectSnapshot& snap = self.latest_snap;
        if (polyscope::pick::haveSelection()) {
            auto [struct_ptr, local_idx] = polyscope::pick::getSelection();
            if (polyscope::hasPointCloud("MAT Verts") &&
                struct_ptr == polyscope::getPointCloud("MAT Verts") &&
                local_idx < snap.vertices.size())
            {
                int vid = (int)local_idx;
                if (vid != vs.selected_vid) {
                    vs.selected_vid = vid;
                    ShowVdeAncestorPoints(snap, (unsigned)vid, vs.enabled_snapshot);
                }
            }
        }

        // Auto-clear selection if the picked vertex no longer exists in snap.
        if (vs.selected_vid >= 0 &&
            (size_t)vs.selected_vid >= snap.vertices.size())
        {
            vs.selected_vid = -1;
            polyscope::pick::resetSelection();
            if (polyscope::hasPointCloud("unsimp_mat_crspnd_points"))
                polyscope::getPointCloud("unsimp_mat_crspnd_points")->setEnabled(false);
            if (polyscope::hasPointCloud("MAT Vert Selected"))
                polyscope::getPointCloud("MAT Vert Selected")->setEnabled(false);
            RestoreEnabledPolyscopeStructures(vs.enabled_snapshot);
            vs.enabled_snapshot.clear();
        }

        // Selection info + clear button.
        if (vs.selected_vid >= 0 &&
            (size_t)vs.selected_vid < snap.vertices.size())
        {
            const int vid = vs.selected_vid;
            ImGui::Text("Selected vertex: %d", vid);
            if (vid < (int)snap.vertex_original_ancestors.size()) {
                const auto& anc = snap.vertex_original_ancestors[vid];
                ImGui::Text("  original_ancestors: %d", (int)anc.size());
                if (!anc.empty()) {
                    constexpr size_t kMaxShown = 32;
                    std::string s;
                    size_t shown = 0;
                    for (unsigned id : anc) {
                        if (shown >= kMaxShown) { s += ", ..."; break; }
                        if (!s.empty()) s += ", ";
                        s += std::to_string(id);
                        ++shown;
                    }
                    ImGui::TextWrapped("    {%s}", s.c_str());
                }
            }
            if (ImGui::Button("Clear selection")) {
                vs.selected_vid = -1;
                polyscope::pick::resetSelection();
                if (polyscope::hasPointCloud("unsimp_mat_crspnd_points"))
                    polyscope::getPointCloud("unsimp_mat_crspnd_points")->setEnabled(false);
                if (polyscope::hasPointCloud("MAT Vert Selected"))
                    polyscope::getPointCloud("MAT Vert Selected")->setEnabled(false);
                RestoreEnabledPolyscopeStructures(vs.enabled_snapshot);
                vs.enabled_snapshot.clear();
            }
            ImGui::Separator();
        }

        if (ImGui::Button("Export MAT as OFF")) {
            std::string path = vs.outputPrefix
                + "_snapshot_" + std::to_string(vs.collapse_count) + ".off";
            ExportMatAsOff(sm, path);
        }

        ImGui::PopItemWidth();
    };
}

} // namespace

void VdeVisualizer::Setup(SlabMesh& sm)
{
    namespace ps = polyscope;

    ps::init();
    ps::options::programName = "QMAT vcg-direct Simplification Viewer";
    ps::view::bgColor = {0.10f, 0.10f, 0.14f, 1.0f};
    ps::options::groundPlaneMode = ps::GroundPlaneMode::None;

    RegisterInputMesh(sm);
    UpdateMatStructures(BuildMatArrays(sm), vs_);

    InstallVdePanel(sm, vs_, *this);
    sm.on_collapse_cb = nullptr;

    for (int i = 0; i < 5; ++i)
        ps::frameTick();
}

LiveUpdateCallback VdeVisualizer::MakeLiveCallback()
{
    return [this](const VcgDirectSnapshot& snap) {
        latest_snap = snap;   // cache for the panel's pick handler
        RenderVcgDirectSnapshot(snap, vs_.show_struct_colors);
        vs_.collapse_count++;
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
    latest_snap = snap;
    RenderVcgDirectSnapshot(snap, vs_.show_struct_colors);
}

#endif  // QMAT_WITH_POLYSCOPE && QMAT_WITH_VCGLIB
