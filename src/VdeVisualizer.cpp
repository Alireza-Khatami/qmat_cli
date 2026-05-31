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

// Registers Initial MAT Faces/Edges/Struct Edges/Struct Verts from the slab
// state at Setup time (the slab is unsimplified — vcg-direct never mutates it).
// All coords are multiplied by sm.pmesh->bb_diagonal_length so they sit in the
// same world space as the scaled snapshot vcg-direct renders.
// All four layers start disabled; the "Initial MAT" toggle flips them on.
void RegisterVdeInitialMatScene(const SlabMesh& sm, ViewerState& vs)
{
    namespace ps = polyscope;
    using E2 = std::array<size_t, 2>;
    using C3 = std::array<float,  3>;

    const double scale = sm.pmesh ? sm.pmesh->bb_diagonal_length : 1.0;

    // Initial MAT Faces + Edges (geometry + face Struct-ID quantity).
    MatArrays init = BuildMatArrays(sm);
    vs.init_idx_to_fid          = init.idx_to_fid;
    vs.init_mat_face_vert_count = init.verts.size();
    std::vector<std::array<double,3>> verts_scaled;
    verts_scaled.reserve(init.verts.size());
    for (const auto& v : init.verts)
        verts_scaled.push_back({ v[0]*scale, v[1]*scale, v[2]*scale });

    if (!init.faces.empty()) {
        auto* mm = ps::registerSurfaceMesh("Initial MAT Faces", verts_scaled, init.faces);
        mm->setSurfaceColor(glm::vec3(0.55f, 0.55f, 0.55f));
        mm->setTransparency(1.0f);
        mm->setEnabled(false);
        if (!init.face_struct_id_colors.empty())
            mm->addFaceColorQuantity("Struct ID", init.face_struct_id_colors)->setEnabled(false);
    }
    if (!init.edges.empty()) {
        auto* cn = ps::registerCurveNetwork("Initial MAT Edges", verts_scaled, init.edges);
        cn->setColor(glm::vec3(0.6f, 0.6f, 0.6f));
        cn->setRadius(0.0005f, true);
        cn->setEnabled(false);
    }

    // Initial MAT Struct Edges — struct-bearing edges, coloured by struct_id.
    {
        std::vector<std::array<double,3>> nodes;
        std::vector<E2>                   segs;
        std::vector<C3>                   edge_colors;
        std::unordered_map<unsigned,size_t> remap;
        auto addV = [&](unsigned vid) -> size_t {
            auto it = remap.find(vid);
            if (it != remap.end()) return it->second;
            size_t idx = nodes.size();
            remap[vid] = idx;
            const auto& c = sm.vertices[vid].second->sphere.center;
            nodes.push_back({c.X()*scale, c.Y()*scale, c.Z()*scale});
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

    // Initial MAT Struct Verts — struct-bearing vertices, coloured by struct_id.
    {
        std::vector<std::array<double,3>> pts;
        std::vector<C3>                   pt_colors;
        for (unsigned i = 0; i < (unsigned)sm.vertices.size(); ++i) {
            if (!sm.vertices[i].first) continue;
            if (sm.vertices[i].second->struct_ids.empty()) continue;
            const auto& c = sm.vertices[i].second->sphere.center;
            pts.push_back({c.X()*scale, c.Y()*scale, c.Z()*scale});
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

// Flip the Initial vs. Current MAT layers.  Mirrors QMAT's ApplyInitialStructToggle.
void ApplyVdeInitialStructToggle(bool show_initial, bool show_struct_colors)
{
    namespace ps = polyscope;

    if (ps::hasSurfaceMesh("Initial MAT Faces")) {
        auto* mm = ps::getSurfaceMesh("Initial MAT Faces");
        mm->setEnabled(show_initial);
        if (auto* q = mm->getQuantity("Struct ID"))
            q->setEnabled(show_initial);
    }
    if (ps::hasCurveNetwork("Initial MAT Edges"))
        ps::getCurveNetwork("Initial MAT Edges")->setEnabled(show_initial);
    if (ps::hasCurveNetwork("Initial MAT Struct Edges"))
        ps::getCurveNetwork("Initial MAT Struct Edges")->setEnabled(show_initial);
    if (ps::hasPointCloud("Initial MAT Struct Verts"))
        ps::getPointCloud("Initial MAT Struct Verts")->setEnabled(show_initial);

    if (ps::hasSurfaceMesh("MAT Faces")) {
        auto* mm = ps::getSurfaceMesh("MAT Faces");
        mm->setEnabled(!show_initial);
        if (auto* q = mm->getQuantity("Struct ID"))
            q->setEnabled(!show_initial && show_struct_colors);
    }
    if (ps::hasCurveNetwork("MAT Edges"))
        ps::getCurveNetwork("MAT Edges")->setEnabled(!show_initial);
    bool cur = !show_initial && show_struct_colors;
    if (ps::hasCurveNetwork("MAT Struct Edges"))
        ps::getCurveNetwork("MAT Struct Edges")->setEnabled(cur);
    if (ps::hasPointCloud("MAT Struct Verts"))
        ps::getPointCloud("MAT Struct Verts")->setEnabled(cur);
}

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
            vs.show_initial_struct = false;
            ApplyVdeInitialStructToggle(false, vs.show_struct_colors);
        }
        // Initial ↔ Current MAT toggle, only meaningful when struct colors are on.
        if (vs.show_struct_colors || vs.show_initial_struct) {
            ImGui::SameLine();
            const char* lbl = vs.show_initial_struct ? "Current MAT##structtog" : "Initial MAT##structtog";
            if (ImGui::Button(lbl)) {
                vs.show_initial_struct = !vs.show_initial_struct;
                ApplyVdeInitialStructToggle(vs.show_initial_struct, vs.show_struct_colors);
            }
        }

        ImGui::Separator();

        // Pick handling — "MAT Verts" click shows initial-MAT ancestors,
        // "MAT Edges" click selects an edge slot for the info panel below.
        // local_idx is compact (matches snap arrays); for curve networks the
        // pick space is [nodes, edges), edge slot = local_idx - vertices.size().
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
                    vs.selected_eid = -1;
                    ShowVdeAncestorPoints(snap, (unsigned)vid, vs.enabled_snapshot);
                }
            }
            else if (polyscope::hasCurveNetwork("MAT Edges") &&
                     struct_ptr == polyscope::getCurveNetwork("MAT Edges"))
            {
                if (local_idx >= snap.vertices.size()) {
                    size_t edge_slot = local_idx - snap.vertices.size();
                    if (edge_slot < snap.edges.size()) {
                        vs.selected_eid = (int)edge_slot;
                        vs.selected_vid = -1;
                    }
                }
                polyscope::pick::resetSelection();
            }
            // Initial MAT Faces pick → show the imported .ma face elem_id.
            // Surface-mesh pick space: [verts, faces, edges, halfedges).
            else if (polyscope::hasSurfaceMesh("Initial MAT Faces") &&
                     struct_ptr == polyscope::getSurfaceMesh("Initial MAT Faces"))
            {
                if (local_idx >= vs.init_mat_face_vert_count &&
                    local_idx <  vs.init_mat_face_vert_count + vs.init_idx_to_fid.size())
                {
                    size_t face_slot = local_idx - vs.init_mat_face_vert_count;
                    vs.selected_init_fid = (int)vs.init_idx_to_fid[face_slot];
                    vs.selected_vid = -1;
                    vs.selected_eid = -1;
                }
                polyscope::pick::resetSelection();
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
            // T-type (cluster type) — shared MS_*/T0..T4 naming from MatVisualizerCommon.
            if (vid < (int)snap.vertex_cluster_type.size())
                ImGui::Text("  T-type: %s", ClusterTypeName(snap.vertex_cluster_type[vid]));
            // Full struct_ids set (empty == not part of any named struct).
            if (vid < (int)snap.vertex_struct_ids.size()) {
                const auto& sids = snap.vertex_struct_ids[vid];
                if (sids.empty()) {
                    ImGui::Text("  struct_ids: (none)");
                } else {
                    std::string s;
                    for (int id : sids) { if (!s.empty()) s += ", "; s += std::to_string(id); }
                    ImGui::Text("  struct_ids: {%s}", s.c_str());
                }
            }
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

        // Edge selection info — topo_type name + full struct_ids.
        // Auto-clear if the selected edge slot no longer exists in the snapshot.
        if (vs.selected_eid >= 0 && (size_t)vs.selected_eid >= snap.edges.size())
            vs.selected_eid = -1;
        if (vs.selected_eid >= 0) {
            const int eid = vs.selected_eid;
            const auto& e = snap.edges[eid];
            ImGui::Text("Selected edge: %d", eid);
            ImGui::Text("  v0=%d  v1=%d", e[0], e[1]);
            if (eid < (int)snap.edge_topo_type.size()) {
                auto tt = static_cast<SlabEdge::TopoType>(snap.edge_topo_type[eid]);
                ImGui::Text("  topo_type: %s", SlabEdge::TopoTypeName(tt));
            }
            if (eid < (int)snap.edge_struct_ids.size()) {
                const auto& sids = snap.edge_struct_ids[eid];
                if (sids.empty()) {
                    ImGui::Text("  struct_ids: (none)");
                } else {
                    std::string s;
                    for (int id : sids) { if (!s.empty()) s += ", "; s += std::to_string(id); }
                    ImGui::Text("  struct_ids: {%s}", s.c_str());
                }
            }
            if (ImGui::Button("Clear edge selection"))
                vs.selected_eid = -1;
            ImGui::Separator();
        }

        // Initial-MAT face selection info — face elem_id + struct_id.
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
            ImGui::Separator();
        } else if (vs.selected_init_fid >= 0) {
            vs.selected_init_fid = -1;
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
    RegisterVdeInitialMatScene(sm, vs_);

    // Live MAT layers must sit in the SAME scaled space as the snapshot render
    // (vcg-direct multiplies all coords by bb_diagonal_length).  BuildMatArrays
    // returns raw slab coords, so scale them before handing to UpdateMatStructures.
    {
        MatArrays live = BuildMatArrays(sm);
        const double scale = sm.pmesh ? sm.pmesh->bb_diagonal_length : 1.0;
        for (auto& v : live.verts)    { v[0]*=scale; v[1]*=scale; v[2]*=scale; }
        for (auto& v : live.ns_verts) { v[0]*=scale; v[1]*=scale; v[2]*=scale; }
        for (auto& v : live.unknown_ttype_verts) { v[0]*=scale; v[1]*=scale; v[2]*=scale; }
        for (auto& bucket : live.cluster_filter_verts)
            for (auto& v : bucket) { v[0]*=scale; v[1]*=scale; v[2]*=scale; }
        for (auto& v : live.sharp_verts) { v[0]*=scale; v[1]*=scale; v[2]*=scale; }
        UpdateMatStructures(live, vs_);
    }

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
        // Snapshot render re-enables live overlays; re-apply Initial-MAT mode.
        if (vs_.show_initial_struct)
            ApplyVdeInitialStructToggle(true, vs_.show_struct_colors);
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
    if (vs_.show_initial_struct)
        ApplyVdeInitialStructToggle(true, vs_.show_struct_colors);
}

#endif  // QMAT_WITH_POLYSCOPE && QMAT_WITH_VCGLIB
