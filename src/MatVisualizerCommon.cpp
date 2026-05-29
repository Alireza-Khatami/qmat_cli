#include "MatVisualizerCommon.h"

#ifdef QMAT_WITH_POLYSCOPE

#include <cmath>
#include <cstdint>
#include <map>
#include <set>
#include <unordered_map>

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"

#include "SlabMesh.h"

void MatVisualizer::Show()
{
    polyscope::show();
}

std::array<float,3> HsvToRgb(float h, float s, float v)
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

std::array<float,3> StructIdColor(int struct_id)
{
    if (struct_id < 0) return {0.5f, 0.5f, 0.5f};
    const float golden = 0.618033988749895f;
    float hue = std::fmod(struct_id * golden, 1.0f);
    return HsvToRgb(hue, 0.85f, 0.95f);
}

namespace {

// Bounding-box center of the input mesh.
std::array<double,3> InputMeshCenter(const MPMesh* pmesh)
{
    return {
        (pmesh->m_min[0] + pmesh->m_max[0]) * 0.5,
        (pmesh->m_min[1] + pmesh->m_max[1]) * 0.5,
        (pmesh->m_min[2] + pmesh->m_max[2]) * 0.5
    };
}

// Point clouds managed by the click/clear path; excluded from snapshots.
constexpr const char* kClickManagedNames[] = {
    "unsimp_mat_crspnd_points",
    "MAT Vert Selected",
};

} // namespace

MatArrays BuildMatArrays(const SlabMesh& sm)
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

        using CT2 = SlabVertex::ClusterType;
        if (ct == CT2::T0 || ct == CT2::T5)
            out.unknown_ttype_verts.push_back({c.X(), c.Y(), c.Z()});

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
        {
            const SlabEdge* se = sm.edges[i].second;
            if (se->struct_ids.empty()) {
                out.edge_structure_collapsible_colors.push_back({0.55f, 0.55f, 0.55f});
            } else {
                const SlabVertex* va2 = sm.vertices[se->vertices_.first].second;
                const SlabVertex* vb2 = sm.vertices[se->vertices_.second].second;
                bool match = (se->struct_ids == va2->struct_ids) && (se->struct_ids == vb2->struct_ids);
                out.edge_structure_collapsible_colors.push_back(match
                    ? std::array<float,3>{0.1f, 0.9f, 0.1f}
                    : std::array<float,3>{0.9f, 0.1f, 0.1f});
            }
        }
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
            out.ns_faces.push_back({ns_vid(a), ns_vid(b), ns_vid(c)});
    }

    for (unsigned i = 0; i < sm.edges.size(); ++i) {
        if (!sm.edges[i].first) continue;
        const SlabVertex* va = sm.vertices[sm.edges[i].second->vertices_.first].second;
        const SlabVertex* vb = sm.vertices[sm.edges[i].second->vertices_.second].second;
        if (va->nmn_cluster_type == CT::T1_spike || vb->nmn_cluster_type == CT::T1_spike) continue;
        size_t a = vid_map.at(sm.edges[i].second->vertices_.first);
        size_t b = vid_map.at(sm.edges[i].second->vertices_.second);
        auto ia = ns_remap.find(a), ib = ns_remap.find(b);
        if (ia == ns_remap.end() || ib == ns_remap.end()) continue;
        out.ns_edges.push_back({ia->second, ib->second});
    }

    return out;
}

void UpdateMatStructures(const MatArrays& arr, ViewerState& vs)
{
    namespace ps = polyscope;

    vs.idx_to_vid  = arr.idx_to_vid;
    vs.idx_to_eid  = arr.idx_to_eid;
    vs.mat_edge_node_count = arr.verts.size();

    if (!arr.faces.empty()) {
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
        cn->setColor(glm::vec3(0.0f, 0.0f, 0.0f));
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
        mm->setSurfaceColor(glm::vec3(0.0f, 0.0f, 0.0f));
        mm->setTransparency(0.45f);
        mm->setEnabled(en);
    } else if (ps::hasSurfaceMesh("Spike Faces")) {
        ps::removeStructure("Spike Faces");
    }

    if (!arr.ns_faces.empty()) {
        bool en = ps::hasSurfaceMesh("MAT Faces (No Spikes)")
                  ? ps::getSurfaceMesh("MAT Faces (No Spikes)")->isEnabled() : false;
        auto* mm = ps::registerSurfaceMesh("MAT Faces (No Spikes)", arr.ns_verts, arr.ns_faces);
        mm->setSurfaceColor(glm::vec3(0.3f, 0.8f, 1.0f));
        mm->setTransparency(0.25f);
        mm->setEnabled(en);
    } else if (ps::hasSurfaceMesh("MAT Faces (No Spikes)")) {
        ps::removeStructure("MAT Faces (No Spikes)");
    }

    if (ps::hasCurveNetwork("MAT Edges (No Spikes)")) {
        ps::removeStructure("MAT Edges (No Spikes)");
    }

    if (!arr.verts.empty()) {
        bool en = ps::hasPointCloud("MAT Verts")
                  ? ps::getPointCloud("MAT Verts")->isEnabled() : true;
        auto* pc = ps::registerPointCloud("MAT Verts", arr.verts);
        pc->setPointRadius(0.00297, true);
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

    if (!arr.unknown_ttype_verts.empty()) {
        bool en = ps::hasPointCloud("MAT Verts (Unknown TType)")
                  ? ps::getPointCloud("MAT Verts (Unknown TType)")->isEnabled()
                  : (vs.color_mode == ViewerState::ColorMode::UnknownTType);
        auto* upc = ps::registerPointCloud("MAT Verts (Unknown TType)", arr.unknown_ttype_verts);
        upc->setPointColor(glm::vec3(1.0f, 0.0f, 1.0f));
        upc->setPointRadius(0.0015, true);
        upc->setEnabled(en && vs.color_mode == ViewerState::ColorMode::UnknownTType);
    } else if (ps::hasPointCloud("MAT Verts (Unknown TType)")) {
        ps::getPointCloud("MAT Verts (Unknown TType)")->setEnabled(false);
    }

    if (!arr.sharp_verts.empty()) {
        auto* spc = ps::registerPointCloud("Sharp Feature Verts", arr.sharp_verts);
        spc->setPointColor(glm::vec3(1.0f, 0.15f, 0.15f));
        spc->setPointRadius(0.004, true);
        spc->setEnabled(false);
    } else if (ps::hasPointCloud("Sharp Feature Verts")) {
        ps::getPointCloud("Sharp Feature Verts")->setEnabled(false);
    }

    // Per-cluster-type filter clouds.
    static const char* kCFCloudNames[7] = {
        "MAT Verts [MS_Sheet]", "MAT Verts [MS_Seam]",
        "MAT Verts [MS_Boundary]", "MAT Verts [MS_Junction]",
        "MAT Verts [MS_Sheet_Boundary]", "MAT Verts [MS_Seam_Boundary]",
        "MAT Verts [MS_Junction_Boundary]"
    };
    static const int kCFCtIdx[7] = { 8, 9, 10, 11, 12, 13, 14 };
    for (int fi = 0; fi < 7; ++fi) {
        const auto& fv = arr.cluster_filter_verts[fi];
        const char* name = kCFCloudNames[fi];
        bool want_visible = (vs.cluster_filter == kCFCtIdx[fi]);
        if (!fv.empty()) {
            auto* fpc = ps::registerPointCloud(name, fv);
            const auto& col = kClusterTypeColors[kCFCtIdx[fi]];
            fpc->setPointColor(glm::vec3(col[0], col[1], col[2]));
            fpc->setPointRadius(0.003, true);
            fpc->setEnabled(want_visible);
        } else if (ps::hasPointCloud(name)) {
            ps::getPointCloud(name)->setEnabled(false);
        }
    }

    if (kModifyGlobalEdgeThickness)
        ApplyGlobalEdgeThickness(vs.edge_thickness);
}

void ApplyGlobalEdgeThickness(float r)
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

void RegisterInputMesh(const SlabMesh& sm)
{
    namespace ps = polyscope;
    if (!sm.pmesh) return;

    const double inv_diag = 1.0 / sm.pmesh->bb_diagonal_length;
    const auto cm = InputMeshCenter(sm.pmesh);
    std::vector<std::array<double,3>> verts;
    verts.reserve(sm.pmesh->pVertexList.size());
    for (unsigned i = 0; i < sm.pmesh->pVertexList.size(); ++i) {
        const auto& p = sm.pmesh->pVertexList[i]->point();
        verts.push_back({(p[0]-cm[0])*inv_diag, (p[1]-cm[1])*inv_diag, (p[2]-cm[2])*inv_diag});
    }

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

    auto registerFeatureEdges = [&](
        const std::set<std::array<int,2>>& edge_set,
        const char* name,
        glm::vec3 color,
        bool default_enabled)
    {
        std::vector<std::array<double,3>> nodes;
        std::vector<std::array<size_t,2>> segs;
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

void SnapshotEnabledPolyscopeStructures(
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

void DisableAllPolyscopeStructures()
{
    for (const auto& [type, named] : polyscope::state::structures)
        for (const auto& [name, sptr] : named)
            if (sptr) sptr->setEnabled(false);
}

void RestoreEnabledPolyscopeStructures(
    const std::vector<std::pair<std::string,std::string>>& snapshot)
{
    for (const auto& [type, name] : snapshot)
        if (polyscope::hasStructure(type, name))
            polyscope::getStructure(type, name)->setEnabled(true);
}

#endif  // QMAT_WITH_POLYSCOPE
