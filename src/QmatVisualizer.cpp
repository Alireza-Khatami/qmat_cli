#include "QmatVisualizer.h"

#include <array>
#include <climits>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <vector>

#include "SlabMesh.h"

// ── QMAT-only CLI-side exporter (no polyscope dep) ───────────────────────
// Post-simplification visualize_info JSON.  Must be called BEFORE
// SlabMesh::Export() (that calls AdjustStorage(), which compacts storage
// and invalidates the edge_last_rejection map).

void ExportSimpVisualizeInfo(const SlabMesh& sm, const std::string& path)
{
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

#ifdef QMAT_WITH_POLYSCOPE

#include <chrono>
#include <cstdio>
#include <map>
#include <string>
#include <thread>
#include <unordered_map>

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"
#include "polyscope/pick.h"
#include "imgui.h"

#include "QemRejectionViz.h"

namespace {

// kClusterTypeNames / ClusterTypeName live in MatVisualizerCommon.h (shared with VDE).

// ── Current MAT struct overlay ───────────────────────────────────────────

void UpdateStructColorVisualization(const SlabMesh& sm, bool enabled)
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

    // Seam + boundary edges coloured by struct_id.
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
            auto* cn = ps::registerCurveNetwork("MAT Struct Edges", nodes, segs);
            cn->addEdgeColorQuantity("Struct ID", edge_colors)->setEnabled(true);
            cn->setRadius(0.003f, true);
            cn->setEnabled(true);
        }
    }

    // Junction vertices coloured by struct_id.
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

// Captures the initial-MAT struct colorization as persistent (initially
// disabled) Polyscope structures.  Run once after the initial MAT is
// registered.
void RegisterInitialStructViz(const SlabMesh& sm, const MatArrays& init_arrays)
{
    namespace ps = polyscope;

    if (ps::hasSurfaceMesh("Initial MAT Faces") && !init_arrays.face_struct_id_colors.empty()) {
        auto* mm = ps::getSurfaceMesh("Initial MAT Faces");
        mm->addFaceColorQuantity("Struct ID", init_arrays.face_struct_id_colors)->setEnabled(false);
    }

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

// Toggles between initial-MAT struct viz and current-MAT struct viz.
void ApplyInitialStructToggle(bool show_initial, bool show_struct_colors)
{
    namespace ps = polyscope;

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

// ── Rejection edge / cause viz ───────────────────────────────────────────

// uint8 RGB → float [0,1].
std::array<float,3> RejectionReasonColor(SlabMesh::RejectionReason rr)
{
    auto c = SlabMesh::RejectionReasonColorU8(rr);
    return { c[0] / 255.0f, c[1] / 255.0f, c[2] / 255.0f };
}

void UpdateRejectionEdgeColors(const SlabMesh& sm, ViewerState& vs)
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
        vs.rejection_eid_order.push_back(i);

        auto it = sm.edge_last_rejection.find(i);
        edge_colors.push_back(it != sm.edge_last_rejection.end()
            ? RejectionReasonColor(it->second)
            : std::array<float,3>{1.0f, 1.0f, 1.0f});
    }

    if (segs.empty()) return;

    vs.rejection_node_count = nodes.size();

    bool en = ps::hasCurveNetwork("MAT Rejection Edges")
              ? ps::getCurveNetwork("MAT Rejection Edges")->isEnabled() : false;
    auto* cn = ps::registerCurveNetwork("MAT Rejection Edges", nodes, segs);
    cn->setRadius(0.0010f, true);
    cn->setEnabled(en);
    cn->addEdgeColorQuantity("Rejection Reason", edge_colors)->setEnabled(true);
}

void ClearRejectionPrimitives()
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

void ShowRejectionPrimitives(const SlabMesh& sm, unsigned eid, bool show_spheres = false)
{
    namespace ps = polyscope;

    auto it = sm.edge_reason_primitives.find(eid);
    if (it == sm.edge_reason_primitives.end()) { ClearRejectionPrimitives(); return; }
    const auto& prims = it->second;

    auto vertPos = [&](unsigned vid) -> std::array<double,3> {
        if (vid < sm.vertices.size() && sm.vertices[vid].first) {
            const auto& c = sm.vertices[vid].second->sphere.center;
            return {c.X(), c.Y(), c.Z()};
        }
        return {0.0, 0.0, 0.0};
    };

    if (!prims.vertices.empty()) {
        std::vector<std::array<double,3>> pts;
        pts.reserve(prims.vertices.size());
        for (unsigned vid : prims.vertices)
            pts.push_back(vertPos(vid));
        auto* pc = ps::registerPointCloud("Rejection Verts", pts);
        pc->setPointColor(glm::vec3(1.0f, 1.0f, 0.0f));
        pc->setPointRadius(0.006, true);
        pc->setEnabled(true);
    } else if (ps::hasPointCloud("Rejection Verts")) {
        ps::getPointCloud("Rejection Verts")->setEnabled(false);
    }

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
        cn->setColor(glm::vec3(1.0f, 1.0f, 0.0f));
        cn->setRadius(0.004f, true);
        cn->setEnabled(true);
    } else if (ps::hasCurveNetwork("Rejection Edges")) {
        ps::getCurveNetwork("Rejection Edges")->setEnabled(false);
    }

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
        mm->setSurfaceColor(glm::vec3(1.0f, 1.0f, 0.0f));
        mm->setTransparency(0.35f);
        mm->setEnabled(true);
    } else if (ps::hasSurfaceMesh("Rejection Faces")) {
        ps::getSurfaceMesh("Rejection Faces")->setEnabled(false);
    }

    if (prims.targ_ver.has_value()) {
        std::vector<std::array<double,3>> pts = { *prims.targ_ver };
        auto* pc = ps::registerPointCloud("Rejection Target", pts);
        pc->setPointColor(glm::vec3(0.0f, 1.0f, 0.0f));
        pc->setPointRadius(0.008, true);
        pc->setEnabled(true);
    } else if (ps::hasPointCloud("Rejection Target")) {
        ps::getPointCloud("Rejection Target")->setEnabled(false);
    }

    // Flipped face (InversionWouldOccur).
    if (prims.flipped_face.has_value()) {
        const auto& f = *prims.flipped_face;
        {
            std::vector<std::array<double,3>> nodes = { f[0][0], f[0][1], f[0][2] };
            std::vector<std::array<size_t,3>> tris  = { {0, 1, 2} };
            auto* mm = ps::registerSurfaceMesh("Flipped Face Before", nodes, tris);
            mm->setSurfaceColor(glm::vec3(1.0f, 0.0f, 0.0f));
            mm->setTransparency(1.0f);
            mm->setEnabled(true);
        }
        {
            std::vector<std::array<double,3>> nodes = { f[1][0], f[1][1], f[1][2] };
            std::vector<std::array<size_t,3>> tris  = { {0, 1, 2} };
            auto* mm = ps::registerSurfaceMesh("Flipped Face After", nodes, tris);
            mm->setSurfaceColor(glm::vec3(1.0f, 0.55f, 0.0f));
            mm->setTransparency(1.0f);
            mm->setEnabled(true);
        }
    } else {
        if (ps::hasSurfaceMesh("Flipped Face Before"))
            ps::getSurfaceMesh("Flipped Face Before")->setEnabled(false);
        if (ps::hasSurfaceMesh("Flipped Face After"))
            ps::getSurfaceMesh("Flipped Face After")->setEnabled(false);
    }

    // Endpoint + target spheres (true medial radii).
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
                if (radius_ptr < 1e-6f) radius_ptr = 1e-4f;
                radii.push_back(radius_ptr);
            }
        };
        addSphere(v1);
        addSphere(v2);

        if (prims.targ_ver.has_value()) {
            pts.push_back(*prims.targ_ver);
            radii.push_back(sm.edges[eid].second->sphere.radius);
        }

        if (!pts.empty()) {
            auto* pc = ps::registerPointCloud("Rejection Spheres", pts);
            pc->setPointColor(glm::vec3(0.4f, 0.8f, 1.0f));
            pc->setTransparency(0.5f);
            pc->addScalarQuantity("radius", radii)->setEnabled(false);
            pc->setPointRadiusQuantity("radius", false);
            pc->setEnabled(true);
        }
    } else if (ps::hasPointCloud("Rejection Spheres")) {
        ps::getPointCloud("Rejection Spheres")->setEnabled(false);
    }
}

// ── Pre-simplification "ancestors of this vertex" click viz ──────────────

void ShowUnsimpMatCrspndPoints(
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

    // Snapshot pre-click scene state once per active selection so A→B clicks
    // don't overwrite the true pre-A snapshot with the post-A scene state.
    if (enabled_snapshot.empty())
        SnapshotEnabledPolyscopeStructures(enabled_snapshot);
    DisableAllPolyscopeStructures();

    auto* uc = polyscope::registerPointCloud("unsimp_mat_crspnd_points", pts);
    uc->setPointRadius(0.0020, true);
    uc->setPointColor(glm::vec3(0.0f, 1.0f, 0.85f));
    uc->setEnabled(true);

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

    if (polyscope::hasSurfaceMesh("Initial MAT Faces")) {
        auto ps_mesh = polyscope::getSurfaceMesh("Initial MAT Faces");
        ps_mesh->setEnabled(true);
        ps_mesh->setEdgeWidth(1.24f);
    }
}

// ── InstallQmatScene ─────────────────────────────────────────────────────
//
// QMAT-flavored initial polyscope scene: input mesh, full initial MAT (with
// struct viz), live MAT, QEM rejection seed, per-collapse placeholders, and
// the full QMAT ImGui panel.  Called once from QmatVisualizer::Setup.

void InstallQmatScene(SlabMesh& sm, ViewerState& vs)
{
    namespace ps = polyscope;

    ps::init();
    ps::options::programName = "QMAT Simplification Viewer";
    ps::view::bgColor = {0.10f, 0.10f, 0.14f, 1.0f};
    ps::options::groundPlaneMode = ps::GroundPlaneMode::None;

    RegisterInputMesh(sm);

    // Initial MAT.  init_idx_to_fid is captured once because only the
    // initial-MAT face ids correspond to imported .ma indices; live MAT face
    // ids drift across collapses.
    MatArrays init = BuildMatArrays(sm);
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

    RegisterInitialStructViz(sm, init);
    UpdateMatStructures(init, vs);

#if defined(ONLY_USE_QEM_CONDITION_CHECKS)
    qemviz::UpdateEdgeColors(sm, vs.qem_viz);
#endif

    // Placeholder Collapsed Edge highlight.
    {
        std::vector<std::array<double,3>>  p = {{0,0,0},{0,0,0}};
        std::vector<std::array<size_t,2>>  e = {{0,1}};
        auto* ce = ps::registerCurveNetwork("Collapsed Edge", p, e);
        ce->setColor(glm::vec3(1.0f, 0.15f, 0.15f));
        ce->setRadius(0.0030f, true);
        ce->setEnabled(false);
    }
    // Placeholder v1 / v2 / result points.
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
    // Click-selection placeholders.
    {
        std::vector<std::array<double,3>> p = {{0,0,0}};

        auto* bpSel = ps::registerPointCloud("unsimp_mat_crspnd_points", p);
        bpSel->setPointColor(glm::vec3(0.0f, 1.0f, 0.85f));
        bpSel->setPointRadius(0.0020, true);
        bpSel->setEnabled(false);

        auto* mSel = ps::registerPointCloud("MAT Vert Selected", p);
        mSel->setPointRadius(0.0040, true);
        mSel->setEnabled(false);
    }

    // QMAT ImGui panel.
    polyscope::state::userCallback = [&vs, &sm]() {
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

        // Pick handling.
        if (polyscope::pick::haveSelection()) {
            auto [struct_ptr, local_idx] = polyscope::pick::getSelection();

#if defined(ONLY_USE_QEM_CONDITION_CHECKS)
            if (qemviz::HandlePick(sm, vs.qem_viz, struct_ptr, local_idx)) {
                polyscope::pick::resetSelection();
            } else
#endif
            // MAT Rejection Edges pick → highlight cause.
            if (polyscope::hasCurveNetwork("MAT Rejection Edges") &&
                struct_ptr == polyscope::getCurveNetwork("MAT Rejection Edges"))
            {
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
            // MAT Verts pick → show bplist ancestors.
            else if (struct_ptr == polyscope::getPointCloud("MAT Verts") &&
                local_idx < vs.idx_to_vid.size())
            {
                unsigned vid = vs.idx_to_vid[local_idx];

                using clock = std::chrono::steady_clock;
                auto now = clock::now();
                auto ms  = std::chrono::duration_cast<std::chrono::milliseconds>(
                               now - vs.last_pick_time).count();
                bool is_double_click = ((int)vid == vs.last_picked_vid) &&
                                       (ms < ViewerState::kDoubleClickMs);
                vs.last_picked_vid  = (int)vid;
                vs.last_pick_time   = now;

                if (is_double_click) {
                    // Camera pivots to vertex but keeps its position.
                    const auto& c = sm.vertices[vid].second->sphere.center;
                    glm::vec3 target = {(float)c.X(), (float)c.Y(), (float)c.Z()};
                    glm::vec3 camPos = polyscope::view::getCameraWorldPosition();
                    polyscope::view::lookAt(camPos, target, /*flyTo=*/true);
                    polyscope::pick::resetSelection();
                } else if ((int)vid != vs.selected_vid) {
                    vs.selected_vid           = (int)vid;
                    vs.selected_eid           = -1;
                    vs.selected_init_fid      = -1;
                    vs.selected_rejection_eid = -1;
                    ClearRejectionPrimitives();
                    ShowUnsimpMatCrspndPoints(sm, vid, vs.enabled_snapshot);
                }
            }
            // MAT edge pick → show struct_ids.
            else if (polyscope::hasCurveNetwork("MAT Edges") &&
                     struct_ptr == polyscope::getCurveNetwork("MAT Edges"))
            {
                if (local_idx >= vs.mat_edge_node_count) {
                    size_t edge_slot = local_idx - vs.mat_edge_node_count;
                    if (edge_slot < vs.idx_to_eid.size()) {
                        vs.selected_eid = (int)vs.idx_to_eid[edge_slot];
                        vs.selected_vid = -1;
                        vs.selected_init_fid = -1;
                    }
                }
                polyscope::pick::resetSelection();
            }
            // Initial MAT Faces pick → show imported .ma face idx.
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

        // Auto-clear vertex selection if the vertex was collapsed away.
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
            // Initial MAT Faces forced off on clear regardless of pre-click state.
            if (polyscope::hasSurfaceMesh("Initial MAT Faces"))
                polyscope::getSurfaceMesh("Initial MAT Faces")->setEnabled(false);
        }

        // Vertex selection info.
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
                polyscope::pick::resetSelection();
                if (polyscope::hasPointCloud("unsimp_mat_crspnd_points"))
                    polyscope::getPointCloud("unsimp_mat_crspnd_points")->setEnabled(false);
                if (polyscope::hasPointCloud("MAT Vert Selected"))
                    polyscope::getPointCloud("MAT Vert Selected")->setEnabled(false);
                RestoreEnabledPolyscopeStructures(vs.enabled_snapshot);
                vs.enabled_snapshot.clear();
                if (polyscope::hasSurfaceMesh("Initial MAT Faces"))
                    polyscope::getSurfaceMesh("Initial MAT Faces")->setEnabled(false);
            }
        } else {
            ImGui::TextDisabled("Click a MAT vertex to see its bplist");
        }

        // Edge selection info.
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
            if (vs.selected_eid >= 0) vs.selected_eid = -1;
        }

        // Initial-MAT face selection info.
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
            if (vs.selected_init_fid >= 0) vs.selected_init_fid = -1;
        }

        // Rejection edge selection info.
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
        qemviz::DrawPanel(sm, vs.qem_viz);
#endif

        if (kModifyGlobalEdgeThickness) {
            ImGui::Separator();
            if (ImGui::SliderFloat("Edge Thickness", &vs.edge_thickness, 0.0001f, 0.005f, "%.4f"))
                ApplyGlobalEdgeThickness(vs.edge_thickness);
        }

        // Struct color overlay.
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

        // Vertex color mode.
        ImGui::Separator();
        ImGui::Text("MAT Vertex Coloring:");
        using CM = ViewerState::ColorMode;
        bool use_cluster      = vs.color_mode == CM::ClusterType;
        bool use_unk_ttype    = vs.color_mode == CM::UnknownTType;

        auto setColorMode = [&](CM new_mode) {
            vs.color_mode = new_mode;
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

        if (use_cluster) {
            for (int k = 0; k < 8; ++k) {
                const auto& col = kClusterTypeColors[7+k];
                ImGui::TextColored(ImVec4(col[0], col[1], col[2], 1.0f), "  %s", kClusterTypeNames[k]);
            }
        }

        // Cluster filter.
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
            if (polyscope::hasPointCloud("MAT Verts"))
                polyscope::getPointCloud("MAT Verts")->setEnabled(
                    showing_all && (vs.color_mode == CM::ClusterType));
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

    // Warm-up frames so the GL window is fully up before frameTick() fires
    // from the collapse loop.
    for (int i = 0; i < 5; ++i)
        polyscope::frameTick();
}

} // namespace

// ── QmatVisualizer::Setup ─────────────────────────────────────────────────

void QmatVisualizer::Setup(SlabMesh& sm)
{
    InstallQmatScene(sm, vs_);

    // Per-collapse callback — fired just before each MergeVertices call.
    sm.on_collapse_cb = [vsp = &vs_, smp = &sm](
        unsigned /*raw_v1*/, const Wm4::Vector3d& v1p, double /*v1r*/,
        unsigned /*raw_v2*/, const Wm4::Vector3d& v2p, double /*v2r*/,
        const Sphere& /*result*/)
    {
        namespace ps = polyscope;
        using clock  = std::chrono::steady_clock;
        ViewerState& vs = *vsp;
        SlabMesh&    sm = *smp;

        ++vs.collapse_count;

        if (ps::windowRequestsClose()) {
            sm.on_collapse_cb = nullptr;
            return;
        }

        // prev_was_manual_step forces a rebuild on the next callback so the
        // user-requested step's resulting state is visible (the callback fires
        // BEFORE the merge that produces it).
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
                prev_was_manual_step = true;
            }
        }

        // Auto-pause + per-collapse sleep once vs.collapse_delay_after is
        // reached (so the user can inspect the late-simplification frames).
        if (vs.collapse_delay_after >= 0 && vs.collapse_count == vs.collapse_delay_after) {
            vs.paused = true;
        }
        if (vs.collapse_delay_after >= 0 && vs.collapse_count >= vs.collapse_delay_after) {
            force_update = true;
            std::this_thread::sleep_for(std::chrono::milliseconds(vs.collapse_delay_ms));
            ps::frameTick();
        }

        // Rate-limit: skip rebuild if too soon and nothing structural is due.
        static const auto kMinInterval = std::chrono::milliseconds(33);  // ~30 fps
        bool do_update = (vs.collapse_count % vs.update_every == 0) || force_update;
        bool time_ok   = (clock::now() - vs.last_frame) >= kMinInterval;

        if (!do_update && !time_ok)
            return;

        vs.last_frame = clock::now();

        if (do_update) {
            // v1 and v2 are still active here (pre-merge).
            UpdateRejectionEdgeColors(sm, vs);
#if defined(ONLY_USE_QEM_CONDITION_CHECKS)
            qemviz::UpdateEdgeColors(sm, vs.qem_viz);
#endif
            UpdateMatStructures(BuildMatArrays(sm), vs);
            UpdateStructColorVisualization(sm, !vs.show_initial_struct && vs.show_struct_colors);

            // Highlighted edge: v1 → v2.
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
        }

        ps::frameTick();
    };
}

// ── QmatVisualizer::RenderFinal ──────────────────────────────────────────

void QmatVisualizer::RenderFinal(SlabMesh& sm)
{
    // Stop the live collapse loop from firing.
    sm.on_collapse_cb = nullptr;

    UpdateRejectionEdgeColors(sm, vs_);
#if defined(ONLY_USE_QEM_CONDITION_CHECKS)
    qemviz::UpdateEdgeColors(sm, vs_.qem_viz);
#endif
    UpdateMatStructures(BuildMatArrays(sm), vs_);
    UpdateStructColorVisualization(sm, !vs_.show_initial_struct && vs_.show_struct_colors);
    if (polyscope::hasCurveNetwork("Collapsed Edge"))
        polyscope::getCurveNetwork("Collapsed Edge")->setEnabled(false);
}

#endif  // QMAT_WITH_POLYSCOPE
