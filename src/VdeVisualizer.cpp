#include "VdeVisualizer.h"

#ifdef QMAT_WITH_VCGLIB

#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <string>
#include <vector>

#include "SlabMesh.h"   // for RejectionReason{Name,ColorU8}, SlabMesh::RejectionReason::Count

// ── VDE-only CLI-side exporters (snapshot-driven, no polyscope dep) ──────
//
// All six functions reconstruct the corresponding QMAT-side file from a
// VcgDirectSnapshot alone, so the vcg-direct path can emit the same family
// of artefacts as the QMAT path without ever touching SlabMesh.  The
// snapshot already stores positions, faces, derived edges, per-vertex
// radius/cluster_type/struct_ids/topo_flags/ancestors, per-edge
// topo_type/struct_ids/last_rejection, and per-face struct_id.

namespace {

// Cluster-type names — indexed by SlabVertex::ClusterType uint8_t value.
// Mirrors the table in SlabMesh::ExportTypedMA so the .mat_typed text format
// stays identical.
constexpr std::array<const char*, 15> kSnapClusterNames = {{
    "T0", "T1_spike", "T2", "T3", "T4", "T5", "T1_non_spike",
    "MS_Unknown", "MS_Sheet", "MS_Seam", "MS_Boundary", "MS_Junction",
    "MS_Sheet_Boundary", "MS_Seam_Boundary", "MS_Junction_Boundary",
}};

constexpr std::array<std::array<int,3>, 15> kSnapClusterRgb = {{
    {230,   0, 230}, { 26,  26,  26}, {  0, 217, 255}, {255, 128,   0},
    {255,  26,  26}, {255, 255, 255}, {140, 140, 140}, { 89,  89,  89},
    {  0, 255,  77}, {255, 230,   0}, {  0, 128, 255}, {255,   0, 128},
    {  0, 230, 255}, {255,  89,   0}, {153,   0, 255},
}};

constexpr std::array<const char*, 6> kSnapEdgeTopoNames = {{
    "Unknown", "Sheet", "Seam", "Boundary", "Seam_Boundary", "Orphan",
}};
constexpr std::array<std::array<int,3>, 6> kSnapEdgeTopoRgb = {{
    {140, 140, 140}, {  0, 255,  77}, {255, 230,   0},
    {  0, 128, 255}, {255,  89,   0}, {230,   0, 230},
}};

inline const char* SnapClusterName(uint8_t ct) {
    return (ct < kSnapClusterNames.size()) ? kSnapClusterNames[ct] : "MS_Unknown";
}

inline std::array<int,3> SnapClusterRgb(uint8_t ct) {
    return (ct < kSnapClusterRgb.size()) ? kSnapClusterRgb[ct]
                                         : std::array<int,3>{89, 89, 89};
}

}  // namespace

void ExportSnapshotAsOff(const VcgDirectSnapshot& snap, const std::string& path)
{
    std::ofstream f(path);
    if (!f) {
        std::cerr << "[ExportSnapshotAsOff] cannot open: " << path << "\n";
        return;
    }
    f << "OFF\n";
    f << snap.vertices.size() << " " << snap.faces.size() << " 0\n";
    f << std::fixed << std::setprecision(10);
    for (const auto& v : snap.vertices)
        f << v[0] << " " << v[1] << " " << v[2] << "\n";
    for (const auto& tri : snap.faces)
        f << "3 " << tri[0] << " " << tri[1] << " " << tri[2] << "\n";
    std::cout << "[ExportSnapshotAsOff] wrote "
              << snap.vertices.size() << " verts, "
              << snap.faces.size() << " faces to " << path << "\n";
}

void ExportSnapshotVisualizeInfo(const VcgDirectSnapshot& snap, const std::string& path)
{
    std::ofstream f(path);
    if (!f) {
        std::cerr << "[ExportSnapshotVisualizeInfo] cannot open: " << path << "\n";
        return;
    }

    auto write_rgb_i = [&](const std::array<int,3>& c) {
        f << "[" << c[0] << "," << c[1] << "," << c[2] << "]";
    };
    auto write_rgb_u8 = [&](const std::array<uint8_t,3>& c) {
        f << "[" << (int)c[0] << "," << (int)c[1] << "," << (int)c[2] << "]";
    };
    auto write_intvec = [&](const std::vector<int>& s) {
        f << "[";
        for (size_t i = 0; i < s.size(); ++i) { if (i) f << ","; f << s[i]; }
        f << "]";
    };

    f << std::fixed << std::setprecision(10);
    f << "{\n";

    // ── legends ─────────────────────────────────────────────────────────────
    f << "  \"legends\": {\n";

    f << "    \"cluster_types\": [\n";
    for (size_t i = 0; i < kSnapClusterNames.size(); ++i) {
        f << "      {\"id\": " << i
          << ", \"name\": \"" << kSnapClusterNames[i] << "\""
          << ", \"rgb\": ";
        write_rgb_i(kSnapClusterRgb[i]);
        f << "}";
        if (i + 1 < kSnapClusterNames.size()) f << ",";
        f << "\n";
    }
    f << "    ],\n";

    f << "    \"edge_topo_types\": [\n";
    for (size_t i = 0; i < kSnapEdgeTopoNames.size(); ++i) {
        f << "      {\"id\": " << i
          << ", \"name\": \"" << kSnapEdgeTopoNames[i] << "\""
          << ", \"rgb\": ";
        write_rgb_i(kSnapEdgeTopoRgb[i]);
        f << "}";
        if (i + 1 < kSnapEdgeTopoNames.size()) f << ",";
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
    for (size_t i = 0; i < snap.vertices.size(); ++i) {
        if (i) f << ",\n";
        const auto& p = snap.vertices[i];
        f << "    {\"pos\": [" << p[0] << "," << p[1] << "," << p[2] << "]"
          << ", \"struct_ids\": ";
        if (i < snap.vertex_struct_ids.size()) write_intvec(snap.vertex_struct_ids[i]);
        else                                    f << "[]";
        f << ", \"cluster_type\": "
          << (int)(i < snap.vertex_cluster_type.size() ? snap.vertex_cluster_type[i] : 0);
        f << ", \"original_ancestors\": [";
        if (i < snap.vertex_original_ancestors.size()) {
            const auto& anc = snap.vertex_original_ancestors[i];
            for (size_t j = 0; j < anc.size(); ++j) {
                if (j) f << ",";
                f << anc[j];
            }
        }
        f << "]}";
    }
    f << "\n  ],\n";

    // ── edges ───────────────────────────────────────────────────────────────
    f << "  \"edges\": [\n";
    for (size_t i = 0; i < snap.edges.size(); ++i) {
        if (i) f << ",\n";
        const auto& e = snap.edges[i];
        f << "    {\"v\": [" << e[0] << "," << e[1] << "]"
          << ", \"struct_ids\": ";
        if (i < snap.edge_struct_ids.size()) write_intvec(snap.edge_struct_ids[i]);
        else                                  f << "[]";
        f << ", \"topo_type\": "
          << (int)(i < snap.edge_topo_type.size() ? snap.edge_topo_type[i] : 0);
        if (i < snap.edge_last_rejection.size() && snap.edge_last_rejection[i] != 255)
            f << ", \"rejection_reason\": " << (int)snap.edge_last_rejection[i];
        else
            f << ", \"rejection_reason\": null";
        f << "}";
    }
    f << "\n  ],\n";

    // ── faces ───────────────────────────────────────────────────────────────
    f << "  \"faces\": [\n";
    for (size_t i = 0; i < snap.faces.size(); ++i) {
        if (i) f << ",\n";
        const auto& fc = snap.faces[i];
        f << "    {\"v\": [" << fc[0] << "," << fc[1] << "," << fc[2] << "]"
          << ", \"struct_id\": "
          << (i < snap.face_struct_id.size() ? snap.face_struct_id[i] : -1)
          << "}";
    }
    f << "\n  ],\n";

    // ── original positions ──────────────────────────────────────────────────
    f << "  \"original_positions\": [\n";
    for (size_t i = 0; i < snap.original_positions.size(); ++i) {
        if (i) f << ",\n";
        const auto& p = snap.original_positions[i];
        f << "    [" << p[0] << "," << p[1] << "," << p[2] << "]";
    }
    f << "\n  ]\n";

    f << "}\n";

    std::cout << "[ExportSnapshotVisualizeInfo] wrote "
              << snap.vertices.size() << " verts, "
              << snap.edges.size() << " edges, "
              << snap.faces.size() << " faces to " << path << "\n";
}

void ExportSnapshotMatTyped(const VcgDirectSnapshot& snap, const std::string& path)
{
    std::ofstream f(path);
    if (!f) {
        std::cerr << "[ExportSnapshotMatTyped] cannot open: " << path << "\n";
        return;
    }
    f << std::fixed << std::setprecision(15);
    f << snap.vertices.size() << " " << snap.edges.size() << " " << snap.faces.size() << "\n";
    for (size_t i = 0; i < snap.vertices.size(); ++i) {
        const auto& c = snap.vertices[i];
        const double r = (i < snap.vertex_radius.size()) ? snap.vertex_radius[i] : 0.0;
        const uint8_t t = (i < snap.vertex_cluster_type.size()) ? snap.vertex_cluster_type[i] : 7;
        f << "v " << c[0] << " " << c[1] << " " << c[2]
          << " " << r << " " << SnapClusterName(t) << "\n";
    }
    for (const auto& e : snap.edges)
        f << "e " << e[0] << " " << e[1] << "\n";
    for (const auto& fc : snap.faces)
        f << "f " << fc[0] << " " << fc[1] << " " << fc[2] << "\n";
    std::cerr << "[ExportSnapshotMatTyped] wrote "
              << snap.vertices.size() << " verts, "
              << snap.edges.size() << " edges, "
              << snap.faces.size() << " faces to " << path << "\n";
}

void ExportSnapshotMa(const VcgDirectSnapshot& snap, const std::string& path_prefix)
{
    std::string fname = path_prefix;
    fname += "___v_" + std::to_string(snap.vertices.size());
    fname += "___e_" + std::to_string(snap.edges.size());
    fname += "___f_" + std::to_string(snap.faces.size());
    fname += ".ma";

    std::ofstream f(fname);
    if (!f) {
        std::cerr << "[ExportSnapshotMa] cannot open: " << fname << "\n";
        return;
    }
    f << snap.vertices.size() << " " << snap.edges.size() << " " << snap.faces.size() << "\n";
    f << std::fixed << std::setprecision(15);
    for (size_t i = 0; i < snap.vertices.size(); ++i) {
        const auto& c = snap.vertices[i];
        const double r = (i < snap.vertex_radius.size()) ? snap.vertex_radius[i] : 0.0;
        f << "v " << c[0] << " " << c[1] << " " << c[2] << " " << r << "\n";
    }
    for (const auto& e : snap.edges)
        f << "e " << e[0] << " " << e[1] << "\n";
    for (const auto& fc : snap.faces)
        f << "f " << fc[0] << " " << fc[1] << " " << fc[2] << "\n";
    std::cout << "[ExportSnapshotMa] wrote " << fname << "\n";
}

void ExportSnapshotClusterPLY(const VcgDirectSnapshot& snap, const std::string& path)
{
    std::ofstream ply(path);
    if (!ply) {
        std::cerr << "[ExportSnapshotClusterPLY] cannot open: " << path << "\n";
        return;
    }
    ply << "ply\nformat ascii 1.0\n"
        << "element vertex " << snap.vertices.size() << "\n"
        << "property float x\nproperty float y\nproperty float z\n"
        << "property uchar red\nproperty uchar green\nproperty uchar blue\n"
        << "element face " << snap.faces.size() << "\n"
        << "property list uchar int vertex_indices\n"
        << "element edge " << snap.edges.size() << "\n"
        << "property int vertex1\nproperty int vertex2\n"
        << "end_header\n";
    ply << std::fixed << std::setprecision(10);
    for (size_t i = 0; i < snap.vertices.size(); ++i) {
        const auto& c = snap.vertices[i];
        const uint8_t ct = (i < snap.vertex_cluster_type.size()) ? snap.vertex_cluster_type[i] : 7;
        const auto rgb = SnapClusterRgb(ct);
        ply << c[0] << " " << c[1] << " " << c[2] << " "
            << rgb[0] << " " << rgb[1] << " " << rgb[2] << "\n";
    }
    for (const auto& fc : snap.faces)
        ply << "3 " << fc[0] << " " << fc[1] << " " << fc[2] << "\n";
    for (const auto& e : snap.edges)
        ply << e[0] << " " << e[1] << "\n";
    std::cerr << "[ExportSnapshotClusterPLY] wrote "
              << snap.vertices.size() << " verts, "
              << snap.faces.size() << " faces, "
              << snap.edges.size() << " edges to " << path << "\n";
}

void ExportSnapshotRejectionSkeleton(const VcgDirectSnapshot& snap,
                                     const std::string& path,
                                     double radius_frac)
{
    // ── bbox-derived cylinder radius (matches SlabMesh::ExportSkeletonPLY) ──
    double minx =  std::numeric_limits<double>::max();
    double miny =  std::numeric_limits<double>::max();
    double minz =  std::numeric_limits<double>::max();
    double maxx = -std::numeric_limits<double>::max();
    double maxy = -std::numeric_limits<double>::max();
    double maxz = -std::numeric_limits<double>::max();
    for (const auto& p : snap.vertices) {
        minx = std::min(minx, p[0]); maxx = std::max(maxx, p[0]);
        miny = std::min(miny, p[1]); maxy = std::max(maxy, p[1]);
        minz = std::min(minz, p[2]); maxz = std::max(maxz, p[2]);
    }
    const double dx = maxx - minx, dy = maxy - miny, dz = maxz - minz;
    const double diag = std::sqrt(dx*dx + dy*dy + dz*dz);
    const double radius = (diag > 1e-12 ? diag : 1.0) * radius_frac;

    struct V { float x, y, z; uint8_t r, g, b, a; };
    struct T { int a, b, c; };
    std::vector<V> verts;
    std::vector<T> tris;
    verts.reserve(snap.edges.size() * 16);
    tris.reserve(snap.edges.size() * 16);

    const int N = 8;
    const double pi = std::acos(-1.0);

    auto reasonCol = [](uint8_t rr) -> std::array<uint8_t,3> {
        if (rr == 255) return {255, 255, 255};
        return SlabMesh::RejectionReasonColorU8(static_cast<SlabMesh::RejectionReason>(rr));
    };

    // ── cylinder per edge ──────────────────────────────────────────────────
    for (size_t i = 0; i < snap.edges.size(); ++i) {
        const auto& e = snap.edges[i];
        if (e[0] < 0 || e[1] < 0 ||
            (size_t)e[0] >= snap.vertices.size() ||
            (size_t)e[1] >= snap.vertices.size()) continue;
        const auto& p0 = snap.vertices[e[0]];
        const auto& p1 = snap.vertices[e[1]];

        double dxe = p1[0]-p0[0], dye = p1[1]-p0[1], dze = p1[2]-p0[2];
        double len = std::sqrt(dxe*dxe + dye*dye + dze*dze);
        if (len < 1e-10) continue;
        dxe /= len; dye /= len; dze /= len;

        // Reference vector for the perpendicular basis.
        double rx = (std::abs(dxe) < 0.9) ? 1.0 : 0.0;
        double ry = (std::abs(dxe) < 0.9) ? 0.0 : 1.0;
        double rz = 0.0;
        // u = d x ref (then normalize), v = d x u
        double ux = dye*rz - dze*ry;
        double uy = dze*rx - dxe*rz;
        double uz = dxe*ry - dye*rx;
        double ulen = std::sqrt(ux*ux + uy*uy + uz*uz);
        if (ulen < 1e-12) continue;
        ux /= ulen; uy /= ulen; uz /= ulen;
        double vx = dye*uz - dze*uy;
        double vy = dze*ux - dxe*uz;
        double vz = dxe*uy - dye*ux;

        const uint8_t rr = (i < snap.edge_last_rejection.size())
                           ? snap.edge_last_rejection[i] : 255;
        const auto col = reasonCol(rr);

        const int base = (int)verts.size();
        for (int s = 0; s < N; ++s) {
            const double angle = 2.0 * pi * s / N;
            const double cs = radius * std::cos(angle);
            const double sn = radius * std::sin(angle);
            const double ox = ux*cs + vx*sn;
            const double oy = uy*cs + vy*sn;
            const double oz = uz*cs + vz*sn;
            verts.push_back({ (float)(p0[0]+ox), (float)(p0[1]+oy), (float)(p0[2]+oz),
                              col[0], col[1], col[2], 255 });
            verts.push_back({ (float)(p1[0]+ox), (float)(p1[1]+oy), (float)(p1[2]+oz),
                              col[0], col[1], col[2], 255 });
        }
        for (int s = 0; s < N; ++s) {
            const int sn = (s + 1) % N;
            const int a = base + s  * 2;
            const int b = base + s  * 2 + 1;
            const int c = base + sn * 2;
            const int dd = base + sn * 2 + 1;
            tris.push_back({ a, b, dd });
            tris.push_back({ a, dd, c });
        }
    }

    // ── semi-transparent face triangles ────────────────────────────────────
    const uint8_t fr = 180, fg = 180, fb = 200, fa = 60;
    unsigned face_tri_count = 0;
    for (const auto& fc : snap.faces) {
        if (fc[0] < 0 || fc[1] < 0 || fc[2] < 0) continue;
        if ((size_t)fc[0] >= snap.vertices.size() ||
            (size_t)fc[1] >= snap.vertices.size() ||
            (size_t)fc[2] >= snap.vertices.size()) continue;
        const auto& a = snap.vertices[fc[0]];
        const auto& b = snap.vertices[fc[1]];
        const auto& c = snap.vertices[fc[2]];
        const int base = (int)verts.size();
        verts.push_back({ (float)a[0], (float)a[1], (float)a[2], fr, fg, fb, fa });
        verts.push_back({ (float)b[0], (float)b[1], (float)b[2], fr, fg, fb, fa });
        verts.push_back({ (float)c[0], (float)c[1], (float)c[2], fr, fg, fb, fa });
        tris.push_back({ base, base + 1, base + 2 });
        ++face_tri_count;
    }

    std::ofstream ply(path);
    if (!ply) {
        std::cerr << "[ExportSnapshotRejectionSkeleton] cannot open: " << path << "\n";
        return;
    }
    ply << "ply\nformat ascii 1.0\n"
        << "comment VDE rejection skeleton -- vertex colour encodes last collapse rejection reason\n"
        << "comment alpha=255 -> cylinder (opaque)  alpha=60 -> MAT face (semi-transparent)\n"
        << "element vertex " << verts.size() << "\n"
        << "property float x\nproperty float y\nproperty float z\n"
        << "property uchar red\nproperty uchar green\nproperty uchar blue\nproperty uchar alpha\n"
        << "element face " << tris.size() << "\n"
        << "property list uchar int vertex_indices\n"
        << "end_header\n";
    ply << std::fixed << std::setprecision(6);
    for (const auto& vt : verts)
        ply << vt.x << " " << vt.y << " " << vt.z << " "
            << (int)vt.r << " " << (int)vt.g << " " << (int)vt.b << " " << (int)vt.a << "\n";
    for (const auto& tr : tris)
        ply << "3 " << tr.a << " " << tr.b << " " << tr.c << "\n";
    std::cerr << "[ExportSnapshotRejectionSkeleton] "
              << (verts.size() / (N * 2)) << " cylinders, "
              << face_tri_count << " MAT faces to " << path << "\n";
}

#endif  // QMAT_WITH_VCGLIB

#if defined(QMAT_WITH_POLYSCOPE) && defined(QMAT_WITH_VCGLIB)

#include <chrono>
#include <thread>
#include <unordered_map>

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"
#include "polyscope/pick.h"
#include "imgui.h"

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

    // MAT Rejection Edges — every edge coloured by its last rejection reason
    // (255 / "never attempted" → white).  Off by default; user toggles from the
    // layer panel.  Pick space: [nodes, edges), edge slot = local_idx - V.
    if (!snap.edges.empty() && !snap.edge_last_rejection.empty()) {
        std::vector<E2> redges; redges.reserve(snap.edges.size());
        std::vector<C3> rcolors; rcolors.reserve(snap.edges.size());
        for (size_t i = 0; i < snap.edges.size(); ++i) {
            redges.push_back({ (size_t)snap.edges[i][0], (size_t)snap.edges[i][1] });
            uint8_t r = (i < snap.edge_last_rejection.size())
                        ? snap.edge_last_rejection[i] : 255;
            if (r == 255) {
                rcolors.push_back({1.0f, 1.0f, 1.0f});
            } else {
                auto rgb = SlabMesh::RejectionReasonColorU8(
                    static_cast<SlabMesh::RejectionReason>(r));
                rcolors.push_back({rgb[0]/255.0f, rgb[1]/255.0f, rgb[2]/255.0f});
            }
        }
        bool en = ps::hasCurveNetwork("MAT Rejection Edges")
                  ? ps::getCurveNetwork("MAT Rejection Edges")->isEnabled() : false;
        auto* cn = ps::registerCurveNetwork("MAT Rejection Edges", snap.vertices, redges);
        cn->setRadius(0.0010f, true);
        cn->addEdgeColorQuantity("Rejection Reason", rcolors)->setEnabled(true);
        cn->setEnabled(en);
    } else if (ps::hasCurveNetwork("MAT Rejection Edges")) {
        ps::removeStructure("MAT Rejection Edges");
    }

    ps::frameTick();
}

// Disable the three click-driven rejection primitives without removing them.
void ClearVdeRejectionPrimitives()
{
    namespace ps = polyscope;
    if (ps::hasPointCloud("Rejection Verts"))
        ps::getPointCloud("Rejection Verts")->setEnabled(false);
    if (ps::hasCurveNetwork("Rejection Edges"))
        ps::getCurveNetwork("Rejection Edges")->setEnabled(false);
    if (ps::hasPointCloud("Rejection Target"))
        ps::getPointCloud("Rejection Target")->setEnabled(false);
}

// Draw the offending primitives for a selected rejected edge.
// Phase D 1.4 scope: verts (endpoints) + edge + target only — no faces, no
// flipped triangles, no spheres.
void ShowVdeRejectionPrimitives(const VcgDirectSnapshot& snap, unsigned eid)
{
    namespace ps = polyscope;
    if (eid >= snap.edges.size()) { ClearVdeRejectionPrimitives(); return; }

    const int a = snap.edges[eid][0];
    const int b = snap.edges[eid][1];
    if (a < 0 || b < 0 ||
        (size_t)a >= snap.vertices.size() || (size_t)b >= snap.vertices.size())
    {
        ClearVdeRejectionPrimitives();
        return;
    }

    {
        std::vector<std::array<double,3>> pts = { snap.vertices[a], snap.vertices[b] };
        auto* pc = ps::registerPointCloud("Rejection Verts", pts);
        pc->setPointColor(glm::vec3(1.0f, 1.0f, 0.0f));
        pc->setPointRadius(0.006, true);
        pc->setEnabled(true);
    }
    {
        std::vector<std::array<double,3>> nodes = { snap.vertices[a], snap.vertices[b] };
        std::vector<std::array<size_t,2>> segs = {{0, 1}};
        auto* cn = ps::registerCurveNetwork("Rejection Edges", nodes, segs);
        cn->setColor(glm::vec3(1.0f, 1.0f, 0.0f));
        cn->setRadius(0.004f, true);
        cn->setEnabled(true);
    }
    if (eid < snap.edge_rejection_target.size()) {
        std::vector<std::array<double,3>> pts = { snap.edge_rejection_target[eid] };
        auto* pc = ps::registerPointCloud("Rejection Target", pts);
        pc->setPointColor(glm::vec3(0.0f, 1.0f, 0.0f));
        pc->setPointRadius(0.008, true);
        pc->setEnabled(true);
    } else if (ps::hasPointCloud("Rejection Target")) {
        ps::getPointCloud("Rejection Target")->setEnabled(false);
    }
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
            // MAT Rejection Edges pick → highlight the offending cause.
            // Curve-network pick space: [nodes, edges); subtract V for the edge slot.
            else if (polyscope::hasCurveNetwork("MAT Rejection Edges") &&
                     struct_ptr == polyscope::getCurveNetwork("MAT Rejection Edges"))
            {
                if (local_idx >= snap.vertices.size()) {
                    size_t edge_slot = local_idx - snap.vertices.size();
                    if (edge_slot < snap.edges.size() &&
                        (int)edge_slot != vs.selected_rejection_eid)
                    {
                        vs.selected_rejection_eid = (int)edge_slot;
                        vs.selected_vid = -1;
                        vs.selected_eid = -1;
                        vs.selected_init_fid = -1;
                        ClearVdeRejectionPrimitives();
                        ShowVdeRejectionPrimitives(snap, (unsigned)edge_slot);
                    }
                }
                polyscope::pick::resetSelection();
            }
        }

        // Auto-clear rejection selection if the picked edge no longer exists.
        if (vs.selected_rejection_eid >= 0 &&
            (size_t)vs.selected_rejection_eid >= snap.edges.size())
        {
            vs.selected_rejection_eid = -1;
            ClearVdeRejectionPrimitives();
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

        // Rejection-edge selection info — reason name + (for struct_ids /
        // cluster mismatches) the triplet of sets so the cause is visible.
        if (vs.selected_rejection_eid >= 0) {
            const int eid = vs.selected_rejection_eid;
            ImGui::Text("Rejection edge: %d", eid);
            const uint8_t reason = (eid < (int)snap.edge_last_rejection.size())
                                   ? snap.edge_last_rejection[eid] : 255;
            if (reason == 255) {
                ImGui::Text("  Reason: (none / not attempted)");
            } else {
                auto rr = static_cast<SlabMesh::RejectionReason>(reason);
                ImGui::Text("  Reason: %s", SlabMesh::RejectionReasonName(rr));

                auto setStr = [](const std::vector<int>& s) -> std::string {
                    if (s.empty()) return "{}";
                    std::string out = "{";
                    for (int id : s) { if (out.size() > 1) out += ","; out += std::to_string(id); }
                    return out + "}";
                };

                const int a = snap.edges[eid][0];
                const int b = snap.edges[eid][1];
                if (rr == SlabMesh::RejectionReason::struct_ids_sets_different) {
                    if (eid < (int)snap.edge_struct_ids.size())
                        ImGui::Text("    edge  struct_ids: %s",
                            setStr(snap.edge_struct_ids[eid]).c_str());
                    if (a >= 0 && a < (int)snap.vertex_struct_ids.size())
                        ImGui::Text("    v%d  struct_ids: %s", a,
                            setStr(snap.vertex_struct_ids[a]).c_str());
                    if (b >= 0 && b < (int)snap.vertex_struct_ids.size())
                        ImGui::Text("    v%d  struct_ids: %s", b,
                            setStr(snap.vertex_struct_ids[b]).c_str());
                } else if (rr == SlabMesh::RejectionReason::DifferentClusterType) {
                    if (a >= 0 && a < (int)snap.vertex_cluster_type.size())
                        ImGui::Text("    v%d  topo_type: %s", a,
                            ClusterTypeName(snap.vertex_cluster_type[a]));
                    if (b >= 0 && b < (int)snap.vertex_cluster_type.size())
                        ImGui::Text("    v%d  topo_type: %s", b,
                            ClusterTypeName(snap.vertex_cluster_type[b]));
                }
            }
            if (ImGui::Button("Clear rejection selection")) {
                vs.selected_rejection_eid = -1;
                ClearVdeRejectionPrimitives();
            }
            ImGui::Separator();
        }

        if (ImGui::Button("Export MAT as OFF")) {
            std::string path = vs.outputPrefix
                + "_snapshot_" + std::to_string(vs.collapse_count) + ".off";
            ExportSnapshotAsOff(self.latest_snap, path);
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
