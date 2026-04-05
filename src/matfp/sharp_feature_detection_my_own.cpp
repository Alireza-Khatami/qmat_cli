#include "sharp_feature_detection.h"

#include <cmath>
#include <map>

static const double SFD_PI = 3.14159265358979323846;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

using HCH = MPMesh::Halfedge_const_handle;
using FCH = MPMesh::Facet_const_handle;
using Vec = MPMesh::Vector;
using Pt  = MPMesh::Point;

// Compute the unit outward normal of the face that halfedge h belongs to.
// Returns a zero vector when the face is degenerate.
static Vec face_normal(FCH f)
{
    HCH h  = f->halfedge();
    Pt  p0 = h->vertex()->point();
    Pt  p1 = h->next()->vertex()->point();
    Pt  p2 = h->next()->next()->vertex()->point();
    Vec n  = CGAL::cross_product(p1 - p0, p2 - p0);
    double len = std::sqrt(n.squared_length());
    if (len < 1e-15) return Vec(0.0, 0.0, 0.0);
    return n / len;
}

// Compute the centroid of face f.
static Pt face_centroid(FCH f)
{
    HCH h  = f->halfedge();
    Pt  p0 = h->vertex()->point();
    Pt  p1 = h->next()->vertex()->point();
    Pt  p2 = h->next()->next()->vertex()->point();
    return Pt((p0.x() + p1.x() + p2.x()) / 3.0,
              (p0.y() + p1.y() + p2.y()) / 3.0,
              (p0.z() + p1.z() + p2.z()) / 3.0);
}

// ---------------------------------------------------------------------------
// Main function
// ---------------------------------------------------------------------------

void find_feature_edges(
    const MPMesh& mesh,
    double thres_concave,
    double angle_sharp_deg,
    std::set<std::pair<int,int>>& sharp_edges,
    std::set<std::pair<int,int>>& concave_edges,
    std::set<int>& corners)
{
    sharp_edges.clear();
    concave_edges.clear();
    corners.clear();

    const double angle_sharp_rad = angle_sharp_deg * SFD_PI / 180.0;

    // Iterate over all halfedges; process each undirected edge exactly once
    // by only handling the halfedge whose address is smaller than its opposite.
    for (auto h = mesh.halfedges_begin(); h != mesh.halfedges_end(); ++h)
    {
        // Canonical representative: skip the "larger" halfedge of each pair
        if (&*h > &*h->opposite()) continue;

        // Skip border (open boundary) edges
        if (h->is_border() || h->opposite()->is_border()) continue;

        FCH fa = h->facet();
        FCH fb = h->opposite()->facet();

        Vec n_a = face_normal(fa);
        Vec n_b = face_normal(fb);
        Pt  c_a = face_centroid(fa);
        Pt  c_b = face_centroid(fb);

        int id_a = h->vertex()->id;
        int id_b = h->opposite()->vertex()->id;
        auto edge = std::make_pair(std::min(id_a, id_b), std::max(id_a, id_b));

        // ------------------------------------------------------------------
        // Concave detection
        //
        // If (c_b - c_a).normalized · n_a > threshold  (A looking toward B)
        // or (c_a - c_b).normalized · n_b > threshold  (B looking toward A)
        // the edge is concave (dihedral > 180°).
        // ------------------------------------------------------------------
        Vec diff_ab = c_b - c_a;
        double dist = std::sqrt(diff_ab.squared_length());

        double tmp_concave   = 0.0;
        double tmp_concave_2 = 0.0;
        if (dist > 1e-15)
        {
            Vec d = diff_ab / dist;
            tmp_concave   =  d * n_a;   // (c_b-c_a)/|…| · n_a
            tmp_concave_2 = -d * n_b;   // (c_a-c_b)/|…| · n_b
        }

        if (tmp_concave > thres_concave || tmp_concave_2 > thres_concave)
        {
            concave_edges.insert(edge);
        }
        else
        {
            // ------------------------------------------------------------------
            // Sharp (convex) detection
            //
            // Dihedral angle theta = acos(n_a · n_b).
            // acos range is [0, π].  Edge is sharp when theta > angle_sharp_rad.
            // ------------------------------------------------------------------
            double dot = n_a * n_b;
            dot = std::max(-1.0, std::min(1.0, dot));   // numerical clamp
            double theta = std::acos(dot);
            if (theta > angle_sharp_rad && theta < SFD_PI)
            {
                sharp_edges.insert(edge);
            }
        }
    }

    // ------------------------------------------------------------------
    // Corner detection: vertices incident to >= 3 sharp edges
    // ------------------------------------------------------------------
    std::map<int, std::set<int>> neighbor_v;
    for (const auto& e : sharp_edges)
    {
        neighbor_v[e.first].insert(e.second);
        neighbor_v[e.second].insert(e.first);
    }
    for (const auto& kv : neighbor_v)
    {
        if (kv.second.size() > 2)
            corners.insert(kv.first);
    }
}
