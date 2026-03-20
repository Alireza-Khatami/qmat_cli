// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
#include "sharp_feature_detection.h"

#include <cmath>
#include <cstdio>
#include <map>
#include <stdexcept>
#include <unordered_set>

using namespace matfp;

namespace pre_matfp {

void find_feature_edges(const Args& args,
                        const std::vector<Vector3>& input_vertices,
                        const std::vector<Vector3i>& input_faces,
                        std::set<std::array<int, 2>>& s_edges,
                        std::set<std::array<int, 2>>& cc_edges,
                        std::set<int>& corners) {
  fprintf(stderr,
          "[sharp_feature_detection] Detecting sharp/concave edges and corners "
          "using threshold: %f\n",
          args.thres_concave);
  s_edges.clear();
  cc_edges.clear();
  corners.clear();

  std::vector<std::array<int, 2>> edges;
  std::map<int, std::unordered_set<int>> conn_tris;
  for (int i = 0; i < (int)input_faces.size(); i++) {
    const auto& f = input_faces[i];
    for (int j = 0; j < 3; j++) {
      std::array<int, 2> e = {{f[j], f[(j + 1) % 3]}};
      if (e[0] > e[1]) std::swap(e[0], e[1]);
      edges.push_back(e);
      conn_tris[input_faces[i][j]].insert(i);
    }
  }
  vector_unique(edges);

  // find sharp edges and concave edges
  for (const auto& e : edges) {
    std::vector<int> n12_f_ids;
    set_intersection(conn_tris[e[0]], conn_tris[e[1]], n12_f_ids);

    if (n12_f_ids.size() == 1) {  // open boundary
      fprintf(stderr,
              "[sharp_feature_detection] ERROR: Detect open boundary! "
              "edge (%d,%d) has only 1 adjacent face.\n",
              e[0], e[1]);
      throw std::runtime_error(
          "ERROR: we don't know how to handle open boundary!!");
    }
    int f_id = n12_f_ids[0];
    int j = 0;
    for (; j < 3; j++) {
      if ((input_faces[f_id][j] == e[0] &&
           input_faces[f_id][mod3(j + 1)] == e[1]) ||
          (input_faces[f_id][j] == e[1] &&
           input_faces[f_id][mod3(j + 1)] == e[0]))
        break;
    }
    Vector3 n = get_normal(input_vertices[input_faces[f_id][0]],
                           input_vertices[input_faces[f_id][1]],
                           input_vertices[input_faces[f_id][2]]);
    Vector3 c_n = get_triangle_centroid(input_vertices[input_faces[f_id][0]],
                                        input_vertices[input_faces[f_id][1]],
                                        input_vertices[input_faces[f_id][2]]);

    for (int k = 0; k < (int)n12_f_ids.size(); k++) {
      if (n12_f_ids[k] == f_id) continue;
      Vector3 n1 = get_normal(input_vertices[input_faces[n12_f_ids[k]][0]],
                              input_vertices[input_faces[n12_f_ids[k]][1]],
                              input_vertices[input_faces[n12_f_ids[k]][2]]);
      Vector3 c_n1 =
          get_triangle_centroid(input_vertices[input_faces[n12_f_ids[k]][0]],
                                input_vertices[input_faces[n12_f_ids[k]][1]],
                                input_vertices[input_faces[n12_f_ids[k]][2]]);

      std::array<int, 2> ref_fs_pair = {{f_id, n12_f_ids[k]}};
      std::sort(ref_fs_pair.begin(), ref_fs_pair.end());
      std::array<Vector3, 2> ref_fs_normals = {{n, n1}};
      bool is_debug = false;

      //////////
      // Since cosine can only measure dihedral angle from (0, 180)
      // but concave has angle larger than 180
      // therefore we use different measurement for concave, and sharp edges
      //////////
      // Concave edges
      // c_n is a random vertex on plane A, c_n1 is a random vertex on plane B
      // n is normal of A
      //
      // 2021-09-04 ninwang:
      // If na and nb are the normals of the both adjacent faces,
      // and pa and pb vertices of the both faces that are not connected to
      // the edge, wherein na and pa belongs to the face A, and nb and pb to
      // the face B, then ( pb - pa ) . na > 0 => concave edge
      double tmp_concave = (c_n1 - c_n).normalized().dot(n);     // A, B
      double tmp_concave_2 = (c_n - c_n1).normalized().dot(n1);  // B, A
      if (tmp_concave > args.thres_concave ||
          tmp_concave_2 > args.thres_concave) {  // SCALAR_ZERO is too small
        if (is_debug)
          fprintf(stderr,
                  "[sharp_feature_detection] edge (%d,%d) is a concave edge, "
                  "tmp_concave: %f\n",
                  e[0], e[1], tmp_concave);
        cc_edges.insert(e);  // once concave, never sharp
      } else {
        // Sharp edges (when it's convex)
        // angle between two normals of convex faces: theta
        // => cos(theta) = n1.dot(n)
        // acosine() range in [0, pi]
        // sharp edges => theta in (angle_sharp, pi)
        // here angle_sharp = 30
        // Note that, using theta CANNOT differentiate concave or convex,
        // so the concave detection must run first
        double tmp_convex = std::acos(n1.dot(n));
        double angle_sharp = PI * (args.thres_convex / 180.);
        if (angle_sharp < tmp_convex && tmp_convex < PI) {
          // fprintf(stderr, "[sharp_feature_detection] sharp edge: theta is %f\n", tmp_convex);
          s_edges.insert(e);
        }
      }
    }  // for n12_f_ids
  }    // for edges

  // vector_unique(s_edges);

  // find corners
  // connect to at least 3 sharp edges
  // (not including concave edges)
  std::map<int, std::set<int>> neighbor_v;
  for (const auto& e : s_edges) {
    neighbor_v[e[0]].insert(e[1]);
    neighbor_v[e[1]].insert(e[0]);
  }
  for (const auto& pair : neighbor_v) {
    // Found a corner
    if (pair.second.size() > 2) {
      corners.insert(pair.first);
    }
  }

  fprintf(stderr, "[sharp_feature_detection] #concave_edges = %zu\n",
          cc_edges.size());
  fprintf(stderr, "[sharp_feature_detection] #sharp_edges = %zu\n",
          s_edges.size());
  fprintf(stderr, "[sharp_feature_detection] #corners = %zu\n",
          corners.size());
}

}  // namespace pre_matfp
