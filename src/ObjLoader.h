#ifndef OBJ_LOADER_H
#define OBJ_LOADER_H

#include <string>
#include <vector>
#include <iostream>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include "Mesh.h"

// Builder class for constructing CGAL Polyhedron from OBJ data
// Template parameter PointType allows using different point types for different kernels
template <class HDS, class PointType>
class ObjPolyhedronBuilderT : public CGAL::Modifier_base<HDS> {
public:
    std::vector<double> vertices;  // x, y, z triplets
    std::vector<std::vector<int>> faces;  // face vertex indices
    std::string error;

    ObjPolyhedronBuilderT() {}

    void operator()(HDS& hds) {
        CGAL::Polyhedron_incremental_builder_3<HDS> builder(hds, true);

        size_t numVertices = vertices.size() / 3;
        size_t numFaces = faces.size();

        builder.begin_surface(numVertices, numFaces);

        // Add vertices
        for (size_t i = 0; i < numVertices; ++i) {
            builder.add_vertex(PointType(
                vertices[i * 3],
                vertices[i * 3 + 1],
                vertices[i * 3 + 2]
            ));
        }

        // Add faces
        for (size_t i = 0; i < numFaces; ++i) {
            if (builder.test_facet(faces[i].begin(), faces[i].end())) {
                builder.add_facet(faces[i].begin(), faces[i].end());
            } else {
                std::cerr << "Warning: Skipping invalid facet " << i << std::endl;
            }
        }

        builder.end_surface();
    }
};

// Type alias for MPMesh (uses simple_kernel Point)
template <class HDS>
using ObjPolyhedronBuilder = ObjPolyhedronBuilderT<HDS, Point>;

// Load an OBJ file into a CGAL Polyhedron mesh (MPMesh type)
// Returns true on success, false on failure
// Error message is stored in 'error' parameter
bool LoadObjFile(const std::string& filename, Mesh& mesh, std::string& error);

// Load an OBJ file into a basic CGAL Polyhedron (for mesh domain)
// Returns true on success, false on failure
// Error message is stored in 'error' parameter
bool LoadObjFile(const std::string& filename, Polyhedron& mesh, std::string& error);

// Check if a filename has .obj extension (case-insensitive)
bool IsObjFile(const std::string& filename);

// Check if a filename has .off extension (case-insensitive)
bool IsOffFile(const std::string& filename);

// Get file extension in lowercase
std::string GetFileExtension(const std::string& filename);

#endif // OBJ_LOADER_H
