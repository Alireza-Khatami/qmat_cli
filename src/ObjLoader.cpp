#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"

#include "ObjLoader.h"
#include <algorithm>
#include <cctype>

std::string GetFileExtension(const std::string& filename) {
    size_t dotPos = filename.rfind('.');
    if (dotPos == std::string::npos || dotPos == filename.length() - 1) {
        return "";
    }
    std::string ext = filename.substr(dotPos + 1);
    // Convert to lowercase
    std::transform(ext.begin(), ext.end(), ext.begin(),
        [](unsigned char c) { return std::tolower(c); });
    return ext;
}

bool IsObjFile(const std::string& filename) {
    return GetFileExtension(filename) == "obj";
}

bool IsOffFile(const std::string& filename) {
    return GetFileExtension(filename) == "off";
}

// Internal struct to hold parsed OBJ data
struct ObjData {
    std::vector<double> vertices;
    std::vector<std::vector<int>> faces;
};

// Internal function to parse OBJ file into raw data
static bool ParseObjFile(const std::string& filename, ObjData& data, std::string& error) {
    tinyobj::attrib_t attrib;
    std::vector<tinyobj::shape_t> shapes;
    std::vector<tinyobj::material_t> materials;
    std::string warn;
    std::string err;

    // Extract base directory for material files
    std::string baseDir = "";
    size_t lastSlash = filename.find_last_of("/\\");
    if (lastSlash != std::string::npos) {
        baseDir = filename.substr(0, lastSlash + 1);
    }

    bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &err,
                                 filename.c_str(), baseDir.empty() ? nullptr : baseDir.c_str(),
                                 true);  // triangulate = true

    if (!warn.empty()) {
        std::cerr << "TinyObjLoader warning: " << warn << std::endl;
    }

    if (!err.empty()) {
        error = err;
    }

    if (!ret) {
        if (error.empty()) {
            error = "Failed to load OBJ file: " + filename;
        }
        return false;
    }

    if (attrib.vertices.empty()) {
        error = "OBJ file contains no vertices";
        return false;
    }

    // Copy vertices (tinyobj stores as x,y,z triplets)
    // Convert from tinyobj::real_t (float) to double
    data.vertices.reserve(attrib.vertices.size());
    for (const auto& v : attrib.vertices) {
        data.vertices.push_back(static_cast<double>(v));
    }

    // Collect all faces from all shapes
    for (const auto& shape : shapes) {
        size_t indexOffset = 0;
        for (size_t f = 0; f < shape.mesh.num_face_vertices.size(); ++f) {
            int fv = shape.mesh.num_face_vertices[f];

            std::vector<int> faceIndices;
            faceIndices.reserve(fv);

            for (int v = 0; v < fv; ++v) {
                tinyobj::index_t idx = shape.mesh.indices[indexOffset + v];
                faceIndices.push_back(idx.vertex_index);
            }

            data.faces.push_back(faceIndices);
            indexOffset += fv;
        }
    }

    if (data.faces.empty()) {
        error = "OBJ file contains no faces";
        return false;
    }

    std::cout << "  OBJ: " << data.vertices.size() / 3 << " vertices, "
              << data.faces.size() << " faces" << std::endl;

    return true;
}

// Load OBJ into MPMesh (uses simple_kernel::Point_3)
bool LoadObjFile(const std::string& filename, Mesh& mesh, std::string& error) {
    ObjData data;
    if (!ParseObjFile(filename, data, error)) {
        return false;
    }

    // Create the builder with simple_kernel Point type
    ObjPolyhedronBuilderT<Mesh::HalfedgeDS, simple_kernel::Point_3> builder;
    builder.vertices = std::move(data.vertices);
    builder.faces = std::move(data.faces);

    // Build the polyhedron
    try {
        mesh.delegate(builder);
    } catch (const std::exception& e) {
        error = std::string("Failed to build polyhedron: ") + e.what();
        return false;
    }

    if (!mesh.is_valid()) {
        error = "Resulting mesh is not valid";
        return false;
    }

    if (!mesh.is_closed()) {
        std::cerr << "Warning: Mesh is not closed (has boundary edges)" << std::endl;
    }

    return true;
}

// Load OBJ into basic Polyhedron (uses K::Point_3 for mesh domain)
bool LoadObjFile(const std::string& filename, Polyhedron& mesh, std::string& error) {
    ObjData data;
    if (!ParseObjFile(filename, data, error)) {
        return false;
    }

    // Create the builder with K (Exact_predicates_inexact_constructions_kernel) Point type
    ObjPolyhedronBuilderT<Polyhedron::HalfedgeDS, K::Point_3> builder;
    builder.vertices = std::move(data.vertices);
    builder.faces = std::move(data.faces);

    // Build the polyhedron
    try {
        mesh.delegate(builder);
    } catch (const std::exception& e) {
        error = std::string("Failed to build polyhedron: ") + e.what();
        return false;
    }

    if (!mesh.is_valid()) {
        error = "Resulting mesh is not valid";
        return false;
    }

    if (!mesh.is_closed()) {
        std::cerr << "Warning: Mesh is not closed (has boundary edges)" << std::endl;
    }

    return true;
}
