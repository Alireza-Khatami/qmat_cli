/**
 * QMAT Command Line Interface
 *
 * A standalone CLI to run QMAT's core medial axis computation and simplification
 * without Qt GUI dependencies.
 *
 * Supported file formats:
 *   .off  - Object File Format
 *   .obj  - Wavefront OBJ format
 *
 * Usage:
 *   qmat_cli <input.off|input.obj> [options]
 *
 * Options:
 *   --simplify <N>     Simplify to N vertices (default: no simplification)
 *   --k <value>        K factor for slab initialization (default: 0.00001)
 *   --output <prefix>  Output file prefix (default: input filename without extension)
 *   --help             Show this help message
 *
 * Examples:
 *   qmat_cli model.off
 *   qmat_cli model.obj
 *   qmat_cli model.obj --simplify 500 --k 0.0001 --output result
 */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstring>
#include <ctime>
#include <chrono>
#include <thread>
#include <filesystem>
#include <unordered_map>
#include "ThreeDimensionalShape.h"
#include "ObjLoader.h"
#include "MatVisualizerCommon.h"   // ExportMatAsOff, ExportInitialVisualizeInfo
#include "QmatVisualizer.h"        // ExportSimpVisualizeInfo
#ifdef QMAT_WITH_VCGLIB
#  include "VcgDirectSimplifier.h"
#  include "VdeVisualizer.h"       // ExportSnapshot* (CLI-side, no polyscope dep)
#endif

// Bodies of ExportMatAsOff / ExportInitialVisualizeInfo /
// ExportSimpVisualizeInfo were moved out of main_cli.cpp:
//   - shared SlabMesh-side: MatVisualizerCommon.h/.cpp
//   - QMAT post-simp:       QmatVisualizer.h/.cpp
//   - VDE snapshot-side:    VdeVisualizer.h/.cpp (ExportSnapshot* family)
// Headers are already included at the top of this file.

// Simple command line argument parsing
struct CLIOptions {
    std::string inputFile;
    std::string outputPrefix;
    int simplifyTarget = -1;  // --simplify <N>  : vertex target (-1 = no simplification)
    int faceTarget = -1;      // --nf <N>        : face target (-1 = unset)

    // --simplifier picks which decimator runs.  All three read the same target
    // inputs above.  qmat/vcg-port both go through SlabMesh::Simplify and which
    // one runs is fixed at compile time by ONLY_USE_QEM_CONDITION_CHECKS — a
    // mismatched flag at runtime warns.  vcg-direct uses vcglib directly and
    // requires -DQMAT_WITH_VCGLIB=ON.
    enum class SimplifierMode { Qmat, VcgPort, VcgDirect };
#if defined(ONLY_USE_QEM_CONDITION_CHECKS)
    SimplifierMode simplifier = SimplifierMode::VcgPort;
#else
    SimplifierMode simplifier = SimplifierMode::Qmat;
#endif
    double k = 0.00001;
    double featureAngleDeg = 10.0; // turning-angle threshold for sharp feature protection
    // VCG-faithful QEM thresholds (ONLY_USE_QEM_CONDITION_CHECKS path).
    // Defaults mirror VcgQuadricParameter's struct defaults.
    double qemQualityThr = 0.3;     // --qem-quality-thr  → VcgQuadricParameter::QualityThr
    double qemBoundaryWeight = 1.0; // --qem-boundary-weight = MeshLab UI "Boundary
                                    // Preserving Weight" slider; effective weight is
                                    // 0.5 * this (meshfilter.cpp:14).  Default 1.0 → 0.5.
    // VCG-faithful QEM toggles (MeshLab QEM checkboxes), default to MeshLab's GUI
    // defaults; override with --qem-* <0|1> so no recompile is needed.
    bool qemOptimalPlacement = true; // --qem-optimal-placement (MeshLab "Optimal position")
    bool qemPreserveTopology = true; // --qem-preserve-topology (MeshLab "Preserve Topology")
    bool qemNormalCheck      = true; // --qem-normal-check       (MeshLab "Preserve Normal")
    bool visualize = false;
    bool showHelp = false;
    bool valid = true;
    std::string errorMessage;
#ifdef USE_MATSTRUCT_INITIALIZATION
    std::string matstructFile;  // path to typed .ma file from external tool
#endif
};

void printUsage(const char* programName) {
    std::cout << "QMAT Command Line Interface\n"
              << "Compute medial axis and optionally simplify.\n\n"
              << "Usage:\n"
              << "  " << programName << " <input.off|input.obj> [options]\n\n"
              << "Supported formats:\n"
              << "  .off               Object File Format\n"
              << "  .obj               Wavefront OBJ format\n\n"
              << "Options:\n"
              << "  --simplify <N>     Simplify to N vertices (default: no simplification)\n"
              << "  --nf <N>           Stop QEM when MAT face count reaches N (requires --simplify; default: off)\n"
              << "  --qem-quality-thr <v>     QEM QualityThr threshold (default: 0.3)\n"
              << "  --qem-boundary-weight <v> MeshLab UI 'Boundary Preserving Weight' slider;\n"
              << "                            effective weight = 0.5*v (default: 1.0 -> 0.5)\n"
              << "  --qem-optimal-placement <0|1>  Optimal vertex position (default: 1)\n"
              << "  --qem-preserve-topology <0|1>  Preserve topology / link condition (default: 1)\n"
              << "  --qem-normal-check <0|1>       Preserve normal (default: 1)\n"
              << "  --simplifier <mode>  qmat | vcg-port | vcg-direct\n"
              << "                       qmat/vcg-port: SlabMesh::Simplify (compile-time pick).\n"
              << "                       vcg-direct: vcglib TriEdgeCollapseQuadric (needs QMAT_WITH_VCGLIB).\n"
              << "  --k <value>        K factor for slab initialization (default: 0.00001)\n"
              << "  --feature-angle <deg>  Turning-angle threshold for sharp feature protection (default: 30)\n"
              << "  --output <prefix>  Output file prefix (default: input filename)\n"
#ifdef QMAT_WITH_POLYSCOPE
              << "  --visualize        Open live Polyscope viewer during simplification\n"
#endif
#ifdef USE_MATSTRUCT_INITIALIZATION
              << "  --matstruct <file> Path to typed .ma file from external MAT tool\n"
#endif
              << "  --help             Show this help message\n\n"
              << "Examples:\n"
              << "  " << programName << " model.off\n"
              << "  " << programName << " model.obj\n"
              << "  " << programName << " model.obj --simplify 500 --k 0.0001 --output result\n";
}

CLIOptions parseArguments(int argc, char* argv[]) {
    CLIOptions options;

    if (argc < 2) {
        options.valid = false;
        options.errorMessage = "No input file specified.";
        return options;
    }

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];

        if (arg == "--help" || arg == "-h") {
            options.showHelp = true;
            return options;
        }
        else if (arg == "--simplify") {
            if (i + 1 >= argc) {
                options.valid = false;
                options.errorMessage = "--simplify requires a value.";
                return options;
            }
            try {
                options.simplifyTarget = std::stoi(argv[++i]);
                if (options.simplifyTarget <= 0) {
                    options.valid = false;
                    options.errorMessage = "--simplify value must be positive.";
                    return options;
                }
            } catch (...) {
                options.valid = false;
                options.errorMessage = "Invalid value for --simplify.";
                return options;
            }
        }
        else if (arg == "--nf") {
            if (i + 1 >= argc) {
                options.valid = false;
                options.errorMessage = "--nf requires a value.";
                return options;
            }
            try {
                options.faceTarget = std::stoi(argv[++i]);
                if (options.faceTarget <= 0) {
                    options.valid = false;
                    options.errorMessage = "--nf value must be positive.";
                    return options;
                }
            } catch (...) {
                options.valid = false;
                options.errorMessage = "Invalid value for --nf.";
                return options;
            }
        }
        else if (arg == "--qem-quality-thr") {
            if (i + 1 >= argc) {
                options.valid = false;
                options.errorMessage = "--qem-quality-thr requires a value.";
                return options;
            }
            try {
                options.qemQualityThr = std::stod(argv[++i]);
                if (options.qemQualityThr <= 0) {
                    options.valid = false;
                    options.errorMessage = "--qem-quality-thr value must be positive.";
                    return options;
                }
            } catch (...) {
                options.valid = false;
                options.errorMessage = "Invalid value for --qem-quality-thr.";
                return options;
            }
        }
        else if (arg == "--qem-boundary-weight") {
            if (i + 1 >= argc) {
                options.valid = false;
                options.errorMessage = "--qem-boundary-weight requires a value.";
                return options;
            }
            try {
                options.qemBoundaryWeight = std::stod(argv[++i]);
                if (options.qemBoundaryWeight < 0) {
                    options.valid = false;
                    options.errorMessage = "--qem-boundary-weight value must be non-negative.";
                    return options;
                }
            } catch (...) {
                options.valid = false;
                options.errorMessage = "Invalid value for --qem-boundary-weight.";
                return options;
            }
        }
        else if (arg == "--simplifier") {
            if (i + 1 >= argc) {
                options.valid = false;
                options.errorMessage = "--simplifier requires a value (qmat|vcg-port|vcg-direct).";
                return options;
            }
            std::string v = argv[++i];
            if      (v == "qmat")       options.simplifier = CLIOptions::SimplifierMode::Qmat;
            else if (v == "vcg-port")   options.simplifier = CLIOptions::SimplifierMode::VcgPort;
            else if (v == "vcg-direct") options.simplifier = CLIOptions::SimplifierMode::VcgDirect;
            else {
                options.valid = false;
                options.errorMessage = "--simplifier must be qmat, vcg-port, or vcg-direct.";
                return options;
            }
        }
        else if (arg == "--qem-optimal-placement" || arg == "--qem-preserve-topology"
                 || arg == "--qem-normal-check") {
            if (i + 1 >= argc) {
                options.valid = false;
                options.errorMessage = arg + " requires a value (0|1).";
                return options;
            }
            std::string v = argv[++i];
            bool on;
            if (v == "1" || v == "true" || v == "on")        on = true;
            else if (v == "0" || v == "false" || v == "off") on = false;
            else {
                options.valid = false;
                options.errorMessage = arg + " must be 0 or 1.";
                return options;
            }
            if      (arg == "--qem-optimal-placement") options.qemOptimalPlacement = on;
            else if (arg == "--qem-preserve-topology") options.qemPreserveTopology = on;
            else                                       options.qemNormalCheck      = on;
        }
        else if (arg == "--k") {
            if (i + 1 >= argc) {
                options.valid = false;
                options.errorMessage = "--k requires a value.";
                return options;
            }
            try {
                options.k = std::stod(argv[++i]);
                if (options.k <= 0) {
                    options.valid = false;
                    options.errorMessage = "--k value must be positive.";
                    return options;
                }
            } catch (...) {
                options.valid = false;
                options.errorMessage = "Invalid value for --k.";
                return options;
            }
        }
        else if (arg == "--output") {
            if (i + 1 >= argc) {
                options.valid = false;
                options.errorMessage = "--output requires a value.";
                return options;
            }
            options.outputPrefix = argv[++i];
        }
#ifdef USE_MATSTRUCT_INITIALIZATION
        else if (arg == "--matstruct") {
            if (i + 1 >= argc) {
                options.valid = false;
                options.errorMessage = "--matstruct requires a file path.";
                return options;
            }
            options.matstructFile = argv[++i];
        }
#endif
        else if (arg == "--feature-angle") {
            if (i + 1 >= argc) {
                options.valid = false;
                options.errorMessage = "--feature-angle requires a value in degrees.";
                return options;
            }
            try { options.featureAngleDeg = std::stod(argv[++i]); std::cout << "Feature angle threshold set to " << options.featureAngleDeg << " degrees.\n"; }
            catch (...) {
                options.valid = false;
                options.errorMessage = "Invalid value for --feature-angle.";
                return options;
            }
        }
        else if (arg == "--visualize") {
#ifdef QMAT_WITH_POLYSCOPE
            options.visualize = true;
#else
            std::cerr << "Warning: --visualize ignored (built without QMAT_WITH_POLYSCOPE)\n";
#endif
        }
        else if (arg[0] == '-') {
            options.valid = false;
            options.errorMessage = "Unknown option: " + arg;
            return options;
        }
        else {
            // Positional argument - input file
            if (options.inputFile.empty()) {
                options.inputFile = arg;
            } else {
                options.valid = false;
                options.errorMessage = "Multiple input files specified.";
                return options;
            }
        }
    }

    if (options.inputFile.empty()) {
        options.valid = false;
        options.errorMessage = "No input file specified.";
        return options;
    }

    // --nf only refines an existing --simplify request; ignore it otherwise.
    if (options.faceTarget > 0 && options.simplifyTarget <= 0) {
        std::cerr << "Warning: --nf ignored (requires --simplify to also be set).\n";
        options.faceTarget = -1;
    }

    // Set default output prefix from input filename
    if (options.outputPrefix.empty()) {
        options.outputPrefix = options.inputFile;
        // Remove .off or .obj extension if present
        size_t dotPos = options.outputPrefix.rfind('.');
        if (dotPos != std::string::npos) {
            std::string ext = options.outputPrefix.substr(dotPos);
            if (ext == ".off" || ext == ".OFF" || ext == ".obj" || ext == ".OBJ") {
                options.outputPrefix = options.outputPrefix.substr(0, dotPos);
            }
        }
    }

    // ── redirect all output into a subfolder named after the mesh stem ────────
    // e.g.  path/to/cube_subdiv_500  →  create  path/to/cube_subdiv_500/
    //                                   prefix  path/to/cube_subdiv_500/cube_subdiv_500
    {
        namespace fs = std::filesystem;
        fs::path p(options.outputPrefix);
        std::string stem = p.filename().string();   // final component only
        fs::path outDir  = p.parent_path() / stem;  // sibling folder named after stem
        std::error_code ec;
        fs::create_directories(outDir, ec);
        if (ec)
            std::cerr << "Warning: could not create output directory "
                      << outDir << ": " << ec.message() << "\n";
        options.outputPrefix = (outDir / stem).string();
    }

    return options;
}

// RenderVcgDirectSnapshot moved to VdeVisualizer.cpp (Phase 2).

int main(int argc, char* argv[]) {

    std::cout << "argc = " << argc << "\n";
    // for (int i = 0; i < argc; ++i)
    //     std::cout << "argv[" << i << "] = [" << argv[i] << "]\n";
    // Parse command line arguments
    CLIOptions options = parseArguments(argc, argv);

    if (options.showHelp) {
        printUsage(argv[0]);
        return 0;
    }

    if (!options.valid) {
        std::cerr << "Error: " << options.errorMessage << std::endl;
        std::cerr << "Use --help for usage information." << std::endl;
        return 1;
    }

    std::cout << "QMAT CLI - Medial Axis Computation" << std::endl;
    std::cout << "===================================" << std::endl;
    std::cout << "Input file: " << options.inputFile << std::endl;
    std::cout << "Output prefix: " << options.outputPrefix << std::endl;
    std::cout << "K value: " << options.k << std::endl;
    if (options.simplifyTarget > 0) {
        std::cout << "Simplify target: " << options.simplifyTarget << " vertices" << std::endl;
        if (options.faceTarget > 0)
            std::cout << "Face target (--nf): stop QEM at " << options.faceTarget << " faces" << std::endl;
        std::cout << "QEM QualityThr: " << options.qemQualityThr
                  << ", QEM BoundaryWeight: " << options.qemBoundaryWeight << std::endl;
    }
    std::cout << std::endl;

    // Create ThreeDimensionalShape object
    ThreeDimensionalShape shape;

    // Step 1: Load the mesh file (OFF or OBJ)
    std::cout << "Loading mesh from " << options.inputFile << "..." << std::endl;
    long startTime = clock();

    bool loadSuccess = false;
    if (IsObjFile(options.inputFile)) {
        // Load OBJ file using tinyobjloader
        std::string objError;
        loadSuccess = LoadObjFile(options.inputFile, shape.input, objError);
        if (!loadSuccess) {
            std::cerr << "Error loading OBJ file: " << objError << std::endl;
            return 1;
        }
    } else if (IsOffFile(options.inputFile)) {
        // Load OFF file using CGAL
        std::ifstream stream(options.inputFile.c_str());
        if (!stream) {
            std::cerr << "Error: Could not open file " << options.inputFile << std::endl;
            return 1;
        }
        stream >> shape.input;
        stream.close();
        loadSuccess = true;
    } else {
        std::cerr << "Error: Unsupported file format. Use .off or .obj files." << std::endl;
        return 1;
    }

    // Compute mesh properties
    shape.input.computebb();
    shape.input.GenerateList();
    shape.input.GenerateRandomColor();
    shape.input.compute_normals();

    long loadTime = clock() - startTime;
    std::cout << "  Loaded mesh with " << shape.input.size_of_vertices() << " vertices, "
              << shape.input.size_of_facets() << " faces" << std::endl;
    std::cout << "  Load time: " << loadTime << " ms" << std::endl;

    // Step 2: Create CGAL mesh domain for inside/outside queries
    std::cout << "Creating mesh domain..." << std::endl;
    Polyhedron pol;
    if (IsObjFile(options.inputFile)) {
        // Load OBJ file for mesh domain
        std::string objError;
        if (!LoadObjFile(options.inputFile, pol, objError)) {
            std::cerr << "Error loading OBJ file for mesh domain: " << objError << std::endl;
            return 1;
        }
    } else {
        // Load OFF file for mesh domain
        std::ifstream streamPol(options.inputFile.c_str());
        streamPol >> pol;
        streamPol.close();
    }

    Mesh_domain* domain = new Mesh_domain(pol);
    shape.input.domain = domain;
    shape.input_nmm.domain = domain;
    shape.input_nmm.pmesh = &shape.input;
    shape.input_nmm.meshname = options.outputPrefix;

#ifndef USE_MATSTRUCT_INITIALIZATION
    // Step 3: Compute Delaunay Triangulation and Medial Axis
    // If the .ma file already exists on disk, skip the expensive DT + MAT
    // computation and reuse it directly.
    const std::string maFile = options.outputPrefix + ".ma";
    if (std::filesystem::exists(maFile)) {
        std::cout << "Found existing MA file: " << maFile << std::endl;
        std::cout << "  Skipping DT and MAT computation." << std::endl;
    } else {
        std::cout << "Computing Delaunay Triangulation..." << std::endl;
        startTime = clock();
        shape.input.computedt();
        long dtTime = clock() - startTime;
        std::cout << "  DT computation time: " << dtTime << " ms" << std::endl;

        // Cluster boundary sample points (input mesh vertices) by position + normal.
        // Output: <outputPrefix>_boundary_clusters.txt




        std::cout << "Computing Medial Axis..." << std::endl;
        startTime = clock();
        shape.ComputeInputNMM();
        long maTime = clock() - startTime;
        std::cout << "  MA computation time: " << maTime << " ms" << std::endl;
        std::cout << "  Raw MA exported to: " << maFile << std::endl;
    }
#endif // USE_MATSTRUCT_INITIALIZATION

#ifdef USE_MATSTRUCT_INITIALIZATION
    // Step 3b: Also compute the Voronoi-based initial MAT for comparison.
    // This is the raw MAT derived from the Delaunay triangulation of the input
    // surface samples — the same one produced without USE_MATSTRUCT.
    // All output goes into a dedicated init_vor_mat/ subfolder:
    //   init_vor_mat/<stem>.ma                        — raw .ma format
    //   init_vor_mat/<stem>.ply                       — ASCII PLY (vertex/face/edge)
    //   init_vor_mat/<stem>_sampledpoints.txt
    //   init_vor_mat/<stem>_surf_2_vor_mat_map.txt
    //   init_vor_mat/<stem>_mat_vertices_with_fields.txt
    {
        namespace fs = std::filesystem;
        fs::path pp(options.outputPrefix);

        // Create the dedicated subfolder.
        fs::path vorDir = pp.parent_path() / "init_vor_mat";
        fs::create_directories(vorDir);
        const std::string vorPrefix = (vorDir / pp.filename()).string();

        std::cout << std::endl << "Computing Voronoi MAT for comparison..." << std::endl;
        shape.input_nmm.meshname = vorPrefix;   // Export() and sidecar writers use this

        startTime = clock();
        shape.input.computedt();
        shape.ComputeInputNMM();                // writes .ma and all sidecar files
        long vorTime = clock() - startTime;
        std::cout << "  Voronoi MAT computation time: " << vorTime << " ms" << std::endl;
        std::cout << "  Voronoi MAT written to: " << vorPrefix << ".ma" << std::endl;

        // ── Export as ASCII PLY for MeshLab visualization ─────────────────────
        // PLY natively supports vertex, face, and edge as first-class elements.
        //
        // Bug fix: faces in the MAT are stored as std::set<unsigned> of vertex
        // indices. In degenerate cases two of the three circumcenter tags can be
        // equal, leaving the set with fewer than 3 entries. Writing the actual
        // set size (not a hardcoded 3) keeps the face list correctly aligned so
        // MeshLab does not lose sync and drop subsequent faces/edges.
        const std::string vorPly = vorPrefix + ".ply";
        {
            std::ofstream ply(vorPly);
            if (ply.is_open()) {
                const auto& nmverts = shape.input_nmm.vertices;
                const auto& nmfaces = shape.input_nmm.faces;
                const auto& nmedges = shape.input_nmm.edges;

                // ── Header ────────────────────────────────────────────────────
                ply << "ply\n"
                    << "format ascii 1.0\n"
                    << "element vertex " << nmverts.size() << "\n"
                    << "property double x\n"
                    << "property double y\n"
                    << "property double z\n";

                if (!nmfaces.empty())
                    ply << "element face " << nmfaces.size() << "\n"
                        << "property list uchar int vertex_indices\n";

                if (!nmedges.empty())
                    ply << "element edge " << nmedges.size() << "\n"
                        << "property int vertex1\n"
                        << "property int vertex2\n";

                ply << "end_header\n";

                // ── Vertices ──────────────────────────────────────────────────
                ply << std::fixed << std::setprecision(15);
                for (const auto& bvp : nmverts) {
                    const auto& c = bvp.second->sphere.center;
                    ply << c[0] << " " << c[1] << " " << c[2] << "\n";
                }

                // ── Faces ─────────────────────────────────────────────────────
                if (!nmfaces.empty()) {
                    for (const auto& bfp : nmfaces) {
                        ply << bfp.second->vertices_.size();
                        for (unsigned vi : bfp.second->vertices_)
                            ply << " " << vi;
                        ply << "\n";
                    }
                }

                // ── Edges ─────────────────────────────────────────────────────
                if (!nmedges.empty()) {
                    for (const auto& bep : nmedges)
                        ply << bep.second->vertices_.first
                            << " " << bep.second->vertices_.second << "\n";
                }

                std::cout << "  Voronoi MAT PLY written to: " << vorPly
                          << "  (v=" << nmverts.size()
                          << " f=" << nmfaces.size()
                          << " e=" << nmedges.size() << ")" << std::endl;
            } else {
                std::cerr << "  Warning: could not open " << vorPly << " for writing\n";
            }
        }

        // Restore meshname so the rest of the pipeline uses the normal prefix.
        shape.input_nmm.meshname = options.outputPrefix;
    }
#endif // USE_MATSTRUCT_INITIALIZATION


    // Step 4: If simplification requested, load into slab mesh and simplify
    if (options.simplifyTarget > 0) {
        std::cout << std::endl << "Loading MA for simplification..." << std::endl;

        // Setup slab mesh
        shape.slab_mesh.pmesh = &shape.input;
        shape.slab_mesh.type = 1;
        shape.slab_mesh.k = options.k;
        shape.slab_mesh.bound_weight = 1.0;
        shape.slab_mesh.export_prefix = options.outputPrefix;

        // VCG-faithful QEM thresholds + optional face-count stop (--nf).
        // Consumed only by the ONLY_USE_QEM_CONDITION_CHECKS path in Simplify().
        shape.slab_mesh.qem_quality_thr     = options.qemQualityThr;
        shape.slab_mesh.qem_boundary_weight = options.qemBoundaryWeight;
        shape.slab_mesh.qem_face_target     = options.faceTarget;
        shape.slab_mesh.qem_optimal_placement = options.qemOptimalPlacement;
        shape.slab_mesh.qem_preserve_topology = options.qemPreserveTopology;
        shape.slab_mesh.qem_normal_check      = options.qemNormalCheck;

        // Initialize slab mesh settings (same as GUI initialize())
        shape.slab_mesh.preserve_boundary_method = 0;
        shape.slab_mesh.hyperbolic_weight_type = 3;
        shape.slab_mesh.compute_hausdorff = false;
        shape.slab_mesh.boundary_compute_scale = 0;
        shape.slab_mesh.prevent_inversion = true;

#ifdef USE_MATSTRUCT_INITIALIZATION
        // Load the typed .ma file from the external MAT tool.
        shape.LoadMatstructMA(options.matstructFile);
        // struct_ids sets are populated during LoadMatstructMA; collapsibility is
        // checked inline in CanMerge — no separate pre-computation needed.
        shape.ComputeFeatureEdges();          // detect sharp/concave edges on input mesh
        shape.ExportSurfacemeshFeatureEdges(options.outputPrefix);
        // ClusterNMNBplist() is skipped — T-types are already loaded from file.
#else
        // Load the MA file into the slab mesh (computed above or loaded from cache)
        shape.LoadInputNMM(maFile);
        shape.ComputeFeatureEdges();          // detect sharp/concave edges on input mesh
        shape.ExportSurfacemeshFeatureEdges(options.outputPrefix);
        shape.slab_mesh.ClusterNMNBplist();   // must run before DetermineTopology (T1 filtering)
#endif
        // shape.slab_mesh.DetermineTopology();  // uses T-types to correctly classify topology
        shape.slab_mesh.feature_angle_threshold = options.featureAngleDeg;
        shape.slab_mesh.MarkSharpFeatureVertices(options.featureAngleDeg);

        std::cout << "  Loaded slab mesh with " << shape.slab_mesh.numVertices << " vertices" << std::endl;

        // Initialize slab mesh for simplification
        std::cout << "Initializing slab mesh..." << std::endl;
        startTime = clock();
        long initTime = shape.LoadSlabMesh();
        std::cout << "  Initialization time: " << initTime << " ms" << std::endl;

        // Simplify
        int currentVertices = shape.slab_mesh.numVertices;

        // ── vcg-direct branch: vcg's own TriEdgeCollapseQuadric via vcglib ──
        // Bypasses SlabMesh::Simplify entirely.  Takes --simplify (vertex) and/or
        // --nf (face) as targets; --nf wins if both are set.  Output is shown in
        // polyscope as "VCG Direct Simplified MAT"; QMAT post-export steps are
        // skipped because the slab mesh itself is left untouched.
        if (options.simplifier == CLIOptions::SimplifierMode::VcgDirect) {
#ifdef QMAT_WITH_VCGLIB
            if (options.simplifyTarget <= 0 && options.faceTarget <= 0) {
                std::cout << "Warning: --simplifier vcg-direct requires --simplify or --nf.\n";
            } else {
                startTime = clock();
                shape.slab_mesh.CleanIsolatedVertices();
                ExportMatAsOff(shape.slab_mesh, options.outputPrefix + "_mat_initial.off");
                ExportInitialVisualizeInfo(shape.slab_mesh,
                                           options.outputPrefix + "_initial_visualize_info.json");

                VcgDirectParams p;
                p.TargetFaceNum         = options.faceTarget;
                p.TargetVertexNum       = options.simplifyTarget;
                p.QualityThr            = options.qemQualityThr;
                p.BoundaryQuadricWeight = options.qemBoundaryWeight;
                p.OptimalPlacement      = options.qemOptimalPlacement;
                p.PreserveTopology      = options.qemPreserveTopology;
                p.NormalCheck           = options.qemNormalCheck;

                LiveUpdateCallback live_cb;
#ifdef QMAT_WITH_POLYSCOPE
                VdeVisualizer vde;
                ViewerState& vs = vde.State();
                if (options.visualize) {
                    vs.outputPrefix = options.outputPrefix;
                    vde.Setup(shape.slab_mesh);
                    live_cb = vde.MakeLiveCallback();
                }
#endif

                VcgDirectResult res;
                bool ok = RunVcgDirectSimplify(shape.slab_mesh, p, res, live_cb);
                std::cout << "  vcg-direct time: " << (clock() - startTime) << " ms\n";

                // Snapshot-driven post-simp exports — the SlabMesh is unchanged
                // by vcg-direct, so we cannot reuse QMAT's SlabMesh-based exporters.
                if (ok) {
                    const VcgDirectSnapshot& snap = res.snapshot;
                    ExportSnapshotAsOff(snap, options.outputPrefix + "_mat_simplified.off");
                    ExportSnapshotVisualizeInfo(snap,
                        options.outputPrefix + "_simp_visualize_info.json");
                    ExportSnapshotRejectionSkeleton(snap,
                        options.outputPrefix + "_rejection_skeleton.ply");
                    ExportSnapshotMatTyped(snap, options.outputPrefix + ".mat_typed");
                    ExportSnapshotClusterPLY(snap, options.outputPrefix + "_cluster.ply");
                    ExportSnapshotMa(snap, options.outputPrefix);
                }

#ifdef QMAT_WITH_POLYSCOPE
                if (ok && options.visualize) {
                    vde.Render(res.snapshot);
                    std::cout << "vcg-direct done. Close the viewer to exit.\n";
                    vde.Show();
                }
#endif
            }
#else
            std::cerr << "--simplifier vcg-direct requires -DQMAT_WITH_VCGLIB=ON.\n";
#endif
        }
        else if (options.simplifyTarget >= currentVertices) {
            std::cout << "Warning: Target vertex count (" << options.simplifyTarget
                      << ") >= current count (" << currentVertices << "). Skipping simplification." << std::endl;
        } else {
            int reductionCount = currentVertices - options.simplifyTarget;
            std::cout << "Simplifying from " << currentVertices << " to " << options.simplifyTarget
                      << " vertices (removing " << reductionCount << ")..." << std::endl;

            // Label every vertex with its topological class before simplification.
            // Labels (mirrors check_non_manifold in read_qmat_output.py):
            //   REGULAR    – interior vertex (not on any boundary or non-manifold edge)
            //   BOUNDARY   – on a boundary edge (edge shared by exactly 1 face)
            //   NM_EDGE    – on a non-manifold edge (>2 faces), but not a corner
            //   NM_CORNER  – on both a boundary and a non-manifold edge
            // Simplify will only collapse an edge whose two endpoints share the same label.

            startTime = clock();
            shape.slab_mesh.CleanIsolatedVertices();

            // Export MAT before simplification starts
            ExportMatAsOff(shape.slab_mesh, options.outputPrefix + "_mat_initial.off");

            // JSON snapshot of the initial MAT: vertex positions, edge/face
            // connectivity, and per-primitive struct_ids.  Companion to the
            // post-simplification _simp_visualize_info.json so external tools
            // can compare structure assignments before vs. after simplifying.
            ExportInitialVisualizeInfo(shape.slab_mesh,
                                       options.outputPrefix + "_initial_visualize_info.json");

#ifdef QMAT_WITH_POLYSCOPE
            QmatVisualizer qmat;
            ViewerState& vs = qmat.State();
            vs.outputPrefix = options.outputPrefix;
            if (options.visualize)
                qmat.Setup(shape.slab_mesh);
#endif
            // qmat/vcg-port: which path runs is fixed at compile time by
            // ONLY_USE_QEM_CONDITION_CHECKS — warn if --simplifier disagrees.
#if defined(ONLY_USE_QEM_CONDITION_CHECKS)
            if (options.simplifier == CLIOptions::SimplifierMode::Qmat)
                std::cerr << "[warn] --simplifier qmat ignored: binary built with ONLY_USE_QEM_CONDITION_CHECKS → vcg-port runs.\n";
#else
            if (options.simplifier == CLIOptions::SimplifierMode::VcgPort)
                std::cerr << "[warn] --simplifier vcg-port ignored: binary built without ONLY_USE_QEM_CONDITION_CHECKS → qmat runs.\n";
#endif
            shape.slab_mesh.allow_steep_collapse = true;
            shape.slab_mesh.Simplify(reductionCount);
            long simplifyTime = clock() - startTime;

            std::cout << "  Simplification time: " << simplifyTime << " ms" << std::endl;
            std::cout << "  Final vertex count: " << shape.slab_mesh.numVertices << std::endl;

#ifdef QMAT_WITH_POLYSCOPE
            if (options.visualize) {
                qmat.RenderFinal(shape.slab_mesh);
                std::cout << "Simplification done. Close the viewer window to exit.\n";
                qmat.Show();
            }
#endif


            // Must run before Export(): Export() calls AdjustStorage() which
            // remaps all edge IDs, invalidating edge_last_rejection.
            // radius is a fraction of the geometry's bbox diagonal (computed inside
            // ExportSkeletonPLY); the old 1.0f was an absolute value that, against the
            // normalised ≈[-1,1] coordinates, made every cylinder model-sized.
            shape.slab_mesh.ExportSkeletonPLY(options.outputPrefix + "_rejection_skeleton.ply");
            std::cout << "  Skeleton PLY exported: " << options.outputPrefix << "_rejection_skeleton.ply\n";

            // JSON with per-primitive struct_ids / cluster-type / topo-type /
            // rejection-reason plus id->color legends, for external
            // re-visualization tools.  Must run before Export() (same
            // AdjustStorage dependency as the skeleton PLY above).
            ExportSimpVisualizeInfo(shape.slab_mesh,
                                    options.outputPrefix + "_simp_visualize_info.json");

            // .mat_typed: same as .ma but each vertex line has its ClusterType
            // name appended, e.g. "v x y z r MS_Seam_Boundary".
            shape.slab_mesh.ExportTypedMA(options.outputPrefix + ".mat_typed");
            std::cout << "  Typed MA exported: " << options.outputPrefix << ".mat_typed\n";

            // PLY with per-vertex colours matching the cluster-type palette.
            shape.slab_mesh.ExportClusterPLY(options.outputPrefix + "_cluster.ply");
            std::cout << "  Cluster PLY exported: " << options.outputPrefix << "_cluster.ply\n";

            // Compute final mesh properties
            shape.slab_mesh.ComputeFacesNormal();
            shape.slab_mesh.ComputeVerticesNormal();
            shape.slab_mesh.ComputeEdgesCone();
            shape.slab_mesh.ComputeFacesSimpleTriangles();

            // Export simplified mesh
            std::cout << "Exporting simplified MA..." << std::endl;
            shape.slab_mesh.Export(options.outputPrefix);
            std::cout << "  Simplified MA exported with prefix: " << options.outputPrefix << std::endl;
            ExportMatAsOff(shape.slab_mesh, options.outputPrefix + "_mat_simplified.off");


        }
    }

    std::cout << std::endl << "Done!" << std::endl;

    return 0;
}
