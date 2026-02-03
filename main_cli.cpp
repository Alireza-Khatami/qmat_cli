/**
 * QMAT Command Line Interface
 *
 * A standalone CLI to run QMAT's core medial axis computation and simplification
 * without Qt GUI dependencies.
 *
 * Usage:
 *   qmat_cli <input.off> [options]
 *
 * Options:
 *   --simplify <N>     Simplify to N vertices (default: no simplification)
 *   --k <value>        K factor for slab initialization (default: 0.00001)
 *   --output <prefix>  Output file prefix (default: input filename without extension)
 *   --help             Show this help message
 *
 * Examples:
 *   qmat_cli model.off
 *   qmat_cli model.off --simplify 1000
 *   qmat_cli model.off --simplify 500 --k 0.0001 --output result
 */

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <ctime>
#include "tiny_obj_loader.h"
#include "ThreeDimensionalShape.h"

// Simple command line argument parsing
struct CLIOptions {
    std::string inputFile;
    std::string outputPrefix;
    int simplifyTarget = -1;  // -1 means no simplification
    double k = 0.00001;
    bool showHelp = false;
    bool valid = true;
    std::string errorMessage;
};

void printUsage(const char* programName) {
    std::cout << "QMAT Command Line Interface\n"
              << "Compute medial axis and optionally simplify.\n\n"
              << "Usage:\n"
              << "  " << programName << " <input.off> [options]\n\n"
              << "Options:\n"
              << "  --simplify <N>     Simplify to N vertices (default: no simplification)\n"
              << "  --k <value>        K factor for slab initialization (default: 0.00001)\n"
              << "  --output <prefix>  Output file prefix (default: input filename)\n"
              << "  --help             Show this help message\n\n"
              << "Examples:\n"
              << "  " << programName << " model.off\n"
              << "  " << programName << " model.off --simplify 1000\n"
              << "  " << programName << " model.off --simplify 500 --k 0.0001 --output result\n";
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

    // Set default output prefix from input filename
    if (options.outputPrefix.empty()) {
        options.outputPrefix = options.inputFile;
        // Remove .off extension if present
        size_t dotPos = options.outputPrefix.rfind(".off");
        if (dotPos != std::string::npos && dotPos == options.outputPrefix.length() - 4) {
            options.outputPrefix = options.outputPrefix.substr(0, dotPos);
        }
    }

    return options;
}

int main(int argc, char* argv[]) {

    std::cout << "argc = " << argc << "\n";
    for (int i = 0; i < argc; ++i)
        std::cout << "argv[" << i << "] = [" << argv[i] << "]\n";
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
    }
    std::cout << std::endl;

    // Create ThreeDimensionalShape object
    ThreeDimensionalShape shape;

    // Step 1: Load the OFF file
    std::cout << "Loading mesh from " << options.inputFile << "..." << std::endl;
    long startTime = clock();

    std::ifstream stream(options.inputFile.c_str());
    if (!stream) {
        std::cerr << "Error: Could not open file " << options.inputFile << std::endl;
        return 1;
    }

    // Read the mesh
    stream >> shape.input;
    stream.close();

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
    std::ifstream streamPol(options.inputFile.c_str());
    Polyhedron pol;
    streamPol >> pol;
    streamPol.close();

    Mesh_domain* domain = new Mesh_domain(pol);
    shape.input.domain = domain;
    shape.input_nmm.domain = domain;
    shape.input_nmm.pmesh = &shape.input;
    shape.input_nmm.meshname = options.outputPrefix;

    // Step 3: Compute Delaunay Triangulation and Medial Axis
    std::cout << "Computing Delaunay Triangulation..." << std::endl;
    startTime = clock();
    shape.input.computedt();
    long dtTime = clock() - startTime;
    std::cout << "  DT computation time: " << dtTime << " ms" << std::endl;

    std::cout << "Computing Medial Axis..." << std::endl;
    startTime = clock();
    shape.ComputeInputNMM();
    long maTime = clock() - startTime;
    std::cout << "  MA computation time: " << maTime << " ms" << std::endl;
    std::cout << "  Raw MA exported to: " << options.outputPrefix << ".ma" << std::endl;

    // Step 4: If simplification requested, load into slab mesh and simplify
    if (options.simplifyTarget > 0) {
        std::cout << std::endl << "Loading MA for simplification..." << std::endl;

        // Setup slab mesh
        shape.slab_mesh.pmesh = &shape.input;
        shape.slab_mesh.type = 1;
        shape.slab_mesh.k = options.k;
        shape.slab_mesh.bound_weight = 1.0;

        // Initialize slab mesh settings (same as GUI initialize())
        shape.slab_mesh.preserve_boundary_method = 0;
        shape.slab_mesh.hyperbolic_weight_type = 3;
        shape.slab_mesh.compute_hausdorff = false;
        shape.slab_mesh.boundary_compute_scale = 0;
        shape.slab_mesh.prevent_inversion = false;

        // Load the MA file we just exported into the slab mesh
        std::string maFile =  options.outputPrefix + ".ma";
        shape.LoadInputNMM(maFile);

        std::cout << "  Loaded slab mesh with " << shape.slab_mesh.numVertices << " vertices" << std::endl;

        // Initialize slab mesh for simplification
        std::cout << "Initializing slab mesh..." << std::endl;
        startTime = clock();
        long initTime = shape.LoadSlabMesh();
        std::cout << "  Initialization time: " << initTime << " ms" << std::endl;

        // Simplify
        int currentVertices = shape.slab_mesh.numVertices;
        if (options.simplifyTarget >= currentVertices) {
            std::cout << "Warning: Target vertex count (" << options.simplifyTarget
                      << ") >= current count (" << currentVertices << "). Skipping simplification." << std::endl;
        } else {
            int reductionCount = currentVertices - options.simplifyTarget;
            std::cout << "Simplifying from " << currentVertices << " to " << options.simplifyTarget
                      << " vertices (removing " << reductionCount << ")..." << std::endl;

            startTime = clock();
            shape.slab_mesh.CleanIsolatedVertices();
            shape.slab_mesh.Simplify(reductionCount);
            long simplifyTime = clock() - startTime;

            std::cout << "  Simplification time: " << simplifyTime << " ms" << std::endl;
            std::cout << "  Final vertex count: " << shape.slab_mesh.numVertices << std::endl;

            // Compute final mesh properties
            shape.slab_mesh.ComputeFacesNormal();
            shape.slab_mesh.ComputeVerticesNormal();
            shape.slab_mesh.ComputeEdgesCone();
            shape.slab_mesh.ComputeFacesSimpleTriangles();

            // Export simplified mesh
            std::cout << "Exporting simplified MA..." << std::endl;
            shape.slab_mesh.Export(options.outputPrefix);
            std::cout << "  Simplified MA exported with prefix: " << options.outputPrefix << std::endl;
        }
    }

    std::cout << std::endl << "Done!" << std::endl;

    return 0;
}
