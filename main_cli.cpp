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
#include <string>
#include <cstring>
#include <ctime>
#include <chrono>
#include <thread>
#include <filesystem>
#include "ThreeDimensionalShape.h"
#include "ObjLoader.h"

#ifdef QMAT_WITH_POLYSCOPE
#  include "polyscope/polyscope.h"
#  include "polyscope/surface_mesh.h"
#  include "polyscope/curve_network.h"
#  include "polyscope/point_cloud.h"
#  include "polyscope/pick.h"
#  include "imgui.h"

// ── Polyscope live-simplification viewer ─────────────────────────────────────

// Builds contiguous arrays from the currently active MAT in `sm`.
// Also fills idx_to_vid: contiguous point-cloud index → original vertex id.
struct MatArrays {
    std::vector<std::array<double,3>>  verts;
    std::vector<std::array<size_t,2>>  edges;
    std::vector<std::array<size_t,3>>  faces;
    std::vector<unsigned>              idx_to_vid;   // index in verts → slab vertex id
    std::vector<std::array<float,3>>   vert_colors;  // per-vertex color (black = steep)
};

static MatArrays BuildMatArrays(const SlabMesh& sm)
{
    MatArrays out;
    std::map<unsigned,size_t> vid_map;
    for (unsigned i = 0; i < sm.vertices.size(); ++i) {
        if (!sm.vertices[i].first) continue;
        vid_map[i] = out.verts.size();
        const auto& c = sm.vertices[i].second->sphere.center;
        out.verts.push_back({c.X(), c.Y(), c.Z()});
        out.idx_to_vid.push_back(i);
        if (sm.vertices[i].second->is_steep_tetrahedron)
            out.vert_colors.push_back({0.0f, 0.0f, 0.0f});
        else
            out.vert_colors.push_back({1.0f, 1.0f, 0.4f});
    }
    for (unsigned i = 0; i < sm.edges.size(); ++i) {
        if (!sm.edges[i].first) continue;
        size_t a = vid_map.at(sm.edges[i].second->vertices_.first);
        size_t b = vid_map.at(sm.edges[i].second->vertices_.second);
        out.edges.push_back({a, b});
    }
    for (unsigned i = 0; i < sm.faces.size(); ++i) {
        if (!sm.faces[i].first) continue;
        auto it = sm.faces[i].second->vertices_.begin();
        size_t a = vid_map.at(*it++);
        size_t b = vid_map.at(*it++);
        size_t c = vid_map.at(*it);
        out.faces.push_back({a, b, c});
    }
    return out;
}

// Get the 3-D positions of all nmn_bplist entries for a slab vertex.
// Points come from the input mesh (sm.pmesh->pVertexList).
static std::vector<std::array<double,3>> BplistPositions(
    const SlabMesh& sm, unsigned vid)
{
    std::vector<std::array<double,3>> pts;
    if (!sm.vertices[vid].first) return pts;
    const auto& bps = sm.vertices[vid].second->nmn_bplist;
    pts.reserve(bps.size());
    // MAT vertices are stored normalized by bb_diagonal_length (see LoadInputNMM).
    // Divide bp positions by the same factor so they land in the same space.
    const double inv_diag = 1.0 / sm.pmesh->bb_diagonal_length;
    for (unsigned bp : bps) {
        const auto& p = sm.pmesh->pVertexList[bp]->point();
        pts.push_back({p[0] * inv_diag, p[1] * inv_diag, p[2] * inv_diag});
    }
    return pts;
}

// Registers the input surface mesh (pmesh) into Polyscope.
// Uses the CGAL Polyhedron face/vertex iterators directly so it works
// regardless of whether GenerateFaceList() was called.
static void RegisterInputMesh(const SlabMesh& sm)
{
    namespace ps = polyscope;
    if (!sm.pmesh) return;

    // Vertices — normalize by bb_diagonal_length to match the MAT coordinate space.
    const double inv_diag = 1.0 / sm.pmesh->bb_diagonal_length;
    std::vector<std::array<double,3>> verts;
    verts.reserve(sm.pmesh->pVertexList.size());
    for (unsigned i = 0; i < sm.pmesh->pVertexList.size(); ++i) {
        const auto& p = sm.pmesh->pVertexList[i]->point();
        verts.push_back({p[0] * inv_diag, p[1] * inv_diag, p[2] * inv_diag});
    }

    // Faces — iterate directly via CGAL halfedge structure.
    // Each facet_begin() halfedge circulates the face vertices in order.
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

    auto* mm = ps::registerSurfaceMesh("Input Mesh", verts, faces);
    mm->setSurfaceColor(glm::vec3(0.55f, 0.70f, 0.85f));
    mm->setTransparency(0.55f);
    mm->setEnabled(true);
}

struct ViewerState {
    int  collapse_count = 0;
    bool paused         = false;
    bool step_once      = false;
    int  update_every   = 1;
    std::chrono::steady_clock::time_point last_frame =
        std::chrono::steady_clock::now();

    // Maps the contiguous "MAT Verts" point cloud index → original slab vertex id.
    // Rebuilt every time BuildMatArrays is called.
    std::vector<unsigned> idx_to_vid;

    // Currently selected MAT vertex (-1 = none).
    int selected_vid = -1;
};

// Re-registers (or updates) the live MAT structures in Polyscope.
// Also updates vs.idx_to_vid from arr.
static void UpdateMatStructures(const MatArrays& arr, ViewerState& vs)
{
    namespace ps = polyscope;

    vs.idx_to_vid = arr.idx_to_vid;

    if (!arr.faces.empty()) {
        auto* mm = ps::registerSurfaceMesh("MAT Faces", arr.verts, arr.faces);
        mm->setSurfaceColor(glm::vec3(0.9f, 0.6f, 0.2f));
        mm->setTransparency(0.35f);
    } else if (ps::hasSurfaceMesh("MAT Faces")) {
        ps::removeStructure("MAT Faces");
    }

    if (!arr.edges.empty()) {
        auto* cn = ps::registerCurveNetwork("MAT Edges", arr.verts, arr.edges);
        cn->setColor(glm::vec3(1.0f, 0.80f, 0.30f));
        cn->setRadius(0.0008f, true);
    }

    if (!arr.verts.empty()) {
        auto* pc = ps::registerPointCloud("MAT Verts", arr.verts);
        pc->setPointRadius(0.0015, true);
        pc->setEnabled(false);
        if (!arr.vert_colors.empty())
            pc->addColorQuantity("type", arr.vert_colors)->setEnabled(true);
        else
            pc->setPointColor(glm::vec3(1.0f, 1.0f, 0.4f));
    }
}

// Call once before Simplify(). Initialises Polyscope, registers the initial
// MAT, installs the ImGui callback, and wires up sm.on_collapse_cb.
static void SetupSimplificationViewer(SlabMesh& sm, ViewerState& vs)
{
    namespace ps = polyscope;

    ps::init();
    ps::options::programName = "QMAT Simplification Viewer";
    ps::view::bgColor = {0.10f, 0.10f, 0.14f, 1.0f};

    // Register the input surface mesh so bplist points can be verified
    // visually against it.
    RegisterInputMesh(sm);

    // Register initial MAT (faded, hidden by default — shown via layer panel)
    MatArrays init = BuildMatArrays(sm);
    if (!init.faces.empty()) {
        auto* mm = ps::registerSurfaceMesh("Initial MAT Faces", init.verts, init.faces);
        mm->setSurfaceColor(glm::vec3(0.55f, 0.55f, 0.55f));
        mm->setTransparency(0.70f);
        mm->setEnabled(false);
    }
    if (!init.edges.empty()) {
        auto* cn = ps::registerCurveNetwork("Initial MAT Edges", init.verts, init.edges);
        cn->setColor(glm::vec3(0.6f, 0.6f, 0.6f));
        cn->setRadius(0.0005f, true);
        cn->setEnabled(false);
    }

    // Live MAT — updated every N collapses
    UpdateMatStructures(init, vs);

    // Placeholder: highlighted edge (v1 → v2), shown/updated per collapse
    {
        std::vector<std::array<double,3>>  p = {{0,0,0},{0,0,0}};
        std::vector<std::array<size_t,2>>  e = {{0,1}};
        auto* ce = ps::registerCurveNetwork("Collapsed Edge", p, e);
        ce->setColor(glm::vec3(1.0f, 0.15f, 0.15f));
        ce->setRadius(0.0030f, true);
        ce->setEnabled(false);
    }
    // Placeholder point clouds: v1 (red), v2 (orange), result (green)
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
    // Placeholder bplist clouds (per-collapse: v1 bps blue, v2 bps orange)
    {
        std::vector<std::array<double,3>> p = {{0,0,0}};

        auto* bp1 = ps::registerPointCloud("BPList v1", p);
        bp1->setPointColor(glm::vec3(0.30f, 0.50f, 1.0f));
        bp1->setPointRadius(0.0020, true);
        bp1->setEnabled(false);

        auto* bp2 = ps::registerPointCloud("BPList v2", p);
        bp2->setPointColor(glm::vec3(1.0f, 0.55f, 0.10f));
        bp2->setPointRadius(0.0020, true);
        bp2->setEnabled(false);

        // For manual pick selection
        auto* bpSel = ps::registerPointCloud("BPList selected", p);
        bpSel->setPointColor(glm::vec3(0.0f, 1.0f, 0.85f));
        bpSel->setPointRadius(0.0020, true);
        bpSel->setEnabled(false);
    }

    // ImGui panel — assigned to polyscope::state::userCallback
    polyscope::state::userCallback = [&vs, &sm]() {
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
        ImGui::SliderInt("Update every N", &vs.update_every, 1, 500);
        ImGui::Separator();

        // ── pick handling: click a MAT vertex to see its nmn_bplist ──────────
        if (polyscope::pick::haveSelection()) {
            auto [struct_ptr, local_idx] = polyscope::pick::getSelection();
            // Only react when the user clicked the live MAT vertex cloud
            if (struct_ptr == polyscope::getPointCloud("MAT Verts") &&
                local_idx < vs.idx_to_vid.size())
            {
                unsigned vid = vs.idx_to_vid[local_idx];
                if ((int)vid != vs.selected_vid) {
                    vs.selected_vid = (int)vid;
                    auto pts = BplistPositions(sm, vid);
                    if (!pts.empty()) {
                        auto* bpSel = polyscope::registerPointCloud("BPList selected", pts);
                        bpSel->setPointColor(glm::vec3(0.0f, 1.0f, 0.85f));
                        bpSel->setPointRadius(0.0020, true);
                        bpSel->setEnabled(true);
                    } else {
                        polyscope::getPointCloud("BPList selected")->setEnabled(false);
                    }
                }
            }
        }

        // Show info about selected vertex
        if (vs.selected_vid >= 0 &&
            (unsigned)vs.selected_vid < sm.vertices.size() &&
            sm.vertices[vs.selected_vid].first)
        {
            const auto& sv = *sm.vertices[vs.selected_vid].second;
            ImGui::Text("Selected vertex: %d", vs.selected_vid);
            ImGui::Text("  nmn_bplist size: %d", (int)sv.nmn_bplist.size());
            ImGui::Text("  (cyan points on mesh surface)");
            if (ImGui::Button("Clear selection")) {
                vs.selected_vid = -1;
                polyscope::getPointCloud("BPList selected")->setEnabled(false);
            }
        } else {
            ImGui::TextDisabled("Click a MAT vertex to see its bplist");
        }

        ImGui::Separator();
        ImGui::TextDisabled("Left panel layers:");
        ImGui::TextDisabled("  MAT Faces/Edges – orange/yellow");
        ImGui::TextDisabled("  Collapsed Edge  – red");
        ImGui::TextDisabled("  v1/v2 bplist    – blue/orange dots");
        ImGui::TextDisabled("  BPList selected – cyan dots");
        ImGui::PopItemWidth();
    };

    // Warm up the GLFW/WGL context: render a few frames so the window is fully
    // initialized before we start calling frameTick() from the collapse loop.
    for (int i = 0; i < 5; ++i)
        polyscope::frameTick();

    // Per-collapse callback — fired just before each MergeVertices call
    sm.on_collapse_cb = [&vs, &sm](
        unsigned raw_v1, const Wm4::Vector3d& v1p, double /*v1r*/,
        unsigned raw_v2, const Wm4::Vector3d& v2p, double /*v2r*/,
        const Sphere& result)
    {
        namespace ps = polyscope;
        using clock  = std::chrono::steady_clock;

        ++vs.collapse_count;

        // ── stop if the user closed the window ───────────────────────────────
        if (ps::windowRequestsClose()) {
            sm.on_collapse_cb = nullptr;   // disable further callbacks
            return;
        }

        // ── pause / step ─────────────────────────────────────────────────────
        if (vs.paused) {
            while (vs.paused && !vs.step_once && !ps::windowRequestsClose()) {
                std::this_thread::sleep_for(std::chrono::milliseconds(16));
                ps::frameTick();
            }
            vs.step_once = false;
        }

        // ── rate-limit: skip render if too soon and no structural update due ─
        static const auto kMinInterval = std::chrono::milliseconds(33); // ~30 fps
        bool do_update = (vs.collapse_count % vs.update_every == 0);
        bool time_ok   = (clock::now() - vs.last_frame) >= kMinInterval;

        if (!do_update && !time_ok)
            return;

        vs.last_frame = clock::now();

        if (do_update) {
            // Rebuild live MAT (v1 + v2 still active at this point — pre-merge)
            UpdateMatStructures(BuildMatArrays(sm), vs);

            // Highlighted edge: v1 → v2
            std::vector<std::array<double,3>>  ep = {
                {v1p.X(), v1p.Y(), v1p.Z()},
                {v2p.X(), v2p.Y(), v2p.Z()}
            };
            std::vector<std::array<size_t,2>> ec = {{0, 1}};
            auto* ce = ps::registerCurveNetwork("Collapsed Edge", ep, ec);
            ce->setColor(glm::vec3(1.0f, 0.15f, 0.15f));
            ce->setRadius(0.0030f, true);
            ce->setEnabled(true);

            // v1, v2, result sphere-centre markers
            std::vector<std::array<double,3>> p1 = {{v1p.X(), v1p.Y(), v1p.Z()}};
            std::vector<std::array<double,3>> p2 = {{v2p.X(), v2p.Y(), v2p.Z()}};
            std::vector<std::array<double,3>> pr = {{result.center.X(),
                                                     result.center.Y(),
                                                     result.center.Z()}};
            ps::getPointCloud("v1")->updatePointPositions(p1);
            ps::getPointCloud("v1")->setEnabled(true);
            ps::getPointCloud("v2")->updatePointPositions(p2);
            ps::getPointCloud("v2")->setEnabled(true);
            ps::getPointCloud("result")->updatePointPositions(pr);
            ps::getPointCloud("result")->setEnabled(true);

            // ── bplist of v1 and v2 on the input mesh surface ────────────────
            auto bp1_pts = BplistPositions(sm, raw_v1);
            auto bp2_pts = BplistPositions(sm, raw_v2);

            if (!bp1_pts.empty()) {
                auto* bp1 = ps::registerPointCloud("BPList v1", bp1_pts);
                bp1->setPointColor(glm::vec3(0.30f, 0.50f, 1.0f));
                bp1->setPointRadius(0.0020, true);
                bp1->setEnabled(true);
            }
            if (!bp2_pts.empty()) {
                auto* bp2 = ps::registerPointCloud("BPList v2", bp2_pts);
                bp2->setPointColor(glm::vec3(1.0f, 0.55f, 0.10f));
                bp2->setPointRadius(0.0020, true);
                bp2->setEnabled(true);
            }
        }

        ps::frameTick();
    };
}
#endif  // QMAT_WITH_POLYSCOPE

// Simple command line argument parsing
struct CLIOptions {
    std::string inputFile;
    std::string outputPrefix;
    int simplifyTarget = -1;  // -1 means no simplification
    double k = 0.00001;
    bool visualize = false;
    bool showHelp = false;
    bool valid = true;
    std::string errorMessage;
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
              << "  --k <value>        K factor for slab initialization (default: 0.00001)\n"
              << "  --output <prefix>  Output file prefix (default: input filename)\n"
#ifdef QMAT_WITH_POLYSCOPE
              << "  --visualize        Open live Polyscope viewer during simplification\n"
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

    // Step 3: Compute Delaunay Triangulation and Medial Axis
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

            // Label every vertex with its topological class before simplification.
            // Labels (mirrors check_non_manifold in read_qmat_output.py):
            //   REGULAR    – interior vertex (not on any boundary or non-manifold edge)
            //   BOUNDARY   – on a boundary edge (edge shared by exactly 1 face)
            //   NM_EDGE    – on a non-manifold edge (>2 faces), but not a corner
            //   NM_CORNER  – on both a boundary and a non-manifold edge
            // Simplify will only collapse an edge whose two endpoints share the same label.

            startTime = clock();
            shape.slab_mesh.CleanIsolatedVertices();

#ifdef QMAT_WITH_POLYSCOPE
            ViewerState vs;
            if (options.visualize)
                SetupSimplificationViewer(shape.slab_mesh, vs);
#endif

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

#ifdef QMAT_WITH_POLYSCOPE
            if (options.visualize) {
                // Clear the per-collapse callback so frameTick is no longer
                // driven by a running Simplify loop.
                shape.slab_mesh.on_collapse_cb = nullptr;
                // Register the final simplified MAT and hand control to the
                // Polyscope window for interactive inspection.
                UpdateMatStructures(BuildMatArrays(shape.slab_mesh), vs);
                if (polyscope::hasCurveNetwork("Collapsed Edge"))
                    polyscope::getCurveNetwork("Collapsed Edge")->setEnabled(false);
                polyscope::state::userCallback = []() {
                    ImGui::Text("Simplification complete.");
                    ImGui::Text("Close window to exit.");
                };
                std::cout << "Simplification done. Close the viewer window to exit.\n";
                polyscope::show();
            }
#endif
        }
    }

    std::cout << std::endl << "Done!" << std::endl;

    return 0;
}
