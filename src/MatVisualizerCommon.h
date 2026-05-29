// Shared types/palettes for the QMAT and VDE Polyscope visualizers.

#pragma once

#ifdef QMAT_WITH_POLYSCOPE

#include <array>
#include <chrono>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include "QemRejectionViz.h"   // for qemviz::State (macro-guarded inside)

class SlabMesh;

// ── MatArrays ────────────────────────────────────────────────────────────

// Contiguous arrays for the live MAT, fed to polyscope by UpdateMatStructures.
struct MatArrays {
    std::vector<std::array<double,3>>  verts;
    std::vector<std::array<size_t,2>>  edges;
    std::vector<std::array<size_t,2>>  boundary_edges; // both endpoints boundary
    std::vector<std::array<size_t,2>>  spike_edges;    // at least one T1_spike endpoint
    std::vector<std::array<size_t,3>>  faces;
    std::vector<std::array<size_t,3>>  spike_faces;    // at least one T1_spike vertex
    // Spike-free temp mesh: spike faces removed first, then spike edges.
    std::vector<std::array<double,3>>  ns_verts;
    std::vector<std::array<size_t,2>>  ns_edges;
    std::vector<std::array<size_t,3>>  ns_faces;
    std::vector<unsigned>              idx_to_vid;
    std::vector<unsigned>              idx_to_eid;
    std::vector<unsigned>              idx_to_fid;
    std::vector<std::array<float,3>>   vert_colors;          // by cluster type
    std::vector<std::array<double,3>>  unknown_ttype_verts;  // T0/T5 only
    // [0]=MS_Sheet [1]=MS_Seam [2]=MS_Boundary [3]=MS_Junction
    // [4]=MS_Sheet_Boundary [5]=MS_Seam_Boundary [6]=MS_Junction_Boundary
    std::array<std::vector<std::array<double,3>>, 7> cluster_filter_verts;
    std::vector<std::array<double,3>>  sharp_verts;
    // green=collapsible, red=mismatch, grey=not a struct edge
    std::vector<std::array<float,3>>   edge_structure_collapsible_colors;
    // by SlabEdge::topo_type (kEdgeTopoTypeColors)
    std::vector<std::array<float,3>>   edge_topo_type_colors;
    std::vector<std::array<float,3>>   face_struct_id_colors;
};

// ImGui edge-thickness slider overrides all curve network radii when true.
inline constexpr bool kModifyGlobalEdgeThickness = true;

// ── Palettes ─────────────────────────────────────────────────────────────

// Indexed by ClusterType enum value (T-types 0-6, MS_* types 7-14).
inline constexpr std::array<std::array<float,3>, 15> kClusterTypeColors = {{
    {0.9f, 0.0f, 0.9f},    // 0  T0                  invalid
    {0.1f, 0.1f, 0.1f},    // 1  T1_spike
    {0.0f, 0.85f, 1.0f},   // 2  T2
    {1.0f, 0.5f, 0.0f},    // 3  T3
    {1.0f, 0.1f, 0.1f},    // 4  T4
    {1.0f, 1.0f, 1.0f},    // 5  T5                  invalid
    {0.55f, 0.55f, 0.55f}, // 6  T1_non_spike
    {0.35f, 0.35f, 0.35f}, // 7  MS_Unknown
    {0.0f, 1.0f, 0.3f},    // 8  MS_Sheet
    {1.0f, 0.9f, 0.0f},    // 9  MS_Seam
    {0.0f, 0.5f, 1.0f},    // 10 MS_Boundary
    {1.0f, 0.0f, 0.5f},    // 11 MS_Junction
    {0.0f, 0.9f, 1.0f},    // 12 MS_Sheet_Boundary
    {1.0f, 0.35f, 0.0f},   // 13 MS_Seam_Boundary
    {0.6f, 0.0f, 1.0f},    // 14 MS_Junction_Boundary
}};

// Indexed by SlabEdge::TopoType enum value.  Matches MS_* hues above.
inline constexpr std::array<std::array<float,3>, 6> kEdgeTopoTypeColors = {{
    {0.55f, 0.55f, 0.55f}, // 0 Unknown
    {0.0f,  1.0f,  0.3f }, // 1 Sheet
    {1.0f,  0.9f,  0.0f }, // 2 Seam
    {0.0f,  0.5f,  1.0f }, // 3 Boundary
    {1.0f,  0.35f, 0.0f }, // 4 Seam_Boundary
    {0.9f,  0.0f,  0.9f }, // 5 Orphan
}};

// ── Color helpers ────────────────────────────────────────────────────────

std::array<float,3> HsvToRgb(float h, float s, float v);

// Distinct color per struct_id via golden-ratio hue.  Grey for sid < 0.
std::array<float,3> StructIdColor(int struct_id);

// ── ViewerState ──────────────────────────────────────────────────────────

struct ViewerState {
    int  collapse_count = 0;
    bool paused         = true;
    bool step_once      = false;
    int  update_every   = 1;
    // Adds collapse_delay_ms per collapse after this count.  -1 disables.
    int  collapse_delay_after = 3000;
    int  collapse_delay_ms    = 1;
    std::chrono::steady_clock::time_point last_frame =
        std::chrono::steady_clock::now();

    // "MAT Verts" cloud index → slab vertex id.  Rebuilt by BuildMatArrays.
    std::vector<unsigned> idx_to_vid;

    // "MAT Rejection Edges" pick index → slab edge id.  Curve-network pick
    // space is flat [nodes, edges); subtract rejection_node_count for edges.
    std::vector<unsigned> rejection_eid_order;
    size_t rejection_node_count = 0;

    int  selected_rejection_eid = -1;
    bool show_rejection_spheres = false;

#if defined(ONLY_USE_QEM_CONDITION_CHECKS)
    qemviz::State qem_viz;
#endif

    // Global radius for all curve networks when kModifyGlobalEdgeThickness.
    float edge_thickness = 0.0008f;

    bool show_struct_colors  = false;
    bool show_initial_struct = false;

    bool struct_collapsible_quantity_enabled = true;

    int selected_vid = -1;
    int selected_eid = -1;

    // "MAT Edges" edge slot → slab edge id.  Rebuilt by UpdateMatStructures.
    std::vector<unsigned> idx_to_eid;
    size_t mat_edge_node_count = 0;

    // Selection on "Initial MAT Faces" only (live face ids are unstable post
    // collapse).  Same id as the .ma file's face elem_id.
    int selected_init_fid = -1;

    // Surface-mesh pick order is [verts, faces, edges, halfedges]; face_slot
    // = local_idx - init_mat_face_vert_count.  Set once in Setup.
    std::vector<unsigned> init_idx_to_fid;
    size_t init_mat_face_vert_count = 0;

    std::string outputPrefix;

    enum class ColorMode { ClusterType, UnknownTType } color_mode = ColorMode::ClusterType;

    // -1 = show all; 8..14 = isolate that MS_* cluster type as own cloud.
    int cluster_filter = -1;

    int  last_picked_vid  = -1;
    std::chrono::steady_clock::time_point last_pick_time = {};
    static constexpr int kDoubleClickMs = 300;

    // (type, name) of structures enabled at click time, for restore on clear.
    std::vector<std::pair<std::string,std::string>> enabled_snapshot;
};

// ── Abstract visualizer base ─────────────────────────────────────────────
//
// QmatVisualizer and VdeVisualizer derive from this.  Both call Setup once
// to initialize their polyscope state and (where applicable) install
// per-collapse callbacks.  Show() runs polyscope::show().
class MatVisualizer {
public:
    virtual ~MatVisualizer() = default;

    virtual void Setup(SlabMesh& sm) = 0;
    void Show();   // polyscope::show()

    ViewerState&       State()       { return vs_; }
    const ViewerState& State() const { return vs_; }

protected:
    ViewerState vs_;
};

// ── Shared render helpers ────────────────────────────────────────────────

MatArrays BuildMatArrays(const SlabMesh& sm);
void UpdateMatStructures(const MatArrays& arr, ViewerState& vs);
void ApplyGlobalEdgeThickness(float r);
void RegisterInputMesh(const SlabMesh& sm);

// Scene-snapshot helpers used by the click/clear path.
void SnapshotEnabledPolyscopeStructures(
    std::vector<std::pair<std::string,std::string>>& snapshot);
void DisableAllPolyscopeStructures();
void RestoreEnabledPolyscopeStructures(
    const std::vector<std::pair<std::string,std::string>>& snapshot);

#endif  // QMAT_WITH_POLYSCOPE
