#pragma once
#include <vector>
#include <array>
#include <unordered_map>
#include <string>
#include <cstdint>

// Forward-declare to avoid circular dependency (SlabMesh.h includes this header).
class SlabMesh;

// ─────────────────────────────────────────────────────────────────────────────
//  MatCollapseHistory
//
//  Tracks the full history of edge collapses during MAT simplification.
//
//  Two complementary views are supported:
//
//  1. Global scrubbing (keyframes)
//     Full mesh snapshots are taken every `keyframe_interval` collapses.
//     GetKeyframeBefore(step) returns the nearest snapshot at or before that step.
//
//  2. Per-vertex lineage (merge tree)
//     MergeVertices() creates a fresh vid_tgt and deletes both vid_src1 and
//     vid_src2, so each CollapseRecord stores both removed vertices plus their
//     bplists before the merge and the combined bplist after.
//     GetAncestors(vid) returns the full subtree of original vertices that
//     were eventually merged into vid.
//     GetLineage(vid) returns those CollapseRecords in chronological order,
//     letting the viewer replay how the bplist grew at each step.
// ─────────────────────────────────────────────────────────────────────────────

// One record per edge collapse.
struct CollapseRecord {
    unsigned step;                       // global collapse index (0-based)
    unsigned vid_src1;                   // first vertex removed
    unsigned vid_src2;                   // second vertex removed
    unsigned vid_tgt;                    // new merged vertex
    std::array<double,3>  pos_src1;      // sphere centre of vid_src1 BEFORE merge
    std::array<double,3>  pos_src2;      // sphere centre of vid_src2 BEFORE merge
    std::vector<unsigned> bplist_src1;                    // nmn_bplist of vid_src1 BEFORE merge
    std::vector<unsigned> bplist_src2;                    // nmn_bplist of vid_src2 BEFORE merge
    std::vector<unsigned> bplist_after;                   // nmn_bplist of vid_tgt AFTER merge
    std::vector<std::vector<unsigned>> clusters_after;    // nmn_bplist_clusters of vid_tgt AFTER merge
};

// Lightweight snapshot of the active MAT at a particular step.
// Only active (non-tombstoned) vertices are stored.
struct MeshSnapshot {
    unsigned step;

    struct VertSnap {
        unsigned vid;
        std::array<double,3>  pos;    // normalized sphere centre
        uint8_t               type;   // SlabVertex::ClusterType cast to uint8_t
        std::vector<unsigned> bplist; // nmn_bplist at this step
    };

    std::vector<VertSnap> verts;
};

// ─────────────────────────────────────────────────────────────────────────────

class MatCollapseHistory {
public:
    // How often to take a full mesh snapshot (every N collapses).
    int keyframe_interval = 500;

    // ── Recording API (called from MinCostEdgeCollapse / MinCostBoundaryEdgeCollapse) ──

    // Record one edge collapse.
    // pos_src1/src2 and bplist_src1/src2 must be captured BEFORE MergeVertices().
    // bplist_after must be captured AFTER MergeVertices() from vid_tgt.
    void Record(unsigned step,
                unsigned vid_src1,
                unsigned vid_src2,
                unsigned vid_tgt,
                std::array<double,3>  pos_src1,
                std::array<double,3>  pos_src2,
                std::vector<unsigned> bplist_src1,
                std::vector<unsigned> bplist_src2,
                std::vector<unsigned> bplist_after,
                std::vector<std::vector<unsigned>> clusters_after);

    // Capture a full mesh snapshot at the given step.
    // Called automatically inside Record() every keyframe_interval collapses.
    void TakeKeyframe(unsigned step, const SlabMesh& sm);

    // ── Query API (called from main_cli.cpp / Polyscope UI) ──────────────────

    // All vertices (including vid itself) whose collapse chain eventually
    // reaches `vid`.  These are the leaves + internal nodes of the merge subtree.
    std::vector<unsigned> GetAncestors(unsigned vid) const;

    // All CollapseRecords along the lineage of `vid`, in chronological order.
    // A record is included when vid_tgt is in the ancestor subtree of `vid`.
    std::vector<CollapseRecord> GetLineage(unsigned vid) const;

    // Nearest keyframe snapshot at or before `step`.
    // Returns nullptr if no keyframe has been taken yet.
    const MeshSnapshot* GetKeyframeBefore(unsigned step) const;

    // Total number of collapses recorded so far.
    unsigned TotalSteps() const { return static_cast<unsigned>(records_.size()); }

    // Read-only access to the full collapse log.
    const std::vector<CollapseRecord>& Records() const { return records_; }

    // Read-only access to all keyframe snapshots.
    const std::vector<MeshSnapshot>& Keyframes() const { return keyframes_; }

    // Clear all recorded data (call before a new simplification run).
    void Clear();

private:
    std::vector<CollapseRecord> records_;   // one entry per collapse, in order
    std::vector<MeshSnapshot>   keyframes_; // sorted by step (ascending)

    // Merge forest: parent_[vid] = vid_tgt that absorbed vid.
    // Both vid_src1 and vid_src2 are registered here on every collapse.
    std::unordered_map<unsigned,unsigned> parent_;
};
