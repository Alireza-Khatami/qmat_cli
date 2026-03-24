#include "MatCollapseHistory.h"
#include "SlabMesh.h"
#include <algorithm>

// ─────────────────────────────────────────────────────────────────────────────
//  Recording
// ─────────────────────────────────────────────────────────────────────────────

void MatCollapseHistory::Record(unsigned step,
                                unsigned vid_src1,
                                unsigned vid_src2,
                                unsigned vid_tgt,
                                std::array<double,3>  pos_src1,
                                std::array<double,3>  pos_src2,
                                std::vector<unsigned> bplist_src1,
                                std::vector<unsigned> bplist_src2,
                                std::vector<unsigned> bplist_after)
{
    records_.push_back({
        step,
        vid_src1, vid_src2, vid_tgt,
        pos_src1, pos_src2,
        std::move(bplist_src1),
        std::move(bplist_src2),
        std::move(bplist_after)
    });

    // Register both sources in the merge forest.
    parent_[vid_src1] = vid_tgt;
    parent_[vid_src2] = vid_tgt;

    // Auto-keyframe every keyframe_interval collapses.
    // TakeKeyframe() must be called explicitly by the caller because it
    // needs a reference to the full SlabMesh.
}

void MatCollapseHistory::TakeKeyframe(unsigned step, const SlabMesh& sm)
{
    MeshSnapshot snap;
    snap.step = step;

    for (unsigned i = 0; i < (unsigned)sm.vertices.size(); ++i) {
        if (!sm.vertices[i].first) continue;
        const SlabVertex* sv = sm.vertices[i].second;
        const auto& c = sv->sphere.center;

        MeshSnapshot::VertSnap vs;
        vs.vid  = i;
        vs.pos  = {c.X(), c.Y(), c.Z()};
        vs.type = static_cast<uint8_t>(sv->nmn_cluster_type);
        vs.bplist.assign(sv->nmn_bplist.begin(), sv->nmn_bplist.end());
        snap.verts.push_back(std::move(vs));
    }

    keyframes_.push_back(std::move(snap));
}

// ─────────────────────────────────────────────────────────────────────────────
//  Querying
// ─────────────────────────────────────────────────────────────────────────────

std::vector<unsigned> MatCollapseHistory::GetAncestors(unsigned vid) const
{
    // Build a children map: vid_tgt → {vid_src1, vid_src2, ...}
    std::unordered_map<unsigned, std::vector<unsigned>> children;
    for (const auto& kv : parent_)
        children[kv.second].push_back(kv.first);

    // DFS from vid downward toward the original leaf vertices.
    std::vector<unsigned> result;
    std::vector<unsigned> stack = {vid};
    while (!stack.empty()) {
        unsigned cur = stack.back(); stack.pop_back();
        result.push_back(cur);
        auto it = children.find(cur);
        if (it != children.end())
            for (unsigned child : it->second)
                stack.push_back(child);
    }
    return result;
}

std::vector<CollapseRecord> MatCollapseHistory::GetLineage(unsigned vid) const
{
    // Collect the full ancestor subtree for this vertex.
    const std::vector<unsigned> ancestors = GetAncestors(vid);
    std::unordered_map<unsigned,bool> in_subtree;
    in_subtree.reserve(ancestors.size());
    for (unsigned v : ancestors) in_subtree[v] = true;

    // Walk the global log and keep records where vid_tgt is in the subtree.
    // Records are already in chronological order (append-only log).
    std::vector<CollapseRecord> lineage;
    for (const auto& rec : records_)
        if (in_subtree.count(rec.vid_tgt))
            lineage.push_back(rec);

    return lineage;
}

const MeshSnapshot* MatCollapseHistory::GetKeyframeBefore(unsigned step) const
{
    if (keyframes_.empty()) return nullptr;

    // Keyframes are stored in ascending step order.
    // Return the last one whose step <= requested step.
    const MeshSnapshot* best = nullptr;
    for (const auto& kf : keyframes_) {
        if (kf.step <= step) best = &kf;
        else break;
    }
    return best ? best : &keyframes_.front();
}

void MatCollapseHistory::Clear()
{
    records_.clear();
    keyframes_.clear();
    parent_.clear();
}
