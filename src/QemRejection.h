#ifndef _QEM_REJECTION_H
#define _QEM_REJECTION_H

// ─────────────────────────────────────────────────────────────────────────────
// QemRejection — a SEPARATE rejection-record system for the VCG-faithful QEM
// edge-collapse simplifier (VcgQuadricSimplifier).
//
// This is deliberately independent of QMAT's own SlabMesh::RejectionReason /
// ReasonPrimitives machinery: it has its own enum, its own colour/name tables,
// and its own primitive struct, so the two paths never share state and the QEM
// viewer panel is wholly separate from the QMAT one.  (Per the user's choice of
// a "separate QemRejection system".)
//
// On the VCG path an edge collapse is only ever *hard*-rejected for one of three
// reasons, and only when the corresponding gate is enabled (which the Simplify
// dispatch turns on for this path):
//   • HardQualityCheckFailed  — the collapse would create a sliver triangle
//   • NormalFlipped           — the collapse would flip/fold a surviving face
//   • NonManifoldLinkCondition— the collapse would break manifoldness
// See md_files/qem/hard_rejection_causes.md.
//
// The whole header compiles to nothing unless ONLY_USE_QEM_CONDITION_CHECKS is
// defined, so non-flag builds are byte-identical.
// ─────────────────────────────────────────────────────────────────────────────

#if defined(ONLY_USE_QEM_CONDITION_CHECKS)

#include <array>
#include <vector>
#include <string>
#include <optional>
#include <utility>
#include <cstdint>

// One reason a VCG-path collapse was hard-rejected.  `None` is the sentinel for
// "never attempted / no rejection recorded" (rendered white in the viewer).
enum class QemRejectionReason : uint8_t {
	None = 0,
	HardQualityCheckFailed,    // FailsHardQualityCheck — post-collapse sliver
	NormalFlipped,             // CheckForFlip — face flip / >150° fold
	NonManifoldLinkCondition,  // IsFeasible (PreserveTopology) — link condition
	Count                      // sentinel — keep last
};

// RGB (0-255 per channel) for each reason.  Single source of truth, paired with
// QemRejectionReasonName below — both the viewer and any export use these.
inline std::array<uint8_t,3> QemRejectionReasonColorU8(QemRejectionReason rr)
{
	switch (rr) {
		case QemRejectionReason::HardQualityCheckFailed:   return { 128,   0, 128 }; // PURPLE
		case QemRejectionReason::NormalFlipped:            return { 255,   0,   0 }; // RED
		case QemRejectionReason::NonManifoldLinkCondition: return {   0, 255, 255 }; // CYAN
		case QemRejectionReason::None:
		case QemRejectionReason::Count:
		default:                                           return { 255, 255, 255 }; // WHITE
	}
}

// Human-readable name.  A switch (not an indexed array) so a new enum value can
// never read out of bounds; unhandled values return "???".
inline const char* QemRejectionReasonName(QemRejectionReason rr)
{
	switch (rr) {
		case QemRejectionReason::None:                     return "None";
		case QemRejectionReason::HardQualityCheckFailed:   return "HardQualityCheckFailed";
		case QemRejectionReason::NormalFlipped:            return "NormalFlipped";
		case QemRejectionReason::NonManifoldLinkCondition: return "NonManifoldLinkCondition";
		case QemRejectionReason::Count:                    break;
	}
	return "???";
}

// Geometry that explains a single rejection, so the viewer can highlight exactly
// what went wrong.  Mirrors the shape of QMAT's ReasonPrimitives but is its own
// type.  Slab vertex IDs index into SlabMesh::vertices; explicit world-space
// coordinates are used where the offending geometry no longer exists at current
// vertex positions (the post-collapse triangles).
struct QemReasonPrimitives {
	// Highlight by current position (drawn at sm.vertices[id].center).
	std::vector<unsigned>               vertices;  // slab vertex IDs
	std::vector<std::array<unsigned,2>> edges;     // [v0,v1] slab vertex-ID pairs
	std::vector<std::array<unsigned,3>> faces;     // [v0,v1,v2] slab vertex-ID triples

	// World-space position the collapse would have targeted (optPos).
	std::optional<std::array<double,3>> targ_ver;

	// NormalFlipped: the surviving face whose normal changed most, as explicit
	// world coords.  [0] = BEFORE collapse, [1] = AFTER collapse (one endpoint
	// moved to optPos).  Drawn as two coloured triangles (before/after).
	std::optional<std::array<std::array<std::array<double,3>,3>,2>> flipped_face;

	// Explicit post-collapse triangles in world coords (e.g. the sliver that
	// triggered HardQualityCheckFailed), drawn at their AFTER positions.
	std::vector<std::array<std::array<double,3>,3>> tris_after;

	// Named scalar metrics ("<label> = <value>" lines), e.g.
	// {"quality (post-collapse)", newQual}, {"threshold", HardQualityThr}.
	std::vector<std::pair<std::string,double>> metrics;

	// Free-text explanation of why this collapse was rejected (shown in the panel).
	std::string message;
};

#endif  // ONLY_USE_QEM_CONDITION_CHECKS
#endif  // _QEM_REJECTION_H
