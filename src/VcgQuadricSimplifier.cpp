#include "VcgQuadricSimplifier.h"

#if defined(ONLY_USE_QEM_CONDITION_CHECKS)

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <set>
#include <string>
#include <utility>
#include <vector>
#include <iostream>
#include <Eigen/Dense>

// ─────────────────────────────────────────────────────────────────────────────
// Local geometry helpers — faithful ports of the vcg free functions used by the
// quadric collapse (vcglib/triangle3.h).
// ─────────────────────────────────────────────────────────────────────────────
namespace {

// vcg Quality(p0,p1,p2) == QualityFace = 2*Area / max(edge^2).  Range [0,0.866].
double QualityFace(const Wm4::Vector3d& p0,
                   const Wm4::Vector3d& p1,
                   const Wm4::Vector3d& p2)
{
	Wm4::Vector3d d10 = p1 - p0;
	Wm4::Vector3d d20 = p2 - p0;
	Wm4::Vector3d d12 = p1 - p2;
	double a = d10.Cross(d20).Length();   // 2*area
	if (a == 0.0) return 0.0;
	double b = d10.SquaredLength();
	if (b == 0.0) return 0.0;
	double t = d20.SquaredLength(); if (b < t) b = t;
	t = d12.SquaredLength();        if (b < t) b = t;
	return a / b;
}

// vcg NormalizedTriangleNormal — (p1-p0) x (p2-p0), normalised (zero on degenerate).
Wm4::Vector3d NormalizedTriangleNormal(const Wm4::Vector3d& p0,
                                       const Wm4::Vector3d& p1,
                                       const Wm4::Vector3d& p2)
{
	Wm4::Vector3d n = (p1 - p0).Cross(p2 - p0);
	double len = n.Length();
	if (len < 1e-30) return Wm4::Vector3d(0.0, 0.0, 0.0);
	return n / len;
}

// vcg math::Quadric<double>::Minimum — solve 2Ax+b=0 (i.e. Ax = -b/2) with
// Eigen fullPivLu and the same relative-error acceptance test (RelativeErrorThr
// = 1e-6).  Returns false (caller falls back to the midpoint) if the solution
// does not fit the system.
bool QuadricMinimum(const VcgQuadric& q, Wm4::Vector3d& x)
{
	Eigen::Matrix3d A;
	Eigen::Vector3d be;
	A << q.a[0], q.a[1], q.a[2],
	     q.a[1], q.a[3], q.a[4],
	     q.a[2], q.a[4], q.a[5];
	be << -q.b[0] / 2.0, -q.b[1] / 2.0, -q.b[2] / 2.0;

	Eigen::Vector3d xe = A.fullPivLu().solve(be);
	double error = (A * xe - be).norm();
	if (error > be.norm() * 1e-6) return false;

	x = Wm4::Vector3d(xe[0], xe[1], xe[2]);
	return true;
}

} // anonymous namespace

// ─────────────────────────────────────────────────────────────────────────────
// Small adjacency helpers.
// ─────────────────────────────────────────────────────────────────────────────
void VcgQuadricSimplifier::EnsureSized()
{
	if (vq.size() < m.vertices.size())    vq.resize(m.vertices.size());
	if (imark.size() < m.vertices.size()) imark.resize(m.vertices.size(), 0);
}

bool VcgQuadricSimplifier::FaceVerts(unsigned fid, unsigned out[3]) const
{
	if (fid >= m.faces.size() || !m.faces[fid].first) return false;
	const auto& vs = m.faces[fid].second->vertices_;
	if (vs.size() != 3) return false;
	int i = 0;
	for (unsigned v : vs) out[i++] = v;   // sorted ascending (std::set)
	return true;
}

double VcgQuadricSimplifier::EdgeLen(unsigned v0, unsigned v1) const
{
	return (Pos(v0) - Pos(v1)).Length();
}

// ─────────────────────────────────────────────────────────────────────────────
// InitQuadric (tri_edge_collapse_quadric.h:556).
//   • Zero every per-vertex quadric.
//   • For each face: build the (area-weighted) plane quadric and add it to the
//     three corner vertices.  For each face edge that is a border (exactly one
//     incident face) — or always, if QualityQuadric — add the perpendicular
//     border quadric scaled by BoundaryQuadricWeight (resp. QualityQuadricWeight).
//   • ScaleIndependent: ScaleFactor = 1e8 * (1/bboxDiag)^6 over all vertices.
// ─────────────────────────────────────────────────────────────────────────────
void VcgQuadricSimplifier::InitQuadrics()
{
	EnsureSized();

	for (size_t i = 0; i < m.vertices.size(); ++i)
		if (m.vertices[i].first)
			vq[i].SetZero();

	for (size_t fid = 0; fid < m.faces.size(); ++fid)
	{
		unsigned fv[3];
		if (!FaceVerts(static_cast<unsigned>(fid), fv)) continue;
		if (!m.vertices[fv[0]].first || !m.vertices[fv[1]].first || !m.vertices[fv[2]].first)
			continue;

		const Wm4::Vector3d& v0 = Pos(fv[0]);
		const Wm4::Vector3d& v1 = Pos(fv[1]);
		const Wm4::Vector3d& v2 = Pos(fv[2]);

		// vcg: dirArea = (v1-v0) x (v2-v0); facePlane direction = normalize(dirArea);
		//      area = |dirArea| (== 2*triangle area); offset = dir·v0.
		Wm4::Vector3d dirArea = (v1 - v0).Cross(v2 - v0);
		double area = dirArea.Length();
		if (area < 1e-30) continue;       // degenerate face contributes nothing
		Wm4::Vector3d nrm = dirArea / area;
		double off = nrm.Dot(v0);

		VcgQuadric q;
		q.ByPlane(nrm, off);
		if (params.UseArea) q *= area;

		for (int j = 0; j < 3; ++j)
			vq[fv[j]] += q;

		// Border / quality quadrics, one per face edge.
		for (int j = 0; j < 3; ++j)
		{
			unsigned a = fv[j];
			unsigned b = fv[(j + 1) % 3];

			unsigned eid;
			bool isB = false;
			if (m.Edge(a, b, eid) && m.edges[eid].first)
				isB = (m.edges[eid].second->faces_.size() == 1);

			if (!isB && !params.QualityQuadric) continue;

			// vcg: borderDir = faceNormal x normalize(edge); then amplify by the
			// weight (the *direction* is scaled, so the quadric scales by w²).
			Wm4::Vector3d edge = Pos(b) - Pos(a);
			double elen = edge.Length();
			if (elen < 1e-30) continue;
			Wm4::Vector3d borderDir = nrm.Cross(edge / elen);
			double w = isB ? params.BoundaryQuadricWeight : params.QualityQuadricWeight;
			borderDir = borderDir * w;
			double boff = borderDir.Dot(Pos(a));

			VcgQuadric bq;
			bq.ByPlane(borderDir, boff);
			vq[a] += bq;
			vq[b] += bq;
		}
	}

	// ScaleIndependent — bbox over all active vertices (vcg UpdateBounding::Box).
	if (params.ScaleIndependent)
	{
		Wm4::Vector3d mn( std::numeric_limits<double>::infinity(),
		                  std::numeric_limits<double>::infinity(),
		                  std::numeric_limits<double>::infinity());
		Wm4::Vector3d mx(-std::numeric_limits<double>::infinity(),
		                 -std::numeric_limits<double>::infinity(),
		                 -std::numeric_limits<double>::infinity());
		bool any = false;
		for (size_t i = 0; i < m.vertices.size(); ++i)
			if (m.vertices[i].first)
			{
				const Wm4::Vector3d& p = Pos(static_cast<unsigned>(i));
				mn.X() = std::min(mn.X(), p.X()); mn.Y() = std::min(mn.Y(), p.Y()); mn.Z() = std::min(mn.Z(), p.Z());
				mx.X() = std::max(mx.X(), p.X()); mx.Y() = std::max(mx.Y(), p.Y()); mx.Z() = std::max(mx.Z(), p.Z());
				any = true;
			}
		double diag = any ? (mx - mn).Length() : 0.0;
		if (diag > 1e-30) params.ScaleFactor = 1e8 * std::pow(1.0 / diag, 6.0);
		else              params.ScaleFactor = 1.0;
	}

	std::cerr << "[VcgQuadricSimplifier] InitQuadrics: ScaleFactor=" << params.ScaleFactor
	          << "  BoundaryQuadricWeight=" << params.BoundaryQuadricWeight
	          << "  UseArea=" << params.UseArea << "\n";
}

// ─────────────────────────────────────────────────────────────────────────────
// ComputePosition (tri_edge_collapse_quadric.h:158).
//   midpoint, unless OptimalPlacement and the summed quadric is non-flat at the
//   midpoint (Apply(mid) > 2*QuadricEpsilon), in which case the quadric minimum
//   is used.  If the minimum solve does not fit the system, fall back to the
//   midpoint (vcg leaves the position unchanged; the midpoint is the safe value).
// ─────────────────────────────────────────────────────────────────────────────
Wm4::Vector3d VcgQuadricSimplifier::ComputePosition(unsigned v0, unsigned v1)
{
	Wm4::Vector3d mid = (Pos(v0) + Pos(v1)) * 0.5;
	if (!params.OptimalPlacement) return Pos(v1);

	if ((vq[v0].Apply(mid) + vq[v1].Apply(mid)) > 2.0 * params.QuadricEpsilon)
	{
		VcgQuadric q = vq[v0];
		q += vq[v1];
		Wm4::Vector3d x;
		if (QuadricMinimum(q, x)) return x;
	}
	return mid;
}

// ─────────────────────────────────────────────────────────────────────────────
// ComputePriority (tri_edge_collapse_quadric.h:300).  Computes the heap error
// and the cached optimal position.  Reproduces all branches; with vcg/MeshLab
// defaults (QualityCheck on, everything else off) it reduces to
//   error = ScaleFactor * (q0+q1).Apply(optPos) / newQual.
// ─────────────────────────────────────────────────────────────────────────────
double VcgQuadricSimplifier::ComputePriority(unsigned v0, unsigned v1, Wm4::Vector3d& outPos,
                                             QemRejectionReason*  outReason,
                                             QemReasonPrimitives* outPrims)
{
	const Wm4::Vector3d oldP0 = Pos(v0);
	const Wm4::Vector3d oldP1 = Pos(v1);

	// ── Collect ORIGINAL normals / area / quality, exactly as vcg does. ──────
	std::vector<Wm4::Vector3d> origNormals;
	auto collectOrigNormals = [&](unsigned va, unsigned skip) {
		for (unsigned fid : m.vertices[va].second->faces_)
		{
			unsigned fv[3];
			if (!FaceVerts(fid, fv)) continue;
			if (fv[0] == skip || fv[1] == skip || fv[2] == skip) continue;  // skip faces with the partner
			origNormals.push_back(NormalizedTriangleNormal(Pos(fv[0]), Pos(fv[1]), Pos(fv[2])));
		}
	};
	if (params.NormalCheck) { collectOrigNormals(v0, v1); collectOrigNormals(v1, v0); }

	double origArea = 0.0;
	auto collectArea = [&](unsigned va, bool skipShared, unsigned partner) {
		for (unsigned fid : m.vertices[va].second->faces_)
		{
			unsigned fv[3];
			if (!FaceVerts(fid, fv)) continue;
			if (skipShared && (fv[0] == partner || fv[1] == partner || fv[2] == partner)) continue;
			origArea += (Pos(fv[1]) - Pos(fv[0])).Cross(Pos(fv[2]) - Pos(fv[0])).Length();  // DoubleArea
		}
	};
	if (params.AreaCheck) { collectArea(v0, false, v1); collectArea(v1, true, v0); }

	double origQual = std::numeric_limits<double>::max();
	auto collectOrigQual = [&](unsigned va, bool skipShared, unsigned partner) {
		for (unsigned fid : m.vertices[va].second->faces_)
		{
			unsigned fv[3];
			if (!FaceVerts(fid, fv)) continue;
			if (skipShared && (fv[0] == partner || fv[1] == partner || fv[2] == partner)) continue;
			origQual = std::min(origQual, QualityFace(Pos(fv[0]), Pos(fv[1]), Pos(fv[2])));
		}
	};
	if (params.HardQualityCheck) { collectOrigQual(v0, false, v1); collectOrigQual(v1, true, v0); }

	// ── Simulate the collapse: both endpoints move to the optimal position. ──
	outPos = ComputePosition(v0, v1);

	// Helper: quality / normal of a surviving face with one endpoint substituted.
	auto faceWithSubst = [&](unsigned fid, unsigned moved, Wm4::Vector3d p[3]) -> bool {
		unsigned fv[3];
		if (!FaceVerts(fid, fv)) return false;
		for (int k = 0; k < 3; ++k)
			p[k] = (fv[k] == moved) ? outPos : Pos(fv[k]);
		return true;
	};

	double MinCos = std::numeric_limits<double>::max();
	if (params.NormalCheck)
	{
		int i = 0;
		auto scanNormals = [&](unsigned va, unsigned skip) {
			for (unsigned fid : m.vertices[va].second->faces_)
			{
				unsigned fv[3];
				if (!FaceVerts(fid, fv)) continue;
				if (fv[0] == skip || fv[1] == skip || fv[2] == skip) continue;
				Wm4::Vector3d p[3];
				faceWithSubst(fid, va, p);
				Wm4::Vector3d nn = NormalizedTriangleNormal(p[0], p[1], p[2]);
				if (i < (int)origNormals.size())
					MinCos = std::min(MinCos, nn.Dot(origNormals[i++]));
			}
		};
		scanNormals(v0, v1);
		scanNormals(v1, v0);
	}

	double newQual = std::numeric_limits<double>::max();
	std::array<std::array<double,3>,3> worstQualTri{};
	bool haveWorstQualTri = false;
	if (params.QualityCheck)
	{
		auto scanQual = [&](unsigned va, unsigned skip) {
			for (unsigned fid : m.vertices[va].second->faces_)
			{
				unsigned fv[3];
				if (!FaceVerts(fid, fv)) continue;
				if (fv[0] == skip || fv[1] == skip || fv[2] == skip) continue;  // destroyed shared face
				Wm4::Vector3d p[3];
				faceWithSubst(fid, va, p);
				double q = QualityFace(p[0], p[1], p[2]);
				if (q < newQual) {
					newQual = q;
					worstQualTri = { { {p[0].X(), p[0].Y(), p[0].Z()},
					                   {p[1].X(), p[1].Y(), p[1].Z()},
					                   {p[2].X(), p[2].Y(), p[2].Z()} } };
					haveWorstQualTri = true;
				}
			}
		};
		scanQual(v0, v1);
		scanQual(v1, v0);
	}

	double newArea = 0.0;
	if (params.AreaCheck)
	{
		auto scanArea = [&](unsigned va, bool skipShared, unsigned partner) {
			for (unsigned fid : m.vertices[va].second->faces_)
			{
				unsigned fv[3];
				if (!FaceVerts(fid, fv)) continue;
				if (skipShared && (fv[0] == partner || fv[1] == partner || fv[2] == partner)) continue;
				Wm4::Vector3d p[3];
				faceWithSubst(fid, va, p);
				newArea += (p[1] - p[0]).Cross(p[2] - p[0]).Length();
			}
		};
		scanArea(v0, false, v1);
		scanArea(v1, true, v0);
	}

	// ── QuadErr = ScaleFactor * (q0+q1).Apply(optPos). ───────────────────────
	VcgQuadric qq = vq[v0];
	qq += vq[v1];
	double QuadErr = params.ScaleFactor * qq.Apply(outPos);

	if (newQual > params.QualityThr) newQual = params.QualityThr;

	if (params.NormalCheck)
	{
		if (MinCos > params.CosineThr) MinCos = params.CosineThr;
		MinCos = std::fabs((MinCos + 1.0) / 2.0);   // → [0,1], 0 = worst
	}

	QuadErr = std::max(QuadErr, params.QuadricEpsilon);
	if (QuadErr <= params.QuadricEpsilon)
		QuadErr *= EdgeLen(v0, v1);

	// UseVertexWeight: W==1 for all vertices in our port → no-op (kept faithful).
	if (params.UseVertexWeight) QuadErr *= 1.0;

	double error;
	if (!params.QualityCheck && !params.NormalCheck) error = QuadErr;
	else if ( params.QualityCheck && !params.NormalCheck) error = QuadErr / newQual;
	else if (!params.QualityCheck &&  params.NormalCheck) error = QuadErr / MinCos;
	else                                                  error = QuadErr / (newQual * MinCos);

	const double MAXVAL = std::numeric_limits<double>::max();
	const bool   capture = (outReason && outPrims);

	if (params.AreaCheck && (std::fabs(origArea - newArea) / (origArea + newArea) > 0.01))
		error = MAXVAL;

	if (params.HardQualityCheck && (newQual < params.HardQualityThr && newQual < origQual * 0.9))
	{
		if (!params.DiagnoseOnly) error = MAXVAL;   // DiagnoseOnly: colour but don't veto
		if (capture) {
			*outReason = QemRejectionReason::HardQualityCheckFailed;
			outPrims->vertices = { v0, v1 };
			outPrims->edges    = { { v0, v1 } };
			outPrims->targ_ver = { outPos.X(), outPos.Y(), outPos.Z() };
			if (haveWorstQualTri) outPrims->tris_after.push_back(worstQualTri);
			outPrims->metrics  = { {"quality (post-collapse)", newQual},
			                       {"hard quality threshold",  params.HardQualityThr},
			                       {"orig quality x0.9",       origQual * 0.9} };
			outPrims->message  =
				"Collapse would create a sliver triangle: the worst surviving face's "
				"quality drops to " + std::to_string(newQual) + ", below the hard threshold "
				+ std::to_string(params.HardQualityThr) + " and below 90% of the original "
				"worst quality (" + std::to_string(origQual * 0.9) + "). Yellow = the sliver.";
		}
	}

	if (params.HardNormalCheck && CheckForFlip(v0, v1, outPos))   // fast path (no capture)
	{
		if (!params.DiagnoseOnly) error = MAXVAL;   // DiagnoseOnly: colour but don't veto
		{
			if (capture) {
				// Second pass only on a real flip: capture the offending face geometry.
				std::array<std::array<std::array<double,3>,3>,2> flipped;
				CheckForFlip(v0, v1, outPos, &flipped);
				*outReason = QemRejectionReason::NormalFlipped;
				outPrims->vertices     = { v0, v1 };
				outPrims->edges        = { { v0, v1 } };
				outPrims->targ_ver     = { outPos.X(), outPos.Y(), outPos.Z() };
				outPrims->flipped_face = flipped;
				outPrims->tris_after.clear();   // flip viz supersedes any sliver viz
				outPrims->metrics.clear();
				outPrims->message =
					"Collapse would flip a surviving face: moving the merged vertex to the "
					"optimal position turns a face inside-out (dihedral > 150 deg) or creates "
					"a near-zero-quality sliver. Red = face before, orange = face after.";
			}
		}
	}

	return error;
}

// ─────────────────────────────────────────────────────────────────────────────
// CheckForFlip (tri_edge_collapse_quadric.h:459).  Returns true if, after the
// collapse, two surviving faces sharing an outgoing edge form a dihedral angle
// over 150° (or any surviving face becomes a near-zero-quality sliver).
// Only invoked when HardNormalCheck is on (off by default).
// ─────────────────────────────────────────────────────────────────────────────
bool VcgQuadricSimplifier::CheckForFlip(unsigned v0, unsigned v1, const Wm4::Vector3d& newPos,
                                        std::array<std::array<std::array<double,3>,3>,2>* outFlipped)
{
	const double angleThrRad = 150.0 * 3.14159265358979323846 / 180.0;

	// For the optional flip visualization: track the surviving face whose normal
	// changes most (smallest dot of before-normal with after-normal) and store
	// its before/after triangle into *outFlipped.
	double worstDot = 2.0;
	auto considerFlipViz = [&](const unsigned fv[3], unsigned va) {
		if (!outFlipped) return;
		Wm4::Vector3d bp[3], ap[3];
		for (int k = 0; k < 3; ++k) {
			bp[k] = Pos(fv[k]);
			ap[k] = (fv[k] == va) ? newPos : Pos(fv[k]);
		}
		Wm4::Vector3d bn = NormalizedTriangleNormal(bp[0], bp[1], bp[2]);
		Wm4::Vector3d an = NormalizedTriangleNormal(ap[0], ap[1], ap[2]);
		double d = bn.Dot(an);
		if (d < worstDot) {
			worstDot = d;
			(*outFlipped)[0] = { { {bp[0].X(), bp[0].Y(), bp[0].Z()},
			                       {bp[1].X(), bp[1].Y(), bp[1].Z()},
			                       {bp[2].X(), bp[2].Y(), bp[2].Z()} } };
			(*outFlipped)[1] = { { {ap[0].X(), ap[0].Y(), ap[0].Z()},
			                       {ap[1].X(), ap[1].Y(), ap[1].Z()},
			                       {ap[2].X(), ap[2].Y(), ap[2].Z()} } };
		}
	};

	auto scan = [&](unsigned va, unsigned skip) -> bool {
		std::vector<std::pair<unsigned, Wm4::Vector3d>> edgeNorm;  // (other-vertex, face-normal)
		double maxAngle = 0.0;
		bool   sliver = false;
		for (unsigned fid : m.vertices[va].second->faces_)
		{
			unsigned fv[3];
			if (!FaceVerts(fid, fv)) continue;
			if (fv[0] == skip || fv[1] == skip || fv[2] == skip) continue;

			considerFlipViz(fv, va);

			Wm4::Vector3d p[3];
			for (int k = 0; k < 3; ++k) p[k] = (fv[k] == va) ? newPos : Pos(fv[k]);
			if (QualityFace(p[0], p[1], p[2]) < 0.01) { sliver = true; if (!outFlipped) return true; }

			Wm4::Vector3d n = NormalizedTriangleNormal(p[0], p[1], p[2]);
			for (int k = 0; k < 3; ++k)
			{
				if (fv[k] == va) continue;
				unsigned other = fv[k];
				bool found = false;
				for (auto& en : edgeNorm)
					if (en.first == other)
					{
						double d = std::max(-1.0, std::min(1.0, n.Dot(en.second)));
						maxAngle = std::max(maxAngle, std::acos(d));
						found = true;
						break;
					}
				if (!found) edgeNorm.push_back({other, n});
			}
		}
		return sliver || (maxAngle > angleThrRad);
	};

	// Fast path (no capture): short-circuit exactly like the original.
	if (!outFlipped)
	{
		if (scan(v0, v1)) return true;
		if (scan(v1, v0)) return true;
		return false;
	}
	// Capture path: run both scans so the flip viz reflects the worst face.
	bool flipped = false;
	if (scan(v0, v1)) flipped = true;
	if (scan(v1, v0)) flipped = true;
	return flipped;
}

// ─────────────────────────────────────────────────────────────────────────────
// IsUpToDate (tri_edge_collapse.h:198) — an element is stale if either endpoint
// has been deleted or has been touched (IMark bumped) since the element was made.
// ─────────────────────────────────────────────────────────────────────────────
bool VcgQuadricSimplifier::IsUpToDate(const HeapElem& e) const
{
	if (e.v0 >= m.vertices.size() || e.v1 >= m.vertices.size()) return false;
	if (!m.vertices[e.v0].first || !m.vertices[e.v1].first)     return false;
	if (e.localMark < imark[e.v0] || e.localMark < imark[e.v1]) return false;
	return true;
}

// ─────────────────────────────────────────────────────────────────────────────
// IsFeasible (tri_edge_collapse_quadric.h:149) — link condition, only when
// PreserveTopology is set (off in MeshLab defaults).  Maps onto SlabMesh's own
// link-condition port WouldCreateNonManifold.
// ─────────────────────────────────────────────────────────────────────────────
bool VcgQuadricSimplifier::IsFeasible(unsigned v0, unsigned v1, QemReasonPrimitives* outPrims) const
{
	if (!params.PreserveTopology) return true;
	if (!m.WouldCreateNonManifold(v0, v1)) return true;

	if (outPrims) {
		outPrims->vertices = { v0, v1 };
		outPrims->edges    = { { v0, v1 } };
		// Offending neighbourhood: vertices shared by the one-rings of v0 and v1
		// (beyond the partner itself) — the link-condition violation.
		std::set<unsigned> r0, r1;
		auto ring = [&](unsigned va, std::set<unsigned>& out) {
			if (va >= m.vertices.size() || !m.vertices[va].first) return;
			for (unsigned fid : m.vertices[va].second->faces_) {
				unsigned fv[3];
				if (!FaceVerts(fid, fv)) continue;
				for (int k = 0; k < 3; ++k) if (fv[k] != va) out.insert(fv[k]);
			}
		};
		ring(v0, r0); ring(v1, r1);
		for (unsigned s : r0)
			if (s != v1 && s != v0 && r1.count(s)) outPrims->vertices.push_back(s);
		outPrims->message =
			"Collapse would break manifoldness (link condition): the one-rings of the two "
			"endpoints share vertices other than the two opposite the collapsed edge, so "
			"merging them would create a non-manifold edge/vertex. Highlighted vertices are "
			"the shared neighbours.";
	}
	return false;
}

// ─────────────────────────────────────────────────────────────────────────────
// AddCollapseToHeap (tri_edge_collapse_quadric.h:500) — compute the priority for
// (v0,v1), and push it unless the error exceeds the admissible max.  Symmetric
// collapse (OptimalPlacement) ⇒ a single directed entry per pair.
// ─────────────────────────────────────────────────────────────────────────────
void VcgQuadricSimplifier::AddCollapseToHeap(unsigned v0, unsigned v1)
{
	if (v0 == v1) return;
	if (v0 >= m.vertices.size() || v1 >= m.vertices.size()) return;
	if (!m.vertices[v0].first || !m.vertices[v1].first) return;

	HeapElem e;
	e.v0 = v0; e.v1 = v1;
	e.localMark = globalMark;

	QemRejectionReason  reason = QemRejectionReason::None;
	QemReasonPrimitives prims;
	e.pri = ComputePriority(v0, v1, e.optPos, &reason, &prims);

	if (e.pri >= std::numeric_limits<double>::max()) {     // maxAdmitErr — hard veto
		if (reason != QemRejectionReason::None)
			RecordRejection(v0, v1, reason, std::move(prims));
		return;                                            // dropped (not DiagnoseOnly)
	}

	// Finite priority → this candidate WILL go on the heap.  In DiagnoseOnly mode a
	// failing check still set `reason` (but not +inf): colour the edge yet keep it.
	// Otherwise the edge is genuinely fine — clear any stale rejection mark.
	if (reason != QemRejectionReason::None)
		RecordRejection(v0, v1, reason, std::move(prims));
	else
		ClearRejection(v0, v1);

	heap.push_back(e);
	std::push_heap(heap.begin(), heap.end(), HeapLess);
}

// ─────────────────────────────────────────────────────────────────────────────
// RecordRejection / ClearRejection — resolve the slab edge id for (v0,v1) and
// store / erase its QEM rejection in SlabMesh's qem_edge_* maps (read by the
// QEM rejection viewer in main_cli).
// ─────────────────────────────────────────────────────────────────────────────
void VcgQuadricSimplifier::RecordRejection(unsigned v0, unsigned v1,
                                           QemRejectionReason reason, QemReasonPrimitives prims)
{
	unsigned eid;
	if (!m.Edge(v0, v1, eid)) return;   // no slab edge → nothing to colour
	m.qem_edge_last_rejection[eid]    = reason;
	m.qem_edge_reason_primitives[eid] = std::move(prims);
}

void VcgQuadricSimplifier::ClearRejection(unsigned v0, unsigned v1)
{
	unsigned eid;
	if (!m.Edge(v0, v1, eid)) return;
	m.qem_edge_last_rejection.erase(eid);
	m.qem_edge_reason_primitives.erase(eid);
}

// ─────────────────────────────────────────────────────────────────────────────
// Execute (tri_edge_collapse_quadric.h:182) — q(v1) += q(v0); collapse v0→v1 to
// the optimal position.  We accumulate the merged quadric, run SlabMesh's
// topological merge (which allocates a fresh target id), then move the merged
// vertex onto optPos.
// ─────────────────────────────────────────────────────────────────────────────
unsigned VcgQuadricSimplifier::Execute(unsigned v0, unsigned v1, const Wm4::Vector3d& optPos)
{
	VcgQuadric merged = vq[v0];
	merged += vq[v1];

	const double r0 = m.vertices[v0].second->sphere.radius;
	const double r1 = m.vertices[v1].second->sphere.radius;
	const double rt = 0.5 * (r0 + r1);

#ifdef QMAT_WITH_POLYSCOPE
	if (m.on_collapse_cb)
	{
		Sphere s; s.center = optPos; s.radius = rt;
		m.on_collapse_cb(v0, Pos(v0), r0, v1, Pos(v1), r1, s);
	}
#endif

	unsigned vid_tgt = (unsigned)-1;
	if (!m.MergeVertices(v0, v1, vid_tgt)) return (unsigned)-1;

	// NOTE: v0 and v1 are now deleted (MergeVertices freed them) — do not touch
	// them past this point.
	EnsureSized();
	vq[vid_tgt] = merged;
	m.vertices[vid_tgt].second->sphere.center = optPos;
	m.vertices[vid_tgt].second->sphere.radius = rt;
	return vid_tgt;
}

// ─────────────────────────────────────────────────────────────────────────────
// UpdateHeap (tri_edge_collapse_quadric.h:523) — bump GlobalMark and the IMark
// of the surviving vertex and its ring, then re-add every collapse incident to
// the surviving vertex plus the outer-ring edges (pairs of ring vertices that
// share a face with it).  Older entries touching the ring are now stale.
// ─────────────────────────────────────────────────────────────────────────────
void VcgQuadricSimplifier::UpdateHeap(unsigned vid_tgt)
{
	++globalMark;
	EnsureSized();
	imark[vid_tgt] = globalMark;

	// Ring = vertices sharing a face with vid_tgt (vcg iterates faces via VFIterator).
	std::set<unsigned> ring;
	for (unsigned fid : m.vertices[vid_tgt].second->faces_)
	{
		unsigned fv[3];
		if (!FaceVerts(fid, fv)) continue;
		for (int k = 0; k < 3; ++k)
			if (fv[k] != vid_tgt) ring.insert(fv[k]);
	}
	for (unsigned r : ring) imark[r] = globalMark;

	// Re-add edges from vid_tgt to each ring vertex.
	for (unsigned r : ring)
		AddCollapseToHeap(vid_tgt, r);

	// Re-add the outer-ring edges (the two non-target vertices of each incident
	// face).  Dedup so a shared outer edge is not pushed twice.
	std::set<std::pair<unsigned, unsigned>> seen;
	for (unsigned fid : m.vertices[vid_tgt].second->faces_)
	{
		unsigned fv[3];
		if (!FaceVerts(fid, fv)) continue;
		unsigned a = (unsigned)-1, b = (unsigned)-1;
		for (int k = 0; k < 3; ++k)
			if (fv[k] != vid_tgt) { (a == (unsigned)-1 ? a : b) = fv[k]; }
		if (a == (unsigned)-1 || b == (unsigned)-1) continue;
		auto key = std::minmax(a, b);
		if (seen.insert({key.first, key.second}).second)
			AddCollapseToHeap(a, b);
	}
}

// ─────────────────────────────────────────────────────────────────────────────
// ClearHeap (local_optimization.h:237) — drop stale elements when the heap grows
// past SimplexNumber * HeapSimplexRatio, then rebuild.
// ─────────────────────────────────────────────────────────────────────────────
void VcgQuadricSimplifier::ClearHeap()
{
	std::vector<HeapElem> kept;
	kept.reserve(heap.size());
	for (const HeapElem& e : heap)
		if (IsUpToDate(e)) kept.push_back(e);
	heap.swap(kept);
	std::make_heap(heap.begin(), heap.end(), HeapLess);
}

// ─────────────────────────────────────────────────────────────────────────────
// Run — vcg LocalOptimization::Init + DoOptimization.
//   • seed the heap with one entry per face edge (vcg seeds via faces, so edges
//     with no incident face are never collapsed — matching vcg on a tri mesh);
//   • pop the cheapest, skip if stale, skip if infeasible, else Execute +
//     UpdateHeap.  Stop after maxCollapses collapses (Simplify's threshold) or
//     when only one vertex remains.
// ─────────────────────────────────────────────────────────────────────────────
unsigned VcgQuadricSimplifier::Run(int maxCollapses)
{
	params.CosineThr = std::cos(params.NormalThrRad);

	globalMark = 0;
	imark.assign(m.vertices.size(), 0);
	vq.assign(m.vertices.size(), VcgQuadric());

	InitQuadrics();

	// Seed: every distinct edge that belongs to at least one face.
	heap.clear();
	{
		std::set<std::pair<unsigned, unsigned>> seen;
		for (size_t fid = 0; fid < m.faces.size(); ++fid)
		{
			unsigned fv[3];
			if (!FaceVerts(static_cast<unsigned>(fid), fv)) continue;
			for (int j = 0; j < 3; ++j)
			{
				unsigned a = fv[j], b = fv[(j + 1) % 3];
				auto key = std::minmax(a, b);
				if (seen.insert({key.first, key.second}).second)
					AddCollapseToHeap(key.first, key.second);
			}
		}
	}
	std::make_heap(heap.begin(), heap.end(), HeapLess);

	const float ratio = params.OptimalPlacement ? 4.0f : 8.0f;  // HeapSimplexRatio

	std::cerr << "[VcgQuadricSimplifier] seeded heap with " << heap.size()
	          << " candidate collapses; target collapses = " << maxCollapses
	          << ", start MAT vertices = " << m.numVertices << "\n";

	unsigned collapses = 0;
	while (collapses < (unsigned)maxCollapses && m.numVertices > 1 && !heap.empty())
	{
		if (heap.size() > (size_t)(m.numFaces * ratio)) ClearHeap();
		if (heap.empty()) break;

		std::pop_heap(heap.begin(), heap.end(), HeapLess);
		HeapElem e = heap.back();
		heap.pop_back();

		if (!IsUpToDate(e)) continue;
		{
			QemReasonPrimitives prims;
			if (!IsFeasible(e.v0, e.v1, &prims)) {
				RecordRejection(e.v0, e.v1, QemRejectionReason::NonManifoldLinkCondition,
				                std::move(prims));
				// Veto only when not diagnosing; DiagnoseOnly colours it yet collapses.
				if (!params.DiagnoseOnly) continue;
			}
		}

		unsigned vid_tgt = Execute(e.v0, e.v1, e.optPos);
		if (vid_tgt == (unsigned)-1) continue;

		UpdateHeap(vid_tgt);
		++collapses;
	}

	std::cerr << "[VcgQuadricSimplifier] done: " << collapses
	          << " collapses, MAT vertices = " << m.numVertices << "\n";
	return collapses;
}

#endif  // ONLY_USE_QEM_CONDITION_CHECKS
