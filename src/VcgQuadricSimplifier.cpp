#include "VcgQuadricSimplifier.h"

#if defined(ONLY_USE_QEM_CONDITION_CHECKS)

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
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
double VcgQuadricSimplifier::ComputePriority(unsigned v0, unsigned v1, Wm4::Vector3d& outPos)
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
				newQual = std::min(newQual, QualityFace(p[0], p[1], p[2]));
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
	if (params.AreaCheck && (std::fabs(origArea - newArea) / (origArea + newArea) > 0.01))
		error = MAXVAL;
	if (params.HardQualityCheck && (newQual < params.HardQualityThr && newQual < origQual * 0.9))
		error = MAXVAL;
	if (params.HardNormalCheck && CheckForFlip(v0, v1, outPos))
		error = MAXVAL;

	return error;
}

// ─────────────────────────────────────────────────────────────────────────────
// CheckForFlip (tri_edge_collapse_quadric.h:459).  Returns true if, after the
// collapse, two surviving faces sharing an outgoing edge form a dihedral angle
// over 150° (or any surviving face becomes a near-zero-quality sliver).
// Only invoked when HardNormalCheck is on (off by default).
// ─────────────────────────────────────────────────────────────────────────────
bool VcgQuadricSimplifier::CheckForFlip(unsigned v0, unsigned v1, const Wm4::Vector3d& newPos)
{
	const double angleThrRad = 150.0 * 3.14159265358979323846 / 180.0;

	auto scan = [&](unsigned va, unsigned skip) -> bool {
		std::vector<std::pair<unsigned, Wm4::Vector3d>> edgeNorm;  // (other-vertex, face-normal)
		double maxAngle = 0.0;
		for (unsigned fid : m.vertices[va].second->faces_)
		{
			unsigned fv[3];
			if (!FaceVerts(fid, fv)) continue;
			if (fv[0] == skip || fv[1] == skip || fv[2] == skip) continue;

			Wm4::Vector3d p[3];
			for (int k = 0; k < 3; ++k) p[k] = (fv[k] == va) ? newPos : Pos(fv[k]);
			if (QualityFace(p[0], p[1], p[2]) < 0.01) return true;

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
		return maxAngle > angleThrRad;
	};

	if (scan(v0, v1)) return true;
	if (scan(v1, v0)) return true;
	return false;
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
bool VcgQuadricSimplifier::IsFeasible(unsigned v0, unsigned v1) const
{
	if (!params.PreserveTopology) return true;
	return !m.WouldCreateNonManifold(v0, v1);
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
	e.pri = ComputePriority(v0, v1, e.optPos);

	if (e.pri >= std::numeric_limits<double>::max()) return;  // maxAdmitErr

	heap.push_back(e);
	std::push_heap(heap.begin(), heap.end(), HeapLess);
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
		if (!IsFeasible(e.v0, e.v1)) continue;

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
