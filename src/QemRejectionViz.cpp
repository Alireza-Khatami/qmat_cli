#include "QemRejectionViz.h"

#if defined(ONLY_USE_QEM_CONDITION_CHECKS) && defined(QMAT_WITH_POLYSCOPE)

#include "SlabMesh.h"          // SlabMesh, QemRejectionReason, QemReasonPrimitives

#include "polyscope/polyscope.h"
#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"
#include "polyscope/surface_mesh.h"
#include "imgui.h"

#include <array>
#include <map>
#include <unordered_map>

namespace qemviz {

namespace ps = polyscope;

// QemRejectionReason colour (0-1 floats) — wraps the 0-255 table in QemRejection.h.
static std::array<float,3> ReasonColorF(QemRejectionReason rr)
{
	auto c = QemRejectionReasonColorU8(rr);
	return { c[0] / 255.0f, c[1] / 255.0f, c[2] / 255.0f };
}

// World position of a slab vertex (zero if deleted/out of range).
static std::array<double,3> VertPos(const SlabMesh& sm, unsigned vid)
{
	if (vid < sm.vertices.size() && sm.vertices[vid].first) {
		const auto& c = sm.vertices[vid].second->sphere.center;
		return {c.X(), c.Y(), c.Z()};
	}
	return {0.0, 0.0, 0.0};
}

// ─────────────────────────────────────────────────────────────────────────────
void UpdateEdgeColors(const SlabMesh& sm, State& st)
{
	std::map<unsigned, size_t> vid_map;
	std::vector<std::array<double,3>> nodes;
	std::vector<std::array<size_t,2>> segs;
	std::vector<std::array<float,3>>  edge_colors;
	st.eid_order.clear();

	for (unsigned i = 0; i < (unsigned)sm.vertices.size(); ++i) {
		if (!sm.vertices[i].first) continue;
		vid_map[i] = nodes.size();
		const auto& c = sm.vertices[i].second->sphere.center;
		nodes.push_back({c.X(), c.Y(), c.Z()});
	}

	for (unsigned i = 0; i < (unsigned)sm.edges.size(); ++i) {
		if (!sm.edges[i].first) continue;
		size_t a = vid_map.at(sm.edges[i].second->vertices_.first);
		size_t b = vid_map.at(sm.edges[i].second->vertices_.second);
		segs.push_back({a, b});
		st.eid_order.push_back(i);

		auto it = sm.qem_edge_last_rejection.find(i);
		edge_colors.push_back(it != sm.qem_edge_last_rejection.end()
			? ReasonColorF(it->second)
			: std::array<float,3>{1.0f, 1.0f, 1.0f});  // white = never rejected
	}

	if (segs.empty()) return;

	st.node_count = nodes.size();

	bool en = ps::hasCurveNetwork("MAT QEM Rejection Edges")
	          ? ps::getCurveNetwork("MAT QEM Rejection Edges")->isEnabled() : false;
	auto* cn = ps::registerCurveNetwork("MAT QEM Rejection Edges", nodes, segs);
	cn->setRadius(0.0012f, true);
	cn->setEnabled(en);
	cn->addEdgeColorQuantity("QEM Rejection Reason", edge_colors)->setEnabled(true);
}

// ─────────────────────────────────────────────────────────────────────────────
void ClearPrimitives()
{
	const char* pcs[] = { "QEM Rejection Verts", "QEM Rejection Target", "QEM Rejection Spheres" };
	for (const char* n : pcs) if (ps::hasPointCloud(n)) ps::getPointCloud(n)->setEnabled(false);
	if (ps::hasCurveNetwork("QEM Rejection Edges")) ps::getCurveNetwork("QEM Rejection Edges")->setEnabled(false);
	const char* sms[] = { "QEM Rejection Faces", "QEM Sliver Face",
	                      "QEM Flipped Face Before", "QEM Flipped Face After" };
	for (const char* n : sms) if (ps::hasSurfaceMesh(n)) ps::getSurfaceMesh(n)->setEnabled(false);
}

// ─────────────────────────────────────────────────────────────────────────────
void ShowPrimitives(const SlabMesh& sm, unsigned eid, bool show_spheres)
{
	auto it = sm.qem_edge_reason_primitives.find(eid);
	if (it == sm.qem_edge_reason_primitives.end()) { ClearPrimitives(); return; }
	const auto& prims = it->second;

	// ── Vertices (yellow) ────────────────────────────────────────────────────
	if (!prims.vertices.empty()) {
		std::vector<std::array<double,3>> pts;
		for (unsigned vid : prims.vertices) pts.push_back(VertPos(sm, vid));
		auto* pc = ps::registerPointCloud("QEM Rejection Verts", pts);
		pc->setPointColor(glm::vec3(1.0f, 1.0f, 0.0f));
		pc->setPointRadius(0.006, true);
		pc->setEnabled(true);
	} else if (ps::hasPointCloud("QEM Rejection Verts")) {
		ps::getPointCloud("QEM Rejection Verts")->setEnabled(false);
	}

	// ── Edges (yellow) ───────────────────────────────────────────────────────
	if (!prims.edges.empty()) {
		std::unordered_map<unsigned,size_t> remap;
		std::vector<std::array<double,3>> nodes;
		std::vector<std::array<size_t,2>> segs;
		auto addV = [&](unsigned vid) -> size_t {
			auto jt = remap.find(vid);
			if (jt != remap.end()) return jt->second;
			size_t idx = nodes.size(); remap[vid] = idx; nodes.push_back(VertPos(sm, vid)); return idx;
		};
		for (const auto& e : prims.edges) segs.push_back({addV(e[0]), addV(e[1])});
		auto* cn = ps::registerCurveNetwork("QEM Rejection Edges", nodes, segs);
		cn->setColor(glm::vec3(1.0f, 1.0f, 0.0f));
		cn->setRadius(0.004f, true);
		cn->setEnabled(true);
	} else if (ps::hasCurveNetwork("QEM Rejection Edges")) {
		ps::getCurveNetwork("QEM Rejection Edges")->setEnabled(false);
	}

	// ── Faces by id (yellow translucent) ─────────────────────────────────────
	if (!prims.faces.empty()) {
		std::unordered_map<unsigned,size_t> remap;
		std::vector<std::array<double,3>> nodes;
		std::vector<std::array<size_t,3>> tris;
		auto addV = [&](unsigned vid) -> size_t {
			auto jt = remap.find(vid);
			if (jt != remap.end()) return jt->second;
			size_t idx = nodes.size(); remap[vid] = idx; nodes.push_back(VertPos(sm, vid)); return idx;
		};
		for (const auto& f : prims.faces) tris.push_back({addV(f[0]), addV(f[1]), addV(f[2])});
		auto* mm = ps::registerSurfaceMesh("QEM Rejection Faces", nodes, tris);
		mm->setSurfaceColor(glm::vec3(1.0f, 1.0f, 0.0f));
		mm->setTransparency(0.35f);
		mm->setEnabled(true);
	} else if (ps::hasSurfaceMesh("QEM Rejection Faces")) {
		ps::getSurfaceMesh("QEM Rejection Faces")->setEnabled(false);
	}

	// ── Target vertex (green) ────────────────────────────────────────────────
	if (prims.targ_ver.has_value()) {
		std::vector<std::array<double,3>> pts = { *prims.targ_ver };
		auto* pc = ps::registerPointCloud("QEM Rejection Target", pts);
		pc->setPointColor(glm::vec3(0.0f, 1.0f, 0.0f));
		pc->setPointRadius(0.008, true);
		pc->setEnabled(true);
	} else if (ps::hasPointCloud("QEM Rejection Target")) {
		ps::getPointCloud("QEM Rejection Target")->setEnabled(false);
	}

	// ── Sliver triangle(s) at post-collapse positions — HardQualityCheckFailed ─
	if (!prims.tris_after.empty()) {
		std::vector<std::array<double,3>> nodes;
		std::vector<std::array<size_t,3>> tris;
		for (const auto& t : prims.tris_after) {
			size_t base = nodes.size();
			nodes.push_back(t[0]); nodes.push_back(t[1]); nodes.push_back(t[2]);
			tris.push_back({base, base + 1, base + 2});
		}
		auto* mm = ps::registerSurfaceMesh("QEM Sliver Face", nodes, tris);
		mm->setSurfaceColor(glm::vec3(1.0f, 1.0f, 0.0f));
		mm->setEnabled(true);
	} else if (ps::hasSurfaceMesh("QEM Sliver Face")) {
		ps::getSurfaceMesh("QEM Sliver Face")->setEnabled(false);
	}

	// ── Flipped face before(red)/after(orange) — NormalFlipped ───────────────
	if (prims.flipped_face.has_value()) {
		const auto& f = *prims.flipped_face;
		{
			std::vector<std::array<double,3>> nodes = { f[0][0], f[0][1], f[0][2] };
			std::vector<std::array<size_t,3>> tris  = { {0, 1, 2} };
			auto* mm = ps::registerSurfaceMesh("QEM Flipped Face Before", nodes, tris);
			mm->setSurfaceColor(glm::vec3(1.0f, 0.0f, 0.0f));
			mm->setEnabled(true);
		}
		{
			std::vector<std::array<double,3>> nodes = { f[1][0], f[1][1], f[1][2] };
			std::vector<std::array<size_t,3>> tris  = { {0, 1, 2} };
			auto* mm = ps::registerSurfaceMesh("QEM Flipped Face After", nodes, tris);
			mm->setSurfaceColor(glm::vec3(1.0f, 0.55f, 0.0f));
			mm->setEnabled(true);
		}
	} else {
		if (ps::hasSurfaceMesh("QEM Flipped Face Before")) ps::getSurfaceMesh("QEM Flipped Face Before")->setEnabled(false);
		if (ps::hasSurfaceMesh("QEM Flipped Face After"))  ps::getSurfaceMesh("QEM Flipped Face After")->setEnabled(false);
	}

	// ── Endpoint & target spheres (optional) ─────────────────────────────────
	if (show_spheres && eid < sm.edges.size() && sm.edges[eid].first) {
		std::vector<std::array<double,3>> pts;
		std::vector<double> radii;
		auto addSphere = [&](unsigned vid) {
			if (vid < sm.vertices.size() && sm.vertices[vid].first) {
				const auto& c = sm.vertices[vid].second->sphere.center;
				pts.push_back({c.X(), c.Y(), c.Z()});
				float r = sm.vertices[vid].second->sphere.radius;
				if (r < 1e-6f) r = 1e-4f;
				radii.push_back(r);
			}
		};
		addSphere(sm.edges[eid].second->vertices_.first);
		addSphere(sm.edges[eid].second->vertices_.second);
		if (prims.targ_ver.has_value()) {
			pts.push_back(*prims.targ_ver);
			radii.push_back(sm.edges[eid].second->sphere.radius);
		}
		if (!pts.empty()) {
			auto* pc = ps::registerPointCloud("QEM Rejection Spheres", pts);
			pc->setPointColor(glm::vec3(0.4f, 0.8f, 1.0f));
			pc->setTransparency(0.5f);
			pc->addScalarQuantity("radius", radii)->setEnabled(false);
			pc->setPointRadiusQuantity("radius", false);
			pc->setEnabled(true);
		}
	} else if (ps::hasPointCloud("QEM Rejection Spheres")) {
		ps::getPointCloud("QEM Rejection Spheres")->setEnabled(false);
	}
}

// ─────────────────────────────────────────────────────────────────────────────
bool HandlePick(const SlabMesh& sm, State& st, polyscope::Structure* struct_ptr, std::size_t local_idx)
{
	if (!ps::hasCurveNetwork("MAT QEM Rejection Edges")) return false;
	// Typed comparison (CurveNetwork* → Structure*) so the pointer is adjusted
	// correctly — matches how main_cli compares the QMAT rejection network.
	if (ps::getCurveNetwork("MAT QEM Rejection Edges") != struct_ptr)
		return false;

	// Flat pick index: [0, node_count) = node, [node_count, ...) = edge.
	if (local_idx >= st.node_count) {
		std::size_t edge_slot = local_idx - st.node_count;
		if (edge_slot < st.eid_order.size()) {
			unsigned eid = st.eid_order[edge_slot];
			if ((int)eid != st.selected_eid) {
				st.selected_eid = (int)eid;
				ClearPrimitives();
				ShowPrimitives(sm, eid, st.show_spheres);
			}
		}
	}
	return true;
}

// ─────────────────────────────────────────────────────────────────────────────
void DrawPanel(const SlabMesh& sm, State& st)
{
	ImGui::Separator();
	ImGui::Text("QEM Rejection (VCG path)");
	if (st.selected_eid < 0) {
		ImGui::TextDisabled("Click a 'MAT QEM Rejection Edges' edge to see why it was rejected");
		return;
	}

	unsigned eid = (unsigned)st.selected_eid;
	ImGui::Text("Rejection edge: %u", eid);

	auto rit = sm.qem_edge_last_rejection.find(eid);
	if (rit != sm.qem_edge_last_rejection.end())
		ImGui::Text("  Reason: %s", QemRejectionReasonName(rit->second));

	auto pit = sm.qem_edge_reason_primitives.find(eid);
	if (pit != sm.qem_edge_reason_primitives.end()) {
		const auto& p = pit->second;
		if (!p.message.empty()) ImGui::TextWrapped("  %s", p.message.c_str());
		for (const auto& mm : p.metrics)
			ImGui::Text("    %s = %.6f", mm.first.c_str(), mm.second);
	}

	if (ImGui::Checkbox("Show endpoint/target spheres##qem", &st.show_spheres))
		ShowPrimitives(sm, eid, st.show_spheres);
	if (ImGui::Button("Clear QEM rejection selection")) {
		st.selected_eid = -1;
		ClearPrimitives();
	}
}

} // namespace qemviz

#endif  // ONLY_USE_QEM_CONDITION_CHECKS && QMAT_WITH_POLYSCOPE
