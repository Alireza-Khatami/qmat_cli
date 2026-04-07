#include "SlabMesh.h"
#include <omp.h>
#include <unordered_map>
#include <map>
#include <functional>
#include <fstream>
#include <iomanip>

// Forward declarations — defined later in this file.
static std::vector<std::set<unsigned>> ComputeBpClusters(
	const std::set<unsigned>& bps,
	const std::vector<Mesh::Vertex_iterator>& vlist,
	unsigned n_mv,
	const std::set<std::array<int,2>>& sharp_edges,
	const std::set<std::array<int,2>>& concave_edges);
static SlabVertex::ClusterType ClusterTypeFromCount(unsigned n);
static SlabVertex::ClusterType ClusterTypeFromCountAndBplist(unsigned n_clusters, unsigned n_bplist);

// Compacts vertices/edges/faces by removing deleted slots and remapping all
// indices to a dense 0-based numbering.  Called by Export() before writing .ma.
// WARNING: invalidates any map keyed on edge/vertex/face IDs (e.g. edge_last_rejection).
// Always call ExportSkeletonPLY before Export().
void SlabMesh::AdjustStorage()
{
	std::vector<unsigned> newv;
	std::vector<unsigned> newe;
	std::vector<unsigned> newf;

	newv.resize(vertices.size());
	newe.resize(edges.size());
	newf.resize(faces.size());

	unsigned count = 0;
	for(unsigned i = 0; i < vertices.size(); i ++)
		if(vertices[i].first)
			newv[i] = count ++;

	count = 0;
	for(unsigned i = 0; i < edges.size(); i ++)
		if(edges[i].first)
			newe[i] = count ++;

	count = 0;
	for(unsigned i = 0; i < faces.size(); i ++)
		if(faces[i].first)
			newf[i] = count ++;

	std::vector<Bool_SlabVertexPointer> new_vertices;
	std::vector<Bool_SlabEdgePointer> new_edges;
	std::vector<Bool_SlabFacePointer> new_faces;

	for(unsigned i = 0; i < vertices.size(); i ++)
		if(vertices[i].first)
		{
			Bool_SlabVertexPointer bvp;
			bvp = vertices[i];
			std::set<unsigned> neweset;
			std::set<unsigned> newfset;

			for(std::set<unsigned>::iterator si = bvp.second->edges_.begin();
				si != bvp.second->edges_.end(); si ++)
				neweset.insert(newe[*si]);
			for(std::set<unsigned>::iterator si = bvp.second->faces_.begin();
				si != bvp.second->faces_.end(); si ++)
				newfset.insert(newf[*si]);

			bvp.second->edges_ = neweset;
			bvp.second->faces_ = newfset;
			new_vertices.push_back(bvp);
		}

		for(unsigned i = 0; i < edges.size(); i ++)
			if(edges[i].first)
			{
				Bool_SlabEdgePointer bep;
				bep = edges[i];
				std::pair<unsigned,unsigned> newvpair;
				std::set<unsigned> newfset;

				newvpair.first = newv[bep.second->vertices_.first];
				newvpair.second = newv[bep.second->vertices_.second];

				for(std::set<unsigned>::iterator si = bep.second->faces_.begin();
					si != bep.second->faces_.end(); si ++)
					newfset.insert(newf[*si]);

				bep.second->vertices_ = newvpair;
				bep.second->faces_ = newfset;
				new_edges.push_back(bep);
			}

			for(unsigned i = 0; i < faces.size(); i ++)
				if(faces[i].first)
				{
					Bool_SlabFacePointer bfp;
					bfp = faces[i];
					std::set<unsigned> newvset;
					std::set<unsigned> neweset;

					for(std::set<unsigned>::iterator si = bfp.second->vertices_.begin();
						si != bfp.second->vertices_.end(); si ++)
						newvset.insert(newv[*si]);

					for(std::set<unsigned>::iterator si = bfp.second->edges_.begin();
						si != bfp.second->edges_.end(); si ++)
						neweset.insert(newe[*si]);

					bfp.second->vertices_ = newvset;
					bfp.second->edges_ = neweset;
					new_faces.push_back(bfp);
				}

				vertices = new_vertices;
				edges = new_edges;
				faces = new_faces;			
}

bool SlabMesh::ValidVertex(unsigned vid){
	if(vid > vertices.size())
		return false;

	return vertices[vid].first;
}

bool SlabMesh::Edge(unsigned vid0, unsigned vid1, unsigned & eid)
{
	if(!ValidVertex(vid0) || !ValidVertex(vid1))
		return false;

	for(std::set<unsigned>::iterator si = (*vertices[vid0].second).edges_.begin(); si != (*vertices[vid0].second).edges_.end(); si ++)
	{
		if(edges[*si].first)
		{
			if(edges[*si].second->HasVertex(vid1))
			{
				eid = *si;
				return true;
			}
		}
	}

	return false;
}

bool SlabMesh::Face(const std::set<unsigned> & vset, unsigned & fid)
{
	if(vset.size() <= 0)
		return false;

	for(std::set<unsigned>::iterator si = vset.begin(); si != vset.end(); si ++)
		if(!ValidVertex(*si))
			return false;

	unsigned vid0 = *(vset.begin());
	for(std::set<unsigned>::iterator si = vertices[vid0].second->faces_.begin();
		si != vertices[vid0].second->faces_.end(); si ++)
		if(faces[*si].first)
		{
			if(faces[*si].second->vertices_ == vset)
			{
				fid = *si;
				return true;
			}
		}

		return false;
}

void SlabMesh::UpdateCentroid(unsigned fid)
{
	if(!faces[fid].first)
		return;

	faces[fid].second->centroid.center = Wm4::Vector3d::ZERO;
	faces[fid].second->centroid.radius = 0.0;

	unsigned count(0);
	for(std::set<unsigned>::iterator si = faces[fid].second->vertices_.begin(); si != faces[fid].second->vertices_.end(); si ++, count ++)
	{
		faces[fid].second->centroid.center += vertices[*si].second->sphere.center;
		faces[fid].second->centroid.radius += vertices[*si].second->sphere.radius;
	}
	faces[fid].second->centroid.center /= count;
	faces[fid].second->centroid.radius /= count;
}

void SlabMesh::ComputeFacesCentroid()
{
	for(unsigned i = 0; i < faces.size(); i ++)
		if(faces[i].first)
			UpdateCentroid(i);
}

void SlabMesh::UpdateNormal(unsigned fid)
{
	if(!faces[fid].first)
		return;

	Vector3d v[3];
	std::set<unsigned>::iterator si = faces[fid].second->vertices_.begin();
	v[0] = vertices[*si].second->sphere.center;
	si ++;
	v[1] = vertices[*si].second->sphere.center;
	si ++;
	v[2] = vertices[*si].second->sphere.center;
	faces[fid].second->normal = (v[1]-v[0]).Cross(v[2]-v[0]);
	faces[fid].second->normal.Normalize();
}

void SlabMesh::ComputeFacesNormal()
{
	for(unsigned i = 0; i < faces.size(); i ++)
		if(faces[i].first)
			UpdateNormal(i);
}

void SlabMesh::UpdateVertexNormal(unsigned vid)
{
	if(!vertices[vid].first)
		return;

	Vector3d v[3];
	std::set<unsigned> fs = vertices[vid].second->faces_;

	Vector3d vnormal;

	for (set<unsigned>::iterator si = fs.begin(); si != fs.end(); si++)
		vnormal += faces[*si].second->normal;

	vertices[vid].second->normal = vnormal;
	vertices[vid].second->normal.Normalize();
}

void SlabMesh::ComputeVerticesNormal()
{
	for(unsigned i = 0; i < vertices.size(); i ++)
		if(vertices[i].first)
			UpdateVertexNormal(i);
}

void SlabMesh::GetNeighborVertices(unsigned vid, std::set<unsigned> & neighborvertices)
{
	if(!vertices[vid].first)
		return;
	for(std::set<unsigned>::iterator si = vertices[vid].second->edges_.begin();
		si != vertices[vid].second->edges_.end(); si ++)
	{
		if(edges[*si].first)
		{
			if(edges[*si].second->vertices_.first == vid)
				neighborvertices.insert(edges[*si].second->vertices_.second);
			else
				neighborvertices.insert(edges[*si].second->vertices_.first);
		}
	}
}

void SlabMesh::GetLinkedEdges(unsigned eid, std::set<unsigned> & neighboredges)
{
	if(!edges[eid].first)
		return;
	neighboredges.clear();
	unsigned vid[2];
	vid[0] = edges[eid].second->vertices_.first;
	vid[1] = edges[eid].second->vertices_.second;
	for(unsigned k = 0; k < 2; k ++)
	{
		if(vertices[vid[k]].first)
		{
			for(std::set<unsigned>::iterator si = vertices[vid[k]].second->edges_.begin();
				si != vertices[vid[k]].second->edges_.end(); si ++)
				if(edges[*si].first)
					neighboredges.insert(*si);
		}
	}
	neighboredges.erase(eid);
}

void SlabMesh::GetAdjacentFaces(unsigned fid, std::set<unsigned> & neighborfaces)
{
	if(!faces[fid].first)
		return;

	neighborfaces.clear();
	for(std::set<unsigned>::iterator si = faces[fid].second->edges_.begin();
		si != faces[fid].second->edges_.end(); si ++)
	{
		if(edges[*si].first)
		{
			for(std::set<unsigned>::iterator si2 = edges[*si].second->faces_.begin();
				si2 != edges[*si].second->faces_.end(); si2 ++)
				neighborfaces.insert(*si2);
		}
	}
	neighborfaces.erase(fid);
}

bool SlabMesh::Contractible(unsigned vid_src, unsigned vid_tgt)
{
	if( !vertices[vid_src].first || !vertices[vid_tgt].first )
		return false;

	set<unsigned> ns;
	set<unsigned> nt;
	set<unsigned> ni;
	GetNeighborVertices(vid_src, ns);
	GetNeighborVertices(vid_tgt, nt);
	for(set<unsigned>::iterator si = ns.begin(); si != ns.end(); si ++)
		if(nt.find(*si) != nt.end())
			ni.insert(*si);
	unsigned eid;
	bool foundedge = Edge(vid_src, vid_tgt, eid);
	if(!foundedge)
		return false;
	if(edges[eid].second->faces_.size() != ni.size())
		return false;

	for(std::set<unsigned>::iterator si = vertices[vid_src].second->faces_.begin();
		si != vertices[vid_src].second->faces_.end(); si ++)
	{
		if(faces[*si].first)
		{
			if(!faces[*si].second->HasVertex(vid_tgt))
			{
				Vector3d pp[3], pa[3];
				unsigned count = 0;
				for(std::set<unsigned>::iterator si2 = faces[*si].second->vertices_.begin();
					si2 != faces[*si].second->vertices_.end(); si2 ++)
				{
					pp[count] = vertices[*si2].second->sphere.center;
					pa[count++] = (*si2 != vid_src)?vertices[*si2].second->sphere.center:vertices[vid_tgt].second->sphere.center;
				}
				Vector3d pnorm = TriangleNormal(pp[0],pp[1],pp[2]);
				Vector3d anorm = TriangleNormal(pa[0],pa[1],pa[2]);
				if(pnorm.Dot(anorm) < 0)
					return false;
			}
		}
	}
	return true;
}

bool SlabMesh::MergeVertices(unsigned vid_src1, unsigned vid_src2, unsigned &vid_tgt)
{
	if(vid_src1 == vid_src2)
		return false;

	unsigned eid;
	InsertVertex(new SlabVertex, vid_tgt);

	if (vertices[vid_src1].second->saved_vertex || vertices[vid_src2].second->saved_vertex)
		vertices[vid_tgt].second->saved_vertex = true;

	//if (vertices[vid_src1].second->fake_boundary_vertex || vertices[vid_src2].second->fake_boundary_vertex)
	//	vertices[vid_tgt].second->fake_boundary_vertex = true;

	//if (vertices[vid_src1].second->boundary_vertex || vertices[vid_src2].second->boundary_vertex)
	//	vertices[vid_tgt].second->boundary_vertex = true;

	// ── merge boundary-point and topology data ───────────────────────────────
	SlabVertex* sv1 = vertices[vid_src1].second;
	SlabVertex* sv2 = vertices[vid_src2].second;
	SlabVertex* svt = vertices[vid_tgt].second;

	// nmn_bplist: union of both sources
	svt->nmn_bplist = sv1->nmn_bplist;
	svt->nmn_bplist.insert(sv2->nmn_bplist.begin(), sv2->nmn_bplist.end());

	// Merge clusters by connectivity: treat each existing cluster as an atomic
	// node and union only clusters that share at least one cross mesh-edge.
	// Intra-cluster edges (already known) are never re-examined.
	if (pmesh)
	{
		// Take local copies before svt (which may alias sv1 or sv2) is mutated.
		const std::vector<std::set<unsigned>> cl1 = sv1->nmn_bplist_clusters;
		const std::vector<std::set<unsigned>> cl2 = sv2->nmn_bplist_clusters;
		const unsigned n1 = (unsigned)cl1.size();
		const unsigned n2 = (unsigned)cl2.size();
		const unsigned n_mv = (unsigned)pmesh->pVertexList.size();
		const auto& vlist = pmesh->pVertexList;

		// Union-find over cluster indices:
		// indices 0..n1-1  → sv1 clusters
		// indices n1..n1+n2-1 → sv2 clusters
		std::vector<unsigned> parent(n1 + n2);
		for (unsigned k = 0; k < n1 + n2; ++k) parent[k] = k;

		std::function<unsigned(unsigned)> find = [&](unsigned x) -> unsigned {
			if (parent[x] != x) parent[x] = find(parent[x]);
			return parent[x];
		};
		auto unite = [&](unsigned a, unsigned b) {
			a = find(a); b = find(b);
			if (a != b) parent[a] = b;
		};

		// For each (sv1-cluster, sv2-cluster) pair: unite if they share an
		// identical bp OR if any bp in cl1[i] has a surface mesh edge to any
		// bp in cl2[j].  Shared-bp check must come first so that a bp that
		// appears in both sources is not duplicated across two separate clusters.
		for (unsigned i = 0; i < n1; ++i)
		{
			for (unsigned j = 0; j < n2; ++j)
			{
				bool found = false;
				for (unsigned bp : cl1[i])
				{
					if (found) break;
					// Identity: same bp appears in both clusters.
					if (cl2[j].count(bp)) { found = true; break; }
					// Adjacency: surface mesh edge between different bps.
					if (bp >= n_mv) continue;
					auto circ = vlist[bp]->vertex_begin();
					auto done = circ;
					do {
						unsigned nbr = (unsigned)circ->opposite()->vertex()->id;
						if (cl2[j].count(nbr)) { found = true; break; }
					} while (++circ != done);
				}
				if (found) unite(i, n1 + j);
			}
		}

		// Collect the resulting merged clusters.
		std::unordered_map<unsigned, std::set<unsigned>> comp_map;
		for (unsigned i = 0; i < n1; ++i)
			for (unsigned bp : cl1[i])
				comp_map[find(i)].insert(bp);
		for (unsigned j = 0; j < n2; ++j)
			for (unsigned bp : cl2[j])
				comp_map[find(n1 + j)].insert(bp);

		svt->nmn_bplist_clusters.clear();
		svt->nmn_bplist_clusters.reserve(comp_map.size());
		for (auto& [root, members] : comp_map)
			svt->nmn_bplist_clusters.push_back(std::move(members));

		// Inherit cluster type from predecessors (both are same type, enforced by CanMerge).
		// For MS_* types the old recomputation from bp count would produce a wrong T0-T5 value.
		svt->nmn_cluster_type = sv1->nmn_cluster_type;
		// [OLD METHOD] recompute from bp cluster count — use for Voronoi/DT path:
		// svt->nmn_cluster_type = ClusterTypeFromCountAndBplist(
		//     (unsigned)svt->nmn_bplist_clusters.size(),
		//     (unsigned)svt->nmn_bplist.size());
	}

	// merged vertex is never steep — it now represents a larger surface region
	svt->is_spike = false;

	// topology flags: conservative union
	// OR for seam/junction/boundary; AND for sheet (only if both were sheets
	// and no seam/boundary edges are introduced)
	svt->topo_is_seam     = sv1->topo_is_seam     || sv2->topo_is_seam;
	svt->topo_is_junction = sv1->topo_is_junction  || sv2->topo_is_junction;
	svt->topo_is_boundary = sv1->topo_is_boundary  || sv2->topo_is_boundary;
	svt->topo_is_sheet    = sv1->topo_is_sheet     && sv2->topo_is_sheet
	                        && !svt->topo_is_seam  && !svt->topo_is_boundary;

	std::vector< std::set<unsigned> > tri_vec;
	for(std::set<unsigned>::iterator si = vertices[vid_src1].second->faces_.begin();
		si != vertices[vid_src1].second->faces_.end(); si ++)
		if(!faces[*si].second->HasVertex(vid_tgt))
		{
			std::set<unsigned> vset = faces[*si].second->vertices_;
			vset.erase(vid_src1);
			vset.insert(vid_tgt);
			tri_vec.push_back(vset);
		}

		for(std::set<unsigned>::iterator si = vertices[vid_src2].second->faces_.begin();
			si != vertices[vid_src2].second->faces_.end(); si ++)
			if(!faces[*si].second->HasVertex(vid_tgt))
			{
				std::set<unsigned> vset = faces[*si].second->vertices_;
				vset.erase(vid_src2);
				vset.insert(vid_tgt);
				tri_vec.push_back(vset);
			}

			std::vector< std::pair<unsigned,unsigned> > edge_vec;
			for(std::set<unsigned>::iterator si = vertices[vid_src1].second->edges_.begin();
				si != vertices[vid_src1].second->edges_.end(); si ++)
				if(!edges[*si].second->HasVertex(vid_tgt))
				{
					std::pair<unsigned, unsigned> vp = edges[*si].second->vertices_;
					if(vp.first == vid_src1)
						vp.first = vid_tgt;
					if(vp.second == vid_src1)
						vp.second = vid_tgt;
					edge_vec.push_back(vp);
				}

				for(std::set<unsigned>::iterator si = vertices[vid_src2].second->edges_.begin();
					si != vertices[vid_src2].second->edges_.end(); si ++)
					if(!edges[*si].second->HasVertex(vid_tgt))
					{
						std::pair<unsigned, unsigned> vp = edges[*si].second->vertices_;
						if(vp.first == vid_src2)
							vp.first = vid_tgt;
						if(vp.second == vid_src2)
							vp.second = vid_tgt;
						edge_vec.push_back(vp);
					}

					DeleteVertex(vid_src1);
					DeleteVertex(vid_src2);

					for(unsigned i = 0; i < tri_vec.size(); i ++)
						InsertFace(tri_vec[i]);

					for(unsigned i = 0; i < edge_vec.size(); i ++)
					{
						unsigned neweid;
						InsertEdge(edge_vec[i].first, edge_vec[i].second, neweid);
					}

					return true;

}

unsigned SlabMesh::VertexIncidentEdgeCount(unsigned vid)
{
	if(!vertices[vid].first)
		return 0;
	return (unsigned)vertices[vid].second->edges_.size();
}

unsigned SlabMesh::VertexIncidentFaceCount(unsigned vid)
{
	if(!vertices[vid].first)
		return 0;
	return (unsigned)vertices[vid].second->faces_.size();
}

unsigned SlabMesh::EdgeIncidentFaceCount(unsigned eid)
{
	if(!edges[eid].first)
		return 0;
	return (unsigned)edges[eid].second->faces_.size();
}

void SlabMesh::DeleteFace(unsigned fid)
{
	if(!faces[fid].first)
		return;

	for(std::set<unsigned>::iterator si = faces[fid].second->vertices_.begin();
		si != faces[fid].second->vertices_.end(); si ++)
		vertices[*si].second->faces_.erase(fid);

	for(std::set<unsigned>::iterator si = faces[fid].second->edges_.begin();
		si != faces[fid].second->edges_.end(); si ++)
		edges[*si].second->faces_.erase(fid);

	delete faces[fid].second;
	faces[fid].first = false;
	numFaces --;
}

void SlabMesh::DeleteEdge(unsigned eid)
{
	if(!edges[eid].first)
		return;

	// �����boundary_edge����������Եĸ���
	if (edges[eid].second->fake_boundary_edge)
	{
		if(vertices[edges[eid].second->vertices_.first].first)
		{
			vertices[edges[eid].second->vertices_.first].second->boundary_edge_vec.erase(eid);
			vertices[edges[eid].second->vertices_.first].second->fake_boundary_vertex = 
				vertices[edges[eid].second->vertices_.first].second->boundary_edge_vec.size() > 0 ? true : false;
		}
		if(vertices[edges[eid].second->vertices_.second].first)
		{	
			vertices[edges[eid].second->vertices_.second].second->boundary_edge_vec.erase(eid);
			vertices[edges[eid].second->vertices_.second].second->fake_boundary_vertex = 
				vertices[edges[eid].second->vertices_.second].second->boundary_edge_vec.size() > 0 ? true : false;
		}
	}

	if(vertices[edges[eid].second->vertices_.first].first)
		vertices[edges[eid].second->vertices_.first].second->edges_.erase(eid);
	if(vertices[edges[eid].second->vertices_.second].first)
		vertices[edges[eid].second->vertices_.second].second->edges_.erase(eid);
	std::set<unsigned> faces_del;
	for(std::set<unsigned>::iterator si = edges[eid].second->faces_.begin();
		si != edges[eid].second->faces_.end(); si ++)
		faces_del.insert(*si);
	for(std::set<unsigned>::iterator si = faces_del.begin(); si != faces_del.end(); si ++)
		DeleteFace(*si);

	delete edges[eid].second;
	edges[eid].first = false;
	numEdges --;
}

void SlabMesh::DeleteVertex(unsigned vid)
{
	if(!vertices[vid].first)
		return;

	std::set<unsigned> edges_del;
	for(std::set<unsigned>::iterator si = vertices[vid].second->edges_.begin();
		si != vertices[vid].second->edges_.end(); si ++)
		edges_del.insert(*si);

	std::set<unsigned> faces_del;
	for(std::set<unsigned>::iterator si = vertices[vid].second->faces_.begin();
		si != vertices[vid].second->faces_.end(); si ++)
		faces_del.insert(*si);

	for(std::set<unsigned>::iterator si = edges_del.begin(); si != edges_del.end(); si ++)
		DeleteEdge(*si);

	for(std::set<unsigned>::iterator si = faces_del.begin(); si != faces_del.end(); si ++)
		DeleteFace(*si);

	delete vertices[vid].second;
	vertices[vid].first = false;
	numVertices --;
}

void SlabMesh::InsertVertex(SlabVertex *vertex, unsigned &vid){
	Bool_SlabVertexPointer bvp;
	bvp.first = true;
	bvp.second = vertex;
	vid = (unsigned)vertices.size();
	bvp.second->index = vid;
	vertices.push_back(bvp);
	numVertices ++;

}

void SlabMesh::InsertEdge(unsigned vid0, unsigned vid1, unsigned & eid)
{
	if(Edge(vid0,vid1,eid))
		return;
	if (!ValidVertex(vid0) || !ValidVertex(vid1))
	{
		return;
	}
	Bool_SlabEdgePointer bep;
	bep.first = true;
	bep.second = new SlabEdge;
	bep.second->vertices_.first = vid0;
	bep.second->vertices_.second = vid1;
	vertices[vid0].second->edges_.insert((unsigned)edges.size());
	vertices[vid1].second->edges_.insert((unsigned)edges.size());
	eid = (unsigned)edges.size();
	bep.second->index = eid;
	edges.push_back(bep);
	ComputeEdgeCone(eid);
	numEdges ++;
}

void SlabMesh::InsertFace(std::set<unsigned> vset)
{
	unsigned fid;
	if(Face(vset,fid))
		return;

	unsigned vid[3];
	std::set<unsigned>::iterator si = vset.begin();
	vid[0] = *si;
	si ++;
	vid[1] = *si;
	si ++;
	vid[2] = *si;

	if(!vertices[vid[0]].first || !vertices[vid[1]].first || !vertices[vid[2]].first)
		return;

	for(std::set<unsigned>::iterator si = vertices[vid[0]].second->faces_.begin(); 
		si != vertices[vid[0]].second->faces_.end(); si ++)
		if(faces[*si].second->vertices_ == vset) // duplicate
			return;

	unsigned eid[3];
	InsertEdge(vid[0],vid[1],eid[0]);
	InsertEdge(vid[0],vid[2],eid[1]);
	InsertEdge(vid[1],vid[2],eid[2]);

	Bool_SlabFacePointer bfp;
	bfp.first = true;
	bfp.second = new SlabFace;
	bfp.second->vertices_.insert(vid[0]);
	bfp.second->vertices_.insert(vid[1]);
	bfp.second->vertices_.insert(vid[2]);
	bfp.second->edges_.insert(eid[0]);
	bfp.second->edges_.insert(eid[1]);
	bfp.second->edges_.insert(eid[2]);
	bfp.second->index = faces.size();
	vertices[vid[0]].second->faces_.insert(faces.size());
	vertices[vid[1]].second->faces_.insert(faces.size());
	vertices[vid[2]].second->faces_.insert(faces.size());
	edges[eid[0]].second->faces_.insert(faces.size());
	edges[eid[1]].second->faces_.insert(faces.size());
	edges[eid[2]].second->faces_.insert(faces.size());
	faces.push_back(bfp);
	UpdateCentroid((unsigned)faces.size()-1);
	UpdateNormal((unsigned)faces.size()-1);
	ComputeFaceSimpleTriangles((unsigned)faces.size() - 1);
	numFaces ++;
}

void SlabMesh::CleanIsolatedVertices()
{
	for(unsigned i = 0; i < vertices.size(); i ++)
		if(vertices[i].first)
		{
			if( (vertices[i].second->edges_.size() == 0) && (vertices[i].second->faces_.size() == 0) )
			{
				delete vertices[i].second;
				vertices[i].first = false;
				numVertices --;
			}
		}
}

void SlabMesh::ComputeVertexProperty(unsigned vid)
{
	if(!vertices[vid].first)
		return;
	std::set<unsigned> neighbor_vertices;
	GetNeighborVertices(vid,neighbor_vertices);
	for(std::set<unsigned>::iterator si = neighbor_vertices.begin();
		si != neighbor_vertices.end(); si ++)
		vertices[*si].second->tag = 0;

	for(std::set<unsigned>::iterator si = vertices[vid].second->faces_.begin();
		si != vertices[vid].second->faces_.end(); si ++)
		for(std::set<unsigned>::iterator ssi = faces[*si].second->vertices_.begin();
			ssi != faces[*si].second->vertices_.end(); ssi ++)
			vertices[*ssi].second->tag ++;

	bool has_one_tag(false);
	bool has_two_plus_tag(false);

	for(std::set<unsigned>::iterator si = neighbor_vertices.begin(); si != neighbor_vertices.end(); si ++)
		if( vertices[*si].second->tag > 2)
			has_two_plus_tag = true; // non-manifold
		else if( vertices[*si].second->tag < 2)
			has_one_tag = true;

	vertices[vid].second->is_boundary = has_one_tag?true:false;
	vertices[vid].second->is_disk = (has_one_tag||has_two_plus_tag)?false:true;
	vertices[vid].second->is_non_manifold = has_two_plus_tag?true:false;
}

void SlabMesh::ComputeVerticesProperty()
{
	for(unsigned i = 0; i < vertices.size(); i ++)
		if(vertices[i].first)
			ComputeVertexProperty(i);
}

void SlabMesh::ComputeEdgeCone(unsigned eid)
{
	if(!edges[eid].first)
		return;

	// test validation
	Vector3d c0 = vertices[edges[eid].second->vertices_.first].second->sphere.center;
	Vector3d c1 = vertices[edges[eid].second->vertices_.second].second->sphere.center;
	double r0 = vertices[edges[eid].second->vertices_.first].second->sphere.radius;
	double r1 = vertices[edges[eid].second->vertices_.second].second->sphere.radius;
	Vector3d c0c1 = c1-c0;
	double templeng = c0c1.Length() - abs(r1 - r0);


	Cone newc(vertices[edges[eid].second->vertices_.first].second->sphere.center, vertices[edges[eid].second->vertices_.first].second->sphere.radius,
		vertices[edges[eid].second->vertices_.second].second->sphere.center, vertices[edges[eid].second->vertices_.second].second->sphere.radius);
	edges[eid].second->cone = newc;
	if(newc.type == 1)
		edges[eid].second->valid_cone = false;
	else
		edges[eid].second->valid_cone = true;
}

void SlabMesh::ComputeEdgesCone()
{
	for(unsigned i = 0; i < edges.size(); i ++)
		if(edges[i].first)
			ComputeEdgeCone(i);
}

void SlabMesh::ComputeFaceSimpleTriangles(unsigned fid)
{
	if(!faces[fid].first)
		return;
	SimpleTriangle st0,st1;
	Wm4::Vector3d pos[3];
	double radius[3];
	unsigned count = 0;
	for(std::set<unsigned>::iterator si = faces[fid].second->vertices_.begin();
		si != faces[fid].second->vertices_.end(); si ++, count ++)
	{
		pos[count] = vertices[*si].second->sphere.center;
		radius[count] = vertices[*si].second->sphere.radius;
	}
	if(TriangleFromThreeSpheres(pos[0],radius[0],pos[1],radius[1],pos[2],radius[2],st0,st1))
	{
		faces[fid].second->st[0] = st0;
		faces[fid].second->st[1] = st1;
		faces[fid].second->valid_st = true;
	}
	else
		faces[fid].second->valid_st = false;

}

void SlabMesh::ComputeFacesSimpleTriangles()
{
	for(unsigned i = 0; i < faces.size(); i ++)
		if(faces[i].first)
			ComputeFaceSimpleTriangles(i);
}

void SlabMesh::DistinguishVertexType()
{
	for (unsigned i = 0; i != edges.size(); i++)
	{
		if (edges[i].first && (edges[i].second->faces_.size() == 1 || edges[i].second->faces_.size() == 0))
			//if (edges[i].first && edges[i].second->faces_.size() == 1)
		{
			// fake boundary edge and fake boundary vertex
			edges[i].second->fake_boundary_edge = true;
			vertices[edges[i].second->vertices_.first].second->fake_boundary_vertex = true;
			vertices[edges[i].second->vertices_.second].second->fake_boundary_vertex = true;
			vertices[edges[i].second->vertices_.first].second->boundary_edge_vec.insert(i);
			vertices[edges[i].second->vertices_.second].second->boundary_edge_vec.insert(i);
			//boundary_vertexes.insert(edges[i].second->vertices_.first);
			//boundary_vertexes.insert(edges[i].second->vertices_.second);
		}
		else if (edges[i].first && edges[i].second->faces_.size() >= 3)
		{
			// non manifold edge
			edges[i].second->non_manifold_edge = true;
			vertices[edges[i].second->vertices_.first].second->non_manifold_vertex = true;
			vertices[edges[i].second->vertices_.second].second->non_manifold_vertex = true;
			//vertices[edges[i].second->vertices_.first].second->collaspe_weight += edges[i].second->faces_.size() - 2;
			//vertices[edges[i].second->vertices_.second].second->collaspe_weight += edges[i].second->faces_.size() - 2;
		}
	}

#if 0
	// real boundary edge and vertex, saved vertex
	for (unsigned i = 0; i != vertices.size(); i++)
	{
		if (vertices[i].first && vertices[i].second->fake_boundary_vertex)
		{
			//vertices[i].second->collaspe_weight += 1.0;

			set<unsigned>::iterator si = vertices[i].second->edges_.begin();
			vector<unsigned> fake_boundary_edge_number;
			for (; si != vertices[i].second->edges_.end(); si++)
				if (edges[*si].second->fake_boundary_edge)
					fake_boundary_edge_number.push_back(*si);

			if (fake_boundary_edge_number.size() >= 2)
				vertices[i].second->boundary_vertex = true;

			set<unsigned> face_set;
			unsigned face_number = 0;
			for (unsigned j = 0; j < fake_boundary_edge_number.size(); j++)
			{
				edges[fake_boundary_edge_number[j]].second->boundary_edge = true;

				face_number += edges[fake_boundary_edge_number[j]].second->faces_.size();
				face_set.insert(edges[fake_boundary_edge_number[j]].second->faces_.begin(), 
					edges[fake_boundary_edge_number[j]].second->faces_.end());
			}

			if (face_number != face_set.size())
			{
				vertices[i].second->collaspe_weight += 10.0 * edges[fake_boundary_edge_number[0]].second->sphere.center.Dot(edges[fake_boundary_edge_number[1]].second->sphere.center) 
					/ edges[fake_boundary_edge_number[0]].second->sphere.center.Length() / edges[fake_boundary_edge_number[1]].second->sphere.center.Length();

				vertices[i].second->saved_vertex = true;
			}
		}
	}

	// judge the bounding points of each boundary point
	for(std::set<unsigned>::iterator sit = boundary_vertexes.begin(); sit != boundary_vertexes.end(); sit++)
		for(std::set<unsigned>::iterator sit2 = boundary_vertexes.begin(); sit2 != boundary_vertexes.end(); sit2++)
		{
			if (*sit == *sit2)
				continue;
			Vector3d tempvec = vertices[*sit].second->sphere.center - vertices[*sit2].second->sphere.center;

			//if (tempvec.Length() <= vertices[*sit].second->sphere.radius) 
			if (tempvec.Length() <= vertices[*sit].second->sphere.radius) 
			{
				vertices[*sit].second->bound_point_vec.push_back(*sit2);
				vertices[*sit].second->boundVec += tempvec;
			}
		}

		for(std::set<unsigned>::iterator sit = boundary_vertexes.begin(); sit != boundary_vertexes.end(); sit++)
		{
			bound_vector.push_back(pair<unsigned, unsigned>(*sit, vertices[*sit].second->bound_point_vec.size()));
			if (vertices[*sit].second->bound_point_vec.size() == 0)
				boundVec_vector.push_back(pair<unsigned, double>(*sit, 0));
			else
				boundVec_vector.push_back(pair<unsigned, double>(*sit, vertices[*sit].second->boundVec.Length() / vertices[*sit].second->bound_point_vec.size()));
			//boundVec_vector.push_back(pair<unsigned, unsigned>(*sit, vertices[*sit].second->boundVec.Length()));
		}

		sort(bound_vector.begin(), bound_vector.end(), cmpByValue<unsigned>());
		sort(boundVec_vector.begin(), boundVec_vector.end(), cmpByBiggerValue<double>());

#endif
}

unsigned SlabMesh::GetSavedPointNumber()
{
	//unsigned count = vertices.size();
	//for (unsigned i = 0; i < count; i++)
	//{
	//	if (vertices[i].first)
	//	{
	//		if (vertices[i].second->edges_.size() == 1 && vertices[i].second->faces_.size() == 0)
	//		{
	//			if (vertices[i].second->sphere.radius > 1.0e-3)
	//			{
	//				InsertSavedPoint(vertices[i].second->index);
	//			}

	//		}
	//	}
	//}
	//return count;

	unsigned count = edges.size();
	for (unsigned i = 0; i < count; i++)
	{
		if (edges[i].first)
		{
			if (edges[i].second->faces_.size() == 0)
			{
				if (vertices[edges[i].second->vertices_.first].second->edges_.size() == 1)
				{
					InsertSavedPoint(edges[i].second->vertices_.first);
					//if (vertices[edges[i].second->vertices_.first].second->sphere.radius > 1.0e-3)
					//	InsertSavedPoint(edges[i].second->vertices_.first);
					//else if (vertices[edges[i].second->vertices_.second].second->sphere.radius > 1.0e-3)
					//	InsertSavedPoint(edges[i].second->vertices_.second);
				}else if (vertices[edges[i].second->vertices_.second].second->edges_.size() == 1)
				{
					InsertSavedPoint(edges[i].second->vertices_.second);
					//if (vertices[edges[i].second->vertices_.second].second->sphere.radius > 1.0e-3)
					//	InsertSavedPoint(edges[i].second->vertices_.second);
					//else if (vertices[edges[i].second->vertices_.first].second->sphere.radius > 1.0e-3)
					//	InsertSavedPoint(edges[i].second->vertices_.first);
				}
			}
		}
	}
	return count;
}

unsigned SlabMesh::GetConnectPointNumber()
{
	unsigned count = 0;
	for (unsigned i = 0; i < edges.size(); i++)
	{
		if (edges[i].first && edges[i].second->faces_.size() == 0)
		{
			if (vertices[edges[i].second->vertices_.first].second->edges_.size() <= 3)
				count++;

			if (vertices[edges[i].second->vertices_.second].second->edges_.size() <= 3)
				count++;
		}
	}
	return count;
}

void SlabMesh::InsertSavedPoint(unsigned vid)
{
	if (vertices[vid].second->saved_vertex)
		return;

	unsigned vid_tgt;
	InsertVertex(new SlabVertex, vid_tgt);

	std::vector< std::set<unsigned> > tri_vec;
	for(std::set<unsigned>::iterator si = vertices[vid].second->faces_.begin();
		si != vertices[vid].second->faces_.end(); si ++)
		if(!faces[*si].second->HasVertex(vid_tgt))
		{
			std::set<unsigned> vset = faces[*si].second->vertices_;
			vset.erase(vid);
			vset.insert(vid_tgt);
			tri_vec.push_back(vset);
		}

		std::vector< std::pair<unsigned,unsigned> > edge_vec;
		for(std::set<unsigned>::iterator si = vertices[vid].second->edges_.begin();
			si != vertices[vid].second->edges_.end(); si ++)
			if(!edges[*si].second->HasVertex(vid_tgt))
			{
				std::pair<unsigned, unsigned> vp = edges[*si].second->vertices_;
				if(vp.first == vid)
					vp.first = vid_tgt;
				if(vp.second == vid)
					vp.second = vid_tgt;
				edge_vec.push_back(vp);
			}

			vertices[vid_tgt].second->saved_vertex = true;
			vertices[vid_tgt].second->sphere = vertices[vid].second->sphere;
			vertices[vid_tgt].second->bplist          = vertices[vid].second->bplist;
			vertices[vid_tgt].second->nmn_bplist      = vertices[vid].second->nmn_bplist;
			vertices[vid_tgt].second->is_spike = vertices[vid].second->is_spike;
			vertices[vid_tgt].second->topo_is_sheet    = vertices[vid].second->topo_is_sheet;
			vertices[vid_tgt].second->topo_is_seam     = vertices[vid].second->topo_is_seam;
			vertices[vid_tgt].second->topo_is_junction = vertices[vid].second->topo_is_junction;
			vertices[vid_tgt].second->topo_is_boundary = vertices[vid].second->topo_is_boundary;

			vertices[vid_tgt].second->slab_A = vertices[vid].second->slab_A;
			vertices[vid_tgt].second->slab_b = vertices[vid].second->slab_b;
			vertices[vid_tgt].second->slab_c = vertices[vid].second->slab_c;

			vertices[vid_tgt].second->mean_square_error = vertices[vid].second->mean_square_error;
			vertices[vid_tgt].second->related_face = vertices[vid].second->related_face;

			DeleteVertex(vid);

			for(unsigned i = 0; i < tri_vec.size(); i ++)
				InsertFace(tri_vec[i]);

			for(unsigned i = 0; i < edge_vec.size(); i ++)
			{
				unsigned neweid;
				InsertEdge(edge_vec[i].first, edge_vec[i].second, neweid);
			}

			for (std::set<unsigned>::iterator si = vertices[vid_tgt].second->edges_.begin(); si != vertices[vid_tgt].second->edges_.end(); si ++)
			{
				EvaluateEdgeCollapseCost(*si);
				if (edges[*si].second->collapse_cost != DBL_MAX)
				{
					edge_collapses_queue.push(EdgeInfo(*si, edges[*si].second->collapse_cost));
					ComputeEdgeCone(*si);
				}
			}

			return;
}

// �ж��Ƿ����������η�ת���
bool SlabMesh::Contractible(unsigned vid_src1, unsigned vid_src2, const Vector3d &v_tgt)
{
	if( !vertices[vid_src1].first || !vertices[vid_src2].first )
		return false;

	set<unsigned> fs1;
	set<unsigned> fs2;
	fs1 = vertices[vid_src1].second->faces_;
	fs2 = vertices[vid_src2].second->faces_;

	for (std::set<unsigned>::iterator si = fs1.begin(); si != fs1.end(); si++)
	{
		if(faces[*si].first)
		{
			if ( faces[*si].second->vertices_.find(vid_src1) != faces[*si].second->vertices_.end() 
				&& faces[*si].second->vertices_.find(vid_src2) != faces[*si].second->vertices_.end())
				continue;

			Vector3d pp[3], pa[3];
			unsigned count = 0;
			for(std::set<unsigned>::iterator si2 = faces[*si].second->vertices_.begin();
				si2 != faces[*si].second->vertices_.end(); si2 ++)
			{
				pp[count] = vertices[*si2].second->sphere.center;
				pa[count++] = (*si2 != vid_src1) ? vertices[*si2].second->sphere.center : v_tgt;
			}
			Vector3d pnorm = TriangleNormal(pp[0],pp[1],pp[2]);
			Vector3d anorm = TriangleNormal(pa[0],pa[1],pa[2]);

			//double angle = VectorAngle(pnorm, anorm);
			//if (angle > Wm4::Math<double>::PI * 2.0 / 3.0)
			//	return false;

			if(pnorm.Dot(anorm) < 0)
				return false;
		}
	}

	for (std::set<unsigned>::iterator si = fs2.begin(); si != fs2.end(); si++)
	{
		if(faces[*si].first)
		{
			if ( faces[*si].second->vertices_.find(vid_src1) != faces[*si].second->vertices_.end() 
				&& faces[*si].second->vertices_.find(vid_src2) != faces[*si].second->vertices_.end())
				continue;

			Vector3d pp[3], pa[3];
			unsigned count = 0;
			for(std::set<unsigned>::iterator si2 = faces[*si].second->vertices_.begin();
				si2 != faces[*si].second->vertices_.end(); si2 ++)
			{
				pp[count] = vertices[*si2].second->sphere.center;
				pa[count++] = (*si2 != vid_src2) ? vertices[*si2].second->sphere.center : v_tgt;
			}
			Vector3d pnorm = TriangleNormal(pp[0],pp[1],pp[2]);
			Vector3d anorm = TriangleNormal(pa[0],pa[1],pa[2]);

			//double angle = VectorAngle(pnorm, anorm);
			//if (angle > Wm4::Math<double>::PI * 2.0 / 3.0)
			//	return false;

			if(pnorm.Dot(anorm) < 0)
				return false;
		}
	}

	return true;
}

bool SlabMesh::MinCostBoundaryEdgeCollapse(unsigned & eid)
{
	//merge 2 vertices of the edge first, then move the combined vertex to the preferred point and resize it.
	unsigned v1, v2;
	v1 = edges[eid].second->vertices_.first;
	v2 = edges[eid].second->vertices_.second;

	{
		RejectionReason reason;
		ReasonPrimitives prims;
		if (!CanMerge(v1, v2, &reason, &prims))
		{ LogCollapseRejection("boundary", eid, v1, v2, edges[eid].second->collapse_cost, reason, std::move(prims)); return false; }
	}

	Wm4::Matrix4d A = edges[eid].second->slab_A;
	Wm4::Vector4d b = edges[eid].second->slab_b;
	double c = edges[eid].second->slab_c;
	Sphere sphere = edges[eid].second->sphere;

	if (prevent_inversion == true)
	{
		if (!Contractible(v1, v2, sphere.center))
		{
			ReasonPrimitives prims; prims.vertices = { v1, v2 }; prims.edges = { {v1, v2} };
			LogCollapseRejection("boundary", eid, v1, v2, edges[eid].second->collapse_cost,
			                     RejectionReason::InversionWouldOccur, std::move(prims));
			return false;
		}
	}


	set<unsigned> temp_bplist;
	for (set<unsigned>::iterator it = vertices[v1].second->bplist.begin(); it != vertices[v1].second->bplist.end(); it++)
		temp_bplist.insert(*it);
	for (set<unsigned>::iterator it = vertices[v2].second->bplist.begin(); it != vertices[v2].second->bplist.end(); it++)
		temp_bplist.insert(*it);	

	unsigned temp_related_face = vertices[v1].second->related_face + vertices[v2].second->related_face;
	double temp_mean_squre_error = edges[eid].second->collapse_cost < 0 ? 0 : edges[eid].second->collapse_cost / temp_related_face;
	max_mean_squre_error = max(temp_mean_squre_error, max_mean_squre_error);

#ifdef QMAT_WITH_POLYSCOPE
	if (on_collapse_cb) {
		on_collapse_cb(
			v1, vertices[v1].second->sphere.center, vertices[v1].second->sphere.radius,
			v2, vertices[v2].second->sphere.center, vertices[v2].second->sphere.radius,
			sphere);
	}
#endif

	// Capture positions and bplists before merge for history recording.
	const auto& c1 = vertices[v1].second->sphere.center;
	const auto& c2 = vertices[v2].second->sphere.center;
	std::array<double,3> hist_pos1 = {c1.X(), c1.Y(), c1.Z()};
	std::array<double,3> hist_pos2 = {c2.X(), c2.Y(), c2.Z()};
	std::vector<unsigned> hist_bp1(vertices[v1].second->nmn_bplist.begin(),
	                               vertices[v1].second->nmn_bplist.end());
	std::vector<unsigned> hist_bp2(vertices[v2].second->nmn_bplist.begin(),
	                               vertices[v2].second->nmn_bplist.end());

	unsigned former_edge_number = edges.size();
	unsigned vid_tgt;
	if(MergeVertices(v1, v2, vid_tgt)){
		vertices[vid_tgt].second->slab_A = A;
		vertices[vid_tgt].second->slab_b = b;
		vertices[vid_tgt].second->slab_c = c;
		vertices[vid_tgt].second->sphere = sphere;
		vertices[vid_tgt].second->related_face = temp_related_face;
		vertices[vid_tgt].second->mean_square_error = temp_mean_squre_error;
		vertices[vid_tgt].second->bplist = temp_bplist;

		// Record collapse in history (bplist_after and clusters_after captured post-merge).
		{
			std::vector<unsigned> hist_bpt(vertices[vid_tgt].second->nmn_bplist.begin(),
			                               vertices[vid_tgt].second->nmn_bplist.end());
			std::vector<std::vector<unsigned>> hist_clusters;
			for (const auto& cl : vertices[vid_tgt].second->nmn_bplist_clusters)
				hist_clusters.push_back(std::vector<unsigned>(cl.begin(), cl.end()));
			const unsigned step = history.TotalSteps();
			history.Record(step, v1, v2, vid_tgt,
			               hist_pos1, hist_pos2,
			               std::move(hist_bp1), std::move(hist_bp2),
			               std::move(hist_bpt), std::move(hist_clusters));
			if (history.keyframe_interval > 0 &&
			    (int)history.TotalSteps() % history.keyframe_interval == 0)
				history.TakeKeyframe(history.TotalSteps(), *this);
		}

		// Refresh topology for the new merged vertex.
		RecomputeVertexTopology(vid_tgt);

		switch(boundary_compute_scale)
		{
		case 1:			
			// �������ӵı��ж�����
			for (unsigned i = former_edge_number; i < edges.size(); i++)
			{
				if (edges[i].second->faces_.size() <= 1)
				{
					edges[i].second->fake_boundary_edge = true;
					vertices[edges[i].second->vertices_.first].second->fake_boundary_vertex = true;
					vertices[edges[i].second->vertices_.second].second->fake_boundary_vertex = true;
					vertices[edges[i].second->vertices_.first].second->boundary_edge_vec.insert(i);
					vertices[edges[i].second->vertices_.second].second->boundary_edge_vec.insert(i);
				}
			}
			//for (std::set<unsigned>::iterator si = vertices[vid_tgt].second->edges_.begin(); si != vertices[vid_tgt].second->edges_.end(); si ++)
			//{
			//	if (edges[*si].second->faces_.size() <= 1)
			//	{
			//		unsigned fir = edges[*si].second->vertices_.first;
			//		unsigned sec = edges[*si].second->vertices_.second;
			//		vertices[fir].second->fake_boundary_vertex = true;
			//		vertices[sec].second->fake_boundary_vertex = true;
			//	}
			//}
			break;
		case 2:
			// �������ӵı��ж�����
			for (unsigned i = former_edge_number; i < edges.size(); i++)
			{
				if (edges[i].second->faces_.size() <= 1)
				{
					edges[i].second->fake_boundary_edge = true;
					vertices[edges[i].second->vertices_.first].second->fake_boundary_vertex = true;
					vertices[edges[i].second->vertices_.second].second->fake_boundary_vertex = true;
					vertices[edges[i].second->vertices_.first].second->boundary_edge_vec.insert(i);
					vertices[edges[i].second->vertices_.second].second->boundary_edge_vec.insert(i);
				}
			}
			//for (std::set<unsigned>::iterator si = vertices[vid_tgt].second->edges_.begin(); si != vertices[vid_tgt].second->edges_.end(); si ++)
			//{
			//	if (edges[*si].second->faces_.size() <= 1)
			//	{
			//		unsigned fir = edges[*si].second->vertices_.first;
			//		unsigned sec = edges[*si].second->vertices_.second;
			//		vertices[fir].second->fake_boundary_vertex = true;
			//		vertices[sec].second->fake_boundary_vertex = true;
			//	}
			//}
			break;
		case 3:
			break;
		default:
			break;
		}

		if (compute_hausdorff == true)
		{
			double temp_sum_haus_dis = meanhausdorff_distance * pmesh->pVertexList.size();
			for (set<unsigned>::iterator it = temp_bplist.begin(); it != temp_bplist.end(); it++)
			{
				unsigned temp_ind = *it;
				Vector3d bou_ver(pmesh->pVertexList[temp_ind]->point()[0], pmesh->pVertexList[temp_ind]->point()[1], pmesh->pVertexList[temp_ind]->point()[2]);

				temp_sum_haus_dis -= pmesh->pVertexList[temp_ind]->slab_hausdorff_dist;

				double min_dis = DBL_MAX;
				unsigned min_index = -1;
				for (unsigned j = 0; j < vertices.size(); j++)
				{
					if (vertices[j].first)
					{
						Sphere ma_ver = vertices[j].second->sphere;
						double temp_length = abs((bou_ver - ma_ver.center).Length() - ma_ver.radius);
						if (temp_length >= 0 && temp_length < min_dis)
						{
							min_dis = temp_length;
							min_index = j;
						}
					}
				}

				if (min_index != -1)
				{	
					double temp_near_dis = NearestPoint(bou_ver, min_index);
					min_dis = min(temp_near_dis, min_dis);

					maxhausdorff_distance = max(maxhausdorff_distance, min_dis);
					pmesh->pVertexList[temp_ind]->slab_hansdorff_index = min_index;
					pmesh->pVertexList[temp_ind]->slab_hausdorff_dist = min_dis;

					temp_sum_haus_dis += min_dis;
				}
				else
				{
					pmesh->pVertexList[temp_ind]->slab_hansdorff_index = vid_tgt;
					double temp_len = abs((vertices[vid_tgt].second->sphere.center - bou_ver).Length() - vertices[vid_tgt].second->sphere.radius);
					if (temp_len >= 0)
					{
						pmesh->pVertexList[temp_ind]->slab_hausdorff_dist = temp_len;
						maxhausdorff_distance = max(maxhausdorff_distance,temp_len);
						temp_sum_haus_dis += temp_len;
					}
				}
			}
			meanhausdorff_distance = temp_sum_haus_dis / pmesh->pVertexList.size();
		}

		for (std::set<unsigned>::iterator si = vertices[vid_tgt].second->edges_.begin(); si != vertices[vid_tgt].second->edges_.end(); si ++)
		{
			unsigned fir = edges[*si].second->vertices_.first;
			unsigned sec = edges[*si].second->vertices_.second;

			switch(boundary_compute_scale)
			{
			case 1:
				if (!vertices[fir].second->fake_boundary_vertex || !vertices[sec].second->fake_boundary_vertex)
					continue;
				break;
			case 2:
				if (!vertices[fir].second->fake_boundary_vertex && !vertices[sec].second->fake_boundary_vertex)
					//if (vertices[fir].second->boundary_edge_vec.size() < 2 && vertices[sec].second->boundary_edge_vec.size() < 2)
					continue;
				break;
			case 3:
				if (!vertices[fir].second->fake_boundary_vertex && !vertices[sec].second->fake_boundary_vertex)
					continue;
				break;
			default:
				break;
			}
			EvaluateEdgeHausdorffCost(*si);
			//ReEvaluateEdgeHausdorffCost(*si);
			boundary_edge_collapses_queue.push(EdgeInfo(*si, edges[*si].second->collapse_cost));
		}
	}

	return true;
}

bool SlabMesh::MinCostEdgeCollapse(unsigned& eid, CollapseContext ctx){
	//merge 2 vertices of the edge first, then move the combined vertex to the preferred point and resize it.
	unsigned v1, v2;
	v1 = edges[eid].second->vertices_.first;
	v2 = edges[eid].second->vertices_.second;

	const char* q_name = (ctx == CollapseContext::Spike)    ? "spike"
	                   : (ctx == CollapseContext::Boundary) ? "boundary"
	                   :                                      "edge";

	// Spike collapses bypass cluster/topo CanMerge checks intentionally,
	// but non-manifold check always applies.
	if (ctx != CollapseContext::Spike) {
		RejectionReason reason;
		ReasonPrimitives prims;
		if (!CanMerge(v1, v2, &reason, &prims))
		{ LogCollapseRejection(q_name, eid, v1, v2, edges[eid].second->collapse_cost, reason, std::move(prims)); return false; }
	} else {
		RejectionReason nm_reason = RejectionReason::NonManifold_LinkCondition;
		ReasonPrimitives prims;
		if (WouldCreateNonManifold(v1, v2, &nm_reason, &prims))
		{ LogCollapseRejection(q_name, eid, v1, v2, edges[eid].second->collapse_cost, nm_reason, std::move(prims)); return false; }
	}

	Wm4::Matrix4d A = edges[eid].second->slab_A;
	Wm4::Vector4d b = edges[eid].second->slab_b;
	double c = edges[eid].second->slab_c;
	Sphere sphere = edges[eid].second->sphere;
	double hyperbolic_weight = vertices[v1].second->hyperbolic_weight + vertices[v2].second->hyperbolic_weight;

	//// ���ںϲ��ᷢ�����˸ı�ıߣ����������кϲ�
	if (!edges[eid].second->topo_contractable)
	{
		ReasonPrimitives prims; prims.vertices = { v1, v2 }; prims.edges = { {v1, v2} };
		LogCollapseRejection(q_name, eid, v1, v2, edges[eid].second->collapse_cost,
		                     RejectionReason::TopoNotContractable, std::move(prims));
		return false;
	}

	// ��������˷�ת�Ĵ�����ʽ
	if (prevent_inversion == true)
	{
		// ��������תʱ��ѡȡû������ת�ķ�ʽ���кϲ�
		if (!Contractible(v1, v2, sphere.center))
		{
			Wm4::Vector4d lamdar;
			double coll_cost = 0.0;

			int count = 0;		
			double *collapse_costs = new double[3];
			Sphere *min_sphere = new Sphere[3];
			Vector4d min_vertex;
			int min_index = 0;
			if (Contractible(v1, v2, vertices[v1].second->sphere.center))
			{
				min_sphere[count] = vertices[v1].second->sphere;
				min_vertex = Vector4d(min_sphere[count].center.X(), min_sphere[count].center.Y(), min_sphere[count].center.Z(), min_sphere[count].radius);
				collapse_costs[count] = 0.5 * (min_vertex * A).Dot(min_vertex) - b.Dot(min_vertex) + c;
				count++;
			}
			if (Contractible(v1, v2, vertices[v2].second->sphere.center))
			{
				min_sphere[count] = vertices[v2].second->sphere;
				min_vertex = Vector4d(min_sphere[count].center.X(), min_sphere[count].center.Y(), min_sphere[count].center.Z(), min_sphere[count].radius);
				collapse_costs[count] = 0.5 * (min_vertex * A).Dot(min_vertex) - b.Dot(min_vertex) + c;
				count++;
			}
			if (Contractible(v1, v2, (vertices[v1].second->sphere.center + vertices[v2].second->sphere.center) / 2.0))
			{
				min_sphere[count] = (vertices[v1].second->sphere + vertices[v2].second->sphere) * 0.5;
				min_vertex = Vector4d(min_sphere[count].center.X(), min_sphere[count].center.Y(), min_sphere[count].center.Z(), min_sphere[count].radius);
				collapse_costs[count] = 0.5 * (min_vertex * A).Dot(min_vertex) - b.Dot(min_vertex) + c;
				count++;
			}

			if (count == 1)
			{
				lamdar = Vector4d(min_sphere[0].center.X(), min_sphere[0].center.Y(), min_sphere[0].center.Z(), min_sphere[0].radius);
				coll_cost = collapse_costs[0];
			}else if (count == 2)
			{
				min_index = collapse_costs[0] > collapse_costs[1] ? 1 : 0;
				lamdar = Vector4d(min_sphere[min_index].center.X(), min_sphere[min_index].center.Y(), min_sphere[min_index].center.Z(), min_sphere[min_index].radius);
				coll_cost = collapse_costs[min_index];
			}else if (count == 3)
			{
				if (collapse_costs[0] >= collapse_costs[1]) min_index = 1;
				min_index = collapse_costs[min_index] > collapse_costs[2] ? 2 : min_index;
				lamdar = Vector4d(min_sphere[min_index].center.X(), min_sphere[min_index].center.Y(), min_sphere[min_index].center.Z(), min_sphere[min_index].radius);
				coll_cost = collapse_costs[min_index];
			}else{
				coll_cost += 1e9;
				ReasonPrimitives prims; prims.vertices = { v1, v2 }; prims.edges = { {v1, v2} };
				LogCollapseRejection(q_name, eid, v1, v2, edges[eid].second->collapse_cost,
				                     RejectionReason::InversionWouldOccur, std::move(prims));
			}
			delete [] collapse_costs;
			delete [] min_sphere;

			edges[eid].second->qem_error = coll_cost;
			//coll_cost = (coll_cost + k) * edges[eid].second->hyperbolic_weight * edges[eid].second->hyperbolic_weight 
			//	* edges[eid].second->hyperbolic_weight	* edges[eid].second->hyperbolic_weight
			//	* edges[eid].second->hyperbolic_weight	* edges[eid].second->hyperbolic_weight;

			coll_cost = (coll_cost + k) * edges[eid].second->hyperbolic_weight * edges[eid].second->hyperbolic_weight;

			edges[eid].second->collapse_cost = coll_cost;

			edges[eid].second->sphere.center = Wm4::Vector3d(lamdar.X(), lamdar.Y(), lamdar.Z());
			edges[eid].second->sphere.radius = lamdar.W();

			edge_collapses_queue.push(EdgeInfo(eid, coll_cost));

			return false;
		}
	}

	// �����Ǳ߽�߽��м򻯻����ڲ��߽��м�
	if (edges[eid].second->faces_.size() <= 1)
		simplified_boundary_edges++;
	else
		simplified_inside_edges++;
	// ÿ�μ���1000����֮������򻯽��
	if (simplified_boundary_edges + simplified_inside_edges == 1000)
	{
		ExportSimplifyResult();
		simplified_inside_edges = simplified_boundary_edges = 0;
	}

	set<unsigned> temp_bplist;
	if (compute_hausdorff)
	{
		for (set<unsigned>::iterator it = vertices[v1].second->bplist.begin(); it != vertices[v1].second->bplist.end(); it++)
			temp_bplist.insert(*it);
		for (set<unsigned>::iterator it = vertices[v2].second->bplist.begin(); it != vertices[v2].second->bplist.end(); it++)
			temp_bplist.insert(*it);
	}

	unsigned temp_related_face = vertices[v1].second->related_face + vertices[v2].second->related_face;
	double temp_mean_squre_error = edges[eid].second->collapse_cost < 0 ? 0 : edges[eid].second->collapse_cost / temp_related_face;
	max_mean_squre_error = max(temp_mean_squre_error, max_mean_squre_error);

#ifdef QMAT_WITH_POLYSCOPE
	if (on_collapse_cb) {
		on_collapse_cb(
			v1, vertices[v1].second->sphere.center, vertices[v1].second->sphere.radius,
			v2, vertices[v2].second->sphere.center, vertices[v2].second->sphere.radius,
			sphere);
	}
#endif

	// Capture positions and bplists before merge for history recording.
	const auto& hc1 = vertices[v1].second->sphere.center;
	const auto& hc2 = vertices[v2].second->sphere.center;
	std::array<double,3> hist_pos1 = {hc1.X(), hc1.Y(), hc1.Z()};
	std::array<double,3> hist_pos2 = {hc2.X(), hc2.Y(), hc2.Z()};
	std::vector<unsigned> hist_bp1(vertices[v1].second->nmn_bplist.begin(),
	                               vertices[v1].second->nmn_bplist.end());
	std::vector<unsigned> hist_bp2(vertices[v2].second->nmn_bplist.begin(),
	                               vertices[v2].second->nmn_bplist.end());

	// ── Fold-over check ───────────────────────────────────────────────────────
	// sphere.center is the finalised target position.  Check that none of the
	// new edges produced by this collapse would geometrically cross an existing
	// edge.  This is the check that WouldCreateNonManifold (topology-only) and
	// Contractible (face-normal-only) both miss for 1-D boundary loops.
	{
		ReasonPrimitives prims;
		if (WouldCreateFoldOver(v1, v2, sphere.center, &prims))
		{ LogCollapseRejection(q_name, eid, v1, v2, edges[eid].second->collapse_cost, RejectionReason::WouldCreateFoldOver, std::move(prims)); return false; }
	}

	{
		using CT = SlabVertex::ClusterType;
		auto isChainType = [](CT c) {
			return c == CT::MS_Boundary      || c == CT::MS_Seam ||
			       c == CT::MS_Seam_Boundary || c == CT::MS_Sheet_Boundary;
		};
		const CT c1 = vertices[v1].second->nmn_cluster_type;
		const CT c2 = vertices[v2].second->nmn_cluster_type;
		if (isChainType(c1) || isChainType(c2)) {
			ReasonPrimitives prims;
			if (WouldExceedCurvatureThreshold(v1, v2, &prims))
			{ LogCollapseRejection(q_name, eid, v1, v2, edges[eid].second->collapse_cost, RejectionReason::WouldExceedCurvatureThreshold, std::move(prims)); return false; }
		}
	}

	unsigned vid_tgt;
	if(MergeVertices(v1, v2, vid_tgt)){
		vertices[vid_tgt].second->slab_A = A;
		vertices[vid_tgt].second->slab_b = b;
		vertices[vid_tgt].second->slab_c = c;
		vertices[vid_tgt].second->sphere = sphere;
		vertices[vid_tgt].second->related_face = temp_related_face;
		vertices[vid_tgt].second->mean_square_error = temp_mean_squre_error;
		vertices[vid_tgt].second->hyperbolic_weight = hyperbolic_weight;

		// Record collapse in history (bplist_after and clusters_after captured post-merge).
		{
			std::vector<unsigned> hist_bpt(vertices[vid_tgt].second->nmn_bplist.begin(),
			                               vertices[vid_tgt].second->nmn_bplist.end());
			std::vector<std::vector<unsigned>> hist_clusters;
			for (const auto& cl : vertices[vid_tgt].second->nmn_bplist_clusters)
				hist_clusters.push_back(std::vector<unsigned>(cl.begin(), cl.end()));
			const unsigned step = history.TotalSteps();
			history.Record(step, v1, v2, vid_tgt,
			               hist_pos1, hist_pos2,
			               std::move(hist_bp1), std::move(hist_bp2),
			               std::move(hist_bpt), std::move(hist_clusters));
			if (history.keyframe_interval > 0 &&
			    (int)history.TotalSteps() % history.keyframe_interval == 0)
				history.TakeKeyframe(history.TotalSteps(), *this);
		}

		// Refresh topology for the new merged vertex.
		RecomputeVertexTopology(vid_tgt);

		// ����������Ϣ
		InitialTopologyProperty(vid_tgt);

		for (std::set<unsigned>::iterator si = vertices[vid_tgt].second->edges_.begin(); si != vertices[vid_tgt].second->edges_.end(); si ++)
		{
			unsigned fir = edges[*si].second->vertices_.first;
			unsigned sec = edges[*si].second->vertices_.second;

			EvaluateEdgeCollapseCost(*si);
			ComputeEdgeCone(*si);
			edge_collapses_queue.push(EdgeInfo(*si, edges[*si].second->collapse_cost));
		}

		if (compute_hausdorff)
		{
			double temp_sum_haus_dis = meanhausdorff_distance * pmesh->pVertexList.size();
			for (set<unsigned>::iterator it = temp_bplist.begin(); it != temp_bplist.end(); it++)
			{
				unsigned temp_ind = *it;
				Vector3d bou_ver(pmesh->pVertexList[temp_ind]->point()[0], pmesh->pVertexList[temp_ind]->point()[1], pmesh->pVertexList[temp_ind]->point()[2]);
				bou_ver /= pmesh->bb_diagonal_length;

				temp_sum_haus_dis -= pmesh->pVertexList[temp_ind]->slab_hausdorff_dist;

				double min_dis = DBL_MAX;
				unsigned min_index = -1;
				for (unsigned j = 0; j < vertices.size(); j++)
				{
					if (vertices[j].first)
					{
						Sphere ma_ver = vertices[j].second->sphere;
						double temp_length = abs((bou_ver - ma_ver.center).Length() - ma_ver.radius);
						if (temp_length >= 0 && temp_length < min_dis)
						{
							min_dis = temp_length;
							min_index = j;
						}
					}
				}

				if (min_index != -1)
				{	
					double temp_near_dis = NearestPoint(bou_ver, min_index);
					min_dis = min(temp_near_dis, min_dis);

					vertices[min_index].second->bplist.insert(temp_ind);
					maxhausdorff_distance = max(maxhausdorff_distance, min_dis);
					pmesh->pVertexList[temp_ind]->slab_hansdorff_index = min_index;
					pmesh->pVertexList[temp_ind]->slab_hausdorff_dist = min_dis;

					temp_sum_haus_dis += min_dis;
				}
				else
				{
					pmesh->pVertexList[temp_ind]->slab_hansdorff_index = vid_tgt;
					double temp_len = abs((vertices[vid_tgt].second->sphere.center - bou_ver).Length() - vertices[vid_tgt].second->sphere.radius);
					if (temp_len >= 0)
					{
						pmesh->pVertexList[temp_ind]->slab_hausdorff_dist = temp_len;
						maxhausdorff_distance = max(maxhausdorff_distance,temp_len);
						temp_sum_haus_dis += temp_len;
					}
				}
			}
			meanhausdorff_distance = temp_sum_haus_dis / pmesh->pVertexList.size();
		}
	}

	return true;
}

void SlabMesh::EvaluateEdgeCollapseCost(unsigned eid){
	if (!edges[eid].first)
		return ;

	unsigned v1, v2;
	v1 = edges[eid].second->vertices_.first;
	v2 = edges[eid].second->vertices_.second; 

	double weight = vertices[v1].second->hyperbolic_weight + vertices[v2].second->hyperbolic_weight;

	// set the hyperbolic weight to the related edge
	switch(hyperbolic_weight_type)
	{
	case 1:
		edges[eid].second->hyperbolic_weight = GetHyperbolicLength(eid);
		break;
	case 2:
		edges[eid].second->hyperbolic_weight = GetHyperbolicLength(eid);
		break;
	case 3:
		edges[eid].second->hyperbolic_weight = GetRatioHyperbolicEuclid(eid);
		break;
	default:
		break;
	}


	//double w1 = 1e-5, w2 = 1e-5;
	double w1 = 1.0, w2 = 1.0;
	//// �Բ�ͬratio�ı߽���ӳ�䴦����С��0.2�Ĳ���������(0.2,1)ӳ�䵽(2, 10)
	//w1 = edges[eid].second->hyperbolic_weight * edges[eid].second->hyperbolic_weight * edges[eid].second->hyperbolic_weight;
	//w2 = w1;
	//if (edges[eid].second->hyperbolic_weight >= 0.2)
	//{
	//	w1 = (edges[eid].second->hyperbolic_weight - 0.2) / 0.8 * (10 - 2) + 2;
	//	w2 = w1;
	//}

	Wm4::Matrix4d A1 = vertices[v1].second->slab_A;
	Wm4::Matrix4d A2 = vertices[v2].second->slab_A; 
	Wm4::Matrix4d add_A1 = vertices[v1].second->add_A * w1;
	Wm4::Matrix4d add_A2 = vertices[v2].second->add_A * w2;

	Wm4::Vector4d b1 = vertices[v1].second->slab_b;
	Wm4::Vector4d b2 = vertices[v2].second->slab_b;	
	Wm4::Vector4d add_b1 = vertices[v1].second->add_b * w1;
	Wm4::Vector4d add_b2 = vertices[v2].second->add_b * w2;

	double c1 = vertices[v1].second->slab_c;
	double c2 = vertices[v2].second->slab_c;
	double add_c1 = vertices[v1].second->add_c * w1;
	double add_c2 = vertices[v2].second->add_c * w2;

	edges[eid].second->slab_A = A1 + A2;
	edges[eid].second->slab_b = b1 + b2;
	edges[eid].second->slab_c = c1 + c2;

	Matrix4d inverse_A_matrix = edges[eid].second->slab_A.Inverse();
	Wm4::Vector4d lamdar;
	double coll_cost = 0.0;

	if ((vertices[v1].second->saved_vertex && !vertices[v2].second->saved_vertex) ||
		(vertices[v2].second->saved_vertex && !vertices[v1].second->saved_vertex))
	{
		Sphere mid_sphere; 
		if (vertices[v1].second->saved_vertex)
			mid_sphere = vertices[v1].second->sphere;
		else
			mid_sphere = vertices[v2].second->sphere;

		lamdar = Vector4d(mid_sphere.center.X(), mid_sphere.center.Y(), mid_sphere.center.Z(), mid_sphere.radius);
	}	
	else if ((vertices[v1].second->saved_vertex && vertices[v2].second->saved_vertex))
	{
		//edges[eid].second->collapse_cost = DBL_MAX;
		//GetBestBoundaryPoint(eid);
		double collapse_costs[3];
		Sphere min_sphere[3];
		Vector4d min_vertex;
		int min_index = 0;

		min_sphere[0] = vertices[v1].second->sphere;
		min_vertex = Vector4d(min_sphere[0].center.X(), min_sphere[0].center.Y(), min_sphere[0].center.Z(), min_sphere[0].radius);
		collapse_costs[0] = 0.5 * (min_vertex * edges[eid].second->slab_A).Dot(min_vertex) 
			- edges[eid].second->slab_b.Dot(min_vertex) + edges[eid].second->slab_c;
		min_sphere[1] = vertices[v2].second->sphere;
		min_vertex = Vector4d(min_sphere[1].center.X(), min_sphere[1].center.Y(), min_sphere[1].center.Z(), min_sphere[1].radius);
		collapse_costs[1] = 0.5 * (min_vertex * edges[eid].second->slab_A).Dot(min_vertex) 
			- edges[eid].second->slab_b.Dot(min_vertex) + edges[eid].second->slab_c;
		min_sphere[2] = (vertices[v1].second->sphere + vertices[v2].second->sphere) * 0.5;
		min_vertex = Vector4d(min_sphere[2].center.X(), min_sphere[2].center.Y(), min_sphere[2].center.Z(), min_sphere[2].radius);
		collapse_costs[2] = 0.5 * (min_vertex * edges[eid].second->slab_A).Dot(min_vertex) 
			- edges[eid].second->slab_b.Dot(min_vertex) + edges[eid].second->slab_c;

		if (collapse_costs[0] >= collapse_costs[1]) min_index = 1;

		min_index = collapse_costs[min_index] > collapse_costs[2] ? 2 : min_index;

		edges[eid].second->collapse_cost = collapse_costs[min_index];
		edges[eid].second->sphere.center = min_sphere[min_index].center;
		edges[eid].second->sphere.radius = min_sphere[min_index].radius;

		return;
	}
	else
	{
		if (inverse_A_matrix != Matrix4d() || edges[eid].second->faces_.size() == 0)
		{
			// add the boundary preserving.
			//if (edges[eid].second->hyperbolic_weight >= 0.1)
			//{
			edges[eid].second->slab_A = A1 + A2 + add_A1 + add_A2;
			edges[eid].second->slab_b = b1 + b2 + add_b1 + add_b2;
			edges[eid].second->slab_c = c1 + c2 + add_c1 + add_c2;
			//}
			inverse_A_matrix = edges[eid].second->slab_A.Inverse();
			if (inverse_A_matrix != Matrix4d())
			{
				lamdar = inverse_A_matrix * edges[eid].second->slab_b;
				if (lamdar.W() < 0)
				{
					Sphere mid_sphere = (vertices[v1].second->sphere + vertices[v2].second->sphere) * 0.5;
					lamdar = Vector4d(mid_sphere.center.X(), mid_sphere.center.Y(), mid_sphere.center.Z(), mid_sphere.radius);
				}
			}
			else
			{
				// it's now calculate as the middle of the spheres
				Sphere mid_sphere = (vertices[v1].second->sphere + vertices[v2].second->sphere) * 0.5;
				lamdar = Vector4d(mid_sphere.center.X(), mid_sphere.center.Y(), mid_sphere.center.Z(), mid_sphere.radius);
			}
		}
		else
		{
			//// it's now calculate as the middle of the spheres
			//Sphere mid_sphere = (vertices[v1].second->sphere + vertices[v2].second->sphere) * 0.5;
			//lamdar = Vector4d(mid_sphere.center.X(), mid_sphere.center.Y(), mid_sphere.center.Z(), mid_sphere.radius);

			double collapse_costs[3];
			Sphere min_sphere[3];
			Vector4d min_vertex;
			int min_index = 0;

			min_sphere[0] = vertices[v1].second->sphere;
			min_vertex = Vector4d(min_sphere[0].center.X(), min_sphere[0].center.Y(), min_sphere[0].center.Z(), min_sphere[0].radius);
			collapse_costs[0] = 0.5 * (min_vertex * edges[eid].second->slab_A).Dot(min_vertex) 
				- edges[eid].second->slab_b.Dot(min_vertex) + edges[eid].second->slab_c;
			min_sphere[1] = vertices[v2].second->sphere;
			min_vertex = Vector4d(min_sphere[1].center.X(), min_sphere[1].center.Y(), min_sphere[1].center.Z(), min_sphere[1].radius);
			collapse_costs[1] = 0.5 * (min_vertex * edges[eid].second->slab_A).Dot(min_vertex) 
				- edges[eid].second->slab_b.Dot(min_vertex) + edges[eid].second->slab_c;
			min_sphere[2] = (vertices[v1].second->sphere + vertices[v2].second->sphere) * 0.5;
			min_vertex = Vector4d(min_sphere[2].center.X(), min_sphere[2].center.Y(), min_sphere[2].center.Z(), min_sphere[2].radius);
			collapse_costs[2] = 0.5 * (min_vertex * edges[eid].second->slab_A).Dot(min_vertex) 
				- edges[eid].second->slab_b.Dot(min_vertex) + edges[eid].second->slab_c;

			if (collapse_costs[0] >= collapse_costs[1]) min_index = 1;
			min_index = collapse_costs[min_index] > collapse_costs[2] ? 2 : min_index;

			coll_cost = collapse_costs[min_index];
			switch(hyperbolic_weight_type)
			{
			case 1:
				coll_cost = coll_cost * edges[eid].second->hyperbolic_weight;
				break;
			case 2:
				if (weight <= 1e-12)
					coll_cost = 0.0;
				else
					coll_cost = collapse_costs[min_index] / weight;
				break;
			case 3:
				//coll_cost = edges[eid].second->hyperbolic_weight;
				edges[eid].second->qem_error = coll_cost;
				coll_cost = (coll_cost + k) * edges[eid].second->hyperbolic_weight * edges[eid].second->hyperbolic_weight;

				//coll_cost = (coll_cost + k) * edges[eid].second->hyperbolic_weight * edges[eid].second->hyperbolic_weight 
				//	* edges[eid].second->hyperbolic_weight	* edges[eid].second->hyperbolic_weight
				//	* edges[eid].second->hyperbolic_weight	* edges[eid].second->hyperbolic_weight;

				//coll_cost = (coll_cost + k) * edges[eid].second->hyperbolic_weight;
				//coll_cost = (coll_cost + edges[eid].second->hyperbolic_weight) * edges[eid].second->hyperbolic_weight;
				break;
			default:
				break;
			}

			edges[eid].second->collapse_cost = coll_cost;
			edges[eid].second->sphere.center = min_sphere[min_index].center;
			edges[eid].second->sphere.radius = min_sphere[min_index].radius;

			return;
		}
	}

	coll_cost = 0.5 * (lamdar * edges[eid].second->slab_A).Dot(lamdar) 
		- edges[eid].second->slab_b.Dot(lamdar) + edges[eid].second->slab_c;

	// ��������תʱ��ѡȡû������ת�ķ�ʽ���кϲ�
	if (!Contractible(v1, v2, Wm4::Vector3d(lamdar.X(), lamdar.Y(), lamdar.Z())))
	{
		int count = 0;		
		double *collapse_costs = new double[3];
		Sphere *min_sphere = new Sphere[3];
		Vector4d min_vertex;
		int min_index = 0;
		if (!Contractible(v1, v2, vertices[v1].second->sphere.center))
		{
			min_sphere[count] = vertices[v1].second->sphere;
			min_vertex = Vector4d(min_sphere[count].center.X(), min_sphere[count].center.Y(), min_sphere[count].center.Z(), min_sphere[count].radius);
			collapse_costs[count] = 0.5 * (min_vertex * edges[eid].second->slab_A).Dot(min_vertex) 
				- edges[eid].second->slab_b.Dot(min_vertex) + edges[eid].second->slab_c;
			count++;
		}
		if (!Contractible(v1, v2, vertices[v2].second->sphere.center))
		{
			min_sphere[count] = vertices[v2].second->sphere;
			min_vertex = Vector4d(min_sphere[count].center.X(), min_sphere[count].center.Y(), min_sphere[count].center.Z(), min_sphere[count].radius);
			collapse_costs[count] = 0.5 * (min_vertex * edges[eid].second->slab_A).Dot(min_vertex) 
				- edges[eid].second->slab_b.Dot(min_vertex) + edges[eid].second->slab_c;
			count++;
		}
		if (!Contractible(v1, v2, (vertices[v1].second->sphere.center + vertices[v2].second->sphere.center) / 2.0))
		{
			min_sphere[count] = (vertices[v1].second->sphere + vertices[v2].second->sphere) * 0.5;
			min_vertex = Vector4d(min_sphere[count].center.X(), min_sphere[count].center.Y(), min_sphere[count].center.Z(), min_sphere[count].radius);
			collapse_costs[count] = 0.5 * (min_vertex * edges[eid].second->slab_A).Dot(min_vertex) 
				- edges[eid].second->slab_b.Dot(min_vertex) + edges[eid].second->slab_c;
			count++;
		}

		if (count == 1)
		{
			lamdar = Vector4d(min_sphere[0].center.X(), min_sphere[0].center.Y(), min_sphere[0].center.Z(), min_sphere[0].radius);
			//coll_cost = collapse_costs[0];
			//coll_cost = 0.0;
		}else if (count == 2)
		{
			min_index = collapse_costs[0] > collapse_costs[1] ? 1 : 0;
			//min_index = min_sphere[0].radius > min_sphere[1].radius ? 0 : 1;
			lamdar = Vector4d(min_sphere[min_index].center.X(), min_sphere[min_index].center.Y(), min_sphere[min_index].center.Z(), min_sphere[min_index].radius);
			//coll_cost = collapse_costs[min_index];
			//coll_cost = 0.0;
		}else if (count == 3)
		{
			if (collapse_costs[0] >= collapse_costs[1]) min_index = 1;
			min_index = collapse_costs[min_index] > collapse_costs[2] ? 2 : min_index;
			//if (min_sphere[0].radius >= min_sphere[1].radius) min_index = 0;
			//min_index = min_sphere[min_index].radius > min_sphere[2].radius ? min_index : 2;
			lamdar = Vector4d(min_sphere[min_index].center.X(), min_sphere[min_index].center.Y(), min_sphere[min_index].center.Z(), min_sphere[min_index].radius);
			//coll_cost = collapse_costs[min_index];
			//coll_cost = 0.0;
		}else
			coll_cost += 1e9;
		delete [] collapse_costs;
		delete [] min_sphere;
	}

	switch(hyperbolic_weight_type)
	{
	case 1:
		coll_cost = coll_cost * edges[eid].second->hyperbolic_weight;
		break;
	case 2:
		if (weight <= 1e-12)
			coll_cost = 0.0;
		else
			coll_cost = coll_cost / weight;
		break;
	case 3:
		//coll_cost = coll_cost * edges[eid].second->hyperbolic_weight;
		//coll_cost = edges[eid].second->hyperbolic_weight;
		edges[eid].second->qem_error = coll_cost;
		coll_cost = (coll_cost + k) * edges[eid].second->hyperbolic_weight * edges[eid].second->hyperbolic_weight;

		//coll_cost = (coll_cost + k) * edges[eid].second->hyperbolic_weight * edges[eid].second->hyperbolic_weight 
		//	* edges[eid].second->hyperbolic_weight * edges[eid].second->hyperbolic_weight
		//	* edges[eid].second->hyperbolic_weight * edges[eid].second->hyperbolic_weight;

		//coll_cost = (coll_cost + k) * edges[eid].second->hyperbolic_weight;
		//coll_cost = (coll_cost + edges[eid].second->hyperbolic_weight) * edges[eid].second->hyperbolic_weight;
		break;
	default:
		break; 
	}

	edges[eid].second->collapse_cost = coll_cost;
	edges[eid].second->sphere.center = Wm4::Vector3d(lamdar.X(), lamdar.Y(), lamdar.Z());
	edges[eid].second->sphere.radius = lamdar.W();
}

void SlabMesh::EvaluateEdgeHausdorffCost(unsigned eid)
{
	if (!edges[eid].first)
		return ;

	unsigned v1, v2;
	v1 = edges[eid].second->vertices_.first;
	v2 = edges[eid].second->vertices_.second;

	double w1 = 1.0, w2 = 1.0;

	Wm4::Matrix4d A1 = vertices[v1].second->slab_A;
	Wm4::Matrix4d A2 = vertices[v2].second->slab_A;
	Wm4::Matrix4d add_A1 = vertices[v1].second->add_A * w1;
	Wm4::Matrix4d add_A2 = vertices[v2].second->add_A * w2;

	Wm4::Vector4d b1 = vertices[v1].second->slab_b;
	Wm4::Vector4d b2 = vertices[v2].second->slab_b;	
	Wm4::Vector4d add_b1 = vertices[v1].second->add_b * w1;
	Wm4::Vector4d add_b2 = vertices[v2].second->add_b * w2;

	double c1 = vertices[v1].second->slab_c;
	double c2 = vertices[v2].second->slab_c;
	double add_c1 = vertices[v1].second->add_c * w1;
	double add_c2 = vertices[v2].second->add_c * w2;

	edges[eid].second->slab_A = A1 + A2;
	edges[eid].second->slab_b = b1 + b2;
	edges[eid].second->slab_c = c1 + c2;

	Matrix4d inverse_A_matrix = edges[eid].second->slab_A.Inverse();
	Wm4::Vector4d lamdar;
	double coll_cost;

	if (inverse_A_matrix != Matrix4d())
	{
		// add the boundary preserving.
		edges[eid].second->slab_A = A1 + A2 + add_A1 + add_A2;
		edges[eid].second->slab_b = b1 + b2 + add_b1 + add_b2;
		edges[eid].second->slab_c = c1 + c2 + add_c1 + add_c2;

		inverse_A_matrix = edges[eid].second->slab_A.Inverse();
		if (inverse_A_matrix != Matrix4d())
		{
			lamdar = inverse_A_matrix * edges[eid].second->slab_b;
			if (lamdar.W() < 0)
			{
				Sphere mid_sphere = (vertices[v1].second->sphere + vertices[v2].second->sphere) * 0.5;
				lamdar = Vector4d(mid_sphere.center.X(), mid_sphere.center.Y(), mid_sphere.center.Z(), mid_sphere.radius);
			}
		}
		else
		{
			// it's now calculate as the middle of the spheres
			Sphere mid_sphere = (vertices[v1].second->sphere + vertices[v2].second->sphere) * 0.5;
			lamdar = Vector4d(mid_sphere.center.X(), mid_sphere.center.Y(), mid_sphere.center.Z(), mid_sphere.radius);
		}
	}
	else
	{
		// it's now calculate as the middle of the spheres
		Sphere mid_sphere = (vertices[v1].second->sphere + vertices[v2].second->sphere) * 0.5;
		lamdar = Vector4d(mid_sphere.center.X(), mid_sphere.center.Y(), mid_sphere.center.Z(), mid_sphere.radius);
	}

	set<unsigned> temp_bplist;
	for (set<unsigned>::iterator it = vertices[v1].second->bplist.begin(); it != vertices[v1].second->bplist.end(); it++)
		temp_bplist.insert(*it);
	for (set<unsigned>::iterator it = vertices[v2].second->bplist.begin(); it != vertices[v2].second->bplist.end(); it++)
		temp_bplist.insert(*it);

	double max_hausdorff = 0;
	for (set<unsigned>::iterator it = temp_bplist.begin(); it != temp_bplist.end(); it++)
	{
		unsigned temp_ind = *it;
		Vector3d bou_ver(pmesh->pVertexList[temp_ind]->point()[0], pmesh->pVertexList[temp_ind]->point()[1], pmesh->pVertexList[temp_ind]->point()[2]);

		//double min_dis = DBL_MIN;
		//unsigned min_index = -1;
		//for (unsigned j = 0; j < vertices.size(); j++)
		//{
		//	if (vertices[j].first && (j != v1 && j != v2))
		//	{
		//		Sphere ma_ver = vertices[j].second->sphere;
		//		double temp_length = abs((bou_ver - ma_ver.center).Length() - ma_ver.radius);
		//		if (temp_length >= 0 && temp_length < min_dis)
		//		{
		//			min_dis = temp_length;
		//			min_index = j;
		//		}
		//	}
		//}

		//double temp_near_dis = NearestPoint(bou_ver, min_index);
		//min_dis = min(temp_near_dis, min_dis);

		// ��������С���Ƿ��������ɵĵ�
		double temp_length = abs((bou_ver - Wm4::Vector3d(lamdar.X(), lamdar.Y(), lamdar.Z())).Length() - lamdar.W());
		//min_dis = max(temp_length, min_dis);

		max_hausdorff = max(temp_length, max_hausdorff);
	}
	//set<unsigned> neighbors_v, tdneighbors_v;
	//GetNeighborVertices(v1, neighbors_v);
	//GetNeighborVertices(v2, tdneighbors_v);
	//neighbors_v.insert(tdneighbors_v.begin(), tdneighbors_v.end());
	//neighbors_v.erase(v1);
	//neighbors_v.erase(v2);

	//set< set<unsigned> > neighbors_f;
	//for(set<unsigned>::iterator si = vertices[v1].second->faces_.begin();
	//	si != vertices[v1].second->faces_.end(); si ++)
	//{
	//	if(!faces[*si].second->HasVertex(v2))
	//	{
	//		set<unsigned> vset = faces[*si].second->vertices_;
	//		vset.erase(v1);
	//		neighbors_f.insert(vset);
	//	}
	//}
	//for(set<unsigned>::iterator si = vertices[v2].second->faces_.begin();
	//	si != vertices[v2].second->faces_.end(); si ++)
	//{
	//	if(!faces[*si].second->HasVertex(v1))
	//	{
	//		set<unsigned> vset = faces[*si].second->vertices_;
	//		vset.erase(v2);
	//		neighbors_f.insert(vset);
	//	}
	//}
	//max_hausdorff = EvaluateVertexDistanceErrorEnvelope(lamdar, neighbors_v, neighbors_f, temp_bplist);

	coll_cost = max_hausdorff;

	double temp_coll_cost = 0.5 * (lamdar * edges[eid].second->slab_A).Dot(lamdar) 
		- edges[eid].second->slab_b.Dot(lamdar) + edges[eid].second->slab_c;

	edges[eid].second->collapse_cost = coll_cost;
	edges[eid].second->sphere.center = Wm4::Vector3d(lamdar.X(), lamdar.Y(), lamdar.Z());
	edges[eid].second->sphere.radius = lamdar.W();
}

void SlabMesh::ReEvaluateEdgeHausdorffCost(unsigned eid)
{
	if (!edges[eid].first)
		return ;

	unsigned v[2];
	v[0] = edges[eid].second->vertices_.first;
	v[1] = edges[eid].second->vertices_.second;

	Wm4::Matrix4d A[2];
	Wm4::Vector4d b[2];
	double c[2]  = {0, 0};

	for(unsigned i = 0; i < 2; i++)
	{
		SlabVertex sv = *vertices[v[i]].second;
		std::set<unsigned> fset = sv.faces_;
		Vector4d C1(sv.sphere.center.X(), sv.sphere.center.Y(), sv.sphere.center.Z(), sv.sphere.radius);

		for (set<unsigned>::iterator si = fset.begin(); si != fset.end(); si++)
		{
			SlabFace sf = *faces[*si].second;

			if (sf.valid_st == false || sf.st[0].normal == Vector3d(0., 0., 0.) || 
				sf.st[1].normal == Vector3d(0., 0., 0.))
				continue;

			Vector4d normal1(sf.st[0].normal.X(), sf.st[0].normal.Y(), sf.st[0].normal.Z(), 1.0);
			Vector4d normal2(sf.st[1].normal.X(), sf.st[1].normal.Y(), sf.st[1].normal.Z(), 1.0);

			// compute the matrix of A
			Matrix4d temp_A1, temp_A2;
			temp_A1.MakeTensorProduct(normal1, normal1);
			temp_A2.MakeTensorProduct(normal2, normal2);
			temp_A1 *= 2.0;
			temp_A2 *= 2.0;

			// compute the matrix of b
			double normal_mul_point1 = normal1.Dot(C1);
			double normal_mul_point2 = normal2.Dot(C1);
			Wm4::Vector4d temp_b1 = normal1 * 2 * normal_mul_point1;
			Wm4::Vector4d temp_b2 = normal2 * 2 * normal_mul_point2;

			//compute c
			double temp_c1 = normal_mul_point1 * normal_mul_point1;
			double temp_c2 = normal_mul_point2 * normal_mul_point2;

			A[i] += temp_A1;
			A[i] += temp_A2;
			b[i] += temp_b1;
			b[i] += temp_b2;
			c[i] += temp_c1;
			c[i] += temp_c2;
		}
	}

	edges[eid].second->slab_A = A[0] + A[1];
	edges[eid].second->slab_b = b[0] + b[1];
	edges[eid].second->slab_c = c[0] + c[1];

	Matrix4d inverse_A_matrix = edges[eid].second->slab_A.Inverse();
	Wm4::Vector4d lamdar;
	double coll_cost;

	if (inverse_A_matrix != Matrix4d())
	{
		lamdar = inverse_A_matrix * edges[eid].second->slab_b;
		if (lamdar.W() < 0)
		{
			Sphere mid_sphere = (vertices[v[0]].second->sphere + vertices[v[1]].second->sphere) * 0.5;
			lamdar = Vector4d(mid_sphere.center.X(), mid_sphere.center.Y(), mid_sphere.center.Z(), mid_sphere.radius);
		}
	}
	else
	{
		// it's now calculate as the middle of the spheres
		Sphere mid_sphere = (vertices[v[0]].second->sphere + vertices[v[1]].second->sphere) * 0.5;
		lamdar = Vector4d(mid_sphere.center.X(), mid_sphere.center.Y(), mid_sphere.center.Z(), mid_sphere.radius);
	}

	set<unsigned> temp_bplist;
	for (set<unsigned>::iterator it = vertices[v[0]].second->bplist.begin(); it != vertices[v[0]].second->bplist.end(); it++)
		temp_bplist.insert(*it);
	for (set<unsigned>::iterator it = vertices[v[1]].second->bplist.begin(); it != vertices[v[1]].second->bplist.end(); it++)
		temp_bplist.insert(*it);

	double max_hausdorff = 0;
	for (set<unsigned>::iterator it = temp_bplist.begin(); it != temp_bplist.end(); it++)
	{
		unsigned temp_ind = *it;
		Vector3d bou_ver(pmesh->pVertexList[temp_ind]->point()[0], pmesh->pVertexList[temp_ind]->point()[1], pmesh->pVertexList[temp_ind]->point()[2]);

		// ��������С���Ƿ��������ɵĵ�
		double temp_length = abs((bou_ver - Wm4::Vector3d(lamdar.X(), lamdar.Y(), lamdar.Z())).Length() - lamdar.W());
		//min_dis = max(temp_length, min_dis);

		max_hausdorff = max(temp_length, max_hausdorff);
	}

	coll_cost = max_hausdorff;

	double temp_coll_cost = 0.5 * (lamdar * edges[eid].second->slab_A).Dot(lamdar) 
		- edges[eid].second->slab_b.Dot(lamdar) + edges[eid].second->slab_c;

	edges[eid].second->collapse_cost = coll_cost;
	edges[eid].second->sphere.center = Wm4::Vector3d(lamdar.X(), lamdar.Y(), lamdar.Z());
	edges[eid].second->sphere.radius = lamdar.W();
}

void SlabMesh::Simplify(int threshold){

	// ���򻯵�С��50������ʱ�������������˵�ı߽��кϲ�
	if (numVertices <= 100)
	{
		if (initial_boundary_preserve == false)
		{
			initial_boundary_preserve = true;
			InitialTopologyProperty();
			for (int i = 0; i < vertices.size(); i++)
			{
				if (vertices[i].first)
				{
					set<unsigned> fir_edges = vertices[i].second->edges_;
					for (set<unsigned>::iterator si = fir_edges.begin(); si != fir_edges.end(); si++)
					{
						int index = edges[*si].second->vertices_.first == i ? 
							edges[*si].second->vertices_.second : edges[*si].second->vertices_.first;

						if (vertices[index].second->edges_.size() == 1 && vertices[index].second->faces_.size() == 0)
						{
							edges[*si].second->topo_contractable = false;
						}
					}
				}
			}
		}
	}

	int deleteSphereNum = 0;

	std::cerr << "[Simplify] Per-phase rejection logs: " << export_prefix << "_rejection_log_{phase}.txt\n";

	// Helper: print a reason-count summary by scanning the phase rejection log.
	auto printPhaseSummary = [&](const std::string& phase, unsigned attempted, unsigned collapsed) {
		unsigned rejected = attempted - collapsed;
		std::cerr << "[Simplify] " << phase << " summary: "
		          << attempted << " attempted, " << collapsed << " collapsed, "
		          << rejected << " rejected\n";
		// Count reasons from log file.
		std::ifstream log(export_prefix + "_rejection_log_" + phase + ".txt");
		if (!log) return;
		std::map<std::string, unsigned> counts;
		std::string line;
		while (std::getline(log, line)) {
			auto pos = line.find("reason=");
			if (pos != std::string::npos)
				counts[line.substr(pos + 7)]++;
		}
		for (auto& kv : counts)
			std::cerr << "    " << kv.first << ": " << kv.second << "\n";
	};

	// Truncates the rejection log for a phase and writes a header with queue
	// size, total MAT vertex count, and per-cluster-type vertex breakdown.
	// Call once before the collapse loop for each phase.
	static const char* ct_names_phase[] = {
		"T0","T1_spike","T2","T3","T4","T5","T1_non_spike",
		"MS_Unknown","MS_Sheet","MS_Seam","MS_Boundary","MS_Junction",
		"MS_Sheet_Boundary","MS_Seam_Boundary","MS_Junction_Boundary"
	};
	auto startPhaseLog = [&](const std::string& phase, unsigned queue_size) {
		std::ofstream log(export_prefix + "_rejection_log_" + phase + ".txt", std::ios::out);
		if (!log) return;
		// Tally active vertices by cluster type.
		std::map<uint8_t, unsigned> ct_counts;
		for (unsigned i = 0; i < (unsigned)vertices.size(); ++i)
			if (vertices[i].first)
				++ct_counts[static_cast<uint8_t>(vertices[i].second->nmn_cluster_type)];
		log << "==============================\n"
		    << " Phase      : " << phase << "\n"
		    << " Queue size : " << queue_size << "\n"
		    << " MAT verts  : " << numVertices << "\n"
		    << " Vertex cluster type breakdown:\n";
		for (auto& kv : ct_counts) {
			const char* name = (kv.first < 15) ? ct_names_phase[kv.first] : "???";
			log << "   " << name << ": " << kv.second << "\n";
		}
		log << "==============================\n";
	};

	// --- Phase 0: collapse all spike edges first ---
	// Spike edges connect to T1 (spike) vertices and should be removed before
	// the main simplification pass so they don't interfere with topology.
	// Phase 0: collapse spike edges in a loop until no T1 vertices remain.
	// Each merge may produce a new T1 vertex, so we iterate to convergence.
	const std::string prefix = export_prefix.empty() ? "mat" : export_prefix;
	int spikeCollapsed = 0;
	int spikePass = 0;
	current_phase = "spike";
	initSpikeCollapseQueue();
	startPhaseLog("spike", (unsigned)spike_collapse_queue.size());
	while (!spike_collapse_queue.empty())
	{
		++spikePass;
		std::cerr << "[Simplify] Phase 0 pass " << spikePass
		          << ": queue size = " << spike_collapse_queue.size()
		          << "  MAT vertices = " << numVertices << "\n";

		while (!spike_collapse_queue.empty())
		{
			EdgeInfo topEdge = spike_collapse_queue.top();
			spike_collapse_queue.pop();
			unsigned eid = topEdge.edge_num;
			if (!edges[eid].first)
			{ LogCollapseRejection("spike", eid, UINT_MAX, UINT_MAX, topEdge.collapse_cost, RejectionReason::StaleEdge, {}); continue; }
			unsigned v1 = edges[eid].second->vertices_.first;
			unsigned v2 = edges[eid].second->vertices_.second;
			if (!ValidVertex(v1) || !ValidVertex(v2))
			{
				ReasonPrimitives prims;
				if (v1 != UINT_MAX && v1 < vertices.size()) prims.vertices.push_back(v1);
				if (v2 != UINT_MAX && v2 < vertices.size()) prims.vertices.push_back(v2);
				LogCollapseRejection("spike", eid, v1, v2, topEdge.collapse_cost,
				                     RejectionReason::InvalidVertex, std::move(prims));
				continue;
			}
			if (MinCostEdgeCollapse(eid, CollapseContext::Spike))
			{
				deleteSphereNum++;
				spikeCollapsed++;
			}
		}
		// Rebuild the queue — some merged vertices may still be T1.
		initSpikeCollapseQueue();
	}

	std::cerr << "[Simplify] Phase 0 done: " << spikePass << " pass(es), "
	          << spikeCollapsed << " total spike collapses, "
	          << "MAT vertices remaining = " << numVertices << "\n";

	// Re-determine topology now that spikes are gone.  DetermineTopology() uses
	// only MAT mesh structure (edges/faces/adjacency) — no ClusterType needed.
	std::cerr << "[Simplify] Re-running DetermineTopology after Phase 0...\n";
	DetermineTopology();

	// Export MAT state immediately after all spike edges have been collapsed.
	ExportOff(prefix + "_post_spike.off");

	// ── Phase 1: Main simplification (type-independent, all edges) ──────────
	current_phase = "main";
	int mainPass = 0;
	initTopoCollapseQueue();
	startPhaseLog("main", (unsigned)topo_collapse_queue.size());
	while (deleteSphereNum < threshold && numVertices > 1 && !topo_collapse_queue.empty())
	{
		++mainPass;
		unsigned attempted = 0, collapsed = 0;
		std::cerr << "[Simplify] Phase 1 pass " << mainPass << " (main): queue size = "
		          << topo_collapse_queue.size() << "  MAT vertices = " << numVertices << "\n";
		while (deleteSphereNum < threshold && numVertices > 1 && !topo_collapse_queue.empty())
		{
			EdgeInfo topEdge = topo_collapse_queue.top(); topo_collapse_queue.pop();
			unsigned eid = topEdge.edge_num;
			if (edges[eid].first && ValidVertex(edges[eid].second->vertices_.first) && ValidVertex(edges[eid].second->vertices_.second))
			{ ++attempted; if (MinCostEdgeCollapse(eid)) { ++collapsed; deleteSphereNum++; } }
		}
		printPhaseSummary("main", attempted, collapsed);
		if (collapsed == 0) break;
		initTopoCollapseQueue();
	}
	std::cerr << "[Simplify] Phase 1 done: " << mainPass << " pass(es), MAT vertices = " << numVertices << "\n";
	std::cerr << "[Simplify] edge_last_rejection has " << edge_last_rejection.size() << " entries for ExportSkeletonPLY.\n";

}

void SlabMesh::initCollapseQueue(){

	static const char* ct_names[] = {
		"T0","T1_spike","T2","T3","T4","T5","T1_non_spike",
		"MS_Unknown","MS_Sheet","MS_Seam","MS_Boundary","MS_Junction",
		"MS_Sheet_Boundary","MS_Seam_Boundary","MS_Junction_Boundary"
	};

	std::map<uint8_t, unsigned> type_counts;
	unsigned total_queued = 0;
	for (int i = 0; i < numEdges; i++)
	{
		if (edges[i].first)
		{
			EvaluateEdgeCollapseCost(i);
			edge_collapses_queue.push(EdgeInfo(i, edges[i].second->collapse_cost));

			const unsigned fir = edges[i].second->vertices_.first;
			const unsigned sec = edges[i].second->vertices_.second;
			if (vertices[fir].first) ++type_counts[static_cast<uint8_t>(vertices[fir].second->nmn_cluster_type)];
			if (vertices[sec].first) ++type_counts[static_cast<uint8_t>(vertices[sec].second->nmn_cluster_type)];
			++total_queued;
		}
	}

	std::cerr << "[initCollapseQueue] edges queued = " << total_queued << "\n";
}

void SlabMesh::initBoundaryCollapseQueue()
{
	std::map<uint8_t, unsigned> bq_type_counts;
	unsigned bq_total = 0;

	for (int i = 0; i < edges.size(); i ++)
	{
		if (edges[i].first)
		{
			unsigned fir = edges[i].second->vertices_.first;
			unsigned sec = edges[i].second->vertices_.second;

			switch(boundary_compute_scale)
			{
			case 1:
				if (!vertices[fir].second->fake_boundary_vertex || !vertices[sec].second->fake_boundary_vertex)
					continue;
				break;
			case 2:
				if (!vertices[fir].second->fake_boundary_vertex && !vertices[sec].second->fake_boundary_vertex)
					//if (vertices[fir].second->boundary_edge_vec.size() < 2 && vertices[sec].second->boundary_edge_vec.size() < 2)
					continue;
				break;
			case 3:
				if (!vertices[fir].second->fake_boundary_vertex && !vertices[sec].second->fake_boundary_vertex)
					continue;
				break;
			default:
				break;
			}

			//EvaluateEdgeCollapseCost(i);
			EvaluateEdgeHausdorffCost(i);
			boundary_edge_collapses_queue.push(EdgeInfo(i, edges[i].second->collapse_cost));

			const unsigned fir2 = edges[i].second->vertices_.first;
			const unsigned sec2 = edges[i].second->vertices_.second;
			if (vertices[fir2].first) ++bq_type_counts[static_cast<uint8_t>(vertices[fir2].second->nmn_cluster_type)];
			if (vertices[sec2].first) ++bq_type_counts[static_cast<uint8_t>(vertices[sec2].second->nmn_cluster_type)];
			++bq_total;
		}
	}

	static const char* ct_names[] = {
		"T0","T1_spike","T2","T3","T4","T5","T1_non_spike",
		"MS_Unknown","MS_Sheet","MS_Seam","MS_Boundary","MS_Junction",
		"MS_Sheet_Boundary","MS_Seam_Boundary","MS_Junction_Boundary"
	};

	std::cerr << "[initBoundaryCollapseQueue] edges queued = " << bq_total << "\n";
}

// ─────────────────────────────────────────────────────────────────────────────
// initTopoCollapseQueue  — all active edges, type-independent
// ─────────────────────────────────────────────────────────────────────────────
void SlabMesh::initTopoCollapseQueue()
{
	while (!topo_collapse_queue.empty()) topo_collapse_queue.pop();

	unsigned total = 0;
	for (unsigned i = 0; i < (unsigned)edges.size(); ++i)
	{
		if (!edges[i].first) continue;
		const unsigned fir = edges[i].second->vertices_.first;
		const unsigned sec = edges[i].second->vertices_.second;
		if (!vertices[fir].first || !vertices[sec].first) continue;
		EvaluateEdgeCollapseCost(i);
		topo_collapse_queue.push(EdgeInfo(i, edges[i].second->collapse_cost));
		++total;
	}
	std::cerr << "[initTopoCollapseQueue] queued " << total << " edges\n";
}

void SlabMesh::initSpikeCollapseQueue()
{
	// Drain any stale entries first.
	while (!spike_collapse_queue.empty())
		spike_collapse_queue.pop();

	for (unsigned i = 0; i < (unsigned)edges.size(); ++i)
	{
		if (!edges[i].first) continue;
		const unsigned fir = edges[i].second->vertices_.first;
		const unsigned sec = edges[i].second->vertices_.second;
		// Include the edge if either endpoint is a T1_spike (true spike) vertex.
		// T1_non_spike vertices are sheet boundaries and use the regular collapse queue.
		if (vertices[fir].second->nmn_cluster_type == SlabVertex::ClusterType::T1_spike ||
		    vertices[sec].second->nmn_cluster_type == SlabVertex::ClusterType::T1_spike)
		{
			// EvaluateEdgeCollapseCost(i);
			// spike_collapse_queue.push(EdgeInfo(i, edges[i].second->collapse_cost));
			spike_collapse_queue.push(EdgeInfo(i, 0.0f));// We want to collapse all spike edges regardless of cost, so we push a dummy cost of 0.0f.
		}
	}

	std::cerr << "[initSpikeCollapseQueue] spike edges queued = " << spike_collapse_queue.size() << "\n";
}

void SlabMesh::ExportOff(const std::string& path) const
{
	std::unordered_map<unsigned, unsigned> vid_to_idx;
	std::vector<unsigned> active_vids;
	for (unsigned i = 0; i < vertices.size(); ++i)
	{
		if (!vertices[i].first) continue;
		vid_to_idx[i] = (unsigned)active_vids.size();
		active_vids.push_back(i);
	}

	unsigned n_faces = 0;
	for (unsigned i = 0; i < faces.size(); ++i)
		if (faces[i].first) ++n_faces;

	std::ofstream f(path);
	if (!f) { std::cerr << "[ExportOff] Cannot open " << path << "\n"; return; }

	const double scale = pmesh ? pmesh->bb_diagonal_length : 1.0;
	f << "OFF\n";
	f << active_vids.size() << " " << n_faces << " 0\n";
	f << std::fixed << std::setprecision(10);

	for (unsigned vid : active_vids)
	{
		const auto& c = vertices[vid].second->sphere.center;
		f << c[0]*scale << " " << c[1]*scale << " " << c[2]*scale << "\n";
	}

	for (unsigned i = 0; i < faces.size(); ++i)
	{
		if (!faces[i].first) continue;
		const auto& vset = faces[i].second->vertices_;
		f << vset.size();
		for (unsigned v : vset)
			f << " " << vid_to_idx.at(v);
		f << "\n";
	}
	f.close();
	std::cout << "[ExportOff] Written to: " << path << "\n";
}

double SlabMesh::NearestPoint(Vector3d point, unsigned vid)
{
	set<unsigned> near_faces = vertices[vid].second->faces_;
	set<unsigned> near_edges = vertices[vid].second->edges_;

	double mind = DBL_MAX;
	// calculation of related faces
	for (set<unsigned>::iterator si = near_faces.begin(); si != near_faces.end(); si++)
	{
		if (!faces[*si].first)
			continue;
		SlabFace sf = *faces[*si].second;
		if (sf.valid_st == false || sf.st[0].normal == Vector3d(0., 0., 0.) || 
			sf.st[1].normal == Vector3d(0., 0., 0.))
			continue;

		Vector3d v[3];
		Vector3d tfp;
		double td;
		for (int i = 0; i < 2; i++)
		{
			v[0] = sf.st[i].v[0];
			v[1] = sf.st[i].v[1];
			v[2] = sf.st[i].v[2];
			ProjectOntoTriangle(point, v[0], v[1], v[2], tfp, td);

			if(td < mind)	mind = td;
		}
	}

	// calculation of related edges
	for (set<unsigned>::iterator si = near_edges.begin(); si != near_edges.end(); si++)
	{
		if (!edges[*si].first)
			continue;
		SlabEdge se = *edges[*si].second;
		if (se.valid_cone == false)
			continue;

		SlabVertex v[2];
		Vector3d tfp;
		double td;
		double tr;
		v[0] = *vertices[se.vertices_.first].second;
		v[1] = *vertices[se.vertices_.second].second;
		Vector3d v0 = v[0].sphere.center;
		Vector3d v1 = v[1].sphere.center;
		double t((point-v0).Dot(v1-v0) / (v1-v0).SquaredLength());
		if( (t >= 0.0) && (t <= 1.0) )
		{
			tfp = (1.0-t)*v0 + t*v1;
			td = (point-tfp).Length();
			tr = (1.0-t)*v[0].sphere.radius + t*v[1].sphere.radius;

			if (abs(td - tr) < mind) mind = abs(td - tr);
		}
	}

	return mind;
}

// simple method, do not add any plane to preserve the boundary 
void SlabMesh::PreservBoundaryMethodOne()
{
	// �����е�boundary_edge���ϱ߽籣��
	for (unsigned i = 0; i < edges.size(); i++)
	{
		if (!edges[i].first || !edges[i].second->fake_boundary_edge)
			continue;

		SlabEdge se = *(edges[i].second);
		if (edges[i].second->faces_.size() == 0)
		{
			unsigned ver_index[2];
			ver_index[0] = se.vertices_.first;
			ver_index[1] = se.vertices_.second;
			Sphere ver[3];
			ver[0]= vertices[ver_index[0]].second->sphere;
			ver[1] = vertices[ver_index[1]].second->sphere;
			ver[2].center = (ver[0].center + ver[1].center);
			ver[2].radius = (ver[0].radius + ver[1].radius) / 2.0;

			SimpleTriangle st[2];
			Wm4::Vector3d pos[3];
			double radius[3];
			for(int count = 0; count < 3; count++)
			{
				pos[count] = ver[count].center;
				radius[count] = ver[count].radius;
			}
			if(TriangleFromThreeSpheres(pos[0],radius[0],pos[1],radius[1],pos[2],radius[2],st[0],st[1]))
			{
				if (st[0].normal == Vector3d(0., 0., 0.) || st[1].normal == Vector3d(0., 0., 0.))
					continue;

				// �ӵ�һ��slab�е�����ƽ��
				Vector4d normal1(st[0].normal.X(), st[0].normal.Y(), st[0].normal.Z(), 1.0);
				Vector4d normal2(st[1].normal.X(), st[1].normal.Y(), st[1].normal.Z(), 1.0);
				// compute the matrix of A
				Matrix4d temp_A1, temp_A2;
				temp_A1.MakeTensorProduct(normal1, normal1);
				temp_A2.MakeTensorProduct(normal2, normal2);
				temp_A1 *= 2.0;
				temp_A2 *= 2.0;
				for (int i = 0; i < 2; i++)
				{
					Vector4d C1(ver[i].center.X(), ver[i].center.Y(), ver[i].center.Z(), ver[i].radius);
					// compute the matrix of b
					double normal_mul_point1 = normal1.Dot(C1);
					double normal_mul_point2 = normal2.Dot(C1);
					Wm4::Vector4d temp_b1 = normal1 * 2 * normal_mul_point1;
					Wm4::Vector4d temp_b2 = normal2 * 2 * normal_mul_point2;
					//compute c
					double temp_c1 = normal_mul_point1 * normal_mul_point1;
					double temp_c2 = normal_mul_point2 * normal_mul_point2;
					vertices[ver_index[i]].second->add_A += temp_A1;
					vertices[ver_index[i]].second->add_A += temp_A2;
					vertices[ver_index[i]].second->add_b += temp_b1;
					vertices[ver_index[i]].second->add_b += temp_b2;
					vertices[ver_index[i]].second->add_c += temp_c1;
					vertices[ver_index[i]].second->add_c += temp_c2;
				}

				// �ӵڶ���slab�е�����ƽ��
				Vector3d ver1_to_ver2 = ver[0].center - ver[1].center;
				Vector3d t1 = ver1_to_ver2.Cross(st[0].normal);
				Vector3d t2 = ver1_to_ver2.Cross(st[1].normal);
				Vector4d tnormal1(t1.X(), t1.Y(), t1.Z(), 1.0);
				Vector4d tnormal2(t2.X(), t2.Y(), t2.Z(), 1.0);
				// compute the matrix of A
				temp_A1.MakeTensorProduct(tnormal1, tnormal1);
				temp_A2.MakeTensorProduct(tnormal2, tnormal2);
				temp_A1 *= 2.0;
				temp_A2 *= 2.0;
				for (int i = 0; i < 2; i++)
				{
					Vector4d C1(ver[i].center.X(), ver[i].center.Y(), ver[i].center.Z(), ver[i].radius);
					// compute the matrix of b
					double normal_mul_point1 = tnormal1.Dot(C1);
					double normal_mul_point2 = tnormal2.Dot(C1);
					Wm4::Vector4d temp_b1 = tnormal1 * 2 * normal_mul_point1;
					Wm4::Vector4d temp_b2 = tnormal2 * 2 * normal_mul_point2;
					//compute c
					double temp_c1 = normal_mul_point1 * normal_mul_point1;
					double temp_c2 = normal_mul_point2 * normal_mul_point2;
					vertices[ver_index[i]].second->add_A += temp_A1;
					vertices[ver_index[i]].second->add_A += temp_A2;
					vertices[ver_index[i]].second->add_b += temp_b1;
					vertices[ver_index[i]].second->add_b += temp_b2;
					vertices[ver_index[i]].second->add_c += temp_c1;
					vertices[ver_index[i]].second->add_c += temp_c2;
				}
			}

			//// ���ڱ�¶�����ĵ��ټ�һ������ƽ��
			//if (vertices[ver_index[0]].second->edges_.size() == 1)
			//{
			//	Vector3d ver1_to_ver2 = ver[0].center - ver[1].center;
			//	Vector4d normal1(ver1_to_ver2.X(), ver1_to_ver2.Y(), ver1_to_ver2.Z(), 1.0);
			//	// compute the matrix of A
			//	Matrix4d temp_A1;
			//	temp_A1.MakeTensorProduct(normal1, normal1);
			//	temp_A1 *= 2.0;
			//	Vector4d C1(ver[0].center.X(), ver[0].center.Y(), ver[0].center.Z(), ver[0].radius);
			//	// compute the matrix of b
			//	double normal_mul_point1 = normal1.Dot(C1);
			//	Wm4::Vector4d temp_b1 = normal1 * 2 * normal_mul_point1;
			//	//compute c
			//	double temp_c1 = normal_mul_point1 * normal_mul_point1;
			//	vertices[ver_index[0]].second->add_A += temp_A1 * 100.0;
			//	vertices[ver_index[0]].second->add_b += temp_b1 * 100.0;
			//	vertices[ver_index[0]].second->add_c += temp_c1 * 100.0;
			//}
			//if (vertices[ver_index[1]].second->edges_.size() == 1)
			//{
			//	Vector3d ver2_to_ver1 = ver[1].center - ver[0].center;
			//	Vector4d normal1(ver2_to_ver1.X(), ver2_to_ver1.Y(), ver2_to_ver1.Z(), 1.0);
			//	// compute the matrix of A
			//	Matrix4d temp_A1;
			//	temp_A1.MakeTensorProduct(normal1, normal1);
			//	temp_A1 *= 2.0;
			//	Vector4d C1(ver[1].center.X(), ver[1].center.Y(), ver[1].center.Z(), ver[1].radius);
			//	// compute the matrix of b
			//	double normal_mul_point1 = normal1.Dot(C1);
			//	Wm4::Vector4d temp_b1 = normal1 * 2 * normal_mul_point1;
			//	//compute c
			//	double temp_c1 = normal_mul_point1 * normal_mul_point1;
			//	vertices[ver_index[1]].second->add_A += temp_A1 * 100.0;
			//	vertices[ver_index[1]].second->add_b += temp_b1 * 100.0;
			//	vertices[ver_index[1]].second->add_c += temp_c1 * 100.0;
			//}
		}
	}
}

// just add one normal to all the Boundary Vertex Sphere with two Boundary edges.
void SlabMesh::PreservBoundaryMethodTwo()
{
	// for each boundary vertex sphere, add a normal
	for (unsigned i = 0; i < vertices.size(); i++)
	{
		if (vertices[i].first && vertices[i].second->boundary_edge_vec.size() == 2)
		{
			Vector3d add_normal(0.0, 0.0, 0.0);
			int valid_edge = 0;
			for (set<unsigned>::iterator vi = vertices[i].second->boundary_edge_vec.begin(); vi != vertices[i].second->boundary_edge_vec.end(); vi++)
			{
				if (edges[*vi].second->faces_.size() == 0)
					continue;
				unsigned face_num = *(edges[*vi].second->faces_.begin());			
				SlabFace sf = *faces[face_num].second;
				if (sf.valid_st == false || sf.st[0].normal == Vector3d(0., 0., 0.) || 
					sf.st[1].normal == Vector3d(0., 0., 0.))
					continue;

				valid_edge++;
				unsigned ver_index = edges[*vi].second->vertices_.first == i ? edges[*vi].second->vertices_.second 
					: edges[*vi].second->vertices_.first;
				Vector3d temp_norm = vertices[i].second->sphere.center - vertices[ver_index].second->sphere.center;
				temp_norm.Normalize();
				add_normal += temp_norm;
			}
			if (valid_edge == 0)
				continue;

			add_normal /= valid_edge;
			add_normal.Normalize();

			vertices[i].second->boundVec = add_normal;

			// ȷ�����ӵ����Ȩ�ش�С
			//Vector3d boudary_vec[2];
			//for (int index = 0; index < 2; index++)
			//{
			//	set<unsigned>::iterator vi = vertices[i].second->boundary_edge_vec.begin() + index;
			//	unsigned ver_index = edges[*vi].second->vertices_.first == i ? edges[*vi].second->vertices_.second 
			//		: edges[*vi].second->vertices_.first;
			//	Vector3d temp_norm = vertices[i].second->sphere.center - vertices[ver_index].second->sphere.center;
			//	boudary_vec[index] = temp_norm;
			//}
			//vertices[i].second->collaspe_weight = sin(VectorAngle(boudary_vec[0], boudary_vec[1]));

			// �ж�����߽��������͹�����ĵ㻹�ǰ���ȥ�ĵ�
			bool boundary_vertex = false;
			for (auto si = vertices[i].second->edges_.begin(); si != vertices[i].second->edges_.end(); si++)
			{
				unsigned temp_ind = edges[*si].second->vertices_.first == i ? edges[*si].second->vertices_.second 
					: edges[*si].second->vertices_.first;
				Vector3d temp_vec = vertices[temp_ind].second->sphere.center - vertices[i].second->sphere.center;
				double temp_angle = acos(temp_vec.Dot(add_normal) / temp_vec.Length());
				boundary_vertex = temp_angle < Wm4::Math<double>::PI / 2.0 ? true : false;
				if (boundary_vertex == true)
				{
					vertices[i].second->boundary_vertex = true;
					break;
				}
			}

			if (boundary_vertex == false)
			{
				Vector4d normal(add_normal.X(), add_normal.Y(), add_normal.Z(), 1.0);
				Matrix4d temp_A;
				temp_A.MakeTensorProduct(normal, normal);
				temp_A *= 2.0;
				Sphere ve = vertices[i].second->sphere;
				Vector4d C1(ve.center.X(), ve.center.Y(), ve.center.Z(), ve.radius);
				double normal_mul_point = normal.Dot(C1);
				Wm4::Vector4d temp_b = normal * 2 * normal_mul_point;
				double temp_c = normal_mul_point * normal_mul_point;

				vertices[i].second->add_A += temp_A;
				vertices[i].second->add_b += temp_b;
				vertices[i].second->add_c += temp_c;

				vertices[i].second->related_face ++;
			}
			else
			{
				add_normal = add_normal * (-1.0);
				Vector4d normal(add_normal.X(), add_normal.Y(), add_normal.Z(), 1.0);
				Matrix4d temp_A;
				temp_A.MakeTensorProduct(normal, normal);
				temp_A *= 2.0;
				Sphere ve = vertices[i].second->sphere;
				Vector4d C1(ve.center.X(), ve.center.Y(), ve.center.Z(), ve.radius);
				double normal_mul_point = normal.Dot(C1);
				Wm4::Vector4d temp_b = normal * 2 * normal_mul_point;
				double temp_c = normal_mul_point * normal_mul_point;

				vertices[i].second->add_A += temp_A;
				vertices[i].second->add_b += temp_b;
				vertices[i].second->add_c += temp_c;

				vertices[i].second->related_face ++;
			}
		}
	}
}

// add a plane to the Possible Spike Vertex, and add a slab to the Possible Boundary Vertex.
void SlabMesh::PreservBoundaryMethodThree()
{
	// for each boundary vertex sphere, add a normal
	for (unsigned i = 0; i < vertices.size(); i++)
	{
		if (vertices[i].first && vertices[i].second->boundary_edge_vec.size() == 2)
		{
			Vector3d add_normal(0.0, 0.0, 0.0);
			int valid_edge = 0;
			for (set<unsigned>::iterator vi = vertices[i].second->boundary_edge_vec.begin(); vi != vertices[i].second->boundary_edge_vec.end(); vi++)
			{
				if (edges[*vi].second->faces_.size() == 0)
					continue;
				unsigned face_num = *(edges[*vi].second->faces_.begin());			
				SlabFace sf = *faces[face_num].second;
				if (sf.valid_st == false || sf.st[0].normal == Vector3d(0., 0., 0.) || 
					sf.st[1].normal == Vector3d(0., 0., 0.))
					continue;

				valid_edge++;
				unsigned ver_index = edges[*vi].second->vertices_.first == i ? edges[*vi].second->vertices_.second 
					: edges[*vi].second->vertices_.first;
				Vector3d temp_norm = vertices[i].second->sphere.center - vertices[ver_index].second->sphere.center;
				temp_norm.Normalize();
				add_normal += temp_norm;
			}
			if (valid_edge == 0)
				continue;

			add_normal /= valid_edge;
			add_normal.Normalize();

			vertices[i].second->boundVec = add_normal;

			// ȷ�����ӵ����Ȩ�ش�С
			//Vector3d boudary_vec[2];
			//for (int index = 0; index < 2; index++)
			//{
			//	auto vi = vertices[i].second->boundary_edge_vec.begin() + index;
			//	unsigned ver_index = edges[*vi].second->vertices_.first == i ? edges[*vi].second->vertices_.second 
			//		: edges[*vi].second->vertices_.first;
			//	Vector3d temp_norm = vertices[i].second->sphere.center - vertices[ver_index].second->sphere.center;
			//	boudary_vec[index] = temp_norm;
			//}
			//vertices[i].second->collaspe_weight = sin(VectorAngle(boudary_vec[0], boudary_vec[1]));

			// �ж�����߽��������͹�����ĵ㻹�ǰ���ȥ�ĵ�
			bool boundary_vertex = false;
			for (auto si = vertices[i].second->edges_.begin(); si != vertices[i].second->edges_.end(); si++)
			{
				unsigned temp_ind = edges[*si].second->vertices_.first == i ? edges[*si].second->vertices_.second 
					: edges[*si].second->vertices_.first;
				Vector3d temp_vec = vertices[temp_ind].second->sphere.center - vertices[i].second->sphere.center;
				double temp_angle = acos(temp_vec.Dot(add_normal) / temp_vec.Length());
				boundary_vertex = temp_angle < Wm4::Math<double>::PI / 2.0 ? true : false;
				if (boundary_vertex == true)
				{
					vertices[i].second->boundary_vertex = true;
					break;
				}
			}

			if (boundary_vertex == false)
			{
				Vector4d normal(add_normal.X(), add_normal.Y(), add_normal.Z(), 1.0);
				Matrix4d temp_A;
				temp_A.MakeTensorProduct(normal, normal);
				temp_A *= 2.0;
				Sphere ve = vertices[i].second->sphere;
				Vector4d C1(ve.center.X(), ve.center.Y(), ve.center.Z(), ve.radius);
				double normal_mul_point = normal.Dot(C1);
				Wm4::Vector4d temp_b = normal * 2 * normal_mul_point;
				double temp_c = normal_mul_point * normal_mul_point;

				vertices[i].second->add_A += temp_A;
				vertices[i].second->add_b += temp_b;
				vertices[i].second->add_c += temp_c;

				vertices[i].second->related_face ++;
			}
		}
	}

	// �����е�boundary_edge���ϱ߽籣��
	for (unsigned i = 0; i < edges.size(); i++)
	{
		if (edges[i].first && edges[i].second->faces_.size() == 1)
		{
			unsigned face_num = *(edges[i].second->faces_.begin());

			unsigned ver_index[2];
			ver_index[0] = edges[i].second->vertices_.first;
			ver_index[1] = edges[i].second->vertices_.second;
			Sphere ver[2];
			ver[0]= vertices[ver_index[0]].second->sphere;
			ver[1] = vertices[ver_index[1]].second->sphere;
			Vector3d ver1_to_ver2 = ver[0].center - ver[1].center;

			SlabFace sf = *faces[face_num].second;
			if (sf.valid_st == false || sf.st[0].normal == Vector3d(0., 0., 0.) || 
				sf.st[1].normal == Vector3d(0., 0., 0.))
				continue;

			Vector3d temp_normal1(sf.st[0].normal.X(), sf.st[0].normal.Y(), sf.st[0].normal.Z());
			Vector3d temp_normal2(sf.st[1].normal.X(), sf.st[1].normal.Y(), sf.st[1].normal.Z());
			Vector3d t1 = ver1_to_ver2.Cross(temp_normal1);
			Vector3d t2 = ver1_to_ver2.Cross(temp_normal2);

			// ����boundary_edge����һ��slab���б߽籣��
			Vector4d normal1(t1.X(), t1.Y(), t1.Z(), 1.0);
			Vector4d normal2(t2.X(), t2.Y(), t2.Z(), 1.0);

			// compute the matrix of A
			Matrix4d temp_A1, temp_A2;
			temp_A1.MakeTensorProduct(normal1, normal1);
			temp_A2.MakeTensorProduct(normal2, normal2);
			temp_A1 *= 2.0;
			temp_A2 *= 2.0;

			for (int i = 0; i < 2; i++)
			{	
				//if (vertices[ver_index[i]].second->boundary_vertex == false)
				//	continue;

				Vector4d C1(ver[i].center.X(), ver[i].center.Y(), ver[i].center.Z(), ver[i].radius);

				// compute the matrix of b
				double normal_mul_point1 = normal1.Dot(C1);
				double normal_mul_point2 = normal2.Dot(C1);
				Wm4::Vector4d temp_b1 = normal1 * 2 * normal_mul_point1;
				Wm4::Vector4d temp_b2 = normal2 * 2 * normal_mul_point2;

				//compute c
				double temp_c1 = normal_mul_point1 * normal_mul_point1;
				double temp_c2 = normal_mul_point2 * normal_mul_point2;

				vertices[ver_index[i]].second->add_A += temp_A1;
				vertices[ver_index[i]].second->add_A += temp_A2;
				vertices[ver_index[i]].second->add_b += temp_b1;
				vertices[ver_index[i]].second->add_b += temp_b2;
				vertices[ver_index[i]].second->add_c += temp_c1;
				vertices[ver_index[i]].second->add_c += temp_c2;
			}
		}
	}
}

void SlabMesh::PreservBoundaryMethodFour()
{
	// �����е�boundary_edge���ϱ߽籣��
	for (unsigned i = 0; i < edges.size(); i++)
	{
		if (!edges[i].first || !edges[i].second->fake_boundary_edge)
			continue;

		SlabEdge se = *(edges[i].second);
		if (se.faces_.size() == 1)
		{
			unsigned face_num = *(se.faces_.begin());

			set<unsigned> sv = faces[face_num].second->vertices_;
			Vector3d face_normal = faces[face_num].second->normal;

			unsigned ver_index[3];
			ver_index[0] = se.vertices_.first;
			sv.erase(ver_index[0]);
			ver_index[1] = se.vertices_.second;
			sv.erase(ver_index[1]);
			ver_index[2] = *(sv.begin());
			Sphere ver[2];
			ver[0]= vertices[ver_index[0]].second->sphere;
			ver[1] = vertices[ver_index[1]].second->sphere;

			Vector3d v1v2 = vertices[ver_index[0]].second->sphere.center - vertices[ver_index[1]].second->sphere.center;
			Vector3d v1v3 = vertices[ver_index[2]].second->sphere.center - vertices[ver_index[0]].second->sphere.center;
			Vector3d temp_nor = face_normal.Cross(v1v2);

			double temp_angle = acos(temp_nor.Dot(v1v3) / temp_nor.Length() / v1v3.Length());
			bool dir = temp_angle > Wm4::Math<double>::PI / 2.0 ? true : false;

			if (dir == false)
				temp_nor *= -1;

			// ����boundary_edge����һ��ƽ����б߽籣��
			Vector4d normal1(temp_nor.X(), temp_nor.Y(), temp_nor.Z(), 1.0);
			// compute the matrix of A
			Matrix4d temp_A1;
			temp_A1.MakeTensorProduct(normal1, normal1);
			temp_A1 *= 2.0;

			// �Բ�ͬratio�ı߽���ӳ�䴦����С��0.2�Ĳ���������(0.2,1)ӳ�䵽(2, 10)
			double ratio = GetRatioHyperbolicEuclid(i);
			double w1 = 1.0;
			//w1 = 0.02 * ratio * ratio * ratio * ratio * ratio * ratio;
			w1 =  0.1 * ratio * ratio;
			//w1 = 3 * ratio * ratio;
			//w1 = 1 * ratio * ratio * ratio;
			//if (ratio >= 0.3)
			//{
			//	w1 = (ratio - 0.3) / 0.7 * (3 - 1) + 1;
			//}

			for (int i = 0; i < 2; i++)
			{	
				Vector4d C1(ver[i].center.X(), ver[i].center.Y(), ver[i].center.Z(), ver[i].radius);

				// compute the matrix of b
				double normal_mul_point1 = normal1.Dot(C1);
				Wm4::Vector4d temp_b1 = normal1 * 2 * normal_mul_point1;

				//compute c
				double temp_c1 = normal_mul_point1 * normal_mul_point1;

				vertices[ver_index[i]].second->add_A += temp_A1 * bound_weight * w1;
				vertices[ver_index[i]].second->add_b += temp_b1 * bound_weight * w1;
				vertices[ver_index[i]].second->add_c += temp_c1 * bound_weight * w1;
			}
		}else if (edges[i].second->faces_.size() == 0)
		{
			unsigned ver_index[2];
			ver_index[0] = se.vertices_.first;
			ver_index[1] = se.vertices_.second;
			Sphere ver[3];
			ver[0]= vertices[ver_index[0]].second->sphere;
			ver[1] = vertices[ver_index[1]].second->sphere;
			ver[2].center = (ver[0].center + ver[1].center);
			ver[2].radius = (ver[0].radius + ver[1].radius) / 2.0;

			SimpleTriangle st[2];
			Wm4::Vector3d pos[3];
			double radius[3];
			for(int count = 0; count < 3; count++)
			{
				pos[count] = ver[count].center;
				radius[count] = ver[count].radius;
			}
			if(TriangleFromThreeSpheres(pos[0],radius[0],pos[1],radius[1],pos[2],radius[2],st[0],st[1]))
			{
				if (st[0].normal == Vector3d(0., 0., 0.) || st[1].normal == Vector3d(0., 0., 0.))
					continue;

				// �ӵ�һ��slab�е�����ƽ��
				Vector4d normal1(st[0].normal.X(), st[0].normal.Y(), st[0].normal.Z(), 1.0);
				Vector4d normal2(st[1].normal.X(), st[1].normal.Y(), st[1].normal.Z(), 1.0);
				// compute the matrix of A
				Matrix4d temp_A1, temp_A2;
				temp_A1.MakeTensorProduct(normal1, normal1);
				temp_A2.MakeTensorProduct(normal2, normal2);
				temp_A1 *= 2.0;
				temp_A2 *= 2.0;
				for (int i = 0; i < 2; i++)
				{
					Vector4d C1(ver[i].center.X(), ver[i].center.Y(), ver[i].center.Z(), ver[i].radius);
					// compute the matrix of b
					double normal_mul_point1 = normal1.Dot(C1);
					double normal_mul_point2 = normal2.Dot(C1);
					Wm4::Vector4d temp_b1 = normal1 * 2 * normal_mul_point1;
					Wm4::Vector4d temp_b2 = normal2 * 2 * normal_mul_point2;
					//compute c
					double temp_c1 = normal_mul_point1 * normal_mul_point1;
					double temp_c2 = normal_mul_point2 * normal_mul_point2;
					vertices[ver_index[i]].second->add_A += temp_A1;
					vertices[ver_index[i]].second->add_A += temp_A2;
					vertices[ver_index[i]].second->add_b += temp_b1;
					vertices[ver_index[i]].second->add_b += temp_b2;
					vertices[ver_index[i]].second->add_c += temp_c1;
					vertices[ver_index[i]].second->add_c += temp_c2;
				}

				// �ӵڶ���slab�е�����ƽ��
				Vector3d ver1_to_ver2 = ver[0].center - ver[1].center;
				Vector3d t1 = ver1_to_ver2.Cross(st[0].normal);
				Vector3d t2 = ver1_to_ver2.Cross(st[1].normal);
				Vector4d tnormal1(t1.X(), t1.Y(), t1.Z(), 1.0);
				Vector4d tnormal2(t2.X(), t2.Y(), t2.Z(), 1.0);
				// compute the matrix of A
				temp_A1.MakeTensorProduct(tnormal1, tnormal1);
				temp_A2.MakeTensorProduct(tnormal2, tnormal2);
				temp_A1 *= 2.0;
				temp_A2 *= 2.0;
				for (int i = 0; i < 2; i++)
				{
					Vector4d C1(ver[i].center.X(), ver[i].center.Y(), ver[i].center.Z(), ver[i].radius);
					// compute the matrix of b
					double normal_mul_point1 = tnormal1.Dot(C1);
					double normal_mul_point2 = tnormal2.Dot(C1);
					Wm4::Vector4d temp_b1 = tnormal1 * 2 * normal_mul_point1;
					Wm4::Vector4d temp_b2 = tnormal2 * 2 * normal_mul_point2;
					//compute c
					double temp_c1 = normal_mul_point1 * normal_mul_point1;
					double temp_c2 = normal_mul_point2 * normal_mul_point2;
					vertices[ver_index[i]].second->add_A += temp_A1;
					vertices[ver_index[i]].second->add_A += temp_A2;
					vertices[ver_index[i]].second->add_b += temp_b1;
					vertices[ver_index[i]].second->add_b += temp_b2;
					vertices[ver_index[i]].second->add_c += temp_c1;
					vertices[ver_index[i]].second->add_c += temp_c2;
				}
			}

			double ratio = GetRatioHyperbolicEuclid(i);
			double w1 = 1.0;
			//w1 = 0.02 * ratio * ratio * ratio * ratio * ratio * ratio;
			w1 = 0.1 * ratio * ratio;
			//w1 = 1 * ratio * ratio * ratio * ratio * ratio * ratio;
			//w1 = 1 * ratio * ratio * ratio;
			// ���ڱ�¶�����ĵ��ټ�һ������ƽ��
			if (vertices[ver_index[0]].second->edges_.size() == 1)
			{

				Vector3d ver1_to_ver2 = ver[0].center - ver[1].center;
				Vector4d normal1(ver1_to_ver2.X(), ver1_to_ver2.Y(), ver1_to_ver2.Z(), 1.0);
				// compute the matrix of A
				Matrix4d temp_A1;
				temp_A1.MakeTensorProduct(normal1, normal1);
				temp_A1 *= 2.0;
				Vector4d C1(ver[0].center.X(), ver[0].center.Y(), ver[0].center.Z(), ver[0].radius);
				// compute the matrix of b
				double normal_mul_point1 = normal1.Dot(C1);
				Wm4::Vector4d temp_b1 = normal1 * 2 * normal_mul_point1;
				//compute c
				double temp_c1 = normal_mul_point1 * normal_mul_point1;
				//vertices[ver_index[0]].second->add_A += (temp_A1 * 100.0);
				//vertices[ver_index[0]].second->add_b += (temp_b1 * 100.0);
				//vertices[ver_index[0]].second->add_c += (temp_c1 * 100.0);
				vertices[ver_index[0]].second->add_A += temp_A1 * bound_weight * w1;
				vertices[ver_index[0]].second->add_b += temp_b1 * bound_weight * w1;
				vertices[ver_index[0]].second->add_c += temp_c1 * bound_weight * w1;
			}
			if (vertices[ver_index[1]].second->edges_.size() == 1)
			{
				Vector3d ver2_to_ver1 = ver[1].center - ver[0].center;
				Vector4d normal1(ver2_to_ver1.X(), ver2_to_ver1.Y(), ver2_to_ver1.Z(), 1.0);
				// compute the matrix of A
				Matrix4d temp_A1;
				temp_A1.MakeTensorProduct(normal1, normal1);
				temp_A1 *= 2.0;
				Vector4d C1(ver[1].center.X(), ver[1].center.Y(), ver[1].center.Z(), ver[1].radius);
				// compute the matrix of b
				double normal_mul_point1 = normal1.Dot(C1);
				Wm4::Vector4d temp_b1 = normal1 * 2 * normal_mul_point1;
				//compute c
				double temp_c1 = normal_mul_point1 * normal_mul_point1;
				//vertices[ver_index[1]].second->add_A += (temp_A1 * 100.0);
				//vertices[ver_index[1]].second->add_b += (temp_b1 * 100.0);
				//vertices[ver_index[1]].second->add_c += (temp_c1 * 100.0);
				vertices[ver_index[1]].second->add_A += temp_A1 * bound_weight * w1;
				vertices[ver_index[1]].second->add_b += temp_b1 * bound_weight * w1;
				vertices[ver_index[1]].second->add_c += temp_c1 * bound_weight * w1;
			}
		}
	} 
}

void SlabMesh::clear()
{
	for (unsigned i = 0; i < vertices.size(); i++)
	{
		if (vertices[i].first)
		{
			vertices[i].second->slab_A.MakeZero();
			vertices[i].second->add_A.MakeZero();
			vertices[i].second->slab_b = Wm4::Vector4<double>::ZERO;
			vertices[i].second->add_b = Wm4::Vector4<double>::ZERO;
			vertices[i].second->slab_c = 0.0;
			vertices[i].second->add_c = 0.0;
			vertices[i].second->fake_boundary_vertex = false;
			vertices[i].second->boundary_vertex = false;
			vertices[i].second->boundary_edge_vec.clear();
			vertices[i].second->mean_square_error = 0.0;
			vertices[i].second->related_face = 0;
		}
	}

	for (unsigned i = 0; i < edges.size(); i++)
	{
		if (edges[i].first)
		{
			edges[i].second->slab_A.MakeZero();
			edges[i].second->slab_b = Wm4::Vector4<double>::ZERO;
			edges[i].second->slab_c = 0.0;
			edges[i].second->fake_boundary_edge = false;
			edges[i].second->non_manifold_edge = false;
		}
	}

	max_mean_squre_error = 0.0;

	DistinguishVertexType();
}

void SlabMesh::RecomputerVertexType()
{
	for (unsigned i = 0; i < vertices.size(); i++)
	{
		if (vertices[i].first)
		{
			vertices[i].second->fake_boundary_vertex = false;
			vertices[i].second->boundary_vertex = false;
			vertices[i].second->boundary_edge_vec.clear();
		}
	}

	for (unsigned i = 0; i < edges.size(); i++)
	{
		if (edges[i].first)
		{
			edges[i].second->fake_boundary_edge = false;
			edges[i].second->non_manifold_edge = false;
		}
	}

	DistinguishVertexType();
}

void SlabMesh::computebb()
{
	m_min[0] = 1e20;
	m_min[1] = 1e20;
	m_min[2] = 1e20;
	m_max[0] = -1e20;
	m_max[1] = -1e20;
	m_max[2] = -1e20;

	for (unsigned i = 0; i < vertices.size(); i++)
	{
		if (!vertices[i].first)
			continue;

		Vector3d ver = vertices[i].second->sphere.center;
		if (ver[0] < m_min[0])
			m_min[0] = ver[0]; 
		if (ver[1] < m_min[1])
			m_min[1] = ver[1]; 
		if (ver[2] < m_min[2])
			m_min[2] = ver[2];

		if (ver[0] > m_max[0])
			m_max[0] = ver[0]; 
		if (ver[1] > m_max[1])
			m_max[1] = ver[1]; 
		if (ver[2] > m_max[2])
			m_max[2] = ver[2];
	}
}

void SlabMesh::GetEnvelopeSet(const Vector4d & lamder, const set<unsigned> & neighbor_v, const set< std::set<unsigned> > & adj_faces, vector<Sphere> & sph_vec, vector<Cone> & con_vec, vector<SimpleTriangle> & st_vec)
{
	Vector3d ps(lamder.X(), lamder.Y(), lamder.Z());
	double rs = lamder.W();

	sph_vec.push_back(Sphere(ps,rs));

	for(std::set<unsigned>::iterator si = neighbor_v.begin(); si != neighbor_v.end(); si ++)
	{
		//sph_vec.push_back(Sphere(vertices[*si].second->pos, vertices[*si].second->radius));
		Cone newc = Cone(ps,rs,vertices[*si].second->sphere.center, vertices[*si].second->sphere.radius);
		if(newc.type != 1)
			con_vec.push_back(newc);
	}

	for(std::set< std::set<unsigned> >::iterator si = adj_faces.begin(); si != adj_faces.end(); si ++)
	{
		Vector3d cen[2];
		double rad[2];
		unsigned count = 0;
		for(std::set<unsigned>::iterator si2 = (*si).begin(); si2 != (*si).end(); si2 ++, count ++)
		{
			cen[count] = vertices[*si2].second->sphere.center;
			rad[count] = vertices[*si2].second->sphere.radius;
		}

		SimpleTriangle st[2];
		if(TriangleFromThreeSpheres(cen[0], rad[0], cen[1], rad[1], ps, rs, st[0], st[1]))
		{
			st_vec.push_back(st[0]);
			st_vec.push_back(st[1]);
		}

	}

	return;
}

double SlabMesh::EvaluateVertexDistanceErrorEnvelope(Vector4d & lamdar, set<unsigned> & neighbor_vertices, set< set<unsigned> > & neighbor_faces, set<unsigned> & bplist)
{
	bool valid_cone = true;
	vector<Sphere> sph_vec;
	vector<Cone> con_vec;
	vector<SimpleTriangle> st_vec;
	GetEnvelopeSet(lamdar, neighbor_vertices, neighbor_faces, sph_vec, con_vec, st_vec);

	double maxerror(0.0);
	for(set<unsigned>::iterator si = bplist.begin(); si != bplist.end(); si ++)
	{
		Vector3d p(pmesh->pVertexList[*si]->point()[0], pmesh->pVertexList[*si]->point()[1], pmesh->pVertexList[*si]->point()[2]);
		double tempdist;
		Vector3d tempfp;

		double mindist(1e20);
		for(unsigned i = 0; i < sph_vec.size(); i ++)
		{
			sph_vec[i].ProjectOntoSphere(p,tempfp,tempdist);
			tempdist = fabs(tempdist);
			mindist = min(tempdist,mindist);
		}

		tempdist = abs((p - Wm4::Vector3d(lamdar.X(), lamdar.Y(), lamdar.Z())).Length() - lamdar.W());
		mindist = min(tempdist,mindist);

		for(unsigned i = 0; i < con_vec.size(); i ++)
		{
			con_vec[i].ProjectOntoCone(p,tempfp,tempdist);
			//if(tempdist < -2.*error_threshold)
			//	valid_cone = false;
			if(tempdist < 0)
				tempdist = -tempdist;
			mindist = min(tempdist,mindist);
		}

		for(unsigned i = 0; i < st_vec.size(); i ++)
		{
			st_vec[i].ProjectOntoSimpleTriangle(p,tempfp,tempdist);
			mindist = min(tempdist, mindist);
		}
		maxerror = max(maxerror, mindist);
	}
	//vertices[vid].second->v_evaluated_distance_error_envelope = maxerror;
	return maxerror;
}

double SlabMesh::GetHyperbolicLength(unsigned eid)
{
	double hyperbolic_weight;
	unsigned v1 = edges[eid].second->vertices_.first;
	unsigned v2 = edges[eid].second->vertices_.second;
	double edge_length = Vector3d(vertices[v1].second->sphere.center - vertices[v2].second->sphere.center).Length();
	double r1 = vertices[v1].second->sphere.radius;
	double r2 = vertices[v2].second->sphere.radius;
	if (r1 <= r2)
		hyperbolic_weight = edge_length - (r2 - r1);
	else
		hyperbolic_weight = edge_length - (r1 - r2);
	hyperbolic_weight = max(hyperbolic_weight, 0.0); 
	return hyperbolic_weight;
}

double SlabMesh::GetRatioHyperbolicEuclid(unsigned eid)
{
	double hyperbolic_distance;
	unsigned v1 = edges[eid].second->vertices_.first;
	unsigned v2 = edges[eid].second->vertices_.second;
	double edge_length = Vector3d(vertices[v1].second->sphere.center - vertices[v2].second->sphere.center).Length();
	double r1 = vertices[v1].second->sphere.radius;
	double r2 = vertices[v2].second->sphere.radius;
	if (r1 <= r2)
		hyperbolic_distance = edge_length - (r2 - r1);
	else
		hyperbolic_distance = edge_length - (r1 - r2);
	hyperbolic_distance = max(hyperbolic_distance, 0.0); 

	if (edge_length == 0.0)
		return 0.0;

	return hyperbolic_distance / edge_length;
}

void SlabMesh::ExportSimplifyResult()
{
	//std::ofstream f_result_out;
	//f_result_out.open("Result.txt", ios::app);

	//f_result_out << simplified_boundary_edges << "\t" << simplified_inside_edges << "\t" << maxhausdorff_distance << endl;
}

void SlabMesh::Export(std::string fname){
	fname += "___v_";
	fname += std::to_string(static_cast<long long>(numVertices));
	fname += "___e_";
	fname += std::to_string(static_cast<long long>(numEdges));
	fname += "___f_";
	fname += std::to_string(static_cast<long long>(numFaces));

	AdjustStorage();

	std::string maname = fname;
	maname += ".ma";

	std::ofstream fout(maname);

	//	GraphVertexIterator gvi,gvi_end;

	fout << numVertices << " " << numEdges << " " << numFaces << std::endl;

	//fout << num_vertices(*g) << " " << num_edges(*g) << " " << g->tris.size() << std::endl;

	for(unsigned i = 0; i < vertices.size(); i ++)
		//fout << "v " << vertices[i].second->sphere.center << " " << vertices[i].second->sphere.radius << std::endl;
		fout << "v " << setiosflags(ios::fixed) << setprecision(15) << (vertices[i].second->sphere.center * pmesh->bb_diagonal_length) << " " << (vertices[i].second->sphere.radius * pmesh->bb_diagonal_length) << std::endl;

	for(unsigned i = 0; i < edges.size(); i ++)
		fout << "e " << edges[i].second->vertices_.first << " " << edges[i].second->vertices_.second << std::endl;
	for(unsigned i = 0; i < faces.size(); i ++)
	{
		fout << "f";
		for(std::set<unsigned>::iterator si = faces[i].second->vertices_.begin();
			si != faces[i].second->vertices_.end(); si ++)
			fout << " " << *si;
		fout << std::endl;
	}
	fout.close();
}


void SlabMesh::InitialTopologyProperty(unsigned vid) {
	if (numVertices > 100)
		return;

	set<unsigned> fir_faces = vertices[vid].second->faces_;
	set<unsigned> fir_edges = vertices[vid].second->edges_;
	for (set<unsigned>::iterator si = fir_edges.begin(); si != fir_edges.end(); si++)
	{
		int index = edges[*si].second->vertices_.first == vid ? 
			edges[*si].second->vertices_.second : edges[*si].second->vertices_.first;

		set<unsigned> sec_edges = vertices[index].second->edges_;
		sec_edges.erase(*si);
		for (set<unsigned>::iterator si2 = sec_edges.begin(); si2 != sec_edges.end(); si2++)
		{
			int index2 = edges[*si2].second->vertices_.first == index ? 
				edges[*si2].second->vertices_.second : edges[*si2].second->vertices_.first;

			set<unsigned> third_edges = vertices[index2].second->edges_;
			third_edges.erase(*si2);
			for (set<unsigned>::iterator si3 = third_edges.begin(); si3 != third_edges.end(); si3++)
			{
				int index3 = edges[*si3].second->vertices_.first == index2 ? 
					edges[*si3].second->vertices_.second : edges[*si3].second->vertices_.first;

				if (index3 == vid)
				{
					// �������γ��˻�·��������滹��hole
					bool is_hole = true;
					for (set<unsigned>::iterator fi = fir_faces.begin(); fi != fir_faces.end(); fi++) 
					{
						set<unsigned> ver = faces[*fi].second->vertices_;
						if (ver.find(index) != ver.end() && ver.find(index2) != ver.end())
						{
							is_hole = false;
							break;
						}
					}

					if (is_hole)
					{
						if (edges[*si].second->faces_.size() <= 1 && edges[*si2].second->faces_.size() <= 1 
							&& edges[*si3].second->faces_.size() <= 1)
						{
							edges[*si].second->topo_contractable = false;
							edges[*si2].second->topo_contractable = false;
							edges[*si3].second->topo_contractable = false; 
						}
					}
				}
			}
		}
	}

	// ���򻯵�С��50������ʱ�������������˵�ı߽��кϲ�
	if (numVertices <= 100)
	{
		for (set<unsigned>::iterator si = fir_edges.begin(); si != fir_edges.end(); si++)
		{
			int index = edges[*si].second->vertices_.first == vid ? 
				edges[*si].second->vertices_.second : edges[*si].second->vertices_.first;

			if (vertices[index].second->edges_.size() == 1 && vertices[index].second->faces_.size() == 0)
			{
				edges[*si].second->topo_contractable = false;
			}
		}
	}
}

void SlabMesh::InitialTopologyProperty() {
	for (int i = 0 ;i < vertices.size(); i++)
	{
		if(vertices[i].first)
		{
			InitialTopologyProperty(i);
		}
	}
}

// Assign a topological label to every active vertex using the flags set by
// DistinguishVertexType().  Must be called after LoadSlabMesh().
//
//   NM_CORNER  : fake_boundary_vertex && non_manifold_vertex
//   NM_EDGE    : !fake_boundary_vertex && non_manifold_vertex
//   BOUNDARY   : fake_boundary_vertex && !non_manifold_vertex
//   REGULAR    : neither
//
// Collapse is only allowed between vertices that share the same label.

// ─────────────────────────────────────────────────────────────────────────────
// SlabMesh::DetermineTopology
// ─────────────────────────────────────────────────────────────────────────────
//
// Recomputes topological classification flags on every active SlabVertex by
// inspecting the face-valence of each incident edge:
//
//   edge face-count == 1  →  boundary edge
//   edge face-count == 2  →  manifold (sheet) edge
//   edge face-count  > 2  →  non-manifold (seam) edge
//
// Per-vertex flags set (a vertex can appear in multiple sets):
//   topo_is_sheet    – appears in at least one 2-manifold (sheet) edge.
//   topo_is_boundary – appears in at least one boundary edge (nf == 1).
//   topo_is_seam     – appears in at least one non-manifold edge (nf > 2).
//   topo_is_junction – appears in BOTH a seam AND a boundary edge
//                      (intersection of those two vertex sets).
//
// topo_type priority: junction > seam > boundary > sheet.
//
// is_spike is loaded from the sidecar (set by ComputeInputNMM)
// and just counted here — not recomputed.
//
// Call this after simplification to refresh flags that MergeVertices left
// as conservative approximations.

void SlabMesh::DetermineTopology()
{
	// ── Pass 1: iterate edges, populate vertex sets by edge type ─────────────
	// An edge's type is determined by the number of active faces incident to it:
	//   nf == 1  →  boundary edge
	//   nf == 2  →  sheet edge
	//   nf  > 2  →  seam edge
	// Junction: vertex with >= 3 seam edges incident.

	// Reset nf on all active vertices before the edge loop writes into them.
	for (unsigned i = 0; i < vertices.size(); ++i)
		if (vertices[i].first) vertices[i].second->nf = 0;

	std::set<unsigned> sheet_verts, boundary_verts, seam_verts;
	std::unordered_map<unsigned, unsigned> seam_edge_count; // vid → # of seam edges

	for (unsigned eid = 0; eid < edges.size(); ++eid)
	{
		if (!edges[eid].first) continue;

		const unsigned va = edges[eid].second->vertices_.first;
		const unsigned vb = edges[eid].second->vertices_.second;

		if (va >= vertices.size() || !vertices[va].first) continue;
		if (vb >= vertices.size() || !vertices[vb].first) continue;

		// Count active incident faces.
		unsigned nf = 0;
		for (unsigned fid : edges[eid].second->faces_)
			if (fid < faces.size() && faces[fid].first) ++nf;

		// nf == 0: dangling edge with no faces — skip (spike remnant or isolated).
		if (nf == 0) continue;

		// Record max nf on each endpoint for GUI debugging.
		// if (nf > vertices[va].second->nf) vertices[va].second->nf = nf;
		// if (nf > vertices[vb].second->nf) vertices[vb].second->nf = nf;

		if (nf == 1) { boundary_verts.insert(va); boundary_verts.insert(vb); }
		if (nf == 2) { sheet_verts.insert(va);    sheet_verts.insert(vb);    }
		if (nf > 2)  {
			seam_verts.insert(va); seam_verts.insert(vb);
			++seam_edge_count[va]; ++seam_edge_count[vb];
		}
	}

	// Junction: vertex with >= 3 seam edges.
	std::set<unsigned> junction_verts;
	for (const auto& kv : seam_edge_count)
		if (kv.second >= 3) junction_verts.insert(kv.first);

	// ── Pass 2: assign per-vertex flags and topo_type ────────────────────────
	// Base type from strongest edge; _Boundary suffix if any boundary edge present.

	unsigned n_sheet = 0, n_sheet_b = 0, n_seam = 0, n_seam_b = 0,
	         n_junction = 0, n_junction_b = 0, n_steep = 0, unknown = 0;

	using TT = SlabVertex::TopoType;

	for (unsigned i = 0; i < vertices.size(); ++i)
	{
		if (!vertices[i].first) continue;
		SlabVertex* v = vertices[i].second;

		v->topo_is_sheet    = sheet_verts.count(i)    > 0;
		v->topo_is_boundary = boundary_verts.count(i) > 0;
		v->topo_is_seam     = seam_verts.count(i)     > 0;
		v->topo_is_junction = junction_verts.count(i) > 0;

		if      (v->topo_is_junction && v->topo_is_boundary) { v->topo_type = TT::Junction_Boundary; ++n_junction_b; v->nf = 3; }
		else if (v->topo_is_junction)                        { v->topo_type = TT::Junction;          ++n_junction;   v->nf = 3; }
		else if (v->topo_is_seam     && v->topo_is_boundary) { v->topo_type = TT::Seam_Boundary;     ++n_seam_b;     v->nf = 4; }
		else if (v->topo_is_seam)                            { v->topo_type = TT::Seam;              ++n_seam;       v->nf = 4; }
		else if (v->topo_is_boundary)                        { v->topo_type = TT::Sheet_Boundary;    ++n_sheet_b;    v->nf = 1; }
		else if (v->topo_is_sheet)                           { v->topo_type = TT::Sheet;             ++n_sheet;      v->nf = 2; }
		else                                                 { v->topo_type = TT::Unknown;            ++unknown;      v->nf = -1; }

		if (v->is_spike) ++n_steep;
	}

	std::cout << "[SlabMesh::DetermineTopology]"
	          << "  sheet="            << n_sheet
	          << "  sheet_boundary="   << n_sheet_b
	          << "  seam="             << n_seam
	          << "  seam_boundary="    << n_seam_b
	          << "  junction="         << n_junction
	          << "  junction_boundary="<< n_junction_b
	          << "  steep="            << n_steep
	          << "  unknown="          << unknown << "\n";
}

// ─────────────────────────────────────────────────────────────────────────────
// SlabMesh::RecomputeVertexTopology
// ─────────────────────────────────────────────────────────────────────────────
//
// Recomputes topology flags and topo_type for a single vertex by examining
// only its incident edges.  Mirrors the per-vertex logic of DetermineTopology()
// but operates on one vertex — called on vid_tgt after each MergeVertices().

void SlabMesh::RecomputeVertexTopology(unsigned vid)
{
	if (vid >= vertices.size() || !vertices[vid].first) return;
	SlabVertex* v = vertices[vid].second;

	v->nf               = 0;
	v->topo_is_sheet    = false;
	v->topo_is_boundary = false;
	v->topo_is_seam     = false;
	v->topo_is_junction = false;

	unsigned n_seam_edges = 0;

	for (unsigned eid : v->edges_)
	{
		if (eid >= edges.size() || !edges[eid].first) continue;

		// Count active incident faces on this edge.
		unsigned nf = 0;
		for (unsigned fid : edges[eid].second->faces_)
			if (fid < faces.size() && faces[fid].first) ++nf;

		if (nf == 0) continue;

		if (nf > v->nf) v->nf = nf;

		if (nf == 1) v->topo_is_boundary = true;
		if (nf == 2) v->topo_is_sheet    = true;
		if (nf > 2)  { v->topo_is_seam = true; ++n_seam_edges; }
	}

	if (n_seam_edges >= 3) v->topo_is_junction = true;

	using TT = SlabVertex::TopoType;
	if      (v->topo_is_junction && v->topo_is_boundary) v->topo_type = TT::Junction_Boundary;
	else if (v->topo_is_junction)                        v->topo_type = TT::Junction;
	else if (v->topo_is_seam     && v->topo_is_boundary) v->topo_type = TT::Seam_Boundary;
	else if (v->topo_is_seam)                            v->topo_type = TT::Seam;
	else if (v->topo_is_boundary)                        v->topo_type = TT::Sheet_Boundary;
	else if (v->topo_is_sheet)                           v->topo_type = TT::Sheet;
	else                                                 v->topo_type = TT::Unknown;
}

// ─────────────────────────────────────────────────────────────────────────────
// SlabMesh::CanMerge
// ─────────────────────────────────────────────────────────────────────────────
//
// Three conditions must all hold (Main / Boundary context):
//
// Condition 1 — same TopoType:
//   v1->topo_type == v2->topo_type
//
// Condition 2 — same ClusterType (T-type):
//   v1->nmn_cluster_type == v2->nmn_cluster_type
//
// Condition 3 — bplists are surface-mesh-edge neighbours:
//   At least one bp in v1->nmn_bplist shares a surface mesh edge with at least
//   one bp in v2->nmn_bplist (checked via CGAL halfedge circulator).
//   If no such pair exists, the collapse is rejected.

// Computes connected components of `bps` using input mesh edge connectivity.
// Returns one set per component (union-find over mesh-adjacent bp pairs).
static std::vector<std::set<unsigned>> ComputeBpClusters(
	const std::set<unsigned>& bps,
	const std::vector<Mesh::Vertex_iterator>& vlist,
	unsigned n_mv)
{
	if (bps.empty()) return {};

	std::unordered_map<unsigned, unsigned> parent;
	for (unsigned bp : bps) parent[bp] = bp;

	std::function<unsigned(unsigned)> find = [&](unsigned x) -> unsigned {
		if (parent[x] != x) parent[x] = find(parent[x]);
		return parent[x];
	};
	auto unite = [&](unsigned a, unsigned b) {
		a = find(a); b = find(b);
		if (a != b) parent[a] = b;
	};

	for (unsigned bp : bps)
	{
		if (bp >= n_mv) continue;
		auto circ = vlist[bp]->vertex_begin();
		auto done = circ;
		do {
			unsigned nbr = static_cast<unsigned>(circ->opposite()->vertex()->id);
			if (bps.count(nbr))
				unite(bp, nbr);
		} while (++circ != done);
	}

	std::unordered_map<unsigned, std::set<unsigned>> comp_map;
	for (unsigned bp : bps)
		comp_map[find(bp)].insert(bp);

	std::vector<std::set<unsigned>> result;
	result.reserve(comp_map.size());
	for (auto& [root, members] : comp_map)
		result.push_back(std::move(members));
	return result;
}

static SlabVertex::ClusterType ClusterTypeFromCount(unsigned n)
{
	switch (n) {
		case 0:  return SlabVertex::ClusterType::T0;
		case 1:  return SlabVertex::ClusterType::T1_spike;  // refined by bplist size below
		case 2:  return SlabVertex::ClusterType::T2;
		case 3:  return SlabVertex::ClusterType::T3;
		case 4:  return SlabVertex::ClusterType::T4;
		default: return SlabVertex::ClusterType::T5;
	}
}

// Refines T1 into T1_spike (exactly 4 bpoints) or T1_non_spike (>4 bpoints).
// A true spike sphere is tangent to exactly 4 surface points (Delaunay tetrahedron).
static SlabVertex::ClusterType ClusterTypeFromCountAndBplist(
	unsigned n_clusters, unsigned n_bplist)
{
	SlabVertex::ClusterType ct = ClusterTypeFromCount(n_clusters);
	//uncomment me 
	// if (ct == SlabVertex::ClusterType::T1_spike && n_bplist != 4)
	// 	ct = SlabVertex::ClusterType::T1_non_spike;
	return ct;
}

void SlabMesh::ClusterNMNBplist()
{
	if (!pmesh) return;
	const auto& vlist = pmesh->pVertexList;
	const unsigned n_mv = static_cast<unsigned>(vlist.size());

	// Build the set of all feature vertex IDs (endpoints of sharp/concave edges
	// and corner vertices). These are excluded from bplist before clustering so
	// that sharp-feature boundaries act as hard separators between patches.
	std::set<unsigned> feature_verts;
	for (const auto& e : sharp_edges) {
		feature_verts.insert((unsigned)e[0]);
		feature_verts.insert((unsigned)e[1]);
	}

	for (int v : feature_corners)
		feature_verts.insert((unsigned)v);

	for (unsigned i = 0; i < vertices.size(); ++i)
	{
		if (!vertices[i].first) continue;
		SlabVertex* sv = vertices[i].second;

		// Filter out feature vertices from the bplist before clustering.
		// Feature vertices sit on geometric discontinuities and would
		// incorrectly bridge separate surface patches into one cluster.
		std::set<unsigned> filtered_bps;
		for (unsigned bp : sv->nmn_bplist)
			if (!feature_verts.count(bp))
				filtered_bps.insert(bp);

		// If filtering removes the entire bplist (all bps are feature vertices),
		// fall back to the unfiltered bplist so the vertex is not incorrectly
		// classified as T0 (0 clusters) when it should be T1.
		const std::set<unsigned>* bps_to_cluster = &filtered_bps;
		std::set<unsigned> full_bps;
		if (filtered_bps.empty() && !sv->nmn_bplist.empty()) {
			full_bps.insert(sv->nmn_bplist.begin(), sv->nmn_bplist.end());
			bps_to_cluster = &full_bps;
		}

		sv->nmn_bplist_clusters = ComputeBpClusters(*bps_to_cluster, vlist, n_mv);
		sv->nmn_cluster_type    = ClusterTypeFromCountAndBplist(
			(unsigned)sv->nmn_bplist_clusters.size(),
			(unsigned)sv->nmn_bplist.size());
	}
}

bool SlabMesh::WouldCreateNonManifold(unsigned vid0, unsigned vid1,
                                      RejectionReason*  out_reason,
                                      ReasonPrimitives* out_prims) const
{
	// Translated from PMP is_collapse_ok().  Returns true if collapsing the
	// MAT edge (vid0,vid1) would produce a non-manifold result.

	// ── helpers ──────────────────────────────────────────────────────────────

	// Count active faces on an edge.
	auto activeNF = [&](unsigned eid) -> unsigned {
		if (eid >= edges.size() || !edges[eid].first) return 0;
		unsigned n = 0;
		for (unsigned fid : edges[eid].second->faces_)
			if (fid < faces.size() && faces[fid].first) ++n;
		return n;
	};

	auto isBoundaryEdge = [&](unsigned eid) -> bool {
		return activeNF(eid) == 1;
	};

	// Find the edge ID connecting two vertices (UINT_MAX if none).
	auto findEdge = [&](unsigned va, unsigned vb) -> unsigned {
		if (va >= vertices.size() || !vertices[va].first) return UINT_MAX;
		for (unsigned eid : vertices[va].second->edges_) {
			if (eid >= edges.size() || !edges[eid].first) continue;
			if (edges[eid].second->vertices_.first  == vb ||
				edges[eid].second->vertices_.second == vb)
				return eid;
		}
		return UINT_MAX;
	};

	auto isBoundaryVertex = [&](unsigned v) -> bool {
		if (v >= vertices.size() || !vertices[v].first) return false;
		for (unsigned eid : vertices[v].second->edges_)
			if (isBoundaryEdge(eid)) return true;
		return false;
	};

	// ── find shared edge and its incident faces ───────────────────────────────

	unsigned shared_eid = findEdge(vid0, vid1);
	if (shared_eid == UINT_MAX) return false; // no edge between them

	// Collect ALL third vertices of incident faces into a set.
	// Using a set (not just vl/vr) correctly handles non-manifold seam edges
	// that have 3+ incident faces — all their third vertices are "allowed"
	// shared neighbours in the link condition test.
	std::set<unsigned> face_third_verts;
	unsigned edge_nf = 0;

	for (unsigned fid : edges[shared_eid].second->faces_) {
		if (fid >= faces.size() || !faces[fid].first) continue;
		++edge_nf;
		for (unsigned v : faces[fid].second->vertices_)
			if (v != vid0 && v != vid1)
				face_third_verts.insert(v);
	}

	bool edge_is_boundary = (edge_nf == 1);

	// ── PMP test 1 & 2: per-face boundary-edge pair check ────────────────────
	// For each incident face, the other two edges of that triangle must not
	// both be boundary edges.
	for (unsigned ft : face_third_verts) {
		unsigned e1t = findEdge(vid1, ft);
		unsigned et0 = findEdge(ft,   vid0);
		if (e1t != UINT_MAX && et0 != UINT_MAX &&
			isBoundaryEdge(e1t) && isBoundaryEdge(et0))
		{
			if (out_reason) *out_reason = RejectionReason::NonManifold_BoundaryEdgePair;
			if (out_prims) {
				out_prims->vertices = { ft };
				out_prims->edges    = { {vid1, ft}, {ft, vid0} };
				out_prims->faces    = { {vid0, vid1, ft} };
			}
			return true;
		}
	}

	// ── PMP test 3: two incident faces share the same third vertex ────────────
	// In the manifold case this is vl == vr.  In the non-manifold case it
	// can't happen (each face has a distinct third vertex), so we only trigger
	// if there are exactly 2 incident faces and their third vertex is the same.
	if (edge_nf == 2 && face_third_verts.size() == 1)
	{
		if (out_reason) *out_reason = RejectionReason::NonManifold_SharedThirdVert;
		if (out_prims) {
			unsigned ft = *face_third_verts.begin();
			out_prims->vertices = { ft };
			// Both incident faces share the same three vertices.
			out_prims->faces    = { {vid0, vid1, ft}, {vid0, vid1, ft} };
		}
		return true;
	}

	// ── PMP test 4: boundary-vertex / boundary-edge consistency ──────────────
	if (isBoundaryVertex(vid0) && isBoundaryVertex(vid1) && !edge_is_boundary)
	{
		if (out_reason) *out_reason = RejectionReason::NonManifold_BoundaryVertEdge;
		if (out_prims) {
			out_prims->vertices = { vid0, vid1 };
			out_prims->edges    = { {vid0, vid1} };
		}
		return true;
	}

	// ── PMP test 5: link condition (one-ring intersection) ───────────────────
	// The one-rings of vid0 and vid1 must only intersect at the set of exempt
	// third vertices.  Any other shared neighbour would cause a duplicate face
	// or an over-valent edge after collapse.
	//
	// Special case: edge_nf == 0 means the collapsed edge has no incident faces
	// (e.g. MS_Seam edges from the MatStruct .ma file that live only on the
	// 1-skeleton with no MAT triangles).  There are no faces to violate, so the
	// link condition is trivially satisfied — skip it entirely.
	if (edge_nf == 0) return false;

	// Exemption rule: ALL face_third_verts are exempt.
	// Any vertex V that appears in a face {vid0, vid1, V} is connected to both
	// endpoints by definition.  During MergeVertices that face will be deleted,
	// collapsing the two edges vid0-V and vid1-V into one.  The shared connection
	// is therefore safely resolved and does NOT produce a non-manifold result.
	// This holds regardless of how many incident faces V appears in (manifold or
	// non-manifold seam edge).  Only a shared neighbour that is NOT a face third
	// vertex would create an over-valent edge after collapse.
	const std::set<unsigned>& exempt_verts = face_third_verts;

	std::set<unsigned> nbrs1;
	for (unsigned eid : vertices[vid1].second->edges_) {
		if (eid >= edges.size() || !edges[eid].first) continue;
		unsigned a = edges[eid].second->vertices_.first;
		unsigned b = edges[eid].second->vertices_.second;
		nbrs1.insert(a == vid1 ? b : a);
	}

	for (unsigned eid : vertices[vid0].second->edges_) {
		if (eid >= edges.size() || !edges[eid].first) continue;
		unsigned a = edges[eid].second->vertices_.first;
		unsigned b = edges[eid].second->vertices_.second;
		unsigned nbr = (a == vid0) ? b : a;
		if (nbr == vid1 || exempt_verts.count(nbr)) continue;
		if (nbrs1.count(nbr))
		{
			if (out_reason) *out_reason = RejectionReason::NonManifold_LinkCondition;
			if (out_prims) {
				// nbr is a shared neighbour of both vid0 and vid1 that is not an
				// exempt face-third vertex — it would create a duplicate edge.
				out_prims->vertices = { nbr };
				out_prims->edges    = { {vid0, nbr}, {vid1, nbr} };
			}
			return true;
		}
	}

	return false; // all tests passed — collapse is topologically safe
}

// ── 3D segment crossing test ──────────────────────────────────────────────────
// Returns true if segment (A,B) and segment (C,D) cross in 3D.
//
// NOTE: The signed-volume straddling approach was attempted but is mathematically
// broken in 3D — signedVol(A,B,D,C) = -signedVol(A,B,C,D) always, making the
// product always ≤ 0 and causing every non-coplanar pair to be flagged. It is
// left below as a comment for reference.
//
// Current approach: minimum distance between the two segments.
// A crossing is detected when:
//   1. The closest points on each segment are strictly interior (s,t ∈ (eps,1-eps))
//   2. The distance between them is below a threshold relative to the shorter edge.
// This correctly handles skew 3D segments — skew segments have a non-zero minimum
// distance and will not be flagged as crossings.
static bool SegmentsCross3D(const Wm4::Vector3d& A, const Wm4::Vector3d& B,
                             const Wm4::Vector3d& C, const Wm4::Vector3d& D)
{
	// --- Commented-out signed-volume approach (broken in 3D) ---
	// auto signedVol = [](const Wm4::Vector3d& P, const Wm4::Vector3d& Q,
	//                     const Wm4::Vector3d& R, const Wm4::Vector3d& S) -> double {
	//     return (Q-P).Dot((R-P).Cross(S-P));
	// };
	// double d1 = signedVol(A,B,C,D); double d2 = signedVol(A,B,D,C); // d2 == -d1 always
	// double d3 = signedVol(C,D,A,B); double d4 = signedVol(C,D,B,A);
	// return (d1*d2 < 0.0 && d3*d4 < 0.0); // always true when non-coplanar → wrong

	// Minimum distance between segments (A,B) and (C,D).
	// Parametric form: P(s) = A + s*(B-A), Q(t) = C + t*(D-C), s,t ∈ [0,1].
	const Wm4::Vector3d d1 = B - A;
	const Wm4::Vector3d d2 = D - C;
	const Wm4::Vector3d r  = A - C;

	const double a = d1.Dot(d1); // squared length of AB
	const double e = d2.Dot(d2); // squared length of CD
	const double f = d2.Dot(r);

	static const double kDegen = 1e-14;
	double s, t;

	if (a <= kDegen && e <= kDegen) {
		// Both segments degenerate to points.
		s = t = 0.0;
	} else if (a <= kDegen) {
		s = 0.0;
		t = f / e;
		t = std::max(0.0, std::min(1.0, t));
	} else {
		const double c = d1.Dot(r);
		if (e <= kDegen) {
			t = 0.0;
			s = std::max(0.0, std::min(1.0, -c / a));
		} else {
			const double b     = d1.Dot(d2);
			const double denom = a * e - b * b;
			if (std::abs(denom) > kDegen) {
				s = std::max(0.0, std::min(1.0, (b * f - c * e) / denom));
			} else {
				s = 0.0; // parallel — pick arbitrary s
			}
			t = (b * s + f) / e;
			// Clamp t then recompute s.
			if (t < 0.0) {
				t = 0.0;
				s = std::max(0.0, std::min(1.0, -c / a));
			} else if (t > 1.0) {
				t = 1.0;
				s = std::max(0.0, std::min(1.0, (b - c) / a));
			}
		}
	}

	// Closest point distance.
	const Wm4::Vector3d closest = (A + d1 * s) - (C + d2 * t);
	const double dist = closest.Length();

	// Threshold: fraction of the shorter segment length.
	// 1e-3 means the segments must come within 0.1% of the shorter edge length
	// to be considered crossing — tight enough to miss genuine skew edges but
	// catch coplanar crossings.
	const double shorter = std::sqrt(std::min(a, e));
	const double threshold = shorter * 1e-3;

	// Only flag as crossing when both closest points are strictly interior
	// (not at endpoints) — endpoint sharing is legal in a connected graph.
	const double kEndPt = 1e-4;
	const bool s_interior = (s > kEndPt && s < 1.0 - kEndPt);
	const bool t_interior = (t > kEndPt && t < 1.0 - kEndPt);

	return (dist < threshold && s_interior && t_interior);
}

// Returns true if collapsing edge (vid0,vid1) to v_tgt would cause any of the
// resulting new edges to geometrically cross an existing edge — the fold-over /
// polyline self-intersection problem not caught by WouldCreateNonManifold or
// Contractible (both blind to 1-D loop geometry).
bool SlabMesh::WouldCreateFoldOver(unsigned vid0, unsigned vid1,
                                    const Wm4::Vector3d& v_tgt,
                                    ReasonPrimitives* out_prims) const
{
	if (!vertices[vid0].first || !vertices[vid1].first) return false;

	// ── Step 1: collect neighbors of the merged vertex ────────────────────────
	// After the collapse, the new vertex connects to every neighbor of vid0 and
	// vid1 except each other.  Shared neighbors (via a common face) appear once.
	std::set<unsigned> new_nbrs;
	for (unsigned eid : vertices[vid0].second->edges_) {
		if (!edges[eid].first) continue;
		unsigned nbr = (edges[eid].second->vertices_.first == vid0)
		               ? edges[eid].second->vertices_.second
		               : edges[eid].second->vertices_.first;
		if (nbr != vid1) new_nbrs.insert(nbr);
	}
	for (unsigned eid : vertices[vid1].second->edges_) {
		if (!edges[eid].first) continue;
		unsigned nbr = (edges[eid].second->vertices_.first == vid1)
		               ? edges[eid].second->vertices_.second
		               : edges[eid].second->vertices_.first;
		if (nbr != vid0) new_nbrs.insert(nbr);
	}

	// ── Step 2: one-ring exclusion set ───────────────────────────────────────
	// Any edge (Y,Z) that shares an endpoint with the new vertex cannot properly
	// cross one of its edges (they share a point).  Skip all such edges.
	std::set<unsigned> one_ring(new_nbrs);
	one_ring.insert(vid0);
	one_ring.insert(vid1);

	// ── Step 3: for each new edge (v_tgt → X), test against all other edges ──
	// We check all active edges in the mesh (not just the 2-hop neighbourhood)
	// so that long-range crossings — typical in boundary-loop fold-overs — are
	// also detected.
	for (unsigned X : new_nbrs) {
		if (!vertices[X].first) continue;
		const Wm4::Vector3d& pos_X = vertices[X].second->sphere.center;

		for (unsigned eid_c = 0; eid_c < (unsigned)edges.size(); ++eid_c) {
			if (!edges[eid_c].first) continue;
			unsigned Y = edges[eid_c].second->vertices_.first;
			unsigned Z = edges[eid_c].second->vertices_.second;
			// Skip edges that touch the one-ring — shared endpoint → cannot cross.
			if (one_ring.count(Y) || one_ring.count(Z)) continue;

			const Wm4::Vector3d& pos_Y = vertices[Y].second->sphere.center;
			const Wm4::Vector3d& pos_Z = vertices[Z].second->sphere.center;

			if (SegmentsCross3D(v_tgt, pos_X, pos_Y, pos_Z)) {
				if (out_prims) {
					// X is the new neighbor whose edge to v_tgt would cross [Y,Z].
					out_prims->vertices = { X, Y, Z };
					out_prims->edges    = { {Y, Z} };
				}
				return true;
			}
		}
	}

	return false; // no crossing found — collapse is geometrically safe
}

// ─────────────────────────────────────────────────────────────────────────────
// WouldExceedCurvatureThreshold
// Runtime check: collapse edge (vid0,vid1) only if BOTH endpoints sit on a
// locally straight portion of their boundary/seam chain.
//
// For each endpoint V, find its one same-type chain neighbour on the far side
// (i.e. not the other collapsing vertex).  Then compute the turning angle at V
// between that far neighbour, V itself, and the other endpoint:
//
//   chain: ... A — vid0 — vid1 — D ...
//   turning at vid0 = acos(-normalize(A-vid0) · normalize(vid1-vid0))
//   turning at vid1 = acos(-normalize(vid0-vid1) · normalize(D-vid1))
//
// Return true (reject) if either angle > feature_angle_threshold.
// Also reject if an endpoint has ≥ 2 other same-type neighbours (junction-like).
// Skip the check if an endpoint has 0 other same-type neighbours (chain end).
// ─────────────────────────────────────────────────────────────────────────────
bool SlabMesh::WouldExceedCurvatureThreshold(unsigned vid0, unsigned vid1,
                                              ReasonPrimitives* out_prims) const
{
	using CT = SlabVertex::ClusterType;
	if (!vertices[vid0].first || !vertices[vid1].first) return false;

	const CT ct_0 = vertices[vid0].second->nmn_cluster_type;
	const CT ct_1 = vertices[vid1].second->nmn_cluster_type;
	const double threshold_rad = feature_angle_threshold * M_PI / 180.0;

	// Find the single same-type neighbour of 'vid' (matching 'ct') excluding 'exclude'.
	// Returns UINT_MAX = chain end, UINT_MAX-1 = junction-like (multiple same-type neighbours).
	auto farNeighbour = [&](unsigned vid, unsigned exclude, CT ct) -> unsigned {
		unsigned found = UINT_MAX;
		for (unsigned eid : vertices[vid].second->edges_) {
			if (eid >= edges.size() || !edges[eid].first) continue;
			unsigned nbr = (edges[eid].second->vertices_.first == vid)
			               ? edges[eid].second->vertices_.second
			               : edges[eid].second->vertices_.first;
			if (nbr == exclude) continue;
			if (!vertices[nbr].first) continue;
			if (vertices[nbr].second->nmn_cluster_type != ct) continue;
			if (found != UINT_MAX) return UINT_MAX - 1;
			found = nbr;
		}
		return found;
	};

	// Returns true if endpoint 'vid' is the cause; also populates out_prims.
	auto checkEndpoint = [&](unsigned vid, unsigned partner, CT ct) -> bool {
		unsigned far = farNeighbour(vid, partner, ct);
		if (far == UINT_MAX)     return false; // chain end — no angle to check
		if (far == UINT_MAX - 1) {             // junction-like — always reject
			if (out_prims) {
				out_prims->vertices = { vid };
			}
			return true;
		}

		const Wm4::Vector3d& pV   = vertices[vid].second->sphere.center;
		const Wm4::Vector3d& pFar = vertices[far].second->sphere.center;
		const Wm4::Vector3d& pPrt = vertices[partner].second->sphere.center;
		Wm4::Vector3d u = pFar - pV;
		Wm4::Vector3d w = pPrt - pV;
		double lu = u.Length(), lw = w.Length();
		if (lu < 1e-14 || lw < 1e-14) return false;
		u /= lu; w /= lw;
		double dot = -(u.Dot(w));
		dot = std::max(-1.0, std::min(1.0, dot));
		if (std::acos(dot) > threshold_rad) {
			if (out_prims) {
				// vid is the bent vertex; far is its chain neighbour on the far side.
				out_prims->vertices = { vid, far };
				out_prims->edges    = { {vid, far}, {vid, partner} };
			}
			return true;
		}
		return false;
	};

	return checkEndpoint(vid0, vid1, ct_0) || checkEndpoint(vid1, vid0, ct_1);
}

// ─────────────────────────────────────────────────────────────────────────────
// MarkSharpFeatureVertices
// Pre-simplification pass over MS_Boundary, MS_Seam, and MS_Junction vertices.
//
//   MS_Junction: always marked sharp — by definition they are branch points
//     where multiple feature chains meet; collapsing them would destroy topology.
//
//   MS_Boundary / MS_Seam: marked sharp if the vertex has >= 2 same-type
//     neighbours AND the minimum turning angle across all neighbour pairs
//     exceeds angle_deg_threshold.
//     Turning angle = deviation from straight when traversing A → V → B:
//       = acos( -normalize(A-V) · normalize(B-V) )  (0° = straight, π = U-turn)
// ─────────────────────────────────────────────────────────────────────────────
void SlabMesh::MarkSharpFeatureVertices(double angle_deg_threshold)
{
	using CT = SlabVertex::ClusterType;
	const double threshold_rad = angle_deg_threshold * M_PI / 180.0;

	unsigned marked = 0;

	for (unsigned vid = 0; vid < (unsigned)vertices.size(); ++vid) {
		if (!vertices[vid].first) continue;
		SlabVertex* sv = vertices[vid].second;

		const CT ct = sv->nmn_cluster_type;

		// Junctions (and junction+boundary) are always sharp — branch points by definition.
		if (ct == CT::MS_Junction || ct == CT::MS_Junction_Boundary) {
			sv->sharpNotContractable = true;
			++marked;
			continue;
		}

		// Only seam / boundary / their compound variants participate in the angle test.
		if (ct != CT::MS_Boundary && ct != CT::MS_Seam &&
		    ct != CT::MS_Seam_Boundary && ct != CT::MS_Sheet_Boundary) continue;

		// Collect same-type neighbours.
		std::vector<unsigned> same_nbrs;
		for (unsigned eid : sv->edges_) {
			if (eid >= edges.size() || !edges[eid].first) continue;
			unsigned nbr = (edges[eid].second->vertices_.first == vid)
			               ? edges[eid].second->vertices_.second
			               : edges[eid].second->vertices_.first;
			if (!vertices[nbr].first) continue;
			if (vertices[nbr].second->nmn_cluster_type == ct)
				same_nbrs.push_back(nbr);
		}

		// Need >= 2 same-type neighbours to form a chain angle.
		if (same_nbrs.size() < 2) continue;

		// For each pair of same-type neighbours, compute the turning angle at V.
		// We take the minimum turning angle across all pairs — if any pair is
		// nearly straight the vertex is collapsible; only mark sharp when ALL
		// pairs exceed the threshold (i.e. every direction is a sharp turn).
		const Wm4::Vector3d& pV = sv->sphere.center;
		double min_turning = M_PI; // start high

		for (size_t i = 0; i < same_nbrs.size(); ++i) {
			for (size_t j = i + 1; j < same_nbrs.size(); ++j) {
				const Wm4::Vector3d& pA = vertices[same_nbrs[i]].second->sphere.center;
				const Wm4::Vector3d& pB = vertices[same_nbrs[j]].second->sphere.center;
				Wm4::Vector3d u = pA - pV;
				Wm4::Vector3d w = pB - pV;
				double lu = u.Length(), lw = w.Length();
				if (lu < 1e-14 || lw < 1e-14) continue;
				u /= lu;  w /= lw;
				// turning angle: deviation from straight (0 = collinear, π = U-turn)
				double dot = -(u.Dot(w));  // = cos(turning_angle)
				dot = std::max(-1.0, std::min(1.0, dot));
				double turning = std::acos(dot);
				if (turning < min_turning ||turning >= 90.0f  ) min_turning = turning;
			}
		}

		if (min_turning > threshold_rad) {
			sv->sharpNotContractable = true;
			++marked;
		}
	}

	std::cerr << "[MarkSharpFeatureVertices] threshold=" << angle_deg_threshold
	          << " deg  marked=" << marked << " vertices as sharpNotContractable\n";
}

bool SlabMesh::CanMerge(unsigned vid1, unsigned vid2,
                        RejectionReason*  out_reason,
                        ReasonPrimitives* out_prims) const
{
	const SlabVertex* v1 = vertices[vid1].second;
	const SlabVertex* v2 = vertices[vid2].second;

	// Condition 0: sharp feature protection.
	{
		using CT = SlabVertex::ClusterType;
		if (v1->nmn_cluster_type == CT::MS_Junction          ||
		    v1->nmn_cluster_type == CT::MS_Junction_Boundary  ||
		    v2->nmn_cluster_type == CT::MS_Junction          ||
		    v2->nmn_cluster_type == CT::MS_Junction_Boundary)
		{
			if (out_reason) *out_reason = RejectionReason::SharpNotContractable;
			if (out_prims) {
				using CT2 = SlabVertex::ClusterType;
				if (v1->nmn_cluster_type == CT2::MS_Junction ||
				    v1->nmn_cluster_type == CT2::MS_Junction_Boundary)
					out_prims->vertices.push_back(vid1);
				if (v2->nmn_cluster_type == CT2::MS_Junction ||
				    v2->nmn_cluster_type == CT2::MS_Junction_Boundary)
					out_prims->vertices.push_back(vid2);
			}
			return false;
		}
	}
	// Same boundary-related type pairs → honour sharpNotContractable.
	{
		using CT = SlabVertex::ClusterType;
		const CT c1 = v1->nmn_cluster_type, c2 = v2->nmn_cluster_type;
		bool same_boundary_pair = (c1 == c2) &&
		    (c1 == CT::MS_Boundary       ||
		     c1 == CT::MS_Sheet_Boundary ||
		     c1 == CT::MS_Seam_Boundary  ||
		     c1 == CT::MS_Seam);
		if (same_boundary_pair && (v1->sharpNotContractable || v2->sharpNotContractable))
		{
			if (out_reason) *out_reason = RejectionReason::SharpNotContractable;
			if (out_prims) {
				if (v1->sharpNotContractable) out_prims->vertices.push_back(vid1);
				if (v2->sharpNotContractable) out_prims->vertices.push_back(vid2);
			}
			return false;
		}
	}

	// Condition 2: same cluster type (T-type).
	if (v1->nmn_cluster_type != v2->nmn_cluster_type)
	{
		if (out_reason) *out_reason = RejectionReason::DifferentClusterType;
		if (out_prims)  out_prims->vertices = { vid1, vid2 };
		return false;
	}

	// Condition 4: link condition — thread out_prims so the sub-reason geometry is captured.
	if (WouldCreateNonManifold(vid1, vid2, out_reason, out_prims))
		return false;

	return true;
}

bool SlabMesh::CanMerge(unsigned vid1, unsigned vid2, CollapseContext ctx) const
{
	switch (ctx)
	{
	case CollapseContext::Spike:
		// Spike edges are force-collapsed — bypass all boundary and topology
		// checks.  topo_contractable is still enforced inside MinCostEdgeCollapse.
		return true;

	case CollapseContext::Boundary:
	case CollapseContext::Main:
	default:
		return CanMerge(vid1, vid2);
	}
}

// ─────────────────────────────────────────────────────────────────────────────
// SlabMesh::LogCollapseRejection
// ─────────────────────────────────────────────────────────────────────────────
void SlabMesh::LogCollapseRejection(const char* queue_name,
                                    unsigned eid, unsigned v1, unsigned v2,
                                    double cost, RejectionReason reason,
                                    ReasonPrimitives prims) const
{
	// Record the latest rejection reason and primitives per edge.
	if (eid < (unsigned)edges.size()) {
		edge_last_rejection[eid]    = reason;
		edge_reason_primitives[eid] = std::move(prims);
	}
	static const char* ct_names[] = {
		"T0","T1_spike","T2","T3","T4","T5","T1_non_spike",
		"MS_Unknown","MS_Sheet","MS_Seam","MS_Boundary","MS_Junction",
		"MS_Sheet_Boundary","MS_Seam_Boundary","MS_Junction_Boundary"
	};
	static const char* tt_names[] = {
		"Unknown","Sheet","Sheet_Boundary","Seam","Seam_Boundary","Junction","Junction_Boundary"
	};
	static const char* reason_names[] = {
		"StaleEdge","InvalidVertex",
		"DifferentTopoType","DifferentClusterType",
		"BplistNotNeighbors","NoPmesh",
		"TopoNotContractable","InversionWouldOccur",
		"NonManifold_BoundaryEdgePair",
		"NonManifold_SharedThirdVert",
		"NonManifold_BoundaryVertEdge",
		"NonManifold_LinkCondition",
		"WouldCreateFoldOver",
		"SharpNotContractable",
		"WouldExceedCurvatureThreshold",
	};

	auto ct_name = [&](unsigned vid) -> const char* {
		if (vid >= vertices.size() || !vertices[vid].first) return "deleted";
		uint8_t idx = static_cast<uint8_t>(vertices[vid].second->nmn_cluster_type);
		return idx < 15 ? ct_names[idx] : "???";
	};
	auto tt_name = [&](unsigned vid) -> const char* {
		if (vid >= vertices.size() || !vertices[vid].first) return "deleted";
		uint8_t idx = static_cast<uint8_t>(vertices[vid].second->topo_type);
		return idx < 7 ? tt_names[idx] : "???";
	};

	const std::string phase_tag = current_phase.empty() ? queue_name : current_phase;
	std::ofstream log(export_prefix + "_rejection_log_" + phase_tag + ".txt", std::ios::app);
	if (!log) return;

	uint8_t r = static_cast<uint8_t>(reason);
	log << "[" << queue_name << "] REJECTED"
	    << "  edge=" << eid
	    << "  cost=" << cost
	    << "  reason=" << (r < std::size(reason_names) ? reason_names[r] : "???") << "\n"
	    << "    v1=" << v1 << "  cluster=" << ct_name(v1) << "  topo=" << tt_name(v1) << "\n"
	    << "    v2=" << v2 << "  cluster=" << ct_name(v2) << "  topo=" << tt_name(v2) << "\n";
}

// ─────────────────────────────────────────────────────────────────────────────
// SlabMesh::ExportSkeletonPLY
// Exports remaining active edges as PLY cylinders with per-vertex RGB colour
// based on the last collapse-rejection reason recorded for each edge.
// ─────────────────────────────────────────────────────────────────────────────
void SlabMesh::ExportTypedMA(const std::string& path) const
{
	static const char* ct_names[] = {
		"T0","T1_spike","T2","T3","T4","T5","T1_non_spike",
		"MS_Unknown","MS_Sheet","MS_Seam","MS_Boundary","MS_Junction",
		"MS_Sheet_Boundary","MS_Seam_Boundary","MS_Junction_Boundary"
	};

	// Build compact old→new index maps without modifying the mesh.
	std::vector<unsigned> newv(vertices.size(), UINT_MAX);
	std::vector<unsigned> newe(edges.size(),    UINT_MAX);
	std::vector<unsigned> newf(faces.size(),    UINT_MAX);
	unsigned cv = 0, ce = 0, cf = 0;
	for (unsigned i = 0; i < (unsigned)vertices.size(); ++i) if (vertices[i].first) newv[i] = cv++;
	for (unsigned i = 0; i < (unsigned)edges.size();    ++i) if (edges[i].first)    newe[i] = ce++;
	for (unsigned i = 0; i < (unsigned)faces.size();    ++i) if (faces[i].first)    newf[i] = cf++;

	std::ofstream f(path);
	if (!f) { std::cerr << "[ExportTypedMA] cannot open: " << path << "\n"; return; }

	const double scale = pmesh ? pmesh->bb_diagonal_length : 1.0;
	f << std::fixed << std::setprecision(15);
	f << cv << " " << ce << " " << cf << "\n";

	for (unsigned i = 0; i < (unsigned)vertices.size(); ++i) {
		if (!vertices[i].first) continue;
		const auto& c = vertices[i].second->sphere.center;
		const double  r = vertices[i].second->sphere.radius;
		const uint8_t t = static_cast<uint8_t>(vertices[i].second->nmn_cluster_type);
		const char* name = (t < 15) ? ct_names[t] : "MS_Unknown";
		f << "v " << c.X()*scale << " " << c.Y()*scale << " " << c.Z()*scale
		  << " " << r*scale << " " << name << "\n";
	}
	for (unsigned i = 0; i < (unsigned)edges.size(); ++i) {
		if (!edges[i].first) continue;
		f << "e " << newv[edges[i].second->vertices_.first]
		  << " "  << newv[edges[i].second->vertices_.second] << "\n";
	}
	for (unsigned i = 0; i < (unsigned)faces.size(); ++i) {
		if (!faces[i].first) continue;
		f << "f";
		for (unsigned vid : faces[i].second->vertices_)
			f << " " << newv[vid];
		f << "\n";
	}
	std::cerr << "[ExportTypedMA] wrote " << cv << " verts, " << ce << " edges, "
	          << cf << " faces to " << path << "\n";
}

void SlabMesh::ExportClusterPLY(const std::string& path) const
{
	// Cluster-type colour table — mirrors kClusterTypeColors in main_cli.cpp.
	// Index = uint8_t value of SlabVertex::ClusterType.
	struct Col { uint8_t r, g, b; };
	static const Col kCTCol[15] = {
		{ 229,   0, 229 }, // 0  T0                   magenta
		{  25,  25,  25 }, // 1  T1_spike             near-black
		{   0, 216, 255 }, // 2  T2                   bright cyan
		{ 255, 127,   0 }, // 3  T3                   vivid orange
		{ 255,  25,  25 }, // 4  T4                   vivid red
		{ 255, 255, 255 }, // 5  T5                   white
		{ 140, 140, 140 }, // 6  T1_non_spike         grey
		{  89,  89,  89 }, // 7  MS_Unknown           dark grey
		{   0, 255,  76 }, // 8  MS_Sheet             vivid green
		{ 255, 229,   0 }, // 9  MS_Seam              vivid yellow
		{   0, 127, 255 }, // 10 MS_Boundary          vivid blue
		{ 255,   0, 127 }, // 11 MS_Junction          vivid pink
		{   0, 229, 255 }, // 12 MS_Sheet_Boundary    bright cyan
		{ 255,  89,   0 }, // 13 MS_Seam_Boundary     vivid orange-red
		{ 153,   0, 255 }, // 14 MS_Junction_Boundary vivid purple
	};
	auto vertCol = [&](unsigned vid) -> Col {
		if (!vertices[vid].first) return { 89, 89, 89 };
		uint8_t t = static_cast<uint8_t>(vertices[vid].second->nmn_cluster_type);
		return (t < 15) ? kCTCol[t] : kCTCol[7];
	};

	// Compact index maps (no AdjustStorage side-effects)
	const double scale = pmesh ? pmesh->bb_diagonal_length : 1.0;
	std::vector<unsigned> newv(vertices.size(), UINT_MAX);
	unsigned cv = 0, ce = 0, cf = 0;
	for (unsigned i = 0; i < (unsigned)vertices.size(); ++i) if (vertices[i].first) newv[i] = cv++;
	for (unsigned i = 0; i < (unsigned)edges.size();    ++i) if (edges[i].first)    ++ce;
	for (unsigned i = 0; i < (unsigned)faces.size();    ++i) if (faces[i].first)    ++cf;

	std::ofstream ply(path);
	if (!ply) { std::cerr << "[ExportClusterPLY] cannot open: " << path << "\n"; return; }

	// PLY header — vertices, faces, edges
	ply << "ply\nformat ascii 1.0\n"
	    << "element vertex " << cv << "\n"
	    << "property float x\nproperty float y\nproperty float z\n"
	    << "property uchar red\nproperty uchar green\nproperty uchar blue\n"
	    << "element face " << cf << "\n"
	    << "property list uchar int vertex_indices\n"
	    << "element edge " << ce << "\n"
	    << "property int vertex1\nproperty int vertex2\n"
	    << "end_header\n";

	ply << std::fixed << std::setprecision(10);

	// Vertices
	for (unsigned i = 0; i < (unsigned)vertices.size(); ++i) {
		if (!vertices[i].first) continue;
		const auto& c = vertices[i].second->sphere.center;
		Col col = vertCol(i);
		ply << c.X()*scale << " " << c.Y()*scale << " " << c.Z()*scale
		    << " " << (int)col.r << " " << (int)col.g << " " << (int)col.b << "\n";
	}
	// Faces
	for (unsigned i = 0; i < (unsigned)faces.size(); ++i) {
		if (!faces[i].first) continue;
		ply << faces[i].second->vertices_.size();
		for (unsigned vid : faces[i].second->vertices_)
			ply << " " << newv[vid];
		ply << "\n";
	}
	// Edges
	for (unsigned i = 0; i < (unsigned)edges.size(); ++i) {
		if (!edges[i].first) continue;
		ply << newv[edges[i].second->vertices_.first]
		    << " " << newv[edges[i].second->vertices_.second] << "\n";
	}

	std::cerr << "[ExportClusterPLY] wrote " << cv << " verts, " << cf << " faces, "
	          << ce << " edges to " << path << "\n";
}

void SlabMesh::ExportSkeletonPLY(const std::string& path, double radius) const
{
	// Colour helper — use the shared table in SlabMesh.h.
	struct Col { uint8_t r, g, b; };
	auto reasonToColor = [](RejectionReason rr) -> Col {
		auto c = RejectionReasonColorU8(rr);
		return { c[0], c[1], c[2] };
	};



//   ┌─────────────────────────────────────────────┬──────────────────────┐
//   │                   Reason                    │        Color         │
//   ├─────────────────────────────────────────────┼──────────────────────┤
//   │ TopoNotContractable                         │ 🔴 RED 255,0,0       │
//   ├─────────────────────────────────────────────┼──────────────────────┤
//   │ NonManifold_BoundaryEdgePair                │ 🟠 ORANGE 255,165,0  │
//   ├─────────────────────────────────────────────┼──────────────────────┤
//   │ NonManifold_SharedThirdVert                 │ 🟡 YELLOW 255,255,0  │
//   ├─────────────────────────────────────────────┼──────────────────────┤
//   │ NonManifold_BoundaryVertEdge                │ 🟢 GREEN 0,255,0     │
//   ├─────────────────────────────────────────────┼──────────────────────┤
//   │ NonManifold_LinkCondition                   │ 🩵 CYAN 0,255,255    │
//   ├─────────────────────────────────────────────┼──────────────────────┤
//   │ WouldCreateFoldOver                         │ 🔵 BLUE 0,0,255      │
//   ├─────────────────────────────────────────────┼──────────────────────┤
//   │ SharpNotContractable                        │ 🟣 VIOLET 148,0,211  │
//   ├─────────────────────────────────────────────┼──────────────────────┤
//   │ WouldExceedCurvatureThreshold               │ 🟣 MAGENTA 255,0,255 │
//   ├─────────────────────────────────────────────┼──────────────────────┤
//   │ Less critical (stale/invalid/type mismatch) │ Greys                │
//   ├─────────────────────────────────────────────┼──────────────────────┤
//   │ Never attempted                             │ ⬜ WHITE             │
//   └─────────────────────────────────────────────┴──────────────────────┘

	// // ── Debug: dump map state and active edge IDs to file ───────────────────
	// {
	// 	std::ofstream dbg(export_prefix + "_skeleton_debug.txt");
	// 	unsigned activeEdges = 0;
	// 	for (unsigned i = 0; i < (unsigned)edges.size(); ++i)
	// 		if (edges[i].first) ++activeEdges;
	// 	dbg << "active_edges=" << activeEdges
	// 	    << "  edge_last_rejection.size()=" << edge_last_rejection.size() << "\n\n";
	// 	// All entries in the map
	// 	dbg << "-- map contents --\n";
	// 	for (auto& kv : edge_last_rejection)
	// 		dbg << "  eid=" << kv.first << "  reason=" << (int)kv.second
	// 		    << "  edge_active=" << (kv.first < edges.size() && edges[kv.first].first ? "yes" : "NO") << "\n";
	// 	// All active edges and whether they appear in the map
	// 	dbg << "\n-- active edges vs map --\n";
	// 	for (unsigned i = 0; i < (unsigned)edges.size(); ++i) {
	// 		if (!edges[i].first) continue;
	// 		dbg << "  eid=" << i
	// 		    << (edge_last_rejection.count(i) ? "  IN_MAP" : "  NOT_IN_MAP") << "\n";
	// 	}
	// }

	// ── Collect geometry (PLY header needs total counts upfront) ──────────────
	struct Vert { float x, y, z; uint8_t r, g, b; };
	struct Tri  { int a, b, c; };

	std::vector<Vert> verts;
	std::vector<Tri>  tris;
	verts.reserve(edges.size() * 16);
	tris.reserve(edges.size() * 16);

	const int    N  = 8;
	const double pi = std::acos(-1.0);

	for (unsigned i = 0; i < (unsigned)edges.size(); ++i)
	{
		if (!edges[i].first) continue;
		const unsigned vid0 = edges[i].second->vertices_.first;
		const unsigned vid1 = edges[i].second->vertices_.second;
		if (!vertices[vid0].first || !vertices[vid1].first) continue;

		const Wm4::Vector3d p0 = vertices[vid0].second->sphere.center;
		const Wm4::Vector3d p1 = vertices[vid1].second->sphere.center;

		Wm4::Vector3d d = p1 - p0;
		const double len = d.Length();
		if (len < 1e-10) continue;
		d /= len;

		// Orthonormal basis perpendicular to the edge direction
		Wm4::Vector3d ref = (std::abs(d.X()) < 0.9)
		                    ? Wm4::Vector3d(1.0, 0.0, 0.0)
		                    : Wm4::Vector3d(0.0, 1.0, 0.0);
		Wm4::Vector3d u = d.Cross(ref); u.Normalize();
		Wm4::Vector3d v = d.Cross(u);

		// Colour from last rejection reason (white if never attempted)
		auto it = edge_last_rejection.find(i);
		const Col col = (it != edge_last_rejection.end())
		                ? reasonToColor(it->second)
		                : Col{ 255, 255, 255 };

		const int base = (int)verts.size();

		// 2*N vertices: alternating p0/p1 rings
		for (int s = 0; s < N; ++s)
		{
			const double angle = 2.0 * pi * s / N;
			const Wm4::Vector3d off = u * (radius * std::cos(angle))
			                        + v * (radius * std::sin(angle));
			const Wm4::Vector3d q0 = p0 + off;
			const Wm4::Vector3d q1 = p1 + off;
			verts.push_back({ (float)q0.X(), (float)q0.Y(), (float)q0.Z(), col.r, col.g, col.b });
			verts.push_back({ (float)q1.X(), (float)q1.Y(), (float)q1.Z(), col.r, col.g, col.b });
		}

		// 2*N triangles forming the tube sides
		for (int s = 0; s < N; ++s)
		{
			const int sn = (s + 1) % N;
			const int a  = base + s  * 2;
			const int b  = base + s  * 2 + 1;
			const int c  = base + sn * 2;
			const int dd = base + sn * 2 + 1;
			tris.push_back({ a, b, dd });
			tris.push_back({ a, dd, c });
		}
	}

	// ── Write PLY ASCII ───────────────────────────────────────────────────────
	std::ofstream ply(path);
	if (!ply.is_open()) {
		std::cerr << "[ExportSkeletonPLY] ERROR: cannot open: " << path << "\n";
		return;
	}

	ply << "ply\n"
	    << "format ascii 1.0\n"
	    << "comment MAT edge skeleton -- vertex colour encodes last collapse rejection reason\n"
	    << "comment white=never_attempted red=DifferentTopoType orange=DifferentClusterType\n"
	    << "comment yellow=BplistNotNeighbors dark_blue=TopoNotContractable purple=InversionWouldOccur\n"
	    << "comment cyan/blue=NonManifold magenta=FoldOver green=SharpNotContractable lime=CurvatureThreshold\n"
	    << "element vertex " << verts.size() << "\n"
	    << "property float x\n"
	    << "property float y\n"
	    << "property float z\n"
	    << "property uchar red\n"
	    << "property uchar green\n"
	    << "property uchar blue\n"
	    << "element face " << tris.size() << "\n"
	    << "property list uchar int vertex_indices\n"
	    << "end_header\n";

	ply << std::fixed << std::setprecision(6);
	for (const auto& vt : verts)
		ply << vt.x << " " << vt.y << " " << vt.z << " "
		    << (int)vt.r << " " << (int)vt.g << " " << (int)vt.b << "\n";

	for (const auto& tr : tris)
		ply << "3 " << tr.a << " " << tr.b << " " << tr.c << "\n";

	std::cerr << "[ExportSkeletonPLY] "
	          << (verts.size() / (N * 2)) << " cylinders ("
	          << verts.size() << " verts, " << tris.size() << " tris) -> "
	          << path << "\n";
}