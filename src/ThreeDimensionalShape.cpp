#include "ThreeDimensionalShape.h"

// Note: QString include removed - was unused and prevents CLI build without Qt

#include <deque>
#include <fstream>
#include "matfp/sharp_feature_detection.h"
#include "matfp/Args.h"

void ThreeDimensionalShape::ComputeInputNMM()
{
	input_nmm.numVertices = 0;
	input_nmm.numEdges = 0;
	input_nmm.numFaces = 0;

	input_nmm.vertices.clear();
	input_nmm.edges.clear();
	input_nmm.faces.clear();

	Triangulation * pt = &(input.dt);

	
	num_vor_v = 0;
	num_vor_e = 0;
	num_vor_f = 0;

	double len[4];
	len[0] = input.m_max[0] - input.m_min[0];
	len[1] = input.m_max[1] - input.m_min[1];
	len[2] = input.m_max[2] - input.m_min[2];
	len[3] = sqrt(len[0]*len[0]+len[1]*len[1]+len[2]*len[2]);
	input_nmm.diameter = len[3];

	// Recover exact original vertex positions via the ID stored during computedt().
	// Using fvi->point() here would give slightly jittered coordinates (the jitter
	// is applied inside computedt() to break degenerate DT configurations).
	for(Finite_vertices_iterator_t fvi = pt->finite_vertices_begin(); fvi != pt->finite_vertices_end(); fvi ++)
	{
		int orig_id = fvi->info().id;
		const auto& orig_pt = input.pVertexList[orig_id]->point();
		input_nmm.BoundaryPoints.push_back(SamplePoint(orig_pt[0], orig_pt[1], orig_pt[2]));
	}

	// Export all sample points to <meshname>_sampledpoints.txt
	{
		std::string sp_fname = input_nmm.meshname + "_sampledpoints.txt";
		std::ofstream sp_ofs(sp_fname.c_str());
		if (sp_ofs.is_open())
		{
			sp_ofs << "# id x y z\n";
			sp_ofs << input_nmm.BoundaryPoints.size() << "\n";
			for (unsigned i = 0; i < input_nmm.BoundaryPoints.size(); i++)
			{
				const SamplePoint& sp = input_nmm.BoundaryPoints[i];
				sp_ofs << i
					<< " " << sp.X()
					<< " " << sp.Y()
					<< " " << sp.Z() << "\n";
			}
			sp_ofs.close();
		}
	}

	int mas_vertex_count(0);
	for(Finite_cells_iterator_t fci = pt->finite_cells_begin(); fci != pt->finite_cells_end(); fci ++)
	{
		if(fci->info().inside == false)
		{
			fci->info().tag = -1;
			continue;
		}
		fci->info().tag = mas_vertex_count ++;
		

		Bool_VertexPointer bvp;
		bvp.first = true;
		bvp.second = new NonManifoldMesh_Vertex;
		(*bvp.second).sphere.center = to_wm4(CGAL::circumcenter(pt->tetrahedron(fci)));
		(*bvp.second).is_pole = fci->info().is_pole;
		(*bvp.second).is_spike = false;
		(*bvp.second).topo_is_sheet    = false;
		(*bvp.second).topo_is_seam     = false;
		(*bvp.second).topo_is_junction = false;
		(*bvp.second).topo_is_boundary = false;
		for(unsigned k = 0; k < 4; k ++)
			(*bvp.second).bplist.insert(fci->vertex(k)->info().id);
		(*bvp.second).sphere.radius = pt->TetCircumRadius(pt->tetrahedron(fci));
		input_nmm.vertices.push_back(bvp);
		input_nmm.numVertices ++;
		num_vor_v ++;


	}
	
	for(Finite_facets_iterator_t ffi = pt->finite_facets_begin(); ffi != pt->finite_facets_end(); ffi ++)
	{
		Triangulation::Object o = pt->dual(*ffi);
		if(const Triangulation::Segment *s = CGAL::object_cast<Triangulation::Segment>(&o))
		{
			if( (ffi->first->info().inside == false) || (pt->mirror_facet(*ffi).first->info().inside == false) )
				continue;
			Bool_EdgePointer bep;
			bep.first = true;
			bep.second = new NonManifoldMesh_Edge;
			(*bep.second).vertices_.first = ffi->first->info().tag;
			(*bep.second).vertices_.second = pt->mirror_facet(*ffi).first->info().tag;
			(*input_nmm.vertices[ffi->first->info().tag].second).edges_.insert(input_nmm.edges.size());
			(*input_nmm.vertices[pt->mirror_facet(*ffi).first->info().tag].second).edges_.insert(input_nmm.edges.size());
			input_nmm.edges.push_back(bep);
			input_nmm.numEdges ++;
			num_vor_e ++;

		}
	}
	
	for(Finite_edges_iterator_t fei = pt->finite_edges_begin(); fei != pt->finite_edges_end(); fei ++)
	{
		bool all_finite_inside = true;
		std::vector<Cell_handle_t> vec_ch;
		Cell_circulator_t cc = pt->incident_cells(*fei);
		do
		{
			if(pt->is_infinite(cc))
				all_finite_inside = false;
			else if(cc->info().inside == false)
				all_finite_inside = false;
			vec_ch.push_back(cc++);
		}while(cc != pt->incident_cells(*fei));
		if(!all_finite_inside)
			continue;

		for(unsigned k = 2; k < vec_ch.size() - 1; k ++)
		{
			Bool_EdgePointer bep;
			bep.first = true;
			bep.second = new NonManifoldMesh_Edge;
			(*bep.second).vertices_.first = vec_ch[0]->info().tag;
			(*bep.second).vertices_.second = vec_ch[k]->info().tag;
			(*input_nmm.vertices[vec_ch[0]->info().tag].second).edges_.insert(input_nmm.edges.size());
			(*input_nmm.vertices[vec_ch[k]->info().tag].second).edges_.insert(input_nmm.edges.size());
			input_nmm.edges.push_back(bep);
			input_nmm.numEdges ++;

		}

		for(unsigned k = 1; k < vec_ch.size() - 1; k ++)
		{
			Bool_FacePointer bfp;
			bfp.first = true;
			bfp.second = new NonManifoldMesh_Face;
			unsigned vid[3];
			vid[0] = vec_ch[0]->info().tag;
			vid[1] = vec_ch[k]->info().tag;
			vid[2] = vec_ch[k+1]->info().tag;
			(*bfp.second).vertices_.insert(vec_ch[0]->info().tag);
			(*bfp.second).vertices_.insert(vec_ch[k]->info().tag);
			(*bfp.second).vertices_.insert(vec_ch[k+1]->info().tag);
			unsigned eid[3];
			if(input_nmm.Edge(vid[0],vid[1],eid[0]))
				(*bfp.second).edges_.insert(eid[0]);
			if(input_nmm.Edge(vid[0],vid[2],eid[1]))
				(*bfp.second).edges_.insert(eid[1]);
			if(input_nmm.Edge(vid[1],vid[2],eid[2]))
				(*bfp.second).edges_.insert(eid[2]);
			input_nmm.vertices[vid[0]].second->faces_.insert(input_nmm.faces.size());
			input_nmm.vertices[vid[1]].second->faces_.insert(input_nmm.faces.size());
			input_nmm.vertices[vid[2]].second->faces_.insert(input_nmm.faces.size());
			input_nmm.edges[eid[0]].second->faces_.insert(input_nmm.faces.size());
			input_nmm.edges[eid[1]].second->faces_.insert(input_nmm.faces.size());
			input_nmm.edges[eid[2]].second->faces_.insert(input_nmm.faces.size());
			input_nmm.faces.push_back(bfp);
			input_nmm.numFaces ++;
			num_vor_f ++;
			
		}
	}
	input_nmm.Export(input_nmm.meshname);
	// ── topology (must run before Export so sidecar can be written) ──────────
	// DetermineTopology();

	// ── tag steep tetrahedra: all 4 bplist points mutually mesh-edge connected ─
	// ── tag steep tetrahedra ─────────────────────────────────────────────────
	// A MAT vertex is steep when its 4 responsible boundary points form a quad
	// on the input mesh: two adjacent triangles sharing an edge, where all four
	// vertices of the quad are in the bplist.
	// For each edge (bp1,bp2) between bplist points, check if the two faces
	// sharing that edge each have their third vertex also in the bplist.
	// {
	// 	const unsigned n_mv = static_cast<unsigned>(input.pVertexList.size());
	// 	for (unsigned i = 0; i < input_nmm.vertices.size(); ++i)
	// 	{
	// 		if (!input_nmm.vertices[i].first) continue;
	// 		NonManifoldMesh_Vertex* mv = input_nmm.vertices[i].second;
	// 		mv->is_spike = false;
	// 		if (mv->bplist.size() != 4) continue;
	// 		const std::set<unsigned>& bps = mv->bplist;
	// 		bool is_quad = false;
	// 		for (unsigned bp1 : bps) {
	// 			if (is_quad) break;
	// 			if (bp1 >= n_mv) continue;
	// 			auto circ = input.pVertexList[bp1]->vertex_begin();
	// 			auto done = circ;
	// 			do {
	// 				unsigned bp2 = static_cast<unsigned>(circ->opposite()->vertex()->id);
	// 				if (bps.count(bp2) &&
	// 				    !circ->is_border() && !circ->opposite()->is_border())
	// 				{
	// 					unsigned va = static_cast<unsigned>(circ->next()->vertex()->id);
	// 					unsigned vb = static_cast<unsigned>(circ->opposite()->next()->vertex()->id);
	// 					if (bps.count(va) && bps.count(vb) && va != vb)
	// 					{
	// 						is_quad = true;
	// 						break;
	// 					}
	// 				}
	// 			} while (++circ != done);
	// 		}
	// 		mv->is_spike = is_quad;
	// 	}
	// }


	// ── write sidecar topology file ───────────────────────────────────────────
	// Format:
	//   # comments
	//   <N>
	//   <index> <steep> <sheet> <seam> <junction> <topo_boundary>
	// //           <bp_count> <bp0> <bp1> ...
	// {
	// 	std::string topo_fname = input_nmm.meshname + "_mat_topo.txt";
	// 	std::ofstream tof(topo_fname.c_str());
	// 	if (tof.is_open())
	// 	{
	// 		tof << "# QMAT MAT vertex topology + boundary-point sidecar\n";
	// 		tof << "# index  steep  sheet  seam  junction  topo_boundary"
	// 		       "  bp_count  bp_id x bp_count\n";

	// 		// count active vertices (all should be active here)
	// 		unsigned active_v = 0;
	// 		for (unsigned i = 0; i < input_nmm.vertices.size(); ++i)
	// 			if (input_nmm.vertices[i].first) ++active_v;
	// 		tof << active_v << "\n";

	// 		for (unsigned i = 0; i < input_nmm.vertices.size(); ++i)
	// 		{
	// 			if (!input_nmm.vertices[i].first) continue;
	// 			NonManifoldMesh_Vertex* v = input_nmm.vertices[i].second;

	// 			tof << i
	// 			    << " " << (int)v->is_spike
	// 			    << " " << (int)v->topo_is_sheet
	// 			    << " " << (int)v->topo_is_seam
	// 			    << " " << (int)v->topo_is_junction
	// 			    << " " << (int)v->topo_is_boundary
	// 			    << " " << v->bplist.size();
	// 			for (unsigned bp : v->bplist)
	// 				tof << " " << bp;
	// 			tof << "\n";
	// 		}
	// 		tof.close();
	// 		std::cout << "[ComputeInputNMM]  written sidecar -> " << topo_fname << "\n";
	// 	}
	// 	else
	// 	{
	// 		std::cerr << "[ComputeInputNMM]  WARNING: could not write sidecar "
	// 		          << topo_fname << "\n";
	// 	}
	// }

	// Save vertex-to-sample-point associations:
	// For each Voronoi vertex (circumcenter), record its sphere and the IDs + coordinates
	// of the 4 Delaunay vertices (original sample points) that define its tetrahedron.
	{
		std::string sample_fname = input_nmm.meshname + "_surf_2_vor_mat_map.txt";
		std::ofstream ofs(sample_fname.c_str());
		if (ofs.is_open())
		{
			ofs << "# vertex_index cx cy cz radius num_sample_points\n";
			ofs << "# s sample_point_id px py pz\n";
			ofs << input_nmm.numVertices << "\n";
			for (unsigned i = 0; i < input_nmm.vertices.size(); i++)
			{
				if (!input_nmm.vertices[i].first)
					continue;
				NonManifoldMesh_Vertex* v = input_nmm.vertices[i].second;
				ofs << "v " << i
					<< " " << v->sphere.center.X()
					<< " " << v->sphere.center.Y()
					<< " " << v->sphere.center.Z()
					<< " " << v->sphere.radius
					<< " " << v->bplist.size() << "\n";
				for (std::set<unsigned>::iterator it = v->bplist.begin(); it != v->bplist.end(); ++it)
				{
					unsigned pid = *it;
					const SamplePoint& sp = input_nmm.BoundaryPoints[pid];
					ofs << "s " << pid
						<< " " << sp.X()
						<< " " << sp.Y()
						<< " " << sp.Z() << "\n";
				}
			}
			ofs.close();
		}
	}

	

	input_nmm.numVertices = 0;
	input_nmm.numEdges = 0;
	input_nmm.numFaces = 0;
	//input_nmm.ComputeFacesNormal();
	//input_nmm.ComputeFacesCentroid();
	//input_nmm.ComputeFacesSimpleTriangles();
	//input_nmm.ComputeEdgesCone();

	//slab_mesh.ComputeFacesCentroid();
	//slab_mesh.ComputeFacesNormal();
	//slab_mesh.ComputeVerticesNormal();

}

void ThreeDimensionalShape::LoadInputNMM(std::string fname){
	std::ifstream mastream(fname.c_str());
	NonManifoldMesh newinputnmm;
	newinputnmm.numVertices = 0;
	newinputnmm.numEdges = 0;
	newinputnmm.numFaces = 0;
	int nv, ne, nf;
	mastream >> nv >> ne >> nf;

	// slab mesh
	slab_mesh.numVertices = 0;
	slab_mesh.numEdges = 0;
	slab_mesh.numFaces = 0;

	double len[4];
	len[0] = input.m_max[0] - input.m_min[0];
	len[1] = input.m_max[1] - input.m_min[1];
	len[2] = input.m_max[2] - input.m_min[2];
	len[3] = sqrt(len[0]*len[0]+len[1]*len[1]+len[2]*len[2]);
	newinputnmm.diameter = len[3];
	slab_mesh.bound_weight = 0.1; 

	for(unsigned i = 0; i < input.pVertexList.size(); i ++)
		newinputnmm.BoundaryPoints.push_back(SamplePoint(
		input.pVertexList[i]->point()[0],
		input.pVertexList[i]->point()[1],
		input.pVertexList[i]->point()[2]
	));

	for(unsigned i = 0; i < nv; i ++)
	{
		char ch;
		double x,y,z,r;
		mastream >> ch >> x >> y >> z >> r;

		//Bool_VertexPointer bvp;
		//bvp.first = true;
		//bvp.second = new NonManifoldMesh_Vertex;
		//(*bvp.second).sphere.center[0] = x / input.bb_diagonal_length;
		//(*bvp.second).sphere.center[1] = y / input.bb_diagonal_length;
		//(*bvp.second).sphere.center[2] = z / input.bb_diagonal_length;
		//(*bvp.second).sphere.radius = r / input.bb_diagonal_length;
		//newinputnmm.vertices.push_back(bvp);
		//newinputnmm.numVertices ++;

		// handle the slab mesh
		Bool_SlabVertexPointer bsvp2;
		bsvp2.first = true;
		bsvp2.second = new SlabVertex;
		(*bsvp2.second).sphere.center[0] = x / input.bb_diagonal_length;
		(*bsvp2.second).sphere.center[1] = y / input.bb_diagonal_length;
		(*bsvp2.second).sphere.center[2] = z / input.bb_diagonal_length;
		(*bsvp2.second).sphere.radius = r / input.bb_diagonal_length;
		(*bsvp2.second).index = slab_mesh.vertices.size();
		slab_mesh.vertices.push_back(bsvp2);
		slab_mesh.numVertices ++;
	}

	for(unsigned i = 0; i < ne; i ++)
	{
		char ch;
		unsigned ver[2];
		mastream >> ch;
		mastream >> ver[0];
		mastream >> ver[1];

		//Bool_EdgePointer bep;
		//bep.first = true;
		//bep.second = new NonManifoldMesh_Edge;
		//(*bep.second).vertices_.first = ver[0];
		//(*bep.second).vertices_.second = ver[1];
		//(*newinputnmm.vertices[(*bep.second).vertices_.first].second).edges_.insert(newinputnmm.edges.size());
		//(*newinputnmm.vertices[(*bep.second).vertices_.second].second).edges_.insert(newinputnmm.edges.size());
		//newinputnmm.edges.push_back(bep);
		//newinputnmm.numEdges ++;

		// handle the slab mesh
		Bool_SlabEdgePointer bsep2;
		bsep2.first = true;
		bsep2.second = new SlabEdge;
		(*bsep2.second).vertices_.first = ver[0];
		(*bsep2.second).vertices_.second = ver[1];
		(*slab_mesh.vertices[(*bsep2.second).vertices_.first].second).edges_.insert(slab_mesh.edges.size());
		(*slab_mesh.vertices[(*bsep2.second).vertices_.second].second).edges_.insert(slab_mesh.edges.size());
		(*bsep2.second).index = slab_mesh.edges.size();
		slab_mesh.edges.push_back(bsep2);
		slab_mesh.numEdges ++;
	}

	for(unsigned i = 0; i < nf; i ++)
	{
		char ch;
		unsigned vid[3];
		unsigned eid[3];
		mastream >> ch >> vid[0] >> vid[1] >> vid[2];

		//Bool_FacePointer bfp;
		//bfp.first = true;
		//bfp.second = new NonManifoldMesh_Face;
		//(*bfp.second).vertices_.insert(vid[0]);
		//(*bfp.second).vertices_.insert(vid[1]);
		//(*bfp.second).vertices_.insert(vid[2]);
		//if(newinputnmm.Edge(vid[0],vid[1],eid[0]))
		//	(*bfp.second).edges_.insert(eid[0]);
		//if(newinputnmm.Edge(vid[0],vid[2],eid[1]))
		//	(*bfp.second).edges_.insert(eid[1]);
		//if(newinputnmm.Edge(vid[1],vid[2],eid[2]))
		//	(*bfp.second).edges_.insert(eid[2]);
		//newinputnmm.vertices[vid[0]].second->faces_.insert(newinputnmm.faces.size());
		//newinputnmm.vertices[vid[1]].second->faces_.insert(newinputnmm.faces.size());
		//newinputnmm.vertices[vid[2]].second->faces_.insert(newinputnmm.faces.size());
		//newinputnmm.edges[eid[0]].second->faces_.insert(newinputnmm.faces.size());
		//newinputnmm.edges[eid[1]].second->faces_.insert(newinputnmm.faces.size());
		//newinputnmm.edges[eid[2]].second->faces_.insert(newinputnmm.faces.size());
		//newinputnmm.faces.push_back(bfp);
		//newinputnmm.numFaces ++;

		// handle the slab mesh	
		Bool_SlabFacePointer bsfp2;
		bsfp2.first = true;
		bsfp2.second = new SlabFace;
		(*bsfp2.second).vertices_.insert(vid[0]);
		(*bsfp2.second).vertices_.insert(vid[1]);
		(*bsfp2.second).vertices_.insert(vid[2]);
		if(slab_mesh.Edge(vid[0],vid[1],eid[0]))
			(*bsfp2.second).edges_.insert(eid[0]);
		if(slab_mesh.Edge(vid[0],vid[2],eid[1]))
			(*bsfp2.second).edges_.insert(eid[1]);
		if(slab_mesh.Edge(vid[1],vid[2],eid[2]))
			(*bsfp2.second).edges_.insert(eid[2]);
		(*bsfp2.second).index = slab_mesh.faces.size();
		slab_mesh.vertices[vid[0]].second->faces_.insert(slab_mesh.faces.size());
		//slab_mesh.vertices[vid[0]].second->related_face += 2;
		slab_mesh.vertices[vid[1]].second->faces_.insert(slab_mesh.faces.size());
		//slab_mesh.vertices[vid[1]].second->related_face += 2;
		slab_mesh.vertices[vid[2]].second->faces_.insert(slab_mesh.faces.size());
		//slab_mesh.vertices[vid[2]].second->related_face += 2;
		slab_mesh.edges[eid[0]].second->faces_.insert(slab_mesh.faces.size());
		slab_mesh.edges[eid[1]].second->faces_.insert(slab_mesh.faces.size());
		slab_mesh.edges[eid[2]].second->faces_.insert(slab_mesh.faces.size());
		slab_mesh.faces.push_back(bsfp2);
		slab_mesh.numFaces++;
	}

	//newinputnmm.ComputeFacesNormal();
	//newinputnmm.ComputeFacesCentroid();
	//newinputnmm.ComputeFacesSimpleTriangles();
	//newinputnmm.ComputeEdgesCone();
	//input_nmm = newinputnmm;

	slab_mesh.iniNumVertices = slab_mesh.numVertices;
	slab_mesh.iniNumEdges = slab_mesh.numEdges;
	slab_mesh.iniNumFaces = slab_mesh.numFaces;

	slab_mesh.CleanIsolatedVertices();
	slab_mesh.computebb();
	slab_mesh.ComputeFacesCentroid();
	slab_mesh.ComputeFacesNormal();
	slab_mesh.ComputeVerticesNormal();
	slab_mesh.ComputeEdgesCone();
	slab_mesh.ComputeFacesSimpleTriangles();
	slab_mesh.DistinguishVertexType();

	// ── load sidecar topology file ────────────────────────────────────────────
	// Derive sidecar path from .ma file path (strip ".ma", append "_mat_topo.txt")
	{
		std::string sidecar = fname;
		const std::string ma_ext = ".ma";
		if (sidecar.size() >= ma_ext.size() &&
		    sidecar.compare(sidecar.size() - ma_ext.size(), ma_ext.size(), ma_ext) == 0)
			sidecar = sidecar.substr(0, sidecar.size() - ma_ext.size());
		sidecar += "_mat_topo.txt";

		std::ifstream sf(sidecar.c_str());
		if (sf.is_open())
		{
			// skip comment lines, read record count
			std::string line;
			unsigned n_rec = 0;
			while (std::getline(sf, line))
			{
				if (line.empty() || line[0] == '#') continue;
				n_rec = static_cast<unsigned>(std::stoul(line));
				break;
			}

			for (unsigned r = 0; r < n_rec; ++r)
			{
				unsigned idx;
				int steep, sheet, seam, junc, tbound;
				unsigned bp_count;
				sf >> idx >> steep >> sheet >> seam >> junc >> tbound >> bp_count;

				std::set<unsigned> bplist;
				for (unsigned b = 0; b < bp_count; ++b)
				{
					unsigned bp;
					sf >> bp;
					bplist.insert(bp);
				}

				if (idx < slab_mesh.vertices.size() && slab_mesh.vertices[idx].first)
				{
					SlabVertex* sv = slab_mesh.vertices[idx].second;
					sv->nmn_bplist = bplist;
					// is_spike and topo_is_* are recomputed by
					// ClusterNMNBplist() + SlabMesh::DetermineTopology().
				}
			}
			sf.close();
			std::cout << "[LoadInputNMM]  loaded sidecar -> " << sidecar << "\n";
		}
		else
		{
			std::cerr << "[LoadInputNMM]  WARNING: sidecar not found: "
			          << sidecar << "\n";
		}
	}

	// ── load Voronoi neighbor sidecar ─────────────────────────────────────────
	{
		std::string vn_sidecar = fname;
		const std::string ma_ext = ".ma";
		if (vn_sidecar.size() >= ma_ext.size() &&
		    vn_sidecar.compare(vn_sidecar.size() - ma_ext.size(),
		                       ma_ext.size(), ma_ext) == 0)
			vn_sidecar = vn_sidecar.substr(0, vn_sidecar.size() - ma_ext.size());
		vn_sidecar += "_voronoi_neighbors.txt";

		std::ifstream vf(vn_sidecar.c_str());
		if (vf.is_open())
		{
			const unsigned n_bp =
			    static_cast<unsigned>(input.pVertexList.size());
			slab_mesh.voronoi_neighbors.assign(n_bp, std::set<unsigned>());

			// skip comment lines, read record count
			std::string line;
			unsigned n_rec = 0;
			while (std::getline(vf, line))
			{
				if (line.empty() || line[0] == '#') continue;
				n_rec = static_cast<unsigned>(std::stoul(line));
				break;
			}

			for (unsigned r = 0; r < n_rec; ++r)
			{
				unsigned bp_id, nb_count;
				vf >> bp_id >> nb_count;
				for (unsigned n = 0; n < nb_count; ++n)
				{
					unsigned nb;
					vf >> nb;
					if (bp_id < n_bp && nb < n_bp)
						slab_mesh.voronoi_neighbors[bp_id].insert(nb);
				}
			}
			vf.close();
			std::cout << "[LoadInputNMM]  loaded Voronoi neighbors -> "
			          << vn_sidecar << "\n";
		}
		else
		{
			std::cerr << "[LoadInputNMM]  WARNING: Voronoi neighbor sidecar not found: "
			          << vn_sidecar << "\n";
		}
	}
}

#ifdef USE_MATSTRUCT_INITIALIZATION
void ThreeDimensionalShape::LoadMatstructMA(std::string fname)
{
	// ── Open geometry file (.ma / .mat) ──────────────────────────────────────
	std::ifstream mastream(fname.c_str());
	if (!mastream.is_open()) {
		std::cerr << "[LoadMatstructMA] ERROR: cannot open file: " << fname << "\n";
		return;
	}

	int nv, ne, nf;
	if (!(mastream >> nv >> ne >> nf) || nv <= 0) {
		std::cerr << "[LoadMatstructMA] ERROR: bad header in: " << fname << "\n";
		std::exit(1);
	}

	slab_mesh.numVertices = 0;
	slab_mesh.numEdges    = 0;
	slab_mesh.numFaces    = 0;
	slab_mesh.bound_weight = 0.1;

	// ── Vertices: format is "v x y z r" (no inline type field) ───────────────
	for (int i = 0; i < nv; ++i)
	{
		char ch;
		double x, y, z, r;
		mastream >> ch >> x >> y >> z >> r;

		Bool_SlabVertexPointer bsvp;
		bsvp.first  = true;
		bsvp.second = new SlabVertex;
		bsvp.second->sphere.center[0] = x / input.bb_diagonal_length;
		bsvp.second->sphere.center[1] = y / input.bb_diagonal_length;
		bsvp.second->sphere.center[2] = z / input.bb_diagonal_length;
		bsvp.second->sphere.radius     = r / input.bb_diagonal_length;
		bsvp.second->index             = slab_mesh.vertices.size();
		bsvp.second->nmn_cluster_type  = SlabVertex::ClusterType::MS_Unknown;
		slab_mesh.vertices.push_back(bsvp);
		slab_mesh.numVertices++;
	}

	// ── Edges: format is "e v0 v1" ────────────────────────────────────────────
	for (int i = 0; i < ne; ++i)
	{
		char ch;
		unsigned v0, v1;
		mastream >> ch >> v0 >> v1;

		Bool_SlabEdgePointer bsep;
		bsep.first  = true;
		bsep.second = new SlabEdge;
		bsep.second->vertices_.first  = v0;
		bsep.second->vertices_.second = v1;
		slab_mesh.vertices[v0].second->edges_.insert(slab_mesh.edges.size());
		slab_mesh.vertices[v1].second->edges_.insert(slab_mesh.edges.size());
		bsep.second->index = slab_mesh.edges.size();
		slab_mesh.edges.push_back(bsep);
		slab_mesh.numEdges++;
	}

	// ── Faces: format is "f v0 v1 v2" ────────────────────────────────────────
	for (int i = 0; i < nf; ++i)
	{
		char ch;
		unsigned vid[3];
		unsigned eid[3];
		mastream >> ch >> vid[0] >> vid[1] >> vid[2];

		Bool_SlabFacePointer bsfp;
		bsfp.first  = true;
		bsfp.second = new SlabFace;
		bsfp.second->vertices_.insert(vid[0]);
		bsfp.second->vertices_.insert(vid[1]);
		bsfp.second->vertices_.insert(vid[2]);
		if (slab_mesh.Edge(vid[0], vid[1], eid[0])) bsfp.second->edges_.insert(eid[0]);
		if (slab_mesh.Edge(vid[0], vid[2], eid[1])) bsfp.second->edges_.insert(eid[1]);
		if (slab_mesh.Edge(vid[1], vid[2], eid[2])) bsfp.second->edges_.insert(eid[2]);
		bsfp.second->index = slab_mesh.faces.size();
		slab_mesh.vertices[vid[0]].second->faces_.insert(slab_mesh.faces.size());
		slab_mesh.vertices[vid[1]].second->faces_.insert(slab_mesh.faces.size());
		slab_mesh.vertices[vid[2]].second->faces_.insert(slab_mesh.faces.size());
		slab_mesh.edges[eid[0]].second->faces_.insert(slab_mesh.faces.size());
		slab_mesh.edges[eid[1]].second->faces_.insert(slab_mesh.faces.size());
		slab_mesh.edges[eid[2]].second->faces_.insert(slab_mesh.faces.size());
		slab_mesh.faces.push_back(bsfp);
		slab_mesh.numFaces++;
	}

	slab_mesh.iniNumVertices = slab_mesh.numVertices;
	slab_mesh.iniNumEdges    = slab_mesh.numEdges;
	slab_mesh.iniNumFaces    = slab_mesh.numFaces;

	std::cout << "[LoadMatstructMA] loaded " << nv << " verts, "
	          << ne << " edges, " << nf << " faces from " << fname << "\n";

	// ── Load struct classification from the same stream ─────────────────────
	// The struct blocks are appended to the same file immediately after the
	// face list (see MAT_struct_format.md).  Keep reading from mastream.
	{
		std::istream& sf = mastream;
		using CT = SlabVertex::ClusterType;
		const unsigned nv_total = (unsigned)slab_mesh.vertices.size();

		// Track the "main" type (sheet/seam/junction) and boundary flag separately
		// so we can produce compound types at the end.
		// Main type priority: Junction(3) > Seam(2) > Sheet(1) > Unknown(0).
		// Boundary is independent — it combines with the main type.
		std::vector<int>  main_prio(nv_total, 0); // 0=unknown,1=sheet,2=seam,3=junction
		std::vector<bool> has_boundary(nv_total, false);

		auto applyMain = [&](unsigned vid, int prio) {
			if (vid < nv_total && slab_mesh.vertices[vid].first)
				if (prio > main_prio[vid]) main_prio[vid] = prio;
		};
		auto applyBoundary = [&](unsigned vid) {
			if (vid < nv_total && slab_mesh.vertices[vid].first)
				has_boundary[vid] = true;
		};

		int num_structs = 0;
		if (!(sf >> num_structs)) {
			std::cerr << "[LoadMatstructMA] WARNING: no struct data found in " << fname
			          << " — all vertices will be typed MS_Unknown.\n";
			num_structs = 0;
		}
		int sheet_structs = 0, seam_structs = 0, boundary_structs = 0, junction_structs = 0;

		// Buffer all struct blocks so both passes can iterate them.
		struct StructBlock { int struct_id, type_id; std::vector<unsigned> elems; };
		std::vector<StructBlock> struct_blocks;
		struct_blocks.reserve(num_structs);
		for (int s = 0; s < num_structs; ++s)
		{
			StructBlock blk;
			int count;
			sf >> blk.struct_id >> blk.type_id >> count;
			blk.elems.resize(count);
			for (int j = 0; j < count; ++j) sf >> blk.elems[j];
			switch (blk.type_id) {
				case 0: ++sheet_structs;    break;
				case 1: ++seam_structs;     break;
				case 2: ++boundary_structs; break;
				case 3: ++junction_structs; break;
			}
			struct_blocks.push_back(std::move(blk));
		}

		// ── PASS 1: struct_id stamping + vertex type computation ─────────────
		// Process in priority order (SEAM > BOUNDARY > SHEET for edges,
		// JUNCTION > all for vertices) and assign nmn_cluster_type.
		std::unordered_set<unsigned> seam_stamped;     // edge ids claimed by SEAM blocks
		std::unordered_set<unsigned> boundary_stamped; // edge ids claimed by BOUNDARY blocks
		std::unordered_set<unsigned> junction_stamped; // vertex ids claimed by JUNCTION blocks

		for (int pass_type : {3, 2, 1, 0})
		{
			for (const auto& blk : struct_blocks)
			{
				if (blk.type_id != pass_type) continue;
				for (unsigned elem_id : blk.elems)
				{
					if (pass_type == 3) // JUNCTION — stamp vertex struct_id
					{
						if (elem_id < nv_total && slab_mesh.vertices[elem_id].first)
						{
							slab_mesh.vertices[elem_id].second->struct_id = blk.struct_id;
							junction_stamped.insert(elem_id);
						}
						applyMain(elem_id, 3);
					}
					else if (pass_type == 2) // BOUNDARY
					{
						if (elem_id < slab_mesh.edges.size() && slab_mesh.edges[elem_id].first)
						{
							auto* e = slab_mesh.edges[elem_id].second;
							e->struct_id = blk.struct_id;
							boundary_stamped.insert(elem_id);
							applyBoundary(e->vertices_.first);
							applyBoundary(e->vertices_.second);
						}
					}
					else if (pass_type == 1) // SEAM (highest edge priority)
					{
						if (elem_id < slab_mesh.edges.size() && slab_mesh.edges[elem_id].first)
						{
							auto* e = slab_mesh.edges[elem_id].second;
							e->struct_id = blk.struct_id;
							seam_stamped.insert(elem_id);
							applyMain(e->vertices_.first,  2);
							applyMain(e->vertices_.second, 2);
						}
					}
					else // SHEET — face edges, skip if already claimed by seam/boundary
					{
						if (elem_id < slab_mesh.faces.size() && slab_mesh.faces[elem_id].first)
						{
							auto* face = slab_mesh.faces[elem_id].second;
							face->struct_id = blk.struct_id;
							for (unsigned vid : face->vertices_)
								applyMain(vid, 1);
							for (unsigned eid : face->edges_)
							{
								if (eid >= slab_mesh.edges.size() || !slab_mesh.edges[eid].first) continue;
								if (seam_stamped.count(eid) || boundary_stamped.count(eid)) continue;
								slab_mesh.edges[eid].second->struct_id = blk.struct_id;
							}
						}
					}
				}
			}
		}

		// Assign nmn_cluster_type from main_prio + has_boundary.
		unsigned n_sheet=0, n_seam=0, n_boundary=0, n_junction=0,
		         n_sheet_b=0, n_seam_b=0, n_junction_b=0, n_unknown=0;
		for (unsigned i = 0; i < nv_total; ++i)
		{
			if (!slab_mesh.vertices[i].first) continue;
			CT ct = CT::MS_Unknown;
			switch (main_prio[i]) {
				case 1: ct = has_boundary[i] ? CT::MS_Sheet_Boundary    : CT::MS_Sheet;    break;
				case 2: ct = has_boundary[i] ? CT::MS_Seam_Boundary     : CT::MS_Seam;     break;
				case 3: ct = has_boundary[i] ? CT::MS_Junction_Boundary : CT::MS_Junction; break;
				default: ct = has_boundary[i] ? CT::MS_Boundary : CT::MS_Unknown; break;
			}
			slab_mesh.vertices[i].second->nmn_cluster_type = ct;
			switch (ct) {
				case CT::MS_Sheet:            ++n_sheet;     break;
				case CT::MS_Seam:             ++n_seam;      break;
				case CT::MS_Boundary:         ++n_boundary;  break;
				case CT::MS_Junction:         ++n_junction;  break;
				case CT::MS_Sheet_Boundary:   ++n_sheet_b;   break;
				case CT::MS_Seam_Boundary:    ++n_seam_b;    break;
				case CT::MS_Junction_Boundary:++n_junction_b;break;
				default:                      ++n_unknown;   break;
			}
		}

		// ── PASS 2: matStruc_struct_collapsible ──────────────────────────────
		// Vertex types are now known.  Rules differ by edge type:
		//   Sheet edges:         collapsible only if BOTH endpoints are MS_Sheet
		//   Seam/boundary edges: collapsible only if NEITHER endpoint is
		//                        MS_Seam_Boundary, MS_Junction, or MS_Junction_Boundary
		//   Non-struct edges:    always false
		for (unsigned i = 0; i < slab_mesh.edges.size(); ++i)
		{
			if (!slab_mesh.edges[i].first) continue;
			auto* e = slab_mesh.edges[i].second;
			if (e->struct_id < 0) { e->matStruc_struct_collapsible = false; continue; }

			const CT va = slab_mesh.vertices[e->vertices_.first].second->nmn_cluster_type;
			const CT vb = slab_mesh.vertices[e->vertices_.second].second->nmn_cluster_type;

			if (seam_stamped.count(i) || boundary_stamped.count(i))
			{
				// Seam or boundary edge: blocked by junction or MS_Seam_Boundary endpoints
				e->matStruc_struct_collapsible =
					(va != CT::MS_Seam_Boundary) && (va != CT::MS_Junction) && (va != CT::MS_Junction_Boundary) &&
					(vb != CT::MS_Seam_Boundary) && (vb != CT::MS_Junction) && (vb != CT::MS_Junction_Boundary);
			}
			else
			{
				// Sheet edge: both endpoints must be pure MS_Sheet
				e->matStruc_struct_collapsible =
					(va == CT::MS_Sheet) &&
					(vb == CT::MS_Sheet);
			}
		}

		std::cout << "[LoadMatstructMA] struct data read from: " << fname << "\n"
		          << "  structs: " << num_structs
		          << " (sheet=" << sheet_structs << " seam=" << seam_structs
		          << " boundary=" << boundary_structs << " junction=" << junction_structs << ")\n"
		          << "  vertex types:"
		          << " sheet=" << n_sheet
		          << " sheet_b=" << n_sheet_b
		          << " seam=" << n_seam
		          << " seam_b=" << n_seam_b
		          << " boundary=" << n_boundary
		          << " junction=" << n_junction
		          << " junction_b=" << n_junction_b
		          << " unknown=" << n_unknown << "\n";
	}

	// ── Finalize slab mesh ────────────────────────────────────────────────────
	slab_mesh.CleanIsolatedVertices();
	slab_mesh.computebb();
	slab_mesh.ComputeFacesCentroid();
	slab_mesh.ComputeFacesNormal();
	slab_mesh.ComputeVerticesNormal();
	slab_mesh.ComputeEdgesCone();
	slab_mesh.ComputeFacesSimpleTriangles();
	slab_mesh.DistinguishVertexType();
}
#endif // USE_MATSTRUCT_INITIALIZATION

void ThreeDimensionalShape::ComputeFeatureEdges()
{
	// Build flat vertex position list from input mesh
	std::vector<pre_matfp::Vector3> verts;
	verts.reserve(input.pVertexList.size());
	for (unsigned i = 0; i < input.pVertexList.size(); ++i)
	{
		const auto& p = input.pVertexList[i]->point();
		verts.emplace_back(p[0], p[1], p[2]);
	}

	// Build flat face list using stored vertex IDs
	std::vector<pre_matfp::Vector3i> faces;
	for (auto fi = input.facets_begin(); fi != input.facets_end(); ++fi)
	{
		auto h = fi->facet_begin();
		int v0 = h->vertex()->id;
		int v1 = h->next()->vertex()->id;
		int v2 = h->next()->next()->vertex()->id;
		faces.emplace_back(v0, v1, v2);
	}

	matfp::Args args;  // default thres_concave=0.18, thres_convex=30.0
	slab_mesh.sharp_edges.clear();
	slab_mesh.concave_edges.clear();
	slab_mesh.feature_corners.clear();
	pre_matfp::find_feature_edges(args, verts, faces,
	                               slab_mesh.sharp_edges,
	                               slab_mesh.concave_edges,
	                               slab_mesh.feature_corners);
}

long ThreeDimensionalShape::LoadSlabMesh()
{
	slab_mesh.clear();
	long startt = clock();
	InitialSlabMesh();
	slab_mesh.initCollapseQueue();
	long endt = clock();
	return endt - startt;

	//if (slab_mesh.clear_error)
	//{
	//	slab_mesh.clear();
	//	long startt = clock();
	//	while(!slab_mesh.boundary_edge_collapses_queue.empty())
	//		slab_mesh.boundary_edge_collapses_queue.pop();

	//	InitialSlabMesh();
	//	slab_mesh.initCollapseQueue();
	//	long endt = clock();
	//	return endt - startt;
	//}
	//else
	//{

	//	slab_mesh.RecomputerVertexType();	
	//	long startt = clock();
	//	while(!slab_mesh.boundary_edge_collapses_queue.empty())
	//		slab_mesh.boundary_edge_collapses_queue.pop();
	//	if (slab_initial == false)
	//	{
	//		switch(slab_mesh.hyperbolic_weight_type)
	//		{
	//		case 1:
	//			InitialSlabMesh();
	//			break;
	//		case 2:
	//			InitialWeightedSlabMesh();
	//			break;
	//		case 3:
	//			InitialSlabMesh();
	//			break;
	//		default:
	//			InitialSlabMesh();
	//			break;
	//		}
	//		slab_initial = true;
	//	}
	//	slab_mesh.initCollapseQueue();
	//	long endt = clock();
	//	return endt - startt;
	//}
}

void ThreeDimensionalShape::InitialSlabMesh()
{
	// handle each face
	for(unsigned i = 0; i < slab_mesh.vertices.size(); i++)
	{ 
		if(!slab_mesh.vertices[i].first)
			continue;

		SlabVertex sv = *slab_mesh.vertices[i].second;
		std::set<unsigned> fset = sv.faces_;
		Vector4d C1(sv.sphere.center.X(), sv.sphere.center.Y(), sv.sphere.center.Z(), sv.sphere.radius);

		for (set<unsigned>::iterator si = fset.begin(); si != fset.end(); si++)
		{
			SlabFace sf = *slab_mesh.faces[*si].second;

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

			slab_mesh.vertices[i].second->slab_A += temp_A1;
			slab_mesh.vertices[i].second->slab_A += temp_A2;
			slab_mesh.vertices[i].second->slab_b += temp_b1;
			slab_mesh.vertices[i].second->slab_b += temp_b2;
			slab_mesh.vertices[i].second->slab_c += temp_c1;
			slab_mesh.vertices[i].second->slab_c += temp_c2;

			slab_mesh.vertices[i].second->related_face += 2;
		}
	}

	switch(slab_mesh.preserve_boundary_method)
	{
	case 1 :
		slab_mesh.PreservBoundaryMethodOne();
		break;
	case 2 :
		//slab_mesh.PreservBoundaryMethodTwo();
		break;
	case 3 :
		slab_mesh.PreservBoundaryMethodThree();
		break;
	default:
		slab_mesh.PreservBoundaryMethodFour();
		break;
	}

}

double ThreeDimensionalShape::NearestPoint(Vector3d point, unsigned vid)
{
	set<unsigned> faces = slab_mesh.vertices[vid].second->faces_;
	set<unsigned> edges = slab_mesh.vertices[vid].second->edges_;
	double mind = DBL_MAX;

	// calculation of related faces
	for (set<unsigned>::iterator si = faces.begin(); si != faces.end(); si++)
	{
		if (!slab_mesh.faces[*si].first)
			continue;
		SlabFace sf = *slab_mesh.faces[*si].second;
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
	for (set<unsigned>::iterator si = edges.begin(); si != edges.end(); si++)
	{
		if (!slab_mesh.edges[*si].first)
			continue;
		SlabEdge se = *slab_mesh.edges[*si].second;
		if (se.valid_cone == false)
			continue;

		Vector3d tfp;
		double td;
		se.cone.ProjectOntoCone(point, tfp, td);
		if (td < 0)
			td = -td;

		if(td < mind)	mind = td;
	}

	return mind;
}

void ThreeDimensionalShape::ComputeHausdorffDistance()
{
#if 0
	long start_time = clock();

	ma_qem_mesh.maxhausdorff_distance = 0.;
	for(unsigned i = 0; i < ma_qem_mesh.faces.size(); i ++)
	{
		Vector3d ve[8];
		if (ma_qem_mesh.faces[i].second->valid_st == false ||
			ma_qem_mesh.faces[i].second->st[0].normal == Vector3d(0., 0., 0.) || 
			ma_qem_mesh.faces[i].second->st[1].normal == Vector3d(0., 0., 0.))
			continue;

		ve[0] = ma_qem_mesh.faces[i].second->st[0].v[0];
		ve[1] = ma_qem_mesh.faces[i].second->st[0].v[1];
		ve[2] = ma_qem_mesh.faces[i].second->st[0].v[2];
		ve[3] = ma_qem_mesh.faces[i].second->st[1].v[0];
		ve[4] = ma_qem_mesh.faces[i].second->st[1].v[1];
		ve[5] = ma_qem_mesh.faces[i].second->st[1].v[2];
		ve[6] = (ve[0] + ve[1] +ve[2]) / 3.0;
		ve[7] = (ve[3] + ve[4] +ve[5]) / 3.0;
		double face_haus = 0.;
		for (int j = 0; j < 8; j++)
		{
			Vector3d fp = input.NearestVertex(ve[j]);
			double len = (ve[j] - fp).Length();
			face_haus = max(len, face_haus);
		}
		ma_qem_mesh.faces[i].second->hausdorff_dist = face_haus;

		ma_qem_mesh.maxhausdorff_distance = max(ma_qem_mesh.maxhausdorff_distance,face_haus);
	}

	long end_time = clock();
	long result = end_time - start_time;
#endif

	//ma_qem_mesh.maxhausdorff_distance = 0.;
	slab_mesh.maxhausdorff_distance = 0;
	double sumhausdorff_distance = 0;
	for (unsigned i = 0; i < input.pVertexList.size(); i++)
	{
		double min_dis = DBL_MAX;
		unsigned min_index = -1;
		Vector3d bou_ver(input.pVertexList[i]->point()[0], input.pVertexList[i]->point()[1], input.pVertexList[i]->point()[2]);
		bou_ver /=  input.bb_diagonal_length; 

		for (unsigned j = 0; j < slab_mesh.numVertices; j++)
		{
			Sphere ma_ver = slab_mesh.vertices[j].second->sphere;
			double temp_length = abs((bou_ver - ma_ver.center).Length() - ma_ver.radius);
			//if (temp_length >= 0 && temp_length < min_dis)
			if (temp_length < min_dis)
			{
				min_dis = temp_length;
				min_index = j;
			}

			//double temp_near_dis = slab_mesh.NearestPoint(bou_ver, min_index);
			//if (temp_near_dis < min_dis)
			//{
			//	min_dis = temp_near_dis;
			//	min_index = j;
			//}

		}

		//// Ϊ������ó��Ľ����ǰ��ó��Ľ����ҪС��
		//double nearest_dis = NearestPoint(bou_ver, min_index);
		//min_dis = min(nearest_dis, min_dis);

		//sumhausdorff_distance += min_dis;


		if (min_index != -1)
		{	
			double temp_near_dis = slab_mesh.NearestPoint(bou_ver, min_index);
			min_dis = min(temp_near_dis, min_dis);

			sumhausdorff_distance += min_dis;

			//ma_qem_mesh.vertices[min_index].second->bplist.push_back(i);
			//ma_qem_mesh.maxhausdorff_distance = max(ma_qem_mesh.maxhausdorff_distance, min_dis);

			slab_mesh.vertices[min_index].second->bplist.insert(i);
			slab_mesh.maxhausdorff_distance = max(slab_mesh.maxhausdorff_distance, min_dis);

			//input.pVertexList[i]->vqem_hausdorff_dist = min_dis / input.bb_diagonal_length;
			input.pVertexList[i]->vqem_hausdorff_dist = min_dis;
			input.pVertexList[i]->vqem_hansdorff_index = min_index;

			//input.pVertexList[i]->slab_hausdorff_dist = min_dis / input.bb_diagonal_length;
			input.pVertexList[i]->slab_hausdorff_dist = min_dis;
			input.pVertexList[i]->slab_hansdorff_index = min_index;
		}
	}
	
	//ma_qem_mesh.meanhausdorff_distance = sumhausdorff_distance / input.pVertexList.size();
	slab_mesh.meanhausdorff_distance = sumhausdorff_distance / input.pVertexList.size();
	//ma_qem_mesh.initialhausdorff_distance = ma_qem_mesh.maxhausdorff_distance;
	slab_mesh.initialhausdorff_distance = slab_mesh.maxhausdorff_distance;
}

void ThreeDimensionalShape::PruningSlabMesh()
{
	slab_mesh.ComputeVerticesProperty();

	bool has_boundary_non_pole;
	do
	{
		has_boundary_non_pole = false;
		unsigned vid;
		for(unsigned i = 0; i < slab_mesh.vertices.size(); i ++)
			if(slab_mesh.vertices[i].first)
				if((slab_mesh.vertices[i].second->is_boundary) &&
					(!slab_mesh.vertices[i].second->is_non_manifold) &&
					(!slab_mesh.vertices[i].second->is_pole) &&
					(slab_mesh.vertices[i].second->edges_.size() == 2))
				{
					vid = i;
					has_boundary_non_pole = true;
					break;
				}
				if(has_boundary_non_pole)
					slab_mesh.DeleteVertex(vid);
	}while(has_boundary_non_pole);


	slab_mesh.CleanIsolatedVertices();
	slab_mesh.ComputeEdgesCone();
	slab_mesh.ComputeFacesSimpleTriangles();
	slab_mesh.DistinguishVertexType();
	slab_mesh.computebb();
}

// ─────────────────────────────────────────────────────────────────────────────
// DetermineTopology
// ─────────────────────────────────────────────────────────────────────────────
//
// Classifies every active MAT vertex by its local topological role, determined
// by inspecting the face-valence of each incident edge:
//
//   edge face-count == 1  →  boundary edge
//   edge face-count == 2  →  manifold (sheet) edge
//   edge face-count  > 2  →  non-manifold (seam) edge
//
// Per-vertex flags set:
//   topo_is_sheet    – has edges, ALL of which are manifold (face-count == 2).
//   topo_is_seam     – has at least one non-manifold edge (face-count > 2).
//   topo_is_junction – has more than 2 seam edges (seams converge here).
//   topo_is_boundary – has at least one boundary edge (face-count == 1).

void ThreeDimensionalShape::DetermineTopology()
{
	unsigned n_sheet = 0, n_seam = 0, n_junction = 0, n_boundary = 0;

	for (unsigned i = 0; i < input_nmm.vertices.size(); ++i)
	{
		if (!input_nmm.vertices[i].first) continue;
		NonManifoldMesh_Vertex* v = input_nmm.vertices[i].second;

		// initialise
		v->topo_is_sheet    = false;
		v->topo_is_seam     = false;
		v->topo_is_junction = false;
		v->topo_is_boundary = false;

		unsigned seam_edge_count = 0;
		bool     has_any_edge    = false;
		bool     all_manifold    = true;

		for (unsigned eid : v->edges_)
		{
			if (eid >= input_nmm.edges.size() || !input_nmm.edges[eid].first)
				continue;

			has_any_edge = true;
			unsigned nf = static_cast<unsigned>(input_nmm.edges[eid].second->faces_.size());

			if (nf == 1)
			{
				v->topo_is_boundary = true;
				all_manifold = false;
			}
			else if (nf > 2)
			{
				v->topo_is_seam = true;
				all_manifold = false;
				++seam_edge_count;
			}
			// nf == 2: manifold edge – no flag needed
		}

		if (seam_edge_count > 2)
			v->topo_is_junction = true;

		// Sheet: has edges and every one of them is 2-manifold
		if (has_any_edge && all_manifold)
			v->topo_is_sheet = true;

		if (v->topo_is_sheet)    ++n_sheet;
		if (v->topo_is_seam)     ++n_seam;
		if (v->topo_is_junction) ++n_junction;
		if (v->topo_is_boundary) ++n_boundary;
	}

	std::cout << "[DetermineTopology]"
	          << "  sheet="    << n_sheet
	          << "  seam="     << n_seam
	          << "  junction=" << n_junction
	          << "  boundary=" << n_boundary << "\n";
}

void ThreeDimensionalShape::ExportSurfacemeshFeatureEdges(const std::string& path) const
{
	std::string outpath = path + ".mesh_features";
	std::ofstream f(outpath);
	if (!f.is_open()) {
		std::cerr << "[ExportSurfacemeshFeatureEdges] Cannot open " << outpath << "\n";
		return;
	}

	f << "[convex]\n";
	for (const auto& e : slab_mesh.sharp_edges)
		f << e[0] << " " << e[1] << "\n";

	f << "[concave]\n";
	for (const auto& e : slab_mesh.concave_edges)
		f << e[0] << " " << e[1] << "\n";

	std::cout << "[ExportSurfacemeshFeatureEdges] Written to " << outpath
	          << "  convex=" << slab_mesh.sharp_edges.size()
	          << "  concave=" << slab_mesh.concave_edges.size() << "\n";
}

