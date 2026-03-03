#include "ThreeDimensionalShape.h"

// Note: QString include removed - was unused and prevents CLI build without Qt

#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/property_map.h>
#include <numeric>
#include <unordered_map>

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

	// Save vertex-to-sample-point associations:
	// For each Voronoi vertex (circumcenter), record its sphere and the IDs + coordinates
	// of the 4 Delaunay vertices (original sample points) that define its tetrahedron.
	{
		std::string sample_fname = input_nmm.meshname + "_vertex_samples.txt";
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
// ClusterBoundaryPoints
// ─────────────────────────────────────────────────────────────────────────────
//
// Clusters input mesh vertices (boundary sample points) by surface orientation
// using a two-pass union-find on a k-NN spatial graph.
//
// Pass 1  – strict angle threshold (angle_threshold_deg, default 25°)
//   Connects neighbours whose normals differ by less than the threshold.
//   Reliably separates cube faces (adjacent face normals differ by 90°).
//   Produces clean face clusters + many small edge/corner fragments.
//
// Pass 2  – loose angle threshold (edge_merge_deg, default 60°)
//   Merges the small fragments from pass 1 together, but NEVER into a
//   large (face) cluster.  A cluster is considered "large" if it has at
//   least  max(3, n/40)  vertices after pass 1.
//   Reason for two passes: edge vertices may have slightly varying normals
//   (due to unequal triangle counts on each side of a mesh edge), causing
//   them to exceed the strict threshold even though they belong to the
//   same geometric edge.  The loose threshold re-connects these fragments
//   without merging them into the face clusters.
//
// Prerequisites
//   • input.compute_normals() must have been called (done in main_cli.cpp).
//   • input_nmm.meshname must be set so the output file can be named.
//
// Parameters
//   k_neighbors         – spatial kd-tree neighbourhood size (default 12).
//   angle_threshold_deg – strict normal angle for pass-1 face clustering (25°).
//   edge_merge_deg      – loose normal angle for pass-2 edge merging (60°).
//
// Result
//   Vector-of-vectors sorted largest-first.
//   Indices refer to input.pVertexList (same order as BoundaryPoints).
//
// Output file  <meshname>_boundary_clusters.txt
//   line 1-3 : comments
//   line 4   : <num_clusters>  <num_points>
//   remaining: <cluster_id>  <vertex_index>  <x>  <y>  <z>  <nx>  <ny>  <nz>

std::vector<std::vector<unsigned>>
ThreeDimensionalShape::ClusterBoundaryPoints(int    k_neighbors,
                                             double angle_threshold_deg,
                                             double edge_merge_deg)
{
	typedef simple_kernel::Point_3  Pt;
	typedef simple_kernel::Vector_3 Vec;

	// ── kd-tree with per-point index ──────────────────────────────────────────
	typedef std::pair<Pt, unsigned>                               PtId;
	typedef CGAL::First_of_pair_property_map<PtId>               PtMap;
	typedef CGAL::Search_traits_3<simple_kernel>                 BaseTraits;
	typedef CGAL::Search_traits_adapter<PtId, PtMap, BaseTraits> SearchTraits;
	typedef CGAL::Orthogonal_k_neighbor_search<SearchTraits>     KNN;
	typedef KNN::Tree                                             KDTree;

	const unsigned n = static_cast<unsigned>(input.pVertexList.size());
	if (n == 0) return {};

	std::vector<PtId> pts(n);
	for (unsigned i = 0; i < n; ++i)
		pts[i] = { input.pVertexList[i]->point(), i };

	KDTree tree(pts.begin(), pts.end());

	// ── union-find helpers (path-halving) ─────────────────────────────────────
	std::vector<unsigned> parent(n);
	std::iota(parent.begin(), parent.end(), 0u);

	auto find = [&](unsigned x) -> unsigned {
		while (parent[x] != x) {
			parent[x] = parent[parent[x]];
			x = parent[x];
		}
		return x;
	};
	auto unite = [&](unsigned a, unsigned b) {
		a = find(a); b = find(b);
		if (a != b) parent[a] = b;
	};

	// ── pass 1: strict threshold – separates face patches ─────────────────────
	const double cos_strict = std::cos(angle_threshold_deg * CGAL_PI / 180.0);

	for (unsigned i = 0; i < n; ++i)
	{
		const Vec& ni = input.pVertexList[i]->normal;
		if (ni.squared_length() < 1e-20) continue;

		KNN search(tree, input.pVertexList[i]->point(), k_neighbors + 1);
		for (auto it = search.begin(); it != search.end(); ++it)
		{
			unsigned j = it->first.second;
			if (j == i) continue;
			const Vec& nj = input.pVertexList[j]->normal;
			if (nj.squared_length() < 1e-20) continue;
			if ((ni * nj) >= cos_strict)
				unite(i, j);
		}
	}

	// ── pass 2: loose threshold – merge edge/corner fragments ─────────────────
	// A cluster is "large" (a face patch) if it has >= min_size vertices.
	// We only merge pairs of vertices that BOTH belong to small clusters,
	// so face patches are never altered in this pass.
	const unsigned min_size   = static_cast<unsigned>(std::max(3, (int)n / 40));
	const double   cos_loose  = std::cos(edge_merge_deg * CGAL_PI / 180.0);

	// count cluster sizes after pass 1
	std::unordered_map<unsigned, unsigned> root_count;
	for (unsigned i = 0; i < n; ++i)
		root_count[find(i)]++;

	// mark which vertices belong to a large (face) cluster
	std::vector<bool> in_large(n);
	for (unsigned i = 0; i < n; ++i)
		in_large[i] = (root_count[find(i)] >= min_size);

	for (unsigned i = 0; i < n; ++i)
	{
		if (in_large[i]) continue;                      // skip face-cluster vertices
		const Vec& ni = input.pVertexList[i]->normal;
		if (ni.squared_length() < 1e-20) continue;

		KNN search(tree, input.pVertexList[i]->point(), k_neighbors + 1);
		for (auto it = search.begin(); it != search.end(); ++it)
		{
			unsigned j = it->first.second;
			if (j == i || in_large[j]) continue;        // never merge into a face cluster
			const Vec& nj = input.pVertexList[j]->normal;
			if (nj.squared_length() < 1e-20) continue;
			if ((ni * nj) >= cos_loose)
				unite(i, j);
		}
	}

	// ── collect final clusters ────────────────────────────────────────────────
	std::unordered_map<unsigned, std::vector<unsigned>> cmap;
	for (unsigned i = 0; i < n; ++i)
		cmap[find(i)].push_back(i);

	std::vector<std::vector<unsigned>> result;
	result.reserve(cmap.size());
	for (auto& kv : cmap)
		result.push_back(std::move(kv.second));

	std::sort(result.begin(), result.end(),
	          [](const auto& a, const auto& b){ return a.size() > b.size(); });

	std::cout << "[ClusterBoundaryPoints]  " << result.size() << " clusters"
	          << "  (k=" << k_neighbors
	          << ", pass1=" << angle_threshold_deg << "deg"
	          << ", pass2=" << edge_merge_deg << "deg"
	          << ", min_size=" << min_size
	          << ", points=" << n << ")\n";

	// ── write file ────────────────────────────────────────────────────────────
	std::string fname = input_nmm.meshname + "_boundary_clusters.txt";
	std::ofstream ofs(fname.c_str());
	if (ofs.is_open())
	{
		ofs << "# QMAT boundary point clusters\n";
		ofs << "# k_neighbors=" << k_neighbors
		    << "  pass1_angle=" << angle_threshold_deg
		    << "  pass2_edge_merge=" << edge_merge_deg
		    << "  min_cluster_size=" << min_size << "\n";
		ofs << "# cluster_id  vertex_index  x  y  z  nx  ny  nz\n";
		ofs << result.size() << "  " << n << "\n";

		for (std::size_t cid = 0; cid < result.size(); ++cid)
		{
			for (unsigned vi : result[cid])
			{
				auto p  = input.pVertexList[vi]->point();
				auto nm = input.pVertexList[vi]->normal;
				ofs << cid << "  " << vi
				    << "  " << p.x()  << "  " << p.y()  << "  " << p.z()
				    << "  " << nm.x() << "  " << nm.y() << "  " << nm.z()
				    << "\n";
			}
		}
		ofs.close();
		std::cout << "[ClusterBoundaryPoints]  written -> " << fname << "\n";
	}
	else
	{
		std::cerr << "[ClusterBoundaryPoints]  WARNING: could not write "
		          << fname << "\n";
	}

	return result;
}