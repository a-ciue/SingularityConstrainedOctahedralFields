/*
 * MatchingsAndAlignmentExtractionT.cc
 *
 *  Created on: Feb 28, 2018
 *      Author: hliu
 */
#define MATCHINGSANDALIGNMENTEXTRACTIONT_C


#include "MatchingsAndAlignmentExtractionT.hh"

namespace SCOF
{

template <class MeshT1, class MeshT2>
void
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
determineMatchings()
{
	if(!check_input())
		return;

	StopWatch<> time, time_total;
	time_total.start();
	time.start();

	//0. initialize properties
	initialize_output_properties();
	initialize_trimesh_and_properties();
	set_singular_arc_alignment_constraints();
	//1.
	std::cerr << "\n######Merge singular arc charts...";
	merge_singular_arc_charts();
	std::cerr << "\nSingular arc merge time:    " <<time.stop()/1000<<" s";
	time.start();
	//2.
	std::cerr << "\n######Merge boundary charts...";
	merge_boundary_charts();
	std::cerr << "\nBoundary merge time:    " <<time.stop()/1000<<" s";
	time.start();
	//3.
	std::cerr << "\n######Merge volume charts...";
	merge_volume_charts();
	std::cerr << "\nVolume merge time:    " <<time.stop()/1000<<" s";
	std::cout << "\n######Total time:    " <<time_total.stop()/1000<<" s";
	//4.
	std::cerr<<"\n######Check transition and axis direction...";
	check_transitions();
	std::cerr<<"\n######Check done!"<< std::endl;
}


template <class MeshT1, class MeshT2>
void
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
initialize_singularity_graph()
{
	sg_.update_node_type_property();
	sg_.set_singularity_graph_component_property();
	sg_.update_singular_vertex_property();

	max_label_ = sg_.max_label();
	max_comp_i_= *std::max_element(sg_comp_i_.begin(), sg_comp_i_.end());
	max_comp_b_= *std::max_element(sg_comp_b_.begin(), sg_comp_b_.end());
}


template <class MeshT1, class MeshT2>
void
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
initialize_chart_index()
{
	for(int i=1; i<= max_label_; ++i)
		real_chart_id_[i] = i;
	real_chart_id_[-1] = -1;
	for(int i=1; i<= max_comp_b_; ++i)
		real_chart_id_[-i] = -i;
}


template <class MeshT1, class MeshT2>
void
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
initialize_output_properties()
{
	//initialize transition property
	for(auto hfi : mesh_.halffaces())
		trans_prop_[hfi] = -1;
	//initialize normal property(+x)
//	for(auto hfi : mesh_.halffaces())
//		if(mesh_.is_boundary(hfi))
//			normal_prop_[hfi] = 0;

	//initilaize axes direction property
	for(auto hei : mesh_.halfedges())
		for(auto hec_it = mesh_.hec_iter(hei); hec_it.valid(); ++hec_it)
			axes_prop_[*hec_it][hei] = -1;
}


template <class MeshT1, class MeshT2>
void
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
initialize_trimesh_and_properties()
{
	create_trimesh();

	trimesh_.add_property(tri_val_);
	trimesh_.add_property(tri_node_type_);
	trimesh_.add_property(tri_vh_trans_);

	set_trimesh_property();

	trimesh_.add_property(tri_trans_);
	trimesh_.add_property(tri_axis_);
}


template <class MeshT1, class MeshT2>
void
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
set_singular_arc_alignment_constraints()
{
	//interior singular edges +x
	for(auto ehi : mesh_.edges())
		if(!mesh_.is_boundary(ehi) && valance_[ehi] != 0)
		{
			auto he = mesh_.halfedge_handle(ehi, 0);
			for(auto hec_it = mesh_.hec_iter(he); hec_it.valid(); ++hec_it)
					axes_prop_[*hec_it][he] = 0;
		}
	//at interior singular nodes
	for(auto vhi : mesh_.vertices())
		if(node_type_[vhi] > 0)
		{
			for(auto voh_it = mesh_.voh_iter(vhi); voh_it.valid(); ++voh_it)
			{
				auto eho = mesh_.edge_handle(*voh_it);
				if(!mesh_.is_boundary(eho) && valance_[eho] != 0)
				{
					for(auto hec_it = mesh_.hec_iter(*voh_it); hec_it.valid(); ++hec_it)
						axes_prop_[*hec_it][*voh_it] = 0;

					auto he_opp = mesh_.opposite_halfedge_handle(*voh_it);
					for(auto hec_it = mesh_.hec_iter(he_opp); hec_it.valid(); ++hec_it)
						axes_prop_[*hec_it][he_opp] = -1;
				}
			}
		}

	//boundary singular edges +y
	for(auto ehi : trimesh_.edges())
	{
		auto he = trimesh_.halfedge_handle(ehi, 0);
		auto he_opp = trimesh_.halfedge_handle(ehi, 1);
		if(trimesh_.property(tri_val_, ehi) != 0)
		{
			trimesh_.property(tri_axis_, he) = 2;
			trimesh_.property(tri_axis_, he_opp) = -1;
		}else
		{
			trimesh_.property(tri_axis_, he) = -1;
			trimesh_.property(tri_axis_, he_opp) = -1;
		}
	}

	//at boundary singular nodes
	for(auto vhi : trimesh_.vertices())
		if(trimesh_.property(tri_node_type_, vhi) < 0)
		{
			for(auto voh_it = trimesh_.voh_iter(vhi); voh_it.is_valid(); ++voh_it)
			{
				int val = trimesh_.property(tri_val_, trimesh_.edge_handle(*voh_it));
				if(val != 0)
				{
					trimesh_.property(tri_axis_, *voh_it) = 2;
					trimesh_.property(tri_axis_, trimesh_.opposite_halfedge_handle(*voh_it)) = -1;
				}
			}
		}
}


template <class MeshT1, class MeshT2>
void
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
add_corner_constraints(std::map<EH, std::vector<int> >& _eh_to_corners)
{
	MP<std::vector<VH> > corner_vt_prop = mesh_.template request_mesh_property<std::vector<VH> >("corner_vertices");
	MP<std::vector<std::vector<HFH> > > dp_prop = mesh_.template request_mesh_property<std::vector<std::vector<HFH> > >("corner_dual_paths");
	mesh_.set_persistent(corner_vt_prop, false);
	mesh_.set_persistent(dp_prop, false);

	unsigned int n_corner = corner_vt_prop[0].size()/4;
	if(n_corner == 0)
		return;
	std::cerr<<"\nInitialize Corner Constraints...";

	for(unsigned int i=0; i<n_corner; ++i)
	{
		std::shared_ptr<CornerT<TetMesh> > corner(new CornerT<TetMesh>());
		for(int j=0; j<3; ++j)
		{
			HEH hej = mesh_.halfedge(corner_vt_prop[0][4*i], corner_vt_prop[0][4*i+j+1]);
			corner.get()->set_ith_halfedge(j, hej);
			corner.get()->set_ith_dual_path(j, dp_prop[0][3*i+j]);
		}

		v_corners_.push_back(corner);
	}
	//set up mapping: edge handle to corners' indices
	std::vector<EH> ehs;
	for(auto i : v_corners_)
		for(int j=0; j<3; ++j)
			ehs.push_back(mesh_.edge_handle(i.get()->halfedge(j)));

	std::sort(ehs.begin(), ehs.end());
	std::unique(ehs.begin(), ehs.end());
	for(auto ehi : ehs)
	{
		std::vector<int> vec;
		_eh_to_corners[ehi] = vec;
	}

	for(unsigned int i=0; i<n_corner; ++i)
		for(int j=0; j<3; ++j)
		{
			EH ehj = mesh_.edge_handle(v_corners_[i].get()->halfedge(j));
			_eh_to_corners[ehj].push_back(i);
		}

	//check input corners
	for(auto i : v_corners_)
	{
		std::vector<CH> chs;
		for(int j=0; j<3; ++j)
		{
			std::vector<HFH> dp;
			 i.get()->dual_path(j, dp);
			 CH chj = mesh_.incident_cell(dp[dp.size()-1]);
			 if(!chj.is_valid())
				 std::cerr<<"\nError: invalid corner dual path! Halfface: "<<dp[dp.size()-1];
			 chs.push_back(chj);
		}

		if(chs[0] != chs[1] || chs[0] != chs[2])
		{
			std::cerr<<"\nError: corner dual path is wrong! Node: "<<i.get()->node(mesh_)
					<<" chs: "<<chs[0]<<" "<<chs[1]<<" "<<chs[2];
		}
	}

	std::cerr<<"\nInitialize corners' dual path properties...";
	for(auto i : v_corners_)
	{
		if(!i.get()->initially_set_transitions(mesh_, cell_tag_, arc_type_, trans_prop_, chart_id_))
		{
			std::cerr<<"\nNode: "<<i.get()->node(mesh_);
			std::cerr<<"\nError: couldn't correctly initialize corners' dual paths properties!";
			return;
		}
	}
}


template <class MeshT1, class MeshT2>
void
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
add_tangent_continuity_constraints()
{
	MP<std::vector<VH> > tc_vt_prop = mesh_.template request_mesh_property<std::vector<VH> >("tangent_continuity_vertices");
	MP<std::vector<std::vector<HFH> > > tc_dp_prop = mesh_.template request_mesh_property<std::vector<std::vector<HFH> > >("tangent_continuity_dual_paths");
	mesh_.set_persistent(tc_vt_prop, false);
	mesh_.set_persistent(tc_dp_prop, false);

	unsigned int n_tc = tc_vt_prop[0].size()/3;
	if(n_tc == 0)
		return;
	std::cerr<<"\nInitialize Tangent Continuity Constraints of Singular Nodes...";

	for(unsigned int i=0; i<n_tc; ++i)
	{
		if(node_type_[tc_vt_prop[0][3*i]] == 4)
		{
			std::shared_ptr<TangentContinuityT<TetMesh> > tc(new TangentContinuityT<TetMesh>());
			for(int j=0; j<2; ++j)
			{
				HEH hej = mesh_.halfedge(tc_vt_prop[0][3*i], tc_vt_prop[0][3*i+j+1]);
				tc.get()->set_ith_halfedge(j, hej);
			}
			tc.get()->set_dual_path(tc_dp_prop[0][i]);

			v_tc_.push_back(tc);
		}
	}


	for(auto i : v_tc_)
	{
		if(!i.get()->initially_set_transitions(mesh_, cell_tag_, trans_prop_, chart_id_))
		{
			std::cerr<<"\nNode: "<<i.get()->node(mesh_);
			std::cerr<<"\nError: couldn't correctly initialize dual paths properties!";
			return;
		}
	}
}


template <class MeshT1, class MeshT2>
void
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
constrained_chart_mergingB()
{
	std::vector<bool> bnd_sgl_vh(trimesh_.n_vertices(), false);
	for(auto ehi : trimesh_.edges())
		if(trimesh_.property(tri_val_, ehi) != 0)
		{
			auto tri_he0 = trimesh_.halfedge_handle(ehi, 0);
			bnd_sgl_vh[trimesh_.from_vertex_handle(tri_he0).idx()] = true;
			bnd_sgl_vh[trimesh_.to_vertex_handle(tri_he0).idx()] = true;
		}
	for(auto vhi : trimesh_.vertices())
	{
		if(bnd_sgl_vh[(vhi).idx()] && trimesh_.property(tri_node_type_, vhi) == 0)
			add_flat_sector_constraints(vhi);
	}


	MP<std::vector<VH> > sector_vt_prop = mesh_.template request_mesh_property<std::vector<VH> >("sector_vertices");
	MP<std::vector<double> > sector_angle_prop = mesh_.template request_mesh_property<std::vector<double> >("sector_angles");
	mesh_.set_persistent(sector_vt_prop, false);
	mesh_.set_persistent(sector_angle_prop, false);

	unsigned int n_sector = sector_vt_prop[0].size()/3;
	for(unsigned int i=0; i<n_sector; ++i)
	{
		std::vector<TriVH> tri_vhs;
		for(int j=0; j<3; ++j)
			tri_vhs.push_back(vh_to_trivh_[sector_vt_prop[0][3*i+j]]);

		add_irregular_sector_constraints(tri_vhs, sector_angle_prop[0][i]);
	}
}


template <class MeshT1, class MeshT2>
void
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
merge_singular_arc_charts()
{
	merge_and_zip_charts();

	constrained_chart_mergingT();

	constrained_chart_mergingB();
}


template <class MeshT1, class MeshT2>
void
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
merge_and_zip_charts()
{
//interior singular edges
	for(auto ehi : mesh_.edges())
		if(!mesh_.is_boundary(ehi) && valance_[ehi] != 0)
			interior_singular_edge_closing(ehi);

	//set chart id
	for(auto ehi : mesh_.edges())
		if(!mesh_.is_boundary(ehi) && valance_[ehi] != 0)
		{
			auto he = mesh_.halfedge_handle(ehi, 0);
			for(auto hec_it = mesh_.hec_iter(he); hec_it.valid(); ++hec_it)
			{
				chart_id_[*hec_it] = label_[ehi];
				cell_tag_[*hec_it] = arc_type_[ehi];
			}
			edge_tag_[ehi] = true;
		}

//boundary singular edges
	for(auto ehi : trimesh_.edges())
	{
		auto he = trimesh_.halfedge_handle(ehi, 0);
		auto he_opp = trimesh_.halfedge_handle(ehi, 1);
		if(trimesh_.property(tri_val_, ehi) == -1)
		{
			trimesh_.property(tri_trans_, he) = 6;
			trimesh_.property(tri_trans_, he_opp) = 7;
		}else if(trimesh_.property(tri_val_, ehi) == 1)
		{
			trimesh_.property(tri_trans_, he) = 7;
			trimesh_.property(tri_trans_, he_opp) = 6;
		}else
		{
			trimesh_.property(tri_trans_, he) = -1;
			trimesh_.property(tri_trans_, he_opp) = -1;
		}
	}

	//outgoing edges of singular nodes
	for(auto vhi : trimesh_.vertices())
		if(trimesh_.property(tri_node_type_, vhi) < 0)
		{
			for(auto voh_it = trimesh_.voh_iter(vhi); voh_it.is_valid(); ++voh_it)
			{
				int val = trimesh_.property(tri_val_, trimesh_.edge_handle(*voh_it));
				if(val == -1)
				{
					trimesh_.property(tri_trans_, *voh_it) = 6;
					trimesh_.property(tri_trans_, trimesh_.opposite_halfedge_handle(*voh_it)) = 7;
				}else if(val == 1)
				{
					trimesh_.property(tri_trans_, *voh_it) = 7;
					trimesh_.property(tri_trans_, trimesh_.opposite_halfedge_handle(*voh_it)) = 6;
				}
			}
		}
}


template <class MeshT1, class MeshT2>
void
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
merge_boundary_charts()
{
	//initialize chart index
	std::vector<int> chart_index(trimesh_.n_faces(), 0);
	std::vector<int> bnd_sgl_vh(trimesh_.n_vertices(), 0);
	for(auto ehi : trimesh_.edges())
		if(trimesh_.property(tri_val_, ehi) != 0)
		{
			auto tri_he0 = trimesh_.halfedge_handle(ehi, 0);
			auto he0 = trimesh_.property(tet_he_, tri_he0);
			bnd_sgl_vh[trimesh_.from_vertex_handle(tri_he0).idx()] = -sg_comp_b_[mesh_.edge_handle(he0)];
			bnd_sgl_vh[trimesh_.to_vertex_handle(tri_he0).idx()] = -sg_comp_b_[mesh_.edge_handle(he0)];
		}
	for(auto vhi : trimesh_.vertices())
		if(bnd_sgl_vh[vhi.idx()])
		{
			//get component number
			for(auto voh_it = trimesh_.voh_iter(vhi); voh_it.is_valid(); ++voh_it)
				chart_index[trimesh_.face_handle(*voh_it).idx()] = bnd_sgl_vh[vhi.idx()];
		}

	//push triangles incident to boundary singularities
	std::vector<TriFH> fhs;
	for(auto fhi : trimesh_.faces())
		if(chart_index[(fhi).idx()])
			fhs.push_back(fhi);
	//if no boundary singularity, randomly choose seed(s)
	if(fhs.empty())
	{
		get_seeds_of_boundaries(fhs);
		for(auto fhi : fhs)
			chart_index[fhi.idx()] = -1;
	}

	dual_spanning_forest_on_trimesh(fhs, chart_index);
	solve_transitions_on_trimesh();

//	std::cerr<<"\nTransfer axes direction... ";
	for(auto hei : trimesh_.halfedges())
		if(trimesh_.property(tri_axis_, hei) == 2)
		{
			auto he = trimesh_.property(tet_he_, hei);
			for(auto hec_it = mesh_.hec_iter(he); hec_it.valid(); ++hec_it)
				axes_prop_[*hec_it][he] = 2;
		}
	//connect boundary singular components(optional)
	if(max_comp_b_ > 1)
	{
		std::cerr<<"\nMax comp_b: "<<max_comp_b_;

		connect_boundary_graph_components(chart_index);
		solve_transitions_on_trimesh();
	}

	//check
	check_transitions_of_trimesh();

	//transfer property to tetmesh
	for(auto ehi : trimesh_.edges())
	{
		auto he0 = trimesh_.halfedge_handle(ehi, 0);
		auto he1 = trimesh_.halfedge_handle(ehi, 1);
		int trans0 = trimesh_.property(tri_trans_, he0);
		int trans1 = trimesh_.property(tri_trans_, he1);

		auto t_he = trimesh_.property(tet_he_, he0);
		auto hf_last = HFH(-1);
		for(auto hehf_it = mesh_.hehf_iter(t_he); hehf_it.valid(); ++hehf_it)
			if(!mesh_.is_boundary(mesh_.face_handle(*hehf_it)))
			{
				if(trans0 == -1)
				{
					trans_prop_[*hehf_it] = -1;
					trans_prop_[mesh_.opposite_halfface_handle(*hehf_it)] = -1;
				}
				else
				{
					trans_prop_[*hehf_it] = 0;
					trans_prop_[mesh_.opposite_halfface_handle(*hehf_it)] = 0;
				}
				hf_last = *hehf_it;
			}

		if(!hf_last.is_valid())
			std::cerr<<"\nError: cell has two boundary faces!";
		trans_prop_[hf_last] = trans0;
		trans_prop_[mesh_.opposite_halfface_handle(hf_last)] = trans1;

		if(trans0 == -1)
		{
			open_edges_.push_back(mesh_.edge_handle(t_he));
//			std::cerr<<" "<<mesh_.edge_handle(t_he);
		}else
			edge_tag_[mesh_.edge_handle(t_he)] = true;


//		if(trans0 !=0 )
//			std::cerr<<" "<<mesh_.edge_handle(t_he);
	}

	//conquer tets
	for(auto ehi : trimesh_.edges())
	{
		auto tri_he = trimesh_.halfedge_handle(ehi, 0);
		auto idx = chart_index[trimesh_.face_handle(tri_he).idx()];
		auto he = trimesh_.property(tet_he_, tri_he);
		if(trimesh_.property(tri_trans_, tri_he) != -1)
		{
			for(auto hec_it = mesh_.hec_iter(he); hec_it.valid(); ++hec_it)
			{
				cell_tag_[*hec_it] = 1;
				chart_id_[*hec_it] = idx;
			}
		}
	}

	zip_edges();
//	check_solved_regular_edges_transitions();
}


template <class MeshT1, class MeshT2>
bool
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
check_transitions_of_trimesh()
{
	for(auto vhi: trimesh_.vertices())
	{
		int idx = 0;
		bool all_trans_known = true;
		for(auto voh_it = trimesh_.voh_ccwiter(vhi); voh_it.is_valid(); ++voh_it)
		{
			int trans = trimesh_.property(tri_trans_, *voh_it);
			if(trans == -1)
			{
				all_trans_known = false;
				break;
			}

			idx = tq_.mult_transitions_idx(trans, idx);
		}
		if(all_trans_known)
		{
			int idx_ref = trimesh_.property(tri_vh_trans_, vhi);
			if(idx != idx_ref)
			{
				std::cerr<<"\nError: boundary shell transition is wrong! Transition product on edge: ";
				for(auto voh_it = trimesh_.voh_iter(vhi); voh_it.is_valid(); ++voh_it)
					std::cerr<<" "<<mesh_.edge_handle(trimesh_.property(tet_he_, *voh_it));
				std::cerr<<" should be: "<<idx_ref<<", but is: "<<idx;

				return false;
			}
		}
	}

	return true;
}


template <class MeshT1, class MeshT2>
void
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
solve_transitions_on_trimesh()
{
	//solve transition
	std::queue<TriVH> que_vhs;
	for(auto vhi : trimesh_.vertices())
		que_vhs.push(vhi);
	while(!que_vhs.empty())
	{
		TriHEH he_cur;
		auto vh_cur = que_vhs.front();
		que_vhs.pop();

		if(is_removable(vh_cur, he_cur))
		{
			int idx_r = trimesh_.property(tri_vh_trans_, vh_cur);
			int idx = 0;
			auto he_it = he_cur;
			while(1)
			{
				//ccw circulator
				he_it = trimesh_.prev_halfedge_handle(he_it);
				he_it = trimesh_.opposite_halfedge_handle(he_it);

				if(he_it == he_cur)
					break;

				int idxi = trimesh_.property(tri_trans_, he_it);
				idx = tq_.mult_transitions_idx(idxi, idx);
			}

			int idx_cur = tq_.mult_transitions_idx(tq_.inverse_transition_idx(idx), idx_r);

			trimesh_.property(tri_trans_, he_cur) = idx_cur;
			trimesh_.property(tri_trans_, trimesh_.opposite_halfedge_handle(he_cur)) = tq_.inverse_transition_idx(idx_cur);

			//push neighbour
			auto vh_t = trimesh_.to_vertex_handle(he_cur);
			que_vhs.push(vh_t);
		}
	}
}


template <class MeshT1, class MeshT2>
bool
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
is_removable(const TriVH _vh, TriHEH &_he)const
{
	int n_unknown = 0;
	for(auto vohccw_it = trimesh_.cvoh_ccwbegin(_vh); vohccw_it.is_valid(); ++vohccw_it)
		if(trimesh_.property(tri_trans_, *vohccw_it) == -1)
		{
			_he = *vohccw_it;
			n_unknown++;
		}

	if(n_unknown == 1)
		return true;
	return false;
}


template <class MeshT1, class MeshT2>
void
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
dual_spanning_forest_on_trimesh(const std::vector<TriFH>& _fhs, std::vector<int>& _chart_index)
{
	std::queue<TriFH> qu;
	for(auto fhi : _fhs)
		qu.push(fhi);

	while(!qu.empty())
	{
		auto cur_face = qu.front();
		qu.pop();

	    for(auto fhccw_it = trimesh_.fh_ccwbegin(cur_face); fhccw_it.is_valid(); ++fhccw_it)
	    {
	    		auto  neigh_fh  = trimesh_.face_handle(trimesh_.opposite_halfedge_handle(*fhccw_it));

	    		// face already part of spanning tree
	    		if( _chart_index[neigh_fh.idx()] )
	    			continue;
	    		else
	    		{
	    			trimesh_.property(tri_trans_, *fhccw_it) = 0;
	    			trimesh_.property(tri_trans_, trimesh_.opposite_halfedge_handle(*fhccw_it)) = 0;

	    			_chart_index[neigh_fh.idx()] = _chart_index[cur_face.idx()];
	    			qu.push(neigh_fh);
	    		}
	    }
	}
}


template <class MeshT1, class MeshT2>
void
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
get_seeds_of_boundaries(std::vector<TriFH>& _seeds)
{
	std::vector<bool> is_visited(trimesh_.n_faces(), false);
	for(auto fhi : trimesh_.faces())
	{
		if(!is_visited[fhi.idx()])
		{
			_seeds.push_back(fhi);

			std::queue<TriFH> que;
			is_visited[fhi.idx()] = true;
			que.push(fhi);
			while(!que.empty())
			{
				auto fh_cur = que.front();
				que.pop();

			    for(auto fh_it = trimesh_.fh_begin(fh_cur); fh_it.is_valid(); ++fh_it)
			    {
			    		auto neigh_fh = trimesh_.face_handle(trimesh_.opposite_halfedge_handle(*fh_it));

			    		// face already visited
			    		if( is_visited[neigh_fh.idx()] )
			    			continue;
			    		else
			    		{
			    			is_visited[neigh_fh.idx()] = true;
			    			que.push(neigh_fh);
			    		}
			    }
			}
		}
	}
}


template <class MeshT1, class MeshT2>
void
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
add_flat_sector_constraints(const TriVH _vh)
{
	std::vector<TriHEH> hes, hes_dir;
	for(auto voh_it = trimesh_.voh_begin(_vh); voh_it.is_valid(); ++voh_it)
	{
		auto eho = trimesh_.edge_handle(*voh_it);
		if(trimesh_.property(tri_val_, eho) != 0)
		{
			if(trimesh_.property(tri_axis_, *voh_it) == 2)
				hes_dir.push_back(*voh_it);
			else if(trimesh_.property(tri_axis_, trimesh_.opposite_halfedge_handle(*voh_it)) == 2)
				hes_dir.push_back(trimesh_.opposite_halfedge_handle(*voh_it));
			else
				std::cerr<<"\nError: no axis direction of the singular edge!";
			hes.push_back(*voh_it);
		}
	}

	if(hes.size() != 2)
		std::cerr<<"\nError: vertex has "<<hes.size()<<" outgoing singular edges!";


	bool same_dir = are_halfedges_the_same_direction(hes_dir[0], hes_dir[1]);
	int idx = 0;
	if(!same_dir)
		idx = 1;

	for(int i=0; i<2; ++i)
	{
		auto he_it = hes[i];
		TriHEH he_last(-1);
		while(1)
		{
			//ccw
			he_it = trimesh_.prev_halfedge_handle(he_it);
			he_it = trimesh_.opposite_halfedge_handle(he_it);

			if(he_it == hes[1-i])
				break;

			trimesh_.property(tri_trans_, he_it) = 0;
			trimesh_.property(tri_trans_, trimesh_.opposite_halfedge_handle(he_it)) = 0;
			he_last = he_it;
		}
		trimesh_.property(tri_trans_, he_last) = idx;
		trimesh_.property(tri_trans_, trimesh_.opposite_halfedge_handle(he_last)) = tq_.inverse_transition_idx(idx);
	}
}


template <class MeshT1, class MeshT2>
void
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
add_irregular_sector_constraints(const std::vector<TriVH>& _vhs, const double _angle)
{
	auto he0 = trimesh_.find_halfedge(_vhs[0], _vhs[1]);
	auto he1 = trimesh_.find_halfedge(_vhs[0], _vhs[2]);

	int trans = sector_transition(_angle);
//
	//assign transition
	auto he_it = he0;
	TriHEH he_last(-1);

	while(1)
	{
		//ccw
		he_it = trimesh_.prev_halfedge_handle(he_it);
		he_it = trimesh_.opposite_halfedge_handle(he_it);

		if(he_it == he1)
			break;

		trimesh_.property(tri_trans_, he_it) = 0;
		trimesh_.property(tri_trans_, trimesh_.opposite_halfedge_handle(he_it)) = 0;
		he_last = he_it;
	}
	trimesh_.property(tri_trans_, he_last) =  tq_.inverse_transition_idx(trans);
	trimesh_.property(tri_trans_, trimesh_.opposite_halfedge_handle(he_last)) = trans;
}


template <class MeshT1, class MeshT2>
int
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
sector_transition(const double _angle) const
{
	int angle = (int)_angle;
	int trans_idx = 0;
	if(angle == 0)
		trans_idx = 0;
	else if(angle == 90)
		trans_idx = 4;
	else if(angle == 180)
		trans_idx = 1;
	else if(angle == 270)
		trans_idx = 5;
	else if(angle == 360)
		trans_idx = 0;

	return trans_idx;
}


template <class MeshT1, class MeshT2>
bool
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
are_halfedges_the_same_direction(const TriHEH _he0, const TriHEH _he1)const
{
	auto vh00 = trimesh_.from_vertex_handle(_he0);
	auto vh01 = trimesh_.to_vertex_handle(_he0);

	auto vh10 = trimesh_.from_vertex_handle(_he1);
	auto vh11 = trimesh_.to_vertex_handle(_he1);

	if(vh01 == vh10 || vh11 == vh00)
		return true;

	return false;
}


template <class MeshT1, class MeshT2>
bool
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
are_halfedges_the_same_direction(const HEH _he0, const HEH _he1)const
{
	auto vh00 = mesh_.halfedge(_he0).from_vertex();
	auto vh01 = mesh_.halfedge(_he0).to_vertex();

	auto vh10 = mesh_.halfedge(_he1).from_vertex();
	auto vh11 = mesh_.halfedge(_he1).to_vertex();

	if(vh01 == vh10 || vh11 == vh00)
		return true;

	return false;
}


template <class MeshT1, class MeshT2>
double
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
distance(const VH _vh0, const VH _vh1)const
{
	auto vec = mesh_.vertex(_vh1) - mesh_.vertex(_vh0);
	return vec.norm();
}


template <class MeshT1, class MeshT2>
void
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
set_trimesh_property()
{
	for(auto vhi : mesh_.vertices())
		if(mesh_.is_boundary(vhi))
		{
			//get singular vertex transition
			int idx = 0;
			for(auto voh_it = mesh_.voh_iter(vhi); voh_it.valid(); ++voh_it)
			{
				auto eho = mesh_.edge_handle(*voh_it);
				if(!mesh_.is_boundary(eho))
				{
					if(valance_[eho] == -1)
						idx = 4;
					else if(valance_[eho] == 1)
						idx = 5;
				}
			}

			auto tri_vh = vh_to_trivh_[vhi];
			trimesh_.property(tri_vh_trans_, tri_vh) = idx;
		}

	for(auto ehi : mesh_.edges())
		if(mesh_.is_boundary(ehi))
		{
			auto he0 = mesh_.halfedge_handle(ehi, 0);

			auto vh_f = mesh_.halfedge(he0).from_vertex();
			auto vh_t = mesh_.halfedge(he0).to_vertex();

			//edge valance property
			auto tri_he0 = trimesh_.find_halfedge(vh_to_trivh_[vh_f], vh_to_trivh_[vh_t]);
			auto tri_eh = trimesh_.edge_handle(tri_he0);
			trimesh_.property(tri_val_, tri_eh) = valance_[ehi];
		}

	//singular node property
	for(auto vhi : mesh_.vertices())
		if(mesh_.is_boundary(vhi))
		{
			auto tri_vh = vh_to_trivh_[vhi];
			trimesh_.property(tri_node_type_, tri_vh) =  node_type_[vhi];
		}
}


template <class MeshT1, class MeshT2>
void
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
create_trimesh()
{
	trimesh_.add_property(tet_he_);
	//get surface mesh of tetmesh
	for(auto vhi : mesh_.vertices())
		if(mesh_.is_boundary(vhi))
		{
			auto tri_vh = trimesh_.add_vertex(mesh_.vertex(vhi));
			vh_to_trivh_[vhi] = tri_vh;
		}

	//add trimesh faces
	for(auto hfi : mesh_.halffaces())
		if(mesh_.is_boundary(hfi))
		{
			std::vector<TriVH> tri_vhs;
			for(auto hfv_it = mesh_.hfv_iter(hfi); hfv_it.valid(); ++hfv_it)
				tri_vhs.push_back(vh_to_trivh_[*hfv_it]);

			trimesh_.add_face(tri_vhs);
		}

	for(auto ehi : mesh_.edges())
	{
		if(mesh_.is_boundary(ehi))
		{
			auto he0 = mesh_.halfedge_handle(ehi, 0);
			auto he1 = mesh_.halfedge_handle(ehi, 1);

			auto vh_f = mesh_.halfedge(he0).from_vertex();
			auto vh_t = mesh_.halfedge(he0).to_vertex();

			//store tet halfedges
			auto tri_he0 = trimesh_.find_halfedge(vh_to_trivh_[vh_f], vh_to_trivh_[vh_t]);
			trimesh_.property(tet_he_, tri_he0) = he0;
			trimesh_.property(tet_he_, trimesh_.opposite_halfedge_handle(tri_he0)) = he1;
		}
	}
}


template <class MeshT1, class MeshT2>
void
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
constrained_chart_mergingT()
{
	VP<bool> unconnected_vh_prop =
			mesh_.template request_vertex_property<bool>("unconnected_vh");
	mesh_.set_persistent(unconnected_vh_prop, true);
	//two edges of the vertex are not connected
	for(auto ehi : mesh_.edges())
		if(!mesh_.is_boundary(ehi) && valance_[ehi] != 0)
		{
			unconnected_vh_prop[mesh_.edge(ehi).from_vertex()] = true;
			unconnected_vh_prop[mesh_.edge(ehi).to_vertex()] = true;
		}

	//for all tangent continuity constraints
	for(int i=1; i<=max_label_; ++i)
	{
		std::vector<EH> ehs = sg_.get_singular_edges_of_label(i);
		if(mesh_.is_boundary(ehs[0]))
			continue;
		for(unsigned int j=0; j<ehs.size()-1; ++j)
		{
			tangent_continuity(ehs[j], ehs[j+1]);

			VH vh_cm = common_vertex_handle(ehs[j], ehs[j+1]);
			if(vh_cm.is_valid())
				unconnected_vh_prop[vh_cm] = false;
		}
	}
}


template <class MeshT1, class MeshT2>
void
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
merge_volume_charts()
{
	constrained_chart_mergingTB();

	//add corner constraints
	std::map<EH, std::vector<int> > eh_to_corners;
	add_corner_constraints(eh_to_corners);
	add_tangent_continuity_constraints();

	//connect independent sg components
	if(max_comp_i_ > 1)
	{
		std::cerr<<"\n#######Connect interior graph components";
		std::cerr<<"\nmax comp_i: "<<max_comp_i_;
		connect_interior_graph_components();

		//check
		if(!check_solved_regular_edges_transitions())
			return;
	}

	//higher genus, and sg is not closed
	add_constraint_to_avoid_twist();

	std::cerr<<"\nMerge as many charts as possible...";
	int tags[4] = {1, 2, -1, -2};
	for(int i=0; i<4; ++i)
		dual_spanning_forest(tags[i]);
	zip_edges();

//	check_solved_regular_edges_transitions();

//	StopWatch<> time;
//	time.start();
	std::cerr<<"\nProcess locally decidable corners and patches...";
	local_decidable_chart_merging(eh_to_corners);

//	std::cerr << "\ntime:    " <<time.stop()/1000<<" s";
}


template <class MeshT1, class MeshT2>
void
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
constrained_chart_mergingTB()
{
	VP<bool> unconnected_vh_prop =
			mesh_.template request_vertex_property<bool>("unconnected_vh");

	//for all interior singular tubers that are incident to the boundary, merge to boundary charts
	for(auto vhi : mesh_.vertices())
		if(mesh_.is_boundary(vhi) && unconnected_vh_prop[vhi])
		{
			HEH he_sg(-1);
			for(auto voh_it = mesh_.voh_iter(vhi); voh_it.valid(); ++voh_it)
			{
				auto eho = mesh_.edge_handle(*voh_it);
				if(!mesh_.is_boundary(eho) && valance_[eho] != 0)
				{
					he_sg = *voh_it;
					break;
				}
			}
			auto ch_isg = *mesh_.hec_iter(he_sg);
			int tube_chart_id = real_chart_id_[chart_id_[ch_isg]];
			//if already connected
			if(tube_chart_id < 0)
				continue;

			HEH he_bnd(-1);
			for(auto voh_it = mesh_.voh_iter(vhi); voh_it.valid(); ++voh_it)
				if(mesh_.is_boundary(*voh_it))
				{
					he_bnd = *voh_it;
					break;
				}
			auto ch_bnd = *mesh_.hec_iter(he_bnd);
			int bc_id = real_chart_id_[chart_id_[ch_bnd]];

			//find dual path to boundary shell
			std::vector<HFH> d_path;
			get_path_from_fan_to_boundary_restricting_to_the_vertex(mesh_.edge_handle(he_sg), d_path);

			int axis = axes_prop_[*mesh_.hec_iter(he_sg)][he_sg] == 0 ? 0 : 1;

			if(axis == 1)
			{
				for(unsigned int i=0; i<d_path.size(); ++i)
				{
					trans_prop_[d_path[i]] = 0;
					trans_prop_[mesh_.opposite_halfface_handle(d_path[i])] = 0;
				}
			}else
			{
				trans_prop_[d_path[0]] = 3;
				trans_prop_[mesh_.opposite_halfface_handle(d_path[0])] = 3;
				for(unsigned int i=1; i<d_path.size(); ++i)
				{
					trans_prop_[d_path[i]] = 0;
					trans_prop_[mesh_.opposite_halfface_handle(d_path[i])] = 0;
				}
			}

			//set cell properties
			for(unsigned int i=0; i<d_path.size(); ++i)
			{
				auto ch_path = mesh_.incident_cell(d_path[i]);
				cell_tag_[ch_path] = 1;
				chart_id_[ch_path] = chart_id_[ch_isg];
			}

			update_chart_index_mapping(tube_chart_id, bc_id);
			unconnected_vh_prop[vhi] = false;
		}

	mesh_.set_persistent(unconnected_vh_prop, false);
}


template <class MeshT1, class MeshT2>
void
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
interior_singular_edge_closing(const EH _eh)
{
	HEH he_dir;
	auto he = mesh_.halfedge_handle(_eh, 0);
	auto ch = *mesh_.hec_iter(he);
	if(axes_prop_[ch][he] != -1)
		he_dir = he;
	else
		he_dir = mesh_.opposite_halfedge_handle(he);

	//assign transition
	auto val = valance_[_eh];
	auto hehf_it = mesh_.hehf_iter(he);
	if(trans_prop_[*hehf_it] != -1)
		std::cerr<<"\nError: edge closing is wrong!";
	if(val == -1)
	{
		if(he_dir == he)
		{
			trans_prop_[*hehf_it] = 4;
			trans_prop_[mesh_.opposite_halfface_handle(*hehf_it)] = 5;
		}else
		{
			trans_prop_[*hehf_it] = 5;
			trans_prop_[mesh_.opposite_halfface_handle(*hehf_it)] = 4;
		}
	}else if(val == 1)
	{
		if(he_dir == he)
		{
			trans_prop_[*hehf_it] = 5;
			trans_prop_[mesh_.opposite_halfface_handle(*hehf_it)] = 4;
		}else
		{
			trans_prop_[*hehf_it] = 4;
			trans_prop_[mesh_.opposite_halfface_handle(*hehf_it)] = 5;
		}
	}


	hehf_it++;
	for(;hehf_it.valid(); ++hehf_it)
	{
		if(trans_prop_[*hehf_it] != -1)
			std::cerr<<"\nError: edge closing is wrong!";
		trans_prop_[*hehf_it] = 0;
		trans_prop_[mesh_.opposite_halfface_handle(*hehf_it)] = 0;
	}
}


template <class MeshT1, class MeshT2>
void
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
tangent_continuity(const EH _eh, const EH _eh_next)
{
	auto he0 = mesh_.halfedge_handle(_eh, 0);
	auto ch = *mesh_.hec_iter(he0);
	if(axes_prop_[ch][he0] == -1)
		he0 = mesh_.opposite_halfedge_handle(he0);
	auto he1 = mesh_.halfedge_handle(_eh_next, 0);
	ch = *mesh_.hec_iter(he1);
	if(axes_prop_[ch][he1] == -1)
		he1 = mesh_.opposite_halfedge_handle(he1);

	std::vector<HFH> d_path;
	auto ch_s = *mesh_.hec_iter(he0);

	shortest_dual_path_from_cell_to_edge_restricting_to_the_vertex(ch_s, _eh_next, d_path);

	if(are_halfedges_the_same_direction(he0, he1))
	{
		for(unsigned int i=0; i<d_path.size(); ++i)
		{
			trans_prop_[d_path[i]] = 0;
			trans_prop_[mesh_.opposite_halfface_handle(d_path[i])] = 0;
		}
	}else
	{
		trans_prop_[d_path[0]] = 3;
		trans_prop_[mesh_.opposite_halfface_handle(d_path[0])] = 3;
		for(unsigned int i=1; i<d_path.size(); ++i)
		{
			trans_prop_[d_path[i]] = 0;
			trans_prop_[mesh_.opposite_halfface_handle(d_path[i])] = 0;
		}
	}

	//mark visited cells
	for(unsigned int i=0; i<d_path.size(); ++i)
	{
		chart_id_[mesh_.incident_cell(d_path[i])] = label_[_eh];
		cell_tag_[mesh_.incident_cell(d_path[i])] = arc_type_[_eh];
	}
}


template <class MeshT1, class MeshT2>
void
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
get_path_from_fan_to_boundary_restricting_to_the_vertex(const EH _eh, std::vector<HFH>& _hfs)const
{
	//get the common vertex
	auto vh_f = mesh_.edge(_eh).from_vertex();
	auto vh_t = mesh_.edge(_eh).to_vertex();

	auto vh_bnd = mesh_.is_boundary(vh_f) ? vh_f : vh_t;

	std::map<CH, bool> cell_free, bnd_shell_cells;
	//previous halfface on the path
	std::map<CH, HFH> pre_dp;
	for(auto vc_it = mesh_.vc_iter(vh_bnd); vc_it.valid(); ++vc_it)
	{
		if(cell_tag_[*vc_it] == 1)
			bnd_shell_cells[*vc_it] = true;
		cell_free[*vc_it] = true;
		pre_dp[*vc_it] = HFH(-1);
	}
	for(auto hec_it = mesh_.hec_iter(mesh_.halfedge_handle(_eh, 0)); hec_it.valid(); ++hec_it)
	{
		bnd_shell_cells[*hec_it] = false;
		cell_free[*hec_it] = false;
	}

	//search for dual path
	auto ch_s = *mesh_.hec_iter(mesh_.halfedge_handle(_eh, 0));
	std::queue<CH> que;
	que.push(ch_s);
	CH ch_e(-1);
	bool found = false;
	while(!que.empty() && !found)
	{
		auto ch = que.front();
		que.pop();

		auto hfs = mesh_.cell(ch).halffaces();
		for(int ii=0; ii<4; ++ii)
		{
			HFH hf_opp = mesh_.opposite_halfface_handle(hfs[ii]);

			if(mesh_.is_boundary(hf_opp))
				continue;
			auto ch_opp = mesh_.incident_cell(hf_opp);

			if(bnd_shell_cells[ch_opp])
			{
				found = true;
				ch_e = ch_opp;
				pre_dp[ch_opp] = hf_opp;
				break;
			}

			if(!cell_free[ch_opp])
				continue;
			else
			{
				cell_free[ch_opp] = false;
				pre_dp[ch_opp] = hf_opp;

				que.push(ch_opp);
			}
		}
	}


	//find the path
	if(found)
	{
		std::vector<HFH> rhfs;
		std::queue<CH> que_dp;
		que_dp.push(ch_e);

		while(!que_dp.empty())
		{
			CH ch = que_dp.front();
			que_dp.pop();

			HFH hf = pre_dp[ch];

			if(hf != HFH(-1))
			{
				rhfs.push_back(hf);
				que_dp.push(mesh_.incident_cell(mesh_.opposite_halfface_handle(hf)));
			}
		}

		for(auto riter = rhfs.rbegin(); riter != rhfs.rend(); ++riter)
			_hfs.push_back(*riter);
	}else
		std::cerr<<"\nError: cannot find path from boundary to singular edge "<<_eh;
}


template <class MeshT1, class MeshT2>
void
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
shortest_dual_path_from_cell_to_edge_restricting_to_the_vertex(const CH _ch, const EH _eh, std::vector<HFH>& _hfs)const
{
	//get the common vertex
	auto cell_vhs = mesh_.get_cell_vertices(_ch);
	std::vector<VH> edge_vhs;
	edge_vhs.push_back(mesh_.edge(_eh).from_vertex());
	edge_vhs.push_back(mesh_.edge(_eh).to_vertex());

	std::sort(cell_vhs.begin(), cell_vhs.end());
	std::sort(edge_vhs.begin(), edge_vhs.end());

	std::vector<VH> common_vh;
	std::set_intersection(cell_vhs.begin(), cell_vhs.end(), edge_vhs.begin(), edge_vhs.end(), std::back_inserter(common_vh));

	if(common_vh.empty())
	{
		std::cerr<<"\nError: cannot find common vertex!";
		return;
	}

	std::map<CH, bool> cell_free;
	//previous halfface on the path
	std::map<CH, HFH> pre_dp;
	for(auto vc_it = mesh_.vc_iter(common_vh[0]); vc_it.valid(); ++vc_it)
	{
		if(cell_tag_[*vc_it] == 0)
			cell_free[*vc_it] = true;
		pre_dp[*vc_it] = HFH(-1);
	}

	_hfs.clear();

	std::queue<CH> que;
	que.push(_ch);
	cell_free[_ch] = false;
	CH ch_e(-1);
	bool found = false;
	while(!que.empty() && !found)
	{
		auto ch = que.front();
		que.pop();

		auto hfs = mesh_.cell(ch).halffaces();
		for(int ii=0; ii<4; ++ii)
		{
			auto hf_opp = mesh_.opposite_halfface_handle(hfs[ii]);

			if(mesh_.is_boundary(hf_opp))
				continue;
			auto ch_opp = mesh_.incident_cell(hf_opp);

			if(cell_has_edge(ch_opp, _eh))
			{
				found = true;
				ch_e = ch_opp;
				pre_dp[ch_opp] = hf_opp;
				break;
			}

			if(cell_tag_[ch_opp] == 0 && cell_free[ch_opp])
			{
				cell_free[ch_opp] = false;
				pre_dp[ch_opp] = hf_opp;

				que.push(ch_opp);
			}
		}
	}


	//find the path
	if(found)
	{
		std::vector<HFH> rhfs;
		std::queue<CH> que_dp;
		que_dp.push(ch_e);

		while(!que_dp.empty())
		{
			auto ch = que_dp.front();
			que_dp.pop();

			auto hf = pre_dp[ch];

			if(hf != HFH(-1))
			{
				rhfs.push_back(hf);
				que_dp.push(mesh_.incident_cell(mesh_.opposite_halfface_handle(hf)));
			}
		}

		for(auto riter = rhfs.rbegin(); riter != rhfs.rend(); ++riter)
		{
			_hfs.push_back(*riter);
//			bool good_path = false;
//			auto hf_vhs = mesh_.get_halfface_vertices(*riter);
//			for(auto j : hf_vhs)
//				if(j == common_vh[0])
//				{
//					good_path = true;
//					break;
//				}
//			if(!good_path)
//				std::cerr<<"\nWarning: searched dual path might be invalid!";
		}
	}else
		std::cerr<<"\nError: cannot find path from cell "<<_ch<<" to singular edge "<<_eh;
}


template <class MeshT1, class MeshT2>
bool
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
cell_has_edge(const CH _ch, const EH _eh)const
{
	std::vector<VH> vhs = mesh_.get_cell_vertices(_ch);
	auto vhe0 = mesh_.edge(_eh).from_vertex();
	auto vhe1 = mesh_.edge(_eh).to_vertex();
	bool has_vhe0 = false;
	for(int i=0; i<4; ++i)
		if(vhs[i] == vhe0)
			has_vhe0 = true;

	if(!has_vhe0)
		return false;
	else
	{
		for(int i=0; i<4; ++i)
			if(vhs[i] == vhe1)
				return true;
	}

	return false;
}

template <class MeshT1, class MeshT2>
VH
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
common_vertex_handle(const EH _eh0, const EH _eh1)const
{
	auto vhe00 = mesh_.edge(_eh0).from_vertex();
	auto vhe01 = mesh_.edge(_eh0).to_vertex();

	auto vhe10 = mesh_.edge(_eh1).from_vertex();
	auto vhe11 = mesh_.edge(_eh1).to_vertex();

	if(vhe10 == vhe00)
		return vhe00;
	if(vhe10 == vhe01)
		return vhe01;
	if(vhe11 == vhe00)
		return vhe00;
	if(vhe11 == vhe01)
		return vhe01;

	return VH(-1);
}


template <class MeshT1, class MeshT2>
HEH
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
common_halfedge_handle(const HFH _hf0, const HFH _hf1)const
{
	auto hes0 = mesh_.halfface(_hf0).halfedges();
	auto hes1 = mesh_.halfface(_hf1).halfedges();
	std::sort(hes0.begin(), hes0.end());
	std::sort(hes1.begin(), hes1.end());

	std::vector<HEH> common_he;
	std::set_intersection(hes0.begin(), hes0.end(), hes1.begin(), hes1.end(), std::back_inserter(common_he));

	if(common_he.size() != 1)
	{
		std::cerr<<"\nError: couldn't find common halfedge!";
		return HEH(-1);
	}
	return common_he[0];
}


template <class MeshT1, class MeshT2>
void
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
cell_edges(const CH _ch, std::vector<EH>& _ehs)const
{
	std::vector<VH> vhs = mesh_.get_cell_vertices(_ch);

	_ehs.push_back(mesh_.edge_handle(mesh_.halfedge(vhs[0], vhs[1])));
	_ehs.push_back(mesh_.edge_handle(mesh_.halfedge(vhs[0], vhs[2])));
	_ehs.push_back(mesh_.edge_handle(mesh_.halfedge(vhs[0], vhs[3])));
	_ehs.push_back(mesh_.edge_handle(mesh_.halfedge(vhs[1], vhs[2])));
	_ehs.push_back(mesh_.edge_handle(mesh_.halfedge(vhs[1], vhs[3])));
	_ehs.push_back(mesh_.edge_handle(mesh_.halfedge(vhs[2], vhs[3])));
}


template <class MeshT1, class MeshT2>
void
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
local_decidable_chart_merging(std::map<EH, std::vector<int> >& _eh_to_corners)
{
	std::queue<std::shared_ptr<CornerT<TetMesh> > >que_cnr;
	bool need_initial_connect = true;
	for(auto i : v_corners_)
		if(i.get()->connectable(mesh_, trans_prop_, chart_id_, real_chart_id_))
		{
			que_cnr.push(i);
			need_initial_connect = false;
			break;
		}

	//if there's no locally decidable corner, e.g. sphere, any solution is globally valid
	if(need_initial_connect)
	{
		std::cerr<<"\nInitial connect...";
		HEH he_connected(-1);
		for(auto i : v_corners_)
			if(i.get()->initial_connectable(mesh_, trans_prop_, chart_id_, real_chart_id_))
			{
				he_connected = i.get()->initial_constrained_chart_mergingC(mesh_, 0, tq_, trans_prop_, chart_id_, real_chart_id_);
				break;
			}

		zip_edges();

		std::vector<EH> ehs;
		if(he_connected.is_valid())
			get_end_edges_of_an_arc(he_connected, ehs);

		for(auto ehi : ehs)
			for(auto i : _eh_to_corners[ehi])
			{
				que_cnr.push(v_corners_[i]);
			}
	}

	//
	do
	{
		//iteratively process local decidable corners
		while(!que_cnr.empty())
		{
			auto cnr_cur = que_cnr.front();
			que_cnr.pop();

			if(cnr_cur.get()->connectable(mesh_, trans_prop_, chart_id_, real_chart_id_))
			{
				auto he_connected = cnr_cur.get()->constrained_chart_mergingC(mesh_, 0, tq_, trans_prop_, chart_id_, real_chart_id_);

//				std::cerr<<"\nconnect edge: "<<mesh_.edge_handle(he_connected)<<" "
//						<<mesh_.halfedge(he_connected).from_vertex()<<" "<<mesh_.halfedge(he_connected).to_vertex();
				zip_edges();

	//			check_solved_regular_edges_transitions();

	//			std::vector<EH> ehs;
	//			if(he_connected.is_valid())
	//				get_end_edges_of_an_arc(he_connected, ehs);
	//			else
	//				std::cerr<<"\nError: connected edge is invalid!";
	//
	//			//push neighbouring corners
	//			for(auto ehi : ehs)
	//				for(auto i : _eh_to_corners[ehi])
	//				{
	//					que_cnr.push(v_corners_[i]);
	//					std::cerr<<" "<<i;
	//				}
				for(auto i : v_corners_)
				{
					if(i.get()->n_connected(trans_prop_) <=2)
						que_cnr.push(i);
				}
			}
		}

		//process locally decidable patch
		constrained_chart_mergingP();
		//in case neither local decidable patch nor local decidable corner exists, pick one valid solution
		apply_exhaustive_search();

		for(auto i : v_corners_)
		{
			if(i.get()->connectable(mesh_, trans_prop_, chart_id_, real_chart_id_))
				que_cnr.push(i);
		}
	}while(!que_cnr.empty());

	//tangent continuity of corners
	for(auto i : v_tc_)
	{
		i.get()->constrained_chart_mergingTC(mesh_,  trans_prop_, chart_id_, real_chart_id_);
		zip_edges();
	}
}


template <class MeshT1, class MeshT2>
void
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
get_end_edges_of_an_arc(const HEH _he_first, std::vector<EH>& _ehs)const
{
	auto eh_first = mesh_.edge_handle(_he_first);
	auto label_first = label_[eh_first];
	_ehs.push_back(eh_first);
	auto node_first = mesh_.halfedge(_he_first).from_vertex();

	auto vhs = sg_.get_singular_vertices_of_label(label_first);
	auto vh_second = (node_first == vhs[0]) ? vhs[vhs.size()-1] : vhs[0];

	if(!mesh_.is_boundary(vh_second))
		for(auto voh_it = mesh_.voh_iter(vh_second); voh_it.valid(); ++voh_it)
		{
			if(label_[mesh_.edge_handle(*voh_it)] == label_[eh_first])
			{
				//push the singular edge of the same singular arc
				_ehs.push_back(mesh_.edge_handle(*voh_it));
				break;
			}
		}
//	std::cerr<<"\ncorresponding edges: ";
//	for(auto ehi : _ehs)
//		std::cerr<<" "<<ehi;
}


template <class MeshT1, class MeshT2>
void
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
update_interior_graph_component_property()
{
	for(auto ehi : mesh_.edges())
		if(sg_comp_i_[ehi] == 2)
			sg_comp_i_[ehi] = 1;

	for(auto ehi : mesh_.edges())
		if(sg_comp_i_[ehi] == max_comp_i_)
			sg_comp_i_[ehi] = 2;

	max_comp_i_--;
}


template <class MeshT1, class MeshT2>
void
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
update_boundary_graph_component_property()
{
	for(auto ehi : mesh_.edges())
		if(sg_comp_b_[ehi] == 2)
			sg_comp_b_[ehi] = 1;

	for(auto ehi : mesh_.edges())
		if(sg_comp_b_[ehi] == max_comp_b_)
			sg_comp_b_[ehi] = 2;

	max_comp_b_--;
}


template <class MeshT1, class MeshT2>
void
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
update_chart_index_mapping(const int _id0, const int _id1)
{
	int id_old = std::abs(_id0) < std::abs(_id1) ? _id1 : _id0;
	int id_new = std::abs(_id0) < std::abs(_id1) ? _id0 : _id1;

	for(auto it = real_chart_id_.begin(); it != real_chart_id_.end(); ++it)
		if(it->second == id_old)
			it->second = id_new;
}


template <class MeshT1, class MeshT2>
void
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
constrained_chart_mergingP()
{
	//find patches
	HFP<int> patch_id = mesh_.template request_halfface_property<int>("patch_id");
	EP<bool> edge_block = mesh_.template request_edge_property<bool>("is_edge_blocked");
	std::vector<std::vector<HFH> > patches;
	get_patches(patch_id, edge_block, patches);

	//search on all patches
	for(auto patch : patches)
	{
		//if the patch is already detemined
		if(trans_prop_[patch[0]] != -1)
			continue;

		std::vector<VH> patch_sg_vhs;
		get_singular_vertices_on_patch(patch, patch_sg_vhs);

		//if no singular vertex on the patch, continue searching next patch
		if(patch_sg_vhs.empty())
			continue;

		std::vector<HFH> hfs;
		std::vector<std::vector<int> > v_trans;
		get_first_constraints_on_patch(patch_id, patch_sg_vhs, hfs, v_trans);

		//find the other orthogonal constraint
		HFH hf_breach(-1);
		int unique_trans = 0;
		bool solved = get_second_constraints_on_patch(patch_id, edge_block, hfs, v_trans, hf_breach, unique_trans);
		if(!solved)
			continue;

		//assign transtion
		trans_prop_[hf_breach] = unique_trans;
		trans_prop_[mesh_.opposite_halfface_handle(hf_breach)] = tq_.inverse_transition_idx(unique_trans);

		zip_edges();

		//update chart id
		int chart0 = real_chart_id_[chart_id_[mesh_.incident_cell(hf_breach)]];
		int chart1 = real_chart_id_[chart_id_[mesh_.incident_cell(mesh_.opposite_halfface_handle(hf_breach))]];
		update_chart_index_mapping(chart0, chart1);
	}
}


template <class MeshT1, class MeshT2>
void
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
get_patches(HFP<int>& _patch_id, EP<bool>& _edge_block, std::vector<std::vector<HFH> >& _patches)const
{
	FP<bool> uk_face = mesh_.template request_face_property<bool>("is_unknown_face");
	for(auto fhi : mesh_.faces())
		if(!mesh_.is_boundary(fhi))
		{
			auto hf0 = mesh_.halfface_handle(fhi, 0);
			if(trans_prop_[hf0] == -1)
				uk_face[fhi] = true;
		}

	//mark the boundary of the patch
	for(auto fhi : mesh_.faces())
	{
		if(uk_face[fhi])
		{
			auto hf0 = mesh_.halfface_handle(fhi, 0);
			auto hes = mesh_.halfface(hf0).halfedges();
			for(auto hei : hes)
				if(!mesh_.is_boundary(hei))
				{
					int count = 0;
					for(auto hehf_it = mesh_.hehf_iter(hei); hehf_it.valid(); ++hehf_it)
					{
						if(trans_prop_[*hehf_it] == -1)
							count++;
					}
					if(count != 2)
						_edge_block[mesh_.edge_handle(hei)] = true;
				}
		}
	}

//	ACG::HaltonColors hcol;
//	for(auto fhi : mesh_.faces())
//		color_[fhi] = ACG::Vec4f(1,1,1,0);

	int n_patch = 0;
	for(auto hfi : mesh_.halffaces())
	{
		auto fhi = mesh_.face_handle(hfi);
		auto hfi_opp = mesh_.opposite_halfface_handle(hfi);
		if(uk_face[fhi] && !_patch_id[hfi] && !_patch_id[hfi_opp])
		{
			n_patch++;
			_patch_id[hfi] = n_patch;
//			ACG::Vec4f cur_col = hcol.get_next_color();
//			color_[fhi] = cur_col;

			std::vector<HFH> patch;
			std::queue<HFH> que;
			patch.push_back(hfi);
			que.push(hfi);
			while(!que.empty())
			{
				auto hf_cur = que.front();
				auto fh_cur = mesh_.face_handle(hf_cur);
				que.pop();

				auto hes = mesh_.halfface(hf_cur).halfedges();
				for(int i=0; i<3; ++i)
				{
					auto he_opp = mesh_.opposite_halfedge_handle(hes[i]);
					if(mesh_.is_boundary(he_opp))
						continue;

					if(_edge_block[mesh_.edge_handle(he_opp)])
						continue;

					HFH hf_next(-1);
					for(auto hehf_it = mesh_.hehf_iter(he_opp); hehf_it.valid(); ++hehf_it)
						if(uk_face[mesh_.face_handle(*hehf_it)] && (mesh_.face_handle(*hehf_it) != fh_cur))
							hf_next = *hehf_it;
					if(!hf_next.is_valid())
						continue;

					if(!_patch_id[hf_next])
					{
						_patch_id[hf_next]  = n_patch;
//						color_[mesh_.face_handle(hf_next)] = cur_col;

						que.push(hf_next);
						patch.push_back(hf_next);
					}
				}
			}

			_patches.push_back(patch);
		}
	}
//	std::cerr<<"\n"<<n_patch<<" patches";
}


template <class MeshT1, class MeshT2>
void
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
get_singular_vertices_on_patch(const std::vector<HFH>& _patch, std::vector<VH>& _patch_sg_vhs)const
{
	std::vector<bool> is_vt_visited(mesh_.n_vertices(), false);
	for(auto hfi : _patch)
	{
		for(auto hfv : mesh_.halfface_vertices(hfi))
		{
			if(is_vt_visited[hfv.idx()])
				continue;
			if((is_singular_vt_[hfv] && mesh_.is_boundary(hfv))
					|| node_type_[hfv] > 0)
			{
				is_vt_visited[hfv.idx()] = true;
				_patch_sg_vhs.push_back(hfv);
			}
		}
	}
}


template <class MeshT1, class MeshT2>
void
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
get_first_constraints_on_patch(const HFP<int>& _patch_id, const std::vector<VH>& _patch_sg_vhs,
		std::vector<HFH>& _hfs, std::vector<std::vector<int> >& _v_trans)const
{
	for(auto vhi : _patch_sg_vhs)
	{
		//if it is singular node. TODO: deal with tangent continuity at singular nodes
		if(node_type_[vhi] > 0)
		{
			for(auto i : v_corners_)
			{
				if(i.get()->node(mesh_) == vhi)
				{
					if(i.get()->n_connected(trans_prop_) != 2)
						continue;

					std::vector<std::vector<HFH> > dpaths(3);
					int pos = 0;
					HFH hf_uk(-1);
					for(int j=0; j<3; ++j)
						i.get()->dual_path(j, dpaths[j]);
					for(int j=0; j<3; ++j)
					{
						HFH hf_last = dpaths[j][dpaths[j].size()-1];
						if((_patch_id[mesh_.opposite_halfface_handle(hf_last)] || _patch_id[hf_last])
								&& trans_prop_[hf_last] == -1)
						{
							hf_uk = hf_last;
							pos = j;
							break;
						}
					}

					if(hf_uk.is_valid())
					{
						int axis0=-1, axis1=-1, axis2=-1;
						if(pos == 0)
						{
							axis0 = axis_after_transition(0, trans_prop_[dpaths[1][dpaths[1].size()-1]]);
							axis1 = axis_after_transition(0, trans_prop_[dpaths[2][dpaths[2].size()-1]]);
							axis2 = the_third_axis(axis0, axis1);
						}else if(pos == 1)
						{
							axis0 = axis_after_transition(0, trans_prop_[dpaths[2][dpaths[2].size()-1]]);
							axis1 = axis_after_transition(0, trans_prop_[dpaths[0][dpaths[0].size()-1]]);
							axis2 = the_third_axis(axis0, axis1);
						}else if(pos == 2)
						{
							axis0 = axis_after_transition(0, trans_prop_[dpaths[0][dpaths[0].size()-1]]);
							axis1 = axis_after_transition(0, trans_prop_[dpaths[1][dpaths[1].size()-1]]);
							axis2 = the_third_axis(axis0, axis1);
						}

						std::vector<int> trans;
						for(int ii=0; ii<24; ++ii)
							if(axis_after_transition(0, ii) == axis2)
								trans.push_back(ii);

						if(_patch_id[hf_uk])
						{
							_hfs.push_back(hf_uk);
						}else if(_patch_id[mesh_.opposite_halfface_handle(hf_uk)])
						{
							_hfs.push_back(mesh_.opposite_halfface_handle(hf_uk));
							for(int ii=0; ii<4; ++ii)
								trans[ii] = tq_.inverse_transition_idx(trans[ii]);
						}
						_v_trans.push_back(trans);
					}
				}
			}
		}

		//if the singular vertex is at the end of a singular arc and on the boundary
		if(mesh_.is_boundary(vhi))
		{
			HEH sg_he(-1);
			for(auto voh_it = mesh_.voh_iter(vhi); voh_it.valid(); ++voh_it)
			{
				EH eho = mesh_.edge_handle(*voh_it);
				if(!mesh_.is_boundary(eho) && valance_[eho] != 0)
				{
					sg_he = *voh_it;
					break;
				}
			}

			auto vh_t = mesh_.halfedge(sg_he).to_vertex();
			//tell if the axis direction is the same as the normal direction
			bool same_dir = true;
			HEH he0 = mesh_.halfedge_handle(mesh_.edge_handle(sg_he), 0);
			if(he0 == sg_he)
				same_dir = false;
			else
				same_dir = true;
			//in case distance between interior singular node and boundary is 1
			if(node_type_[vh_t] > 0)
				same_dir = true;
			int axis = same_dir ? 0 : 1;

			//find dual path to boundary shell
			std::vector<HFH> dpath;
			get_path_from_fan_to_boundary_restricting_to_the_vertex(mesh_.edge_handle(sg_he), dpath);

			HFH hf_tmp(-1);
			int n_unknown = 0;
			for(auto hfi : dpath)
				if(trans_prop_[hfi] == -1)
				{
					hf_tmp = hfi;
					n_unknown++;
				}
			HFH hf_p(-1);
			if((_patch_id[mesh_.opposite_halfface_handle(hf_tmp)] || _patch_id[hf_tmp])
					&& n_unknown == 1)
				hf_p = hf_tmp;

			if(!hf_p.is_valid())
				continue;

			int idx0 = 0, idx1 = 0;
			bool hf_found = false;
			for(auto hfi : dpath)
			{
				if(hfi == hf_p)
				{
					hf_found = true;
					continue;
				}
				if(!hf_found)
					idx0 = tq_.mult_transitions_idx(trans_prop_[hfi], idx0);
				else
					idx1 = tq_.mult_transitions_idx(trans_prop_[hfi], idx1);
			}

			std::vector<int> trans;
			for(int ii=0; ii<24; ++ii)
				if(axis_after_transition(axis, ii) == 0)
				{
					int idx_tmp = tq_.mult_transitions_idx(tq_.inverse_transition_idx(idx1), tq_.mult_transitions_idx(ii, tq_.inverse_transition_idx(idx0)));
					trans.push_back(idx_tmp);
				}

			if(_patch_id[mesh_.opposite_halfface_handle(hf_p)])
			{
				for(int ii=0; ii<4; ++ii)
					trans[ii] = tq_.inverse_transition_idx(trans[ii]);

				hf_p = mesh_.opposite_halfface_handle(hf_p);
			}
			if(hf_p.is_valid())
			{
				_hfs.push_back(hf_p);
				_v_trans.push_back(trans);
			}
		}
	}
}


template <class MeshT1, class MeshT2>
bool
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
get_second_constraints_on_patch(const HFP<int>& _patch_id, const EP<bool>& _edge_block, const std::vector<HFH>& _hfs,
    		const std::vector<std::vector<int> >& _v_trans, HFH& _hf_breach, int& _unique_trans)const
{
	for(size_t i=0; i<_hfs.size(); ++i)
	{
		HFH hf_tb(-1);
		std::vector<HFH> dpath2;
		shortest_dual_path_from_boundary_to_halfface_on_patch(_hfs[i], hf_tb, _patch_id, _edge_block, dpath2);
		//if it is on a invalid patch, continue
		if(!hf_tb.is_valid())
			continue;
		HEH he_bnd(-1);
		for(int ii=0; ii<(int)mesh_.halfface(hf_tb).halfedges().size(); ++ii)
			if(mesh_.is_boundary(mesh_.halfface(hf_tb).halfedges()[ii]))
			{
				he_bnd = mesh_.halfface(hf_tb).halfedges()[ii];
				break;
			}

		//combine two orthogonal constraints and get unique transition
		//express the axis in the chart incident to hf_breach
		int idx = 0;

		HFH hf_bdy(-1);
		for(auto hehf_it = mesh_.hehf_iter(he_bnd); hehf_it.valid(); hehf_it++)
			if(mesh_.is_boundary(*hehf_it))
				hf_bdy = *hehf_it;
		std::vector<HFH> hfs2;
		halffaces_between_two_halffaces(mesh_.opposite_halfedge_handle(he_bnd),
				mesh_.opposite_halfface_handle(hf_bdy), mesh_.opposite_halfface_handle(hf_tb), hfs2);

		for(auto hfi : hfs2)
		{
			if(trans_prop_[hfi] == -1)
				std::cerr<<"\nError: illegal dual path!"<<mesh_.face_handle(hfi);
			idx = tq_.mult_transitions_idx(trans_prop_[hfi], idx);
		}
		for(int k=0; k<(int)dpath2.size()-1; ++k)
		{
			HEH hei = common_halfedge_handle(dpath2[k], mesh_.opposite_halfface_handle(dpath2[k+1]));

			std::vector<HFH> hfs;
			halffaces_between_two_halffaces(hei, dpath2[k], mesh_.opposite_halfface_handle(dpath2[k+1]), hfs);
			for(auto hfi : hfs)
			{
				if(trans_prop_[hfi] == -1)
					std::cerr<<"\nError: illegal dual path!"<<mesh_.face_handle(hfi);
				idx = tq_.mult_transitions_idx(trans_prop_[hfi], idx);
			}
		}
		int axis1 = axis_after_transition(0, idx);

		//express the other axis in the chart incident to hf_breach
		idx = 0;
		HFH hf_bdy2(-1);
		for(auto hehf_it = mesh_.hehf_iter(mesh_.opposite_halfedge_handle(he_bnd)); hehf_it.valid(); hehf_it++)
			if(mesh_.is_boundary(*hehf_it))
				hf_bdy2 = *hehf_it;

		halffaces_between_two_halffaces(he_bnd, mesh_.opposite_halfface_handle(hf_bdy2), hf_tb, hfs2);
		for(auto hfi : hfs2)
		{
			if(trans_prop_[hfi] == -1)
				std::cerr<<"\nError: illegal dual path!"<<mesh_.face_handle(hfi);
			idx = tq_.mult_transitions_idx(trans_prop_[hfi], idx);
		}

		for(int k=0; k<(int)dpath2.size()-1; ++k)
		{
			HEH hei = common_halfedge_handle(mesh_.opposite_halfface_handle(dpath2[k]), dpath2[k+1]);

			std::vector<HFH> hfs;
			halffaces_between_two_halffaces(hei, mesh_.opposite_halfface_handle(dpath2[k]), dpath2[k+1], hfs);
			for(auto hfi : hfs)
			{
				if(trans_prop_[hfi] == -1)
					std::cerr<<"\nError: illegal dual path!"<<mesh_.face_handle(hfi);
				idx = tq_.mult_transitions_idx(trans_prop_[hfi], idx);
			}
		}
		int axis0 = axis_after_transition(0, idx);

		//try to get unique transition
		std::vector<int> unique_trans;
		for(int ii=0; ii<4; ++ii)
			if(axis_after_transition(axis0, _v_trans[i][ii]) == axis1)
				unique_trans.push_back(_v_trans[i][ii]);

		if(unique_trans.size() == 1)
		{
			_hf_breach = _hfs[i];
			_unique_trans = unique_trans[0];
			std::cerr<<"\nDetermine locally decidable patch!";
			return true;
		}else
			continue;
	}

	return false;
}


template <class MeshT1, class MeshT2>
void
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
halffaces_between_two_halffaces(const HEH _he, const HFH _hf0, const HFH _hf1, std::vector<HFH>& _hfs)const
{
	_hfs.clear();

	HFH hf_it = _hf0;
	while(hf_it != _hf1)
	{
		hf_it = mesh_.adjacent_halfface_in_cell(hf_it, _he);
		hf_it = mesh_.opposite_halfface_handle(hf_it);

		if(hf_it != _hf1)
			_hfs.push_back(hf_it);
	}
}


template <class MeshT1, class MeshT2>
void
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
shortest_dual_path_from_boundary_to_halfface_on_patch(const HFH _hf,  HFH& _hf_tb, const HFP<int>& _patch_id,
   		const EP<bool>& _edge_block, std::vector<HFH>& _hfs)const
{
	if(!_hf.is_valid())
	{
		std::cerr<<"\nError: input halfface is invalid!";
		return;
	}
	std::vector<bool> hf_visited(mesh_.n_halffaces(), false);
	std::vector<HFH> pre_hf(mesh_.n_halffaces(), HFH(-1));

	std::queue<HFH> que;
	que.push(_hf);
	hf_visited[_hf.idx()] = true;
	int patch_id = _patch_id[_hf];
	bool hf_found = false;

	while(!que.empty() && !hf_found)
	{
		auto hf_cur = que.front();
		que.pop();

		auto hes = mesh_.halfface(hf_cur).halfedges();
		for(int i=0; i<3; ++i)
		{
			auto he_opp = mesh_.opposite_halfedge_handle(hes[i]);
			if(mesh_.is_boundary(he_opp))
				continue;

			if(_edge_block[mesh_.edge_handle(he_opp)])
				continue;

			HFH hf_next(-1);
			for(auto hehf_it = mesh_.hehf_iter(he_opp); hehf_it.valid(); ++hehf_it)
				if(_patch_id[*hehf_it] == patch_id)
					hf_next = *hehf_it;
			if(!hf_next.is_valid())
				continue;

			if(!hf_visited[hf_next.idx()])
			{
				hf_visited[hf_next.idx()] = true;
				pre_hf[hf_next.idx()] = hf_cur;

				for(int ii=0; ii<(int)mesh_.halfface(hf_next).halfedges().size(); ++ii)
					if(mesh_.is_boundary(mesh_.halfface(hf_next).halfedges()[ii]))
					{
						hf_found = true;
						_hf_tb = hf_next;
						break;
					}

				que.push(hf_next);
			}
		}
	}

	if(!_hf_tb.is_valid())
	{
//		std::cerr<<"\nWarning: no dual path on this patch leads to boundary!";
		return;
	}

	//get the path
	std::vector<HFH> path_dpath;
	std::queue<HFH> que_dpath;
	que_dpath.push(_hf_tb);
	while(!que_dpath.empty())
	{
		auto hf = que_dpath.front();
		que_dpath.pop();
		if(hf == HFH(-1))
			break;

		_hfs.push_back(hf);

		que_dpath.push(pre_hf[hf.idx()]);
	}
}


template <class MeshT1, class MeshT2>
void
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
apply_exhaustive_search()
{
	if(open_edges_.empty())
	{
		std::cerr<<"\nNo need to do exhaustive search!";
		return;
	}

	//get unknown faces
	std::vector<FH> uk_fhs;
	for(auto fhi : mesh_.faces())
		if(!mesh_.is_boundary(fhi))
		{
			auto hf0 = mesh_.halfface_handle(fhi, 0);
			if(trans_prop_[hf0] == -1)
				uk_fhs.push_back(fhi);
		}


	if(uk_fhs.empty())
	{
		std::cerr<<"\nNo need to do exhaustive search!";
		return;
	}

	std::vector<HFH> key_hfs;
	std::vector<std::vector<int> > v_key_trans;
	get_keys(key_hfs, v_key_trans);

	for(auto fhi : uk_fhs)
	{
		auto hf0 = mesh_.halfface_handle(fhi, 0);
		auto hf1 = mesh_.halfface_handle(fhi, 1);

		trans_prop_[hf0] = -1;
		trans_prop_[hf1] = -1;
	}

	exhaustive_search(uk_fhs, key_hfs, v_key_trans);
}


template <class MeshT1, class MeshT2>
void
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
get_keys(std::vector<HFH>& _hfs, std::vector<std::vector<int> >& _v_trans)
{
	EP<bool> ioe_prop = mesh_.template request_edge_property<bool>("is_open_edge");

	update_open_edge_property_on_boundary(ioe_prop);


	//there's only one unknown transition in the half fans of the open edges
	//which means the key halffaces are independent
	for(auto ehi : open_edges_)
	{
		if(ioe_prop[ehi])
		{
			//check
			auto hei = mesh_.halfedge_handle(ehi, 0);
			int n_unknown = 0;
			for(auto hehf_it = mesh_.hehf_iter(hei); hehf_it.valid(); ++hehf_it)
			{
				auto fh = mesh_.face_handle(*hehf_it);
				if(!mesh_.is_boundary(fh) && trans_prop_[*hehf_it] == -1)
					n_unknown++;
			}
			if(n_unknown > 1)
			{
				std::cerr<<"\nWarning: more than one unknown transitions at the open edge: "<<ehi;
				continue;
			}
			if(n_unknown == 0)
			{
				std::cerr<<"\nError: something was wrong in updating open edge: "<<ehi;
				continue;
			}

			//get the four possible transitions that maps axis_0 to axis_0
			std::vector<int> trans;
			for(int i=0; i<24; ++i)
				if(axis_after_transition(0, i) == 0)
					trans.push_back(i);

			int idx_0 = 0,
					idx_1 = 0;
			bool if_found = false;
			HFH hf_uk;
			for(auto hehf_it = mesh_.hehf_iter(hei); hehf_it.valid(); ++hehf_it)
			{
				auto fh = mesh_.face_handle(*hehf_it);
				if(!mesh_.is_boundary(fh))
				{
					if(trans_prop_[*hehf_it] == -1)
					{
						if_found = true;
						hf_uk = *hehf_it;
						continue;
					}
					if(!if_found)
						idx_0 = tq_.mult_transitions_idx(trans_prop_[*hehf_it], idx_0);
					else
						idx_1 = tq_.mult_transitions_idx(trans_prop_[*hehf_it], idx_1);
				}
			}

			//get the transitions of the breach face
			for(int i=0; i<4; ++i)
			{
				trans[i]= tq_.mult_transitions_idx(trans[i], tq_.inverse_transition_idx(idx_0));
				trans[i] = tq_.mult_transitions_idx(tq_.inverse_transition_idx(idx_1), trans[i]);
			}

			if(hf_uk.is_valid())
			{
				_hfs.push_back(hf_uk);
				_v_trans.push_back(trans);
			}
			else
			{
				std::cerr<<"\nError: breach halfface is invalid!";
				return;
			}

			//fake solving
			trans_prop_[hf_uk] = trans[0];
			trans_prop_[mesh_.opposite_halfface_handle(hf_uk)] = tq_.inverse_transition_idx(trans[0]);
			zip_edges(false);

			//update property
			update_open_edge_property_on_boundary(ioe_prop);
		}
	}
}


template <class MeshT1, class MeshT2>
void
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
exhaustive_search(const std::vector<FH>& _uk_fhs, const std::vector<HFH>& _hfs, const std::vector<std::vector<int> >&_v_trans)
{
	if(_hfs.empty())
		return;

	int num = std::pow(4, _hfs.size());
	int count = num-1;
	std::vector<int> trans_indexes(_v_trans.size(), 0);

	while(count >= 0)
	{
		get_indexes(count, trans_indexes);

		for(unsigned int i=0; i<_hfs.size(); ++i)
		{
			trans_prop_[_hfs[i]] = _v_trans[i][trans_indexes[i]];
			trans_prop_[mesh_.opposite_halfface_handle(_hfs[i])] = tq_.inverse_transition_idx(_v_trans[i][trans_indexes[i]]);
		}

		zip_edges(false);
//		std::cerr<<"\nface: "<<mesh_.face_handle(_hfs[0])<<" transiton: "<<trans_prop_[_hfs[0]]<<" count: "<<count;

		if(check_solved_regular_edges_transitions(true))
		{
			//update chart id
			for(auto hfi : _hfs)
			{
				int chart0 = real_chart_id_[chart_id_[mesh_.incident_cell(mesh_.opposite_halfface_handle(hfi))]];
				int chart1 = real_chart_id_[chart_id_[mesh_.incident_cell(hfi)]];
				update_chart_index_mapping(chart0, chart1);
			}

			std::cerr<<"\nExaustive search! One candidate solution found!";
			return;
		}
		for(auto fhi : _uk_fhs)
		{
			auto hf0 = mesh_.halfface_handle(fhi, 0);
			auto hf1 = mesh_.halfface_handle(fhi, 1);

			trans_prop_[hf0] = -1;
			trans_prop_[hf1] = -1;
		}

		count--;
	}

	std::cerr<<"\nError: cannot find a valid solution with ehausive_search!";
}


template <class MeshT1, class MeshT2>
void
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
get_indexes(const int _decimal, std::vector<int>& _indexes)const
{
	if(_decimal == 0)
	{
		_indexes[0] = 0;
		return;
	}
	int dc_num = _decimal;
	int idx_of_indexes = _indexes.size() - 1;
	while(dc_num > 0)
	{
		int mod = dc_num%4;
		dc_num /= 4;
		_indexes[idx_of_indexes] = mod;

		idx_of_indexes--;
	}
}


template <class MeshT1, class MeshT2>
void
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
update_open_edge_property_on_boundary(EP<bool>& _ioe_prop)const
{
	int n_unknown = 0;
	for(auto ehi : open_edges_)
	{
		n_unknown = 0;
		for(auto hehf_it = mesh_.hehf_iter(mesh_.halfedge_handle(ehi, 0)); hehf_it.valid(); ++hehf_it)
		{
			auto fh = mesh_.face_handle(*hehf_it);
			if(!mesh_.is_boundary(fh) && trans_prop_[*hehf_it] == -1)
				n_unknown++;
		}
		if(n_unknown == 0)
			_ioe_prop[ehi] = false;
		else
			_ioe_prop[ehi] = true;
	}
}


template <class MeshT1, class MeshT2>
void
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
dual_spanning_forest(const int _cell_tag)
{
	std::queue<CH> que;
	for(auto chi : mesh_.cells())
		if(cell_tag_[chi] == _cell_tag)
			que.push(chi);
	while(!que.empty())
	{
		auto ch_cur = que.front();
		que.pop();

		auto hfs = mesh_.cell(ch_cur).halffaces();

		for(unsigned int ii=0; ii<4; ++ii)
		{
    			FH fh = mesh_.face_handle(hfs[ii]);
    			if(mesh_.is_boundary(fh))
    				continue;

    			auto hf_opp = mesh_.opposite_halfface_handle(hfs[ii]);
    			auto  neigh_ch  = mesh_.incident_cell(hf_opp);

    			if(cell_tag_[neigh_ch] != 0)
    				continue;

    			trans_prop_[hfs[ii]] = 0;
    			trans_prop_[hf_opp] = 0;
    			cell_tag_[neigh_ch] = _cell_tag;
    			chart_id_[neigh_ch] = chart_id_[ch_cur];
    			que.push(neigh_ch);
		}
	}

}


//iteratively close available edges
template <class MeshT1, class MeshT2>
void
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
zip_edges(const bool _mark_edge)
{
	std::queue<EH> que;
	for(auto ehi : mesh_.edges())
		if(!edge_tag_[ehi])
			que.push(ehi);
	while(!que.empty())
	{
		auto eh = que.front();
		que.pop();

		if(mesh_.is_boundary(eh))
			continue;
		if(edge_tag_[eh])
			continue;

		auto he = mesh_.halfedge_handle(eh, 0);
		HFH hf;
		if(is_removable(he, hf))
		{
			regular_edge_closing(he, hf);
			if(_mark_edge)
				edge_tag_[eh] = true;

			auto v_hes = mesh_.halfface(hf).halfedges();
			for(int ii=0;ii<3;++ii)
			{
				auto ehi = mesh_.edge_handle(v_hes[ii]);
				if(!edge_tag_[ehi])
					que.push(ehi);
			}
		}
	}
}


template <class MeshT1, class MeshT2>
bool
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
is_removable(const HEH _he, HFH &_hf)const
{
	int n_cutfaces = 0;
	for(auto hehf_it = mesh_.hehf_iter(_he); hehf_it.valid();++hehf_it)
		if(trans_prop_[*hehf_it] == -1)
		{
			_hf = *hehf_it;
			n_cutfaces++;
		}

	if(n_cutfaces == 1)
		return true;
	return false;
}


template <class MeshT1, class MeshT2>
void
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
regular_edge_closing(const HEH _he, const HFH _hf)
{
	auto hf_it = mesh_.adjacent_halfface_in_cell(_hf, _he);
	hf_it = mesh_.opposite_halfface_handle(hf_it);
	int idx = 0;
	while(hf_it != _hf)
	{
		idx = tq_.mult_transitions_idx(trans_prop_[hf_it], idx);
		//ccw
		hf_it = mesh_.adjacent_halfface_in_cell(hf_it, _he);
		hf_it = mesh_.opposite_halfface_handle(hf_it);
	}

	trans_prop_[_hf] = tq_.inverse_transition_idx(idx);
	trans_prop_[mesh_.opposite_halfface_handle(_hf)] = idx;
}


template <class MeshT1, class MeshT2>
int
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
axis_after_transition(const int _axis_start, const int _transition)const
{
	if(_transition == -1 || _axis_start > 6 || _axis_start < 0)
		std::cerr<<"\nERROR: input transition or axis is invalid! Axis start: "<<_axis_start<<" transition: "<<_transition;

	auto m1 = tq_.transition_matrix_int(_transition);

	int axes_mapped[3] = {0,0,0};
	for(int j=0; j<3;++j)
		for(int i=0;i<3;++i)
		{
			if(m1(i,j) == 1)
			{
				axes_mapped[j] = 2*i;
				break;
			}
			if(m1(i,j) == -1)
			{
				axes_mapped[j] = 2*i+1;
				break;
			}
		}


	int pos = _axis_start/2;
	int res = _axis_start%2;
	int axis = 0;

	if(res == 0)
		axis = axes_mapped[pos];
	else if(res == 1)
			axis = (axes_mapped[pos] %2 == 0) ? axes_mapped[pos] + 1 : axes_mapped[pos] - 1;

	return axis;
}


//get the third axis in the right hand coordinate system
template <class MeshT1, class MeshT2>
int
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
the_third_axis(const int _axis0, const int _axis1) const
{
	int pos0 = _axis0/2;
	int pos1 = _axis1/2;
	int pos2 = 0;
	if(pos0+pos1 == 1)
		pos2 = 2;
	else if(pos0+pos1 == 2)
		pos2 =1;
	else if(pos0+pos1 == 3)
		pos2 = 0;

	int axis2 = 2*pos2;

	int prod = (_axis1- _axis0)*(axis2 - _axis1)*(_axis0 - axis2);
	if(prod > 0)
		prod = 1;
	else
		prod = 0;
	prod = (_axis0%2 + _axis1%2 + prod)%2;

	if(prod == 1)
		axis2 = axis2+1;

	return axis2;
}


template <class MeshT1, class MeshT2>
bool
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
check_solved_regular_edges_transitions(const bool _mute)
{
	bool smile = true;
	for(auto ehi : mesh_.edges())
		if(valance_[ehi] == 0)
		{
			auto he0 = mesh_.halfedge_handle(ehi, 0);
			if(!edge_tag_[ehi])
			{
				if(!mesh_.is_boundary(ehi))
				{
					int idx = 0;
					for(auto hehf_it =  mesh_.hehf_iter(he0);hehf_it.valid(); ++hehf_it)
					{
						if(trans_prop_[*hehf_it] == -1)
						{
							idx = 0;
							break;
						}
						idx = tq_.mult_transitions_idx(trans_prop_[*hehf_it], idx);
					}

					if(idx != 0)
					{
						smile = false;
						if(!_mute)
						{
							std::cerr<<"\nError: transition around edge "<<ehi<<" is Wrong! Transition is "<<idx;
							status_[ehi].set_selected(true);
						}
						break;
					}
				}else
				{
					int idx = 0;
					for(auto hehf_it =  mesh_.hehf_iter(he0);hehf_it.valid(); ++hehf_it)
						if(!mesh_.is_boundary(mesh_.face_handle(*hehf_it)))
						{
							if(trans_prop_[*hehf_it] == -1)
							{
								idx = 0;
								break;
							}
							idx = tq_.mult_transitions_idx(trans_prop_[*hehf_it], idx);
						}

					if(idx != 0 && idx != 1 && idx != 4 && idx != 5)
					{
						smile = false;
						if(!_mute)
						{
							std::cerr<<"\nError: transition around edge "<<ehi<<" is Wrong! Transition is "<<idx;
							status_[ehi].set_selected(true);
						}
						break;
					}
				}
			}
		}

	return smile;
}


template <class MeshT1, class MeshT2>
void
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
check_transitions()
{
	for(auto fhi : mesh_.faces())
		if(!mesh_.is_boundary(fhi))
		{
			auto hf0 = mesh_.halfface_handle(fhi, 0);
			if(trans_prop_[hf0] == -1)
			{
				std::cerr<<"\nError: solving is not complete!";
				return;
			}
		}
	for(auto ehi : mesh_.edges())
	{
		if(!mesh_.is_boundary(ehi))
		{
			int idx = 0;
			auto he0 = mesh_.halfedge_handle(ehi, 0);
			auto he1 = mesh_.halfedge_handle(ehi, 1);
			auto hehf_it =  mesh_.hehf_iter(he0);
			for(;hehf_it.valid(); ++hehf_it)
				idx = tq_.mult_transitions_idx(trans_prop_[*hehf_it], idx);

			HFH hf_start = *hehf_it;
			HFH hf_start_opp = mesh_.opposite_halfface_handle(hf_start);

			CH ch = mesh_.incident_cell(hf_start_opp);
			if(!mesh_.incident_cell(hf_start_opp).is_valid())
			{
				std::cerr<<"\nError: invalid cell!";
				return;
			}

			if(idx != 0)
			{
				int axis = -1;
				if(axes_prop_[ch][he0] == -1 && axes_prop_[ch][he1] == -1)
				{
					valance_[ehi] = 4;
					std::cerr<<"\nError: wrong holonomy on interior regular edge! Eh: "<<ehi<<" idx: "<<idx;
				}
				else if(axes_prop_[ch][he0] != -1 && axes_prop_[ch][he1] == -1)
					axis = axes_prop_[ch][he0];
				else if(axes_prop_[ch][he0] == -1 && axes_prop_[ch][he1] != -1)
					axis = (axes_prop_[ch][he1]%2 == 0) ? axes_prop_[ch][he1]+1 : axes_prop_[ch][he1]-1;
				else
				{
					valance_[ehi] = 4;
					std::cerr<<"\n\nError: two axis directions assigned to the edge, which should not happen! Eh: "
							<<ehi<<" "<<axes_prop_[ch][he0]<<" " <<axes_prop_[ch][he1];
				}

				if(axis != -1)
				{
					if(axis == (idx-4))
						valance_[ehi] = -1;
					else if(std::abs(axis - (idx-4)) == 1)
						valance_[ehi] = 1;
					else
					{
						std::cerr<<"\nError: axis does not agree with the rotational axis! Eh: "<<ehi<<", axis: "<<axis<<", rotational axis: "<<idx;
						valance_[ehi] = 4;
					}
				}
			}else
				valance_[ehi] = 0;

		}else
		{
			int idx = 0;
			auto he0 = mesh_.halfedge_handle(ehi, 0);
			auto he1 = mesh_.halfedge_handle(ehi, 1);
			auto hf_start = HFH(-1),
					hf_end = HFH(-1);
			//counterclockwise
			for(auto hehf_it =  mesh_.hehf_iter(he0);hehf_it.valid(); ++hehf_it)
			{
				if(!mesh_.is_boundary(mesh_.face_handle(*hehf_it)))
					idx = tq_.mult_transitions_idx(trans_prop_[*hehf_it], idx);

				//store hf_end
				hf_end = *hehf_it;
			}

			hf_start = mesh_.opposite_halfface_handle(*mesh_.hehf_iter(he0));
			int n0 = normal_prop_[hf_start];
			int n0_t = axis_after_transition(n0, idx);

			int n1 = normal_prop_[hf_end];

			if(n0_t == n1)
				valance_[ehi] = 0;
			else if((n0_t/2) == (n1/2))
				valance_[ehi] = 2;
			else
			{
				auto ch_end = mesh_.incident_cell(mesh_.opposite_halfface_handle(hf_end));
				if(!ch_end.is_valid())
				{
					std::cerr<<"\nError: invalid cell ch_end, OVM did not sort the boundary halfedge halffaces! Eh: "
							<<ehi<<" "<<mesh_.halfedge(he0).from_vertex()<<" "<<mesh_.halfedge(he0).to_vertex();

					for(auto hehf_it =  mesh_.hehf_iter(he0);hehf_it.valid(); ++hehf_it)
						std::cerr<<" "<<mesh_.face_handle(*hehf_it);
					return;
				}
				int axis = axes_prop_[ch_end][he0];
				if(axis == -1)
				{
					axis = axes_prop_[ch_end][he1];
					if(axis == -1)
					{
						std::cerr<<"\nError: wrong holonomy on boundary regular edge! Eh: "<<mesh_.edge_handle(he0)<<" "<<idx<<" "<<n0_t;
						valance_[ehi] = 4;
					}
					else
						axis = (axis%2) ? axis-1 : axis+1;
				}

				int axis_rt = the_third_axis(n0_t, n1);
				if(axis_rt == axis)
					valance_[ehi] = 1;
				else if(axis_rt/2 == axis/2)
					valance_[ehi] = -1;
			}
		}
	}
}


template <class MeshT1, class MeshT2>
bool
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
check_input()const
{
	if(!mesh_.template edge_property_exists<int>("edge_valance"))
	{
		std::cerr<<"\nError: no edge_valance property found!";
		return false;
	}
	if(!mesh_.template mesh_property_exists<std::vector<VH> >("corner_vertices") ||
	   !mesh_.template mesh_property_exists<std::vector<std::vector<HFH> > >("corner_dual_paths"))
	{
		std::cerr<<"Warning: no corner constraints found!";
	}
	if(!mesh_.template mesh_property_exists<std::vector<VH> >("sector_vertices") ||
	   !mesh_.template mesh_property_exists<std::vector<double> >("sector_angles"))
	{
		std::cerr<<"Warning: no sector constraints found!";
	}

	if(!check_input_singularity_graph())
		return false;

	return true;
};


template <class MeshT1, class MeshT2>
bool
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
check_input_singularity_graph()const
{
	for(auto chi : mesh_.cells())
	{
		std::vector<EH> ehs;
		cell_edges(chi, ehs);

		int n_sgl = 0;
		EH sg_eh(-1);
		for(int i=0; i<6; ++i)
			if(valance_[ehs[i]] != 0)
			{
				sg_eh = ehs[i];
				n_sgl++;
			}
		if(n_sgl > 1)
		{
			std::cerr<<"\nError: two singular edges in one cell!"<<chi;
			return false;
		}else if(n_sgl == 1 && !mesh_.is_boundary(sg_eh))
		{
			auto hfs = mesh_.cell(chi).halffaces();
			for(auto hfi : hfs)
				if(mesh_.is_boundary(mesh_.face_handle(hfi)))
				{
					std::cerr<<"\nError: singular cell has boundary face!"<<chi;
					return false;
				}
		}
	}

	std::vector<EH> ehs2;
	for(auto ehi : mesh_.edges())
		if(valance_[ehi] == 0)
		{
			HEH he = mesh_.halfedge_handle(ehi, 0);
			VH vh0 = mesh_.halfedge(he).from_vertex();
			VH vh1 = mesh_.halfedge(he).to_vertex();
			if(is_singular_vt_[vh1] && is_singular_vt_[vh0])
			{
				std::cerr<<"\nError: regular edge is incident to two singular tubes!"<<ehi;
				return false;
			}
		}

	return true;
}


template <class MeshT1, class MeshT2>
void
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
connect_interior_graph_components()
{
	if(max_comp_i_ <= 1)
		return;

	//if there's no comp touching boundary, connect first component(circle) to boundary
	bool is_all_circle = true;
	for(auto ehi : mesh_.edges())
		if(sg_comp_i_[ehi] == 1)
		{
			if(arc_type_[ehi] != 2)
				is_all_circle = false;

			break;
		}

	if(is_all_circle)
	{
		//choose a cell
		EH eh_1;
		for(auto ehi : mesh_.edges())
			if(sg_comp_i_[ehi] == 1 && status_[ehi].selected())
			{
				eh_1 = ehi;
				break;
			}
		if(!eh_1.is_valid())
		{
			std::cerr<<"\nError: no additional constraints is found! Solve may get stuck or the solution might not be optimum!";
			return;
		}
		auto ch_1 = *mesh_.hec_iter(mesh_.halfedge_handle(eh_1, 0));

		std::vector<HFH> dpath;
		shortest_dual_path_from_cell_to_boundary(ch_1, dpath);

		connect_via_dual_path(dpath, 2, 0);
	}

//	//add constraints in component one(fertility_rsf)
//	std::vector<EH> ehs2;
//	for(auto ehi : mesh_.edges())
//		if(sg_comp_i_[ehi] == 1 && status_[ehi].selected())
//			ehs2.push_back(ehi);
//
//	if(ehs2.size() == 2)
//	{
//		auto ch_2 = *mesh_.hec_iter(mesh_.halfedge_handle(ehs2[0], 0));
//		auto he2 = (axes_prop_[ch_2][mesh_.halfedge_handle(ehs2[0], 0)] == 0) ?
//				mesh_.halfedge_handle(ehs2[0], 0) : mesh_.halfedge_handle(ehs2[0], 1);
//
//		auto ch_1 = *mesh_.hec_iter(mesh_.halfedge_handle(ehs2[1], 0));
//		auto he1 = (axes_prop_[ch_1][mesh_.halfedge_handle(ehs2[1], 0)] == 0) ?
//				mesh_.halfedge_handle(ehs2[1], 0) : mesh_.halfedge_handle(ehs2[1], 1);
//
//		std::cerr<<"\nch1: "<<ch_1<<" ch_2: "<<ch_2;
//		std::cerr<<"\nhe1: "<<mesh_.halfedge(he1).from_vertex()<<" "<<mesh_.halfedge(he1).to_vertex()
//				<<" he2: "<<mesh_.halfedge(he2).from_vertex()<<" "<<mesh_.halfedge(he2).to_vertex()
//				<<" "<<axes_prop_[ch_1][he1]<<" "<<axes_prop_[ch_2][he2];
//
//		std::vector<HFH> dpath2;
//		shortest_dual_path_between_cells(ch_2, ch_1, dpath2);
//
//		int axis_2 = 0;
//		auto vec0 = mesh_.vertex(mesh_.halfedge(he1).to_vertex()) - mesh_.vertex(mesh_.halfedge(he1).from_vertex());
//		vec0.normalize();
//		auto vec1 = mesh_.vertex(mesh_.halfedge(he2).to_vertex()) - mesh_.vertex(mesh_.halfedge(he2).from_vertex());
//		vec1.normalize();
//
//		if((vec0|vec1) < 0.0)
//			axis_2 = 1;
//
//		connect_via_dual_path(dpath2, 0, axis_2);
//	}


	//connect other components
	//choose a cell
	while(max_comp_i_ > 1)
	{
		std::vector<EH> ehs2;
		for(auto ehi : mesh_.edges())
			if(sg_comp_i_[ehi] == 2 && status_[ehi].selected())
				ehs2.push_back(ehi);

		if(ehs2.empty())
		{
			std::cerr<<"\nError: no additional constraints is found! Solve may get stuck or the solution might not be optimum!";
			update_interior_graph_component_property();
			continue;
		}

		for(auto eh_2 : ehs2)
		{
			auto ch_2 = *mesh_.hec_iter(mesh_.halfedge_handle(eh_2, 0));
			auto he2 = (axes_prop_[ch_2][mesh_.halfedge_handle(eh_2, 0)] == 0) ?
					mesh_.halfedge_handle(eh_2, 0) : mesh_.halfedge_handle(eh_2, 1);

			EH eh_1;
			find_closest_parallel_singular_edges(eh_2, 1, eh_1);
			status_[eh_1].set_selected(true);
			auto ch_1 = *mesh_.hec_iter(mesh_.halfedge_handle(eh_1, 0));
			auto he1 = (axes_prop_[ch_1][mesh_.halfedge_handle(eh_1, 0)] == 0) ?
					mesh_.halfedge_handle(eh_1, 0) : mesh_.halfedge_handle(eh_1, 1);

//			std::cerr<<"\nch1: "<<ch_1<<" ch_2: "<<ch_2;
//			std::cerr<<"\nhe1: "<<mesh_.halfedge(he1).from_vertex()<<" "<<mesh_.halfedge(he1).to_vertex()
//					<<" he2: "<<mesh_.halfedge(he2).from_vertex()<<" "<<mesh_.halfedge(he2).to_vertex()
//					<<" "<<axes_prop_[ch_1][he1]<<" "<<axes_prop_[ch_2][he2];

			std::vector<HFH> dpath2;
			shortest_dual_path_between_cells(ch_2, ch_1, dpath2);

			int axis_2 = 0;
			auto vec0 = mesh_.vertex(mesh_.halfedge(he1).to_vertex()) - mesh_.vertex(mesh_.halfedge(he1).from_vertex());
			vec0.normalize();
			auto vec1 = mesh_.vertex(mesh_.halfedge(he2).to_vertex()) - mesh_.vertex(mesh_.halfedge(he2).from_vertex());
			vec1.normalize();

			if((vec0|vec1) < 0.0)
				axis_2 = 1;

			connect_via_dual_path(dpath2, 0, axis_2);
		}

		//update component property
		update_interior_graph_component_property();
		check_solved_regular_edges_transitions();
	}

	std::cerr<<"\nMax comp_i: "<<max_comp_i_;
}


template <class MeshT1, class MeshT2>
void
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
connect_via_dual_path(const std::vector<HFH>& _dpath, int _axis1, int _axis2)
{
	//find first unknown halfface
	int trans0 = 0, trans1 = 0;
	int pos = -1;
	find_independent_unknown_halfface_transition(_dpath, pos);
	if(pos == -1)
		std::cerr<<"\nError: couldn't find independent unknown halfface!";
//	std::cerr<<"\nfh_uk: "<<mesh_.face_handle(_dpath[pos])<<" trans hf_uk: "<<trans_prop_[_dpath[pos]];

	for(unsigned int i=0; i<_dpath.size(); ++i)
	{
		if((int)i < pos)
		{
			if(trans_prop_[_dpath[i]] == -1)
			{
				trans_prop_[_dpath[i]] = 0;
				trans_prop_[mesh_.opposite_halfface_handle(_dpath[i])] = 0;

				//check if there's edge can be closed
				zip_edges();
			}else
				trans0 = tq_.mult_transitions_idx(trans_prop_[_dpath[i]], trans0);
		}else if((int)i > pos)
		{
			if(trans_prop_[_dpath[i]] == -1)
			{
				trans_prop_[_dpath[i]] = 0;
				trans_prop_[mesh_.opposite_halfface_handle(_dpath[i])] = 0;

				//check if there's edge can be closed
				zip_edges();

			}else
				trans1 = tq_.mult_transitions_idx(trans_prop_[_dpath[i]], trans1);
		}
	}

	//choose one transition for the unknown
	int trans = 0;
	for(int i=0; i<24; ++i)
		if(axis_after_transition(_axis2, i) == _axis1)//x map to +x or -x
		{
			trans = i;
			break;
		}

//	std::cerr<<"\ntransition: "<<trans<<" fh_uk: "<<mesh_.face_handle(_dpath[pos])<<" trans hf_uk: "<<trans_prop_[_dpath[pos]];
	int trans_uk = tq_.mult_transitions_idx(tq_.inverse_transition_idx(trans1), tq_.mult_transitions_idx(trans ,tq_.inverse_transition_idx(trans0)));
	trans_prop_[_dpath[pos]] = trans_uk;
	trans_prop_[mesh_.opposite_halfface_handle(_dpath[pos])] = tq_.inverse_transition_idx(trans_uk);

	zip_edges();


	//update chart id
	for(auto hfi : _dpath)
		if(chart_id_[mesh_.incident_cell(hfi)] == 0)
			chart_id_[mesh_.incident_cell(hfi)] = chart_id_[mesh_.incident_cell(_dpath[_dpath.size()-1])];
	int chart0 = real_chart_id_[chart_id_[mesh_.incident_cell(mesh_.opposite_halfface_handle(_dpath[0]))]];
	int chart1 = real_chart_id_[chart_id_[mesh_.incident_cell(_dpath[_dpath.size()-1])]];
	update_chart_index_mapping(chart0, chart1);
	//update cell tag
	int cell_tag = cell_tag_[mesh_.incident_cell(_dpath[_dpath.size()-1])];
	for(auto hfi : _dpath)
		cell_tag_[mesh_.incident_cell(hfi)] = cell_tag;
}


template <class MeshT1, class MeshT2>
void
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
find_closest_parallel_singular_edges(const EH _eh0, const int _comp, EH &_eh1)const
{
	auto vh_s = mesh_.edge(_eh0).from_vertex();
	//find the other
	double min_length = DBL_MAX;
	if(!mesh_.is_boundary(_eh0))
	{
		for(auto ehi : mesh_.edges())
			if(sg_comp_i_[ehi] == _comp)
			{
				auto vh_e = mesh_.edge(ehi).from_vertex();
				auto dist = distance(vh_s, vh_e);
				if(dist < min_length)
				{
					min_length = dist;
					_eh1 = ehi;
				}
			}
	}else
	{
		for(auto ehi : mesh_.edges())
			if(sg_comp_b_[ehi] == _comp)
			{
				auto vh_e = mesh_.edge(ehi).from_vertex();
				auto dist = distance(vh_s, vh_e);
				if(dist < min_length)
				{
					min_length = dist;
					_eh1 = ehi;
				}
			}
	}


	auto vec0 = mesh_.vertex(mesh_.edge(_eh0).to_vertex()) - mesh_.vertex(vh_s);
	vec0.normalize();
	auto vec1 = mesh_.vertex(mesh_.edge(_eh1).to_vertex()) - mesh_.vertex(mesh_.edge(_eh1).from_vertex());
	vec1.normalize();

	if(std::abs(vec0|vec1) < 0.707)
		std::cerr<<"\nWarning: two edges are not parallel!";
}


template <class MeshT1, class MeshT2>
void
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
find_independent_unknown_halfface_transition(const std::vector<HFH>& _dpath, int& _pos)
{
	HFP<int> trans_prop_t = mesh_.template request_halfface_property<int>("trans_property_temp");
	for(auto hfi : mesh_.halffaces())
		trans_prop_t[hfi] = trans_prop_[hfi];
//	std::cerr<<"\n"<<mesh_.face_handle(_dpath[0])<<" "<<mesh_.face_handle(_dpath[1])
//			<<" "<< trans_prop_[_dpath[1]]<<" "<<trans_prop_t[_dpath[1]]<<"\n";

	for(unsigned int i=0; i<_dpath.size(); ++i)
		if(trans_prop_[_dpath[i]] == -1)
		{
			bool found = true;

			trans_prop_[_dpath[i]] = 0;
			trans_prop_[mesh_.opposite_halfface_handle(_dpath[i])] = 0;
			zip_edges(false);
			for(unsigned int j=0; j<_dpath.size(); ++j)
			{
//				std::cerr<<" "<<mesh_.face_handle(_dpath[j])<<" "<<mesh_.face_handle(_dpath[i])
//						 <<" "<< trans_prop_[_dpath[j]]<<" "<<trans_prop_t[_dpath[j]];

				if(j != i && trans_prop_[_dpath[j]] != trans_prop_t[_dpath[j]])
				{
					found = false;
					break;
				}
			}

			for(auto hfi : mesh_.halffaces())
				trans_prop_[hfi] = trans_prop_t[hfi];

			if(found)
			{
				_pos = i;
				return;
			}
		}
}


template <class MeshT1, class MeshT2>
void
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
shortest_dual_path_from_cell_to_boundary(const CH _ch, std::vector<HFH>& _dpath)const
{
	std::vector<bool> cell_visited(mesh_.n_cells(), false);
	std::vector<HFH> pre_hf(mesh_.n_cells(), HFH(-1));

	std::queue<CH> que;
	que.push(_ch);
	cell_visited[_ch.idx()] = true;
	bool reach_boundary = false;
	CH ch_end(-1);
	while(!que.empty() && !reach_boundary)
	{
		auto ch_cur = que.front();
		que.pop();

		auto hfs = mesh_.cell(ch_cur).halffaces();
		for(int i=0; i<4; ++i)
		{
			auto hf_opp = mesh_.opposite_halfface_handle(hfs[i]);
			if(mesh_.is_boundary(hf_opp))
			{
				reach_boundary = true;
				ch_end = mesh_.incident_cell(hfs[i]);
				break;
			}
			auto ch_opp = mesh_.incident_cell(hf_opp);
			if(!cell_visited[ch_opp.idx()])
			{
				cell_visited[ch_opp.idx()] = true;
				pre_hf[ch_opp.idx()] = hfs[i];

				que.push(ch_opp);
			}
		}
	}

	//get the path
	std::vector<HFH> path_dpath;
	std::queue<HFH> que_dpath;
	que_dpath.push(pre_hf[ch_end.idx()]);
	while(!que_dpath.empty())
	{
		auto hf = que_dpath.front();
		que_dpath.pop();
		if(hf == HFH(-1))
			break;

		path_dpath.push_back(mesh_.opposite_halfface_handle(hf));
		que_dpath.push(pre_hf[mesh_.incident_cell(hf).idx()]);
	}

	for(auto r_it = path_dpath.rbegin(); r_it != path_dpath.rend(); ++r_it)
		_dpath.push_back(*r_it);
}


template <class MeshT1, class MeshT2>
void
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
shortest_dual_path_between_cells(const CH _ch0, const CH _ch1, std::vector<HFH>& _dpath)const
{
	std::vector<bool> cell_visited(mesh_.n_cells(), false);
	std::vector<HFH> pre_hf(mesh_.n_cells(), HFH(-1));

	std::queue<CH> que;
	que.push(_ch0);
	cell_visited[_ch0.idx()] = true;
	bool cell_found = false;

	while(!que.empty() && !cell_found)
	{
		auto ch_cur = que.front();
		que.pop();

		auto hfs = mesh_.cell(ch_cur).halffaces();
		for(int i=0; i<4; ++i)
		{
			auto hf_opp = mesh_.opposite_halfface_handle(hfs[i]);
			if(mesh_.is_boundary(hf_opp))
				continue;
			auto ch_opp = mesh_.incident_cell(hf_opp);
			if(!cell_visited[ch_opp.idx()])
			{
				cell_visited[ch_opp.idx()] = true;
				pre_hf[ch_opp.idx()] = hfs[i];

				if(ch_opp == _ch1)
				{
					cell_found = true;
					break;
				}

				que.push(ch_opp);
			}
		}
	}

	//get the path
	std::vector<HFH> path_dpath;
	std::queue<HFH> que_dpath;
	que_dpath.push(pre_hf[_ch1.idx()]);
	while(!que_dpath.empty())
	{
		auto hf = que_dpath.front();
		que_dpath.pop();
		if(hf == HFH(-1))
			break;

		path_dpath.push_back(mesh_.opposite_halfface_handle(hf));
		que_dpath.push(pre_hf[mesh_.incident_cell(hf).idx()]);
	}

	for(auto r_it = path_dpath.rbegin(); r_it != path_dpath.rend(); ++r_it)
		_dpath.push_back(*r_it);
}




template <class MeshT1, class MeshT2>
void
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
add_constraint_to_avoid_twist()
{
	std::vector<HEH> hes;
	for(auto vhi : mesh_.vertices())
	{
		if(status_[vhi].selected())
			for(auto voh_it = mesh_.voh_iter(vhi); voh_it.valid(); ++voh_it)
			{
				EH eho = mesh_.edge_handle(*voh_it);
				if(!mesh_.is_boundary(eho) && valance_[eho] != 0)
				{
					hes.push_back(*voh_it);
					break;
				}
			}
	}

	if(hes.size() == 2)
	{

		auto ch_2 = *mesh_.hec_iter(hes[0]);
		auto he2 = (axes_prop_[ch_2][hes[0]] == 0) ?
				hes[0] : mesh_.opposite_halfedge_handle(hes[0]);

		auto ch_1 = *mesh_.hec_iter(hes[1]);
		auto he1 = (axes_prop_[ch_1][hes[1]] == 0) ?
				hes[1] : mesh_.opposite_halfedge_handle(hes[1]);

//		std::cerr<<"\nch1: "<<ch_1<<" ch_2: "<<ch_2;
//		std::cerr<<"\nhe1: "<<mesh_.halfedge(he1).from_vertex()<<" "<<mesh_.halfedge(he1).to_vertex()
//				<<" he2: "<<mesh_.halfedge(he2).from_vertex()<<" "<<mesh_.halfedge(he2).to_vertex()
//				<<" "<<axes_prop_[ch_1][he1]<<" "<<axes_prop_[ch_2][he2];

		std::vector<HFH> dpath2;
		shortest_dual_path_between_cells(ch_2, ch_1, dpath2);

		int axis_2 = 0;
		auto vec0 = mesh_.vertex(mesh_.halfedge(he1).to_vertex()) - mesh_.vertex(mesh_.halfedge(he1).from_vertex());
		vec0.normalize();
		auto vec1 = mesh_.vertex(mesh_.halfedge(he2).to_vertex()) - mesh_.vertex(mesh_.halfedge(he2).from_vertex());
		vec1.normalize();

		if((vec0|vec1) < 0.0)
			axis_2 = 1;

		connect_via_dual_path(dpath2, 0, axis_2);
	}
}


template <class MeshT1, class MeshT2>
void
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
connect_boundary_graph_components(std::vector<int>& _chart_index)
{
	//connect other components
	//choose a cell
	while(max_comp_b_ > 1)
	{
		std::vector<EH> ehs2;
		for(auto ehi : mesh_.edges())
			if(sg_comp_b_[ehi] == 2 && status_[ehi].selected())
				ehs2.push_back(ehi);

		if(ehs2.empty())
		{
			std::cerr<<"\nWarning: no additional constraints is found! Solve may get stuck or the solution might not be optimum!";
			update_boundary_graph_component_property();

			continue;
		}


		for(auto eh_2 : ehs2)
		{
			auto ch_2 = *mesh_.hec_iter(mesh_.halfedge_handle(eh_2, 0));
			auto he2 = (axes_prop_[ch_2][mesh_.halfedge_handle(eh_2, 0)] == 2) ?
				mesh_.halfedge_handle(eh_2, 0) : mesh_.halfedge_handle(eh_2, 1);

			auto he2_tri = trimesh_.find_halfedge(vh_to_trivh_[mesh_.halfedge(he2).from_vertex()],
					vh_to_trivh_[mesh_.halfedge(he2).to_vertex()]);

			EH eh_1;
			find_closest_parallel_singular_edges(eh_2, 1, eh_1);
			status_[eh_1].set_selected(true);
			auto ch_1 = *mesh_.hec_iter(mesh_.halfedge_handle(eh_1, 0));
			auto he1 = (axes_prop_[ch_1][mesh_.halfedge_handle(eh_1, 0)] == 2) ?
				mesh_.halfedge_handle(eh_1, 0) : mesh_.halfedge_handle(eh_1, 1);

			auto he1_tri = trimesh_.find_halfedge(vh_to_trivh_[mesh_.halfedge(he1).from_vertex()],
					vh_to_trivh_[mesh_.halfedge(he1).to_vertex()]);

//			std::cerr<<"\nch1: "<<ch_1<<" ch_2: "<<ch_2;
//			std::cerr<<"\nhe1: "<<mesh_.halfedge(he1).from_vertex()<<" "<<mesh_.halfedge(he1).to_vertex()
//						<<" he2: "<<mesh_.halfedge(he2).from_vertex()<<" "<<mesh_.halfedge(he2).to_vertex()
//						<<" "<<axes_prop_[ch_1][he1]<<" "<<axes_prop_[ch_2][he2];

			std::vector<OpenMesh::HalfedgeHandle> dpath2;
			shortest_dual_path_between_faces_tri(trimesh_, trimesh_.face_handle(he2_tri), trimesh_.face_handle(he1_tri), dpath2);

			std::vector<int> v_trans;
			for(int i=0; i<24; ++i)
				if(axis_after_transition(0, i) == 0)//x map to x
					v_trans.push_back(i);

			//get transition
			int axis_2 = 2;
			auto vec0 = mesh_.vertex(mesh_.halfedge(he1).to_vertex()) - mesh_.vertex(mesh_.halfedge(he1).from_vertex());
			vec0.normalize();
			auto vec1 = mesh_.vertex(mesh_.halfedge(he2).to_vertex()) - mesh_.vertex(mesh_.halfedge(he2).from_vertex());
			vec1.normalize();

			if((vec0|vec1) < 0.0)
				axis_2 = 3;
			int trans = 0;
			for(auto i : v_trans)
				if(axis_after_transition(axis_2, i) == 2)//y map to +y or -y
				{
					trans = i;
					break;
				}

			bool breach_found = false;
			int count=0;
			TriHEH he_uk(-1);
			int trans0 = 0, trans1 = 0;
			for(auto hei : dpath2)
			{
				if(trimesh_.property(tri_val_, trimesh_.edge_handle(hei)) != 0)
					continue;
				if(trimesh_.property(tri_trans_, hei) == -1)
				{
					breach_found = true;
					he_uk = hei;
					count++;
					continue;
				}
				if(!breach_found)
				{
					trans0 = tq_.mult_transitions_idx(trimesh_.property(tri_trans_, hei), trans0);
				}else
					trans1 = tq_.mult_transitions_idx(trimesh_.property(tri_trans_, hei), trans1);
			}


			if(!he_uk.is_valid())
			{
				std::cerr<<"\nError: cannot find open edge!";
				return;
			}
			if(count != 1)
				std::cerr<<"\nError: find "<<count<<"open edges!";

			int trans_uk = tq_.mult_transitions_idx(tq_.inverse_transition_idx(trans1), tq_.mult_transitions_idx(trans ,tq_.inverse_transition_idx(trans0)));
			trimesh_.property(tri_trans_,he_uk) = trans_uk;
			trimesh_.property(tri_trans_, trimesh_.opposite_halfedge_handle(he_uk)) = tq_.inverse_transition_idx(trans_uk);

//			std::cerr<<"\ntrans: "<<trans<<" trans0: "<<trans0<<" trans1: "<<trans1<<" trans_uk: "<<trans_uk;
		}

		//maybe not all boundary components are connected
		for(auto fhi : trimesh_.faces())
			if(_chart_index[fhi.idx()] == -2)
				_chart_index[fhi.idx()] = -1;
		for(auto fhi : trimesh_.faces())
			if(_chart_index[fhi.idx()] == -max_comp_b_)
				_chart_index[fhi.idx()] = -2;

		//update component property
		update_boundary_graph_component_property();
	}
	std::cerr<<"\nmax comp: "<<max_comp_b_;
}


template <class MeshT1, class MeshT2>
void
MatchingsAndAlignmentExtractionT<MeshT1, MeshT2>::
shortest_dual_path_between_faces_tri(TriMesh& _trimesh, const TriFH _fh0, const TriFH _fh1,
		std::vector<TriHEH>& _tri_hes)const
{
	_tri_hes.clear();

	//distance to _fh0
	TriFP<int> dist;
	_trimesh.add_property(dist);
	//previous face on the path
	TriFP<TriHEH> pre_dp;
	_trimesh.add_property(pre_dp);

	//initialize face property
	for(auto f_it = _trimesh.faces_begin(); f_it != _trimesh.faces_end(); ++f_it)
		_trimesh.property(pre_dp, *f_it) = TriHEH(-1);

	std::queue<TriFH> que;
	std::vector<bool> fh_visited(_trimesh.n_faces(), false);
	que.push(_fh0);
	fh_visited[_fh0.idx()]=true;
	bool found = false;

	while(!que.empty() && !found)
	{
		auto fh = que.front();
		que.pop();

		for(auto fhe_it = _trimesh.fh_ccwbegin(fh); fhe_it.is_valid(); ++fhe_it)
		{
			auto he_opp = _trimesh.opposite_halfedge_handle(*fhe_it);
			auto fh_opp = _trimesh.face_handle(he_opp);

			if(!fh_visited[fh_opp.idx()])
			{
				_trimesh.property(pre_dp, fh_opp) = *fhe_it;
				fh_visited[fh_opp.idx()] = true;

				if(fh_opp == _fh1)
				{
					found = true;
					break;
				}

				que.push(fh_opp);
			}
		}
	}


	//find the path
	std::vector<TriHEH> rhes;
	std::queue<TriFH> que_dp;
	que_dp.push(_fh1);

	while(!que_dp.empty())
	{
		auto fh = que_dp.front();
		que_dp.pop();

		auto he = _trimesh.property(pre_dp, fh);

		if(he.is_valid())
		{
			rhes.push_back(he);
			que_dp.push(_trimesh.face_handle(he));
		}
	}

	for(auto iter = rhes.rbegin(); iter != rhes.rend(); ++iter)
		_tri_hes.push_back(*iter);
}

}
