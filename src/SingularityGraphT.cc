/*
 * SingularityGraphT.cc
 *
 *  Created on: Dec 31, 2017
 *      Author: hliu
 */
#define SINGULARITYGRAPHT_C

#include <queue>
#include <stack>
#include "SingularityGraphT.hh"

namespace SCOF
{

template <class MeshT>
void
SingularityGraphT<MeshT>::
set_singularity_graph_component_property()
{
	int max_arc_type = *std::max_element(arc_type_.begin(), arc_type_.end());
	if(max_arc_type == 0)
		classify_singular_arc_type();

	//interior components
	int comp = 0;
	//heuristically mark the dominant component as 1
	int n_se = 0, dm_comp = 1,
			n_se_max = 0;

	for(auto hei : mesh_.halfedges())
		if(valance_[mesh_.edge_handle(hei)] != 0 &&
				arc_type_[mesh_.edge_handle(hei)] != -2 &&
				!mesh_.is_boundary(mesh_.edge_handle(hei)) &&
				sg_comp_i_[mesh_.edge_handle(hei)] == 0)
		{
			comp++;
			sg_comp_i_[mesh_.edge_handle(hei)] = comp;
			n_se = 0;

			std::queue<HEH> que;
			que.push(hei);
			que.push(mesh_.opposite_halfedge_handle(hei));
			while(!que.empty())
			{
				auto hei_cur = que.front();
				que.pop();

				auto vht = mesh_.halfedge(hei_cur).to_vertex();
				for(auto voh_it = mesh_.voh_iter(vht); voh_it.valid(); ++voh_it)
					if(valance_[mesh_.edge_handle(*voh_it)] != 0 &&
							arc_type_[mesh_.edge_handle(hei)] != -2 &&
							!mesh_.is_boundary(mesh_.edge_handle(*voh_it)) &&
							sg_comp_i_[mesh_.edge_handle(*voh_it)] == 0)
					{
						sg_comp_i_[mesh_.edge_handle(*voh_it)] = comp;
						que.push(*voh_it);

						//don't let circle to be component 1
						if(arc_type_[mesh_.edge_handle(hei)] != 2)
							n_se++;
					}
			}

			if(n_se > n_se_max)
			{
				dm_comp = comp;
				n_se_max = n_se;
			}
		}

	//swap dominant component with the first one
	std::vector<EH> tmp_ehs;
	for(auto ehi : mesh_.edges())
		if(sg_comp_i_[ehi] == 1)
			tmp_ehs.push_back(ehi);

	for(auto ehi : mesh_.edges())
		if(sg_comp_i_[ehi] == dm_comp)
			sg_comp_i_[ehi] = 1;
	for(auto ehi : tmp_ehs)
		sg_comp_i_[ehi] = dm_comp;


	//boundary components
	comp = 0, n_se = 0, dm_comp = 1, n_se_max = 0;
	for(auto hei : mesh_.halfedges())
		if(valance_[mesh_.edge_handle(hei)] != 0 &&
				mesh_.is_boundary(mesh_.edge_handle(hei)) &&
				sg_comp_b_[mesh_.edge_handle(hei)] == 0)
		{
			comp++;
			n_se = 0;
			sg_comp_b_[mesh_.edge_handle(hei)] = comp;

			std::queue<HEH> que;
			que.push(hei);
			que.push(mesh_.opposite_halfedge_handle(hei));
			while(!que.empty())
			{
				auto hei_cur = que.front();
				que.pop();

				auto vht = mesh_.halfedge(hei_cur).to_vertex();
				for(auto voh_it = mesh_.voh_iter(vht); voh_it.valid(); ++voh_it)
					if(valance_[mesh_.edge_handle(*voh_it)] != 0 &&
							mesh_.is_boundary(mesh_.edge_handle(*voh_it)) &&
							sg_comp_b_[mesh_.edge_handle(*voh_it)] == 0)
					{
						sg_comp_b_[mesh_.edge_handle(*voh_it)] = comp;
						n_se++;
						que.push(*voh_it);
					}
			}

			if(n_se > n_se_max)
			{
				dm_comp = comp;
				n_se_max = n_se;
			}
		}

	//swap dominant component with the first one
	tmp_ehs.clear();
	for(auto ehi : mesh_.edges())
		if(sg_comp_b_[ehi] == 1)
			tmp_ehs.push_back(ehi);

	for(auto ehi : mesh_.edges())
		if(sg_comp_b_[ehi] == dm_comp)
			sg_comp_b_[ehi] = 1;
	for(auto ehi : tmp_ehs)
		sg_comp_b_[ehi] = dm_comp;
}


template <class MeshT>
void
SingularityGraphT<MeshT>::
classify_singular_arc_type()
{
	for(int ll=1; ll<=max_label_; ++ll)
	{
		std::vector<VH> curve_vhs = v_sg_vhs_[ll-1];
		std::vector<EH> curve_ehs = v_sg_ehs_[ll-1];

		if(curve_ehs.empty())
			continue;
		//interior singular curve touching boundary or boundary singular curve: 1
		//interior singular curve not touching boundary: -1
		//boundary singular arc
		if(mesh_.is_boundary(curve_ehs[0]))
		{
			HEH he = mesh_.halfedge(curve_vhs[0], curve_vhs[curve_vhs.size()-1]);
			//boundary circle
			if(he.is_valid())
				for(auto i : curve_ehs)
					arc_type_[i] = 3;
			else//boundary arc that ends at singular nodes
				for(auto i : curve_ehs)
					arc_type_[i] = 4;
		}else
		{
			//interior singular arc touches the boundary
			if(mesh_.is_boundary(curve_vhs[0]) && mesh_.is_boundary(curve_vhs[curve_vhs.size()-1]))
				for(auto i : curve_ehs)
					arc_type_[i] = -2;
			else if(mesh_.is_boundary(curve_vhs[0]) || mesh_.is_boundary(curve_vhs[curve_vhs.size()-1]))
				for(auto i : curve_ehs)
					arc_type_[i] = 1;
			else
			{
				HEH he = mesh_.halfedge(curve_vhs[0], curve_vhs[curve_vhs.size()-1]);
				//interior singular circle
				if(he.is_valid() && curve_vhs.size() > 2)
					for(auto i : curve_ehs)
						arc_type_[i] = 2;
				else//interior arc that ends at singular nodes
					for(auto i : curve_ehs)
						arc_type_[i] = -1;
			}
		}
	}
}


template <class MeshT>
void
SingularityGraphT<MeshT>::
update_label_property()
{
	for(auto e_it = mesh_.e_iter(); e_it.valid(); ++e_it)
		label_[*e_it] = 0;

	int i_label = 1;
	std::vector<bool> edge_conquered(mesh_.n_edges(), false);
	for(auto e_it = mesh_.e_iter(); e_it.valid(); ++e_it)
	{
		if(valance_[*e_it] != 0 && !edge_conquered[(*e_it).idx()])
		{
			auto he0 = mesh_.halfedge_handle(*e_it, 0);
			edge_conquered[(*e_it).idx()] = true;
			label_[*e_it] = i_label;

			auto vh_s = mesh_.halfedge(he0).from_vertex();
			auto vh_t = mesh_.halfedge(he0).to_vertex();

			std::stack<VH> st_vhs;
			if(!is_end_vertex(vh_t))
				st_vhs.push(vh_t);
			if(!is_end_vertex(vh_s))
				st_vhs.push(vh_s);
			//trace back until the boundary or an edge which is conquered, then trace front
			while(!st_vhs.empty())
			{
				VH vh_cur = st_vhs.top();
				st_vhs.pop();

				for(auto voh_it = mesh_.voh_iter(vh_cur); voh_it.valid(); ++voh_it)
				{
					auto eh_og = mesh_.edge_handle(*voh_it);

					if(valance_[eh_og] != 0 && !edge_conquered[eh_og.idx()])
					{
						edge_conquered[eh_og.idx()] = true;
						label_[eh_og] = i_label;

						auto vh_next = mesh_.halfedge(*voh_it).to_vertex();
						if(!is_end_vertex(vh_next))
							st_vhs.push(vh_next);
					}
				}
			}

			i_label++;
		}
	}

	max_label_ = i_label-1;

	for(int i=1; i<=max_label_; ++i)
	{
		std::vector<EH> ehs;
		sort_edges_on_curve(i, ehs);
		v_sg_ehs_.push_back(ehs);

		std::vector<VH> vhs;
		sort_vertices_on_curve(i, vhs);
		v_sg_vhs_.push_back(vhs);
	}
}


//Singular node type(NOT COMPLETE YET)
//interior singularity node property
//TYPE:(b_val-1, b_val0, b_val1, i_val-1, i_val0, i_val1)
//10: turning point(0, 0, 0, 1, 0, 1)
//1: (0, 0, 0, 4, 0, 0)
//2: (0, 0, 0, 2, 2, 2)
//3: (0, 0, 0, 1, 3, 1)
//4: (0, 0, 0, 0, 4, 4)
//boundary singularity node property
//-1: (3, 0, 0, 0, 0, 0)
//-2: (2, 2, 1, 0, 0, 0)
//-3: (2, 2, 2, 0, 0, 0)
//-4: (2, 2, 2, 0, 0, 0)mirror case of -3
//-5: (3, 0, 3, 0, 0, 0)
//-6: (1, 2, 2, 0, 1, 0)
//-7: (0, 0, 3, 0, 3, 0)
//-8: (2, 0, 4, 0, 1, 0)
//-9: (0, 3, 2, 0, 1, 1)
//-10: (0, 3, 0, 1, 0, 0)
//-11: (0, 5, 0, 0, 0, 1)
template <class MeshT>
void
SingularityGraphT<MeshT>::
update_node_type_property()
{
	for(auto v_it = mesh_.v_iter(); v_it.valid(); ++v_it)
	{
		if(!mesh_.is_boundary(*v_it))
		{
			int n_val_ng1_i = 0, n_val1_i = 0, n_invalid = 0;
			for(auto voh_it = mesh_.voh_iter(*v_it); voh_it.valid(); ++voh_it)
			{
				auto eh = mesh_.edge_handle(*voh_it);
				if(valance_[eh] == -1)
					n_val_ng1_i++;
				else if(valance_[eh] == 1)
					n_val1_i++;
				else if(valance_[eh] != 0)
					n_invalid++;
			}

			if((n_val_ng1_i==0 && n_val1_i==0) || (n_val_ng1_i==2 && n_val1_i == 0) || (n_val_ng1_i==0 && n_val1_i == 2))
				node_type_[*v_it] = 0;
			else if(n_val_ng1_i==4 && n_val1_i==0)
				node_type_[*v_it] = 1;
			else if(n_val_ng1_i==2 && n_val1_i == 2)
				node_type_[*v_it] = 2;
			else if(n_val_ng1_i==1 && n_val1_i == 3)
				node_type_[*v_it] = 3;
			else if(n_val_ng1_i==0 && n_val1_i == 4)
				node_type_[*v_it] = 4;
			else if(n_val_ng1_i==1 && n_val1_i==1)
			{
				node_type_[*v_it] = 10;
				std::cerr<<"\nERROR: Wrong Interior Singularity Node Type! Vertex: "<<*v_it
						<<" Valance -1 number: "<<n_val_ng1_i<<" Valance 1 number: "<<n_val1_i;
			}
			else
			{
				node_type_[*v_it] = 20;
				std::cerr<<"\nERROR: Wrong Interior Singularity Node Type! Vertex: "<<*v_it
						<<" Valance -1 number: "<<n_val_ng1_i<<" Valance 1 number: "<<n_val1_i;
			}

			if(n_invalid > 0)
			{
				node_type_[*v_it] = 20;
				std::cerr<<"\nERROR: Wrong Edge Valance Type at Vertex: "<<*v_it<<" Invalid valance number: "<<n_invalid;
			}
		}else
		{
			int n_val_ng1_b = 0, n_val1_b = 0,
					n_val_ng1_i = 0, n_val1_i = 0, n_invalid = 0;
			for(auto voh_it = mesh_.voh_iter(*v_it); voh_it.valid(); ++voh_it)
			{
				EH eh = mesh_.edge_handle(*voh_it);
				if(mesh_.is_boundary(eh))
				{
					if(valance_[eh] == -1)
						n_val_ng1_b++;
					else if(valance_[eh] == 1)
						n_val1_b++;
					else if(valance_[eh] != 0)
						n_invalid++;
				}else
				{
					if(valance_[eh] == -1)
						n_val_ng1_i++;
					else if(valance_[eh] == 1)
						n_val1_i++;
					else if(valance_[eh] != 0)
						n_invalid++;
				}
			}

			if(((n_val_ng1_b==0 && n_val1_b==0 &&n_val_ng1_i ==0 && n_val1_i ==0) ||
					(n_val_ng1_b==2 && n_val1_b == 0 && n_val_ng1_i ==0 && n_val1_i==0) ||
					(n_val_ng1_b==0 && n_val1_b == 2 && n_val_ng1_i ==0 && n_val1_i==0)))
				node_type_[*v_it] = 0;
			else if(n_val_ng1_b==0 && n_val1_b == 0 && n_val_ng1_i ==0 && n_val1_i==1)
				node_type_[*v_it] = -11;
			else if(n_val_ng1_b==0 && n_val1_b == 0 && n_val_ng1_i ==1 && n_val1_i==0)
				node_type_[*v_it] = -10;
			else if((n_val_ng1_b==3 && n_val1_b==0) && n_val_ng1_i ==0 && n_val1_i==0)
				node_type_[*v_it] = -1;
			else if((n_val_ng1_b==2 && n_val1_b==1)&& n_val_ng1_i ==0 && n_val1_i==0)
				node_type_[*v_it] = -2;
			else if((n_val_ng1_b==2 && n_val1_b == 2)&& n_val_ng1_i ==0 && n_val1_i==0)
			{
				//TODO:
				//if val5 -> val3 PI ccw
					node_type_[*v_it] = -3;
				//if val3 -> val5 PI ccw
					node_type_[*v_it] = -4;
					std::cerr<<"\nHas mirror case! To be clarified...";
			}
			else if((n_val_ng1_b==3 && n_val1_b == 3)&& n_val_ng1_i ==0 && n_val1_i==0)
				node_type_[*v_it] = -5;
			else if((n_val_ng1_b==1 && n_val1_b == 2)&& n_val_ng1_i ==0 && n_val1_i==0)
				node_type_[*v_it] = -6;
			else if((n_val_ng1_b==0 && n_val1_b == 3)&& n_val_ng1_i ==0 && n_val1_i==0)
				node_type_[*v_it] = -7;
			else if((n_val_ng1_b==2 && n_val1_b == 4)&& n_val_ng1_i ==0 && n_val1_i==0)
				node_type_[*v_it] = -8;
			else if((n_val_ng1_b==0 && n_val1_b == 2)&& n_val_ng1_i ==0 && n_val1_i==1)
				node_type_[*v_it] = -9;
			else
			{
				node_type_[*v_it] = 20;
				std::cerr<<"\nERROR: Wrong Boundary Singularity Node Type! Vertex: "<<*v_it<<" Valance -1-b: "<<n_val_ng1_b<<" Valance 1-b: "<<n_val1_b
						<<" Valance3: "<<n_val_ng1_i<<" Valance5: "<<n_val1_i;
			}

			if(n_invalid > 0)
			{
				node_type_[*v_it] = 20;
				std::cerr<<"\nERROR: Wrong Edge Valance Type at Vertex: "<<*v_it<<" Invalid valance number: "<<n_invalid;
			}
		}
	}
}


template <class MeshT>
void
SingularityGraphT<MeshT>::
update_singular_vertex_property()
{
	for(auto ehi : mesh_.edges())
		if(valance_[ehi] != 0)
		{
			VH vh0 = mesh_.edge(ehi).from_vertex();
			is_singular_vt_[vh0] = true;
			VH vh1 = mesh_.edge(ehi).to_vertex();
			is_singular_vt_[vh1] = true;
		}
}


template <class MeshT>
void
SingularityGraphT<MeshT>::
measure_the_shortest_curve_length() const
{
	int max_label = *std::max_element(label_.begin(), label_.end());
	double min_length = DBL_MAX;
	int min_label = 0;
	for(int ll=1; ll<=max_label; ++ll)
	{
		std::vector<VH> curve_vertices;
		sort_vertices_on_curve(ll, curve_vertices);

		double length = 0.0;
		for(size_t i=0; i<curve_vertices.size()-1; ++i)
			length += distance(curve_vertices[i], curve_vertices[i+1]);

		if(min_length > length)
		{
			min_length = length;
			min_label = ll;
		}
	}

	std::cerr<<"\nMinimum singular arc length: "<<min_length<<" label: "<<min_label;
}


template <class MeshT>
std::vector<EH>
SingularityGraphT<MeshT>::
get_singular_edges_of_label(const int _label) const
{
	if(_label > max_label_ || _label <= 0)
	{
		std::cerr<<"\nInvalid label!";
		std::vector<EH> ehs(1, EH(0));
		return ehs;
	}

	return v_sg_ehs_[_label-1];
}


template <class MeshT>
std::vector<VH>
SingularityGraphT<MeshT>::
get_singular_vertices_of_label(const int _label) const
{
	if(_label > max_label_ || _label <= 0)
	{
		std::cerr<<"\nInvalid label!";
		std::vector<VH> vhs(1, VH(0));
		return vhs;
	}

	return v_sg_vhs_[_label-1];
}


template <class MeshT>
void
SingularityGraphT<MeshT>::
sort_vertices_on_curve(const int _label, std::vector<VH>& _vhs) const
{
	EH eh_s(-1);
	for(auto e_it = mesh_.e_iter(); e_it.valid(); ++e_it)
		if(label_[*e_it] == _label)
		{
			eh_s = *e_it;
			break;
		}

	VH vh_f = mesh_.edge(eh_s).from_vertex();
	VH vh_t = mesh_.edge(eh_s).to_vertex();
	std::vector<bool> vh_visited(mesh_.n_vertices(), false);
	vh_visited[vh_f.idx()] = true;
	vh_visited[vh_t.idx()] = true;

	std::vector<VH> vhs0, vhs1;
	vhs1.push_back(vh_t);
	vhs0.push_back(vh_f);
	std::queue<VH> que_vhs;
	que_vhs.push(vh_t);

	while(!que_vhs.empty())
	{
		VH vh_cur = que_vhs.front();
		que_vhs.pop();

		for(auto voh_it = mesh_.voh_iter(vh_cur); voh_it.valid(); ++voh_it)
		{
			EH eh_og = mesh_.edge_handle(*voh_it);
			VH vh_next = mesh_.halfedge(*voh_it).to_vertex();
			if(label_[eh_og] == _label && !vh_visited[vh_next.idx()])
			{
				vh_visited[vh_next.idx()] = true;
				vhs1.push_back(vh_next);

				que_vhs.push(vh_next);
			}
		}
	}

	que_vhs.push(vh_f);
	while(!que_vhs.empty())
	{
		VH vh_cur = que_vhs.front();
		que_vhs.pop();

		for(auto voh_it = mesh_.voh_iter(vh_cur); voh_it.valid(); ++voh_it)
		{
			EH eh_og = mesh_.edge_handle(*voh_it);
			VH vh_next = mesh_.halfedge(*voh_it).to_vertex();
			if(label_[eh_og] == _label && !vh_visited[vh_next.idx()])
			{
				vh_visited[vh_next.idx()] = true;
				vhs0.push_back(vh_next);

				que_vhs.push(vh_next);
			}
		}
	}

	for(auto iter = vhs0.rbegin(); iter!= vhs0.rend(); ++iter)
		_vhs.push_back(*iter);
	for(auto iter = vhs1.begin(); iter!= vhs1.end(); ++iter)
		_vhs.push_back(*iter);
}


template <class MeshT>
void
SingularityGraphT<MeshT>::
sort_edges_on_curve(int _label, std::vector<EH>& _ehs) const
{
	std::vector<VH> vhs;
	sort_vertices_on_curve(_label, vhs);

	for(size_t i=0; i<vhs.size()-1; ++i)
		_ehs.push_back(mesh_.edge_handle(mesh_.halfedge(vhs[i], vhs[i+1])));

	auto he = mesh_.halfedge(vhs[vhs.size()-1], vhs[0]);
	if(he.is_valid() && vhs.size() > 2)
		_ehs.push_back(mesh_.edge_handle(he));
}


//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
template <class MeshT>
bool
SingularityGraphT<MeshT>::
is_end_vertex(VH _vh) const
{
	int val = 0;
	for(auto voh_it = mesh_.voh_iter(_vh); voh_it.valid(); ++voh_it)
	{
		auto eh = mesh_.edge_handle(*voh_it);
		if(valance_[eh] != 0)//count the number of singularity line connected to the node
			val++;
	}

	if(val == 2)
		return false;

	return true;
}


template <class MeshT>
double
SingularityGraphT<MeshT>::
distance(VH _vh0, VH _vh1) const
{
	auto vec = mesh_.vertex(_vh1) - mesh_.vertex(_vh0);
	return vec.norm();
}

}
