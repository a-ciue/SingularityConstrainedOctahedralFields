/*
 * CornerT.hh
 *
 *  Created on: Oct 9, 2017
 *      Author: hliu
 */

#ifndef CORNERT_HH
#define CORNERT_HH

#include "TransitionQuaternionEigen.hh"

#include "Typedefs.hh"

namespace SCOF
{

template <class MeshT>
class CornerT
{
public:
	using TetMesh = MeshT;

	CornerT(){}
	~CornerT(){}

	//set the ith dual path from the ith halfedge to the base cell of the corner
	void set_ith_dual_path(const int _i, const std::vector<HFH>& _dpath)
	{
		if(_i==0)
			dpath0_ = _dpath;
		if(_i==1)
			dpath1_ = _dpath;
		if(_i==2)
			dpath2_ = _dpath;
	}

	//set the ith halfedge of the corner
	void set_ith_halfedge(const int i, const HEH _he)
	{
		if(i==0)
			he0_ = _he;
		if(i==1)
			he1_ = _he;
		if(i==2)
			he2_ = _he;
	}

	//return the singular node
	VH node(const TetMesh& _tmesh) const{return _tmesh.halfedge(he0_).from_vertex();}

	//return the base cell
	CH base_cell(const TetMesh& _tmesh) const {return _tmesh.incident_cell(dpath0_[dpath0_.size()-1]);}

	//return the ith dual path
	void dual_path(const int _i, std::vector<HFH>& _dpath)const
	{
		if(_i == 0)
			_dpath = dpath0_;
		else if(_i == 1)
			_dpath = dpath1_;
		else if(_i == 2)
			_dpath = dpath2_;
	}

	//return the ith halfedge
	HEH halfedge(const int _i) const
	{
		if(_i == 0)
			return he0_;
		else if(_i == 1)
			return he1_;
		else if(_i == 2)
			return he2_;

		return HEH(-1);
	}

	//initially set transitions of the dual faces
	bool initially_set_transitions(const TetMesh& _tmesh, CP<int>& _cell_tag,
			const EP<int>& _arc_type, HFP<int>&_trans_prop, CP<int>& _chart_id) const
	{
		if(dpath0_.empty() || dpath1_.empty() || dpath2_.empty())
		{
			std::cerr<<"\nError: dual path is Empty!";
			return false;
		}
		//set cell_tag and chart_id
		CH ch0 = _tmesh.incident_cell(_tmesh.opposite_halfface_handle(dpath0_[0]));
		for(size_t i=0; i<dpath0_.size()-1; ++i)
		{
			CH ch = _tmesh.incident_cell(dpath0_[i]);
			if(!ch.is_valid())
			{
				std::cerr<<"\nError: the cell is invalid! "<<ch<<" hf: "<<dpath0_[i];
				return false;
			}
			if(_cell_tag[ch] != 0)
			{
				std::cerr<<"\nError: the cell is already conquered! "<<ch<<" hf: "<<dpath0_[i];
				for(auto vhi : _tmesh.get_cell_vertices(ch))
					std::cerr<<" "<<vhi;
				return false;
			}

			_cell_tag[ch] = _arc_type[_tmesh.edge_handle(he0_)];
			_chart_id[ch] = _chart_id[ch0];
		}

		CH ch1 = _tmesh.incident_cell(_tmesh.opposite_halfface_handle(dpath1_[0]));
		for(size_t i=0; i<dpath1_.size()-1; ++i)
		{
			CH ch = _tmesh.incident_cell(dpath1_[i]);
			if(!ch.is_valid())
			{
				std::cerr<<"\nError: the cell is invalid! "<<ch<<" hf: "<<dpath1_[i];
				return false;
			}

			if(_cell_tag[ch] != 0)
			{
				std::cerr<<"\nError: the cell is already conquered! "<<ch<<" tag: "<<_cell_tag[ch];
				return false;
			}

			_cell_tag[ch] = _arc_type[_tmesh.edge_handle(he1_)];
			_chart_id[ch] = _chart_id[ch1];
		}

		CH ch2 = _tmesh.incident_cell(_tmesh.opposite_halfface_handle(dpath2_[0]));
		for(size_t i=0; i<dpath2_.size()-1; ++i)
		{
			CH ch = _tmesh.incident_cell(dpath2_[i]);
			if(!ch.is_valid())
			{
				std::cerr<<"\nError: the cell is invalid! "<<ch<<" hf: "<<dpath2_[i];
				return false;
			}

			if(_cell_tag[ch] != 0)
			{
				std::cerr<<"\nError: the cell is already conquered! "<<ch<<" tag: "<<_cell_tag[ch];
				return false;
			}

			_cell_tag[ch] = _arc_type[_tmesh.edge_handle(he2_)];
			_chart_id[ch] = _chart_id[ch2];
		}

		//set transition
		for(size_t i=0; i<dpath0_.size()-1; ++i)
		{
			if(_trans_prop[dpath0_[i]] != -1)
				std::cerr<<"\nError: the face already has transition! "<<_tmesh.face_handle(dpath0_[i]);
			_trans_prop[dpath0_[i]] = 0;
			_trans_prop[_tmesh.opposite_halfface_handle(dpath0_[i])] = 0;
		}

		for(size_t i=0; i<dpath1_.size()-1; ++i)
		{
			if(_trans_prop[dpath1_[i]] != -1)
				std::cerr<<"\nError: the face already has transition! "<<_tmesh.face_handle(dpath1_[i]);
			_trans_prop[dpath1_[i]] = 0;
			_trans_prop[_tmesh.opposite_halfface_handle(dpath1_[i])] = 0;
		}

		for(size_t i=0; i<dpath2_.size()-1; ++i)
		{
			if(_trans_prop[dpath2_[i]] != -1)
				std::cerr<<"\nError: the face already has transition! "<<_tmesh.face_handle(dpath2_[i]);
			_trans_prop[dpath2_[i]] = 0;
			_trans_prop[_tmesh.opposite_halfface_handle(dpath2_[i])] = 0;
		}

		//connect the corner to boundary
		CH ch = _tmesh.incident_cell(dpath2_[dpath2_.size()-1]);
		if(_arc_type[_tmesh.edge_handle(he0_)] == 1)
		{
			_trans_prop[dpath0_[dpath0_.size()-1]] = 0;
			_trans_prop[_tmesh.opposite_halfface_handle(dpath0_[dpath0_.size()-1])] = 0;

			_cell_tag[ch] = 1;
			_chart_id[ch] = _chart_id[ch0];
		}else if(_arc_type[_tmesh.edge_handle(he1_)] == 1)
		{
			_trans_prop[dpath1_[dpath1_.size()-1]] = 0;
			_trans_prop[_tmesh.opposite_halfface_handle(dpath1_[dpath1_.size()-1])] = 0;

			_cell_tag[ch] = 1;
			_chart_id[ch] = _chart_id[ch1];

		}else if(_arc_type[_tmesh.edge_handle(he2_)] == 1)
		{
			_trans_prop[dpath2_[dpath2_.size()-1]] = 0;
			_trans_prop[_tmesh.opposite_halfface_handle(dpath2_[dpath2_.size()-1])] = 0;

			_cell_tag[ch] = 1;
			_chart_id[ch] = _chart_id[ch2];
		}

		return true;
	}

	//initially set the transition for the second dual path if it is initially connectable(like the sphere case)
	HEH initial_constrained_chart_mergingC(const TetMesh& _tmesh, const int _count, const TransitionQuaternion& _tq, HFP<int>&_trans_prop,
			const CP<int>& _chart_id, std::map<int, int>& _real_chart_id) const
	{
		//get transition(any axis different from ±x is feasible)
		int trans = 0, count = 0;
		for(int i=0; i<24; ++i)
		{
			int axis_after = axis_after_transition(_tq, 0, i);
			if(axis_after != 0 && axis_after != 1)
			{
				trans = i;
				if(count == _count)
					break;
				count++;
			}
		}
		//assign transition and update chart id
		if(_trans_prop[dpath0_[dpath0_.size()-1]] != -1)
		{
			_trans_prop[dpath1_[dpath1_.size()-1]] = trans;
			_trans_prop[_tmesh.opposite_halfface_handle(dpath1_[dpath1_.size()-1])] = _tq.inverse_transition_idx(trans);

			//update chart index mapping
			auto idx_old = _real_chart_id[_chart_id[_tmesh.incident_cell(_tmesh.opposite_halfface_handle(dpath1_[0]))]];
			auto idx_new = _real_chart_id[_chart_id[_tmesh.incident_cell(_tmesh.opposite_halfface_handle(dpath0_[0]))]];
			update_chart_index_mapping(_real_chart_id, idx_old, idx_new);
//			std::cerr<<"\nnew char ids: "<<_real_chart_id[_chart_id[_tmesh.incident_cell(_tmesh.opposite_halfface_handle(dpath0_[0]))]]
//                     <<" "<<_real_chart_id[_chart_id[_tmesh.incident_cell(_tmesh.opposite_halfface_handle(dpath1_[0]))]]
//                     <<" "<<_real_chart_id[_chart_id[_tmesh.incident_cell(_tmesh.opposite_halfface_handle(dpath2_[0]))]];
			return he1_;
		}else if(_trans_prop[dpath1_[dpath1_.size()-1]] != -1)
		{
			_trans_prop[dpath2_[dpath2_.size()-1]] = trans;
			_trans_prop[_tmesh.opposite_halfface_handle(dpath2_[dpath2_.size()-1])] = _tq.inverse_transition_idx(trans);

			//update chart index mapping
			auto idx_old = _real_chart_id[_chart_id[_tmesh.incident_cell(_tmesh.opposite_halfface_handle(dpath2_[0]))]];
			auto idx_new = _real_chart_id[_chart_id[_tmesh.incident_cell(_tmesh.opposite_halfface_handle(dpath1_[0]))]];
			update_chart_index_mapping(_real_chart_id, idx_old, idx_new);
//			std::cerr<<"\nnew char ids: "<<_real_chart_id[_chart_id[_tmesh.incident_cell(_tmesh.opposite_halfface_handle(dpath0_[0]))]]
//                     <<" "<<_real_chart_id[_chart_id[_tmesh.incident_cell(_tmesh.opposite_halfface_handle(dpath1_[0]))]]
//                     <<" "<<_real_chart_id[_chart_id[_tmesh.incident_cell(_tmesh.opposite_halfface_handle(dpath2_[0]))]];
			return he2_;
		}else if(_trans_prop[dpath2_[dpath2_.size()-1]] != -1)
		{
			_trans_prop[dpath0_[dpath0_.size()-1]] = trans;
			_trans_prop[_tmesh.opposite_halfface_handle(dpath0_[dpath0_.size()-1])] = _tq.inverse_transition_idx(trans);

			//update chart index mapping
			auto idx_old = _real_chart_id[_chart_id[_tmesh.incident_cell(_tmesh.opposite_halfface_handle(dpath0_[0]))]];
			auto idx_new = _real_chart_id[_chart_id[_tmesh.incident_cell(_tmesh.opposite_halfface_handle(dpath2_[0]))]];
			update_chart_index_mapping(_real_chart_id, idx_old, idx_new);
//			std::cerr<<"\nnew char ids: "<<_real_chart_id[_chart_id[_tmesh.incident_cell(_tmesh.opposite_halfface_handle(dpath0_[0]))]]
//                     <<" "<<_real_chart_id[_chart_id[_tmesh.incident_cell(_tmesh.opposite_halfface_handle(dpath1_[0]))]]
//                     <<" "<<_real_chart_id[_chart_id[_tmesh.incident_cell(_tmesh.opposite_halfface_handle(dpath2_[0]))]];
			return he0_;
		}

		return HEH(-1);
	}


	//if two axes are known, get the unknown axis and assign transition to the dual path
	HEH constrained_chart_mergingC(const TetMesh& _tmesh, int _count, const  TransitionQuaternion& _tq, HFP<int>&_trans_prop,
			const CP<int>& _chart_id, std::map<int, int>& _real_chart_id)const
	{
		if(_trans_prop[dpath0_[dpath0_.size()-1]] == -1)
		{
			int axis1 = axis_after_transition(_tq, 0, _trans_prop[dpath1_[dpath1_.size()-1]]);
			int axis2 = axis_after_transition(_tq, 0, _trans_prop[dpath2_[dpath2_.size()-1]]);
			int axis0 = the_third_axis(axis1, axis2);

			int trans0 = 0;
			int count = 0;
			for(int i=0; i<24; ++i)
				if(axis_after_transition(_tq, 0, i) == axis0)
				{
					trans0 = i;
					if(count == _count)
						break;
					count++;
				}
			_trans_prop[dpath0_[dpath0_.size()-1]] = trans0;
			_trans_prop[_tmesh.opposite_halfface_handle(dpath0_[dpath0_.size()-1])] = _tq.inverse_transition_idx(trans0);

			//update chart index mapping
			auto idx_old = _real_chart_id[_chart_id[_tmesh.incident_cell(_tmesh.opposite_halfface_handle(dpath0_[0]))]];
			auto idx_new = _real_chart_id[_chart_id[_tmesh.incident_cell(_tmesh.opposite_halfface_handle(dpath1_[0]))]];
			update_chart_index_mapping(_real_chart_id, idx_old, idx_new);
//			std::cerr<<"\nnew char ids: "<<_real_chart_id[_chart_id[_tmesh.incident_cell(_tmesh.opposite_halfface_handle(dpath0_[0]))]]
//                     <<" "<<_real_chart_id[_chart_id[_tmesh.incident_cell(_tmesh.opposite_halfface_handle(dpath1_[0]))]]
//                     <<" "<<_real_chart_id[_chart_id[_tmesh.incident_cell(_tmesh.opposite_halfface_handle(dpath2_[0]))]];

			return he0_;
		}else if(_trans_prop[dpath1_[dpath1_.size()-1]] == -1)
		{
			int axis2 = axis_after_transition(_tq, 0, _trans_prop[dpath2_[dpath2_.size()-1]]);
			int axis0 = axis_after_transition(_tq, 0, _trans_prop[dpath0_[dpath0_.size()-1]]);
			int axis1 = the_third_axis(axis2, axis0);

			int trans1 = 0;
			int count = 0;
			for(int i=0; i<24; ++i)
				if(axis_after_transition(_tq, 0, i) == axis1)
				{
					trans1 = i;
					if(count == _count)
						break;
					count++;
				}
			_trans_prop[dpath1_[dpath1_.size()-1]] = trans1;
			_trans_prop[_tmesh.opposite_halfface_handle(dpath1_[dpath1_.size()-1])] = _tq.inverse_transition_idx(trans1);

			//update chart index mapping
			auto idx_old = _real_chart_id[_chart_id[_tmesh.incident_cell(_tmesh.opposite_halfface_handle(dpath1_[0]))]];
			auto idx_new = _real_chart_id[_chart_id[_tmesh.incident_cell(_tmesh.opposite_halfface_handle(dpath2_[0]))]];
			update_chart_index_mapping(_real_chart_id, idx_old, idx_new);
//			std::cerr<<"\nnew char ids: "<<_real_chart_id[_chart_id[_tmesh.incident_cell(_tmesh.opposite_halfface_handle(dpath0_[0]))]]
//                     <<" "<<_real_chart_id[_chart_id[_tmesh.incident_cell(_tmesh.opposite_halfface_handle(dpath1_[0]))]]
//                     <<" "<<_real_chart_id[_chart_id[_tmesh.incident_cell(_tmesh.opposite_halfface_handle(dpath2_[0]))]];

			return he1_;
		}else if(_trans_prop[dpath2_[dpath2_.size()-1]] == -1)
		{
			int axis0 = axis_after_transition(_tq, 0, _trans_prop[dpath0_[dpath0_.size()-1]]);
			int axis1 = axis_after_transition(_tq, 0, _trans_prop[dpath1_[dpath1_.size()-1]]);
			int axis2 = the_third_axis(axis0, axis1);

			int trans2 = 0;
			int count = 0;
			for(int i=0; i<24; ++i)
				if(axis_after_transition(_tq, 0, i) == axis2)
				{
					trans2 = i;
					if(count == _count)
						break;
					count++;
				}

			_trans_prop[dpath2_[dpath2_.size()-1]] = trans2;
			_trans_prop[_tmesh.opposite_halfface_handle(dpath2_[dpath2_.size()-1])] = _tq.inverse_transition_idx(trans2);

			//update chart index mapping
			auto idx_old = _real_chart_id[_chart_id[_tmesh.incident_cell(_tmesh.opposite_halfface_handle(dpath2_[0]))]];
			auto idx_new = _real_chart_id[_chart_id[_tmesh.incident_cell(_tmesh.opposite_halfface_handle(dpath0_[0]))]];
			update_chart_index_mapping(_real_chart_id, idx_old, idx_new);
//			std::cerr<<"\nnew char ids: "<<_real_chart_id[_chart_id[_tmesh.incident_cell(_tmesh.opposite_halfface_handle(dpath0_[0]))]]
//					 <<" "<<_real_chart_id[_chart_id[_tmesh.incident_cell(_tmesh.opposite_halfface_handle(dpath1_[0]))]]
//					 <<" "<<_real_chart_id[_chart_id[_tmesh.incident_cell(_tmesh.opposite_halfface_handle(dpath2_[0]))]];

			return he2_;
		}

		return HEH(-1);
	}

	//check if initial connecting is necessary
	bool initial_connectable(const TetMesh& _tmesh, const HFP<int>&_trans_prop,
			const CP<int>& _chart_id, std::map<int, int>& _real_chart_id)const
	{
		std::vector<int> ids;
		ids.push_back(_real_chart_id[_chart_id[_tmesh.incident_cell(_tmesh.opposite_halfface_handle(dpath0_[0]))]]);
		ids.push_back(_real_chart_id[_chart_id[_tmesh.incident_cell(_tmesh.opposite_halfface_handle(dpath1_[0]))]]);
		ids.push_back(_real_chart_id[_chart_id[_tmesh.incident_cell(_tmesh.opposite_halfface_handle(dpath2_[0]))]]);

		if(n_connected(_trans_prop) == 1)
		{
			std::sort(ids.begin(), ids.end());
			if(ids[0] < 0 && ids[1] > 0)
				return true;
		}

		return false;
	}

	//check if chart merging is legal
	bool connectable(const TetMesh& _tmesh, const HFP<int>&_trans_prop, const CP<int>& _chart_id,
			std::map<int, int>& _real_chart_id)const
	{
		std::vector<int> ids;
		ids.push_back(_real_chart_id[_chart_id[_tmesh.incident_cell(_tmesh.opposite_halfface_handle(dpath0_[0]))]]);
		ids.push_back(_real_chart_id[_chart_id[_tmesh.incident_cell(_tmesh.opposite_halfface_handle(dpath1_[0]))]]);
		ids.push_back(_real_chart_id[_chart_id[_tmesh.incident_cell(_tmesh.opposite_halfface_handle(dpath2_[0]))]]);

		if(n_connected(_trans_prop) == 2)
		{
//			std::cerr<<"\nids: "<<ids[0]<<" "<<ids[1]<<" "<<ids[2];
			if(!dpath0_.empty() && _trans_prop[dpath0_[dpath0_.size()-1]] == -1 && ids[0] > 0)
					return true;
			if(!dpath1_.empty() && _trans_prop[dpath1_[dpath1_.size()-1]] == -1 && ids[1] > 0)
					return true;
			if(!dpath2_.empty() && _trans_prop[dpath2_[dpath2_.size()-1]] == -1 && ids[2] > 0)
					return true;
		}

		return false;
	}

	//return how many dual paths' transition are known
	int n_connected(const HFP<int>&_trans_prop) const
	{
		int n_known = 0;
		if(!dpath0_.empty() && _trans_prop[dpath0_[dpath0_.size()-1]] != -1)
			n_known++;

		if(!dpath1_.empty() && _trans_prop[dpath1_[dpath1_.size()-1]] != -1)
			n_known++;

		if(!dpath2_.empty() && _trans_prop[dpath2_[dpath2_.size()-1]] != -1)
			n_known++;

		return n_known;
	}

	int n_connected_debug(const TetMesh& _tmesh, const CP<int>& _cell_tag, const HFP<int>&_trans_prop)const
	{
		int n_known = 0;
		if(_trans_prop[dpath0_[dpath0_.size()-1]] != -1)
			n_known++;
		if(_trans_prop[dpath1_[dpath1_.size()-1]] != -1)
			n_known++;
		if(_trans_prop[dpath2_[dpath2_.size()-1]] != -1)
			n_known++;
		CH ch = _tmesh.incident_cell(dpath0_[dpath0_.size()-1]);
		std::cerr<<"\nbase cell tag: "<<_cell_tag[ch]<<" ch: "<<ch<<"   ";
		auto hfs = _tmesh.cell(ch).halffaces();
		for(int i=0; i<4; ++i)
		{
			CH ch_opp = _tmesh.incident_cell(_tmesh.opposite_halfface_handle(hfs[i]));
			std::cerr<<" "<<ch_opp<<" ";
		}
		std::cerr<<"   "<<n_known;

		return n_known;
	}

	//return the number of singular arcs that touches the boundary
	int n_touching_boundary(const TetMesh& _tmesh, const EP<int>& _arc_type) const
	{
		int count = 0;
		if(_arc_type[_tmesh.edge_handle(he0_)] == 1)
			count++;
		if(_arc_type[_tmesh.edge_handle(he1_)] == 1)
			count++;
		if(_arc_type[_tmesh.edge_handle(he2_)] == 1)
			count++;
		return count;
	}


private:
	//with the first axis known, counterclockwisely choose the second axis(1 out of 4)
	int the_second_axis(const int _axis0) const
	{
		int pos = _axis0/2;

		int axis1 = 0;
		if(pos == 0)
			axis1 = 2;
		else
			axis1 = 2*(pos-1);

		return axis1;
	}


	//get the third aixs
	int the_third_axis(const int _axis0, const int _axis1) const
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


	//compute the axis after transition
	int axis_after_transition(const  TransitionQuaternion& _tq, const int _axis_start, const int _transition) const
	{
		int axis = 0;

		auto m1 = _tq.transition_matrix_int(_transition);

		int axes_mapped[3] = {0};
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

		if(res == 0)
			axis = axes_mapped[pos];
		else if(res == 1)
				axis = (axes_mapped[pos] %2 == 0) ? axes_mapped[pos] + 1 : axes_mapped[pos] - 1;

		return axis;
	}

	//update the chart index
	void update_chart_index_mapping(std::map<int, int>& _real_chart_id, const int _id0, const int _id1) const
	{
		int id_old = std::abs(_id0) < std::abs(_id1) ? _id1 : _id0;
		int id_new = std::abs(_id0) < std::abs(_id1) ? _id0 : _id1;
		for(auto it = _real_chart_id.begin(); it != _real_chart_id.end(); ++it)
			if(it->second == id_old)
				it->second = id_new;
	}

private:
	std::vector<HFH> dpath0_, dpath1_, dpath2_;
	HEH he0_, he1_, he2_;
};

}

#endif /* CORNERT_HH_ */
