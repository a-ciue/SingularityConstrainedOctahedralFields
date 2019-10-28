/*
 * TangentContinuityT.hh
 *
 *  Created on: May 29, 2018
 *      Author: hliu
 */

#ifndef TANGENTCONTINUITYT_HH
#define TANGENTCONTINUITYT_HH

#include "Typedefs.hh"

namespace SCOF
{

template<class MeshT>
class TangentContinuityT
{
public:
	using TetMesh = MeshT;

	TangentContinuityT(){}
	~TangentContinuityT(){}

	//set the dual path from one halfedge to the other
	void set_dual_path(const std::vector<HFH>& _dpath){dpath_ = _dpath;}

	//set the halfedge
	void set_ith_halfedge(const int _i, const HEH _he)
	{
		if(_i==0)
			he0_ = _he;
		if(_i==1)
			he1_ = _he;
	}

	//return the dual path
	void dual_path(std::vector<HFH>& _dpath) const {_dpath = dpath_;}

	//return the ith halfedge
	HEH halfedge(const int _i) const
	{
		if(_i == 0)
			return he0_;
		else if(_i == 1)
			return he1_;

		return HEH(-1);
	}

	//return the node
	VH node(const TetMesh& _tmesh) const {return _tmesh.halfedge(he0_).from_vertex();}

	//initially set the transition of the dual faces
	bool initially_set_transitions(const TetMesh& _tmesh, CP<int>& _cell_tag,
			HFP<int> &_trans_prop, CP<int>& _chart_id) const
	{
		CH ch0 = _tmesh.incident_cell(_tmesh.opposite_halfface_handle(dpath_[0]));

		for(size_t i=0; i<dpath_.size()-1; ++i)
		{
			CH ch = _tmesh.incident_cell(dpath_[i]);
			if(!ch.is_valid())
			{
				std::cerr<<"\nError: the cell is invalid! "<<ch;
				return false;
			}
			if(_cell_tag[ch] != 0)
			{
				std::cerr<<"\nError: the cell is already conquered! "<<ch;
				return false;
			}

			_cell_tag[ch] = 1;
			_chart_id[ch] = _chart_id[ch0];
		}

		for(size_t i=0; i<dpath_.size()-1; ++i)
		{
			if(_trans_prop[dpath_[i]] != -1)
				std::cerr<<"\nError: the face already has transition!";
			_trans_prop[dpath_[i]] = 0;
			_trans_prop[_tmesh.opposite_halfface_handle(dpath_[i])] = 0;
		}

		return true;
	}

	//merge chart following the tangent continuity constraint
	void constrained_chart_mergingTC(const TetMesh& _tmesh, HFP<int> &_trans_prop, const CP<int>& _chart_id, std::map<int, int>& _real_chart_id) const
	{
		if(_trans_prop[dpath_[dpath_.size()-1]] != -1)
		{
			return;
		}
		_trans_prop[dpath_[dpath_.size()-1]] = 2;
		_trans_prop[_tmesh.opposite_halfface_handle(dpath_[dpath_.size()-1])] = 2;
//		std::cerr<<"\nface: "<<_tmesh.face_handle(dpath_[dpath_.size()-1]);
		//update chart index mapping
		std::vector<int> v_idx;
		v_idx.push_back(_real_chart_id[_chart_id[_tmesh.incident_cell(dpath_[dpath_.size()-1])]]);
		v_idx.push_back(_real_chart_id[_chart_id[_tmesh.incident_cell(_tmesh.opposite_halfface_handle(dpath_[dpath_.size()-1]))]]);
		std::sort(v_idx.begin(), v_idx.end());
		update_chart_index_mapping(_real_chart_id, v_idx[0], v_idx[1]);
//		std::cerr<<"\nnew char ids: "<<_real_chart_id[_chart_id[_tmesh.incident_cell(dpath_[dpath_.size()-1])]]
//													  <<" "<<_real_chart_id[_chart_id[_tmesh.incident_cell(_tmesh.opposite_halfface_handle(dpath_[dpath_.size()-1]))]];
	}


private:
	void update_chart_index_mapping(std::map<int, int>& _real_chart_id, const int _id0, const int _id1) const
	{
		int id_old = std::abs(_id0) < std::abs(_id1) ? _id1 : _id0;
		int id_new = std::abs(_id0) < std::abs(_id1) ? _id0 : _id1;

		for(auto it = _real_chart_id.begin(); it != _real_chart_id.end(); ++it)
			if(it->second == id_old)
				it->second = id_new;
	}

private:
	std::vector<HFH> dpath_;
	HEH he0_, he1_;
};


}

#endif /* TANGENTCONTINUITYT_HH_ */
