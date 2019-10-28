/*
 * MatchingsAndAlignmentExtractionT.hh
 *
 *  Created on: Feb 28, 2018
 *      Author: hliu
 */

#ifndef MATCHINGSANDALIGNMENTEXTRACTIONT_HH
#define MATCHINGSANDALIGNMENTEXTRACTIONT_HH

//#include <ACG/Utils/HaltonColors.hh>
#include <queue>
#include "CornerT.hh"
#include "TangentContinuityT.hh"
#include "TransitionQuaternionEigen.hh"

#include "SingularityGraphT.hh"
#include "Typedefs.hh"
#include "StopWatch.hh"


#define PI 3.14159265358979323846

namespace OVM = OpenVolumeMesh;

namespace SCOF
{

template <class MeshT1, class MeshT2>
class MatchingsAndAlignmentExtractionT
{
public:
	using TetMesh = MeshT1;
	using TriMesh = MeshT2;

	//Constructor
	MatchingsAndAlignmentExtractionT(TetMesh& _mesh, SingularityGraphT<TetMesh>& _sg,
			OVM::StatusAttrib &_status/*, OVM::ColorAttrib<ACG::Vec4f>& _color*/)
	: mesh_(_mesh),
	  valance_(mesh_.template request_edge_property<int>("edge_valance")),

	  node_type_(mesh_.template request_vertex_property<int>("singular_node")),
	  label_(mesh_.template request_edge_property<int>("singularity_line_label")),
	  is_singular_vt_(mesh_.template request_vertex_property<bool>("is_singular_vertex")),

	  axes_prop_(mesh_.template request_cell_property<std::map<HEH, int> >("EdgeAxis")),
	  trans_prop_(mesh_.template request_halfface_property<int>("HalffaceTransiton")),
	  normal_prop_(mesh_.template request_halfface_property<int>("BoundaryFaceNormal")),

	  chart_id_(mesh_.template request_cell_property<int>("chart_index")),
	  cell_tag_(mesh_.template request_cell_property<int>("cell_tag")),
	  edge_tag_(mesh_.template request_edge_property<bool>("edge_tag")),

	  sg_comp_b_(mesh_.template request_edge_property<int>("boundary_singular_graph_component")),
	  sg_comp_i_(mesh_.template request_edge_property<int>("interior_singular_graph_component")),
	  arc_type_(mesh_.template request_edge_property<int>("singular_arc_type")),
	  status_(_status),
//	  color_(_color),
	  sg_(_sg)
	{
		mesh_.set_persistent(axes_prop_, true);
		mesh_.set_persistent(trans_prop_, true);
		mesh_.set_persistent(normal_prop_, true);

		initialize_singularity_graph();
		initialize_chart_index();

		open_edges_.clear();
	}

	void determineMatchings();
	void check_transitions();

	~MatchingsAndAlignmentExtractionT(){}

private:
	bool check_input()const;
    bool check_input_singularity_graph()const;
    //
	void initialize_singularity_graph();
	//initialize chart index mapping
	void initialize_chart_index();
	//
	void initialize_output_properties();
	//create surface mesh and set properties
	void initialize_trimesh_and_properties();
    void create_trimesh();
    void set_trimesh_property();

	//set axis alignments
	void set_singular_arc_alignment_constraints();
	//merge singular arc charts
	void merge_singular_arc_charts();
	void merge_and_zip_charts();
    void constrained_chart_mergingT();
    void tangent_continuity(const EH _eh, const EH _eh_next);
    void interior_singular_edge_closing(const EH _eh);
    void shortest_dual_path_from_cell_to_edge_restricting_to_the_vertex(const CH _ch, const EH _eh, std::vector<HFH>& _hfs, bool _block_free = false)const;
    bool cell_has_edge(const CH _ch, const EH _eh)const;
    VH common_vertex_handle(const EH _eh0, const EH _eh1)const;
    HEH common_halfedge_handle(const HFH _hf0, const HFH _hf1)const;
    void cell_edges(const CH _ch, std::vector<EH>& _ehs)const;
    //add input constraints
	void constrained_chart_mergingB();
    void add_flat_sector_constraints(const TriVH _vh);
    void add_irregular_sector_constraints(const std::vector<TriVH>& _vhs, const double _angle);
    bool are_halfedges_the_same_direction(const TriHEH _he0, const TriHEH _he1) const;
    int sector_transition(const double _angle) const;

	//merge boundary charts
    void merge_boundary_charts();
    //dual spanning forest and zippering on surface mesh
    void get_seeds_of_boundaries(std::vector<TriFH>& _seeds);
    void dual_spanning_forest_on_trimesh(const std::vector<TriFH>& _fhs, std::vector<int>& _chart_index);
    void solve_transitions_on_trimesh();
    bool is_removable(const TriVH _vh, TriHEH &_he)const;
    //check
    bool check_transitions_of_trimesh();

    //merge volume charts
    void merge_volume_charts();
    //add corner constraints
	void add_corner_constraints(std::map<EH, std::vector<int> >& _eh_to_corners);
	void add_tangent_continuity_constraints();
	//merge interior singular tubes to the boundary
    void constrained_chart_mergingTB();
    void get_path_from_fan_to_boundary_restricting_to_the_vertex(const EH _eh, std::vector<HFH>& _hfs)const;
    //dual spanning forest and zippering
   void dual_spanning_forest(const int _cell_tag);
   void zip_edges(const bool _mark_edge = true);
   bool is_removable(const HEH _he, HFH &_hf)const;
   void regular_edge_closing(const HEH _he, const HFH _hf);
    //iteratively process corners and decidable patches
    void local_decidable_chart_merging(	std::map<EH, std::vector<int> >& _eh_to_corners);
    void get_end_edges_of_an_arc(const HEH _he_first, std::vector<EH>& _ehs)const;
    //locally decide patch matchings
    void constrained_chart_mergingP();
    void get_patches(HFP<int>& _patch_id, EP<bool>& _edge_block, std::vector<std::vector<HFH> >& _patches)const;
    void get_singular_vertices_on_patch(const std::vector<HFH>& _patch, std::vector<VH>& _patch_sg_vhs)const;
    void get_first_constraints_on_patch(const HFP<int>& _patch_id,  const std::vector<VH>& _patch_sg_vhs,
    		std::vector<HFH>& _hfs, std::vector<std::vector<int> >& _v_trans)const;
    bool get_second_constraints_on_patch(const HFP<int>& _patch_id, const EP<bool>& _edge_block, const std::vector<HFH>& _hfs,
    		const std::vector<std::vector<int> >& _v_trans, HFH& _hf_breach, int& _unique_trans)const;
    void halffaces_between_two_halffaces(const HEH _he, const HFH _hf0, const HFH _hf1, std::vector<HFH>& _hfs)const;
    void shortest_dual_path_from_boundary_to_halfface_on_patch(const HFH _hf, HFH& _hf_tb, const HFP<int>& _patch_id,
    		const EP<bool>& _edge_block, std::vector<HFH>& _hfs)const;

    //apply exhaustive search if necessary
     void apply_exhaustive_search();
     void get_keys(std::vector<HFH>& _hfs, std::vector<std::vector<int> >& _v_trans);
     void exhaustive_search(const std::vector<FH>& _uk_fhs, const std::vector<HFH>& _hfs, const std::vector<std::vector<int> >&_v_trans);
     void get_indexes(const int _decimal, std::vector<int>& _indexes)const;
     void update_open_edge_property_on_boundary(EP<bool>& _ioe_prop)const;

    int the_third_axis(const int _axis0, const int _axis1) const;
    bool are_halfedges_the_same_direction(const HEH _he0, const HEH _he1)const;
    double distance(const VH _vh0, const VH _vh1)const;
    //check
    bool check_solved_regular_edges_transitions(const bool _mute = false);

    //update real chart id
    void update_chart_index_mapping(int _id0, int _id1);

    //update component property
    void update_interior_graph_component_property(const int _id0, const int _id1);
    void update_boundary_graph_component_property(const int _id0, const int _id1);


    //used for additional constraints
    //connect interior graph components
    void connect_interior_graph_components();
    void connect_via_dual_path(const std::vector<HFH>& _dpath, const int _axis1, const int _axis2);
	int additional_constraint_type(const std::vector<HFH>& _dpath, std::vector<HEH>& _hes);
    void find_independent_unknown_halfface_transition(const std::vector<HFH>& _dpath, int& _pos);


    //connect boundary graph components
    void connect_boundary_graph_components(std::vector<int>& _chart_index);
    void get_dual_paths_on_trimesh(std::vector<std::vector<TriHEH> >& _tri_dps);
    void get_singular_halfedge_in_halfface(const HFH _hf, HEH& _he);
	void update_boundary_chart_index(std::vector<int>& _chart_index, const TriHEH _he0, const TriHEH _he1)const;


private:
	TetMesh& mesh_;
	//input properties
	EP<int> valance_;
	//intermediate properties
	VP<int> node_type_;
	EP<int> label_;
	VP<bool> is_singular_vt_;

	//output properties
	CP<std::map<HEH, int> > axes_prop_;
	HFP<int> trans_prop_;
	HFP<int> normal_prop_;

	//intermediate properties
	//chart index map to new chart index
	std::map<int, int> real_chart_id_;
	CP<int> chart_id_;

	//tag visited cells
	CP<int> cell_tag_;

	//mark zippered edges
	EP<bool> edge_tag_;
	//singularity graph components
	EP<int> sg_comp_b_, sg_comp_i_;
	//singular arc type:
    //-2: interior arc that ends at boundary
    //-1: interior arc that ends at singular nodes
    //1: interior arc which touches boundary
    //2: interior circle
    //3: boundary circle
    //4: boundary arc that ends at singular nodes
	EP<int> arc_type_;

	//maximum singular curve label
	int max_label_;
	int max_comp_i_;
	int max_comp_b_;

	OVM::StatusAttrib& status_;
//	OVM::ColorAttrib<ACG::Vec4f>& color_;

	//trimesh
	TriMesh trimesh_;
	//trimesh properties
	TriEP<int> tri_val_;
	TriVP<int> tri_node_type_;
	TriVP<int> tri_vh_trans_;
	TriHP<int> tri_trans_;
	TriHP<int> tri_axis_;
	TriHP<HEH> tet_he_;
	std::map<VH, TriVH> vh_to_trivh_;


	std::vector<EH> open_edges_;

	std::vector<std::shared_ptr<CornerT<TetMesh> > > v_corners_;
	std::vector<std::shared_ptr<TangentContinuityT<TetMesh> > > v_tc_;

	TransitionQuaternion tq_;

	SingularityGraphT<TetMesh>& sg_;
};

}

//=============================================================================
#if defined(INCLUDE_TEMPLATES) && !defined(MATCHINGSANDALIGNMENTEXTRACTIONT_C)
#define MATCHINGSANDALIGNMENTEXTRACTIONT_TEMPLATES
#include "MatchingsAndAlignmentExtractionT.cc"
#endif
//=============================================================================


#endif /* MATCHINGSANDALIGNMENTEXTRACTIONT_HH */
