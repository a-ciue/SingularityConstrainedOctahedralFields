/*
 * SingularityGraphT.hh
 *
 *  Created on: Dec 31, 2017
 *      Author: hliu
 */

#ifndef SINGULARITYGRAPHT_HH
#define SINGULARITYGRAPHT_HH

#include "Typedefs.hh"

namespace OVM = OpenVolumeMesh;

namespace SCOF
{

template <class MeshT>
class SingularityGraphT
{
public:
	using Mesh = MeshT;

    SingularityGraphT(Mesh& _mesh)
    : mesh_(_mesh),
	  valence_(mesh_.template request_edge_property<int>("edge_valance")),
	  label_(mesh_.template request_edge_property<int>("singularity_line_label")),
	  node_type_(mesh_.template request_vertex_property<int>("singular_node")),
	  is_singular_vt_(mesh_.template request_vertex_property<bool>("is_singular_vertex")),

	  sg_comp_b_(mesh_.template request_edge_property<int>("boundary_singular_graph_component")),
	  sg_comp_i_(mesh_.template request_edge_property<int>("interior_singular_graph_component")),
	  arc_type_(mesh_.template request_edge_property<int>("singular_arc_type"))
    {
    		update_label_property();
    		mesh_.set_persistent(arc_type_, true);
    }
    ~SingularityGraphT(){}

public:
	void set_singularity_graph_component_property();
	void classify_singular_arc_type();

	void update_label_property();
	void update_node_type_property();
	void update_singular_vertex_property();
	void update_valence_property();

	std::vector<EH> get_singular_edges_of_label(const int _label) const;
	std::vector<VH> get_singular_vertices_of_label(const int _label) const;

//	const std::vector<EH>& get_singular_edges_of_label(const int _label) const;
//	const std::vector<VH>& get_singular_vertices_of_label(const int _label) const;

	int max_label() {max_label_ = *std::max_element(label_.begin(), label_.end()); return max_label_;};

	void sort_vertices_on_curve(const int _label, std::vector<VH>& _vhs) const;
	void sort_edges_on_curve(const int _label, std::vector<EH>& _ehs) const;
    void sort_halfedges_on_curve(const int _label, std::vector<HEH>& _hes) const;

    void measure_the_shortest_curve_length() const;

	int arc_type(const int _label) const{ return arc_type_[get_singular_edges_of_label(_label)[0]]; }
	int arc_type(const EH _eh) const {return arc_type_[_eh];}

	int valence(const int _label) const { return valence_[get_singular_edges_of_label(_label)[0]]; }
	int valence(const EH _eh) const { return valence_[_eh]; }

	std::vector<VH> adjacent_singular_nodes(const VH _vh) const;

	VH end_singular_node_on_curve(const int _label, const VH _start_node) const;

    //Heuristically get the node index by counting edge valence
    //Valid singular node type with <= 4 incident singular edges

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
    //-12: (0, 1, 2, 1, 1, 0)
    //-13: (0, 2, 2, 1, 0, 1)
    //-14: (1, 3, 2, 0, 0, 1)
	int node_index(const VH& _vh) const;

private:
    bool is_end_vertex(const VH _vh) const;
    double distance(const VH _vh0, const VH _vh1) const;

private:
	Mesh &mesh_;
    EP<int> valence_;
    EP<int> label_;
    VP<int> node_type_;
	VP<bool> is_singular_vt_;

    //singular graph components
    EP<int> sg_comp_b_, sg_comp_i_;
	//singular arc type:
    //-2: interior arc that ends at boundary
    //-1: interior arc that ends at singular nodes
    //1: interior arc which touches boundary on one end
    //2: interior circle
    //3: boundary circle
    //4: boundary arc that ends at singular nodes
    EP<int> arc_type_;

	//maximum singular curve label
	int max_label_;

	std::vector<std::vector<EH> > v_sg_ehs_;
	std::vector<std::vector<VH> > v_sg_vhs_;
};

}
//=============================================================================
#if defined(INCLUDE_TEMPLATES) && !defined(SINGULARITYGRAPHT_C)
#define SINGULARITYGRAPH_TEMPLATES
#include "SingularityGraphT.cc"
#endif
//=============================================================================

#endif /* SINGULARITYGRAPHT_HH */
