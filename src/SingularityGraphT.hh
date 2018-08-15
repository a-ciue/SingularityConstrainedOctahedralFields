/*
 * SingularityGraphT.hh
 *
 *  Created on: Dec 31, 2017
 *      Author: hliu
 */

#ifndef PLUGIN_SMOOTHFRAMEFIELDOVM_SINGULARITYGRAPHT_HH
#define PLUGIN_SMOOTHFRAMEFIELDOVM_SINGULARITYGRAPHT_HH

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
	  valance_(mesh_.template request_edge_property<int>("edge_valance")),
	  label_(mesh_.template request_edge_property<int>("singularity_line_label")),
	  node_type_(mesh_.template request_vertex_property<int>("singular_node")),
	  is_singular_vt_(mesh_.template request_vertex_property<bool>("is_singular_vertex")),

	  sg_comp_b_(mesh_.template request_edge_property<int>("boundary_singular_graph_component")),
	  sg_comp_i_(mesh_.template request_edge_property<int>("interior_singular_graph_component")),
	  arc_type_(mesh_.template request_edge_property<int>("singular_arc_type"))
    {
    		update_label_property();
    }
    ~SingularityGraphT(){}

public:
    void set_singularity_graph_component_property();
    void classify_singular_arc_type();

    void update_label_property();
    void update_node_type_property();
    void update_singular_vertex_property();

    std::vector<EH> get_singular_edges_of_label(const int _label) const;
    std::vector<VH> get_singular_vertices_of_label(const int _label) const;
    int max_label() const {return max_label_;};

    void sort_vertices_on_curve(const int _label, std::vector<VH>& _vhs) const;
    void sort_edges_on_curve(const int _label, std::vector<EH>& _ehs) const;

    void measure_the_shortest_curve_length() const;

private:
    bool is_end_vertex(const VH _vh) const;
    double distance(const VH _vh0, const VH _vh1) const;

private:
	Mesh &mesh_;
    EP<int> valance_;
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

#endif /* PLUGIN_SMOOTHFRAMEFIELDOVM_SINGULARITYGRAPHT_HH */
