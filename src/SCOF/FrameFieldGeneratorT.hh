//=============================================================================
//
//  CLASS FrameFieldGeneratorT
//
//=============================================================================


#ifndef FRAMEFIELDGENERATORT_HH
#define FRAMEFIELDGENERATORT_HH


//== INCLUDES =================================================================

#include <map>
#include <vector>
#include <queue>
#include <iostream>
#include <iomanip>
#include <float.h>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <OpenVolumeMesh/Mesh/TetrahedralMesh.hh>
#include <OpenVolumeMesh/Core/Properties/PropertyPtr.hh>
#include <OpenVolumeMesh/Mesh/TetrahedralGeometryKernel.hh>

#include "TransitionQuaternionEigen.hh"
#include "AxisAlignment.hh"
 
//== FORWARDDECLARATIONS ======================================================

//== NAMESPACES ===============================================================

namespace SCOF {

//== CLASS DEFINITION =========================================================



	      
/** \class FrameFieldGeneratorT FrameFieldGeneratorT.hh

    Brief Description.
  
    A more elaborate description follows.
*/

template <class TetMeshT>
class FrameFieldGeneratorT
{
public:

  // typedefs for shorter access
  typedef TetMeshT TetMesh;

  typedef OpenVolumeMesh::VertexHandle VH;
  typedef OpenVolumeMesh::EdgeHandle EH;
  typedef OpenVolumeMesh::HalfEdgeHandle HEH;
  typedef OpenVolumeMesh::CellHandle CH;
  typedef OpenVolumeMesh::FaceHandle FH;
  typedef OpenVolumeMesh::HalfFaceHandle HFH;
  typedef OpenVolumeMesh::EdgeIter EIt;
  typedef OpenVolumeMesh::HalfEdgeIter HEIt;
  typedef OpenVolumeMesh::VertexIter VIt;
  typedef OpenVolumeMesh::CellIter CIt;
  typedef OpenVolumeMesh::FaceIter FIt;
  typedef OpenVolumeMesh::BoundaryFaceIter BFIt;
  typedef OpenVolumeMesh::HalfFaceIter HFIt;
  typedef OpenVolumeMesh::HalfEdgeCellIter HECIt;
  typedef OpenVolumeMesh::HalfFaceVertexIter HFVIt;
  typedef OpenVolumeMesh::HalfEdgeHalfFaceIter HEHFIt;
  typedef OpenVolumeMesh::VertexOHalfEdgeIter VOHEIt;
  typedef OpenVolumeMesh::VertexCellIter VCIt;

  typedef typename TetMesh::PointT    Point;
  typedef typename TetMesh::Cell      Cell;
  typedef typename TetMesh::Face      Face;

  // Eigen Types
  typedef Eigen::Triplet<double>      Trip;
  typedef Eigen::SparseMatrix<double> SparseMatrix;
  typedef Eigen::Quaternion<double>   Quaternion;
  typedef Eigen::VectorXd             VectorXd;
  typedef Eigen::Matrix<double, 4, 1> Vec4;
  typedef Eigen::Matrix<double, 3, 1> Vec3;
  typedef Eigen::Matrix<double, 2, 1> Vec2;
  typedef Eigen::Matrix<double, 9, 9> Mat9x9;
  typedef Eigen::Matrix<double, 3, 3> Mat3x3;
  typedef Eigen::Matrix<double, 2, 2> Mat2x2;
  typedef Eigen::Matrix<double, 3, 2> Mat3x2;
  typedef Eigen::Matrix<int, 3, 1>    Vec3i;
  typedef Eigen::Matrix<int, 4, 1>    Vec4i;

//  typedef IntT                             Int;
//  typedef PointT                           Point;
//  typedef TransitionT<IntT, PointT>        Transition;
//  typedef TransitionQuaternion::Matrix     Matrix;

  // constructors
  FrameFieldGeneratorT(TetMesh& _m)
    : mesh_(_m),
      axes_cprop_( mesh_.template request_cell_property<std::map<HEH, int> >("EdgeAxis")),
      trans_hfprop_( mesh_.template request_halfface_property<int>("HalffaceTransiton")),
      normal_hfprop_(mesh_.template request_halfface_property<int>("BoundaryFaceNormal")),
      valance_( mesh_.template request_edge_property<int>("edge_valance")),
      qtrans_hfprop_( mesh_.template request_halfface_property<Quaternion>("HalffaceTransitionQuaternion")),
      quaternion_cprop_( mesh_.template request_cell_property<Quaternion>("FrameFieldQuaternions"))
  {
    mesh_.set_persistent(quaternion_cprop_,true);
  }

  // generate smooth frame field w.r.t. to provided transition functions and partial field alignment
  //  axes_cprop_( mesh_.template request_cell_property<std::map<HEH, int> >("EdgeAxis")),
  //  trans_hfprop_( mesh_.template request_halfface_property<int>("HalffaceTransiton")),
  //  normal_hfprop_(mesh_.template request_halfface_property<int>("BoundaryFaceNormal")),
  //  sgl_eprop_( mesh_.template request_edge_property<int>("edge_valance"))
  void generate_frame_field();

  // check whether the provided input properties are consistent
  void check_input_consistency();

  // save frames to file
  void save_frames(const std::string _filename) const;

private:

  void solve(const std::vector<bool>& face_conquered, std::vector<double>& quaternions);

  void init_consistent_transition_quaternions(std::vector<bool>& _face_conquered);

  void add_smoothness_term(const int _cidx0, const int _cidx1, const Eigen::Quaterniond _trans, std::vector<Trip>& _coeffs) const;

  void add_alignment_coeffs( std::vector<Trip>&     _coeffs,
                             const int              _cidx,
                                   Vec3             _axis,
                             const AxisAlignment    _align,
                             const double           _penalty) const;

  void solve_eigen_problem_inverse_power( const std::vector<Trip>& _coeffs, std::vector<double>& _quaternions) const;

  void solve_unit_length_regularization( const std::vector<Trip>& _coeffs, std::vector<double>& _quaternions) const;

  int number_of_adjacent_boundary_faces(const CH _ch) const;

  int get_boundary_alignment_hfacet(const CH _ch, HFH& _hfh_boundary, AxisAlignment& _axis) const;

  int get_singularity_alignment_hedge(const CH _ch, HEH& _he_singularity, AxisAlignment& _axis) const;

  // collect transitions counter-clockwise around halfedge _he, starting and ending at _ch (or ending at boundary if available)
  int ccw_transition(const CH _ch, const HEH _he) const;

  Eigen::Quaterniond ccw_transition_quaternion(const CH _ch, const HEH _he) const;

  HFH find_adjacent_halfface(const CH _ch, const HEH _heh) const;

  bool is_consistent(const int _transition, const AxisAlignment _axis, const int _valance) const;

  void normalize_rotation(int& _axis, int& _angle) const;

  void adjacent_hehs(const CH _ch, std::vector<HEH>& _adj_heh) const;

  CH get_boundary_start_cell(const HEH _heh) const;

  CH get_boundary_end_cell(const HEH _heh) const;

  Vec3 normal_vector(const HFH _hfh) const;

  double measure_energy( const std::vector<Trip>& _coeffs, std::vector<double>& _quaternions) const;
  double measure_energy_directly( std::vector<double>& _quaternions) const;


  void normalize_quaternions( std::vector<double>& _quaternions) const;

  void project_to_alignment_constraints( std::vector<double>& _quaternions) const;

  void project_to_alignment_constraints( std::vector<double>& _quaternions, const std::vector< std::pair<CH,Eigen::Quaterniond> >& _full_alignment, const std::vector< std::pair<CH, std::pair<AxisAlignment, Vec3> > >& _partial_alignment) const;

  void collect_alignment_conditions( std::vector< std::pair<CH,Eigen::Quaterniond> >& _full_alignment, std::vector< std::pair<CH, std::pair<AxisAlignment, Vec3> > >& _partial_alignment) const;

  Eigen::Quaterniond construct_quaternion(const AxisAlignment _axis0, const AxisAlignment _axis1, const Vec3& _vector0, const Vec3& _vector1) const;

  Eigen::Quaterniond project_quaternion(const Eigen::Quaterniond& _q, const AxisAlignment _axis_label, const Vec3& _axis_vector) const;

  Eigen::Quaterniond get_quaternion(const unsigned int _i, const std::vector<double>& _quaternions) const;

  void set_quaternion(const unsigned int _i, const Eigen::Quaterniond& _q, std::vector<double>& _quaternions) const;

  void init_quaternions_by_propagation(std::vector<double>& _quaternions) const;

//  Eigen::Quaterniond acg2eigen(const QuaternionACG& _q) const {return Eigen::Quaterniond(_q[0], _q[1], _q[2], _q[3]);}
//  QuaternionACG      eigen2acg(const Eigen::Quaterniond& _q) const {return QuaternionACG(_q.w(), _q.x(), _q.y(), _q.z());}

  int n_open_facets( const EH _eh, const std::vector<bool>& _face_conquered, FH& _open_fh) const;
  int n_open_facets( const HEH _he, const std::vector<bool>& _face_conquered, HFH& _open_hf) const;


  HEH find_halfedge_of_halfface(const HFH _hfh, const EH _eh) const;

  Eigen::Quaterniond negate( const Eigen::Quaterniond& _q) const
  {
    return Eigen::Quaterniond(-_q.w(),-_q.x(),-_q.y(),-_q.z());
  }

  void test_loop_conditions() const;

private:
   
  TetMesh& mesh_;

  // #### input
  OpenVolumeMesh::CellPropertyT<std::map<HEH, int> > axes_cprop_;
  OpenVolumeMesh::HalfFacePropertyT<int>             trans_hfprop_;
  OpenVolumeMesh::HalfFacePropertyT<int>             normal_hfprop_;
  OpenVolumeMesh::EdgePropertyT<int>                 valance_;

  // #### intermediate
  OpenVolumeMesh::HalfFacePropertyT<Eigen::Quaterniond> qtrans_hfprop_;

  // #### output
  // field quaternions
  OpenVolumeMesh::CellPropertyT< Eigen::Quaterniond> quaternion_cprop_;

  TransitionQuaternion tq_;
};



//=============================================================================
} // namespace SCOF
//=============================================================================
#if defined(INCLUDE_TEMPLATES) && !defined(FRAMEFIELDGENERATORT_C)
#define FRAMEFIELDGENERATORT_TEMPLATES
#include "FrameFieldGeneratorT.cc"
#endif
//=============================================================================
#endif // FRAMEFIELDGENERATORT_HH defined
//=============================================================================

