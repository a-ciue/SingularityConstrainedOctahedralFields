//=============================================================================
//
//  CLASS FrameFieldGeneratorT - IMPLEMENTATION
//
//=============================================================================

#define FRAMEFIELDGENERATORT_C

//== INCLUDES =================================================================

#include "FrameFieldGeneratorT.hh"
//#include <CoMISo/Config/config.hh>
//
#ifdef ENABLE_SUITESPARSE
  #include <Eigen/CholmodSupport>
#endif

#include "StopWatch.hh"
//== NAMESPACES ===============================================================

namespace SCOF {

//== IMPLEMENTATION ==========================================================

template <class TetMeshT>
void
FrameFieldGeneratorT<TetMeshT>::
generate_frame_field()
{
//  test_loop_conditions();
//  return;

  // 0. check input consistency
  check_input_consistency();
    StopWatch<> time;
    time.start();
  // quaternions are a double-cover of frames, hence the sign has to be chosen carefully in order to obtain a quadratic smoothness energy
  std::vector<bool> face_conquered(mesh_.n_faces(), false);
  init_consistent_transition_quaternions(face_conquered);

  //2. solve
  // init unknown quaternions with identity
  std::vector<double> quaternions(4*mesh_.n_cells(),0.0);
  for(unsigned int i=0; i<mesh_.n_cells(); ++i)
    quaternions[4*i] = 1.0;
  solve(face_conquered, quaternions);

  // 3. store result
  OpenVolumeMesh::CellPropertyT< Quaternion > quaternion_prop = mesh_.template request_cell_property< Quaternion>("CellQuaternion");
  for(CIt c_it = mesh_.c_iter(); c_it.valid(); ++c_it)
  {
    int bidx = 4*c_it->idx();
    Quaternion qc(quaternions[bidx], quaternions[bidx+1], quaternions[bidx+2], quaternions[bidx+3]);
    quaternion_prop[*c_it] = qc;
  }
//
//  //handle cut faces
//  bool has_handle = false;
//  for(auto fhi : mesh_.faces())
//	  if(!mesh_.is_boundary(fhi) && !face_conquered[fhi.idx()])
//	  {
//		  has_handle = true;
//		  break;
//	  }
//  if(has_handle)
//  {
//	  std::cerr<<"\nfaces: ";
//
//	  for(auto fhi : mesh_.faces())
//	  {
//		  if(!mesh_.is_boundary(fhi) && !face_conquered[fhi.idx()])
//		  {
//			  std::cerr<<" "<<fhi;
//		  }
//	  }
//	  //7. handle the cut faces
//	  for(auto fhi : mesh_.faces())
//	  {
//		  if(!mesh_.is_boundary(fhi) && !face_conquered[fhi.idx()])
//		  {
//			  std::vector<HFH> cut_hfs;
//			  HFH hf0 = mesh_.halfface_handle(fhi, 0);
//			  cut_hfs.push_back(hf0);
//			  face_conquered[fhi.idx()] = true;
//
//			  //zipping
//			  // queue on mesh candidate edges (with a single open facet)
//			  std::queue<HEH> heq;
//			  std::vector<HEH> hes = mesh_.halfface(hf0).halfedges();
//			  for(auto hei : hes)
//				  heq.push(mesh_.opposite_halfedge_handle(hei));
//std::cerr<<"\ncells: "<<mesh_.incident_cell(hf0);
//			  while(!heq.empty())
//			  {
//			    // get next candidate
//			    HEH he = heq.front();
//			    heq.pop();
//
//			    if(mesh_.is_boundary(he))
//			    	continue;
//			    HFH hf = HFH(-1);
//
//			    // still a valid candidate?
//			    if(n_open_facets(he, face_conquered, hf) == 1)
//			    {
//			    	cut_hfs.push_back(hf);
////			        Eigen::Quaterniond q = ccw_transition_quaternion(mesh_.incident_cell(hf), he);
////
////			        // check sign and correct if necessary
////			        if(q.w() < 0)
////			        {
////			          HFH hfho = mesh_.opposite_halfface_handle(hf);
////			          qtrans_hfprop_[hf ] = negate(qtrans_hfprop_[hf]);
////			          qtrans_hfprop_[hfho] = negate(qtrans_hfprop_[hfho]);
////			        }
//			      if(!hf.is_valid())
//			    	  	  std::cerr << "Warning: could not find_halfedge_of_halfface, which should never happen here!!!" << std::endl;
//			      // mark as visited
//			      face_conquered[mesh_.face_handle(hf).idx()] = true;
//			      std::cerr<<" "<<mesh_.incident_cell(hf);
//
//			      // add potential new candidates
//			      for(unsigned int i=0; i< mesh_.halfface(hf).halfedges().size(); ++i)
//			        if(mesh_.halfface(hf).halfedges()[i] != he)
//			          heq.push(mesh_.opposite_halfedge_handle(mesh_.halfface(hf).halfedges()[i]));
//			    }
//			  }
//
//
//			  // measure local energy
//			  double local_energy0 = 0.0,
//					  local_energy1 = 0.0;
//			  for(auto hfi : cut_hfs)
//			  {
//				  auto ch1 = mesh_.incident_cell(hfi);
//				  auto ch0 = mesh_.incident_cell(mesh_.opposite_halfface_handle(hfi));
//				  QuaternionACG qi = tq_.transition(trans_hfprop_[hfi]);
//				  QuaternionACG q_diff = qi*quaternion_acg_prop[ch0] - quaternion_acg_prop[ch1];
//				  local_energy0 = q_diff.norm();
//
//				  QuaternionACG qi_ng = eigen2acg(negate(acg2eigen(qi)));
//				  q_diff = qi_ng*quaternion_acg_prop[ch0] - quaternion_acg_prop[ch1];
//				  local_energy1 = q_diff.norm();
//
//				  if(local_energy0 < local_energy1)
//				  {
//					  qtrans_hfprop_[hfi] = acg2eigen(qi);
//					  qtrans_hfprop_[mesh_.opposite_halfface_handle(hfi)] = acg2eigen(qi).inverse();
//				  }else
//				  {
//					  qtrans_hfprop_[hfi] = acg2eigen(qi_ng);
//					  qtrans_hfprop_[mesh_.opposite_halfface_handle(hfi)] = acg2eigen(qi_ng).inverse();
//				  }
//			  }
//
////			  if(local_energy0 < local_energy1)
////			  {
////				  for(auto hfi : cut_hfs)
////				  {
////					  Quaternion qt = acg2eigen(tq_.transition(trans_hfprop_[hfi]));
////					  qtrans_hfprop_[hfi] = qt;
////					  qtrans_hfprop_[mesh_.opposite_halfface_handle(hfi)] = qt.inverse();
////				  }
////			  }else
////			  {
////				  for(auto hfi : cut_hfs)
////				  {
////					  Quaternion qt = acg2eigen(tq_.transition(trans_hfprop_[hfi]));
////					  Quaternion qt_inverse = qt.inverse();
////					  qt = negate(qt);
////					  qtrans_hfprop_[hfi] = qt;
////					  qtrans_hfprop_[mesh_.opposite_halfface_handle(hfi)] = negate(qt_inverse);
////				  }
////			  }
//
//		  }
//	  }
//	  // check if zipping reached all faces
//	  for(unsigned int i=0; i<face_conquered.size(); ++i)
//	    if(!mesh_.is_boundary(FH(i)) && face_conquered[i] == false)
//	      std::cerr << "Warning: face " << i << " was not conquered during zipping!!!" << std::endl;
//
//	  //slove again
//	  	  for(unsigned int i=0; i<mesh_.n_cells(); ++i)
//	  		quaternions[i] = 0.0;
//	    for(unsigned int i=0; i<mesh_.n_cells(); ++i)
//	      quaternions[4*i] = 1.0;
//	    solve(face_conquered, quaternions);
//
//
//	    for(CIt c_it = mesh_.c_iter(); c_it.valid(); ++c_it)
//	    {
//	      int bidx = 4*c_it->idx();
//	      QuaternionACG qc(quaternions[bidx], quaternions[bidx+1], quaternions[bidx+2], quaternions[bidx+3]);
//	      quaternion_acg_prop[*c_it] = qc;
//	    }
//  }

  //  OpenVolumeMesh::CellPropertyT< Eigen::Quaterniond > quaternion_prop = mesh_.template request_cell_property< Eigen::Quaterniond >("CellQuaternion");
  OpenVolumeMesh::CellPropertyT< Point >         frame_x_prop = mesh_.template request_cell_property< Point >("CellFrameX");
  OpenVolumeMesh::CellPropertyT< Point >         frame_y_prop = mesh_.template request_cell_property< Point >("CellFrameY");
  OpenVolumeMesh::CellPropertyT< Point >         frame_z_prop = mesh_.template request_cell_property< Point >("CellFrameZ");

  // make properties persistent
//  mesh_.set_persistent(quaternion_prop);
  mesh_.set_persistent(frame_x_prop);
  mesh_.set_persistent(frame_y_prop);
  mesh_.set_persistent(frame_z_prop);

  for(CIt c_it = mesh_.c_iter(); c_it.valid(); ++c_it)
  {
    int bidx = 4*c_it->idx();
    Quaternion qc(quaternions[bidx], quaternions[bidx+1], quaternions[bidx+2], quaternions[bidx+3]);
    quaternion_prop[*c_it] = qc;
//    quaternion_prop[*c_it] = acg2eigen(qc);
    Quaternion qc2; qc2 = qc;
    qc.normalize();
    if(!std::isfinite(qc.squaredNorm()))
    {
      std::cerr << "Warning: found degenerate quaternion" << std::endl;
      qc = Quaternion(1,0,0,0);
    }

    // create rotation matrices and extract frame axes
    Mat3x3 R = qc.toRotationMatrix();
    frame_x_prop[*c_it] = qc2.norm()*Point(R(0,0), R(1,0), R(2,0));
    frame_y_prop[*c_it] = qc2.norm()*Point(R(0,1), R(1,1), R(2,1));
    frame_z_prop[*c_it] = qc2.norm()*Point(R(0,2), R(1,2), R(2,2));

    // store field quaternion for parametrization
    quaternion_cprop_[*c_it] = qc;
  }
    std::cerr << "\n#####Field Generation time:    " <<time.stop()/1000<<" s";

  save_frames("test.frames");
}


//-----------------------------------------------------------------------------

template <class TetMeshT>
void
FrameFieldGeneratorT<TetMeshT>::
solve(const std::vector<bool>& face_conquered, std::vector<double>& quaternions)
{
	  // 1. set constants
	  double penalty_boundary_alignment = 10.0;
	  double penalty_singularity_alignment = 10.0;

	 // 2. allocate system matrix
	  int n = 4*mesh_.n_cells();
	  std::vector<Trip> A_coeff; A_coeff.reserve(10*n);

	  // 3. setup smoothness energy
	  for(FIt f_it = mesh_.f_iter(); f_it.valid(); ++f_it)
	    if(!mesh_.is_boundary(*f_it) && face_conquered[*f_it])
	    {
	      // get halffaces and cells
	      HFH hfh0 = mesh_.halfface_handle(*f_it, 0);
	      HFH hfh1 = mesh_.halfface_handle(*f_it, 1);
	      CH  ch0  = mesh_.incident_cell(hfh0);
	      CH  ch1  = mesh_.incident_cell(hfh1);

	      // get transition from c0 -> c1
	      Eigen::Quaterniond trans = qtrans_hfprop_[hfh0];

	      add_smoothness_term(ch0.idx(), ch1.idx(), trans, A_coeff);
	    }

	  // 4. setup alignment penalties
	  int nb_complete(0);
	  int na_complete(0);
	  for(CIt c_it = mesh_.c_iter(); c_it.valid(); ++c_it)
	  {
	    //      std::cerr << "#cell alignment conditions: " << axes_cprop_[*c_it].size() << std::endl;

	    // a) align boundary facets
	    // check if boundary alignment facet exists
	    HFH hfh_boundary;
	    AxisAlignment axis;
	    int nba = get_boundary_alignment_hfacet(*c_it, hfh_boundary, axis);

	    if(nba > 1)
	      std::cerr << "ERROR: tetmesh has tetrahedron with " << nba << " boundary half-faces!!!" << std::endl;

	    if(nba == 1)
	    {
	      add_alignment_coeffs( A_coeff,
	          c_it->idx(),
	          normal_vector(hfh_boundary),
	          axis,
	          penalty_boundary_alignment);
	      ++nb_complete;
	    }

	    // b) align singular axes
	    HEH he_singularity;
	    int nsa = get_singularity_alignment_hedge(*c_it, he_singularity, axis);

	    if(nsa > 1)
	    {
	      std::cerr << "ERROR: invalid number of singular edges adjacent at tetrahedron - " << nsa << std::endl;
	    }

	//    if(0)
	    if(nsa == 1) // check singular edge
	    {
	      Point p0 = mesh_.vertex( mesh_.halfedge(he_singularity).from_vertex());
	      Point p1 = mesh_.vertex( mesh_.halfedge(he_singularity).to_vertex());
	      Point e  = p1-p0;
	      e.normalize();
	      if(std::isfinite(e.sqrnorm()))
	      {
	        Vec3 ee(e[0], e[1], e[2]);
	        add_alignment_coeffs( A_coeff,
	            c_it->idx(),
	            ee,
	            axis,
	            penalty_singularity_alignment);
	        ++na_complete;
	      }
	      else std::cerr << "ERROR: input mesh has degenerate edge!!!" << std::endl;
	    }
	  }
	  std::cerr << "#boundary    alignment constraints: " << nb_complete << std::endl;
	  std::cerr << "#singularity alignment constraints: " << na_complete << std::endl;

	  // 5a. solve eigenproblem


	//  // initialize by propagation
	//  init_quaternions_by_propagation(quaternions);
	//  std::cerr << "energy after propagation         : " << measure_energy(A_coeff, quaternions) << std::endl;
	//  std::cerr << "energy after propagation directly: " << measure_energy_directly(quaternions) << std::endl;


	  // regularize
	  std::vector<Trip> A_coeff_reg;
	  A_coeff_reg = A_coeff;
	  const double regularization_eps = 1e-6;
	  for( int i=0; i<n; ++i)
	    A_coeff_reg.push_back(Trip(i,i,regularization_eps));

	  // solve eigenvalue problem
	  solve_eigen_problem_inverse_power( A_coeff_reg, quaternions);

	  std::cerr << "energy after eigen problem: " << measure_energy(A_coeff, quaternions) << std::endl;
	  // correct alignment (due to penalty deviation)
	  normalize_quaternions(quaternions);
	  std::cerr << "energy after normalization: " << measure_energy(A_coeff, quaternions) << std::endl;
	  project_to_alignment_constraints(quaternions);
	  std::cerr << "energy after constraint projection: " << measure_energy(A_coeff, quaternions) << std::endl;

	  // 5b. solve with unit-length regularization
	  solve_unit_length_regularization( A_coeff, quaternions);
	  std::cerr << "energy after unit-lengt opt: " << measure_energy(A_coeff, quaternions) << std::endl;
	  normalize_quaternions(quaternions);
	  std::cerr << "energy after normalization: " << measure_energy(A_coeff, quaternions) << std::endl;
}


//-----------------------------------------------------------------------------

template <class TetMeshT>
void
FrameFieldGeneratorT<TetMeshT>::
init_consistent_transition_quaternions(std::vector<bool>& _face_conquered)
{
  // 1. init transition quaternions
  for(FIt f_it = mesh_.f_iter(); f_it.valid(); ++f_it)
  {
    // get halffaces and cells
    HFH hfh0 = mesh_.halfface_handle(*f_it, 0);
    HFH hfh1 = mesh_.halfface_handle(*f_it, 1);

    Eigen::Quaterniond qt(1,0,0,0);

    if(!mesh_.is_boundary(*f_it))
    {
      // get transition from c0 -> c1
      int trans = trans_hfprop_[hfh0];

      qt = tq_.transition(trans);
    }

    qtrans_hfprop_[hfh0] = qt;
    qtrans_hfprop_[hfh1] = qt.inverse();
  }

  // 2. perform zipping to correct transition signs in one-rings around interior edges

  // 2a. spanning tree
  // spanning tree of cells
  std::vector<bool> cell_visited(mesh_.n_cells(),false);

  // qeueue
  std::queue<HFH> cq;

  // select arbitrary root
  CH ch(0);
  cell_visited[ch.idx()] = true;
  // push hf to neighbors
  for(unsigned int i=0; i<mesh_.cell(ch).halffaces().size(); ++i)
    cq.push(mesh_.cell(ch).halffaces()[i]);

  while(!cq.empty())
  {
    //get next halfface
    HFH hfh = cq.front();
    cq.pop();
    HFH hfh_opp = mesh_.opposite_halfface_handle(hfh);
    if(!mesh_.is_boundary(hfh_opp))
    {
      CH  ch     = mesh_.incident_cell(hfh);
      CH  ch_opp = mesh_.incident_cell(hfh_opp);

      if(!cell_visited[ch_opp.idx()])
      {
        // mark as visited
        cell_visited[ch_opp.idx()] = true;
        // mark face conquered
        _face_conquered[mesh_.face_handle(hfh)] = true;
        // push halffaces to neighbors
        for(unsigned int i=0; i<mesh_.cell(ch_opp).halffaces().size(); ++i)
          cq.push(mesh_.cell(ch_opp).halffaces()[i]);
      }
    }
  }

  // 2b. zipping
  // queue on mesh candidate edges (with a single open facet)
  std::queue<EH> eq;
  FH fh;
  for(EIt e_it = mesh_.e_iter(); e_it.valid(); ++e_it)
    if(!mesh_.is_boundary(*e_it))
      if(n_open_facets(*e_it, _face_conquered, fh) == 1)
        eq.push(*e_it);

  while(!eq.empty())
  {
    // get next candidate
    EH eh = eq.front();
    eq.pop();
    // still a valid candidate?
    if(n_open_facets(eh,_face_conquered, fh) == 1)
    {
      // get halfface and corresponding halfedge
      HFH hfh = mesh_.halfface_handle(fh,0);
      HEH heh = find_halfedge_of_halfface(hfh, eh);
      if(!heh.is_valid()) std::cerr << "Warning: could not find_halfedge_of_halfface, which should never happen here!!!" << std::endl;
      else
      {
        Eigen::Quaterniond q = ccw_transition_quaternion(mesh_.incident_cell(hfh), heh);

        // check sign and correct if necessary
        if(q.w() < 0)
        {
          HFH hfho = mesh_.opposite_halfface_handle(hfh);
          qtrans_hfprop_[hfh ] = negate(qtrans_hfprop_[hfh ]);
          qtrans_hfprop_[hfho] = negate(qtrans_hfprop_[hfho]);
        }

        // check result
        {
          q = ccw_transition_quaternion(mesh_.incident_cell(hfh), heh);
          if(q.w() < 0)
            std::cerr << "Warning: sign correction failed!!!" << std::endl;
        }
      }
      // mark as visited
      _face_conquered[fh.idx()] = true;

      // add potential new candidates
      for(unsigned int i=0; i< mesh_.halfface(hfh).halfedges().size(); ++i)
        if(mesh_.halfface(hfh).halfedges()[i] != heh)
          eq.push(mesh_.edge_handle(mesh_.halfface(hfh).halfedges()[i]));
    }
  }

//  //deal with cut faces
//  for(auto fhi : mesh_.faces())
//	  if(!mesh_.is_boundary(fhi) && !_face_conquered[fhi.idx()])
//    	{
//    		HFH hfh = mesh_.halfface_handle(fhi, 0);
//    		HFH hfho = mesh_.opposite_halfface_handle(hfh);
////    		qtrans_hfprop_[hfh ] = negate(qtrans_hfprop_[hfh ]);
////    		qtrans_hfprop_[hfho] = negate(qtrans_hfprop_[hfho]);
//    		std::cerr<<"\nfhi: "<<fhi;
//    		_face_conquered[fhi.idx()] = true;
//    		for(unsigned int i=0; i< mesh_.halfface(hfh).halfedges().size(); ++i)
//    		{
//    			EH ehi = mesh_.edge_handle(mesh_.halfface(hfh).halfedges()[i]);
//             eq.push(ehi);
//        	 	std::cerr<<"\nehi: "<<ehi;
//         }
//
//    		 while(!eq.empty())
//    		  {
//    		    // get next candidate
//    		    EH eh = eq.front();
//    		    eq.pop();
//
//    		    // still a valid candidate?
//    		    if(n_open_facets(eh,_face_conquered, fh) == 1)
//    		    {
//    		      // get halfface and corresponding halfedge
//    		      HFH hfh = mesh_.halfface_handle(fh,0);
//    		      HEH heh = find_halfedge_of_halfface(hfh, eh);
//    		      if(!heh.is_valid()) std::cerr << "Warning: could not find_halfedge_of_halfface, which should never happen here!!!" << std::endl;
//    		      else
//    		      {
//    		        Eigen::Quaterniond q = ccw_transition_quaternion(mesh_.incident_cell(hfh), heh);
//
//    		        // check sign and correct if necessary
//    		        if(q.w() < 0)
//    		        {
//    		          HFH hfho = mesh_.opposite_halfface_handle(hfh);
//    		          qtrans_hfprop_[hfh ] = negate(qtrans_hfprop_[hfh ]);
//    		          qtrans_hfprop_[hfho] = negate(qtrans_hfprop_[hfho]);
//    		        }
//
//    		        // check result
//    		        {
//    		          q = ccw_transition_quaternion(mesh_.incident_cell(hfh), heh);
//    		          if(q.w() < 0)
//    		            std::cerr << "Warning: sign correction failed!!!" << std::endl;
//    		        }
//    		      }
//    		      // mark as visited
//    		      _face_conquered[fh.idx()] = true;
//
//    		      // add potential new candidates
//    		      for(unsigned int i=0; i< mesh_.halfface(hfh).halfedges().size(); ++i)
//    		        if(mesh_.halfface(hfh).halfedges()[i] != heh)
//    		          eq.push(mesh_.edge_handle(mesh_.halfface(hfh).halfedges()[i]));
//    		    }
//    		  }
//    	}
}

//-----------------------------------------------------------------------------

template <class TetMeshT>
typename FrameFieldGeneratorT<TetMeshT>::HEH
FrameFieldGeneratorT<TetMeshT>::
find_halfedge_of_halfface(const HFH _hfh, const EH _eh) const
{
  for(unsigned int i=0; i<mesh_.halfface(_hfh).halfedges().size(); ++i)
  {
    HEH heh = mesh_.halfface(_hfh).halfedges()[i];
    if(mesh_.edge_handle(heh) == _eh)
      return heh;
  }

  return HEH(-1);
}

//-----------------------------------------------------------------------------

template <class TetMeshT>
int
FrameFieldGeneratorT<TetMeshT>::
n_open_facets( const EH _eh, const std::vector<bool>& _face_conquered, FH& _open_fh) const
{
  int n(0);
  // iterate over all halffaces adjacent to edge
  for( HEHFIt hehf_it = mesh_.hehf_iter(mesh_.halfedge_handle(_eh,0)); hehf_it.valid(); ++hehf_it)
  {
    FH fh = mesh_.face_handle(*hehf_it);
    if(_face_conquered[fh.idx()] == false)
    {
      _open_fh = mesh_.face_handle(*hehf_it);
      ++n;
    }
  }
  return n;
}

//-----------------------------------------------------------------------------
template <class TetMeshT>
int
FrameFieldGeneratorT<TetMeshT>::
n_open_facets( const HEH _he, const std::vector<bool>& _face_conquered, HFH& _open_hf) const
{
	  int n(0);
	  // iterate over all halffaces adjacent to edge
	  for( HEHFIt hehf_it = mesh_.hehf_iter(_he); hehf_it.valid(); ++hehf_it)
	  {
	    FH fh = mesh_.face_handle(*hehf_it);
	    if(_face_conquered[fh.idx()] == false)
	    {
	    	_open_hf = *hehf_it;
	      ++n;
	    }
	  }
	  return n;
}

//-----------------------------------------------------------------------------

template <class TetMeshT>
void
FrameFieldGeneratorT<TetMeshT>::
check_input_consistency()
{
  std::cerr << "check input consistency..." << std::endl;

  for(CIt c_it = mesh_.c_iter(); c_it.valid(); ++c_it)
  {
    // 1. check boundary alignment tags ########################################
    int nb = number_of_adjacent_boundary_faces(*c_it);
    if(nb > 1) std::cerr << "ERROR: input has tetrahedron with multiple boundary facets" << std::endl;

    // check boundary alignment facet
    HFH hfh_boundary;
    AxisAlignment axis;
    int nba = get_boundary_alignment_hfacet(*c_it, hfh_boundary, axis);
    if(nb == 1)
    {
      if(nba == 0) std::cerr << "ERROR: tetrahedron is adjacent to boundary but does not have alignment constraint!!!" << std::endl;
      if(nba >= 1 && !mesh_.is_boundary(hfh_boundary)) std::cerr << "ERROR: get_boundary_alignment_facet returned wrong half-face-handle" << std::endl;
    }
    if(nb == 0 && nba > 0)
      std::cerr << "ERROR: a cell that is not adjacent to the boundary has a boundary alignment axis!!!" << std::endl;

    // 2. check singularity alignment tags ########################################
    HEH he_singularity;
    int nsa = get_singularity_alignment_hedge(*c_it, he_singularity, axis);
    if(nsa > 1)
    {
      std::cerr << "ERROR: invalid number of singular edges adjacent at tetrahedron - " << nsa << std::endl;
      for(std::map<HEH,int>::iterator m_it = axes_cprop_[*c_it].begin(); m_it != axes_cprop_[*c_it].end(); ++m_it)
        std::cerr << m_it->first.idx() << " --> " << m_it->second << std::endl;
    }
    if(nsa == 1) // check singular edge
    {
      if(!mesh_.is_boundary(he_singularity)) // internal singularity --> boundary edges are checked later
      {
        int trans = ccw_transition(*c_it, he_singularity);
        if(!is_consistent(trans, axis, valance_[mesh_.edge_handle(he_singularity)]))
          std::cerr << "ERROR: ccw_transition is not consistent with singularity type" << std::endl;
      }
    }

    // check remaining edges of cell
    std::vector<HEH> adj_heh;
    int n_valid_edges(0);
    adjacent_hehs(*c_it, adj_heh);
    std::vector<HEH>::iterator h_it = adj_heh.begin();
    for(; h_it != adj_heh.end(); ++h_it)
      if(!mesh_.is_boundary(*h_it)) // check only internal --> boundary edges are checked separately
      {
        if(valance_[mesh_.edge_handle(*h_it)] == 0) // check only regular edges
        {
          int trans = ccw_transition(*c_it, *h_it);
          if( trans != 0) // check for identity
            std::cerr << "ERROR: regular edge has transition " << trans << std::endl;
          else
            ++n_valid_edges;
        }
      }
  }

  // check boundary edges
  for(EIt e_it = mesh_.e_iter(); e_it.valid(); ++e_it)
    if(mesh_.is_boundary(*e_it))
    {
      // choose halfedge with correct orientation
      HEH heh = mesh_.halfedge_handle(*e_it,0);

      // get edge type
      int eidx = valance_[mesh_.edge_handle(heh)];

      // get start cell
      CH ch0 = get_boundary_start_cell(heh);
      CH ch1 = get_boundary_end_cell(heh);

      // correct orientation if necessary
      if(eidx != 0)
      {
        HEH he_singularity;
        AxisAlignment axiss;
        int nsa = get_singularity_alignment_hedge(ch0, he_singularity, axiss);

        if(he_singularity.is_valid() && heh != he_singularity)
        {
          heh = mesh_.opposite_halfedge_handle(heh);
          std::swap(ch0, ch1);
        }
      }

      HFH hfh_boundary0, hfh_boundary1;
      AxisAlignment axisb0, axisb1;
      int nba0 = get_boundary_alignment_hfacet(ch0, hfh_boundary0, axisb0);
      int nba1 = get_boundary_alignment_hfacet(ch1, hfh_boundary1, axisb1);

      HEH he_singularity0, he_singularity1;
      AxisAlignment axiss0, axiss1;
      int nsa0 = get_singularity_alignment_hedge(ch0, he_singularity0, axiss0);
      int nsa1 = get_singularity_alignment_hedge(ch1, he_singularity1, axiss1);

      // start singular edge matches target edge?
      if(he_singularity0.is_valid() && mesh_.edge_handle(he_singularity0) != mesh_.edge_handle(heh))
      {
        he_singularity0.reset();
        axiss0 = Align_NONE;
      }

      // start singular edge matches target edge?
      if(he_singularity1.is_valid() && mesh_.edge_handle(he_singularity1) != mesh_.edge_handle(heh))
      {
        he_singularity1.reset();
        axiss1 = Align_NONE;
      }

      if(he_singularity0.is_valid() && he_singularity1.is_valid() && he_singularity0 != he_singularity1)
        std::cerr << "Warning: start and end cell have different alignment edges!!! " << he_singularity0.idx() << " " << he_singularity1.idx() << std::endl;

      // get transition
      int trans = ccw_transition(ch0, heh);

      // check equation
      // get rotation matrix
      Quaternion q = tq_.transition(trans);
      Mat3x3 T = q.toRotationMatrix();

      Mat3x2 F0, F1;
      F0.col(0) = AxisAlignmentHelpers::vector(axisb0);
      F0.col(1) = AxisAlignmentHelpers::vector(axiss0);

      F1.col(0) = AxisAlignmentHelpers::vector(axisb1);
      F1.col(1) = AxisAlignmentHelpers::vector(axiss1);

      Mat3x3 R = Eigen::AngleAxisd( 0.5*M_PI*eidx, F1.col(1)).toRotationMatrix();

      if( (T*F0-R*F1).squaredNorm() > 1e-6)
      {
        std::cerr << "ERROR: found invalid boundary edge transition:" << std::endl;
        std::cerr << "edge_hex_index: " << eidx << std::endl;
        std::cerr << "transition: " << trans << std::endl;
        std::cerr << "start frame normal axis " << axisb0 << " and edge axis " << axiss0 << std::endl;
        std::cerr << "end   frame normal axis " << axisb1 << " and edge axis " << axiss1 << std::endl;
      }
    }

  std::cerr << "done!" << std::endl;
}


//-----------------------------------------------------------------------------

template <class TetMeshT>
void
FrameFieldGeneratorT<TetMeshT>::
add_smoothness_term(const int _cidx0, const int _cidx1, const Eigen::Quaterniond _transition, std::vector<Trip>& _coeffs) const
{
  // the smoothness term is 1/2*||q0*t - q1||^2

  // get transition
  Quaternion t = _transition;

  // setup coefficient matrix corresponding to multiplication with transition function
  // --> 1/2*||A*q0 - q1||^2
  Eigen::Matrix4d A;
  A <<  t.w(), -t.x(), -t.y(), -t.z(),
        t.x(),  t.w(),  t.z(), -t.y(),
        t.y(), -t.z(),  t.w(),  t.x(),
        t.z(),  t.y(), -t.x(),  t.w();

  // get 8 indices
  Vec4i i0;
  i0 << 4*_cidx0, 4*_cidx0+1, 4*_cidx0+2, 4*_cidx0+3;
  Vec4i i1;
  i1 << 4*_cidx1, 4*_cidx1+1, 4*_cidx1+2, 4*_cidx1+3;

  // Hessian is |A^tA   -A^T|
  //            | -A     I_4|
  // with A^TA = I_4
  for(unsigned int i=0; i<4; ++i)
    for(unsigned int j=0; j<4; ++j)
    {
      _coeffs.push_back(Trip(i1[i], i0[j], -A(i,j)));
      _coeffs.push_back(Trip(i0[i], i1[j], -A(j,i)));
    }

  // add diagonal in second quaternion
  for(unsigned int i=0; i<4; ++i)
  {
    _coeffs.push_back(Trip(i0[i], i0[i], 1.0));
    _coeffs.push_back(Trip(i1[i], i1[i], 1.0));
  }
}


//-----------------------------------------------------------------------------

template <class TetMeshT>
void
FrameFieldGeneratorT<TetMeshT>::
add_alignment_coeffs( std::vector<Trip>&     _coeffs,
                           const int              _cidx,
                                 Vec3             _axis,
                           const AxisAlignment    _align,
                           const double           _penalty) const
{
  _axis.normalize();
  Vec3i i0, i1;
  Vec3  c0, c1;

  int bi = 4*_cidx;

  switch(_align)
  {
  case Align_X:
    i0=Vec3i(bi+0,bi+1,bi+2);
    c0=Vec3(_axis[2], -_axis[1], 1.0+_axis[0]);
    i1=Vec3i(bi+0,bi+1,bi+3);
    c1=Vec3(-_axis[1], -_axis[2], 1.0+_axis[0]);
    break;

  case Align_MX:
    i0=Vec3i(bi+0,bi+1,bi+2);
    c0=Vec3(-_axis[2], _axis[1], 1.0-_axis[0]);
    i1=Vec3i(bi+0,bi+1,bi+3);
    c1=Vec3(_axis[1], _axis[2], 1.0-_axis[0]);
    break;

  case Align_Y:
    i0=Vec3i(bi+0,bi+2,bi+3);
    c0=Vec3(_axis[0], -_axis[2], 1.0+_axis[1]);
    i1=Vec3i(bi+0,bi+1,bi+2);
    c1=Vec3(-_axis[2], 1.0+_axis[1], -_axis[0]);
    break;

  case Align_MY:
    i0=Vec3i(bi+0,bi+2,bi+3);
    c0=Vec3(-_axis[0], _axis[2], 1.0-_axis[1]);
    i1=Vec3i(bi+0,bi+1,bi+2);
    c1=Vec3(_axis[2], 1.0-_axis[1], _axis[0]);
    break;

  case Align_Z:
    i0=Vec3i(bi+0,bi+1,bi+3);
    c0=Vec3(_axis[1], 1.0+_axis[2], -_axis[0]);
    i1=Vec3i(bi+0,bi+2,bi+3);
    c1=Vec3(-_axis[0], 1.0+_axis[2], -_axis[1]);
    break;

  case Align_MZ:
    i0=Vec3i(bi+0,bi+1,bi+3);
    c0=Vec3(-_axis[1], 1.0-_axis[2], +_axis[0]);
    i1=Vec3i(bi+0,bi+2,bi+3);
    c1=Vec3(_axis[0], 1.0-_axis[2], +_axis[1]);
    break;

  default:
    return;
  }

  for(unsigned int i=0; i<3; ++i)
    for(unsigned int j=0; j<3; ++j)
    {
      _coeffs.push_back(Trip(i0[i], i0[j], _penalty*c0[i]*c0[j]));
      _coeffs.push_back(Trip(i1[i], i1[j], _penalty*c1[i]*c1[j]));
    }
}


//-----------------------------------------------------------------------------

template <class TetMeshT>
void
FrameFieldGeneratorT<TetMeshT>::
solve_eigen_problem_inverse_power( const std::vector<Trip>& _coeffs, std::vector<double>& _quaternions) const
{
  const double epsilon   = 1e-6;
  const int    max_iters = 200;

  // get vectors in Eigen datastructure
  std::vector<double> q_temp_(_quaternions.size());
  Eigen::Map<Eigen::VectorXd> x(_quaternions.data(), _quaternions.size());
  Eigen::Map<Eigen::VectorXd> y(q_temp_.data()     , q_temp_.size());

  // generate current matrix
  SparseMatrix A(_quaternions.size(),_quaternions.size());
  A.setFromTriplets(_coeffs.begin(),_coeffs.end());
  A.prune(1e-9);

#ifdef ENABLE_SUITESPARSE
  Eigen::CholmodSupernodalLLT<SparseMatrix> chol;
    std::cerr<<"\nsolving with cholmode";
#else
  Eigen::SimplicialLDLT<SparseMatrix> chol;
    std::cerr<<"\nsolving with conjugate";

#endif
  chol.compute(A);

  double residual = DBL_MAX;
  int    iters    = 0;

  do
  {
    // perform next iteration
    y = chol.solve(x);
    y.swap(x);

    // compute relative residual
    x.normalize();
    residual = (x-y).norm();

    std::cerr << "    iter " << iters << " residual " << residual << std::endl;
    ++iters;
  }
  while(iters < max_iters && residual > epsilon);

  // store result
  y=x;
}

//-----------------------------------------------------------------------------

template <class TetMeshT>
void
FrameFieldGeneratorT<TetMeshT>::
solve_unit_length_regularization( const std::vector<Trip>& _coeffs, std::vector<double>& _quaternions) const
{
  std::cerr << "optimize with unit_length regularization..." << std::endl;
  const double w         = 0.5;
  const double epsilon   = 1e-6;
  const int    max_iters = 100;

  // collect constraints for projection
  std::vector< std::pair<CH,Eigen::Quaterniond> > full_alignment;
  std::vector< std::pair<CH, std::pair<AxisAlignment, Vec3> > > partial_alignment;
  collect_alignment_conditions(full_alignment, partial_alignment);

  // get vectors in Eigen datastructure
  std::vector<double> rhs(_quaternions.size());
  Eigen::Map<Eigen::VectorXd> x(_quaternions.data(), _quaternions.size());
  Eigen::Map<Eigen::VectorXd> b(rhs.data()     , rhs.size());
  Eigen::VectorXd x_prev; x_prev = x;

  // update system with regularization terms
  std::vector<Trip> coeffs;
  coeffs = _coeffs;
  for(unsigned int i=0; i<_quaternions.size(); ++i)
    coeffs.push_back(Trip(i,i,w));

  // generate current matrix
  SparseMatrix A(_quaternions.size(),_quaternions.size());
  A.setFromTriplets(coeffs.begin(),coeffs.end());
  A.prune(1e-7);

#ifdef ENABLE_SUITESPARSE
  Eigen::CholmodSupernodalLLT<SparseMatrix> chol;
#else
  Eigen::SimplicialLDLT<SparseMatrix> chol;
#endif
  chol.compute(A);

  double residual = DBL_MAX;
  int    iters    = 0;

  int nq = _quaternions.size()/4;

  do
  {
    // normalize and project to constraints
    normalize_quaternions(_quaternions);
    project_to_alignment_constraints(_quaternions, full_alignment, partial_alignment);

    // setup rhs
    for( int i=0; i<nq; ++i)
    {
      const unsigned int i4 = 4*i;
      Vec4 q(x[i4], x[i4+1], x[i4+2], x[i4+3]);

//        project_to_constrtaints(i, q);

//      q.normalize();
//      if(!std::isfinite(q.squaredNorm()))
//        q = Vec4(0,0,0,0);

      q*=w;

      b[i4  ] = q[0];
      b[i4+1] = q[1];
      b[i4+2] = q[2];
      b[i4+3] = q[3];
    }

    // perform next iteration
    x = chol.solve(b);

    // compute relative residual
    residual = (x-x_prev).norm()/x.norm();
    // store solution for next residual
    x_prev = x;

    std::cerr << "    iter " << iters << " residual " << residual << std::endl;
    ++iters;
  }
  while(iters < max_iters && residual > epsilon);

  // normalize and project to constraints
  normalize_quaternions(_quaternions);
  project_to_alignment_constraints(_quaternions, full_alignment, partial_alignment);
}


//-----------------------------------------------------------------------------

template <class TetMeshT>
double
FrameFieldGeneratorT<TetMeshT>::
measure_energy( const std::vector<Trip>& _coeffs, std::vector<double>& _quaternions) const
{
  Eigen::Map<Eigen::VectorXd> x(_quaternions.data(), _quaternions.size());

  // generate current matrix
  SparseMatrix A(_quaternions.size(),_quaternions.size());
  A.setFromTriplets(_coeffs.begin(),_coeffs.end());

  return x.transpose()*A*x;
}

//-----------------------------------------------------------------------------

template <class TetMeshT>
double
FrameFieldGeneratorT<TetMeshT>::
measure_energy_directly( std::vector<double>& _quaternions) const
{
  double energy(0.0);

  // 3. setup smoothness energy
  for(FIt f_it = mesh_.f_iter(); f_it.valid(); ++f_it)
    if(!mesh_.is_boundary(*f_it))
    {
      // get halffaces and cells
      HFH hfh0 = mesh_.halfface_handle(*f_it, 0);
      HFH hfh1 = mesh_.halfface_handle(*f_it, 1);
      CH  ch0  = mesh_.incident_cell(hfh0);
      CH  ch1  = mesh_.incident_cell(hfh1);

      // get transition from c0 -> c1
      Eigen::Quaterniond trans = qtrans_hfprop_[hfh0];

      Eigen::Quaterniond q0 = get_quaternion(ch0.idx(),_quaternions);
      Eigen::Quaterniond q1 = get_quaternion(ch1.idx(),_quaternions);

      Eigen::Quaterniond q0t = q0*trans;

      double d = (q0t.coeffs()-q1.coeffs()).norm();

//      std::cerr << ch0.idx() << " - " << ch1.idx() << " energy: " << d << std::endl;

      energy += d;
    }

  return energy;
}


//-----------------------------------------------------------------------------

template <class TetMeshT>
void
FrameFieldGeneratorT<TetMeshT>::
normalize_quaternions( std::vector<double>& _quaternions) const
{
  for(unsigned int i=0; i<mesh_.n_cells(); ++i)
  {
    int bidx = 4*i;
    Quaternion qc(_quaternions[bidx], _quaternions[bidx+1], _quaternions[bidx+2], _quaternions[bidx+3]);
    qc.normalize();
    if(!std::isfinite(qc.squaredNorm()))
    {
      std::cerr << "Warning: found degenerate quaternion" << std::endl;
      qc = Quaternion(1,0,0,0);
    }

    _quaternions[bidx  ] = qc.w();
    _quaternions[bidx+1] = qc.x();
    _quaternions[bidx+2] = qc.y();
    _quaternions[bidx+3] = qc.z();
  }
}

//-----------------------------------------------------------------------------


template <class TetMeshT>
void
FrameFieldGeneratorT<TetMeshT>::
project_to_alignment_constraints( std::vector<double>& _quaternions) const
{
  std::vector< std::pair<CH,Eigen::Quaterniond> > full_alignment;
  std::vector< std::pair<CH, std::pair<AxisAlignment, Vec3> > > partial_alignment;
  collect_alignment_conditions(full_alignment, partial_alignment);

  project_to_alignment_constraints(_quaternions, full_alignment, partial_alignment);

}

//-----------------------------------------------------------------------------


template <class TetMeshT>
void
FrameFieldGeneratorT<TetMeshT>::
collect_alignment_conditions( std::vector< std::pair<CH,Eigen::Quaterniond> >& _full_alignment, std::vector< std::pair<CH, std::pair<AxisAlignment, Vec3> > >& _partial_alignment) const
{
  _full_alignment.clear();
  _partial_alignment.clear();

  for(CIt c_it = mesh_.c_iter(); c_it.valid(); ++c_it)
  {
    // get singularity alignment
    HFH hfh_boundary;
    AxisAlignment axis_boundary;
    Vec3 vector_boundary(0,0,0);
    int nba = get_boundary_alignment_hfacet(*c_it, hfh_boundary, axis_boundary);

    if(nba == 1)
      vector_boundary = normal_vector(hfh_boundary);

    // get singularity alignment
    HEH he_singularity;
    AxisAlignment axis_singularity;
    Vec3 vector_singularity(0,0,0);
    int nsa = get_singularity_alignment_hedge(*c_it, he_singularity, axis_singularity);

    if(nsa == 1)
    {
      Point p0 = mesh_.vertex( mesh_.halfedge(he_singularity).from_vertex());
      Point p1 = mesh_.vertex( mesh_.halfedge(he_singularity).to_vertex());
      Point e  = p1-p0;
      e.normalize();
      if(std::isfinite(e.sqrnorm()))
        vector_singularity = Vec3(e[0], e[1], e[2]);
      else
        nsa = 0;
    }

    // partial for boundary
    if(nba == 1 && nsa == 0)
      _partial_alignment.push_back( std::pair<CH, std::pair<AxisAlignment, Vec3> >(*c_it, std::pair<AxisAlignment, Vec3>(axis_boundary, vector_boundary)));

    // partial for singularity
    if(nba == 0 && nsa == 1)
      _partial_alignment.push_back( std::pair<CH, std::pair<AxisAlignment, Vec3> >(*c_it, std::pair<AxisAlignment, Vec3>(axis_singularity, vector_singularity)));

    // full
    if( nba == 1 && nsa == 1)
    {
      Eigen::Quaterniond q = construct_quaternion(axis_boundary, axis_singularity, vector_boundary, vector_singularity);
      _full_alignment.push_back(std::pair<CH,Eigen::Quaterniond>(*c_it,q));
    }
  }
}

//-----------------------------------------------------------------------------


template <class TetMeshT>
Eigen::Quaterniond
FrameFieldGeneratorT<TetMeshT>::
construct_quaternion(const AxisAlignment _axis0, const AxisAlignment _axis1, const Vec3& _vector0, const Vec3& _vector1) const
{
  // only case (X,Y) or (Y,X) implemented so far!!!
  if(_axis0 == Align_X && _axis1 == Align_Y)
  {
    // construct rotation matrix
    Mat3x3 R;
    R.col(0) = _vector0;
    R.col(0).normalize();
    R.col(1) = _vector1;
    R.col(1).normalize();
    R.col(2) = R.col(0).cross(R.col(1));

    if(!std::isfinite(R.norm()))
    {
      std::cerr << "Warning: get_quaternion frame from two axes failed!!!" << std::endl;
      return Eigen::Quaterniond(0,0,0,0);
    }

    Eigen::Quaterniond q = Eigen::Quaterniond(R);

//    std::cerr << "input axes/vectors: " << int(_axis0) << ", " << int(_axis1) << "  /  " << _vector0 << " --- " << _vectors1 << std::endl;
//    std::cerr << R << std::endl;
//    std::cerr << " versus " << std::endl;
//    std::cerr << q.toRotationMatrix() << std::endl;

    // return quaternion of rotation matrix
    return q;
  }
  else if(_axis1 == Align_X && _axis0 == Align_Y)
  {
    // swap input
    return(construct_quaternion(_axis1, _axis0, _vector1, _vector0));
  }
  else
  {
    std::cerr << "Warning: request case of get_quaternion(const AxisAlignment _axis0, const AxisAlignment _axis1, const Vec3& _vector0, const Vec3& _vector1) not implemented yet!!!" << std::endl;
    return Eigen::Quaterniond(1,0,0,0);
  }
}


//-----------------------------------------------------------------------------


template <class TetMeshT>
void
FrameFieldGeneratorT<TetMeshT>::
project_to_alignment_constraints( std::vector<double>& _quaternions, const std::vector< std::pair<CH,Eigen::Quaterniond> >& _full_alignment, const std::vector< std::pair<CH, std::pair<AxisAlignment, Vec3> > >& _partial_alignment) const
{
  // set full alignment conditions
  for(unsigned int i=0; i<_full_alignment.size(); ++i)
  {
    int idx = _full_alignment[i].first.idx();
    Eigen::Quaterniond q  = get_quaternion(idx, _quaternions);
    Eigen::Quaterniond qc = _full_alignment[i].second;
    // choose sign according to current quaternion
    if(q.dot(qc) < 0.0)
      qc = negate(qc);

    set_quaternion(idx, qc, _quaternions);
  }

  // set partial alignment conditions
  for(unsigned int i=0; i<_partial_alignment.size(); ++i)
  {
    int idx = _partial_alignment[i].first.idx();
    Eigen::Quaterniond q = get_quaternion(idx, _quaternions);

    // project to axis constraint
    AxisAlignment axis_label  = _partial_alignment[i].second.first;
    Vec3          axis_vector = _partial_alignment[i].second.second;
    Eigen::Quaterniond qp = project_quaternion(q, axis_label, axis_vector);

//    std::cerr << "original    " << std::endl << q << std::endl;
//    std::cerr << "  projected " << std::endl << qp << std::endl;

    set_quaternion(idx, qp, _quaternions);
  }
}


//-----------------------------------------------------------------------------


template <class TetMeshT>
Eigen::Quaterniond
FrameFieldGeneratorT<TetMeshT>::
project_quaternion(const Eigen::Quaterniond& _q, const AxisAlignment _axis_label, const Vec3& _axis_vector) const
{
  // rotate _q to coordinate system, where _axis_vector is the x-axis
  Eigen::Quaterniond R;
  Vec3 a = AxisAlignmentHelpers::vector(_axis_label);
  R.setFromTwoVectors(a,_axis_vector);
  Eigen::Quaterniond Rip = R.inverse()*_q;

  // project rotated quaternion and normalize
  Eigen::Quaterniond qp;
  switch(_axis_label)
  {
    case Align_X :
    case Align_MX: qp = Eigen::Quaterniond(Rip.w(), Rip.x(), 0, 0); break;
    case Align_Y :
    case Align_MY: qp = Eigen::Quaterniond(Rip.w(), 0, Rip.y(), 0); break;
    case Align_Z :
    case Align_MZ: qp = Eigen::Quaterniond(Rip.w(), 0, 0, Rip.z()); break;
    case Align_NONE: return _q; // no projection
  }

  qp.normalize();

  if(!std::isfinite(qp.squaredNorm()))
  {
    std::cerr << "Warning: projection fo quaternion onto partial constraint failed!!!" << std::endl;
    return _q;
  }

  // express solution in original coordinate system
  qp = R*qp;

  // choose sign according to dot product between quaternions
  if(qp.dot(_q) < 0.0)
    qp = negate(qp);

  return qp;
}


//-----------------------------------------------------------------------------


template <class TetMeshT>
int
FrameFieldGeneratorT<TetMeshT>::
number_of_adjacent_boundary_faces(const CH _ch) const
{
  int abf(0);

  for(unsigned int i=0; i<mesh_.cell(_ch).halffaces().size(); ++i)
  {
    const HFH ho = mesh_.opposite_halfface_handle(mesh_.cell(_ch).halffaces()[i]);
    if(mesh_.is_boundary(ho))
      ++abf;
  }

  return abf;
}


//-----------------------------------------------------------------------------


template <class TetMeshT>
int
FrameFieldGeneratorT<TetMeshT>::
get_boundary_alignment_hfacet(const CH _ch, HFH& _hfh_boundary, AxisAlignment& _axis) const
{
  // return number of faces with boundary alignment
  int nba(0);
  _axis = Align_NONE;
  _hfh_boundary.reset();

  for(unsigned int i=0; i<mesh_.cell(_ch).halffaces().size(); ++i)
  {
    const HFH ho = mesh_.opposite_halfface_handle(mesh_.cell(_ch).halffaces()[i]);
    if(mesh_.is_boundary(ho))
    {
      _hfh_boundary = ho;
      _axis = AxisAlignment(normal_hfprop_[ho]);
      ++nba;
    }
  }

  return nba;
}


//-----------------------------------------------------------------------------


template <class TetMeshT>
int
FrameFieldGeneratorT<TetMeshT>::
get_singularity_alignment_hedge(const CH _ch, HEH& _he_singularity, AxisAlignment& _axis) const
{
  // return number of singular edges
  int nsa(0);
  _axis = Align_NONE;
  _he_singularity.reset();

  for(unsigned int i=0; i<mesh_.cell(_ch).halffaces().size(); ++i)
  {
    const HFH ho = mesh_.cell(_ch).halffaces()[i];
    for(unsigned int j=0; j<mesh_.halfface(ho).halfedges().size(); ++j)
    {
      HEH heh = mesh_.halfface(ho).halfedges()[j];
      std::map<HEH, int>::const_iterator m_it = axes_cprop_[_ch].find(heh);
      if(m_it != axes_cprop_[_ch].end() && AxisAlignment(m_it->second) != Align_NONE)
      {
        _he_singularity = heh;
        _axis = AxisAlignment(m_it->second);
        ++nsa;
      }
    }
  }
  return nsa;
}


//-----------------------------------------------------------------------------


template <class TetMeshT>
int
FrameFieldGeneratorT<TetMeshT>::
ccw_transition(const CH _ch, const HEH _he) const
{
  int trans_idx(0);

  HFH hfh_start = find_adjacent_halfface(_ch, _he);
  if(hfh_start.is_valid())
  {
    HFH hfh_cur   = hfh_start;
    do
    {
      // move to next halfface in same cell
      hfh_cur = mesh_.adjacent_halfface_in_cell (hfh_cur, _he);
      HFH hfh_opp = mesh_.opposite_halfface_handle(hfh_cur);
      // stop at boundary!!!
      if(mesh_.is_boundary(hfh_opp)) break;
      else
      {
        // get new transition and update total transition
//        trans_idx = tq_.mult_transitions_idx(trans_hfprop_[hfh_cur], trans_idx);
        trans_idx = tq_.mult_transitions_idx(trans_idx, trans_hfprop_[hfh_cur]);
        // move to opposite halfface
        hfh_cur = hfh_opp;
      }
    }
    while(hfh_cur != hfh_start);

    // return accumulated transition
    return trans_idx;
  }

  // return identity
  return 0;
}


//-----------------------------------------------------------------------------


template <class TetMeshT>
Eigen::Quaterniond
FrameFieldGeneratorT<TetMeshT>::
ccw_transition_quaternion(const CH _ch, const HEH _he) const
{
  Eigen::Quaterniond trans(1,0,0,0);

  HFH hfh_start = find_adjacent_halfface(_ch, _he);
  if(hfh_start.is_valid())
  {
    HFH hfh_cur   = hfh_start;
    do
    {
      // move to next halfface in same cell
      hfh_cur = mesh_.adjacent_halfface_in_cell (hfh_cur, _he);
      HFH hfh_opp = mesh_.opposite_halfface_handle(hfh_cur);
      // stop at boundary!!!
      if(mesh_.is_boundary(hfh_opp)) break;
      else
      {
        // get new transition and update total transition
        trans = trans * qtrans_hfprop_[hfh_cur];
//        trans_idx = tq_.mult_transitions_idx(trans_idx, trans_hfprop_[hfh_cur]);
        // move to opposite halfface
        hfh_cur = hfh_opp;
      }
    }
    while(hfh_cur != hfh_start);

    // return accumulated transition
    return trans;
  }

  // return identity
  return trans;
}


//-----------------------------------------------------------------------------


template <class TetMeshT>
typename FrameFieldGeneratorT<TetMeshT>::HFH
FrameFieldGeneratorT<TetMeshT>::
find_adjacent_halfface(const CH _ch, const HEH _heh) const
{
  for(unsigned int i=0; i<mesh_.cell(_ch).halffaces().size(); ++i)
  {
    const HFH hfh = mesh_.cell(_ch).halffaces()[i];
    for(unsigned int j=0; j<mesh_.halfface(hfh).halfedges().size(); ++j)
    {
      HEH heh = mesh_.halfface(hfh).halfedges()[j];
      // found ?
      if(heh == _heh)
        return hfh;
    }
  }
  // if not found return invalid handle
  return HFH(-1);
}


//-----------------------------------------------------------------------------


template <class TetMeshT>
bool
FrameFieldGeneratorT<TetMeshT>::
is_consistent(const int _transition, const AxisAlignment _axis, const int _valance) const
{
  if(!tq_.is_2d_rotation(_transition))
  {
    return false;
  }
  else
  {
    int rot_axis  = tq_.rotation_axis_2d (_transition);
    int rot_angle = tq_.rotation_angle_2d(_transition);
    normalize_rotation(rot_axis, rot_angle);

    int e_axis  = int(_axis);
    int e_angle = _valance*90;
    normalize_rotation(e_axis, e_angle);

    // special case of identity
    if(rot_angle == 0 && _valance == 0)
      return true;
    else if(rot_axis == e_axis && rot_angle == e_angle)
      return true;
      else
      {
        std::cerr << "ERROR: inconsistent transition:" << std::endl;
        std::cerr << "transition  axis " << rot_axis   << " angle " << rot_angle << std::endl;
        std::cerr << "edge        axis " << int(_axis) << " index " << _valance << std::endl;
        return false;
      }
  }
}


//-----------------------------------------------------------------------------


template <class TetMeshT>
void
FrameFieldGeneratorT<TetMeshT>::
normalize_rotation(int& _axis, int& _angle) const
{
  // assert that _angle is positive
  if(_angle < 0)
  {
    // negate angle and axis
    _angle *= -1;
    if(_axis % 1)
      --_axis;
    else
      ++_axis;
  }
}


//-----------------------------------------------------------------------------


template <class TetMeshT>
void
FrameFieldGeneratorT<TetMeshT>::
adjacent_hehs(const CH _ch, std::vector<HEH>& _adj_heh) const
{
  // clear set
  _adj_heh.clear();

  for(unsigned int i=0; i<mesh_.cell(_ch).halffaces().size(); ++i)
  {
    const HFH hfh = mesh_.cell(_ch).halffaces()[i];
    for(unsigned int j=0; j<mesh_.halfface(hfh).halfedges().size(); ++j)
      _adj_heh.push_back(mesh_.halfface(hfh).halfedges()[j]);
  }
}


//-----------------------------------------------------------------------------


template <class TetMeshT>
typename FrameFieldGeneratorT<TetMeshT>::CH
FrameFieldGeneratorT<TetMeshT>::
get_boundary_start_cell(const HEH _heh) const
{
  CH ch(-1);

  for(HEHFIt hehf_it = mesh_.hehf_iter(_heh); hehf_it.valid(); ++hehf_it)
    if(mesh_.is_boundary(mesh_.opposite_halfface_handle(*hehf_it)))
    {
      if(ch.is_valid()) std::cerr << "Warning: get_boundary_start_cell found several candidates -> input mesh not 3-manifold!!!" << std::endl;

      ch = mesh_.incident_cell(*hehf_it);
    }
  return ch;
}


//-----------------------------------------------------------------------------


template <class TetMeshT>
typename FrameFieldGeneratorT<TetMeshT>::CH
FrameFieldGeneratorT<TetMeshT>::
get_boundary_end_cell(const HEH _heh) const
{
  CH ch(-1);

  for(HEHFIt hehf_it = mesh_.hehf_iter(_heh); hehf_it.valid(); ++hehf_it)
    if(mesh_.is_boundary(*hehf_it))
    {
      if(ch.is_valid()) std::cerr << "Warning: get_boundary_end_cell found several candidates -> input mesh not 3-manifold!!!" << std::endl;

      ch = mesh_.incident_cell(mesh_.opposite_halfface_handle(*hehf_it));
    }
  return ch;
}

//-----------------------------------------------------------------------------


template <class TetMeshT>
typename FrameFieldGeneratorT<TetMeshT>::Vec3
FrameFieldGeneratorT<TetMeshT>::
normal_vector(const HFH _hfh) const
{
  // get halfface
  Face f = mesh_.halfface(_hfh);
  // get points
  Point p0 = mesh_.vertex( mesh_.halfedge(f.halfedges()[0]).to_vertex());
  Point p1 = mesh_.vertex( mesh_.halfedge(f.halfedges()[1]).to_vertex());
  Point p2 = mesh_.vertex( mesh_.halfedge(f.halfedges()[2]).to_vertex());
  // get outward normal vector
  Point n = (p1-p0)%(p2-p0);
  n.normalize();
  if(std::isfinite(n.sqrnorm()))
    return Vec3(n[0], n[1], n[2]);
  else
    return Vec3(0,0,0);
}


//-----------------------------------------------------------------------------


template <class TetMeshT>
void
FrameFieldGeneratorT<TetMeshT>::
init_quaternions_by_propagation(std::vector<double>& _quaternions) const
{
  // spanning tree of cells
  std::vector<bool> cell_visited(mesh_.n_cells(),false);

  // qeueue
  std::queue<HFH> cq;

  // select root
  CH ch(0);
  cell_visited[ch.idx()] = true;
  // push hf to neighbors
  for(unsigned int i=0; i<mesh_.cell(ch).halffaces().size(); ++i)
    cq.push(mesh_.cell(ch).halffaces()[i]);


  while(!cq.empty())
  {
    //get next halfface
    HFH hfh = cq.front();
    cq.pop();
    HFH hfh_opp = mesh_.opposite_halfface_handle(hfh);
    if(!mesh_.is_boundary(hfh_opp))
    {
      CH  ch     = mesh_.incident_cell(hfh);
      CH  ch_opp = mesh_.incident_cell(hfh_opp);

      if(!cell_visited[ch_opp.idx()])
      {
        // mark as visited
        cell_visited[ch_opp.idx()] = true;
        // propagate
        Eigen::Quaterniond q0 = get_quaternion(ch.idx(), _quaternions);

        // get transition
        Eigen::Quaterniond trans = qtrans_hfprop_[hfh];
        // propagate quaternion
        Eigen::Quaterniond q1 = q0*trans;

//        // get transition
//        int trans = trans_hfprop_[hfh];
//        // propagate quaternion
//        Eigen::Quaterniond q1 = q0*acg2eigen(tq_.transition(trans));

        // store propagated quaternion
        set_quaternion(ch_opp.idx(), q1, _quaternions);
        // push halffaces to neighbors
        for(unsigned int i=0; i<mesh_.cell(ch_opp).halffaces().size(); ++i)
          cq.push(mesh_.cell(ch_opp).halffaces()[i]);

//        // check propagation
//        Eigen::Quaterniond qa = get_quaternion(ch.idx(),_quaternions)*acg2eigen(tq_.transition(trans));
//        Eigen::Quaterniond qb = get_quaternion(ch_opp.idx(),_quaternions);
//        std::cerr << "residual norm after propagation " << qa.angularDistance(qb) << std::endl;
      }
    }
  }

  // check whether all cells have been visited
  for(unsigned int i=0; i<cell_visited.size(); ++i)
    if(cell_visited[i] == false)
      std::cerr << "Warning: cell has not been visited!!!" << std::endl;
}

//-----------------------------------------------------------------------------


template <class TetMeshT>
Eigen::Quaterniond
FrameFieldGeneratorT<TetMeshT>::
get_quaternion(const unsigned int _i, const std::vector<double>& _quaternions) const
{
  unsigned int i4 = 4*_i;
  return Eigen::Quaterniond(_quaternions[i4], _quaternions[i4+1], _quaternions[i4+2], _quaternions[i4+3]);
}


//-----------------------------------------------------------------------------


template <class TetMeshT>
void
FrameFieldGeneratorT<TetMeshT>::
set_quaternion(const unsigned int _i, const Eigen::Quaterniond& _q, std::vector<double>& _quaternions) const
{
  unsigned int i4 = 4*_i;

  _quaternions[i4  ] = _q.w();
  _quaternions[i4+1] = _q.x();
  _quaternions[i4+2] = _q.y();
  _quaternions[i4+3] = _q.z();
}

//-----------------------------------------------------------------------------


template <class TetMeshT>
void
FrameFieldGeneratorT<TetMeshT>::
save_frames(const std::string _filename) const
{
  // quaternion property
  OpenVolumeMesh::CellPropertyT< Eigen::Quaterniond > quaternion_prop = mesh_.template request_cell_property< Eigen::Quaterniond >("CellQuaternion");

  std::ofstream f_write(_filename);

  for(CIt c_it = mesh_.c_iter(); c_it.valid(); ++c_it)
  {
    // get quaternion
    Eigen::Quaterniond q = quaternion_prop[*c_it];
    q.normalize();
    if(!std::isfinite(q.squaredNorm()))
      q = Eigen::Quaterniond(1,0,0,0);

    // get frame matrix
    Mat3x3 R = q.toRotationMatrix();

    for(unsigned int j=0; j<3; ++j)
      for(unsigned int k=0; k<3; ++k)
        f_write << R(k,j) << " ";

    f_write << std::endl;
  }

  f_write.close();
}

//-----------------------------------------------------------------------------


template <class TetMeshT>
void
FrameFieldGeneratorT<TetMeshT>::
test_loop_conditions() const
{
  std::cerr << "test loop conditions..." << std::endl;
  int total_valid(0);
  int total_invalid(0);
  for(unsigned int i=0; i<6; ++i)
    for(unsigned int j=0; j<6; ++j)
      for(unsigned int k=0; k<6; ++k)
      {
        // get transitions around edges
        int ti = i+4;
        int tj = j+4;
        int tk = k+4;

        AxisAlignment ai = AxisAlignment(tq_.rotation_axis_2d(ti));
        AxisAlignment aj = AxisAlignment(tq_.rotation_axis_2d(tj));
        AxisAlignment ak = AxisAlignment(tq_.rotation_axis_2d(tk));

        Mat3x3 F;
        F.col(0) = -AxisAlignmentHelpers::vector(ai);
        F.col(1) = -AxisAlignmentHelpers::vector(aj);
        F.col(2) = -AxisAlignmentHelpers::vector(ak);

        double d = F.determinant();

        bool valid_rhframe = (std::abs(d-1.0) < 1e-4);
        bool valid_lhframe = (std::abs(d+1.0) < 1e-4);

        int triple_prod = tq_.mult_transitions_idx(tk, tq_.mult_transitions_idx(tj,ti));

        bool valid_triple_rh = ( triple_prod >= 10 && triple_prod <=15);
        bool valid_triple_lh = ( triple_prod >= 4  && triple_prod <= 9);

        int pair0 = tq_.mult_transitions_idx(tj,ti);
        int pair1 = tq_.mult_transitions_idx(tk,tj);
        int pair2 = tq_.mult_transitions_idx(ti,tk);

        bool valid_pair0 = ( pair0 >= 16 && pair0 <=23);
        bool valid_pair1 = ( pair1 >= 16 && pair1 <=23);
        bool valid_pair2 = ( pair2 >= 16 && pair2 <=23);

        bool all_classes_correct_rh = (valid_triple_rh && valid_pair0 && valid_pair1 && valid_pair2);
        bool all_classes_correct_lh = (valid_triple_lh && valid_pair0 && valid_pair1 && valid_pair2);

//        bool all_pairs_valid = (valid_pair0 && valid_pair1 && valid_pair2);
//        bool valid_simple_case = !(all_pairs_valid^(valid_rhframe||valid_lhframe));

        // test XOR
        bool valid_case = !(valid_rhframe^all_classes_correct_rh) || !(valid_lhframe^all_classes_correct_lh);
        if(valid_case)
          ++total_valid;
        else
          ++total_invalid;

        std::cerr << ti << "," << tj << "," << tk
                  << std::setw(5) << " determinant " << d
                  << " right-handed " << int(valid_rhframe)
                  << " left-handed  " << int(valid_lhframe)
                  << " triple type " << triple_prod
                  << " valid triple rh " << int(valid_triple_rh)
                  << " valid triple lh " << int(valid_triple_lh)
                  << " valid pairs " << int(valid_pair0) << "," << int(valid_pair1) << "," << int(valid_pair2)
                  << " valid case " << int(valid_case)
//                  << " valid simple case " << valid_simple_case
                  << std::endl;

      }

  std::cerr << "#valid  : " << total_valid << std::endl;
  std::cerr << "#invalid: " << total_invalid << std::endl;
  std::cerr << "end test loop conditions!" << std::endl;
}

//=============================================================================
} // namespace SCOF
//=============================================================================
