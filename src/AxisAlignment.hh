//=============================================================================
//
//  ENUM AxisAlignment
//
//=============================================================================


#ifndef AXISALIGNMENT_HH
#define AXISALIGNMENT_HH


//== INCLUDES =================================================================

#include <Eigen/Dense>

//== NAMESPACES ===============================================================

namespace SCOF {

//== ENUM DEFINITION =========================================================



/** \enum AxisAlignment  AxisAlignment.hh

    Brief Description.
  
    A more elaborate description follows.
*/

  enum AxisAlignment {Align_NONE = -1, Align_X = 0, Align_MX = 1, Align_Y = 2, Align_MY = 3, Align_Z = 4, Align_MZ = 5};

class AxisAlignmentHelpers
{
public:
    static Eigen::Vector3d vector(const AxisAlignment _a)
    {
      switch(_a)
      {
      case Align_NONE: return Eigen::Vector3d( 0, 0, 0);
      case Align_X   : return Eigen::Vector3d( 1, 0, 0);
      case Align_MX  : return Eigen::Vector3d(-1, 0, 0);
      case Align_Y   : return Eigen::Vector3d( 0, 1, 0);
      case Align_MY  : return Eigen::Vector3d( 0,-1, 0);
      case Align_Z   : return Eigen::Vector3d( 0, 0, 1);
      case Align_MZ  : return Eigen::Vector3d( 0, 0,-1);
      }
    }

    static AxisAlignment get_dominant_axis(const Eigen::Vector3d& _v)
    {
      int idx(0);
      if(std::abs(_v[1]) > std::abs(_v[idx]))
        idx = 1;
      if(std::abs(_v[2]) > std::abs(_v[idx]))
        idx = 2;

      if(_v[idx] >= 0.0)
        idx = 2*idx; // positive
      else
        idx = 2*idx + 1; // negative

      return AxisAlignment(idx);
    }
};

  
  //=============================================================================
} // namespace SCOF
//=============================================================================
//=============================================================================
#endif // AXISALIGNMENT_HH defined
//=============================================================================

