#pragma once

#include "G4RotationMatrix.hh"
#include "CLHEP/Matrix/Matrix.h"

class G4HalfSpaceRotation {

public:
  G4HalfSpaceRotation(const G4ThreeVector &tb_xyz) : _tb_xyz(tb_xyz) {}

  G4HalfSpaceRotation(CLHEP::HepMatrix &matrix) {
    auto a_11 = matrix(1,1);
    auto a_12 = matrix(1,2);
    auto a_13 = matrix(1,3);

    auto a_21 = matrix(2,1);
    auto a_22 = matrix(2,2);
    auto a_23 = matrix(2,3);

    auto a_31 = matrix(3,1);
    auto a_32 = matrix(3,2);
    auto a_33 = matrix(3,3);

    if( fabs(a_31-1) < 1e-9 && a_31 > 1.0)
      a_31 = 1.0;
    else if( fabs(a_31+1) < 1e-9 && a_31 < -1.0)
      a_31 = -1.0;

    double x = 0;
    double y = 0;
    double z = 0;

    if( fabs(a_31) != 1) {
      x = atan2(a_32, a_33);
      y = asin(-a_31);
      z = atan2(a_21, a_11);
    }
    else if(a_31 == -1) {
      z = 0.0;
      x = atan2(a_12, a_13) + z;
      y = M_PI/2.0;
    }
    else if(a_31 == 1) {
      z = 0.0;
      x = atan2(-a_12, -a_13) - z;
      y = -M_PI/2.0;
    }

    _tb_xyz.set(x,y,z);
  }

  std::string ToString() {
    std::stringstream ss;
    ss << _tb_xyz;

    return ss.str();
  }

protected:
  G4ThreeVector _tb_xyz= G4ThreeVector();
};