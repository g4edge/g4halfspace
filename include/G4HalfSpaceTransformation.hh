#pragma once

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4AffineTransform.hh"

class G4HalfSpaceTransformation {
public:
  G4HalfSpaceTransformation() = default;

  G4HalfSpaceTransformation(const G4ThreeVector &t,
                            const G4ThreeVector &r) : _t(t), _r(r) {
  }

  G4HalfSpaceTransformation(const G4ThreeVector &v) {

    G4RotationMatrix rmatrix;

    auto vnorm = v / v.mag();
    G4ThreeVector zaxis = G4ThreeVector(0, 0, 1);
    G4ThreeVector raxis = zaxis.cross(vnorm);
    G4double angle = asin(raxis.mag());

    if (angle == 0) {
      raxis = G4ThreeVector(0, 0, 1);
    }
    rmatrix.set(raxis, angle);
    rmatrix.rectify();

    auto a_11 = rmatrix.xx();
    auto a_12 = rmatrix.xy();
    auto a_13 = rmatrix.xz();

    auto a_21 = rmatrix.yx();
    auto a_22 = rmatrix.yy();
    auto a_23 = rmatrix.yz();

    auto a_31 = rmatrix.zx();
    auto a_32 = rmatrix.zy();
    auto a_33 = rmatrix.zz();

    if( fabs(a_31-1) < 1e-9 && a_31 > 1.0)
      a_31 = 1.0;
    else if( fabs(a_31+1) < 1e-9 && a_31 < -1.0)
      a_31 = -1.0;

    double x;
    double y;
    double z;

    if( fabs(a_31) != 1) {
      x = atan2(a_32, a_33);
      y = asin(-a_31);
      z = atan2(a_21, a_11);
    }
    else if(a_31 == -1) {
      z = 0.0;
      x = atan2(a_21, a_11) + z;
      y = M_PI/2.0;
    }
    else if(a_31 == 1) {
      z = 0.0;
      x = atan2(-a_12, -a_13) - z;
      y = -M_PI/2.0;
    }

    _r.set(x,y,z);
  }

  G4HalfSpaceTransformation(const G4AffineTransform &a) {
    auto rmat = a.NetRotation();
    _r = G4ThreeVector(rmat.thetaX(), rmat.thetaY(), rmat.thetaZ());
    _t = a.NetTranslation();
  }

  G4RotationMatrix GetRotationMatrix() {
    G4RotationMatrix rmat;
    rmat.rotateX(_r.x());
    rmat.rotateY(_r.y());
    rmat.rotateZ(_r.z());
    rmat.rectify();
    return rmat;
  }

  G4ThreeVector GetRotation() {
    return _r;
  }

  G4ThreeVector GetTranslation() {
    return _t;
  }

  G4AffineTransform GetAffineTransform() {
    auto rmat = GetRotationMatrix();
    return G4AffineTransform(rmat, _t);
  }

private:
  G4ThreeVector _t = G4ThreeVector();
  G4ThreeVector _r = G4ThreeVector();
};