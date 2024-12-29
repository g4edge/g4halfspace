#pragma once

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4AffineTransform.hh"

class G4HalfSpaceTransformation {
public:
  G4HalfSpaceTransformation() = default;
  G4HalfSpaceTransformation(const G4ThreeVector &t,
                            const G4ThreeVector &r) :
                            _t(t), _r(r) {}

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