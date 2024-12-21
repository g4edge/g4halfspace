#pragma once

#include "G4ThreeVector.hh"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/Vector.h"

#include "G4VHalfSpace.hh"

class G4HalfSpaceQuadric : public G4VHalfSpace {
public:
  G4HalfSpaceQuadric();

  G4HalfSpaceQuadric(double qxx, double qxy, double qxz,
                                 double qyy, double qyz,
                                             double qzz,
                     double px,  double py,  double pz,
                     double r);

  G4HalfSpaceQuadric(CLHEP::HepMatrix &q,
                     CLHEP::HepVector &p,
                     double r);

  virtual G4double Sdf(const G4ThreeVector&p) const override;
  std::vector<G4ThreeVector> Intersection(const G4ThreeVector& p, const G4ThreeVector& v) const override;

  virtual void Translate(const G4ThreeVector& t) override;
  void Rotate(const G4ThreeVector &rv);
  virtual void Rotate(const G4RotationMatrix& r) override;
  virtual void Transform(const G4AffineTransform& a) override;

  virtual G4SurfaceMeshCGAL* GetSurfaceMesh()  override;

protected:
  CLHEP::HepMatrix _q = CLHEP::HepMatrix(3,3);
  CLHEP::HepVector _p = CLHEP::HepVector(3);
  double _r = 0;

};
