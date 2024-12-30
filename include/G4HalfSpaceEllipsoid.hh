#pragma once

#include "G4ThreeVector.hh"
#include "G4VHalfSpace.hh"
#include "G4HalfSpaceZone.hh"
#include "G4HalfSpaceTransformation.hh"

class G4HalfSpaceEllipsoid : public G4VHalfSpace {
public:
  G4HalfSpaceEllipsoid();
  G4HalfSpaceEllipsoid(const G4ThreeVector &f1, const G4ThreeVector &f2, G4double l);
  G4HalfSpaceEllipsoid(const G4ThreeVector &abc, const G4ThreeVector &centre, const G4ThreeVector &rotation);
  void ComputeSurfaces();
  ~G4HalfSpaceEllipsoid() = default;

  virtual G4double Sdf(const G4ThreeVector&p) const override;
  std::vector<G4ThreeVector> Intersection(const G4ThreeVector& p, const G4ThreeVector& v) const override;

  virtual void Translate(const G4ThreeVector& t) override;
  virtual void Rotate(const G4RotationMatrix& r) override;
  virtual void Transform(const G4AffineTransform& a) override;

  virtual G4SurfaceMeshCGAL* GetSurfaceMesh()  override;


protected:
  G4ThreeVector _f1;
  G4ThreeVector _f2;
  G4double _l;

  G4double _c;
  G4double _r;

  /*
  double _a;
  double _b;
  double _c;
  G4ThreeVector _centre;
  G4ThreeVector _rotation;
   */

  G4HalfSpaceZone _hsZone = G4HalfSpaceZone();
};
