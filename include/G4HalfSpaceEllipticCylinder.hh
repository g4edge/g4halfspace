#pragma once

#include "G4ThreeVector.hh"
#include "G4VHalfSpace.hh"
#include "G4HalfSpaceZone.hh"

class G4HalfSpaceEllipticCylinder: public G4VHalfSpace {
public:
  G4HalfSpaceEllipticCylinder();
  G4HalfSpaceEllipticCylinder(const G4ThreeVector& v, const G4ThreeVector& h,
                              const G4ThreeVector& r1, const G4ThreeVector &r2);
  ~G4HalfSpaceEllipticCylinder();

  virtual G4double Sdf(const G4ThreeVector&p) const override;
  virtual std::vector<G4ThreeVector> Intersection(const G4ThreeVector& p, const G4ThreeVector &d) const override;

  virtual void Translate(const G4ThreeVector& t) override;
  virtual void Rotate(const G4RotationMatrix& r) override;
  virtual void Transform(const G4AffineTransform& a) override;

  virtual G4SurfaceMeshCGAL* GetSurfaceMesh()  override;

protected:
  G4ThreeVector _v;
  G4ThreeVector _h;
  G4ThreeVector _r1;
  G4ThreeVector _r2;

  G4HalfSpaceZone _hsZone = G4HalfSpaceZone();
};