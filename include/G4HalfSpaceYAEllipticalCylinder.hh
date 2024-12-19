#pragma once

#include "G4ThreeVector.hh"
#include "G4VHalfSpace.hh"

class G4HalfSpaceYAEllipticalCylinder: public G4VHalfSpace {
public:
  G4HalfSpaceYAEllipticalCylinder();
  G4HalfSpaceYAEllipticalCylinder(G4double x, G4double z, G4double r1, G4double r2);
  ~G4HalfSpaceYAEllipticalCylinder();

  virtual G4double Sdf(const G4ThreeVector&p) const override;
  virtual std::vector<G4ThreeVector> Intersection(const G4ThreeVector& p, const G4ThreeVector &d) const override;

  virtual void Translate(const G4ThreeVector& t) override;
  virtual void Rotate(const G4RotationMatrix& r) override;
  virtual void Transform(const G4AffineTransform& a) override;

  virtual G4SurfaceMeshCGAL* GetSurfaceMesh()  override;

protected:
  G4double _x0 = 0;
  G4double _z0 = 0;
  G4double _r1 = 1.0;
  G4double _r2 = 1.0;
};