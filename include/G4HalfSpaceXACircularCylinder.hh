#pragma once

#include "G4ThreeVector.hh"
#include "G4VHalfSpace.hh"

class G4HalfSpaceXACircularCylinder: public G4VHalfSpace {
public:
  G4HalfSpaceXACircularCylinder();
  G4HalfSpaceXACircularCylinder(G4double y, G4double z, G4double r);
  ~G4HalfSpaceXACircularCylinder();

  virtual G4double Sdf(const G4ThreeVector&p) const override;
  virtual std::vector<G4ThreeVector> Intersection(const G4ThreeVector& p, const G4ThreeVector &d) const override;

  virtual void Translate(const G4ThreeVector& t) override;
  virtual void Rotate(const G4RotationMatrix& r) override;
  virtual void Transform(const G4AffineTransform& a) override;

  virtual G4SurfaceMeshCGAL* GetSurfaceMesh() override;

protected:
  G4double _y0 = 0;
  G4double _z0 = 0;
  G4double _r = 1.0;
};