#pragma once

#include "G4ThreeVector.hh"
#include "G4VHalfSpace.hh"

class G4HalfSpaceZone;

class G4HalfSpaceYZPlane : public G4VHalfSpace {
public:
  G4HalfSpaceYZPlane();
  G4HalfSpaceYZPlane(G4double x);
  ~G4HalfSpaceYZPlane();

  virtual G4double Sdf(const G4ThreeVector&p) const override;
  std::vector<G4ThreeVector> Intersection(const G4ThreeVector& p, const G4ThreeVector& v) const override;

  virtual void Translate(const G4ThreeVector& t) override;
  virtual void Rotate(const G4RotationMatrix& r) override;
  virtual void Transform(const G4AffineTransform& a) override;

  virtual G4SurfaceMeshCGAL* GetSurfaceMesh()  override;

protected:
  G4double _x =  0.0;
};