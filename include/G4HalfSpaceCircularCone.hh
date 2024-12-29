#pragma once

#include "G4ThreeVector.hh"
#include "G4VHalfSpace.hh"
#include "G4HalfSpaceZone.hh"

class G4HalfSpaceZone;

class G4HalfSpaceCircularCone : public G4VHalfSpace {
public:
  G4HalfSpaceCircularCone();
  G4HalfSpaceCircularCone(const G4ThreeVector &v, const G4ThreeVector &h,
                          G4double r1, G4double r2);
  G4HalfSpaceCircularCone(G4double h, G4double r1, G4double r2,
                          const G4ThreeVector &v,
                          const G4ThreeVector &r);
  ~G4HalfSpaceCircularCone();

  virtual G4double Sdf(const G4ThreeVector&p) const override;
  std::vector<G4ThreeVector> Intersection(const G4ThreeVector& p, const G4ThreeVector& v) const override;

  virtual void Translate(const G4ThreeVector& t) override;
  virtual void Rotate(const G4RotationMatrix& r) override;
  virtual void Transform(const G4AffineTransform& a) override;

  virtual G4SurfaceMeshCGAL* GetSurfaceMesh()  override;


protected:
  G4double _h = 25;
  G4double _r1 = 10;
  G4double _r2 = 20;
  G4ThreeVector _v = G4ThreeVector(0,0,0);
  G4ThreeVector _r = G4ThreeVector(0,0,0);

  G4HalfSpaceZone _hsZone = G4HalfSpaceZone();
};