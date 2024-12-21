#pragma once

#include "G4ThreeVector.hh"
#include "G4VHalfSpace.hh"
#include "G4HalfSpaceZone.hh"

class G4HalfSpaceCircularCylinder : public G4VHalfSpace {
public:
  G4HalfSpaceCircularCylinder();
  G4HalfSpaceCircularCylinder(const G4ThreeVector &v,
                              const G4ThreeVector &h,
                              G4double r);
  ~G4HalfSpaceCircularCylinder();

  virtual G4double Sdf(const G4ThreeVector&p) const override;
  std::vector<G4ThreeVector> Intersection(const G4ThreeVector& p, const G4ThreeVector& v) const override;

  virtual void Translate(const G4ThreeVector& t) override;
  virtual void Rotate(const G4RotationMatrix& r) override;
  virtual void Transform(const G4AffineTransform& a) override;

  virtual G4SurfaceMeshCGAL* GetSurfaceMesh()  override;


protected:
  G4ThreeVector _v;
  G4ThreeVector _h;
  G4double _r;

  G4HalfSpaceZone _hsZone = G4HalfSpaceZone();
};