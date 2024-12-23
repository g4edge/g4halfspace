//
// Created by Boogert Stewart on 20/12/2024.
//

#pragma once

#include "G4ThreeVector.hh"
#include "G4VHalfSpace.hh"
#include "G4HalfSpaceZone.hh"

class G4HalfSpaceEllipsoid : public G4VHalfSpace {
public:
  G4HalfSpaceEllipsoid();
  G4HalfSpaceEllipsoid(double a, double b, double c, G4ThreeVector centre, G4ThreeVector rotation);
  ~G4HalfSpaceEllipsoid() = default;

  virtual G4double Sdf(const G4ThreeVector&p) const override;
  std::vector<G4ThreeVector> Intersection(const G4ThreeVector& p, const G4ThreeVector& v) const override;

  virtual void Translate(const G4ThreeVector& t) override;
  virtual void Rotate(const G4RotationMatrix& r) override;
  virtual void Transform(const G4AffineTransform& a) override;

  virtual G4SurfaceMeshCGAL* GetSurfaceMesh()  override;


protected:
  double _a;
  double _b;
  double _c;
  G4ThreeVector _centre;
  G4ThreeVector _rotation;

  G4HalfSpaceZone _hsZone = G4HalfSpaceZone();

};
