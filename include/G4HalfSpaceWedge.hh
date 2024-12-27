#pragma once

#include "G4ThreeVector.hh"
#include "G4VHalfSpace.hh"

class G4HalfSpaceZone;

class G4HalfSpaceWedge : public G4VHalfSpace {
public:
  G4HalfSpaceWedge();
  G4HalfSpaceWedge(std::array<double,3> v,
                  std::array<double,3> h1,
                  std::array<double,3> h2,
                  std::array<double,3> h3);

  G4HalfSpaceWedge(G4ThreeVector v,
                  G4ThreeVector h1,
                  G4ThreeVector h2,
                  G4ThreeVector h3);

  G4HalfSpaceWedge(G4ThreeVector v,
                  G4ThreeVector h,
                  G4ThreeVector r);

  void ComputePlanes();

  ~G4HalfSpaceWedge();

  virtual G4double Sdf(const G4ThreeVector&p) const override;
  std::vector<G4ThreeVector> Intersection(const G4ThreeVector& p, const G4ThreeVector& v) const override;

  virtual void Translate(const G4ThreeVector& t) override;
  virtual void Rotate(const G4RotationMatrix& r) override;
  virtual void Transform(const G4AffineTransform& a) override;

  virtual G4SurfaceMeshCGAL* GetSurfaceMesh()  override;

protected:
  G4ThreeVector _v;
  G4ThreeVector _h1;
  G4ThreeVector _h2;
  G4ThreeVector _h3;

  G4HalfSpaceZone* _hsZone;
};