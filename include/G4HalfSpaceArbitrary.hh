#pragma once

#include "G4ThreeVector.hh"
#include "G4VHalfSpace.hh"
#include "G4HalfSpaceZone.hh"

#include <vector>


class G4HalfSpaceArbitrary : public G4VHalfSpace {
public:
  G4HalfSpaceArbitrary();
  G4HalfSpaceArbitrary(const G4ThreeVector &v1,
                       const G4ThreeVector &v2,
                       const G4ThreeVector &v3,
                       const G4ThreeVector &v4,
                       const G4ThreeVector &v5,
                       const G4ThreeVector &v6,
                       const G4ThreeVector &v7,
                       const G4ThreeVector &v8,
                       G4int fi1,
                       G4int fi2,
                       G4int fi3,
                       G4int fi4,
                       G4int fi5,
                       G4int fi6);

  G4HalfSpaceArbitrary(const G4ThreeVector &c,
                       const G4ThreeVector &v1,
                       const G4ThreeVector &v2,
                       const G4ThreeVector &v3,
                       const G4ThreeVector &v4,
                       const G4ThreeVector &v5,
                       const G4ThreeVector &v6,
                       const G4ThreeVector &v7,
                       const G4ThreeVector &v8,
                       G4int fi1,
                       G4int fi2,
                       G4int fi3,
                       G4int fi4,
                       G4int fi5,
                       G4int fi6);

  void ComputePlanes();

  ~G4HalfSpaceArbitrary();

  virtual G4double Sdf(const G4ThreeVector&p) const override;
  std::vector<G4ThreeVector> Intersection(const G4ThreeVector& p, const G4ThreeVector& v) const override;

  virtual void Translate(const G4ThreeVector& t) override;
  virtual void Rotate(const G4RotationMatrix& r) override;
  virtual void Transform(const G4AffineTransform& a) override;

  virtual G4SurfaceMeshCGAL* GetSurfaceMesh()  override;

  std::string ToString();
protected:

  std::vector<G4ThreeVector> _verts;
  std::vector<G4int> _faceVertCount;
  std::vector<G4int> _faceVertIndex;

  G4HalfSpaceZone _hsZone = G4HalfSpaceZone();

private:
  std::vector<G4int> decodeFaceInteger(G4int flukaFaceDef);
};

