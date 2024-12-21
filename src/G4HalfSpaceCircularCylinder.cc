#include "G4HalfSpaceCircularCylinder.hh"
#include "G4HalfSpaceZone.hh"
#include "G4HalfSpacePlane.hh"
#include "G4HalfSpaceQuadric.hh"
#include "G4SurfaceMeshCGAL.hh"


G4HalfSpaceCircularCylinder::G4HalfSpaceCircularCylinder() {}

G4HalfSpaceCircularCylinder::G4HalfSpaceCircularCylinder(const G4ThreeVector &v,
                                                         const G4ThreeVector &h,
                                                         G4double r) :
                                                         _v(v), _h(h), _r(r) {

  auto hmag = h.mag();

  G4RotationMatrix rmatrix;

  G4ThreeVector zaxis = G4ThreeVector(0,0,1);
  G4ThreeVector hnorm = _h/_h.mag();
  G4ThreeVector raxis = zaxis.cross(hnorm);
  G4double angle = asin(raxis.mag());

  if(angle ==0) {
    raxis = G4ThreeVector(0,0,1);
  }
  rmatrix.set(raxis, angle);

  G4HalfSpaceQuadric *q = new G4HalfSpaceQuadric(1/pow(_r,2), 0, 0,
                                                 1/pow(_r,2), 0,
                                                 0,
                                                 0,0,0,
                                                 -1);
  G4HalfSpacePlane *p1 = new G4HalfSpacePlane(G4ThreeVector(0,0,-1),
                                              hmag/2);
  G4HalfSpacePlane *p2 = new G4HalfSpacePlane(G4ThreeVector(0,0,1),
                                              hmag/2);
  q->Rotate(rmatrix);
  p1->Rotate(rmatrix);
  p2->Rotate(rmatrix);

  _hsZone.AddIntersection(q);
  _hsZone.AddIntersection(p1);
  _hsZone.AddIntersection(p2);
}

G4HalfSpaceCircularCylinder::~G4HalfSpaceCircularCylinder() {
}

G4double G4HalfSpaceCircularCylinder::Sdf(const G4ThreeVector&p) const {
  return _hsZone.Sdf(p);
}

std::vector<G4ThreeVector> G4HalfSpaceCircularCylinder::Intersection(const G4ThreeVector& p, const G4ThreeVector& v) const {
  auto intersections = _hsZone.Intersection(p,v);
  return intersections;
}

void G4HalfSpaceCircularCylinder::Translate(const G4ThreeVector& t) {
  _hsZone.Translate(t);
}

void G4HalfSpaceCircularCylinder::Rotate(const G4RotationMatrix& r) {
  _hsZone.Rotate(r);
}

void G4HalfSpaceCircularCylinder::Transform(const G4AffineTransform& a) {
  _hsZone.Transform(a);
}

G4SurfaceMeshCGAL* G4HalfSpaceCircularCylinder::GetSurfaceMesh()  {
  if(cached_mesh) {
    G4cout << "G4HalfSpaceCircularCylinder::GetSurfaceMesh cached" << G4endl;
    return cached_mesh;
  }

  cached_mesh = _hsZone.GetSurfaceMesh();

  return cached_mesh;
}