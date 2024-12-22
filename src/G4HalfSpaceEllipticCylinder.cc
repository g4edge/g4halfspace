#include "G4HalfSpaceEllipticCylinder.hh"
#include "G4HalfSpacePlane.hh"
#include "G4HalfSpaceQuadric.hh"

G4HalfSpaceEllipticCylinder::G4HalfSpaceEllipticCylinder() {}

G4HalfSpaceEllipticCylinder::G4HalfSpaceEllipticCylinder(const G4ThreeVector& v,
                                               const G4ThreeVector& h,
                                               G4double r1,
                                               G4double r2) :
    _v(h), _h(h), _r1(r1), _r2(r2) {
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

  G4HalfSpaceQuadric *q = new G4HalfSpaceQuadric(1/pow(_r1,2), 0, 0,
                                                 1/pow(_r2,2), 0,
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

G4HalfSpaceEllipticCylinder::~G4HalfSpaceEllipticCylinder() {}

G4double G4HalfSpaceEllipticCylinder::Sdf(const G4ThreeVector&p) const {
  return _hsZone.Sdf(p);
}
std::vector<G4ThreeVector> G4HalfSpaceEllipticCylinder::Intersection(const G4ThreeVector& p, const G4ThreeVector &v) const {
  auto intersections = _hsZone.Intersection(p,v);
  return intersections;
}

void G4HalfSpaceEllipticCylinder::Translate(const G4ThreeVector& t) {
  _hsZone.Translate(t);
}

void G4HalfSpaceEllipticCylinder::Rotate(const G4RotationMatrix& r) {
  _hsZone.Rotate(r);
}

void G4HalfSpaceEllipticCylinder::Transform(const G4AffineTransform& a) {
  _hsZone.Transform(a);
}

G4SurfaceMeshCGAL* G4HalfSpaceEllipticCylinder::GetSurfaceMesh()  {
  if(cached_mesh) {
    G4cout << "G4HalfSpaceEllipticalCylinder::GetSurfaceMesh cached" << G4endl;
    return cached_mesh;
  }

  cached_mesh = _hsZone.GetSurfaceMesh();

  return cached_mesh;}
