#include "G4HalfSpaceAARBox.hh"
#include "G4HalfSpaceZone.hh"
#include "G4HalfSpacePlane.hh"
#include "G4SurfaceMeshCGAL.hh"
#include "G4Box.hh"

G4HalfSpaceAARBox::G4HalfSpaceAARBox() {}
G4HalfSpaceAARBox::G4HalfSpaceAARBox(G4double xmin, G4double xmax,
                                     G4double ymin, G4double ymax,
                                     G4double zmin, G4double zmax) :
                                     _xmin(xmin), _xmax(xmax),
                                     _ymin(ymin), _ymax(ymax),
                                     _zmin(zmin), _zmax(zmax) {

  auto p1 = new G4HalfSpacePlane(G4ThreeVector(1,0,0),G4ThreeVector(xmax,0,0));
  auto p2 = new G4HalfSpacePlane(G4ThreeVector(-1,0,0),G4ThreeVector(xmin,0,0));
  auto p3 = new G4HalfSpacePlane(G4ThreeVector(0,1,0),G4ThreeVector(0,ymax,0));
  auto p4 = new G4HalfSpacePlane(G4ThreeVector(0,-1,0),G4ThreeVector(0,ymin,0));
  auto p5 = new G4HalfSpacePlane(G4ThreeVector(0,0,1),G4ThreeVector(0,0,zmax));
  auto p6 = new G4HalfSpacePlane(G4ThreeVector(0,0,-1),G4ThreeVector(0,0,zmin));

  _hsZone.AddIntersection(p1);
  _hsZone.AddIntersection(p2);
  _hsZone.AddIntersection(p3);
  _hsZone.AddIntersection(p4);
  _hsZone.AddIntersection(p5);
  _hsZone.AddIntersection(p6);
}

G4HalfSpaceAARBox::G4HalfSpaceAARBox(const G4ThreeVector &d, const G4ThreeVector &c) :
  G4HalfSpaceAARBox(-d.x()+c.x(), d.x()+c.x(),
                    -d.y()+c.y(), d.y()+c.y(),
                    -d.z()+c.z(), d.z()+c.z()) {}

G4HalfSpaceAARBox::~G4HalfSpaceAARBox() {
}

G4double G4HalfSpaceAARBox::Sdf(const G4ThreeVector&p) const {
  return _hsZone.Sdf(p);
}

std::vector<G4ThreeVector> G4HalfSpaceAARBox::Intersection(const G4ThreeVector& p, const G4ThreeVector& v) const {
  auto intersections = _hsZone.Intersection(p,v);
  return intersections;
}

void G4HalfSpaceAARBox::Translate(const G4ThreeVector& t) {

  /* translate member variables */
  _xmin += t.x();
  _xmax += t.x();
  _ymin += t.y();
  _ymin += t.y();
  _zmin += t.z();
  _zmin += t.z();

  /* translate hs zone */
  _hsZone.Translate(t);
}

void G4HalfSpaceAARBox::Rotate(const G4RotationMatrix& r) {
  G4cout << "G4HalfSpaceAARBox::Rotate not implemented use RBox" << G4endl;;
}

void G4HalfSpaceAARBox::Transform(const G4AffineTransform& a) {
  G4cout << "G4HalfSpaceAARBox::Transform not implemented use RBox" << G4endl;;
}

G4SurfaceMeshCGAL* G4HalfSpaceAARBox::GetSurfaceMesh()  {
  if(cached_mesh) {
    G4cout << "G4HalfSpaceAARBox::GetSurfaceMesh cached" << G4endl;
    return cached_mesh;
  }

  cached_mesh = _hsZone.GetSurfaceMesh();

  return cached_mesh;
}