#include "G4HalfSpaceWedge.hh"
#include "G4HalfSpaceZone.hh"
#include "G4HalfSpacePlane.hh"
#include "G4SurfaceMeshCGAL.hh"
#include "G4Box.hh"

G4HalfSpaceWedge::G4HalfSpaceWedge() {}

G4HalfSpaceWedge::G4HalfSpaceWedge(std::array<double,3>  v,
                                   std::array<double,3>  h1,
                                   std::array<double,3>  h2,
                                   std::array<double,3>  h3) :
    _v(v[0], v[1], v[2]),
    _h1(h1[0], h1[1], v[2]),
    _h2(h2[0], h2[1], h2[2]),
    _h3(h3[0], h3[1], h3[2])
{
  this->ComputeSurfaces();
}

G4HalfSpaceWedge::G4HalfSpaceWedge(G4ThreeVector  v,
                                   G4ThreeVector  h1,
                                   G4ThreeVector  h2,
                                   G4ThreeVector  h3) :
    _v(v),
    _h1(h1),
    _h2(h2),
    _h3(h3)
{
  this->ComputeSurfaces();
}

G4HalfSpaceWedge::G4HalfSpaceWedge(G4ThreeVector h,
                                   G4ThreeVector c,
                                   G4ThreeVector r) {

  _v = c - h/2;
  _h1 = G4ThreeVector(h[0],0,0);
  _h2 = G4ThreeVector(0,h[1],0);
  _h3 = G4ThreeVector(0,0,h[2]);

  G4RotationMatrix rotation_g4 = G4RotationMatrix();
  rotation_g4.rotateX(r[0]);
  rotation_g4.rotateY(r[1]);
  rotation_g4.rotateZ(r[2]);
  rotation_g4.rectify();

  _h1 = rotation_g4*_h1;
  _h2 = rotation_g4*_h2;
  _h3 = rotation_g4*_h3;

  this->ComputeSurfaces();
}

void G4HalfSpaceWedge::ComputeSurfaces() {
  auto p1 = _v + 0.5*_h1 + 0.5*_h2;
  auto p2 = _v + 0.5*_h1 + 0.5*_h3;
  auto p3 = _v + 0.5*_h2 + 0.5*_h3;
  auto p4 = p2 + _h2;
  auto p5 = _v+_h1;

  auto n1 = _h2.cross(_h1);
  auto n2 = _h1.cross(_h3);
  auto n3 = _h3.cross(_h2);
  auto n4 = _h3.cross(_h1);
  auto n5 = _h2.cross(_h3-_h1);

  n1 = n1/n1.mag();
  n2 = n2/n2.mag();
  n3 = n3/n3.mag();
  n4 = n4/n4.mag();
  n5 = n5/n5.mag();

  _hsZone = new G4HalfSpaceZone();
  auto pl1 = new G4HalfSpacePlane(p1,n1);
  auto pl2 = new G4HalfSpacePlane(p2,n2);
  auto pl3 = new G4HalfSpacePlane(p3,n3);
  auto pl4 = new G4HalfSpacePlane(p4,n4);
  auto pl5 = new G4HalfSpacePlane(p5,n5);

  _hsZone->AddIntersection(pl1);
  _hsZone->AddIntersection(pl2);
  _hsZone->AddIntersection(pl3);
  _hsZone->AddIntersection(pl4);
  _hsZone->AddIntersection(pl5);
}

G4HalfSpaceWedge::~G4HalfSpaceWedge() {
  delete _hsZone;
}

G4double G4HalfSpaceWedge::Sdf(const G4ThreeVector&p) const {
  return _hsZone->Sdf(p);
}

std::vector<G4ThreeVector> G4HalfSpaceWedge::Intersection(const G4ThreeVector& p, const G4ThreeVector& v) const {
  auto intersections = _hsZone->Intersection(p,v);
  return intersections;
}

void G4HalfSpaceWedge::Translate(const G4ThreeVector& t) {
  _v += t;

  _hsZone->Translate(t);
}

void G4HalfSpaceWedge::Rotate(const G4RotationMatrix& r) {
  _h1 = r*_h1;
  _h2 = r*_h2;
  _h3 = r*_h3;

  _hsZone->Rotate(r);
}

void G4HalfSpaceWedge::Transform(const G4AffineTransform& a) {
  _hsZone->Transform(a);
}

G4SurfaceMeshCGAL* G4HalfSpaceWedge::GetSurfaceMesh()  {
  if(cached_mesh) {
    G4cout << "G4HalfSpaceWedge::GetSurfaceMesh cached" << G4endl;
    return cached_mesh;
  }

  cached_mesh = _hsZone->GetSurfaceMesh();

  return cached_mesh;
}