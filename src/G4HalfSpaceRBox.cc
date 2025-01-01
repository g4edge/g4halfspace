#include "G4HalfSpaceRBox.hh"
#include "G4HalfSpaceZone.hh"
#include "G4HalfSpacePlane.hh"
#include "G4SurfaceMeshCGAL.hh"
#include "G4Box.hh"

G4HalfSpaceRBox::G4HalfSpaceRBox() {}

G4HalfSpaceRBox::G4HalfSpaceRBox(std::array<double,3>  v,
                                 std::array<double,3>  h1,
                                 std::array<double,3>  h2,
                                 std::array<double,3>  h3) :
                                 _v(v[0], v[1], v[2]),
                                 _h1(h1[0], h1[1], v[2]),
                                 _h2(h2[0], h2[1], h2[2]),
                                 _h3(h3[0], h3[1], h3[2])
{
  this->ComputePlanes();
}

G4HalfSpaceRBox::G4HalfSpaceRBox(G4ThreeVector  v,
                                 G4ThreeVector  h1,
                                 G4ThreeVector  h2,
                                 G4ThreeVector  h3) :
    _v(v),
    _h1(h1),
    _h2(h2),
    _h3(h3)
{
  this->ComputePlanes();
}

G4HalfSpaceRBox::G4HalfSpaceRBox(G4ThreeVector d,
                                 G4ThreeVector c,
                                 G4ThreeVector r) :
                                 _v(c-d) {
  _h1 = G4ThreeVector(2*d[0],0,0);
  _h2 = G4ThreeVector(0,2*d[1],0);
  _h3 = G4ThreeVector(0,0,2*d[2]);

  G4RotationMatrix rotation_g4 = G4RotationMatrix();
  rotation_g4.rotateX(r[0]);
  rotation_g4.rotateY(r[1]);
  rotation_g4.rotateZ(r[2]);
  rotation_g4.rectify();

  _h1 = rotation_g4*_h1;
  _h2 = rotation_g4*_h2;
  _h3 = rotation_g4*_h3;

  this->ComputePlanes();
}

void G4HalfSpaceRBox::ComputePlanes() {
  auto p1 = _v + 0.5*_h1 + 0.5*_h2;
  auto p2 = _v + 0.5*_h1 + 0.5*_h3;
  auto p3 = _v + 0.5*_h2 + 0.5*_h3;
  auto p4 = p1 + _h3;
  auto p5 = p2 + _h2;
  auto p6 = p3 + _h1;

  auto n1 = _h2.cross(_h1);
  auto n2 = _h1.cross(_h3);
  auto n3 = _h3.cross(_h2);
  auto n4 = _h1.cross(_h2);
  auto n5 = _h3.cross(_h1);
  auto n6 = _h2.cross(_h3);

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
  auto pl6 = new G4HalfSpacePlane(p6,n6);

  _hsZone->AddIntersection(pl1);
  _hsZone->AddIntersection(pl2);
  _hsZone->AddIntersection(pl3);
  _hsZone->AddIntersection(pl4);
  _hsZone->AddIntersection(pl5);
  _hsZone->AddIntersection(pl6);

}

G4HalfSpaceRBox::~G4HalfSpaceRBox() {
  delete _hsZone;
}

G4double G4HalfSpaceRBox::Sdf(const G4ThreeVector&p) const {
  return _hsZone->Sdf(p);
}

std::vector<G4ThreeVector> G4HalfSpaceRBox::Intersection(const G4ThreeVector& p, const G4ThreeVector& v) const {
  auto intersections = _hsZone->Intersection(p,v);
  return intersections;
}

void G4HalfSpaceRBox::Translate(const G4ThreeVector& t) {
  _hsZone->Translate(t);
}

void G4HalfSpaceRBox::Rotate(const G4RotationMatrix& r) {
  _hsZone->Rotate(r);
}

void G4HalfSpaceRBox::Transform(const G4AffineTransform& a) {
  _hsZone->Transform(a);
}

G4SurfaceMeshCGAL* G4HalfSpaceRBox::GetSurfaceMesh()  {
  if(cached_mesh) {
    G4cout << "G4HalfSpaceRBox::GetSurfaceMesh cached" << G4endl;
    return cached_mesh;
  }

  cached_mesh = _hsZone->GetSurfaceMesh();

  return cached_mesh;
}