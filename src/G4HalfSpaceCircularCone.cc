#include "G4HalfSpaceCircularCone.hh"
#include "G4HalfSpaceZone.hh"
#include "G4HalfSpacePlane.hh"
#include "G4HalfSpaceQuadric.hh"
#include "G4SurfaceMeshCGAL.hh"
#include "G4HalfSpaceTransformation.hh"
#include "G4Cons.hh"

G4HalfSpaceCircularCone::G4HalfSpaceCircularCone() {
}

G4HalfSpaceCircularCone::G4HalfSpaceCircularCone(const G4ThreeVector &v,
                                                 const G4ThreeVector &h,
                                                 G4double r1, G4double r2) :
                                                 _v(v), _h(h), _r1(r1), _r2(r2) {
  ComputeSurfaces();
}


G4HalfSpaceCircularCone::G4HalfSpaceCircularCone(G4double h, G4double r1, G4double r2,
                                                 const G4ThreeVector &c, const G4ThreeVector &r)
{
  _r1 = r1;
  _r2 = r2;
  auto t = G4HalfSpaceTransformation(c, r);
  _h = t.GetRotationMatrix()*G4ThreeVector(0,0,h);
  _v = c + t.GetRotationMatrix()*G4ThreeVector(0,0,-h/2);
  ComputeSurfaces();
}

void G4HalfSpaceCircularCone::ComputeSurfaces(){

  auto r1temp = _r1;
  auto r2temp = _r2;

  if (_r1 > _r2) {
    r1temp = _r2;
    r2temp = _r1;
  }

  auto g = (r2temp - r1temp) / _h.mag();

  auto z1 = r1temp / g;
  auto z2 = r2temp / g;

  auto c = 1;
  auto a = (r2temp - r1temp) / _h.mag();

  G4HalfSpaceQuadric *q = new G4HalfSpaceQuadric(1 / pow(a, 2), 0, 0,
                                                 1 / pow(a, 2), 0,
                                                 -1 / pow(c, 2),
                                                 0, 0, 0,
                                                 0);

  G4HalfSpacePlane *p1 = new G4HalfSpacePlane(G4ThreeVector(0, 0, z1),
                                              G4ThreeVector(0, 0, -1));
  G4HalfSpacePlane *p2 = new G4HalfSpacePlane(G4ThreeVector(0, 0, z2),
                                              G4ThreeVector(0, 0, 1));


  _hsZone.AddIntersection(q);
  _hsZone.AddIntersection(p1);
  _hsZone.AddIntersection(p2);

  G4ThreeVector shift = G4ThreeVector(0,0,-(z2+z1)/2.0);
  _hsZone.Translate(shift);

  if(_r1 > _r2) {
    G4RotationMatrix rot;
    rot.rotateX(M_PI);
    rot.rectify();
    _hsZone.Rotate(rot);
  }

  G4HalfSpaceTransformation t = G4HalfSpaceTransformation(_h);
  _hsZone.Rotate(t.GetRotationMatrix());

  _hsZone.Translate(_v+t.GetRotationMatrix()*G4ThreeVector(0,0,(z1+z2)/2-z1));
}

G4HalfSpaceCircularCone::~G4HalfSpaceCircularCone() {
}

G4double G4HalfSpaceCircularCone::Sdf(const G4ThreeVector&p) const {
  return _hsZone.Sdf(p);
}

std::vector<G4ThreeVector> G4HalfSpaceCircularCone::Intersection(const G4ThreeVector& p, const G4ThreeVector& v) const {
  auto intersections = _hsZone.Intersection(p,v);
  return intersections;
}

void G4HalfSpaceCircularCone::Translate(const G4ThreeVector& t) {
  _v += t;
  _hsZone.Translate(t);
}

void G4HalfSpaceCircularCone::Rotate(const G4RotationMatrix& r) {
  _h = r*_h;
  _v = r*_v;
  _hsZone.Rotate(r);
}

void G4HalfSpaceCircularCone::Transform(const G4AffineTransform& a) {
  _hsZone.Transform(a);
}

G4SurfaceMeshCGAL* G4HalfSpaceCircularCone::GetSurfaceMesh()  {


  double angle;
  G4ThreeVector axis;
  G4HalfSpaceTransformation t = G4HalfSpaceTransformation(_h);
  t.GetAxisAngle(axis, angle);


  G4Cons c = G4Cons("temp",0, _r1, 0, _r2, _h.mag()/2, 0, 2*M_PI);
  G4Polyhedron *g4poly = c.GetPolyhedron();
  G4SurfaceMeshCGAL *sm = new G4SurfaceMeshCGAL();
  sm->Fill(g4poly);
  sm->Rotate(axis, -angle);
  sm->Translate(_v+t.GetRotationMatrix()*G4ThreeVector(0,0,_h.mag()/2));

  return sm;
}