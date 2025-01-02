#include "G4HalfSpaceCircularCylinder.hh"
#include "G4HalfSpaceZone.hh"
#include "G4HalfSpacePlane.hh"
#include "G4HalfSpaceQuadric.hh"
#include "G4HalfSpaceTransformation.hh"
#include "G4SurfaceMeshCGAL.hh"
#include "G4Tubs.hh"



G4HalfSpaceCircularCylinder::G4HalfSpaceCircularCylinder() {}

G4HalfSpaceCircularCylinder::G4HalfSpaceCircularCylinder(const G4ThreeVector &v,
                                                         const G4ThreeVector &h,
                                                         G4double r) :
                                                         _v(v), _h(h), _r(r) {

  G4HalfSpaceQuadric *q = new G4HalfSpaceQuadric(1/pow(_r,2), 0, 0,
                                                 1/pow(_r,2), 0,
                                                 0,
                                                 0,0,0,
                                                 -1);
  G4HalfSpacePlane *p1 = new G4HalfSpacePlane(G4ThreeVector(0,0,-1),
                                              h.mag()/2);
  G4HalfSpacePlane *p2 = new G4HalfSpacePlane(G4ThreeVector(0,0,1),
                                              h.mag()/2);
  _hsZone.AddIntersection(q);
  _hsZone.AddIntersection(p1);
  _hsZone.AddIntersection(p2);

  G4HalfSpaceTransformation t = G4HalfSpaceTransformation(_h);
  _hsZone.Rotate(t.GetRotationMatrix());
  _hsZone.Translate(_v + t.GetRotationMatrix()*G4ThreeVector(0,0,_h.mag()/2));

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

  _v += t;

  _hsZone.Translate(t);
}

void G4HalfSpaceCircularCylinder::Rotate(const G4RotationMatrix& r) {

  _v = r*_v;
  _h = r*_h;

  _hsZone.Rotate(r);
}

void G4HalfSpaceCircularCylinder::Transform(const G4AffineTransform& a) {

  Rotate(a.NetRotation());
  Translate(a.NetTranslation());

  _hsZone.Transform(a);
}

G4SurfaceMeshCGAL* G4HalfSpaceCircularCylinder::GetSurfaceMesh()  {

  G4Tubs t = G4Tubs("temp",0, _r, _h.mag()/2, 0, 2*M_PI);
  G4Polyhedron *g4poly = t.GetPolyhedron();
  G4SurfaceMeshCGAL *sm = new G4SurfaceMeshCGAL();
  sm->Fill(g4poly);

  G4HalfSpaceTransformation trans = G4HalfSpaceTransformation(_h);
  G4ThreeVector axis;
  G4double angle;
  trans.GetAxisAngle(axis,angle);

  sm->Rotate(axis,-angle);
  sm->Translate(_v+trans.GetRotationMatrix()*G4ThreeVector(0,0,_h.mag()/2.));

  return sm;
}