#include "G4HalfSpaceEllipticCylinder.hh"
#include "G4HalfSpacePlane.hh"
#include "G4HalfSpaceQuadric.hh"
#include "G4HalfSpaceTransformation.hh"
#include "G4EllipticalTube.hh"

G4HalfSpaceEllipticCylinder::G4HalfSpaceEllipticCylinder() {}

G4HalfSpaceEllipticCylinder::G4HalfSpaceEllipticCylinder(const G4ThreeVector& v,
                                                         const G4ThreeVector& h,
                                                         const G4ThreeVector& r1,
                                                         const G4ThreeVector& r2) :
    _v(v), _h(h), _r1(r1), _r2(r2) {

  auto hmag = h.mag();

  G4HalfSpaceTransformation trans = G4HalfSpaceTransformation(_r1.unit(), _r2.unit(), _h.unit());


  G4HalfSpaceQuadric *q = new G4HalfSpaceQuadric(1/pow(_r1.mag(),2), 0, 0,
                                                 1/pow(_r2.mag(),2), 0,
                                                 0,
                                                 0,0,0,
                                                 -1);
  G4HalfSpacePlane *p1 = new G4HalfSpacePlane(G4ThreeVector(0,0,-1),
                                              hmag/2);
  G4HalfSpacePlane *p2 = new G4HalfSpacePlane(G4ThreeVector(0,0,1),
                                              hmag/2);

  _hsZone.AddIntersection(q);
  _hsZone.AddIntersection(p1);
  _hsZone.AddIntersection(p2);

  _hsZone.Rotate(trans.GetRotationMatrix());
  _hsZone.Translate(_v+trans.GetRotationMatrix()*G4ThreeVector(0,0,hmag/2));
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

  _v += t;

  _hsZone.Translate(t);
}

void G4HalfSpaceEllipticCylinder::Rotate(const G4RotationMatrix& r) {

  _v = r*_v;
  _h = r*_h;
  _r1 = r*_r1;
  _r2 = r*_r2;

  _hsZone.Rotate(r);
}

void G4HalfSpaceEllipticCylinder::Transform(const G4AffineTransform& a) {

  Rotate(a.NetRotation());
  Translate(a.NetTranslation());

  _hsZone.Transform(a);
}

G4SurfaceMeshCGAL* G4HalfSpaceEllipticCylinder::GetSurfaceMesh()  {
  G4EllipticalTube t = G4EllipticalTube("temp",_r1.mag(),_r2.mag(),_h.mag()/2);
  G4Polyhedron *g4poly = t.GetPolyhedron();
  G4SurfaceMeshCGAL *sm = new G4SurfaceMeshCGAL();
  sm->Fill(g4poly);

  G4HalfSpaceTransformation trans = G4HalfSpaceTransformation(_r1.unit(), _r2.unit(), _h.unit());

  sm->Rotate(trans.GetRotationMatrix());
  sm->Translate(_v+trans.GetRotationMatrix()*G4ThreeVector(0,0,_h.mag()/2.));
  return sm;

}
