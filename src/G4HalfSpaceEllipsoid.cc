#include "G4HalfSpaceEllipsoid.hh"
#include "G4HalfSpaceQuadric.hh"
#include "G4HalfSpaceTransformation.hh"

#include "G4SystemOfUnits.hh"
#include "G4Ellipsoid.hh"

G4HalfSpaceEllipsoid::G4HalfSpaceEllipsoid() :
                      _f1(G4ThreeVector(0,0,-25)),
                      _f2(G4ThreeVector(0,0,25)),
                      _l(75)
{
  ComputeSurfaces();
}

G4HalfSpaceEllipsoid::G4HalfSpaceEllipsoid(const G4ThreeVector &f1,
                                           const G4ThreeVector &f2,
                                           G4double l) :
    _f1(f1),
    _f2(f2),
    _l(l)
{
  ComputeSurfaces();
}

G4HalfSpaceEllipsoid::G4HalfSpaceEllipsoid(const G4ThreeVector &abc,
                                           const G4ThreeVector &centre,
                                           const G4ThreeVector &rotation)
{
  auto c = sqrt(pow(abc.z(),2) - pow(abc.x(),2));
  _f1 = G4ThreeVector(0,0,-c);
  _f2 = G4ThreeVector(0,0, c);
  _l = 2*abc.z();

  auto t = G4HalfSpaceTransformation(centre, rotation);
  _f1 = t.GetRotationMatrix()*_f1;
  _f2 = t.GetRotationMatrix()*_f2;

  _f1 += centre;
  _f2 += centre;

  ComputeSurfaces();
}

void G4HalfSpaceEllipsoid::ComputeSurfaces() {

  _c = (_f1 - (_f1 + _f2)/2).mag();
  _r = sqrt(pow(_l/2,2) - pow(_c,2));

  auto transformation = G4HalfSpaceTransformation(_f2-_f1);
  auto centre = (_f2+_f1)/2.0;
  G4HalfSpaceQuadric *q = new G4HalfSpaceQuadric(1/pow(_r,2), 0, 0,
                                                 1/pow(_r,2), 0,
                                                 1/pow(_l/2,2),
                                                 0,0,0,
                                                 -1);
  q->Rotate(transformation.GetRotationMatrix());
  q->Translate(centre);
  _hsZone.AddIntersection(q);
}

G4double G4HalfSpaceEllipsoid::Sdf(const G4ThreeVector&p) const {
  return _hsZone.Sdf(p);
}

std::vector<G4ThreeVector> G4HalfSpaceEllipsoid::Intersection(const G4ThreeVector& p, const G4ThreeVector &v) const {
  auto intersections = _hsZone.Intersection(p,v);
  return intersections;
}

void G4HalfSpaceEllipsoid::Translate(const G4ThreeVector& t) {
  _f1 += t;
  _f2 += t;

  _hsZone.Translate(t);
}

void G4HalfSpaceEllipsoid::Rotate(const G4RotationMatrix& r) {
  _f1 = r*_f1;
  _f2 = r*_f2;

  _hsZone.Rotate(r);
}
void G4HalfSpaceEllipsoid::Transform(const G4AffineTransform& a) {
  _hsZone.Transform(a);
}

G4SurfaceMeshCGAL* G4HalfSpaceEllipsoid::GetSurfaceMesh() {

  G4Ellipsoid t = G4Ellipsoid("test",_r,_r,_l/2);
  G4Polyhedron *g4poly = t.GetPolyhedron();
  G4SurfaceMeshCGAL *sm = new G4SurfaceMeshCGAL();
  sm->Fill(g4poly);

  G4ThreeVector axis;
  G4double angle;

  auto transformation = G4HalfSpaceTransformation(_f2-_f1);
  transformation.GetAxisAngle(axis, angle);
  auto centre = (_f2+_f1)/2.0;

  sm->Rotate(axis,-angle);
  sm->Translate(centre);

  return sm;
}
