#include "G4HalfSpaceEllipsoid.hh"

#include "G4HalfSpaceQuadric.hh"

#include "G4SystemOfUnits.hh"
#include "G4Ellipsoid.hh"

G4HalfSpaceEllipsoid::G4HalfSpaceEllipsoid() : _a(10), _b(20), _c(30), _centre(G4ThreeVector()), _rotation(G4ThreeVector())
{}

G4HalfSpaceEllipsoid::G4HalfSpaceEllipsoid(double a,
                                           double b,
                                           double c,
                                           G4ThreeVector centre,
                                           G4ThreeVector rotation) :
                                           _a(a), _b(b), _c(c),
                                           _centre(centre),
                                           _rotation(rotation)
{
  G4HalfSpaceQuadric *q = new G4HalfSpaceQuadric(1/pow(_a,2), 0, 0,
                                                 1/pow(_b,2), 0,
                                                 1/pow(_c,2),
                                                 0,0,0,
                                                 -1);
  q->Rotate(rotation);
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
  _hsZone.Translate(t);
}

void G4HalfSpaceEllipsoid::Rotate(const G4RotationMatrix& r) {
  _hsZone.Rotate(r);
}
void G4HalfSpaceEllipsoid::Transform(const G4AffineTransform& a) {
  _hsZone.Transform(a);
}

G4SurfaceMeshCGAL* G4HalfSpaceEllipsoid::GetSurfaceMesh() {
  G4Ellipsoid t = G4Ellipsoid("test",_a,_b,_c);
  G4Polyhedron *g4poly = t.GetPolyhedron();
  G4SurfaceMeshCGAL *sm = new G4SurfaceMeshCGAL();
  sm->Fill(g4poly);

  G4RotationMatrix rotationMatrix = G4RotationMatrix();
  rotationMatrix.rotateX(_rotation[0]);
  rotationMatrix.rotateY(_rotation[1]);
  rotationMatrix.rotateZ(_rotation[2]);
  rotationMatrix.rectify();

  double angle;
  G4ThreeVector axis;
  rotationMatrix.getAngleAxis(angle, axis);

  sm->Rotate(axis, angle);
  sm->Translate(_centre[0], _centre[1], _centre[2]);

  return sm;
}
