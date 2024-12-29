#include "G4HalfSpaceSphere.hh"
#include "G4Orb.hh"
#include "G4SurfaceMeshCGAL.hh"

G4HalfSpaceSphere::G4HalfSpaceSphere() {}

G4HalfSpaceSphere::G4HalfSpaceSphere(G4double cx, G4double cy, G4double cz, G4double radius) :
G4HalfSpaceSphere(G4ThreeVector(cx,cy,cz), radius) {}

G4HalfSpaceSphere::G4HalfSpaceSphere(G4ThreeVector centre, G4double radius) :
_centre(centre), _r(radius) {}

G4HalfSpaceSphere::~G4HalfSpaceSphere() {}

G4double G4HalfSpaceSphere::Sdf(const G4ThreeVector&p) const {
  return (p - _centre).mag() - _r;
}

std::vector<G4ThreeVector> G4HalfSpaceSphere::Intersection(const G4ThreeVector& p, const G4ThreeVector& v) const {
  // v = lambda d + p
  // |v-c| = r
  // | lambda d + p - c | = r
  // (lambda d + p - c ) . (lambda d + p -c) = r.r
  // lambda^2 d.d + 2 lambda d.(p - c) + (p - c).(p - c) = r.r
  // a = d.d
  // b = 2 d.(p - c)
  // c = (p - c).(p - c) - r*r

  std::vector<G4ThreeVector> intersections;

  auto a = v.dot(v);
  auto b = 2*v.dot(p - _centre);
  auto c = (p - _centre).dot(p - _centre ) - _r*_r;

  G4int nSoln = 0;
  G4double lambda1 = 0;
  G4double lambda2 = 0;
  G4VHalfSpace::QuadraticSolve(a,b,c,nSoln, lambda1, lambda2);
  if(lambda1 > 0)
    intersections.push_back(lambda1*v + p);
  if(lambda2 > 0)
    intersections.push_back(lambda2*v + p);

  return intersections;
}

void G4HalfSpaceSphere::Translate(const G4ThreeVector& t) {
  _centre = _centre + t;
}

void G4HalfSpaceSphere::Rotate(const G4RotationMatrix& r) {
  G4cout << "G4HalfSpaceSphere::Rotate not implemented as sphere is rotationally invariant" << G4endl;;
}

void G4HalfSpaceSphere::Transform(const G4AffineTransform& a) {
  _centre = a.NetTranslation() + _centre;
}

G4SurfaceMeshCGAL* G4HalfSpaceSphere::GetSurfaceMesh()  {
  G4Orb o = G4Orb("test",_r);
  G4Polyhedron *g4poly = o.GetPolyhedron();
  G4SurfaceMeshCGAL *sm = new G4SurfaceMeshCGAL();
  sm->Fill(g4poly);
  sm->Translate(_centre);
  return sm;
}