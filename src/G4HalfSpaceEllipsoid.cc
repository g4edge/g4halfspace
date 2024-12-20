#include "G4HalfSpaceEllipsoid.hh"

#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/Vector.h"

#include "G4SystemOfUnits.hh"
#include "G4Ellipsoid.hh"

G4HalfSpaceEllipsoid::G4HalfSpaceEllipsoid() : _a(10), _b(20), _c(30), _centre(G4ThreeVector()), _rotation(G4ThreeVector())
{}

G4HalfSpaceEllipsoid::G4HalfSpaceEllipsoid(double a,
                                           double b,
                                           double c,
                                           G4ThreeVector centre,
                                           G4ThreeVector rotation) : _a(a), _b(b), _c(c), _centre(centre), _rotation(rotation)
{}

G4double G4HalfSpaceEllipsoid::Sdf(const G4ThreeVector&p) const {

  auto x = p.x();
  auto y = p.y();
  auto z = p.z();

  auto sdf = pow(x-_centre[0],2)/pow(_a,2) + pow(y-_centre[1],2)/pow(_b,2) + pow(z-_centre[2],2)/pow(_c,2) - 1;
  return sdf;
}

std::vector<G4ThreeVector> G4HalfSpaceEllipsoid::Intersection(const G4ThreeVector& p, const G4ThreeVector &d) const {

  auto q = CLHEP::HepMatrix(3,3,0);
  auto qp = CLHEP::HepVector(3);
  auto n = CLHEP::HepVector(3);
  auto p0 = CLHEP::HepVector(3);

  n(1) = d.x();
  n(2) = d.y();
  n(3) = d.z();

  p0(1) = p.x();
  p0(2) = p.y();
  p0(3) = p.z();

  q(1,1) = 1./pow(_a,2);
  q(2,2) = 1./pow(_b,2);
  q(3,3) = 1./pow(_c,2);

  auto a = (n.T()*q*n)(1,1);
  auto b = (n.T()*q*p0 + qp.T()*n + p0.T()*q*n)(1,1);
  auto c = (qp.T()*p0 + p0.T()*q*p0)(1,1) - 1 ;

  G4int nSoln;
  G4double l1, l2;
  QuadraticSolve(a,b,c,nSoln,l1,l2);

  auto inter = std::vector<G4ThreeVector>();
  if (l1 >= 0 ) {
    inter.push_back(l1*d+p);
  }
  if (l2 >= 0 ) {
    inter.push_back(l2*d+p);
  }

  return inter;
}

void G4HalfSpaceEllipsoid::Translate(const G4ThreeVector& t) {
}

void G4HalfSpaceEllipsoid::Rotate(const G4RotationMatrix& r) {
  G4cout << "G4HalfSpaceEllipsoid::Rotate> Not implemented" << G4endl;
}
void G4HalfSpaceEllipsoid::Transform(const G4AffineTransform& a) {
  G4cout << "G4HalfSpaceEllipsoid::Transform> Not implemented" << G4endl;
}

G4SurfaceMeshCGAL* G4HalfSpaceEllipsoid::GetSurfaceMesh() {
  G4Ellipsoid t = G4Ellipsoid("test",_a,_b,_c);
  G4Polyhedron *g4poly = t.GetPolyhedron();
  G4SurfaceMeshCGAL *sm = new G4SurfaceMeshCGAL();
  sm->Fill(g4poly);
  sm->Translate(_centre[0], _centre[1], _centre[2]);

  return sm;
}
