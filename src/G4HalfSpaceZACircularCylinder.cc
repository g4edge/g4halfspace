#include "G4HalfSpaceZACircularCylinder.hh"

#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/Vector.h"

#include "G4SystemOfUnits.hh"
#include "G4Tubs.hh"

G4HalfSpaceZACircularCylinder::G4HalfSpaceZACircularCylinder() {}

G4HalfSpaceZACircularCylinder::G4HalfSpaceZACircularCylinder(G4double x,
                                                             G4double y,
                                                             G4double r) :
                                                             _x0(x), _y0(y), _r(r) {}

G4HalfSpaceZACircularCylinder::~G4HalfSpaceZACircularCylinder() {}

G4double G4HalfSpaceZACircularCylinder::Sdf(const G4ThreeVector&p) const {
  return sqrt(pow(p.x() - _x0,2) + pow(p.y() - _y0,2) ) - _r;
}

std::vector<G4ThreeVector> G4HalfSpaceZACircularCylinder::Intersection(const G4ThreeVector& p, const G4ThreeVector &d) const {

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

  q(1,1) = 1./pow(_r,2);
  q(2,2) = 1./pow(_r,2);

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

void G4HalfSpaceZACircularCylinder::Translate(const G4ThreeVector& t) {
  _x0 += t.x();
  _y0 += t.y();
}

void G4HalfSpaceZACircularCylinder::Rotate(const G4RotationMatrix& r) {
  G4cout << "G4HalfSpaceZACircularCylinder::Rotate> Not implemented" << G4endl;
}
void G4HalfSpaceZACircularCylinder::Transform(const G4AffineTransform& a) {
  G4cout << "G4HalfSpaceZACircularCylinder::Transform> Not implemented" << G4endl;
}

G4SurfaceMeshCGAL* G4HalfSpaceZACircularCylinder::GetSurfaceMesh() {
  G4Tubs t = G4Tubs("test",0,_r,1000000000,0,2*M_PI*rad);
  G4Polyhedron *g4poly = t.GetPolyhedron();
  G4SurfaceMeshCGAL *sm = new G4SurfaceMeshCGAL();
  sm->Fill(g4poly);
  sm->Translate(_x0, _y0, 0 );

  return sm;
}
