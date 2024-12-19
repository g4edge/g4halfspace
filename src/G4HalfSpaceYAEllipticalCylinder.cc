#include "G4HalfSpaceYAEllipticalCylinder.hh"

#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/Vector.h"

#include "G4SystemOfUnits.hh"
#include "G4Tubs.hh"
#include "G4EllipticalTube.hh"

G4HalfSpaceYAEllipticalCylinder::G4HalfSpaceYAEllipticalCylinder() {}

G4HalfSpaceYAEllipticalCylinder::G4HalfSpaceYAEllipticalCylinder(G4double x,
                                                                 G4double z,
                                                                 G4double r1,
                                                                 G4double r2) :
        _x0(x), _z0(z), _r1(r1), _r2(r2) {}

G4HalfSpaceYAEllipticalCylinder::~G4HalfSpaceYAEllipticalCylinder() {}

G4double G4HalfSpaceYAEllipticalCylinder::Sdf(const G4ThreeVector&p) const {

  auto x = p.x();
  auto y = p.y();
  auto z = p.z();

  auto sdf = pow(x-_x0,2)/pow(_r1,2) + pow(z-_z0,2)/pow(_r2,2) - 1;
  return sdf;
}

std::vector<G4ThreeVector> G4HalfSpaceYAEllipticalCylinder::Intersection(const G4ThreeVector& p, const G4ThreeVector &d) const {

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

  q(1,1) = 1./pow(_r1,2);
  q(3,3) = 1./pow(_r2,2);

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

  G4cout << "G4HalfSpaceYAEllipticalCylinder::Intersection n=" << inter.size() << std::endl;
  return inter;
}

void G4HalfSpaceYAEllipticalCylinder::Translate(const G4ThreeVector& t) {
  _x0 += t.x();
  _z0 += t.z();
}

void G4HalfSpaceYAEllipticalCylinder::Rotate(const G4RotationMatrix& r) {
  G4cout << "G4HalfSpaceYAEllipticalCylinder::Rotate> Not implemented" << G4endl;
}

void G4HalfSpaceYAEllipticalCylinder::Transform(const G4AffineTransform& a) {
  G4cout << "G4HalfSpaceYAEllipticalCylinder::Transform> Not implemented" << G4endl;
}

G4SurfaceMeshCGAL* G4HalfSpaceYAEllipticalCylinder::GetSurfaceMesh()  {
  G4EllipticalTube t = G4EllipticalTube("test", _r1, _r2, 1000000000);
  G4Polyhedron *g4poly = t.GetPolyhedron();
  G4SurfaceMeshCGAL *sm = new G4SurfaceMeshCGAL();
  sm->Fill(g4poly);
  sm->Translate(_x0,0,_z0);
  sm->Rotate(G4ThreeVector(1,0,0), M_PI_2);

  //return sm;
  return sm;
}