#include "G4HalfSpaceQuadric.hh"

#include "G4Ellipsoid.hh"

G4HalfSpaceQuadric::G4HalfSpaceQuadric() {}

G4HalfSpaceQuadric::G4HalfSpaceQuadric(double qxx, double qxy, double qxz,
                                                   double qyy, double qyz,
                                                               double qzz,
                                       double px,  double py,  double pz,
                                       double r) {
  _q(1,1) = qxx;
  _q(1,2) = qxy;
  _q(1,3) = qxz;
  _q(2,2) = qyy;
  _q(2,3) = qyz;
  _q(3,3) = qzz;

  _p(1) = px;
  _p(2) = py;
  _p(3) = pz;

  _r = r;

}

G4HalfSpaceQuadric::G4HalfSpaceQuadric(CLHEP::HepMatrix &q,
                                       CLHEP::HepVector &p,
                                       double r) : _q(q), _p(p), _r(r) {}


G4double G4HalfSpaceQuadric::Sdf(const G4ThreeVector&p) const {

  CLHEP::HepVector pv = CLHEP::HepVector(3);
  pv(1) = p.x();
  pv(2) = p.y();
  pv(3) = p.z();
  CLHEP::HepVector rv = CLHEP::HepVector(1);
  rv(1) = _r;

  auto ls = pv.T() * _q * pv + _p.T() * pv + rv;

  return ls(1);
}

std::vector<G4ThreeVector> G4HalfSpaceQuadric::Intersection(const G4ThreeVector& p, const G4ThreeVector &d) const {

  auto n = CLHEP::HepVector(3);
  auto p0 = CLHEP::HepVector(3);

  n(1) = d.x();
  n(2) = d.y();
  n(3) = d.z();

  p0(1) = p.x();
  p0(2) = p.y();
  p0(3) = p.z();

  auto a = (n.T()*_q*n)(1,1);
  auto b = (n.T()*_q*p0 + _p.T()*n + p0.T()*_q*n)(1,1);
  auto c = (_p.T()*p0 + p0.T()*_q*p0)(1,1) + _r;

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

void G4HalfSpaceQuadric::Translate(const G4ThreeVector& t3v) {
  auto t = CLHEP::HepVector(3);
  t(1) = -t3v[0];
  t(2) = -t3v[1];
  t(3) = -t3v[2];

  auto pp = _q*t + _q.T()*t + _p;
  auto rp = _r + (t.T()*_q*t + _p.T()*t)(1);

  _p = pp;
  _r = rp;


}

void G4HalfSpaceQuadric::Rotate(const G4ThreeVector &rv) {
  G4RotationMatrix rm = G4RotationMatrix();
  rm.rotateZ(-rv[2]);
  rm.rotateY(-rv[1]);
  rm.rotateX(-rv[0]);
  rm.rectify();
  this->Rotate(rm);
}

void G4HalfSpaceQuadric::Rotate(const G4RotationMatrix& r) {
  auto rinv = r.inverse();
  auto m = CLHEP::HepMatrix(3,3);
  m(1,1) = rinv.xx();
  m(1,2) = rinv.xy();
  m(1,3) = rinv.xz();
  m(2,1) = rinv.yx();
  m(2,2) = rinv.yy();
  m(2,3) = rinv.yz();
  m(3,1) = rinv.zx();
  m(3,2) = rinv.zy();
  m(3,3) = rinv.zz();

  auto qp = m.T() * _q * m;
  auto pp = _p.T() * m;
  auto rp = _r;

  _q = qp;
  _p = pp.T();
  _r = rp;
}

void G4HalfSpaceQuadric::Transform(const G4AffineTransform& a) {
  this->Rotate(a.NetRotation());
  this->Translate(a.NetTranslation());
}

G4SurfaceMeshCGAL* G4HalfSpaceQuadric::GetSurfaceMesh() {
  std::cout << "G4HalfSpaceQuadric::GetSurfaceMesh" << std::endl;
  G4Ellipsoid t = G4Ellipsoid("test", 10, 20, 30);
  G4Polyhedron *g4poly = t.GetPolyhedron();
  G4SurfaceMeshCGAL *sm = new G4SurfaceMeshCGAL();
  sm->Fill(g4poly);

  return sm;
}
