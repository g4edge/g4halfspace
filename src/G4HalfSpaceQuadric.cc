#include "G4HalfSpaceQuadric.hh"
#include "G4HalfSpaceRotation.hh"
#include "G4ImplicitSurfaceCGAL.hh"

#include "G4Ellipsoid.hh"

#include <Eigen/Dense>

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

G4HalfSpaceQuadric::G4HalfSpaceQuadric(const G4HalfSpaceQuadric &other) :
  _q(other._q), _p(other._p), _r(other._r) {}

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

#if _DEBUG
  G4cout << "G4HalfSpaceQuadric::Rotate(const G4ThreeVector &rv)" << G4endl;
  G4cout << rm << G4endl;
#endif
  this->Rotate(rm);
}

void G4HalfSpaceQuadric::Rotate(const G4RotationMatrix& rm) {
#if _DEBUG
  G4cout << "G4HalfSpaceQuadric::Rotate(const G4RotationMatrix& r)" << G4endl;
  G4cout << rm << G4endl;
#endif
  auto rminv = rm.inverse();
  auto m = CLHEP::HepMatrix(3,3);
  m(1,1) = rminv.xx();
  m(1,2) = rminv.xy();
  m(1,3) = rminv.xz();
  m(2,1) = rminv.yx();
  m(2,2) = rminv.yy();
  m(2,3) = rminv.yz();
  m(3,1) = rminv.zx();
  m(3,2) = rminv.zy();
  m(3,3) = rminv.zz();

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

#if 0
  // determine rotation of quadric
  G4cout << ToStringEquation() << G4endl;
  G4ThreeVector evals;
  auto evecs = GetRotationFromQuadraticForm(evals);
  G4HalfSpaceRotation rot = G4HalfSpaceRotation(evecs);
#if _DEBUG
  G4cout << rot.ToString() << G4endl;
  G4cout << "G4HalfSpaceQuadric::GetSurfaceMesh" << std::endl;
#endif
  G4Ellipsoid t = G4Ellipsoid("test", 10, 20, 30);
  G4Polyhedron *g4poly = t.GetPolyhedron();
  G4SurfaceMeshCGAL *sm = new G4SurfaceMeshCGAL();
  sm->Fill(g4poly);
  return sm;
#endif

  // mesh quadric mesh
  return make_mesh(*this, 3000*3000);


}

CLHEP::HepMatrix G4HalfSpaceQuadric::GetRotationFromQuadraticForm(G4ThreeVector &g4_eigenvalues) {
  Eigen::Matrix3d mtemp;
  mtemp << _q(1,1), _q(1,2), _q(1,3),
           _q(2,1), _q(2,2), _q(2,3),
           _q(3,1), _q(3,2), _q(3,3);

  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> stemp(mtemp);

  Eigen::Vector3d eigenvalues = stemp.eigenvalues().real();
  Eigen::Matrix3d eigenvectors = stemp.eigenvectors().real();

  //if (eigenvectors.determinant() < 0) eigenvectors.col(0) *= -1;

  if (stemp.info() != Eigen::Success) {
    std::cerr << "Error: Eigenvalue computation failed!" << std::endl;
  }

  Eigen::MatrixXd eigenvectors_reordered(eigenvectors.rows(), eigenvectors.cols());
  Eigen::Vector3d eigenvalues_reordered(eigenvalues.size());
  for (int j = 0; j< eigenvectors.cols(); j++) {
      eigenvectors_reordered.col(2-j) = eigenvectors.col(j);
      eigenvalues_reordered(2-j) = eigenvalues[j];
  }

  G4cout << "Eigenvalues:\n" << eigenvalues_reordered << G4endl;
  G4cout << "Eigenvectors:\n" << eigenvectors_reordered << G4endl;
  G4cout << "Determinant:" << eigenvectors_reordered.determinant() << G4endl;
  g4_eigenvalues.set(eigenvalues[2],eigenvalues[1],eigenvalues[0]);

  CLHEP::HepMatrix g4_eigenvectors = CLHEP::HepMatrix(3,3);

  g4_eigenvectors(1,1) = eigenvectors_reordered(0,0);
  g4_eigenvectors(1,2) = eigenvectors_reordered(0,1);
  g4_eigenvectors(1,3) = eigenvectors_reordered(0,2);

  g4_eigenvectors(2,1) = eigenvectors_reordered(1,0);
  g4_eigenvectors(2,2) = eigenvectors_reordered(1,1);
  g4_eigenvectors(2,3) = eigenvectors_reordered(1,2);

  g4_eigenvectors(3,1) = eigenvectors_reordered(2,0);
  g4_eigenvectors(3,2) = eigenvectors_reordered(2,1);
  g4_eigenvectors(3,3) = eigenvectors_reordered(2,2);

  return g4_eigenvectors;

}

G4ThreeVector G4HalfSpaceQuadric::GetTranslationFromQuadricEqn() {
  return G4ThreeVector();
}

std::string G4HalfSpaceQuadric::ToStringEquation() {
  std::stringstream ss;

  ss << ValueToStringPretty(_q(1,1), "x^2")
     << ValueToStringPretty(_q(2,2), "y^2")
     << ValueToStringPretty(_q(3,3), "z^2")
     << ValueToStringPretty(_q(1,2) + _q(2,1), "xy")
     << ValueToStringPretty(_q(1,3) + _q(3,1), "xz")
     << ValueToStringPretty(_q(2,3) + _q(3,2), "yz")
     << ValueToStringPretty(_p(1), "x")
     << ValueToStringPretty(_p(2), "y")
     << ValueToStringPretty(_p(3), "z")
     << ValueToStringPretty(_r, "");

  return ss.str();
}

std::string G4HalfSpaceQuadric::ValueToStringPretty(double val, std::string term) {
  std::stringstream ss;

  if(val > 0) {
    ss << "+" << val;
    if(term != "")
      ss << "*" << term;
  }
  else if(val < 0) {
    ss << val;
    if(term != "")
      ss << "*" << term;
  }

  return ss.str();
}