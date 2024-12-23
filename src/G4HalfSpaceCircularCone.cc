#include "G4HalfSpaceCircularCone.hh"
#include "G4HalfSpaceZone.hh"
#include "G4HalfSpacePlane.hh"
#include "G4HalfSpaceQuadric.hh"
#include "G4SurfaceMeshCGAL.hh"

G4HalfSpaceCircularCone::G4HalfSpaceCircularCone() {
}

G4HalfSpaceCircularCone::G4HalfSpaceCircularCone(const G4ThreeVector &v,
                                                 const G4ThreeVector &h,
                                                 G4double r1, G4double r2) {
  auto hmag = h.mag();

}


G4HalfSpaceCircularCone::G4HalfSpaceCircularCone(G4double h, G4double r1, G4double r2,
                                                 const G4ThreeVector &v, const G4ThreeVector &r) :
  _h(h), _r1(r1), _r2(r2), _v(v), _r(h)
{
  double g = 0;

  if(r1>r2) {
    auto rtemp1 = _r1;
    auto rtemp2 = _r2;

    r1 = rtemp2;
    r2 = rtemp1;
  }

  g = (r2 - r1)/h;

  auto z1 = r1/g;
  auto z2 = r2/g;

  auto c = 1;
  auto a = (r2-r1)/h;

  G4ThreeVector shift = G4ThreeVector(0,0,-(z2+z1)/2.0);

  G4HalfSpaceQuadric *q = new G4HalfSpaceQuadric(1/pow(a,2), 0, 0,
                                                 1/pow(a,2), 0,
                                                 -1/pow(c,2),
                                                 0,0,0,
                                                 0);

  G4HalfSpacePlane *p1 = new G4HalfSpacePlane(G4ThreeVector(0,0,z1),
                                              G4ThreeVector(0,0,-1));
  G4HalfSpacePlane *p2 = new G4HalfSpacePlane(G4ThreeVector(0,0,z2),
                                              G4ThreeVector(0,0,1));


  _hsZone.AddIntersection(q);
  _hsZone.AddIntersection(p1);
  _hsZone.AddIntersection(p2);

  _hsZone.Translate(shift);

  if(_r1 > _r2) {
    G4RotationMatrix rot;
    rot.rotateX(M_PI);
    rot.rectify();
    _hsZone.Rotate(rot);
  }

}


G4HalfSpaceCircularCone::~G4HalfSpaceCircularCone() {
}

G4double G4HalfSpaceCircularCone::Sdf(const G4ThreeVector&p) const {
  return _hsZone.Sdf(p);
}

std::vector<G4ThreeVector> G4HalfSpaceCircularCone::Intersection(const G4ThreeVector& p, const G4ThreeVector& v) const {
  auto intersections = _hsZone.Intersection(p,v);
  return intersections;
}

void G4HalfSpaceCircularCone::Translate(const G4ThreeVector& t) {
  _hsZone.Translate(t);
}

void G4HalfSpaceCircularCone::Rotate(const G4RotationMatrix& r) {
  _hsZone.Rotate(r);
}

void G4HalfSpaceCircularCone::Transform(const G4AffineTransform& a) {
  _hsZone.Transform(a);
}

G4SurfaceMeshCGAL* G4HalfSpaceCircularCone::GetSurfaceMesh()  {
  if(cached_mesh) {
    G4cout << "G4HalfSpaceCircularCone::GetSurfaceMesh cached" << G4endl;
    return cached_mesh;
  }


  cached_mesh = _hsZone.GetSurfaceMesh();

  return cached_mesh;
}