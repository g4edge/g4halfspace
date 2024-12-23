#include "G4HalfSpaceXZPlane.hh"
#include "G4VSolid.hh"

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_rational.h>
#include <CGAL/Extended_cartesian.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Surface_mesh.h>

typedef CGAL::Exact_rational ER;
typedef CGAL::Extended_cartesian<ER> Kernel_ECER;
typedef CGAL::Nef_polyhedron_3<Kernel_ECER> Nef_polyhedron_3_ECER;
typedef CGAL::Surface_mesh<Kernel_ECER> Surface_mesh_ECER;
typedef Kernel_ECER::Point_3 Point_3_ECER;
typedef Kernel_ECER::Vector_3 Vector_3_ECER;
typedef Kernel_ECER::Plane_3 Plane_3_ECER;
typedef Kernel_ECER::Direction_3 Direction_3_ECER;

G4HalfSpaceXZPlane::G4HalfSpaceXZPlane() : _y(0)
{}

G4HalfSpaceXZPlane::G4HalfSpaceXZPlane(G4double y) : _y(y)
{}

G4HalfSpaceXZPlane::~G4HalfSpaceXZPlane() {};

G4double G4HalfSpaceXZPlane::Sdf(const G4ThreeVector &p) const {
  return p.y() - _y;
}

std::vector<G4ThreeVector> G4HalfSpaceXZPlane::Intersection(const G4ThreeVector& p, const G4ThreeVector &v) const {
  std::vector<G4ThreeVector> intersections;

  auto vNorm = v/v.mag();
  auto lambda = (_y-p.y())/ vNorm.y();

  if(vNorm.y() != 0 && lambda >= 0)
    intersections.push_back(lambda*v+p);

  return intersections;
}

void G4HalfSpaceXZPlane::Translate(const G4ThreeVector& t) {
}

void G4HalfSpaceXZPlane::Rotate(const G4RotationMatrix& r) {
}

void G4HalfSpaceXZPlane::Transform(const G4AffineTransform& a) {
}


G4SurfaceMeshCGAL* G4HalfSpaceXZPlane::GetSurfaceMesh() {
  if(cached_mesh) {
    G4cout << "G4HalfSpaceXZPlane::GetSurfaceMesh cached" << G4endl;
    return cached_mesh;
  }

  Nef_polyhedron_3_ECER nef = Nef_polyhedron_3_ECER(Nef_polyhedron_3_ECER::COMPLETE);
  nef *= Nef_polyhedron_3_ECER(Plane_3_ECER(Point_3_ECER(1000000000,0,0),
                                            Direction_3_ECER(1, 0, 0)));
  nef *= Nef_polyhedron_3_ECER(Plane_3_ECER(Point_3_ECER(-1000000000,0,0),
                                            Direction_3_ECER(-1, 0, 0)));
  nef *= Nef_polyhedron_3_ECER(Plane_3_ECER(Point_3_ECER(0,1000000000,0),
                                            Direction_3_ECER(0, 1, 0)));
  nef *= Nef_polyhedron_3_ECER(Plane_3_ECER(Point_3_ECER(0,-1000000000,0),
                                            Direction_3_ECER(0, -1, 0)));
  nef *= Nef_polyhedron_3_ECER(Plane_3_ECER(Point_3_ECER(0,0,1000000000),
                                            Direction_3_ECER(0, 0, 1)));
  nef *= Nef_polyhedron_3_ECER(Plane_3_ECER(Point_3_ECER(0,0,-1000000000),
                                            Direction_3_ECER(0, 0, -1)));

  nef *= Nef_polyhedron_3_ECER(Plane_3_ECER(Point_3_ECER(0, _y, 0),
                                            Direction_3_ECER(0, 1, 0)));

  cached_mesh = new G4SurfaceMeshCGAL(nef);
  return cached_mesh;
}