#include "G4HalfSpaceYZPlane.hh"
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

G4HalfSpaceYZPlane::G4HalfSpaceYZPlane() : _x(0)
{}

G4HalfSpaceYZPlane::G4HalfSpaceYZPlane(G4double x) : _x(x)
{}

G4HalfSpaceYZPlane::~G4HalfSpaceYZPlane() {};

G4double G4HalfSpaceYZPlane::Sdf(const G4ThreeVector &p) const {
  return p.x() - _x;
}

std::vector<G4ThreeVector> G4HalfSpaceYZPlane::Intersection(const G4ThreeVector& p, const G4ThreeVector &v) const {
  std::vector<G4ThreeVector> intersections;

  auto vNorm = v/v.mag();
  auto lambda = (_x-p.x())/ vNorm.x();

  if(vNorm.x() != 0 && lambda >= 0)
    intersections.push_back(lambda*v+p);

  return intersections;
}

void G4HalfSpaceYZPlane::Translate(const G4ThreeVector& t) {
}

void G4HalfSpaceYZPlane::Rotate(const G4RotationMatrix& r) {
}

void G4HalfSpaceYZPlane::Transform(const G4AffineTransform& a) {
}


G4SurfaceMeshCGAL* G4HalfSpaceYZPlane::GetSurfaceMesh() {
  if(cached_mesh) {
    G4cout << "G4HalfSpaceYZPlane::GetSurfaceMesh cached" << G4endl;
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

  nef *= Nef_polyhedron_3_ECER(Plane_3_ECER(Point_3_ECER(_x, 0, 0),
                                            Direction_3_ECER(1, 0, 0)));

  cached_mesh = new G4SurfaceMeshCGAL(nef);
  return cached_mesh;
}