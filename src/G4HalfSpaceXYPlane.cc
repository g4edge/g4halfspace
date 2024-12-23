#include "G4HalfSpaceXYPlane.hh"
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

G4HalfSpaceXYPlane::G4HalfSpaceXYPlane() : _z(0)
{}

G4HalfSpaceXYPlane::G4HalfSpaceXYPlane(G4double z) : _z(z)
{}

G4HalfSpaceXYPlane::~G4HalfSpaceXYPlane() {};

G4double G4HalfSpaceXYPlane::Sdf(const G4ThreeVector &p) const {
  return p.z() - _z;
}

std::vector<G4ThreeVector> G4HalfSpaceXYPlane::Intersection(const G4ThreeVector& p, const G4ThreeVector &v) const {
  std::vector<G4ThreeVector> intersections;

  auto vNorm = v/v.mag();
  auto lambda = (_z-p.z())/ vNorm.z();

  if(vNorm.z() != 0 && lambda >= 0)
    intersections.push_back(lambda*v+p);

  return intersections;
}

void G4HalfSpaceXYPlane::Translate(const G4ThreeVector& t) {
}

void G4HalfSpaceXYPlane::Rotate(const G4RotationMatrix& r) {
}

void G4HalfSpaceXYPlane::Transform(const G4AffineTransform& a) {
}


G4SurfaceMeshCGAL* G4HalfSpaceXYPlane::GetSurfaceMesh() {
  if(cached_mesh) {
    G4cout << "G4HalfSpaceXYPlane::GetSurfaceMesh cached" << G4endl;
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

  nef *= Nef_polyhedron_3_ECER(Plane_3_ECER(Point_3_ECER(0, 0, _z),
                                            Direction_3_ECER(0, 0, 1)));

  cached_mesh = new G4SurfaceMeshCGAL(nef);
  return cached_mesh;
}