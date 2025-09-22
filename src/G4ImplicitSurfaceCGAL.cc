#include "G4ImplicitSurfaceCGAL.hh"
#include "G4HalfSpaceQuadric.hh"
#include "G4SurfaceMeshCGAL.hh"

FT_GT sphere_function (Point_3_GT p) {
  const FT_GT x2=p.x()*p.x(), y2=p.y()*p.y(), z2=p.z()*p.z();
  return x2+y2+z2-10*10;
}

G4SurfaceMeshCGAL* make_mesh(G4HalfSpaceQuadric &quadric, double sphere_size) {
  std::cout << "G4ImplicitSufaceCGAL make_mesh" << std::endl;
  Tr tr;            // 3D-Delaunay triangulation
  C2t3 c2t3 (tr);   // 2D-complex in 3D-Delaunay triangulation

  std::function<FT_GT(Point_3_GT)> surface_function = [&quadric](Point_3_GT p){
    G4ThreeVector g4_p(p.x(), p.y(), p.z());
    return quadric.Sdf(g4_p);
  };

  // defining the surface
  Surface_3_GT surface(surface_function,             // pointer to function
                       Sphere_3_GT(CGAL::ORIGIN, sphere_size)); // bounding sphere

  // defining meshing criteria
  CGAL::Surface_mesh_default_criteria_3<Tr> criteria(30.,  // angular bound
                                                     10,  // radius bound
                                                     10); // distance bound
  // meshing surface
  CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Non_manifold_tag());

  Surface_mesh_3 *sm = new Surface_mesh_3();
  CGAL::facets_in_complex_2_to_triangle_mesh(c2t3, *sm);
  std::cout << "G4ImplicitSufaceCGAL Final number of points: " << tr.number_of_vertices() << "\n";
  CGAL::Polygon_mesh_processing::reverse_face_orientations(*sm);
  G4SurfaceMeshCGAL *g4_sm = new G4SurfaceMeshCGAL(sm);
  delete sm;

  return g4_sm;
}