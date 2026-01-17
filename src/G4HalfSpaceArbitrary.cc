#include "G4HalfSpaceArbitrary.hh"
#include "G4HalfSpaceZone.hh"
#include "G4HalfSpacePlane.hh"
#include "G4SurfaceMeshCGAL.hh"
#include "G4Box.hh"

#include <vector>

G4HalfSpaceArbitrary::G4HalfSpaceArbitrary() {}

G4HalfSpaceArbitrary::G4HalfSpaceArbitrary(const G4ThreeVector &v1,
                                           const G4ThreeVector &v2,
                                           const G4ThreeVector &v3,
                                           const G4ThreeVector &v4,
                                           const G4ThreeVector &v5,
                                           const G4ThreeVector &v6,
                                           const G4ThreeVector &v7,
                                           const G4ThreeVector &v8,
                                           G4int fi1,
                                           G4int fi2,
                                           G4int fi3,
                                           G4int fi4,
                                           G4int fi5,
                                           G4int fi6) {

  _verts.push_back(v1);
  _verts.push_back(v2);
  _verts.push_back(v3);
  _verts.push_back(v4);
  _verts.push_back(v5);
  _verts.push_back(v6);
  _verts.push_back(v7);
  _verts.push_back(v8);

  auto fi1d = decodeFaceInteger(fi1);
  auto fi2d = decodeFaceInteger(fi2);
  auto fi3d = decodeFaceInteger(fi3);
  auto fi4d = decodeFaceInteger(fi4);
  auto fi5d = decodeFaceInteger(fi5);
  auto fi6d = decodeFaceInteger(fi6);

  _faceVertCount.push_back(fi1d.size());
  _faceVertCount.push_back(fi2d.size());
  _faceVertCount.push_back(fi3d.size());
  _faceVertCount.push_back(fi4d.size());
  _faceVertCount.push_back(fi5d.size());
  _faceVertCount.push_back(fi6d.size());

  _faceVertIndex.insert(_faceVertIndex.end(), fi1d.begin(), fi1d.end());
  _faceVertIndex.insert(_faceVertIndex.end(), fi2d.begin(), fi2d.end());
  _faceVertIndex.insert(_faceVertIndex.end(), fi3d.begin(), fi3d.end());
  _faceVertIndex.insert(_faceVertIndex.end(), fi4d.begin(), fi4d.end());
  _faceVertIndex.insert(_faceVertIndex.end(), fi5d.begin(), fi5d.end());
  _faceVertIndex.insert(_faceVertIndex.end(), fi6d.begin(), fi6d.end());

  G4cout << ToString() << std::endl;

  this->ComputePlanes();
}

G4HalfSpaceArbitrary::G4HalfSpaceArbitrary(const G4ThreeVector &c,
                                           const G4ThreeVector &v1,
                                           const G4ThreeVector &v2,
                                           const G4ThreeVector &v3,
                                           const G4ThreeVector &v4,
                                           const G4ThreeVector &v5,
                                           const G4ThreeVector &v6,
                                           const G4ThreeVector &v7,
                                           const G4ThreeVector &v8,
                                           G4int fi1,
                                           G4int fi2,
                                           G4int fi3,
                                           G4int fi4,
                                           G4int fi5,
                                           G4int fi6) :
G4HalfSpaceArbitrary(v1+c, v2+c, v3+c, v4+c, v5+c, v6+c, v7+c, v8+c,
                     fi1, fi2, fi3, fi4, fi5, fi6) {

}

void G4HalfSpaceArbitrary::ComputePlanes() {

  std::vector<G4int>::iterator fcIter = _faceVertIndex.begin();
  for(auto fc : _faceVertCount) {
    std::vector<G4int> planeVertIndices;

    for(G4int fci = 0; fci < fc; fci++) {
      planeVertIndices.push_back(*fcIter);
      fcIter++;
    }

    if(planeVertIndices.size() < 3)
      continue;

    auto pv1 = _verts[planeVertIndices[0]-1];
    auto pv2 = _verts[planeVertIndices[1]-1];
    auto pv3 = _verts[planeVertIndices[2]-1];

    auto d1 = pv1-pv2;
    auto d2 = pv3-pv2;

    auto n = d2.cross(d1);

    G4cout << "(" << planeVertIndices[0] << "," <<  planeVertIndices[1] << "," << planeVertIndices[2] << ") " << fc << " " << pv2 << " " << n << std::endl;

    n = n/n.mag();


    auto pl = new G4HalfSpacePlane(n, pv2);
    _hsZone.AddIntersection(pl);
  }
}

G4HalfSpaceArbitrary::~G4HalfSpaceArbitrary() {
}

G4double G4HalfSpaceArbitrary::Sdf(const G4ThreeVector&p) const {
  return _hsZone.Sdf(p);
}

std::vector<G4ThreeVector> G4HalfSpaceArbitrary::Intersection(const G4ThreeVector& p, const G4ThreeVector& v) const {
  auto intersections = _hsZone.Intersection(p,v);
  return intersections;
}

void G4HalfSpaceArbitrary::Translate(const G4ThreeVector& t) {
  _hsZone.Translate(t);
}

void G4HalfSpaceArbitrary::Rotate(const G4RotationMatrix& r) {
  _hsZone.Rotate(r);
}

void G4HalfSpaceArbitrary::Transform(const G4AffineTransform& a) {
  _hsZone.Transform(a);
}

G4SurfaceMeshCGAL* G4HalfSpaceArbitrary::GetSurfaceMesh()  {
  if(cached_mesh) {
    G4cout << "G4HalfSpaceArbitrary::GetSurfaceMesh cached" << G4endl;
    return cached_mesh;
  }

  cached_mesh = _hsZone.GetSurfaceMesh();

  return cached_mesh;
}

std::vector<G4int> G4HalfSpaceArbitrary::decodeFaceInteger(G4int flukaFaceDef) {
  auto vertList = std::vector<G4int>();
  auto flukaFaceDefString = std::to_string(flukaFaceDef);

  for(auto c : flukaFaceDefString) {
    if(c != '0') {
      vertList.push_back(std::stoi(std::string(1,c)));
    }
  }

  return vertList;
}

std::string G4HalfSpaceArbitrary::ToString() {
  std::stringstream sstr;

  for(auto v : _verts) {
    sstr << v << std::endl;
  }

  for(auto fc : _faceVertCount) {
    sstr << fc << " ";
  }
  sstr << std::endl;

  for(auto fvi : _faceVertIndex) {
    sstr << fvi << " ";
  }
  sstr << std::endl;

  return sstr.str();
}