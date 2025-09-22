#include "G4HalfSpaceSolid.hh"
#include "G4SurfaceMeshCGAL.hh"
#include "G4VSceneHandler.hh"

G4HalfSpaceSolid::G4HalfSpaceSolid() : G4VSolid("") {};
G4HalfSpaceSolid::G4HalfSpaceSolid(const char *nameIn) : G4VSolid(nameIn) {};
G4HalfSpaceSolid::G4HalfSpaceSolid(const G4String &nameIn) : G4VSolid(nameIn) {};

G4HalfSpaceSolid::~G4HalfSpaceSolid() {};

G4double G4HalfSpaceSolid::Sdf(const G4ThreeVector &p) const {
  G4double sdf = kInfinity;
  for(auto z : _zones) {
    auto d = z->Sdf(p);
    sdf = std::min(d,sdf);
  }

  return sdf;
}

std::vector<G4ThreeVector> G4HalfSpaceSolid::Intersection(const G4ThreeVector& p, const G4ThreeVector &v) const {

  std::vector<G4ThreeVector> trialInter;
  for (auto z : _zones) {
    auto zoneIntersections = z->Intersection(p,v);
    for(auto zi : zoneIntersections) {
      trialInter.push_back(zi);
    }
  }

  // test intersections
  std::vector<G4ThreeVector> inter;
  for(auto i : trialInter) {
    auto sdf = Sdf(i);
    if(fabs(sdf)< 1e-9) {
      inter.push_back(i);
    }
  }

  return inter;
}

void G4HalfSpaceSolid::Translate(const G4ThreeVector& t) {
  for(auto z : _zones) {
    z->Translate(t);
  }
}

void G4HalfSpaceSolid::Rotate(const G4RotationMatrix& r) {
  for(auto z : _zones) {
    z->Rotate(r);
  }
}

void G4HalfSpaceSolid::Transform(const G4AffineTransform& a) {
  for(auto z : _zones) {
    z->Transform(a);
  }
}

G4SurfaceMeshCGAL* G4HalfSpaceSolid::GetSurfaceMesh() {

  G4cout << "G4HalfSpaceSolid::GetSurfaceMesh" << G4endl;

  G4SurfaceMeshCGAL *sm1 = new G4SurfaceMeshCGAL();

  for(auto z : _zones) {
    auto sm2 = z->GetSurfaceMesh();
    G4bool valid;
    sm1 = sm1->Union(sm2, valid);
  }

  return sm1;

}


void G4HalfSpaceSolid::AddZone(G4HalfSpaceZone *zone) {_zones.push_back(zone);}

void G4HalfSpaceSolid::RemoveZone(G4HalfSpaceZone *zone) {_zones.push_back(zone);}

G4int G4HalfSpaceSolid::NumberOfZones() { return _zones.size();}

std::vector<G4HalfSpaceSolid*> G4HalfSpaceSolid::connectedSolids() {
  return std::vector<G4HalfSpaceSolid*>();
}

G4bool G4HalfSpaceSolid::CalculateExtent(const EAxis pAxis,
                                         const G4VoxelLimits& pVoxelLimit,
                                         const G4AffineTransform& pTransform,
                                         G4double& pMin, G4double& pMax) const {
  return false;
}

EInside G4HalfSpaceSolid::Inside(const G4ThreeVector& p) const {

  G4double sdf = Sdf(p);

  G4HalfSpaceTestDataEntry *data = new G4HalfSpaceTestDataEntry();
  data->_type = G4HalfSpaceTestDataEntry::Type::Inside;
  data->_position = p;
  if(_test != nullptr) {
    _test->add_data(data);
  }

  if (sdf < -kCarTolerance/2.0) {
    data->_inside = kInside;
    return kInside;
  }
  else if (fabs(sdf) < kCarTolerance/2.0) {
    data->_inside = kSurface;
    return kSurface;
  }
  else {
    data->_inside = kOutside;
    return kOutside;
  }
}

G4ThreeVector G4HalfSpaceSolid::SurfaceNormal(const G4ThreeVector& p) const {
  G4HalfSpaceTestDataEntry *data = new G4HalfSpaceTestDataEntry();
  data->_type = G4HalfSpaceTestDataEntry::Type::Normal;
  data->_position = p;
  if(_test != nullptr) {
    _test->add_data(data);
  }
  data->_normal = Normal(p,1e-5);

  return Normal(p,1e-5);
}

G4double G4HalfSpaceSolid::DistanceToIn(const G4ThreeVector& p,
                                        const G4ThreeVector& v) const {
  G4HalfSpaceTestDataEntry *data = new G4HalfSpaceTestDataEntry();
  data->_type = G4HalfSpaceTestDataEntry::Type::DistanceToInDir;
  data->_position = p;
  data->_direction = v;
  if(_test != nullptr) {
    _test->add_data(data);
  }

  G4ThreeVector i2;

  if (Inside(p) == EInside::kInside) {
    data->_distance = 0;
    return 0;
  }

  std::vector<G4ThreeVector> trialInter;
  for (auto z: _zones) {
    auto zoneIntersections = z->Intersection(p, v);
    for (auto zi: zoneIntersections) {
      trialInter.push_back(zi);
    }
  }

  // test intersections
  std::vector<G4ThreeVector> inter;
  for (auto i: trialInter) {
    auto sdf = Sdf(i);
    if (fabs(sdf) < 1e-9) {
      inter.push_back(i);
    }
  }

  // order intersections in distance
  auto g4tvSort = ([p](G4ThreeVector v1, G4ThreeVector v2) {
    return ((v1 - p).mag() <= (v2 - p).mag());
  });
  std::sort(inter.begin(), inter.end(), g4tvSort);

  data->_ntrialintersections = trialInter.size();
  data->_nintersections = inter.size();

  for(auto i : inter) {
    if (SurfaceNormal(i).dot(v) <= 0) {
      data->_distance = (i - p).mag();
      return (i - p).mag();
    }
  }

  data->_distance = kInfinity;
  return kInfinity;
}

G4double G4HalfSpaceSolid::DistanceToIn(const G4ThreeVector& p) const {

  G4HalfSpaceTestDataEntry *data = new G4HalfSpaceTestDataEntry();
  data->_type = G4HalfSpaceTestDataEntry::Type::DistanceToIn;
  data->_position = p;
  if(_test != nullptr) {
    _test->add_data(data);
  }

  data->_distance = 0;
  return 0;

  G4double distToIn = Sdf(p);

  if(distToIn >= 0 ) {
    data->_distance = distToIn;
    return distToIn;
  }
  else {
    data->_distance = 0;
    return 0;
  }
}

G4double G4HalfSpaceSolid::DistanceToOut(const G4ThreeVector& p,
                                         const G4ThreeVector& v,
                                         const G4bool calcNorm,
                                         G4bool* validNorm,
                                         G4ThreeVector* n) const {

  G4HalfSpaceTestDataEntry *data = new G4HalfSpaceTestDataEntry();
  data->_type = G4HalfSpaceTestDataEntry::Type::DistanceToOutDir;
  data->_position = p;
  data->_direction = v;
  if(_test != nullptr) {
    _test->add_data(data);
  }

  if(Inside(p) == EInside::kOutside) {
    data->_distance = 0;
    return 0;
  }

  G4ThreeVector i2;

  std::vector<G4ThreeVector> trialInter;
  for (auto z : _zones) {
    auto zoneIntersections = z->Intersection(p,v);
    for(auto zi : zoneIntersections) {
      trialInter.push_back(zi);
    }
  }

  // test intersections
  std::vector<G4ThreeVector> inter;
  for(auto i : trialInter) {
    auto sdf = Sdf(i);
    if(fabs(sdf)< 1e-9) {
      inter.push_back(i);
    }
  }

  // order intersections in distance
  auto g4tvSort = ([p](G4ThreeVector v1, G4ThreeVector v2) {
    return ((v1-p).mag()<=(v2-p).mag());
  });
  std::sort(inter.begin(), inter.end(), g4tvSort);

  data->_ntrialintersections = trialInter.size();
  data->_nintersections = inter.size();

  for(auto i : inter) {
    if (SurfaceNormal(i).dot(v) >= 0) {
      data->_distance = (i-p).mag();
      return (i-p).mag();
    }
  }

  data->_distance = 0;
  return 0;
}

G4double G4HalfSpaceSolid::DistanceToOut(const G4ThreeVector& p) const {

  G4HalfSpaceTestDataEntry *data = new G4HalfSpaceTestDataEntry();
  data->_type = G4HalfSpaceTestDataEntry::Type::DistanceToOut;
  data->_position = p;
  if(_test != nullptr) {
    _test->add_data(data);
  }

  data->_distance = 0;
  return 0;

  G4double distToOut= Sdf(p);

  if(distToOut < 0 ) {
    data->_distance = fabs(distToOut);
    return fabs(distToOut);
  }
  else {
    data->_distance = 0;
    return 0;
  }
}

G4GeometryType G4HalfSpaceSolid::GetEntityType() const {
  return G4String("halfSpaceSolid");
}

std::ostream& G4HalfSpaceSolid::StreamInfo(std::ostream& os) const {
  return os;
}

void G4HalfSpaceSolid::DescribeYourselfTo(G4VGraphicsScene& scene) const {
  G4cout << "G4HalfSpaceSolid::DescribeYourselfTo" << std::endl;

  auto sm = const_cast<G4HalfSpaceSolid*>(this)->GetSurfaceMesh();
  auto ph = sm->GetG4Polyhedron();
  scene.AddPrimitive(*((G4Polyhedron*)ph));

  return;
}

void G4HalfSpaceSolid::AddTestInstrument(G4HalfSpaceTest *test) {
  _test = test;
}