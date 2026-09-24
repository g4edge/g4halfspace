#include "G4LoaderHalfSpace.hh"

#include <algorithm>
#include <cmath>

#include "G4HalfSpacePlane.hh"
#include "G4HalfSpaceQuadric.hh"
#include "G4HalfSpaceSphere.hh"
#include "G4RotationMatrix.hh"

#if defined(G4HALFSPACE_USE_OPENCASCADE)
  #include <BRepAdaptor_Surface.hxx>
  #include <GeomAbs_SurfaceType.hxx>
  #include <IFSelect_ReturnStatus.hxx>
  #include <STEPControl_Reader.hxx>
  #include <TopAbs_Orientation.hxx>
  #include <TopAbs_ShapeEnum.hxx>
  #include <TopExp_Explorer.hxx>
  #include <TopoDS.hxx>
  #include <TopoDS_Face.hxx>
  #include <TopoDS_Shape.hxx>
#endif

namespace {

G4ThreeVector UnitOrDefault(const G4ThreeVector &v, const G4ThreeVector &fallback = G4ThreeVector(0, 0, 1)) {
  if (v.mag2() <= 0.0) {
    return fallback;
  }

  return v.unit();
}

G4RotationMatrix RotationFromZAxis(const G4ThreeVector &direction) {
  auto to = UnitOrDefault(direction);
  const auto from = G4ThreeVector(0, 0, 1);

  auto dot = from.dot(to);
  dot = std::max(-1.0, std::min(1.0, dot));

  G4RotationMatrix rotation;
  if (dot > 1.0 - 1e-12) {
    return rotation;
  }

  if (dot < -1.0 + 1e-12) {
    rotation.rotateX(M_PI);
    rotation.rectify();
    return rotation;
  }

  auto axis = from.cross(to).unit();
  auto angle = std::acos(dot);
  rotation.set(axis, angle);
  rotation.rectify();
  return rotation;
}

} // namespace

const std::map<std::size_t, G4LoaderHalfSpace::SurfaceRecord>& G4LoaderHalfSpace::GetSurfaces() const {
  return surfaces_;
}

void G4LoaderHalfSpace::Clear() {
  surfaces_.clear();
}

bool G4LoaderHalfSpace::LoadStep(const G4String &file_name) {
  surfaces_.clear();

#if defined(G4HALFSPACE_USE_OPENCASCADE)
  STEPControl_Reader step_reader;
  if (step_reader.ReadFile(file_name.c_str()) != IFSelect_RetDone) {
    return false;
  }

  if (step_reader.TransferRoots() == 0) {
    return false;
  }

  TopoDS_Shape shape = step_reader.OneShape();
  std::size_t surface_id = 0;

  for (TopExp_Explorer explorer(shape, TopAbs_FACE); explorer.More(); explorer.Next()) {
    auto face = TopoDS::Face(explorer.Current());
    BRepAdaptor_Surface surface(face, Standard_True);

    SurfaceRecord record;
    record.reversed = face.Orientation() == TopAbs_REVERSED;

    switch (surface.GetType()) {
      case GeomAbs_Plane: {
        auto plane = surface.Plane();
        record.type = SurfaceType::Plane;
        record.location = G4ThreeVector(plane.Location().X(), plane.Location().Y(), plane.Location().Z());
        record.direction = UnitOrDefault(G4ThreeVector(plane.Axis().Direction().X(),
                                                       plane.Axis().Direction().Y(),
                                                       plane.Axis().Direction().Z()));
        break;
      }
      case GeomAbs_Cylinder: {
        auto cylinder = surface.Cylinder();
        record.type = SurfaceType::Cylinder;
        record.location = G4ThreeVector(cylinder.Location().X(), cylinder.Location().Y(), cylinder.Location().Z());
        record.direction = UnitOrDefault(G4ThreeVector(cylinder.Axis().Direction().X(),
                                                       cylinder.Axis().Direction().Y(),
                                                       cylinder.Axis().Direction().Z()));
        record.radius = cylinder.Radius();
        break;
      }
      case GeomAbs_Cone: {
        auto cone = surface.Cone();
        record.type = SurfaceType::Cone;
        record.location = G4ThreeVector(cone.Location().X(), cone.Location().Y(), cone.Location().Z());
        record.direction = UnitOrDefault(G4ThreeVector(cone.Axis().Direction().X(),
                                                       cone.Axis().Direction().Y(),
                                                       cone.Axis().Direction().Z()));
        record.refRadius = cone.RefRadius();
        record.semiAngle = cone.SemiAngle();
        break;
      }
      case GeomAbs_Sphere: {
        auto sphere = surface.Sphere();
        record.type = SurfaceType::Sphere;
        record.location = G4ThreeVector(sphere.Location().X(), sphere.Location().Y(), sphere.Location().Z());
        record.radius = sphere.Radius();
        break;
      }
      default:
        record.type = SurfaceType::Unsupported;
        break;
    }

    if (record.type != SurfaceType::Unsupported) {
      surfaces_[surface_id++] = record;
    }
  }

  return !surfaces_.empty();
#else
  (void)file_name;
  return false;
#endif
}

std::map<std::size_t, G4VHalfSpace*> G4LoaderHalfSpace::CreateHalfSpaces() const {
  std::map<std::size_t, G4VHalfSpace*> half_spaces;

  for (const auto &[surface_id, record] : surfaces_) {
    auto *half_space = CreateHalfSpace(record);
    if (half_space != nullptr) {
      half_spaces[surface_id] = half_space;
    }
  }

  return half_spaces;
}

G4VHalfSpace* G4LoaderHalfSpace::CreateHalfSpace(const SurfaceRecord &record) const {
  switch (record.type) {
    case SurfaceType::Plane: {
      auto normal = record.reversed ? -record.direction : record.direction;
      return new G4HalfSpacePlane(UnitOrDefault(normal), record.location);
    }

    case SurfaceType::Sphere: {
      if (record.radius <= 0.0) {
        return nullptr;
      }
      return new G4HalfSpaceSphere(record.location, record.radius);
    }

    case SurfaceType::Cylinder: {
      if (record.radius <= 0.0) {
        return nullptr;
      }

      auto *quadric = new G4HalfSpaceQuadric(1.0 / std::pow(record.radius, 2), 0, 0,
                                             1.0 / std::pow(record.radius, 2), 0,
                                             0,
                                             0, 0, 0,
                                             -1);
      quadric->Rotate(RotationFromZAxis(record.direction));
      quadric->Translate(record.location);
      return quadric;
    }

    case SurfaceType::Cone: {
      auto tan_angle = std::tan(record.semiAngle);
      if (std::abs(tan_angle) <= 1e-12) {
        return nullptr;
      }

      auto apex = record.location - UnitOrDefault(record.direction) * (record.refRadius / tan_angle);
      auto *quadric = new G4HalfSpaceQuadric(1, 0, 0,
                                             1, 0,
                                             -std::pow(tan_angle, 2),
                                             0, 0, 0,
                                             0);
      quadric->Rotate(RotationFromZAxis(record.direction));
      quadric->Translate(apex);
      return quadric;
    }

    case SurfaceType::Unsupported:
    default:
      return nullptr;
  }
}
