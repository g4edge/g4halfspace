#pragma once

#include <cstddef>
#include <map>

#include "G4String.hh"
#include "G4ThreeVector.hh"

class G4VHalfSpace;

class G4LoaderHalfSpace {
public:
  enum class SurfaceType {
    Plane,
    Cylinder,
    Cone,
    Sphere,
    Unsupported
  };

  struct SurfaceRecord {
    SurfaceType type = SurfaceType::Unsupported;
    G4ThreeVector location = G4ThreeVector();
    G4ThreeVector direction = G4ThreeVector(0, 0, 1);
    G4double radius = 0.0;
    G4double refRadius = 0.0;
    G4double semiAngle = 0.0;
    bool reversed = false;
    bool orientationSupportsHalfSpace = true;
  };

  G4LoaderHalfSpace() = default;
  ~G4LoaderHalfSpace() = default;

  bool LoadStep(const G4String &file_name);
  void Clear();

  const std::map<std::size_t, SurfaceRecord>& GetSurfaces() const;
  std::map<std::size_t, G4VHalfSpace*> CreateHalfSpaces() const;

private:
  G4VHalfSpace* CreateHalfSpace(const SurfaceRecord &record) const;

  std::map<std::size_t, SurfaceRecord> surfaces_;
};
