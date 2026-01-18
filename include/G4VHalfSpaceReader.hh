#pragma once

#include <cstddef>

#include "G4String.hh"

class G4HalfSpaceSolid;

class G4VHalfSpaceReader {
public:
  G4VHalfSpaceReader() = default;
  virtual ~G4VHalfSpaceReader() = default;
  virtual G4HalfSpaceSolid* GetSolid(size_t region) = 0;
  virtual G4HalfSpaceSolid* GetSolid(const G4String &region) = 0;

protected:
  virtual void Load(const G4String &file_name) = 0;
};