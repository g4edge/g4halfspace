#pragma once

#include <vector>
#include <string>
#include <map>

#include "G4VHalfSpaceReader.hh"

class G4VHalfSpace;
class G4HalfSpaceSolid;
class G4HalfSpaceTransformation;
class G4HalfSpaceTest;

class G4HalfSpaceReader : public G4VHalfSpaceReader {
public:
  G4HalfSpaceReader() = default;
  G4HalfSpaceReader(const G4String &file_name);
  ~G4HalfSpaceReader() = default;
  G4HalfSpaceSolid* GetSolid(size_t region) override;
  G4HalfSpaceSolid* GetSolid(const G4String &region) override;
  G4HalfSpaceTest* GetTest(size_t test);

protected:
  void Load(const G4String &file_name) override;

  std::map<size_t, G4HalfSpaceTransformation*> hs_trans_map;
  std::map<size_t, G4VHalfSpace*> hs_surface_map;
  std::map<size_t, G4HalfSpaceSolid*> hs_solid_map;
  std::map<size_t, G4HalfSpaceTest*> hs_test_map;
};