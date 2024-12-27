#pragma once

#include <vector>
#include <string>
#include <map>

class G4VHalfSpace;
class G4HalfSpaceSolid;

std::vector<std::string> split(const std::string &s, char delim);

class G4HalfSpaceReader {
public:
  G4HalfSpaceReader(const G4String &file_name);
  ~G4HalfSpaceReader() = default;
  G4HalfSpaceSolid* GetSolid(size_t region);

protected:
  void Read(const G4String &file_name);

  std::map<size_t, G4VHalfSpace*> hs_surface_map;
  std::map<size_t, G4HalfSpaceSolid*> hs_solid_map;
};