#pragma once

#include "G4VHalfSpaceReader.hh"

#include "utils.hh"

#include <string>
#include <variant>
#include <vector>
#include <map>

struct Card {

};

struct RotoTranslation {
  std::variant<int, std::string> id;
};

struct Body {
  std::string type;
  std::variant<int, std::string> id;
};

struct Zone {
  std::vector<std::variant<int, std::string>> positive_surfaces;
  std::vector<std::variant<int, std::string>> negative_surfaces;
};

struct Region {
  std::variant<int, std::string> id;
  std::vector<Zone> zones;
};

struct Assignma {
  std::variant<int, std::string> region;
  std::string material;
};

G4HalfSpaceSolid* RegionToSolid(const Region& region) {
  return nullptr;
}

class G4FlukaReader : public G4VHalfSpaceReader {
public:
  G4FlukaReader();
  G4FlukaReader(const G4String &file_name);
  ~G4FlukaReader();

  virtual G4HalfSpaceSolid* GetSolid(size_t region) override;
  virtual G4HalfSpaceSolid* GetSolid(const G4String &region) override;


protected:

  void Load(const G4String &file_name) override;

  // states of loader
  bool free = false;
  int transform = -1;
  bool geom = false;
  bool pp_include = true;

  std::map<std::variant<int, std::string>, RotoTranslation> rototranslations;
  std::map<std::variant<int, std::string>, Body> bodies;
  std::map<std::variant<int, std::string>, Region> regions;
  std::map<std::variant<int, std::string>, Assignma> assignmas;

};