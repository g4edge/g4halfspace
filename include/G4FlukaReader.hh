#pragma once

#include "G4VHalfSpaceReader.hh"

#include "utils.hh"

#include <map>
#include <string>
#include <vector>

struct Card {

};

struct RotoTranslation {
  std::string id;
  std::vector<double> parameters;
};

struct Body {
  std::string type;
  std::string id;
  std::vector<double> parameters;
};

struct Zone {
  std::vector<std::string> positive_surfaces;
  std::vector<std::string> negative_surfaces;
};

struct Region {
  std::string id;
  std::string raw_expression;
  std::vector<std::string> expression_tokens;
  std::vector<Zone> zones;
};

struct Assignma {
  std::string region;
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

  const std::map<std::string, RotoTranslation>& GetRotoTranslations() const { return rototranslations; }
  const std::map<std::string, Body>& GetBodies() const { return bodies; }
  const std::map<std::string, Region>& GetRegions() const { return regions; }
  const std::vector<std::string>& GetRegionOrder() const { return region_order; }
  const std::map<std::string, Assignma>& GetAssignmas() const { return assignmas; }


protected:

  void Load(const G4String &file_name) override;

  // states of loader
  bool free = false;
  int transform = -1;
  bool geom = false;
  bool pp_include = true;

  std::map<std::string, RotoTranslation> rototranslations;
  std::map<std::string, Body> bodies;
  std::map<std::string, Region> regions;
  std::vector<std::string> region_order;
  std::map<std::string, Assignma> assignmas;

};