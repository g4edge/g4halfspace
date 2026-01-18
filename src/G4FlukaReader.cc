#include "G4FlukaReader.hh"

#include <string>
#include <fstream>
#include <sstream>
#include <vector>

#include "G4ios.hh"

G4FlukaReader::G4FlukaReader() : G4VHalfSpaceReader() {}

G4FlukaReader::G4FlukaReader(const G4String &file_name) : G4VHalfSpaceReader() {
  this->Load(file_name);
}

G4FlukaReader::~G4FlukaReader() {}

G4HalfSpaceSolid* G4FlukaReader::GetSolid(size_t region) {

  std::size_t currentIndex = 0;

  for (const auto& [key, value] : regions) {
    if (currentIndex == region) {
      return RegionToSolid(value);
    }
    ++currentIndex;
  };

  return nullptr;
}

G4HalfSpaceSolid* G4FlukaReader::GetSolid(const G4String &region) {
  return RegionToSolid(regions[region]);
}

void G4FlukaReader::Load(const G4String &file_name) {
  auto file = std::ifstream(file_name);

  for(std::string line; std::getline(file, line); ) {
    // print line
    G4cout << "G4FlukaReader::Load> line " << line << G4endl;

    // split to find card
    auto line_split = split(line, ' ');

    // empty line
    if(line_split.size() == 0) {
      continue;
    }

    // comment line
    if(line_split[0][0] == '*') {
      continue;
    }

  }
}