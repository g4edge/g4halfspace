#pragma once

#include <cstddef>
#include <map>
#include <string>
#include <vector>

#include "G4String.hh"
#include "G4VHalfSpaceReader.hh"

class G4VHalfSpace;
class G4HalfSpaceSolid;

class G4FlukaReader : public G4VHalfSpaceReader {
public:
  G4FlukaReader();
  G4FlukaReader(const G4String& file_name);
  ~G4FlukaReader() override;
  G4FlukaReader(const G4FlukaReader&) = delete;
  G4FlukaReader& operator=(const G4FlukaReader&) = delete;

  G4HalfSpaceSolid* GetSolid(size_t region) override;
  G4HalfSpaceSolid* GetSolid(const G4String& region) override;

protected:
  void Load(const G4String& file_name) override;

private:
  void ClearOwnedData();
  G4VHalfSpace* BuildBody(const std::string& type,
                          const std::vector<double>& values) const;

  std::map<std::string, G4VHalfSpace*> body_map;
  std::map<std::string, G4HalfSpaceSolid*> region_map;
  std::vector<std::string> region_order;
};
