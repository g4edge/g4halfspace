#include "G4FlukaReader.hh"

#include <algorithm>
#include <cctype>
#include <exception>
#include <fstream>
#include <sstream>
#include <vector>

#include "G4ios.hh"

namespace {

std::string trim(const std::string& input) {
  auto begin = input.begin();
  auto end = input.end();

  while (begin != end && std::isspace(static_cast<unsigned char>(*begin))) {
    ++begin;
  }
  while (begin != end && std::isspace(static_cast<unsigned char>(*(end - 1)))) {
    --end;
  }

  return std::string(begin, end);
}

std::vector<std::string> splitWhitespace(const std::string& line) {
  std::istringstream stream(line);
  std::string token;
  std::vector<std::string> tokens;
  while (stream >> token) {
    tokens.push_back(token);
  }
  return tokens;
}

bool appendNumericParameters(const std::vector<std::string>& tokens, size_t start, std::vector<double>& values) {
  std::vector<double> parsedValues;
  for (size_t index = start; index < tokens.size(); ++index) {
    try {
      parsedValues.push_back(std::stod(tokens[index]));
    }
    catch (const std::exception&) {
      return false;
    }
  }

  values.insert(values.end(), parsedValues.begin(), parsedValues.end());
  return true;
}

std::vector<Zone> expressionTokensToZones(const std::vector<std::string>& expressionTokens) {
  std::vector<Zone> zones;
  Zone currentZone;

  for (const auto& token : expressionTokens) {
    if (token == "|") {
      if (!currentZone.positive_surfaces.empty() || !currentZone.negative_surfaces.empty()) {
        zones.push_back(currentZone);
      }
      currentZone = Zone();
      continue;
    }

    auto cleanedToken = token;
    cleanedToken.erase(std::remove(cleanedToken.begin(), cleanedToken.end(), '('), cleanedToken.end());
    cleanedToken.erase(std::remove(cleanedToken.begin(), cleanedToken.end(), ')'), cleanedToken.end());
    if (cleanedToken.empty()) {
      continue;
    }

    const auto sign = cleanedToken[0];
    if (sign != '+' && sign != '-') {
      continue;
    }
    const auto bodyName = cleanedToken.substr(1);
    if (bodyName.empty()) {
      continue;
    }

    if (sign == '+') {
      currentZone.positive_surfaces.push_back(bodyName);
    }
    else {
      currentZone.negative_surfaces.push_back(bodyName);
    }
  }

  if (!currentZone.positive_surfaces.empty() || !currentZone.negative_surfaces.empty()) {
    zones.push_back(currentZone);
  }

  return zones;
}

std::string joinTokens(const std::vector<std::string>& tokens) {
  std::ostringstream stream;
  for (size_t i = 0; i < tokens.size(); ++i) {
    if (i != 0) {
      stream << ' ';
    }
    stream << tokens[i];
  }
  return stream.str();
}

bool isRegionContinuation(const std::string& token) {
  if (token.empty()) {
    return false;
  }

  const auto first = token[0];
  return first == '+' || first == '-' || first == '|' || first == '(' || first == ')';
}

bool isNumericToken(const std::string& token) {
  if (token.empty()) {
    return false;
  }

  try {
    std::stod(token);
    return true;
  }
  catch (const std::exception&) {
    return false;
  }
}

}

G4FlukaReader::G4FlukaReader() : G4VHalfSpaceReader() {}

G4FlukaReader::G4FlukaReader(const G4String &file_name) : G4VHalfSpaceReader() {
  this->Load(file_name);
}

G4FlukaReader::~G4FlukaReader() {}

G4HalfSpaceSolid* G4FlukaReader::GetSolid(size_t region) {
  if (region >= region_order.size()) {
    return nullptr;
  }

  const auto& regionName = region_order[region];
  auto regionIt = regions.find(regionName);
  if (regionIt == regions.end()) {
    return nullptr;
  }

  return RegionToSolid(regionIt->second);
}

G4HalfSpaceSolid* G4FlukaReader::GetSolid(const G4String &region) {
  auto regionIt = regions.find(region);
  if (regionIt == regions.end()) {
    return nullptr;
  }
  return RegionToSolid(regionIt->second);
}

void G4FlukaReader::Load(const G4String &file_name) {
  enum class GeometrySection {
    kNone,
    kBodies,
    kRegions,
  };

  rototranslations.clear();
  bodies.clear();
  regions.clear();
  assignmas.clear();
  region_order.clear();

  auto file = std::ifstream(file_name);
  if (!file.is_open()) {
    return;
  }

  auto section = GeometrySection::kNone;
  std::string currentRegion;

  for (std::string line; std::getline(file, line); ) {
    auto trimmedLine = trim(line);

    if (trimmedLine.empty()) {
      continue;
    }

    if (trimmedLine[0] == '*') {
      continue;
    }

    const auto lineTokens = splitWhitespace(trimmedLine);
    if (lineTokens.empty()) {
      continue;
    }

    const auto& card = lineTokens[0];

    if (card == "GEOBEGIN") {
      geom = true;
      section = GeometrySection::kBodies;
      currentRegion.clear();
      continue;
    }

    if (card == "GEOEND") {
      geom = false;
      section = GeometrySection::kNone;
      currentRegion.clear();
      continue;
    }

    if (card == "END") {
      if (section == GeometrySection::kBodies) {
        section = GeometrySection::kRegions;
      }
      else if (section == GeometrySection::kRegions) {
        section = GeometrySection::kNone;
      }
      currentRegion.clear();
      continue;
    }

    if (section == GeometrySection::kBodies) {
      if (lineTokens.size() < 2) {
        continue;
      }

      Body body;
      body.type = lineTokens[0];
      body.id = lineTokens[1];
      if (!appendNumericParameters(lineTokens, 2, body.parameters)) {
        continue;
      }
      bodies[body.id] = body;
      continue;
    }

    if (section == GeometrySection::kRegions) {
      const auto expressionStart =
          (lineTokens.size() >= 2 && isNumericToken(lineTokens[1])) ? 2 : 1;
      if (lineTokens.size() > expressionStart) {
        Region region;
        region.id = lineTokens[0];
        region.expression_tokens.assign(lineTokens.begin() + expressionStart, lineTokens.end());
        region.raw_expression = joinTokens(region.expression_tokens);
        region.zones = expressionTokensToZones(region.expression_tokens);
        const auto isNewRegion = regions.find(region.id) == regions.end();
        regions[region.id] = region;
        if (isNewRegion) {
          region_order.push_back(region.id);
        }
        currentRegion = region.id;
        continue;
      }

      if (!currentRegion.empty() && isRegionContinuation(lineTokens[0])) {
        auto regionIt = regions.find(currentRegion);
        if (regionIt == regions.end()) {
          continue;
        }

        regionIt->second.expression_tokens.insert(regionIt->second.expression_tokens.end(),
                                                  lineTokens.begin(),
                                                  lineTokens.end());
        regionIt->second.raw_expression = joinTokens(regionIt->second.expression_tokens);
        regionIt->second.zones = expressionTokensToZones(regionIt->second.expression_tokens);
      }
      continue;
    }

    if (card == "ROT-DEFI") {
      if (lineTokens.size() < 2) {
        continue;
      }

      RotoTranslation rotoTranslation;
      rotoTranslation.id = lineTokens[1];
      if (!appendNumericParameters(lineTokens, 2, rotoTranslation.parameters)) {
        continue;
      }
      rototranslations[rotoTranslation.id] = rotoTranslation;
      continue;
    }

    if (card == "ASSIGNMA" && lineTokens.size() >= 3) {
      Assignma assignma;
      assignma.material = lineTokens[1];
      assignma.region = lineTokens[2];
      assignmas[assignma.region] = assignma;
      continue;
    }
  }
}