#include "G4FlukaReader.hh"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "G4ios.hh"

#include "G4HalfSpaceAARBox.hh"
#include "G4HalfSpaceArbitrary.hh"
#include "G4HalfSpaceCircularCone.hh"
#include "G4HalfSpaceCircularCylinder.hh"
#include "G4HalfSpaceEllipsoid.hh"
#include "G4HalfSpaceEllipticCylinder.hh"
#include "G4HalfSpacePlane.hh"
#include "G4HalfSpaceQuadric.hh"
#include "G4HalfSpaceRBox.hh"
#include "G4HalfSpaceSolid.hh"
#include "G4HalfSpaceWedge.hh"
#include "G4HalfSpaceXACircularCylinder.hh"
#include "G4HalfSpaceXAEllipticalCylinder.hh"
#include "G4HalfSpaceXYPlane.hh"
#include "G4HalfSpaceXZPlane.hh"
#include "G4HalfSpaceYACircularCylinder.hh"
#include "G4HalfSpaceYAEllipticalCylinder.hh"
#include "G4HalfSpaceYZPlane.hh"
#include "G4HalfSpaceZACircularCylinder.hh"
#include "G4HalfSpaceZAEllipticalCylinder.hh"
#include "G4HalfSpaceZone.hh"
#include "G4HalfSpaceSphere.hh"

namespace {
std::string Trim(const std::string& s) {
  const auto begin = s.find_first_not_of(" \t\r\n");
  if (begin == std::string::npos) {
    return "";
  }
  const auto end = s.find_last_not_of(" \t\r\n");
  return s.substr(begin, end - begin + 1);
}

std::string ToUpper(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) {
    return static_cast<char>(std::toupper(c));
  });
  return s;
}

bool IsCommentOrEmpty(const std::string& line) {
  const auto trimmed = Trim(line);
  return trimmed.empty() || trimmed[0] == '*';
}

std::vector<std::string> Tokenize(const std::string& line) {
  std::istringstream stream(line);
  std::vector<std::string> tokens;
  std::string token;
  while (stream >> token) {
    tokens.push_back(token);
  }
  return tokens;
}

std::vector<std::vector<std::pair<char, std::string>>> ParseRegionExpression(
    const std::string& expression) {
  std::vector<std::vector<std::pair<char, std::string>>> zones;
  zones.emplace_back();

  for (std::size_t i = 0; i < expression.size();) {
    const char c = expression[i];
    if (c == '|') {
      zones.emplace_back();
      ++i;
      continue;
    }

    if (c == '+' || c == '-') {
      ++i;
      while (i < expression.size() &&
             std::isspace(static_cast<unsigned char>(expression[i]))) {
        ++i;
      }

      const std::size_t start = i;
      while (i < expression.size()) {
        const char idChar = expression[i];
        if (std::isalnum(static_cast<unsigned char>(idChar)) || idChar == '_' ||
            idChar == '.' || idChar == '$') {
          ++i;
          continue;
        }
        break;
      }

      if (i > start) {
        zones.back().emplace_back(c, expression.substr(start, i - start));
      }
      continue;
    }

    ++i;
  }

  zones.erase(
      std::remove_if(
          zones.begin(), zones.end(),
          [](const std::vector<std::pair<char, std::string>>& zone) {
            return zone.empty();
          }),
      zones.end());

  return zones;
}

std::vector<double> ParseDoubles(const std::vector<std::string>& tokens,
                                 std::size_t startIndex) {
  std::vector<double> values;
  values.reserve(tokens.size() - startIndex);
  for (std::size_t i = startIndex; i < tokens.size(); ++i) {
    values.push_back(std::stod(tokens[i]));
  }
  return values;
}
}  // namespace

G4FlukaReader::G4FlukaReader() : G4VHalfSpaceReader() {}

G4FlukaReader::G4FlukaReader(const G4String& file_name) : G4VHalfSpaceReader() {
  Load(file_name);
}

G4FlukaReader::~G4FlukaReader() = default;

G4HalfSpaceSolid* G4FlukaReader::GetSolid(size_t region) {
  if (region >= region_order.size()) {
    return nullptr;
  }
  return region_map[region_order[region]];
}

G4HalfSpaceSolid* G4FlukaReader::GetSolid(const G4String& region) {
  const auto it = region_map.find(region);
  if (it == region_map.end()) {
    return nullptr;
  }
  return it->second;
}

void G4FlukaReader::Load(const G4String& file_name) {
  std::ifstream file(file_name);
  if (!file.good()) {
    G4cout << "G4FlukaReader::Load could not open file " << file_name << G4endl;
    return;
  }

  enum class ParseState { OutsideGeometry, GeometryHeader, Bodies, Regions };
  ParseState state = ParseState::OutsideGeometry;

  for (std::string rawLine; std::getline(file, rawLine);) {
    if (IsCommentOrEmpty(rawLine)) {
      continue;
    }

    const std::string line = Trim(rawLine);
    const auto tokens = Tokenize(line);
    if (tokens.empty()) {
      continue;
    }

    const std::string card = ToUpper(tokens[0]);

    if (card == "GEOBEGIN") {
      state = ParseState::GeometryHeader;
      continue;
    }

    if (card == "GEOEND") {
      state = ParseState::OutsideGeometry;
      continue;
    }

    if (state == ParseState::OutsideGeometry) {
      continue;
    }

    if (state == ParseState::GeometryHeader) {
      state = ParseState::Bodies;
      continue;
    }

    if (card == "END") {
      if (state == ParseState::Bodies) {
        state = ParseState::Regions;
      } else if (state == ParseState::Regions) {
        state = ParseState::OutsideGeometry;
      }
      continue;
    }

    if (state == ParseState::Bodies) {
      if (tokens.size() < 3) {
        continue;
      }

      try {
        const auto values = ParseDoubles(tokens, 2);
        auto* body = BuildBody(card, values);
        if (body != nullptr) {
          body_map[tokens[1]] = body;
        }
      } catch (const std::exception& error) {
        G4cout << "G4FlukaReader::Load skipping malformed body line: " << line
               << " (" << error.what() << ")" << G4endl;
      }
      continue;
    }

    if (state == ParseState::Regions) {
      if (tokens.size() < 3) {
        continue;
      }

      const std::string regionName = tokens[0];
      auto* solid = new G4HalfSpaceSolid(regionName);
      const auto zoneTerms = ParseRegionExpression(line);
      bool regionValid = true;

      for (const auto& zoneTerm : zoneTerms) {
        auto* zone = new G4HalfSpaceZone();
        bool zoneValid = true;
        for (const auto& [sign, bodyName] : zoneTerm) {
          const auto it = body_map.find(bodyName);
          if (it == body_map.end()) {
            G4cout << "G4FlukaReader::Load unknown body '" << bodyName
                   << "' in region " << regionName << G4endl;
            zoneValid = false;
            regionValid = false;
            break;
          }

          if (sign == '+') {
            zone->AddIntersection(it->second);
          } else if (sign == '-') {
            zone->AddSubtraction(it->second);
          }
        }
        if (zoneValid) {
          solid->AddZone(zone);
        }
      }

      if (!regionValid || zoneTerms.empty()) {
        continue;
      }

      const bool regionExists = region_map.find(regionName) != region_map.end();
      region_map[regionName] = solid;
      if (!regionExists) {
        region_order.push_back(regionName);
      }
    }
  }
}

G4VHalfSpace* G4FlukaReader::BuildBody(const std::string& type,
                                       const std::vector<double>& values) const {
  if (type == "RPP" && values.size() >= 6) {
    return new G4HalfSpaceAARBox(values[0], values[1], values[2], values[3],
                                 values[4], values[5]);
  }

  if (type == "BOX" && values.size() >= 12) {
    return new G4HalfSpaceRBox(G4ThreeVector(values[0], values[1], values[2]),
                               G4ThreeVector(values[3], values[4], values[5]),
                               G4ThreeVector(values[6], values[7], values[8]),
                               G4ThreeVector(values[9], values[10], values[11]));
  }

  if (type == "SPH" && values.size() >= 4) {
    return new G4HalfSpaceSphere(values[0], values[1], values[2], values[3]);
  }

  if (type == "TRC" && values.size() >= 8) {
    return new G4HalfSpaceCircularCone(
        G4ThreeVector(values[0], values[1], values[2]),
        G4ThreeVector(values[3], values[4], values[5]), values[6], values[7]);
  }

  if (type == "ELL" && values.size() >= 7) {
    return new G4HalfSpaceEllipsoid(
        G4ThreeVector(values[0], values[1], values[2]),
        G4ThreeVector(values[3], values[4], values[5]), values[6]);
  }

  if ((type == "WED" || type == "RAW") && values.size() >= 12) {
    return new G4HalfSpaceWedge(G4ThreeVector(values[0], values[1], values[2]),
                                G4ThreeVector(values[3], values[4], values[5]),
                                G4ThreeVector(values[6], values[7], values[8]),
                                G4ThreeVector(values[9], values[10], values[11]));
  }

  if (type == "ARB" && values.size() >= 30) {
    return new G4HalfSpaceArbitrary(
        G4ThreeVector(values[0], values[1], values[2]),
        G4ThreeVector(values[3], values[4], values[5]),
        G4ThreeVector(values[6], values[7], values[8]),
        G4ThreeVector(values[9], values[10], values[11]),
        G4ThreeVector(values[12], values[13], values[14]),
        G4ThreeVector(values[15], values[16], values[17]),
        G4ThreeVector(values[18], values[19], values[20]),
        G4ThreeVector(values[21], values[22], values[23]),
        static_cast<G4int>(values[24]), static_cast<G4int>(values[25]),
        static_cast<G4int>(values[26]), static_cast<G4int>(values[27]),
        static_cast<G4int>(values[28]), static_cast<G4int>(values[29]));
  }

  if (type == "PLA" && values.size() >= 6) {
    return new G4HalfSpacePlane(G4ThreeVector(values[0], values[1], values[2]),
                                G4ThreeVector(values[3], values[4], values[5]));
  }

  if (type == "XYP" && values.size() >= 1) {
    return new G4HalfSpaceXYPlane(values[0]);
  }

  if (type == "XZP" && values.size() >= 1) {
    return new G4HalfSpaceXZPlane(values[0]);
  }

  if (type == "YZP" && values.size() >= 1) {
    return new G4HalfSpaceYZPlane(values[0]);
  }

  if (type == "RCC" && values.size() >= 7) {
    return new G4HalfSpaceCircularCylinder(
        G4ThreeVector(values[0], values[1], values[2]),
        G4ThreeVector(values[3], values[4], values[5]), values[6]);
  }

  if (type == "XCC" && values.size() >= 3) {
    return new G4HalfSpaceXACircularCylinder(values[0], values[1], values[2]);
  }

  if (type == "YCC" && values.size() >= 3) {
    return new G4HalfSpaceYACircularCylinder(values[0], values[1], values[2]);
  }

  if (type == "ZCC" && values.size() >= 3) {
    return new G4HalfSpaceZACircularCylinder(values[0], values[1], values[2]);
  }

  if (type == "REC" && values.size() >= 12) {
    return new G4HalfSpaceEllipticCylinder(
        G4ThreeVector(values[0], values[1], values[2]),
        G4ThreeVector(values[3], values[4], values[5]),
        G4ThreeVector(values[6], values[7], values[8]),
        G4ThreeVector(values[9], values[10], values[11]));
  }

  if (type == "XEC" && values.size() >= 4) {
    return new G4HalfSpaceXAEllipticalCylinder(values[0], values[1], values[2],
                                               values[3]);
  }

  if (type == "YEC" && values.size() >= 4) {
    return new G4HalfSpaceYAEllipticalCylinder(values[0], values[1], values[2],
                                               values[3]);
  }

  if (type == "ZEC" && values.size() >= 4) {
    return new G4HalfSpaceZAEllipticalCylinder(values[0], values[1], values[2],
                                               values[3]);
  }

  if (type == "QUA" && values.size() >= 10) {
    return new G4HalfSpaceQuadric(values[0], values[1], values[2], values[3],
                                  values[4], values[5], values[6], values[7],
                                  values[8], values[9]);
  }

  G4cout << "G4FlukaReader::Load unsupported body type or malformed card: "
         << type << G4endl;
  return nullptr;
}
