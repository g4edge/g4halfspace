//
// Created by Boogert Stewart on 16/12/2024.
//

#include <fstream>
#include <sstream>

#include "G4String.hh"

#include "G4HalfSpaceReader.hh"
#include "G4HalfSpaceSolid.hh"
#include "G4HalfSpaceZone.hh"
#include "G4HalfSpacePlane.hh"
#include "G4HalfSpaceAARBox.hh"
#include "G4HalfSpaceRBox.hh"
#include "G4HalfSpaceSphere.hh"
#include "G4HalfSpaceEllipsoid.hh"
#include "G4HalfSpaceXACircularCylinder.hh"
#include "G4HalfSpaceYACircularCylinder.hh"
#include "G4HalfSpaceZACircularCylinder.hh"
#include "G4HalfSpaceXAEllipticalCylinder.hh"
#include "G4HalfSpaceYAEllipticalCylinder.hh"
#include "G4HalfSpaceZAEllipticalCylinder.hh"


std::vector<std::string> split(const std::string &s, char delim) {
  std::vector<std::string> result;
  std::stringstream ss (s);
  std::string item;

  while (getline (ss, item, delim)) {
    if(item != " ")
      result.push_back (item);
  }

  return result;
}

G4HalfSpaceReader::G4HalfSpaceReader(const G4String &file_name) {
  this->Read(file_name);
}

G4HalfSpaceSolid* G4HalfSpaceReader::GetSolid(size_t region) {
  return hs_solid_map[region];
}

void G4HalfSpaceReader::Read(const G4String &file_name) {
  std::ifstream ifstr(file_name);
  std::stringstream sstr;
  std::string key;

  while(ifstr >> key) {
    std::cout << key << std::endl;

    if(key == "aarbox") {
      size_t surface_id;
      double xmin, xmax, ymin, ymax, zmin, zmax;
      ifstr >> surface_id >> xmin >> xmax >> ymin >> ymax >> zmin >> zmax;
      hs_surface_map[surface_id] = new G4HalfSpaceAARBox(xmin, xmax,
                                                         ymin, ymax,
                                                         zmin, zmax);

    }
    else if(key == "rbox") {
      size_t surface_id;
      double vx, vy, vz, dx, dy, dz, rx, ry, rz;
      ifstr >> surface_id >> vx >> vy >> vz >> dx >> dy >> dz >> rx >> ry >> rz;
      hs_surface_map[surface_id] = new G4HalfSpaceRBox(G4ThreeVector(vx,vy,vz),
                                                       G4ThreeVector(dx,dy,dz),
                                                       G4ThreeVector(rx,ry,rz));
    }
    else if(key == "sphere") {
      size_t surface_id;
      double radius, xcentre, ycentre, zcentre;
      ifstr >> surface_id >> radius >> xcentre >> ycentre >> zcentre;
      hs_surface_map[surface_id] = new G4HalfSpaceSphere(radius, G4ThreeVector(xcentre,
                                                                               ycentre,
                                                                               zcentre));
    }
    else if(key == "ellipsoid") {
      size_t surface_id;
      double a, b, c, xcentre, ycentre, zcentre, xrotation, yrotation, zrotation;
      ifstr >> surface_id >> a >> b >> c >> xcentre >> ycentre >> zcentre >> xrotation >> yrotation >> zrotation;
      hs_surface_map[surface_id] = new G4HalfSpaceEllipsoid(a,b,c,
                                                            G4ThreeVector(xcentre,ycentre,zcentre),
                                                            G4ThreeVector(xrotation, yrotation, zrotation));
    }
    else if(key == "xacc") {
      size_t surface_id;
      double ycentre, zcentre, radius;
      ifstr >> surface_id >> ycentre >> zcentre >> radius;
      hs_surface_map[surface_id] = new G4HalfSpaceXACircularCylinder(ycentre, zcentre, radius);
    }
    else if(key == "yacc") {
      size_t surface_id;
      double xcentre, zcentre, radius;
      ifstr >> surface_id >> xcentre >> zcentre >> radius;
      hs_surface_map[surface_id] = new G4HalfSpaceYACircularCylinder(xcentre, zcentre, radius);
    }
    else if(key == "zacc") {
      size_t surface_id;
      double xcentre, ycentre, radius;
      ifstr >> surface_id >> xcentre >> ycentre >> radius;
      hs_surface_map[surface_id] = new G4HalfSpaceZACircularCylinder(xcentre, ycentre, radius);
    }
    else if(key == "xaec") {
      size_t surface_id;
      double ycentre, zcentre, yradius, zradius;
      ifstr >> surface_id >> ycentre >> zcentre >> yradius >> zradius;
      hs_surface_map[surface_id] = new G4HalfSpaceXAEllipticalCylinder(ycentre, zcentre, yradius, zradius);
    }
    else if(key == "yaec") {
      size_t surface_id;
      double xcentre, zcentre, xradius, zradius;
      ifstr >> surface_id >> xcentre >> zcentre >> xradius >> zradius;
      hs_surface_map[surface_id] = new G4HalfSpaceYAEllipticalCylinder(xcentre, zcentre, xradius, zradius);
    }
    else if(key == "zaec") {
      size_t surface_id;
      double xcentre, ycentre, xradius, yradius;
      ifstr >> surface_id >> xcentre >> ycentre >> xradius >> yradius;
      hs_surface_map[surface_id] = new G4HalfSpaceZAEllipticalCylinder(xcentre, ycentre, xradius, yradius);
    }
    else if(key == "region") {
      size_t region_id;
      std::string region_boolean;
      ifstr >> region_id;

      std::getline(ifstr, region_boolean);

      std::cout << "region_id=" << region_id << " " << region_boolean << std::endl;

      auto zone = split(region_boolean, '|');

      auto hs_solid = new G4HalfSpaceSolid(std::string("solid_")+std::to_string(region_id));
      hs_solid_map[region_id] = hs_solid;

      for (auto z : zone) {
        std::cout << "zone=" << z << std::endl;

        // create zone
        auto zone = new G4HalfSpaceZone();

        // create booleans
        for (auto boo: split(z, ' ')) {
          std::cout << "bool=" << boo << std::endl;

          if (boo[0] == '+') {
            auto surf_id = std::stoi(boo.substr(1));
            zone->AddIntersection(hs_surface_map[surf_id]);
          } else if (boo[0] == '-') {
            auto surf_id = std::stoi(boo.substr(1));
            zone->AddSubtraction(hs_surface_map[surf_id]);
          }
        }

        hs_solid->AddZone(zone);
      }
      std::cout << std::endl;
    }
  }
}
