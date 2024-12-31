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
#include "G4HalfSpaceXYPlane.hh"
#include "G4HalfSpaceXZPlane.hh"
#include "G4HalfSpaceYZPlane.hh"
#include "G4HalfSpaceTransformation.hh"
#include "G4HalfSpaceAARBox.hh"
#include "G4HalfSpaceRBox.hh"
#include "G4HalfSpaceSphere.hh"
#include "G4HalfSpaceCircularCone.hh"
#include "G4HalfSpaceEllipsoid.hh"
#include "G4HalfSpaceWedge.hh"
#include "G4HalfSpaceArbitrary.hh"
#include "G4HalfSpaceCircularCylinder.hh"
#include "G4HalfSpaceXACircularCylinder.hh"
#include "G4HalfSpaceYACircularCylinder.hh"
#include "G4HalfSpaceZACircularCylinder.hh"
#include "G4HalfSpaceEllipticCylinder.hh"
#include "G4HalfSpaceXAEllipticalCylinder.hh"
#include "G4HalfSpaceYAEllipticalCylinder.hh"
#include "G4HalfSpaceZAEllipticalCylinder.hh"
#include "G4HalfSpaceQuadric.hh"
#include "G4HalfSpacePlane.hh"

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

  size_t iLineNumber = 0;
  while(ifstr >> key) {
    iLineNumber++;

    std::cout << key << std::endl;

    if(key[0] == '#') {
      std::string restOfLine;
      std::getline(ifstr, restOfLine);
    }

    if(key == "trans") {
      size_t trans_id;
      double dx, dy, dz, rx, ry, rz;
      ifstr >> trans_id >> dx >> dy >> dz >> rx >> ry >> rz;
      hs_trans_map[trans_id] = new G4HalfSpaceTransformation(G4ThreeVector(dx,dy,dz),
                                                             G4ThreeVector(rx/180.*M_PI,ry/180.*M_PI,rz/180.*M_PI));
    }
    else if(key == "aarbox") {
      size_t surface_id = -1;
      int trans_id = -1;
      double xmin, xmax, ymin, ymax, zmin, zmax;
      ifstr >> surface_id >> xmin >> xmax >> ymin >> ymax >> zmin >> zmax >> trans_id;
      if(trans_id < 0) {
        hs_surface_map[surface_id] = new G4HalfSpaceAARBox(xmin, xmax,
                                                           ymin, ymax,
                                                           zmin, zmax);
      }
      else {
        auto t = hs_trans_map[trans_id];
        auto solid = new G4HalfSpaceRBox(G4ThreeVector((xmax-xmin)/2, (ymax-ymin)/2, (zmax-zmin)/2),
                                         G4ThreeVector((xmax+xmin)/2, (ymax+ymin)/2, (zmax+zmin)/2),
                                         G4ThreeVector(0,0,0));
        solid->Rotate(t->GetRotationMatrix());
        solid->Translate(t->GetTranslation());
        hs_surface_map[surface_id] = solid;
      }
    }
    else if(key == "aarbox_od") {
      size_t surface_id;
      int trans_id;
      double dx, dy, dz, cx, cy, cz;
      ifstr >> surface_id >> dx >> dy >> dz >> cx >> cy >> cz >> trans_id;
      if(trans_id <0) {
        hs_surface_map[surface_id] = new G4HalfSpaceAARBox(G4ThreeVector(dx, dy, dz),
                                                           G4ThreeVector(cx, cy, cz));
      }
      else {
        auto t = hs_trans_map[trans_id];
        auto solid = new G4HalfSpaceRBox(G4ThreeVector(dx,dy,dz),
                                         G4ThreeVector(cx,cy,cz),
                                         G4ThreeVector(0,0,0));
        solid->Rotate(t->GetRotationMatrix());
        solid->Translate(t->GetTranslation());
        hs_surface_map[surface_id] = solid;
      }
    }
    else if(key == "rbox") {
      size_t surface_id;
      int trans_id;
      double vx, vy, vz, h1x, h1y, h1z, h2x, h2y, h2z, h3x, h3y, h3z;
      ifstr >> surface_id >> vx >> vy >> vz
                          >> h1x >> h1y >> h1z
                          >> h2x >> h2y >> h2z
                          >> h3x >> h3y >> h3z
                          >> trans_id;
      auto solid = new G4HalfSpaceRBox(G4ThreeVector(vx, vy, vz),
                                       G4ThreeVector(h1x, h1y, h1z),
                                       G4ThreeVector(h2x, h2y, h2z),
                                       G4ThreeVector(h3x, h3y, h3z));
      if (trans_id > -1) {
        auto t = hs_trans_map[trans_id];
        solid->Rotate(t->GetRotationMatrix());
        solid->Translate(t->GetTranslation());
      }
      hs_surface_map[surface_id] = solid;
    }
    else if(key == "rbox_od") {
      size_t surface_id;
      int trans_id;
      double dx, dy, dz, cx, cy, cz, rx, ry, rz;
      ifstr >> surface_id >> dx >> dy >> dz >> cx >> cy >> cz >> rx >> ry >> rz >> trans_id;
      auto solid = new G4HalfSpaceRBox(G4ThreeVector(dx, dy, dz),
                                       G4ThreeVector(cx, cy, cz),
                                       G4ThreeVector(rx / 180. * M_PI, ry / 180. * M_PI, rz / 180. * M_PI));
      if (trans_id > -1) {
        auto t = hs_trans_map[trans_id];
        solid->Rotate(t->GetRotationMatrix());
        solid->Translate(t->GetTranslation());
      }
      hs_surface_map[surface_id] = solid;
    }
    else if(key == "sphere") {
      size_t surface_id;
      int trans_id;
      double radius, xcentre, ycentre, zcentre;
      ifstr >> surface_id >> xcentre >> ycentre >> zcentre >> radius >> trans_id;
      auto solid = new G4HalfSpaceSphere(G4ThreeVector(xcentre,ycentre,zcentre),radius);
      if (trans_id > -1) {
        auto t = hs_trans_map[trans_id];
        solid->Translate(t->GetTranslation());
      }
      hs_surface_map[surface_id] = solid;
    }
    else if(key == "cone") {
      size_t surface_id;
      int trans_id;
      double vx, vy, vz, hx, hy, hz, r1, r2;
      ifstr >> surface_id >> vx >> vy >> vz >> hx >> hy >> hz >> r1 >> r2 >> trans_id;
      auto solid = new G4HalfSpaceCircularCone(G4ThreeVector(vx,vy,vz),
                                               G4ThreeVector(hx,hy,hz),
                                               r1,r2);
      if (trans_id > -1) {
        auto t = hs_trans_map[trans_id];
        solid->Rotate(t->GetRotationMatrix());
        solid->Translate(t->GetTranslation());
      }
      hs_surface_map[surface_id] = solid;
    }
    else if(key == "cone_od") {
      size_t surface_id;
      int trans_id;
      double h, r1, r2, vx, vy, vz, rx, ry, rz;
      ifstr >> surface_id >> h >> r1 >> r2 >> vx >> vy >> vz >> rx >> ry >> rz >> trans_id;
       auto solid = new G4HalfSpaceCircularCone(h, r1, r2,
                                                               G4ThreeVector(vx,vy,vz),
                                                               G4ThreeVector(rx,ry,rz));
      if (trans_id > -1) {
        auto t = hs_trans_map[trans_id];
        solid->Rotate(t->GetRotationMatrix());
        solid->Translate(t->GetTranslation());
      }
      hs_surface_map[surface_id] = solid;
    }
    else if(key == "ellipsoid") {
      size_t surface_id;
      int trans_id;
      double f1x, f1y, f1z, f2x, f2y, f2z, l;
      ifstr >> surface_id >> f1x >> f1y >> f1z >> f2x >> f2y >> f2z >> l >> trans_id;
      auto solid = new G4HalfSpaceEllipsoid(G4ThreeVector(f1x, f1y, f1z),
                                            G4ThreeVector(f2x, f2y, f2z),
                                            l);
      if (trans_id > -1) {
        auto t = hs_trans_map[trans_id];
        solid->Rotate(t->GetRotationMatrix());
        solid->Translate(t->GetTranslation());
      }
      hs_surface_map[surface_id] = solid;
    }
    else if(key == "ellipsoid_od") {
      size_t surface_id;
      int trans_id;
      double a, b, c, xcentre, ycentre, zcentre, xrotation, yrotation, zrotation;
      ifstr >> surface_id >> a >> b >> c >> xcentre >> ycentre >> zcentre >> xrotation >> yrotation >> zrotation >> trans_id;
      hs_surface_map[surface_id] = new G4HalfSpaceEllipsoid(G4ThreeVector(a,b,c),
                                                            G4ThreeVector(xcentre,ycentre,zcentre),
                                                            G4ThreeVector(xrotation/180*M_PI, yrotation/180*M_PI, zrotation/180*M_PI));
    }
    else if(key == "wedge" || key == "raw") {
      size_t surface_id;
      int trans_id;
      double vx, vy, vz, h1x, h1y, h1z, h2x, h2y, h2z, h3x, h3y, h3z;
      ifstr >> surface_id >> vx >> vy >> vz >> h1x >> h1y >> h1z >> h2x >> h2y >> h2z >> h3x >> h3y >> h3z >> trans_id;
      auto solid = new G4HalfSpaceWedge(G4ThreeVector(vx,vy,vz),
                                        G4ThreeVector(h1x, h1y, h1z),
                                        G4ThreeVector(h2x, h2y, h2z),
                                        G4ThreeVector(h3x, h3y, h3z));
      if (trans_id > -1) {
        auto t = hs_trans_map[trans_id];
        solid->Rotate(t->GetRotationMatrix());
        solid->Translate(t->GetTranslation());
      }
      hs_surface_map[surface_id] = solid;
    }
    else if(key == "wedge_od") {
      size_t surface_id;
      int trans_id;
      double dx, dy, dz, cx, cy, cz, rx, ry, rz;
      ifstr >> surface_id >> dx >> dy >> dz >> cx >> cy >> cz >> rx >> ry >> rz >> trans_id;
      auto solid = new G4HalfSpaceWedge(G4ThreeVector(dx,dy,dz),
                                        G4ThreeVector(cx,cy,cz),
                                        G4ThreeVector(rx/180.*M_PI,ry/180.*M_PI,rz/180.*M_PI));
      if (trans_id > -1) {
        auto t = hs_trans_map[trans_id];
        solid->Rotate(t->GetRotationMatrix());
        solid->Translate(t->GetTranslation());
      }
      hs_surface_map[surface_id] = solid;
    }
    else if(key == "arbitrary") {
      size_t surface_id;
      int trans_id;
      G4double v1x, v1y, v1z,
               v2x, v2y, v2z,
               v3x, v3y, v3z,
               v4x, v4y, v4z,
               v5x, v5y, v5z,
               v6x, v6y, v6z,
               v7x, v7y, v7z,
               v8x, v8y, v8z;
      G4int fi1, fi2, fi3, fi4, fi5, fi6;

      ifstr >> surface_id
            >> v1x >> v1y >> v1z
            >> v2x >> v2y >> v2z
            >> v3x >> v3y >> v3z
            >> v4x >> v4y >> v4z
            >> v5x >> v5y >> v5z
            >> v6x >> v6y >> v6z
            >> v7x >> v7y >> v7z
            >> v8x >> v8y >> v8z
            >> fi1 >> fi2 >> fi3
            >> fi4 >> fi5 >> fi6
            >> trans_id;

      auto solid = new G4HalfSpaceArbitrary(G4ThreeVector(v1x, v1y, v1z),
                                            G4ThreeVector(v2x, v2y, v2z),
                                            G4ThreeVector(v3x, v3y, v3z),
                                            G4ThreeVector(v4x, v4y, v4z),
                                            G4ThreeVector(v5x, v5y, v5z),
                                            G4ThreeVector(v6x, v6y, v6z),
                                            G4ThreeVector(v7x, v7y, v7z),
                                            G4ThreeVector(v8x, v8y, v8z),
                                            fi1, fi2, fi3, fi4, fi5, fi6);

      if (trans_id > -1) {
        auto t = hs_trans_map[trans_id];
        solid->Rotate(t->GetRotationMatrix());
        solid->Translate(t->GetTranslation());
      }
      hs_surface_map[surface_id] = solid;
    }
    else if(key == "arbitrary_od") {
      size_t surface_id;
      int trans_id;
      G4double cx, cy, cz,
          v1x, v1y, v1z,
          v2x, v2y, v2z,
          v3x, v3y, v3z,
          v4x, v4y, v4z,
          v5x, v5y, v5z,
          v6x, v6y, v6z,
          v7x, v7y, v7z,
          v8x, v8y, v8z;
      G4int fi1, fi2, fi3, fi4, fi5, fi6;

      ifstr >> surface_id
            >> cx >> cy >> cz
            >> v1x >> v1y >> v1z
            >> v2x >> v2y >> v2z
            >> v3x >> v3y >> v3z
            >> v4x >> v4y >> v4z
            >> v5x >> v5y >> v5z
            >> v6x >> v6y >> v6z
            >> v7x >> v7y >> v7z
            >> v8x >> v8y >> v8z
            >> fi1 >> fi2 >> fi3
            >> fi4 >> fi5 >> fi6
            >> trans_id;

      auto solid = new G4HalfSpaceArbitrary(G4ThreeVector(cx,cy,cz),
                                            G4ThreeVector(v1x, v1y, v1z),
                                            G4ThreeVector(v2x, v2y, v2z),
                                            G4ThreeVector(v3x, v3y, v3z),
                                            G4ThreeVector(v4x, v4y, v4z),
                                            G4ThreeVector(v5x, v5y, v5z),
                                            G4ThreeVector(v6x, v6y, v6z),
                                            G4ThreeVector(v7x, v7y, v7z),
                                            G4ThreeVector(v8x, v8y, v8z),
                                            fi1, fi2, fi3, fi4, fi5, fi6);

      if (trans_id > -1) {
        auto t = hs_trans_map[trans_id];
        solid->Rotate(t->GetRotationMatrix());
        solid->Translate(t->GetTranslation());
      }
      hs_surface_map[surface_id] = solid;
    }
    else if(key == "plane") {
      size_t surface_id;
      int trans_id;
      double nx, ny, nz, vx, vy, vz;
      ifstr >> surface_id >> nx >> ny >> nz >> vx >> vy >> vz >> trans_id;
      auto solid = new G4HalfSpacePlane(G4ThreeVector(nx,ny,nz), G4ThreeVector(vx,vy,vz));
      if (trans_id > -1) {
        auto t = hs_trans_map[trans_id];
        solid->Rotate(t->GetRotationMatrix());
        solid->Translate(t->GetTranslation());
      }
      hs_surface_map[surface_id] = solid;
    }
    else if(key == "xyplane") {
      size_t surface_id;
      int trans_id;
      double z;
      ifstr >> surface_id >> z >> trans_id;
      if(trans_id < 0) {
        hs_surface_map[surface_id] = new G4HalfSpaceXYPlane(z);
      }
      else {
        auto t = hs_trans_map[trans_id];
        auto solid = new G4HalfSpacePlane(G4ThreeVector(0,0,1),
                                          G4ThreeVector(0,0,z));
        solid->Rotate(t->GetRotationMatrix());
        solid->Translate(t->GetTranslation());
        hs_surface_map[surface_id] = solid;
      }
    }
    else if(key == "xzplane") {
      size_t surface_id;
      int trans_id;
      double y;
      ifstr >> surface_id >> y >> trans_id;
      if(trans_id < 0) {
        hs_surface_map[surface_id] = new G4HalfSpaceXZPlane(y);
      }
      else {
        auto t = hs_trans_map[trans_id];
        auto solid = new G4HalfSpacePlane(G4ThreeVector(0,1,0),
                                          G4ThreeVector(0,y,0));
        solid->Rotate(t->GetRotationMatrix());
        solid->Translate(t->GetTranslation());
        hs_surface_map[surface_id] = solid;
      }
    }
    else if(key == "yzplane") {
      size_t surface_id;
      int trans_id;
      double x;
      ifstr >> surface_id >> x >> trans_id;
      if(trans_id < 0) {
        hs_surface_map[surface_id] = new G4HalfSpaceYZPlane(x);
      }
      else {
        auto t = hs_trans_map[trans_id];
        auto solid = new G4HalfSpacePlane(G4ThreeVector(1,0,0),
                                          G4ThreeVector(x,0,0));
        solid->Rotate(t->GetRotationMatrix());
        solid->Translate(t->GetTranslation());
        hs_surface_map[surface_id] = solid;
      }
    }
    else if(key == "cc") {
      size_t surface_id;
      double vx, vy, vz, hx, hy, hz, r;
      ifstr >> surface_id >> vx >> vy >> vz >> hx >> hy >> hz >> r;
      hs_surface_map[surface_id] = new G4HalfSpaceCircularCylinder(G4ThreeVector(vx,vy,vz),
                                                                   G4ThreeVector(hx,hy,hz),
                                                                   r);
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
    else if(key == "ec") {
      size_t surface_id;
      double vx, vy, vz, hx, hy, hz, r1, r2;
      ifstr >> surface_id >> vx >> vy >> vz >> hx >> hy >> hz >> r1 >> r2;
      hs_surface_map[surface_id] = new G4HalfSpaceEllipticCylinder(G4ThreeVector(vx,vy,vz),
                                                                   G4ThreeVector(hx,hy,hz),r1, r2);
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
    else if(key == "quadric") {
      size_t surface_id;
      double qxx, qxy, qxz, qyy, qyz, qzz, rx, ry, rz, r;
      ifstr >> surface_id >> qxx >> qxy >> qxz >> qyy >> qyz >> qzz >> rx >> ry >> rz >> r;
      auto q =new G4HalfSpaceQuadric(qxx, qxy, qxz,
                                          qyy, qyz,
                                               qzz,
                                     rx, ry, rz,
                                     r);
      //q->Rotate(G4ThreeVector(0,0,0));
      //q->Translate(G4ThreeVector(0,0,0));
      hs_surface_map[surface_id] = q;
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
    else {
      G4cout << "G4HalfSpaceReader::Read key not found " << key << " line " << iLineNumber << G4endl;
      exit(1);
    }
  }
}
