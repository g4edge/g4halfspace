#include "DetectorConstruction.hh"

#include "G4NistManager.hh"
#include "G4Box.hh"

#include "G4HalfSpaceSolid.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4BoxInstrumented.hh"


DetectorConstruction::DetectorConstruction(G4HalfSpaceSolid *hss) {
  _hss = hss;
}

G4VPhysicalVolume* DetectorConstruction::Construct() {

  std::cout << "DetectorConstruction::Construct" << std::endl;

  //////////////////////////////
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  //////////////////////////////
  // Envelope parameters
  G4double env_sizeXY = 20*cm, env_sizeZ = 30*cm;
  G4Material* env_mat = nist->FindOrBuildMaterial("G4_Al");

  //////////////////////////////
  // World
  G4double world_sizeXY = env_sizeXY;
  G4double world_sizeZ  = env_sizeZ;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_Galactic");

  auto solidWorld = new G4Box("World",
                              0.5 * world_sizeXY,
                              0.5 * world_sizeXY,
                              0.5 * world_sizeZ);

  auto logicWorld = new G4LogicalVolume(solidWorld,
                                        world_mat,
                                        "World");

  auto physWorld = new G4PVPlacement(nullptr,
                                     G4ThreeVector(),
                                     logicWorld,
                                     "World",
                                     nullptr,
                                     false,
                                     0,
                                     false);

  std::cout << "DetectorConstruction::Construct" << _hss << std::endl;

  auto logicHSS = new G4LogicalVolume(_hss,
                                      env_mat,
                                      "hss");

  auto physHSS = new G4PVPlacement(nullptr,
                                   G4ThreeVector(),
                                   logicHSS,
                                   "physHSS",
                                    logicWorld,
                                    //nullptr,
                                   false,
                                   0,
                                   false);

  return physWorld;
}



#if 0
auto p1 = new G4HalfSpacePlane(G4ThreeVector(25*mm,0,0),G4ThreeVector(1,0,0));
    auto p2 = new G4HalfSpacePlane(G4ThreeVector(-25*mm,0,0),G4ThreeVector(-1,0,0));
    auto p3 = new G4HalfSpacePlane(G4ThreeVector(0,25*mm,0),G4ThreeVector(0,1,0));
    auto p4 = new G4HalfSpacePlane(G4ThreeVector(0,-25*mm,0),G4ThreeVector(0,-1,0));
    auto p5 = new G4HalfSpacePlane(G4ThreeVector(0,0,25*mm),G4ThreeVector(0,0,1));
    auto p6 = new G4HalfSpacePlane(G4ThreeVector(0,0,-25*mm),G4ThreeVector(0,0,-1));
    auto p7 = new G4HalfSpacePlane(G4ThreeVector(-12.5*mm,-12.5*mm,-12.5*mm), G4ThreeVector(1,1,0));
    auto p8 = new G4HalfSpacePlane(G4ThreeVector(0,-12.5*mm,-12.5*mm), G4ThreeVector(1,0,1));
    auto z = new G4HalfSpaceZone();
    z->AddIntersection(p1);
    z->AddIntersection(p2);
    z->AddIntersection(p3);
    z->AddIntersection(p4);
    z->AddIntersection(p5);
    z->AddIntersection(p6);
    auto hss = new G4HalfSpaceSolid("hsSolid");
    hss->AddZone(z);
#endif

#if 0
auto b1 = new G4HalfSpaceAARBox(-25*mm, 25*mm,-25*mm, 25*mm, -25*mm, 25*mm);
    auto b2 = new G4HalfSpaceAARBox(-0*mm, 50*mm,0*m, 50*mm, 0*mm, 50*mm);
    auto z = new G4HalfSpaceZone();
    z->AddIntersection(b1);
    z->AddSubtraction(b2);
    auto hss = new G4HalfSpaceSolid("hsSolid");
    hss->AddZone(z);
    hss->Rotate(G4RotationMatrix(3.14159/8,0,0));
    //hss->Translate(G4ThreeVector(30,30,30))
#endif

#if 0
auto s1 = new G4HalfSpaceSphere(50*mm, G4ThreeVector(0,0,0));
    auto s2 = new G4HalfSpaceSphere(50*mm, G4ThreeVector(25*mm,25*mm,25*mm));
    auto b2 = new G4HalfSpaceAARBox(0*mm, 50*mm,0*m, 50*mm, 0*mm, 50*mm);
    auto z = new G4HalfSpaceZone();
    z->AddIntersection(s1);
    z->AddSubtraction(b2);
    auto z2 = new G4HalfSpaceZone();
    z2->AddIntersection(b2);
    auto hss = new G4HalfSpaceSolid("hsSolid");
    hss->addZone(z);
#endif

#if 0
auto c1 = new G4HalfSpaceXACircularCylinder(0,0,25*mm);
    auto p1 = new G4HalfSpacePlane(G4ThreeVector(25*mm,0,0),
                                   G4ThreeVector(1,0,0));
    auto p2 = new G4HalfSpacePlane(G4ThreeVector(-25*mm,0,0),
                                   G4ThreeVector(-1,0,0));
    auto z = new G4HalfSpaceZone();
    z->AddIntersection(c1);
    z->AddIntersection(p1);
    z->AddIntersection(p2);
    auto hss = new G4HalfSpaceSolid("hsSolid");
    hss->AddZone(z);
#endif

#if 0
auto c1 = new G4HalfSpaceYACircularCylinder(0,0,25*mm);
    auto p1 = new G4HalfSpacePlane(G4ThreeVector(0, 25*mm,0),
                                   G4ThreeVector(0,1,0));
    auto p2 = new G4HalfSpacePlane(G4ThreeVector(0, -25*mm,0),
                                   G4ThreeVector(0,-1,0));
    auto z = new G4HalfSpaceZone();
    z->AddIntersection(c1);
    z->AddIntersection(p1);
    z->AddIntersection(p2);
    auto hss = new G4HalfSpaceSolid("hsSolid");
    hss->AddZone(z);
#endif

#if 0
auto c1 = new G4HalfSpaceZACircularCylinder(0,0,25*mm);
    auto p1 = new G4HalfSpacePlane(G4ThreeVector(0, 0, 25*mm),
                                   G4ThreeVector(0,0, 1));
    auto p2 = new G4HalfSpacePlane(G4ThreeVector(0, 0,-25*mm),
                                   G4ThreeVector(0,0, -1));
    auto z = new G4HalfSpaceZone();
    z->AddIntersection(c1);
    z->AddIntersection(p1);
    z->AddIntersection(p2);
    auto hss = new G4HalfSpaceSolid("hsSolid");
    hss->AddZone(z);
#endif

#if 0
auto e1 = new G4HalfSpaceXAEllipticalCylinder(0,0,25*mm,50*mm);
    auto p1 = new G4HalfSpacePlane(G4ThreeVector(0, 0, 25*mm),
                                   G4ThreeVector(0,0, 1));
    auto p2 = new G4HalfSpacePlane(G4ThreeVector(0, 0,-25*mm),
                                   G4ThreeVector(0,0, -1));
    auto z = new G4HalfSpaceZone();
    z->AddIntersection(e1);
    z->AddIntersection(p1);
    z->AddIntersection(p2);
    auto hss = new G4HalfSpaceSolid("hsSolid");
    hss->AddZone(z);
#endif