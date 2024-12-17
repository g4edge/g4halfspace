#pragma once

#include "G4VUserDetectorConstruction.hh"

class G4HalfSpaceSolid;

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  DetectorConstruction(G4HalfSpaceSolid *hss);

  virtual G4VPhysicalVolume* Construct();

private:
  G4HalfSpaceSolid *_hss;
  G4VPhysicalVolume* fWorld;
};