#include "G4HalfSpaceSolid.hh"
#include "G4HalfSpaceZone.hh"
#include "G4HalfSpacePlane.hh"
#include "G4ThreeVector.hh"
#include "G4ios.hh"

#include "G4BooleanSolid.hh"
#include "G4GDMLParser.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4RunManagerFactory.hh"
#include "G4TransportationManager.hh"
#include "G4Types.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"

#include "FTFP_BERT.hh"
#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4HalfSpaceReader.hh"

#include <vector>

int main(int argc, char** argv)
{
  auto hsr = G4HalfSpaceReader(argv[1]);
  auto hs_region = std::stoi(argv[2]);
  auto hs_test = std::stoi(argv[3]);

  auto* runManager = G4RunManagerFactory::CreateRunManager();

  auto solid = hsr.GetSolid(hs_region);
  auto test = hsr.GetTest(hs_test);
  solid->AddTestInstrument(test);

  runManager->SetUserInitialization(new DetectorConstruction(solid));
  runManager->SetUserInitialization(new ActionInitialization(test));
  runManager->SetUserInitialization(new FTFP_BERT);

  runManager->Initialize();

  // Initialize visualization
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  runManager->BeamOn(0);

  if (argc == 5)  // batch mode
  {
      G4String command = "/control/execute ";
      G4String fileName = argv[4];
      UImanager->ApplyCommand(command + fileName);
  }
  else  // interactive mode
  {
      G4UIExecutive* ui = new G4UIExecutive(argc, argv);
      UImanager->ApplyCommand("/control/execute vis.mac");
      ui->SessionStart();
      delete ui;
  }

  if(test != nullptr)
    test->print_data();

  delete visManager;
  delete runManager;

  return 0;
}
