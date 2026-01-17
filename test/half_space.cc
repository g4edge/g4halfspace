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
  std::string file_name;
  size_t hs_region = -1;
  size_t hs_test = -1;

  if(argc == 1 || argc > 4) {
    std::cout << "./g4halfspace filename region_number test_number" << std::endl;
    exit(1);
  }

  if(argc >= 2) {
    file_name = argv[1];
  }
  if(argc >= 3) {
    hs_region = std::stoi(argv[2]);
  }
  if(argc == 4) {
    hs_test = std::stoi(argv[3]);
  }

  auto hsr = G4HalfSpaceReader(file_name);

  auto* runManager = G4RunManagerFactory::CreateRunManager();

  auto solid = hsr.GetSolid(hs_region);
  if(!solid) {
    std::cout << "Not a valid solid number" << std::endl;
    exit(1);
  }

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
