#include "G4HalfSpaceSolid.hh"
#include "G4HalfSpaceZone.hh"
#include "G4HalfSpacePlane.hh"
#include "G4ThreeVector.hh"
#include "G4FlukaReader.hh"

#include "G4BooleanSolid.hh"
#include "G4GDMLParser.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4RunManagerFactory.hh"
#include "G4TransportationManager.hh"
#include "G4Types.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4ios.hh"

#include "FTFP_BERT.hh"
#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4HalfSpaceReader.hh"

#include <vector>

int main(int argc, char** argv)
{
  std::string file_name;

  G4HalfSpaceSolid *solid = nullptr;
  G4HalfSpaceTest *test = nullptr;

  if(argc == 1 || argc > 4) {
    std::cout << "./g4halfspace filename region_number test_number" << std::endl;
    exit(1);
  }

  // First argument is file name
  if(argc >= 2) {
    file_name = argv[1];
  }

  // Halfspace file
  if( file_name.find(".inp") == std::string::npos ) {
    G4HalfSpaceReader reader(file_name);

    if (argc >= 3) {
      auto hs_region = std::stoi(argv[2]);
      solid = reader.GetSolid(hs_region);
    }
    if (argc == 4) {
      auto hs_test = std::stoi(argv[3]);
      test = reader.GetTest(hs_test);
    }
  }
  // FLUKA file
  else {
    G4FlukaReader reader(file_name);

    if (argc >= 3) {
      solid = reader.GetSolid(argv[2]);
    }
  }

  // Check solid
  if(solid == nullptr) {
    G4cout << "Solid is nullptr!" << G4endl;
    exit(1);
  }

  auto* runManager = G4RunManagerFactory::CreateRunManager();

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
