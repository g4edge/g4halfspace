#include "PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include "G4HalfSpaceTest.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction()
  : G4VUserPrimaryGeneratorAction(), fParticleGun(0)
{
  G4int n_particle = 1;
  fParticleGun = new G4ParticleGun(n_particle);

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  fParticleGun->SetParticleDefinition(particleTable->FindParticle(particleName = "e-"));
  //fParticleGun->SetParticleDefinition(particleTable->FindParticle(particleName = "geantino"));
  fParticleGun->SetParticleEnergy(1 * GeV);
}

PrimaryGeneratorAction::PrimaryGeneratorAction(G4HalfSpaceTest *test)
  : G4VUserPrimaryGeneratorAction(), fParticleGun(0), _test(test)
{
  G4int n_particle = 1;
  fParticleGun = new G4ParticleGun(n_particle);

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  fParticleGun->SetParticleDefinition(particleTable->FindParticle(particleName = "geantino"));
  fParticleGun->SetParticleEnergy(1 * GeV);
}


PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    if(_test == nullptr) {
      G4ThreeVector v(0, 0, 1);
      fParticleGun->SetParticlePosition(
        G4ThreeVector(900 * (G4UniformRand() - 0.5) * mm,
                      900 * (G4UniformRand() - 0.5) * mm,
                      -900 * mm));
      v = v / v.mag();
      fParticleGun->SetParticleMomentumDirection(v);
    }
    else {
      _test->current_test()++;
      G4ThreeVector dir = _test->direction();
      G4ThreeVector pos = _test->position();

      std::cout << "PrimaryGeneratorAction::GeneratePrimaries> test "
                << pos << " " << dir << std::endl;
      fParticleGun->SetParticlePosition(pos);
      fParticleGun->SetParticleMomentumDirection(dir);
    }
  fParticleGun->GeneratePrimaryVertex(anEvent);
}
