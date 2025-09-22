#pragma once

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4HalfSpaceTest.hh"

class G4Event;
class G4ParticleGun;

/// Minimal primary generator action to demonstrate the use of GDML geometries

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    PrimaryGeneratorAction();
    PrimaryGeneratorAction(G4HalfSpaceTest *test);
    ~PrimaryGeneratorAction();

    virtual void GeneratePrimaries(G4Event* anEvent);

  private:
    G4ParticleGun* fParticleGun;
    G4HalfSpaceTest *_test;
};