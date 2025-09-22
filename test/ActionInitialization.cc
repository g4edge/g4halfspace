#include "ActionInitialization.hh"

#include "PrimaryGeneratorAction.hh"

#include "G4HalfSpaceTest.hh"

ActionInitialization::ActionInitialization() :
  G4VUserActionInitialization(),
  _test(nullptr)
{}

ActionInitialization::ActionInitialization(G4HalfSpaceTest *test) :
  G4VUserActionInitialization(),
  _test(test)
{}

ActionInitialization::~ActionInitialization() {}

void ActionInitialization::BuildForMaster() const {}

void ActionInitialization::Build() const
{
  if(_test == nullptr) {
    SetUserAction(new PrimaryGeneratorAction());
  }
  else {
    SetUserAction(new PrimaryGeneratorAction(_test));
  }
}
