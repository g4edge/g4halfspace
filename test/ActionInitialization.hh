#pragma once

#include "G4VUserActionInitialization.hh"

/// Action initialization class.

class G4HalfSpaceTest;

class ActionInitialization : public G4VUserActionInitialization
{
public:
  ActionInitialization();
  ActionInitialization(G4HalfSpaceTest *test);
  virtual ~ActionInitialization();

  virtual void BuildForMaster() const;
  virtual void Build() const;

protected:
  G4HalfSpaceTest *_test;
};

