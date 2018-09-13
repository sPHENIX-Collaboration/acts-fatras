// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// This code is based on the B1 example of Geant4

#include "Fatras/Plugins/Geant4/B1RunAction.hpp"

#include "G4RunManager.hh"
#include "G4Run.hh"

B1RunAction::B1RunAction()
: G4UserRunAction()
{
}

void B1RunAction::BeginOfRunAction(const G4Run*)
{ 
  // Inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
}

void B1RunAction::EndOfRunAction(const G4Run* run)
{
  if(run->GetNumberOfEvent() == 0) 
	return;
}
