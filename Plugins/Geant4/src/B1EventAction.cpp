// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// This code is based on the B1 example of Geant4

#include "Fatras/Plugins/Geant4/B1EventAction.hpp"
#include "Fatras/Plugins/Geant4/B1RunAction.hpp"

#include <fstream>

B1EventAction::B1EventAction(B1RunAction* runAction)
: G4UserEventAction(),
  fRunAction(runAction)
{
} 

void 
B1EventAction::BeginOfEventAction(const G4Event*)
{
  m_particles.clear();
}

void 
B1EventAction::AddParticle(B1particle& p)
{
	m_particles.push_back(p);
}

std::vector<B1particle> 
B1EventAction::particles()
{
	return m_particles;
}
