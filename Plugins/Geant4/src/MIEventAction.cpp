// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// This code is based on the B1 example of Geant4

#include "Fatras/Plugins/Geant4/MIEventAction.hpp"
#include "Fatras/Plugins/Geant4/MIRunAction.hpp"

#include <fstream>

MIEventAction::MIEventAction(MIRunAction* runAction)
: G4UserEventAction(),
  fRunAction(runAction)
{
} 

void 
MIEventAction::BeginOfEventAction(const G4Event*)
{
  m_particles.clear();
}

void 
MIEventAction::AddParticle(MIparticle& p)
{
	m_particles.push_back(p);
}

std::vector<MIparticle> 
MIEventAction::particles()
{
	return m_particles;
}
