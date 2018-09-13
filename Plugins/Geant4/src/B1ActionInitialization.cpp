// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// This code is based on the B1 example of Geant4

#include "Fatras/Plugins/Geant4/B1ActionInitialization.hpp"
#include "Fatras/Plugins/Geant4/B1PrimaryGeneratorAction.hpp"
#include "Fatras/Plugins/Geant4/B1RunAction.hpp"
#include "Fatras/Plugins/Geant4/B1EventAction.hpp"
#include "Fatras/Plugins/Geant4/B1SteppingAction.hpp"

B1ActionInitialization::B1ActionInitialization(double thickness, G4ParticleGun* pGun)
 : G4VUserActionInitialization(),
 m_thickness(thickness), m_pGun(pGun)
{
	m_runAction = new B1RunAction;
	m_eventAction = new B1EventAction(m_runAction);
}

void B1ActionInitialization::Build() const
{
  SetUserAction(new B1PrimaryGeneratorAction(m_pGun));
  SetUserAction(m_runAction);
  SetUserAction(m_eventAction);
  SetUserAction(new B1SteppingAction(m_eventAction, m_thickness));
}

std::vector<B1particle> 
B1ActionInitialization::particles()
{
	return m_eventAction->particles();
}
	
