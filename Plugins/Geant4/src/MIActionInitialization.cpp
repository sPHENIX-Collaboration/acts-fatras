// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// This code is based on the B1 example of Geant4

#include "Fatras/Plugins/Geant4/MIActionInitialization.hpp"
#include "Fatras/Plugins/Geant4/MIEventAction.hpp"
#include "Fatras/Plugins/Geant4/MIPrimaryGeneratorAction.hpp"
#include "Fatras/Plugins/Geant4/MIRunAction.hpp"
#include "Fatras/Plugins/Geant4/MISteppingAction.hpp"

MIActionInitialization::MIActionInitialization(double thickness,
                                               G4ParticleGun *pGun)
    : G4VUserActionInitialization(), m_thickness(thickness), m_pGun(pGun) {
  m_runAction = new MIRunAction;
  m_eventAction = new MIEventAction(m_runAction);
}

void MIActionInitialization::Build() const {
  SetUserAction(new MIPrimaryGeneratorAction(m_pGun));
  SetUserAction(m_runAction);
  SetUserAction(m_eventAction);
  SetUserAction(new MISteppingAction(m_eventAction, m_thickness));
}

std::vector<MIparticle> MIActionInitialization::particles() {
  return m_eventAction->particles();
}

MIActionInitialization::~MIActionInitialization() { delete (m_pGun); }
