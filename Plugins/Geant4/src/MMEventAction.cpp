// This file is part of the ACTS project.
//
// Copyright (C) 2017 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Fatras/Plugins/Geant4/MMEventAction.hpp"
#include <stdexcept>
#include "Fatras/Plugins/Geant4/MMPrimaryGeneratorAction.hpp"
#include "Fatras/Plugins/Geant4/MMSteppingAction.hpp"
#include "G4Event.hh"
#include "G4RunManager.hh"

FW::G4::MMEventAction* FW::G4::MMEventAction::fgInstance = nullptr;

FW::G4::MMEventAction*
FW::G4::MMEventAction::Instance()
{
  // Static acces function via G4RunManager
  return fgInstance;
}

FW::G4::MMEventAction::MMEventAction() : G4UserEventAction()
{
  if (fgInstance) {
    throw std::logic_error("Attempted to duplicate a singleton");
  } else {
    fgInstance = this;
  }
}

FW::G4::MMEventAction::~MMEventAction()
{
  fgInstance = nullptr;
}

void
FW::G4::MMEventAction::BeginOfEventAction(const G4Event*)
{
  // reset the collection of material steps
  MMSteppingAction::Instance()->Reset();
}

void
FW::G4::MMEventAction::EndOfEventAction(const G4Event* event)
{
  Acts::MaterialStep::Position pos(event->GetPrimaryVertex()->GetX0(),
                                   event->GetPrimaryVertex()->GetY0(),
                                   event->GetPrimaryVertex()->GetZ0());
  // access the initial direction of the track
  G4ThreeVector dir   = MMPrimaryGeneratorAction::Instance()->direction();
  double        theta = dir.theta();
  double        phi   = dir.phi();
  // loop over the material steps and add up the material
  double tX0 = 0;
  double tL0 = 0;
  for (auto& mstep : MMSteppingAction::Instance()->materialSteps()) {
    tX0 += mstep.materialProperties().thicknessInX0();
    tL0 += mstep.materialProperties().thicknessInL0();
  }
  // create the MaterialTrack
  Acts::MaterialTrack mtrecord(
      pos, theta, phi, MMSteppingAction::Instance()->materialSteps(), tX0, tL0);
  // write out the MaterialTrack of one event
  m_records.push_back(mtrecord);
}

void
FW::G4::MMEventAction::Reset()
{
}
