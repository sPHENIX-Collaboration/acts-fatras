// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// This code is based on the B1 example of Geant4

#include "Fatras/Plugins/Geant4/MISteppingAction.hpp"
#include "Fatras/Plugins/Geant4/MIEventAction.hpp"
#include "Fatras/Plugins/Geant4/MIDetectorConstruction.hpp"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"

#include "G4SystemOfUnits.hh"

MISteppingAction::MISteppingAction(MIEventAction* eventAction, double thickness)
: G4UserSteppingAction(),
  fEventAction(eventAction),
  fScoringVolume(0),
  m_thickness(thickness)
{
}

void 
MISteppingAction::UserSteppingAction(const G4Step* step)
{
  // Set scoring volume if not set yet
  if (!fScoringVolume) { 
    const MIDetectorConstruction* detectorConstruction
      = static_cast<const MIDetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    fScoringVolume = detectorConstruction->GetScoringVolume();   
  }

  // Get volume of the current step
  G4LogicalVolume* volume 
    = step->GetPreStepPoint()->GetTouchableHandle()
      ->GetVolume()->GetLogicalVolume();
      
  // Check if we are in scoring volume
  if (volume != fScoringVolume) return;
	// If the particle is at the edge of the world it is stored
	if(fabs((step->GetPostStepPoint()->GetPosition().z() - m_thickness * cm / 10.) / (m_thickness * cm / 10.)) < 1e-15
		|| step->GetPostStepPoint()->GetPosition().z() < 0.
		|| fabs((step->GetPostStepPoint()->GetPosition().x() - 1. * m) / (1. * m)) < 1e-15
		|| fabs((step->GetPostStepPoint()->GetPosition().y() - 1. * m) / (1. * m)) < 1e-15)
	{
		MIparticle p;
		p.momentum[0] = step->GetPostStepPoint()->GetMomentum().x();
		p.momentum[1] = step->GetPostStepPoint()->GetMomentum().y();
		p.momentum[2] = step->GetPostStepPoint()->GetMomentum().z();
		p.pdg = step->GetTrack()->GetDynamicParticle()->GetPDGcode();
		p.mass = step->GetPostStepPoint()->GetMass();
		p.charge = step->GetPostStepPoint()->GetCharge();
		fEventAction->AddParticle(p);
	}
}

