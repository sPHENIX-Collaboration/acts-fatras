// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// This code is based on the B1 example of Geant4

#include "Fatras/Plugins/Geant4/B1PrimaryGeneratorAction.hpp"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

B1PrimaryGeneratorAction::B1PrimaryGeneratorAction(G4ParticleGun* pGun)
: G4VUserPrimaryGeneratorAction()
{
  fParticleGun  = pGun;
}

B1PrimaryGeneratorAction::~B1PrimaryGeneratorAction()
{
  delete fParticleGun;
}

void 
B1PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // This function is called at the begining of each event
  fParticleGun->GeneratePrimaryVertex(anEvent);
}

const G4ParticleGun* 
B1PrimaryGeneratorAction::GetParticleGun() const 
{ 
	return fParticleGun; 
}
