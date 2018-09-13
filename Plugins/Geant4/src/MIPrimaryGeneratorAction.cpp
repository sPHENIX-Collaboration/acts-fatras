// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// This code is based on the B1 example of Geant4

#include "Fatras/Plugins/Geant4/MIPrimaryGeneratorAction.hpp"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

MIPrimaryGeneratorAction::MIPrimaryGeneratorAction(G4ParticleGun *pGun)
    : G4VUserPrimaryGeneratorAction() {
  fParticleGun = pGun;
}

MIPrimaryGeneratorAction::~MIPrimaryGeneratorAction() { delete fParticleGun; }

void MIPrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent) {
  // This function is called at the begining of each event
  fParticleGun->GeneratePrimaryVertex(anEvent);
}

const G4ParticleGun *MIPrimaryGeneratorAction::GetParticleGun() const {
  return fParticleGun;
}
