// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// This code is based on the B1 example of Geant4

#pragma once

#include "G4ParticleGun.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

class G4ParticleGun;
class G4Event;

/// @brief The primary generator action class with particle gun.
class MIPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
  /// @brief Constructor
  ///
  /// @param [in] pGun Particle gun
  MIPrimaryGeneratorAction(G4ParticleGun *pGun);

  /// @brief Destructor
  virtual ~MIPrimaryGeneratorAction();

  /// @brief Method from the base class. Fires the gun.
  ///
  /// @param [out] anEvent Event related to the gun
  virtual void GeneratePrimaries(G4Event *anEvent);

  /// @brief Method to access particle gun.
  ///
  /// @return Pointer to the particle gun
  const G4ParticleGun *GetParticleGun() const;

private:
  // Pointer to the Geant4 particle gun
  G4ParticleGun *fParticleGun;
};
