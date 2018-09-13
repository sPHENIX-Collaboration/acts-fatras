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
#include "G4VUserActionInitialization.hh"
#include "MIEventAction.hpp"
#include "MIRunAction.hpp"

/// @brief Action initialization class.
class MIActionInitialization : public G4VUserActionInitialization {
public:
  /// @brief Constructor
  ///
  /// @param [in] thickness Thickness of the material
  /// @param [in] pGun Particle Gun
  MIActionInitialization(double thickness, G4ParticleGun *pGun);

  /// @brief Destructor
  virtual ~MIActionInitialization();

  /// @brief Sets up all necessacry components
  virtual void Build() const;

  /// @brief Getter of the final state particles
  ///
  /// @return Vector containing the final state particles
  std::vector<MIparticle> particles();

protected:
  // Thickness of the material
  double m_thickness;
  // Pointer to the Particle gun
  G4ParticleGun *m_pGun;
  // Pointer to the Run
  MIRunAction *m_runAction;
  // Pointer to the Event
  MIEventAction *m_eventAction;
};
