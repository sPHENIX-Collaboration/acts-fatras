// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// This code is based on the B1 example of Geant4

#pragma once

#include "Acts/Utilities/Definitions.hpp"
#include "G4UserEventAction.hh"

#include <array>
#include <vector>

class MIRunAction;

/// @brief Storage of the final state of a particle
struct MIparticle {
  // Momentum of a particle
  Acts::Vector3D momentum;
  // PDG code of a particle
  int pdg;
  // Charge of a particle
  int charge;
  // Mass of a particle
  double mass;
};

/// @brief Event action class. Stores the final state particles.
class MIEventAction : public G4UserEventAction {
public:
  /// @brief Constructor
  ///
  /// @param [in, out] runAction Run that corresponds to the event
  MIEventAction(MIRunAction *runAction);

  /// @brief Destructor
  virtual ~MIEventAction() = default;

  /// @brief Initializer of an event. Resets the final state store.
  virtual void BeginOfEventAction(const G4Event *);

  /// @brief Adds a final state particle to the store
  ///
  /// @param [in] p Particle that is added to the store
  void AddParticle(MIparticle &p);

  /// @brief Getter of the final state particles
  ///
  /// @return Vector containing the final state particles
  std::vector<MIparticle> particles();

private:
  // Vector of final state particles
  std::vector<MIparticle> m_particles;
  // Pointer to the run
  MIRunAction *fRunAction;
};
