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

#include <vector>
#include <array>

class B1RunAction;

/// @brief Storage of the final state of a particle
struct B1particle
{
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
class B1EventAction : public G4UserEventAction
{
  public:
	/// @brief Constructor
	///
	/// @param [in, out] runAction Run that corresponds to the event
    B1EventAction(B1RunAction* runAction);
    
    /// @brief Destructor
    virtual ~B1EventAction() = default;

	/// @brief Initializer of an event. Resets the final state store.
    virtual void 
    BeginOfEventAction(const G4Event*);
    
    /// @brief Adds a final state particle to the store
    ///
    /// @param [in] p Particle that is added to the store
    void 
    AddParticle(B1particle& p);
	
	/// @brief Getter of the final state particles
	///
	/// @return Vector containing the final state particles
	std::vector<B1particle> 
	particles();
	
  private:
	// Vector of final state particles
	std::vector<B1particle> m_particles;
	// Pointer to the run
    B1RunAction* fRunAction;
};

    
