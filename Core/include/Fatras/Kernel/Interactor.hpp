// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialProperties.hpp"
#include "Acts/Material/SurfaceMaterial.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Propagator/ActionList.hpp"
#include "Fatras/Kernel/Particle.hpp"
#include "Fatras/Kernel/Definitions.hpp"
#include "Fatras/Kernel/PhysicsList.hpp"
#include "RandomNumberDistributions.hpp"
#include <climits>
#include <cmath>
#include <sstream>
#include <utility>

namespace Fatras {

struct VoidSelector {
  
  bool operator()(const Acts::Surface& ) const { return false;}

};

/// The Fatras Simulator
///
/// This is the Fatras plugin to the ACTS Propagator, it replaces
/// the MaterialInteractor of the reconstruction
///
/// @tparam generator_t Type of the random generator
/// @tparam physics_list_t Type of Extendable physics list that is called
/// @tparam decay_list_t Type of Extendable decay list that is called
///
/// The physics list plays a central role in this DetectorInteractor
/// it is called on each process that is defined at compile time
/// if a process triggers an abort, this will be forwarded to
/// the propagation cache.
template <typename generator_t,
          typename sensitive_selector_t = VoidSelector,
          typename physics_list_t = PhysicsList<>>
struct Interactor {

  /// The random generator to be spawnper event
  generator_t *generator = nullptr;

  /// The slector for sensitive surfaces
  sensitive_selector_t sensitiveSelector;

  /// The physics list provided for this call
  physics_list_t physicsList;

  /// Simple result struct to be returned
  Particle initialParticle;
    
  ///
  /// It mainly acts as an internal state cache which is
  /// created for every propagation/extrapolation step
  struct this_result {

    /// result initialization
    bool initialized = false;

    /// The current particle - updated along the way
    Particle particle;
  
    /// The outgoing particles due to physics processes
    std::vector<Particle> outgoing;
    
    /// The simulated hits created along the way
    std::vector<SensitiveHit> simulatedHits;

  };

  typedef this_result result_type;

  /// Interaction with detector material for the ActionList of the Propagator
  ///
  /// It checks if the cache has a current surface, in which case the action
  /// is performed according to the physics list content.
  ///
  /// Eventual particles produced in electromagnetic or hadronic interactions
  /// are stored in the result struct and can thus be retrieved by the caller
  ///
  /// @tparam  propagator_state_t is the type of Propagtor state
  ///
  /// @param state is the mutable propagator state object
  /// @param result is the mutable result cache object
  ///
  /// return value is void as it is a standard actor in the
  /// propagation
  template <typename propagator_state_t>
  void operator()(propagator_state_t &state, 
                  result_type &result) const {
                    
    // If we are on target, everything should have been done
    if (state.navigation.targetReached)
      return;
    
    // Initialize the result, the state is thread local
    if (!result.initialized){
      // set the initial particle parameters
      result.particle    = initialParticle;
      result.initialized = true;
    }
    // set the stepping position to the particle
    result.particle.position = state.stepping.position();
    result.particle.momentum = state.stepping.momentum();
    result.particle.q        = state.stepping.charge();
      
    // Check if the current surrface a senstive one
    bool isSensitive = state.navigation.currentSurface
                           ? sensitiveSelector(*state.navigation.currentSurface)
                           : false;
    double depositedEnergy = 0.;

    // a current surface has been assigned by the navigator
    if (state.navigation.currentSurface &&
        state.navigation.currentSurface->associatedMaterial()) {
      // get the surface material and the corresponding material properties
      auto sMaterial = state.navigation.currentSurface->associatedMaterial();
      auto mProperties = sMaterial->material(state.stepping.position());
      
      bool breakIndicator = false;
      if (mProperties) {
        // run the Fatras physics list - only when there's material 
        breakIndicator = physicsList(*generator, 
                                     *mProperties,
                                     result.particle, 
                                     result.outgoing);
      }
    }

    // update the stepper cache with the current particle parameters
     state.stepping.update(result.particle.position,
                           result.particle.momentum.unit(),
                           result.particle.momentum.mag());
     
    // create the SensitiveHit and store it
    if (isSensitive) {
      // create and fill the hit
      SensitiveHit senHit;
      senHit.surface   = state.navigation.currentSurface;
      senHit.position  = state.stepping.position();
      senHit.direction = state.stepping.direction();
      senHit.value     = depositedEnergy;
      result.simulatedHits.push_back(std::move(senHit));
    }
  }

  /// Pure observer interface
  /// This does not apply to the Fatras simulator
  template <typename propagator_state_t>
  void operator()(propagator_state_t &) const {}
};
}
