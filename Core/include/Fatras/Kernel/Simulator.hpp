// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_FATRAS_SIMULATOR_H
#define ACTS_FATRAS_SIMULATOR_H

#include <cmath>
#include <climits>
#include <sstream>
#include <utility>
#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialProperties.hpp"
#include "Acts/Material/SurfaceMaterial.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/MaterialInteraction.hpp"
#include "Fatras/Kernel/Particle.hpp"
#include "RandomNumberDistributions.hpp"

#ifndef FATRASSIMULATOR_DEBUG_OUTPUTS
#define FATRASSIMULATOR_DEBUG_OUTPUTS
#define FATLOG(cache, result, message)                              \
  if (debug) {                                                      \
    std::stringstream dstream;                                      \
    dstream << "   " << std::setw(cache.debugPfxWidth);           \
    dstream << "fatras kernel"                                      \
            << " | ";                                               \
    dstream << std::setw(cache.debugMsgWidth) << message << '\n'; \
    cache.debugString += dstream.str();                            \
  }
#endif

namespace Fatras {

  /// Using the Definitions for Detector and Particle
  typedef DetectorInfo Detector;

  /// The Fatras Simulator
  ///   
  /// This is the Fatras plugin to the ACTS Propagator, it replaces
  /// the MaterialInteractor of the reconstruction 
  ///
  /// @tparam RandomGenerator is the given random generator type
  /// @tparam SensitiveSelector is type of selector struct/class
  ///   that detrmines if a state.navigation.currentSurface is sensitive
  /// @tparam PhysicsList is an extendable physics list that is called
  ///  
  /// The physics list plays a central role in this DetectorInteractor
  /// it is called on each process that is defined at compile time
  /// if a process triggers an abort, this will be forwarded to
  /// the propagation cache.
    
  template <typename RandomGenerator, 
            typename SensitiveSelector,
            typename PhysicsList>    
  struct FatrasKernel
  {
  
    /// The random generator to be spawn per event
    RandomGenerator* generator = nullptr;
  
    /// The selector for sensitive surfaces
    SensitiveSelector sensitiveSelector;
    
    /// The physics list provided for this call
    PhysicsList       physicsList;
        
    /// debug output flag
    bool debug                = false;
    
    /// Simple result struct to be returned
    /// 
    /// It mainly acts as an internal state cache which is
    /// created for every propagation/extrapolation step
    struct this_result
    {
      /// The current particle - can be updated of course 
      Particle          particle;
      /// The outgoing particles due to physics processes
      std::vector<Partilce>     outgoing; 
      /// The sensitive hits created along the way
      std::vector<SensitiveHit> sensitiveHits;      
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
    /// @param state is the mutable stepper state object
    /// @param result is the mutable result cache object
    ///
    /// return value is void as it is a standard actor in the 
    /// propagation
    template <typename propagator_state_t>
    void
    operator()(propagator_state_t& state, 
               result_type& result) const
    {
  
      // If we are on target, everything should have been done
      if (state.navigation.targetReached) return;
  
      // Check if the current surrface a senstive one 
      bool isSensitive = state.navigation.currentSurface ?
        sensitiveSelector(*state.navigation.currentSurface) : false;
      double depositedEnergy = 0.;
  
      // a current surface has been assigned by the navigator
      if (state.navigation.currentSurface && state.navigation.currentSurface->associatedMaterial()) {
        // get the surface material and the corresponding material properties
        auto sMaterial   = state.navigation.currentSurface->associatedMaterial();
        auto mProperties = sMaterial->material(state.stepping.position());
        if (mProperties) {
           // run the Fatras physics list
           auto breakIndicator = 
             physicsList(*generator,*mProperties,result.particle,result.outgoing);
         }
      }
      // update the stepper cache with the current particle parameters
      state.stepping.position  = particle.position();
      state.stepping.direction = particle.momentum().unit();
      state.stepping.qop       = particle.charge()/particle.momentum().mag();
      
      // create the SensitiveHit and store it
      if (isSensitive){
        // create and fill the hit 
        SensitiveHit senHit;
        senHit.surface   = state.navigation.currentSurface;
        senHit.position  = state.stepping.position();
        senHit.direction = state.stepping.direction();
        senHit.value     = depositedEnergy;
        result.sensitiveHits.push_back(std::move(senHit));
      }
    }
  
    /// Pure observer interface
    /// This does not apply to the Fatras simulator
    template <typename propagator_state_t>
    void
    operator()(propagator_state_t&) const
    {}
      
  };
}

#endif // ACTS_FATRAS_SIMULATOR_H
