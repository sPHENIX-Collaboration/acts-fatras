// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_FATRAS_FATRASINTERACTOR_H
#define ACTS_FATRAS_FATRASINTERACTOR_H

#include <cmath>
#include <sstream>
#include <utility>
#include "ACTS/EventData/ParticleDefinitions.hpp"
#include "ACTS/Material/Material.hpp"
#include "ACTS/Material/MaterialProperties.hpp"
#include "ACTS/Material/SurfaceMaterial.hpp"
#include "ACTS/Surfaces/Surface.hpp"
#include "ACTS/Utilities/MaterialInteraction.hpp"
#include "RandomNumberDistributions.hpp"
#include "IMultipleScatteringSampler.hpp"
#include "IEnergyLossSampler.hpp"
#include "IHadronicInteractionSampler.hpp"

#ifndef FATRASINTERACTOR_DEBUG_OUTPUTS
#define FATRASINTERACTOR_DEBUG_OUTPUTS
#define FATLOG(cache, result, message)                                         \
  if (debug) {                                                                 \
    std::stringstream dstream;                                                 \
    dstream << "   " << std::setw(cache.debug_pfx_width);                      \
    dstream << "material interaction"                                          \
            << " | ";                                                          \
    dstream << std::setw(cache.debug_msg_width) << message << '\n';            \
    cache.debug_string += dstream.str();                                       \
  }
#endif

namespace Fatras {

  /// The information to written out per hit surface
  struct SensitiveHit
  {
  
    const Acts::Surface* surface = nullptr;
    Acts::Vector3D       position;
    Acts::Vector3D       direction;
    double               value   = 0.;
  
  };

/// The Fatras Material interactor struct
///   
/// This is the Fatras plugin to the ACTS Propagator
/// It replaces the reconstruction material updator
/// for the fast track simulation.
///
/// @tparam RandomGenerator is the given random generator type
/// @tparam SensitvitveSeclector is type of selector struct/class
///   that detrmines if a cache.current_surface is sensitive
/// @tparam PhysicsList is an extendable physics list that is called
///  
/// The physics list plays a central role in this Interactor
/// it is called on each process that is defined at compile time
/// if a process triggers an abort, this will be forwarded to
/// the propagation cache.
  
template <typename RandomGenerator, 
          typename SensitiveSelector,
          typename PhysicsList>    
struct Interactor
{

  /// The random generator to be spawn per event
  RandomGenerator* randomGenerator = nullptr;

  /// The selector for sensitive surfaces
  SensitiveSelector sensitiveSelector;
  
  /// debug output flag
  bool debug = false;
  
  /// The particle masses for lookup
  Acts::ParticleMasses particleMasses;

  /// Simple result struct to be returned
  /// It mainly acts as an interal state cache which is
  /// created for every propagation/extrapolation step
  struct this_result
  {
    // the sensitive hits
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
  /// @tparam cache_t is the type of Stepper cache
  ///
  /// @param cache is the mutable stepper cache object
  /// @param result is the mutable result cache object
  ///
  /// return value is void as it is a standard actor in the 
  /// propagation
  template <typename cache_t>
  void
  operator()(cache_t& cache, result_type& result) const
  {

    // If we are on target, everything should have been done
    if (cache.target_reached) return;

    // Check if the current surrface a senstive one 
    bool isSensitive = cache.current_surface ?
      sensitiveSelector(*cache.current_surface) : false;
    double depositedEnergy = 0.;

    // a current surface has been assigned by the navigator
    if (cache.current_surface && cache.current_surface->associatedMaterial()) {
      // @todo adapt sign & particle type
      Acts::ParticleType pType = Acts::muon;
      // get the surface material and the corresponding material properties
      auto sMaterial   = cache.current_surface->associatedMaterial();
      auto mProperties = sMaterial->material(cache.position());
      if (mProperties) {
        // The: pre - full - post update test
        // check if you have a factor for pre/post/full update to do
        double prepofu = 1.;
        if (cache.start_surface == cache.current_surface) {
          FATLOG(cache, result, "Update on start surface: post-update mode.");
          prepofu = cache.current_surface->associatedMaterial()->factor(
              cache.nav_dir, Acts::postUpdate);
        } else if (cache.target_surface == cache.current_surface) {
          FATLOG(cache, result, "Update on target surface: pre-update mode.");
          prepofu = cache.current_surface->associatedMaterial()->factor(
              cache.nav_dir, Acts::preUpdate);
        } else
          FATLOG(cache, result, "Update while pass through: full mode.");
        
        if (prepofu == 0.) {
          FATLOG(cache, result, "Pre/Post factor set material to zero.");
        } 
        
        
        
      }
    }
    
    // create the SensitiveHit and store it
    if (isSensitive){
      // create and fill the hit 
      SensitiveHit senHit;
      senHit.surface   = cache.current_surface;
      senHit.position  = cache.position();
      senHit.direction = cache.direction();
      senHit.value     = depositedEnergy;
      result.sensitiveHits.push_back(std::move(senHit));
    }
  }

  /// Pure observer interface
  /// This does not apply to the surface collector
  template <typename cache_t>
  void
  operator()(cache_t& cache) const
  {
    (void)cache;
  }
};
}

#endif
