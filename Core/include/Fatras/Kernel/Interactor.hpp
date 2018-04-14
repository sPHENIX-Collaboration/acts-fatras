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
/// - this is the Fatras plugin to the Propagator and will make 
///
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

  /// Interaction with detector material
  /// for the ActionList of the Propagator
  ///
  /// It checks if the cache has a current surface,
  /// in which case the action is performed according to
  /// the configruation, eventual particles produced
  /// in electromagnetic or hadronic interactions
  /// are stored in the result vector 
  ///
  /// @tparam cache_t is the type of Stepper cache
  ///
  /// @param cache is the mutable stepper cache object
  /// @param result is the mutable result cache object
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
        } else {

          // get the path correction due to the incident angle
          double pCorrection = cache.current_surface->pathCorrection(
              cache.position(), cache.direction());
        
          // the momentum at current position
          double p    = std::abs(1. / cache.qop);
          double m    = particleMasses.mass.at(pType);
          double E    = std::sqrt(p * p + m * m);
          double beta = p / E;
        
          // apply the multiple scattering
          if (multipleScatteringSampler && randomGenerator) {
            // the multiple scattering sampler 
            double sTheta = multipleScatteringSampler->simTheta(*randomGenerator,
                                                                *mProperties,
                                                                p, pCorrection,
                                                                pType);
          
            // Create a random uniform distribution between in the intervall [0,1]
            UniformDist uniformDist(0., 1.);
          
            //@todo test this non parametric way - not tested yet
            double psi     = 2. * M_PI * uniformDist(*randomGenerator);
          
            // more complex but "more true"
            Acts::Vector3D pDirection(cache.momentum().unit());
            double         x = -pDirection.y();
            double         y =  pDirection.x();
            double         z = 0.;
            // if it runs along the z axis - no good ==> take the x axis
            if (pDirection.z() * pDirection.z() > 0.999999) { x = 1.;  y = 0.; }
          
            // deflector direction
            Acts::Vector3D deflector(x, y, z);
            // rotate the new direction for scattering using theta and arbitraril in psi
            // create the rotation
            Acts::RotationMatrix3D rotation;
            rotation = Acts::AngleAxis3D(sTheta, deflector)
                * Acts::AngleAxis3D(psi, pDirection);
            // Some screen output
            FATLOG(cache, result, "Multiple scattering deltaPsi / deltaTheta = " 
                                  << psi << " / " << sTheta);         
            // create the transform
            // @todo ? can we not just use rotation ??
            Acts::Transform3D transform(rotation, Acts::Vector3D(0., 0., 0.));
            // get the new direction
            pDirection = transform.linear()*pDirection;
            // assign the new direction - make sure it's the unit vector
            cache.dir = pDirection.unit();
          }
          // apply the energy loss
          if (energyLossSampler && randomGenerator) {
            // energy loss from the eloss sampler
            auto eloss = energyLossSampler->energyLoss(*randomGenerator,
                                                       *mProperties,
                                                       p, pCorrection,
                                                       cache.nav_dir,
                                                       pType);
            depositedEnergy = eloss.deltaE();                                           
            FATLOG(cache, result, "Energy loss sampled " << eloss.deltaE());
            if ((E + eloss.deltaE()) < m) return;
            // indeed possible, apply and update
            double newP = sqrt((E + eloss.deltaE()) * (E + eloss.deltaE()) - m * m);
            double q = cache.qop > 0. ? 1. : -1.;
            cache.qop = q/newP;             
          }
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
