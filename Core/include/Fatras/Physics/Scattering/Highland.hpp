// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef FATRAS_PHYSICS_SCATTERING_HIGHLAND_HPP
#define FATRAS_PHYSICS_SCATTERING_HIGHLAND_HPP

#include "ACTS/Utilities/MaterialInteraction.hpp"
#include "Fatras/Kernel/FatrasDefinitions.hpp"
#include "Fatras/Kernel/RandomNumberDistributions.hpp"

namespace Fatras {
      
  /// The Formula used is the highland formula for the projected scattering angle:
  ///
  /// @f$ \theta_{ms} = \frac{13.6MeV}{p}\cdot\sqrt{t/X_{0}}[1 +
  /// 0.038\ln(t/X_{0})]
  /// @f$
  ///
  /// What is returned is the square of the expectation value of the deflection
  /// @f$ < (\theta_ms)^2 > = \sigma_ms^2 @f$
  ///
  /// Public for the benefit of the Mixture scattererers
  template <typename detector_t, typename particle_t>
  double 
  highlandSigma(const detector_t& detector, 
                particle_t& particle, 
                bool log_include=true) 
  {
    // thickness in X0 
    double tInX0 = detector.thickness/detector.material.X0();    
    // get the sigma from 
    double sigma = Acts::sigmaMS(tInX0, particle.p, particle.beta);

    if (std::abs(particle.pdg) == 11) {        
      // Electron multiple scattering effects 
      /// Source: Highland NIM 129 (1975)  p497-499
      // (Highland extension to the Rossi-Greisen formulation)
      // @note: The formula can be extended by replacing 
      // the momentum^2 term with pi * pf
      double sigma2 = constants::main_RossiGreisen / (particle.beta * particle.p);
      sigma2 *= (sigma2 * tInX0);
      // logarithmic term
      if (log_include) {
        double factor = 1. + constants::log_RossiGreisen * log10(10. * tInX0);
        factor *= factor;
        sigma2 *= factor;
      }
      sigma = std::sqrt(sigma2);
    }
    return sigma;
  }  
  
  // The struct for the Scatterer physics list
  struct Highland {
  
    /// Steering parameter
    bool log_include = true;
  
    /// Call operator to perform this scattering
    /// 
    /// @tparam generator_t is a random number generator type 
    /// @tparam detector_t is the detector information type
    /// @tparam particle_t is the particle information type 
    ///
    /// @param[in] generator is the random number generator
    /// @param[in] detector the detector information 
    /// @param[in] particle the particle which is being scattered
    /// 
    /// @return a scattering angle in 3D 
    template <typename generator_t, typename detector_t, typename particle_t>
    double
    operator()(generator_t& generator,
               const detector_t& detector,
               particle_t& particle) const 
    {      
      // Gauss distribution
      Fatras::GaussDist gaussDist(0., 1.);
      // returned scaled by the projection factor
      return M_SQRT2 * highlandSigma(detector,particle,log_include) * gaussDist(generator);
    }
          
  };

} // end of namespace Fatras

#endif // FATRAS_SAMPLER_HIGHLAND_HPP
