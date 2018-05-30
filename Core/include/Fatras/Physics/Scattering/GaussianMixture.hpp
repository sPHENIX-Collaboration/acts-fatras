// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef FATRAS_SAMPLER_GAUSSIANMIXTURE_HPP
#define FATRAS_SAMPLER_GAUSSIANMIXTURE_HPP

#include "Fatras/Kernel/FatrasDefinitions.hpp"
#include "Fatras/Kernel/RandomNumberDistributions.hpp"
#include "Fatras/Scattering/Highland.hpp"

namespace Fatras {
      
 
  /// The struct to be provided to the Scatterer action
  struct GaussianMixture {
  
    /// Steering parameter
    bool log_include = true;
  
    double gausMixSigma1_a0  = 8.471e-1;
    double gausMixSigma1_a1  = 3.347e-2;
    double gausMixSigma1_a2  = -1.843e-3;
    double gausMixEpsilon_a0 = 4.841e-2;
    double gausMixEpsilon_a1 = 6.348e-3;
    double gausMixEpsilon_a2 = 6.096e-4;
    double gausMixEpsilon_b0 = -1.908e-2;
    double gausMixEpsilon_b1 = 1.106e-1;
    double gausMixEpsilon_b2 = -5.729e-3;
    
    bool optGaussianMixtureG4 = false;
    
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
      Fatras::UniformDist uniformDist(0.,1.);
      
      double sigma = highlandSigma(detector,particle,log_include);
      double sigma2 = sigma*sigma;
      
      // Now correct for the tail fraction 
      // d_0'
      double beta2      = particle.beta*particle.beta;
      double dprime     = detector.thickness/(detector.material.X0() * beta2);
      double log_dprime = std::log(dprime);
      // d_0''
      double log_dprimeprime 
        = std::log(std::pow(detector.material.Z(), 2.0 / 3.0) * dprime);
      // get epsilon
      double epsilon = log_dprimeprime < 0.5
          ? gausMixEpsilon_a0
              + gausMixEpsilon_a1 * log_dprimeprime
              + gausMixEpsilon_a2 * log_dprimeprime * log_dprimeprime
          : gausMixEpsilon_b0
              + gausMixEpsilon_b1 * log_dprimeprime
              + gausMixEpsilon_b2 * log_dprimeprime * log_dprimeprime;
      // the standard sigma
      double sigma1square = gausMixSigma1_a0
          + gausMixSigma1_a1 * log_dprime
          + gausMixSigma1_a2 * log_dprime * log_dprime;
      // G4 optimised / native double Gaussian model
      if (!optGaussianMixtureG4) 
        sigma2 = 225. * dprime / (particle.p * particle.p);
      // throw the random number core/tail
      if (uniformDist(generator) < epsilon) 
        sigma2 *= (1. - (1. - epsilon) * sigma1square) / epsilon;
      
      // return back to the 
      return M_SQRT2 * sqrt(sigma2) * gaussDist(generator);
    }  
    
  };

} // end of namespace Fatras

#endif // FATRAS_SAMPLER_GAUSSIANMIXTURE_HPP
