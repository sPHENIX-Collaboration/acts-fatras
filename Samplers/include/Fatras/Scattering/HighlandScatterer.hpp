// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef FATRAS_SAMPLER_HIGHLAND_HPP
#define FATRAS_SAMPLER_HIGHLAND_HPP

#include "ACTS/Utilities/MaterialInteraction.hpp"

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
  struct HighlandScatterer {
  
    /// Steering parameter
    bool log_include = true;
  
    /// Call operator to perform this scattering
    /// 
    /// @tparam generator_t is a random number generator type 
    /// @tparam detector_t is the detector information type
    /// @tparam particle_t is the particle information type 
    ///
    /// @param[in] gen is the random number generator
    /// @param[in] det the detector information 
    /// @param[in] in the particle which is being scattered
    /// 
    /// @return a scattering angle in 3D 
    template <typename generator_t, typename detector_t, typename particle_t>
    double
    operator()(generator_t& gen,
               const detector_t& det,
               particle_t& particle) const 
    {
      
      // thickness in X0 
      double tInX0 = det.thickness/det.material.X0();    
      // get the sigma from 
      double sigma = Acts::sigmaMS(tInX0, in.p, in.beta);

      if (std::abs(in.pdg) == 11) {        
        // Electron multiple scattering effects 
        /// Source: Highland NIM 129 (1975)  p497-499
        // (Highland extension to the Rossi-Greisen formulation)
        // @note: The formula can be extended by replacing 
        // the momentum^2 term with pi * pf
        double sigma2 = constants::main_RossiGreisen / (in.beta * in.p);
        double sigma2 *= (sigma2 * tInX0);
        // logarithmic term
        if (log_include) {
          double factor = 1. + constants::log_RossiGreisen * log10(10. * tInX0);
          factor *= factor;
          sigma2 *= factor;
        }
      }
      // returned scaled by the projection factor
      return = M_SQRT2 * sigma * gaussDist(gen);
    }
  
  };


} // end of namespace Fatras

#endif // FATRAS_SAMPLER_HIGHLAND_HPP
