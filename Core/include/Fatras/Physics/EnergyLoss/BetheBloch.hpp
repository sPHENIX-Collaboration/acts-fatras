// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef FATRAS_SAMPLER_BETHEBLOCH_HPP
#define FATRAS_SAMPLER_BETHEBLOCH_HPP

#include "ACTS/Utilities/MaterialInteraction.hpp"
#include "Fatras/Kernel/FatrasDefinitions.hpp"
#include "Fatras/Kernel/RandomNumberDistributions.hpp"

namespace Fatras {
      
  // The struct for the EnergyLoss physics list
  struct BetheBloch {
    
    /// MOP / Sigma scaling 
    double scaleFactorMPV   = 1.;
    double scaleFactorSigma = 1.;
  
    /// Call operator
    /// 
    /// @tparam generator_t is a random number generator type 
    /// @tparam detector_t is the detector information type
    /// @tparam particle_t is the particle information type 
    ///
    /// @param[in] generator is the random number generator
    /// @param[in] detector the detector information 
    /// @param[in] particle the particle which is being scattered
    /// 
    /// @return empty vector for BetheBloch
    template <typename generator_t, typename detector_t, typename particle_t>
    std::vector<particle_t>
    operator()(generator_t& generator,
               const detector_t& detector,
               particle_t& particle) const 
    {      

      // Evaluate the energy loss and its sigma
      auto eLoss = Acts::ionizationEnergyLoss_mop(particle.p,
                                                  particle.m,
                                                  detector.material,
                                                  detector.thickness);
      double energyLoss = eLoss.first;
      // the uncertainty of the mean energy loss
      double energyLossSigma = eLoss.second;
      
      // Create a random landau distribution between in the intervall [0,1]
      Fatras::LandauDist landauDist(0., 1.);
      double             landau = landauDist(generator);
      // Simulate the energy loss
      double deltaE = scaleFactorMPV * std::fabs(energyLoss) 
              + scaleFactorSigma * energyLossSigma * landau;
      // protection due to straggling
      // - maximum energy loss is E-m, particle goes to rest
      if (particle.E-deltaE < particle.m){
        particle.E  = particle.m;
        particle.p  = 0.;
        particle.pT = 0.;
        particle.momentum = Acts::Vector3D(0.,0.,0.);
      } else {
        particle.E        -= deltaE;
        particle.p         = std::sqrt(particle.E*particle.E-particle.m*particle.m);
        particle.momentum  = particle.p * particle.momentum.unit();
        particle.pT        = particle.momentum.perp();
      }
      return {};
    }
          
  };

} // end of namespace Fatras

#endif // FATRAS_SAMPLER_BETHEBLOCH_HPP
