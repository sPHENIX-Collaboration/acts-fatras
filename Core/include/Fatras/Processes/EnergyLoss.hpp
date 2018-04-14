// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef FATRAS_PROCESSES_ENERGYLOSS_HPP
#define FATRAS_PROCESSES_ENERGYLOSS_HPP

#include <cmath>

#include "Fatras/Kernel/FatrasDefinitions.hpp"
#include "Fatras/Kernel/RandomNumberDistributions.hpp"
#include "ACTS/Utilities/Definitions.hpp"

namespace Fatras {

  /// @class EnergyLoss
  ///
  /// This is the energy loss plugin to the Physics list.
  /// It takes two selectors, for the incoming particle and a 
  /// second check for the particle on it's way out
      
  template <typename Formula, typename SelectorIn, typename SelectorOut>  
  struct EnergyLoss
  {
      
    /// The scattering formula
    Formula         energyLoss;
    
    /// The selector list
    SelectorIn      selectorIn;
    SelectorOut     selectorOut;
    //SelectorChild   selectorChild;
          
    /// This is the scattering call operator
    template <typename generator_t,
              typename detector_t,
              typename particle_t>
    bool
    operator()(generator_t& gen,
               const detector_t& det,
               particle_t& in,             
               std::vector<particle_t>&) const 
    { 
      // check if scattering applies
      if (selectorIn(in)){
        // get the change in energy
        double eloss = energyLoss(gen,det,in);
        // modify the particle 
        if (in.E + eloss > in.m) {
          // kinematically possible, let's go for it
          in.E = in.E + eloss;
          in.p = std::sqrt(in.E * in.E - in.m * in.m);
          in.momentum = in.p * in.momentum.unit();
          in.pT = in.momentum.perp();
        } else 
          return true;      
      }
      // if the incoming particle does not 
      return (!selectorOut(in)); 
    }
  
  };
  
}  // end of namespace

#endif  // FATRAS_PROCESSES_SCATTERING_HPP
