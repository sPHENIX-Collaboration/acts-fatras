// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef FATRAS_PROCESSES_SCATTERING_HPP
#define FATRAS_PROCESSES_SCATTERING_HPP

#include <memory>

#include <ACTS/EventData/ParticleDefinitions.hpp>
#include <ACTS/Material/MaterialProperties.hpp>
#include <ACTS/Utilities/Definitions.hpp>

#include "Fatras/IMultipleScatteringSampler.hpp"

namespace Fatras {

  /// @class Scatterering
  ///
  /// This is the (multiple) scattering plugin to the
  /// Physics list. It needs a scattering formula in order
  /// to provide the the scattering angle in 3D space.
  /// 
  /// There's two options to apply the scattering
  /// - a parametric action that relates phi and theta (default: off)
  /// - an actuall out of direction scattering applying two random numbers
    
  template <typename Formula>  
  struct Scattering
  {
  
    /// Include the log term
    bool     parametric = false;
    
    /// The scattering formula
    Formula  angle;
    
    /// This is the scattering 
    /// call operator:  
    template <typename cache_t,
              typename generator_t,
              typename detector_t,
              typename particle_t>
    bool
    operator()(cache_t&, 
               generator_t&,
               const detector_t&,
               const particle_t&,             
               std::vector<particle_t>&) const 
    { 
      
      
      

      /// scattering always returns false, it is 
      /// a non-distructive process
      return false; 
    }
    
    
  
  };
  
}  // end of namespace

#endif  // FATRAS_PROCESSES_HIGHLAND_HPP
