// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef FATRAS_SELECTORS_HELPERS_HPP
#define FATRAS_SELECTORS_HELPERS_HPP

#include <climits>

namespace Fatras {
   
  // static selectors 
  template <typename CAST>
  struct Min {
  
    CAST   _cast;      
    double valMin = 0.;
          
    /// Return true for all particles with transverse momentum
    /// bigger than the specified minimum value
    template <typename detector_t, typename particle_t>
    bool
    operator()(const detector_t&, const particle_t& particle) const 
    { 
      double val = _cast(particle);
      return (val >= valMin); 
    }
  };
  
  template <typename CAST>
  struct Max {    

    CAST   _cast;      
    double valMax = std::numeric_limits<double>::max();
        
    /// Return true for all particles with transverse momentum
    /// bigger than the specified minimum value
    template <typename detector_t, typename particle_t>
    bool
    operator()(const detector_t&, const particle_t& particle) const 
    { 
      double val = _cast(particle);
      return (val <= valMax); 
    }
  };

  template <typename CAST>
  struct Range {

    CAST   _cast;      
    double valMin = 0.;
    double valMax = std::numeric_limits<double>::max();
          
    /// Return true for all particles with transverse momentum
    /// within the specified range
    template <typename detector_t, typename particle_t>
    bool
    operator()(const detector_t&, const particle_t& particle) const 
    { 
      double val = _cast(particle);
      return (val >= valMin && val <= valMax); 
    }
  }; 
  
}

#endif // FATRAS_SELECTORS_HELPERS_HPP
