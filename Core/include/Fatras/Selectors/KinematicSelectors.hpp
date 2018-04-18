// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef FATRAS_SELECTORS_KINEAMTICSELECTORS_HPP
#define FATRAS_SELECTORS_KINEAMTICSELECTORS_HPP

namespace Fatras {

  /// The Eta cast operator
  struct Eta {
        
    template <typename particle_t>    
    double
    operator()(const particle_t& particle) const {
      return particle.eta();
    }
  
  };

  /// The Eta cast operator
  struct AbsEta {
        
    template <typename particle_t>    
    double
    operator()(const particle_t& particle) const {
      return std::abs(particle.eta());
    }
  
  };

  /// The Eta cast operator
  struct Pt {
        
    template <typename particle_t>    
    double
    operator()(const particle_t& particle) const {
      return particle.eta();
    }
  
  };

}

#endif // FATRAS_SELECTORS_KINEAMTICSELECTORS_HPP
