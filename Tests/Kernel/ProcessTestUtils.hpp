// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef FATRAS_PROCESSETST_UTILS_HPP
#define FATRAS_PROCESSETST_UTILS_HPP

#include <random>
#include "ACTS/Utilities/Definitions.hpp"
#include "Fatras/Kernel/FatrasDefinitions.hpp"

namespace Fatras {

namespace Test {

  /// the generator
  typedef std::mt19937 Generator;
  
  /// The particle definition
  typedef ParticleInfo Particle;  
  
  /// The selector 
  struct Selector {
    
    /// call operator 
    template <typename particle_t>
    bool
    operator()(const particle_t&) const 
    { return true; }
  
  };
  
  /// The detector
  struct Detector {   
  };

} // namespace Test
} // namespace Fatras

#endif // FATRAS_PROCESSETST_UTILS_HPP
