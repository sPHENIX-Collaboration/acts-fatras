// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_FATRAS_FATRASDEFINITIONS_H
#define ACTS_FATRAS_FATRASDEFINITIONS_H

#include <cmath>
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/Units.hpp"
#include "ACTS/Material/Material.hpp"

namespace Fatras {
  
  
  namespace constants {
  
    // @todo multiply with units -> those should go to the samplers
  
    /// KOverA factor in Bethe-Bloch equation [MeV*cm2/gram]
    constexpr double ka_BetheBloch = 30.7075;
  
    /// Fine structure constexprant
    constexpr double alpha = 1. / 137.;
  
    /// Multiple scattering paramters
    constexpr double main_RutherfordScott = 13.6;
    constexpr double log_RutherfordScott  = 0.038;
  
    constexpr double main_RossiGreisen = 17.5;
    constexpr double log_RossiGreisen  = 0.125;
  
  }  // namespace constants

  // Detector information struct for phsycis process samplers
  struct DetectorInfo {
    
      double         thickness;  //!< the scaled & corrected pathlength 
      Acts::Material material;   //!< the material (potentially with composition)
  };

  /// Particle information struct for phsycis process samplers:
  /// - all quatities are calculated at first construction as they may
  ///   be used by downstream samplers 
  ///
  /// @note if a sampler changes one of the parameters, consistency
  /// can be broken, so it should update the rest (no checking done)
  struct ParticleInfo {
    
    Acts::Vector3D position;     //!< particle position at sampler call
    Acts::Vector3D momentum;     //!< particle momentum at sampler call 
    double         q       = 0.; //!< the charge
    double         m       = 0.; //!< particle mass
    double         E       = 0.; //!< total energy
    double         beta    = 0.; //!< relativistic beta factor
    double         gamma   = 1.; //!< relativistic gamma factor 
    double         p       = 0.; //!< momentum magnitude
    double         pT      = 0.; //!< transverse momentum magnitude
    int            pdg     = 0;  //!< pdg code of the particle
    int            barcode = 0;  //!< barcode of the particle
    
    // Construct a particle consistently
    ParticleInfo(const Acts::Vector3D pposition,
                 const Acts::Vector3D& pmomentum,
                 double pq,
                 double pm,
                 int ppdg = 0, int pbarcode = 0) :
       position(pposition),
       momentum(pmomentum),
       q(pq),
       m(pm),
       p(pmomentum.mag()),
       pT(pmomentum.perp()),
       pdg(ppdg),
       barcode(pbarcode) 
    {
      E    = std::sqrt(p * p + m * m);
      beta = (p/E);
    }
      
    /// boost the particle
    /// @todo implement boost definition
    ///void boost()  
      
  };


}  // namespace Fatras

#endif
