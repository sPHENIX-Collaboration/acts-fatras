// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_FATRAS_PARTICLE_H
#define ACTS_FATRAS_PARTICLE_H

#include <cmath>
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"
#include "Acts/Utilities/GeometryID.hpp"
#include "Acts/Surfaces/Surface.hpp"

namespace Fatras {
  
  /// Typedef the pdg code 
  typedef int pdg_type;
  
  // Typedef the process code
  typedef unsigned int process_code;
  
  /// @brief void Truth particle link
  /// - this can be used to link to any underlying
  ///   truth definition in an application    
  struct VoidTruthLink {};

  /// @Particle information struct for phsycis process samplers:
  /// - all quatities are calculated at first construction as they may
  ///   be used by downstream samplers 
  ///
  /// @note if a sampler changes one of the parameters, consistency
  /// can be broken, so it should update the rest (no checking done)
  template <typename TruthParticleLink_type = VoidTruthLink>
  struct ParticleInfo {
    
    // kinematic section
    Acts::Vector3D position  = Acts::Vector3D(0.,0.,0.);//!< particle position
    Acts::Vector3D momentum  = Acts::Vector3D(0.,0.,0.);//!< particle momentum 
    double         q         = 0.;   //!< the charge
    double         m         = 0.;   //!< particle mass
    double         E         = 0.;   //!< total energy
    double         beta      = 0.;   //!< relativistic beta factor
    double         gamma     = 1.;   //!< relativistic gamma factor 
    double         p         = 0.;   //!< momentum magnitude
    double         pT        = 0.;   //!< transverse momentum magnitude
    int            pdg       = 0;    //!< pdg code of the particle
    int            barcode   = 0;    //!< barcode of the particle
    
    // process/simulation handling section
    double         pathInX0  = 0.;   //!< passed path in X0
    double         limitInX0 = 0.;   //!< path limit in X0
    double         pathInL0  = 0.;   //!< passed path in L0
    double         limitInL0 = 0.;   //!< path limit in X0
    bool           alive     = true; //!< the particle is alive
    
    // truth link to outside thruth collection
    TruthParticleLink_type   truthLink; //!< link to truth particle 
    
    /// Default 
    ParticleInfo(){}
    
    /// @brief construct a particle consistently
    /// 
    /// @param pposition The curren particle position
    /// @param pmomentum The curren particle momentum
    /// @param pq The partilce charge
    /// @param pbarcode The particle barcode
    /// @param ptruthLink The (optional) link to a truth tree
    ParticleInfo(const Acts::Vector3D pposition,
                 const Acts::Vector3D& pmomentum,
                 double pq,
                 double pm,
                 int ppdg = 0, 
                 int pbarcode = 0,
                 TruthParticleLink_type ptruthLink = TruthParticleLink_type()) :
       position(pposition),
       momentum(pmomentum),
       q(pq),
       m(pm),
       p(pmomentum.mag()),
       pT(pmomentum.perp()),
       pdg(ppdg),
       barcode(pbarcode),
       truthLink(std::move(ptruthLink))
    {
      E    = std::sqrt(p * p + m * m);
      beta = (p/E);
    }
      
    /// boost the particle
    /// @todo implement boost definition
    ///void boost()  
      
  };

  /// @Particle information struct for phsycis process samplers:
  /// - all quatities are calculated at first construction as they may
  ///   be used by downstream samplers 
  ///
  /// @note if a sampler changes one of the parameters, consistency
  /// can be broken, so it should update the rest (no checking done)
  template <typename Particle_type>
  struct VertexInfo {
    
    /// The vertex position
    Acts::Vector3D position = Acts::Vector3D(0.,0.,0.);
    /// The ingoing particles in the vertex
    std::vector<Particle_type> ingoingParticles = {};
    /// The outgoing particles from the vertex
    std::vector<Particle_type> outgoingParticles = {};
    /// An optional process code 
    process_code      processCode  = 9;
    /// An optional time stamp
    double            timeStamp = 0.;
      
    /// Default 
    VertexInfo(){}
    
    /// @brief construct a particle consistently
    /// 
    /// @param vertex The vertex position
    /// @param ingoing The ingoing particles - copy
    /// @param outgoing The outgoing particles (std::move!)
    /// @param vprocess The process code 
    /// @param time The time stamp of this vertex
    VertexInfo(const Acts::Vector3D&       vertex,
               const std::vector<Particle_type>& ingoing = {},
               std::vector<Particle_type> outgoing = {},
               process_code process = 0,
               double       time = 0.) :
       position(vertex),
       ingoingParticles(ingoing),
       outgoingParticles(std::move(outgoing)),
       processCode(process),
       timeStamp(time)
    {}
  
  };

  /// typedef for particle w/o external truth link
  typedef ParticleInfo<> Particle;    
  
  /// typef for vertex w/o external truth link
  typedef VertexInfo<Particle> Vertex;


}  // namespace Fatras

#endif // ACTS_FATRAS_PARTICLE_H
