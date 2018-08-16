// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/GeometryID.hpp"
#include "Acts/Utilities/Units.hpp"
#include <cmath>

namespace Fatras {

/// Typedef the pdg code
typedef int pdg_type;

// Typedef the process code
typedef unsigned int process_code;

/// @brief void Truth particle link
/// - this can be used to link to any underlying
///   truth definition in an application
struct VoidTruthLink {};

/// @brief Particle information struct for phsycis process samplers:
/// - all quatities are calculated at first construction as they may
///   be used by downstream samplers
///
/// @note if a sampler changes one of the parameters, consistency
/// can be broken, so it should update the rest (no checking done)
template <typename truth_particle_t = VoidTruthLink> struct ParticleInfo {

  Acts::Vector3D position = Acts::Vector3D(0.,0.,0.); //!< kinematic info
  Acts::Vector3D momentum = Acts::Vector3D(0.,0.,0.); //!< kinematic info

  double m = 0.;     //!< particle mass
  double E = 0.;     //!< total energy
  double q = 0.;     //!< the charge
  double beta = 0.;  //!< relativistic beta factor
  double gamma = 1.; //!< relativistic gamma factor
  double p = 0.;     //!< momentum magnitude
  double pT = 0.;    //!< transverse momentum magnitude
  int pdg = 0;       //!< pdg code of the particle
  int barcode = 0;   //!< barcode of the particle

  double pathInX0 = 0.;  //!< passed path in X0
  double limitInX0 = 0.; //!< path limit in X0
  double pathInL0 = 0.;  //!< passed path in L0
  double limitInL0 = 0.; //!< path limit in X0
  bool   alive = true;   //!< the particle is alive

  // Truth link to external truth tree
  truth_particle_t truthLink; //!< link to truth particle

  /// Default
  ParticleInfo() {}

  /// @brief Construct a particle consistently
  ///
  /// @param pposition The particle position at construction
  /// @param pmomentum The particle momentum at construction
  /// @param pm The particle mass
  /// @param pq The partilce charge
  /// @param pbarcode The particle barcode
  /// @param ptruthLink The (optional) link to external truth tree
  ParticleInfo(const Acts::Vector3D pposition, 
               const Acts::Vector3D &pmomentum,
               double pm, double pq, 
               int ppdg = 0, int pbarcode = 0,
               truth_particle_t ptruthLink = truth_particle_t())
      : position(pposition)
      , momentum(pmomentum)
      , m(pm)
      , q(pq)
      , p(pmomentum.mag())
      , pT(pmomentum.perp())
      , pdg(ppdg)
      , barcode(pbarcode)
      , truthLink(std::move(ptruthLink))
  {
    E = std::sqrt(p * p + m * m);
    beta = (p / E);
  }

  /// @brief Update the particle with a new position and momentum
  /// 
  /// @param pposition New position after update
  /// @param pmomentum New momentum after update
  void 
  update(const Acts::Vector3D pposition, const Acts::Vector3D &pmomentum) 
  {
    position = pposition;
    momentum = pmomentum;
    p    = pmomentum.mag(); 
    pT   = pmomentum.perp();
    E    = std::sqrt(p * p + m * m);
    beta = (p / E);
  }

  /// boost the particle
  /// @todo implement boost definition
  /// void boost()
  
};

/// @brief Vertex information struct for phsycis process samplers:
/// - all quatities are calculated at first construction as they may
///   be used by downstream samplers
///
/// @note if a sampler changes one of the parameters, consistency
/// can be broken, so it should update the rest (no checking done)
template <typename particle_t> struct VertexInfo {

  /// The vertex position
  Acts::Vector3D position = Acts::Vector3D(0., 0., 0.);
  /// The ingoing particles in the vertex
  std::vector<particle_t> in = {};
  /// The outgoing particles from the vertex
  std::vector<particle_t> out = {};
  /// An optional process code
  process_code processCode = 9;
  /// An optional time stamp
  double timeStamp = 0.;

  /// Default
  VertexInfo() {}

  /// @brief construct a particle consistently
  ///
  /// @param vertex The vertex position
  /// @param in The ingoing particles - copy
  /// @param out The outgoing particles (copy - can we do a move ?)
  /// @param vprocess The process code
  /// @param time The time stamp of this vertex
  VertexInfo(const Acts::Vector3D &vertex,
             const std::vector<particle_t> &ingoing = {},
             std::vector<particle_t> outgoing = {}, process_code process = 0,
             double time = 0.)
      : position(vertex), in(ingoing), out(outgoing), processCode(process),
        timeStamp(time) {}
        
  /// Forward the particle access to the outgoing particles: begin
  typename std::vector<particle_t>::iterator outgoing_begin() 
  { 
    return out.begin();
  }      

  /// Forward the particle access to the outgoing particles: end
  typename std::vector<particle_t>::iterator outgoing_end() 
  { 
    return out.end();
  }
  
  /// Forward the particle access to the outgoing particles: insert
  typename std::vector<particle_t>::iterator 
  outgoing_insert(const std::vector<particle_t>& inparticles)
  {     
    return out.insert(out.end(),inparticles.begin(),inparticles.end());
  }
  
  
};

/// typedef for particle w/o external truth link
typedef ParticleInfo<> Particle;

/// typef for vertex w/o external truth link
typedef VertexInfo<Particle> Vertex;

} // namespace Fatras
