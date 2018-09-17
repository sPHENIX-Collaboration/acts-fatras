// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Fatras/Plugins/Geant4/Geant4MaterialInteraction.hpp"

namespace Fatras {
namespace Test {

/// @brief This class serves for testing @c Geant4MaterialInteraction.
/// All functions just call the parent class functions for public calls.
/// The functionality of the methods are described in @c
/// Geant4MaterialInteraction
class Geant4MaterialInteractionStub : public Geant4MaterialInteraction {
public:
  Geant4MaterialInteractionStub() : Geant4MaterialInteraction() {}

  std::pair<double, double>
  angleOfNormalVectorStub(const Acts::Vector3D normalVector,
                          const Acts::Vector3D momentum) const {
    return angleOfNormalVector(normalVector, momentum);
  }

  template <typename particle_t>
  G4ParticleDefinition *
  convertParticleToG4Stub(const particle_t &particle) const {
    return convertParticleToG4(particle);
  }

  template <typename particle_t>
  G4ParticleGun *
  createParticleGunStub(const particle_t &particle,
                        const std::pair<double, double> &angles) const {
    return createParticleGun(particle, angles);
  }

  template <typename particle_t>
  G4ParticleGun *createParticleGunStub(const particle_t &particle) const {
    return createParticleGun(particle);
  }

  template <typename material_t>
  G4Material *convertMaterialToG4Stub(const material_t &material) const {
    return convertMaterialToG4(material);
  }

  template <typename particle_t>
  void convertParticlesFromG4Stub(const std::vector<MIparticle> &particlesG4,
                                  particle_t &particleIn,
                                  const std::pair<double, double> &angles,
                                  std::vector<particle_t> &particles) const {
    convertParticlesFromG4(particlesG4, particleIn, angles, particles);
  }

  template <typename particle_t>
  void convertParticlesFromG4Stub(const std::vector<MIparticle> &particlesG4,
                                  particle_t &particleIn,
                                  std::vector<particle_t> &particles) const {
    convertParticlesFromG4(particlesG4, particleIn, particles);
  }
};
}
}
