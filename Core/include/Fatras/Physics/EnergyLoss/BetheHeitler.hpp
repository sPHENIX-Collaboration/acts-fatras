// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/MaterialInteraction.hpp"
#include "Fatras/Kernel/Definitions.hpp"
#include "Fatras/Kernel/RandomNumberDistributions.hpp"

namespace Fatras {

const double log_2 = std::log(2.);

/// The struct for the EnergyLoss physics list
///
/// Bethe-Heitler for electron brem description as described here:
/// "A Gaussian-mixture approximation of the Bethe–Heitler model of electron
/// energy loss by bremsstrahlung" R. Frühwirth
///
struct BetheHeitler {

  double scaleFactor = 1.;

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
  /// @return eventually produced photons
  template <typename generator_t, typename detector_t, typename particle_t>
  std::vector<particle_t> operator()(generator_t &generator,
                                     const detector_t &detector,
                                     particle_t &particle) const {

    double tInX0 = detector.thickness() / detector.material().X0();

    // Take a random gamma-distributed value
    GammaDist gDist(tInX0 / log_2, 1.);

    double u = gDist(generator);

    double z = std::exp(-1. * u);
    double deltaE = std::abs(scaleFactor * particle.E * (z - 1.));

    // protection due to straggling
    // - maximum energy loss is E-m, particle goes to rest
    if (particle.E - deltaE < particle.m) {
      particle.E = particle.m;
      particle.p = 0.;
      particle.pT = 0.;
      particle.momentum = Acts::Vector3D(0., 0., 0.);
    } else {
      particle.E -= deltaE;
      particle.p = std::sqrt(particle.E * particle.E - particle.m * particle.m);
      particle.momentum = particle.p * particle.momentum.unit();
      particle.pT = particle.momentum.perp();
    }
    // todo return photons
    return {};
  }
};

} // namespace Fatras
