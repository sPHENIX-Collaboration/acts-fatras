// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>

#include "Acts/Utilities/Helpers.hpp"
#include "Fatras/EventData/Particle.hpp"

namespace Fatras {
namespace Casts {

  struct eta
  {
    /// Extract the particle pseudo-rapidity.
    constexpr auto
    operator()(const Particle& particle) const
    {
      return Acts::VectorHelpers::eta(particle.direction());
    }
  };

  struct absEta
  {
    /// Extract the absolute particle pseudo-rapidity.
    constexpr auto
    operator()(const Particle& particle) const
    {
      return std::abs(Acts::VectorHelpers::eta(particle.direction()));
    }
  };

  struct pT
  {
    /// Extract the particle transverse momentum.
    constexpr auto
    operator()(const Particle& particle) const
    {
      return Acts::VectorHelpers::perp(particle.direction())
          * particle.momentum();
    }
  };

  struct p
  {
    /// Extract the particle total momentum.
    constexpr auto
    operator()(const Particle& particle) const
    {
      return particle.momentum();)
    }
  };

  struct E
  {
    /// Extract the particle energy.
    constexpr double
    operator()(const Particle& particle) const
    {
      return particle.energy();
    }
  };

  struct vRho
  {
    /// Extract the distance to the origin in the transverse plane.
    constexpr auto
    operator()(const Particle& particle) const
    {
      return Acts::VectorHelpers::perp(particle.position());
    }
  };

  struct vZ
  {
    /// Extract the z-distance to the origin.
    constexpr auto
    operator()(const Particle& particle) const
    {
      return particle.position().z();
    }
  };

  struct AbsVz
  {
    /// Extract the absolute z-distance to the origin.
    constexpr auto
    operator()(const Particle& particle) const
    {
      return std::abs(particle.position().z());
    }
  };

}  // namespace Casts
}  // namespace Fatras
