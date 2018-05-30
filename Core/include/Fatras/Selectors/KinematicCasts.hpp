// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>

namespace Fatras {

namespace casts {

/// The Eta cast operator
struct eta {

  template <typename particle_t>
  double operator()(const particle_t &particle) const {
    return particle.momentum.eta();
  }
};

/// The Eta cast operator
struct absEta {

  template <typename particle_t>
  double operator()(const particle_t &particle) const {
    return std::abs(particle.momentum.eta());
  }
};

/// The Pt cast operator
struct pT {

  template <typename particle_t>
  double operator()(const particle_t &particle) const {
    return particle.momentum.perp();
  }
};

/// The P cast operator
struct p {

  template <typename particle_t>
  double operator()(const particle_t &particle) const {
    return particle.momentum.mag();
  }
};

/// The E cast operator
struct E {

  template <typename particle_t>
  double operator()(const particle_t &particle) const {
    return particle.E;
  }
};

/// The E cast operator
struct vR {

  template <typename particle_t>
  double operator()(const particle_t &particle) const {
    return particle.position.perp();
  }
};

/// The E cast operator
struct vZ {

  template <typename particle_t>
  double operator()(const particle_t &particle) const {
    return particle.position.z();
  }
};

/// The E cast operator
struct AbsVz {

  template <typename particle_t>
  double operator()(const particle_t &particle) const {
    return std::abs(particle.position.z());
  }
};

} // namespace casts
} // namespace Fatras
