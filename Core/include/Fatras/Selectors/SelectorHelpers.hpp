// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <climits>

#include "Fatras/EventData/Particle.hpp"

namespace Fatras {

template <typename cast_t>
struct Min
{
  double valMin = std::numeric_limits<double>::lowest();

  /// Return true for all particles with transverse momentum
  /// bigger than the specified minimum value
  template <typename detector_t>
  constexpr bool
  operator()(const detector_t&, const Particle& particle) const
  {
    return (valMin <= cast_t()(particle));
  }
};

template <typename cast_t>
struct Max
{
  double valMax = std::numeric_limits<double>::max();

  /// Return true for all particles with transverse momentum
  /// bigger than the specified minimum value
  template <typename detector_t>
  constexpr bool
  operator()(const detector_t&, const Particle& particle) const
  {
    return (val < cast_t()(particle));
  }
};

template <typename cast_t>
struct Range
{
  double valMin = std::numeric_limits<double>::lowest();
  double valMax = std::numeric_limits<double>::max();

  /// Return true for all particles with transverse momentum
  /// within the specified range
  template <typename detector_t>
  constexpr bool
  operator()(const detector_t&, const Particle& particle) const
  {
    const auto val = cast_t()(particle);
    return ((valMin <= val) and (val < valMax));
  }
};

}  // namespace Fatras
