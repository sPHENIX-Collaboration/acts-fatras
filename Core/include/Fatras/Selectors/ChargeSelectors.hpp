// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Fatras/EventData/Particle.hpp"

namespace Fatras {

struct ChargedSelector
{
  /// Return true for all particles with charge != 0.
  template <typename detector_t>
  constexpr bool
  operator()(const detector_t&, const Particle& particle) const
  {
    return (particle.charge() != 0);
  }
};

struct NeutralSelector
{
  /// Return true for all particles with charge == 0.
  template <typename detector_t>
  constexpr bool
  operator()(const detector_t&, const Particle& particle) const
  {
    return (particle.charge() == 0);
  }
};

struct PositiveChargeSelector
{
  /// Return true for all particles with charge > 0.
  template <typename detector_t>
  constexpr bool
  operator()(const detector_t&, const Particle& particle) const
  {
    return (0 < particle.charge());
  }
};

struct NegativeChargeSelector
{
  /// Return true for all particles with charge < 0.
  template <typename detector_t>
  constexpr bool
  operator()(const detector_t&, const Particle& particle) const
  {
    return (particle.charge() < 0);
  }
};

}  // namespace Fatras
