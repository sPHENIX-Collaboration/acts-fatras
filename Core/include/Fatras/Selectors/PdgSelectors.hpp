// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>

namespace Fatras {

template <int Pdg>
struct AbsPdgSelector
{
  /// Return true for particles matching the absolute pdg number.
  template <typename detector_t>
  constexpr bool
  operator()(const detector_t&, const Particle& particle) const
  {
    return (particle.pdg() == std::abs(Pdg));
  }
};

template <int Pdg>
struct PdgSelector
{
  /// Return true for particles matching the pdg number.
  template <typename detector_t>
  constexpr bool
  operator()(const detector_t&, const Particle& particle) const
  {
    return (particle.pdg() == Pdg);
  }
};

template <int Pdg>
struct AbsPdgExcluder
{
  /// Return true for particles not maching the absolute pdg number.
  template <typename detector_t>
  constexpr bool
  operator()(const detector_t&, const Particle& particle) const
  {
    return (particle.pdg() != std::abs(Pdg));
  }
};

template <int Pdg>
struct PdgExcluder
{
  /// Return true for particles not matching the pdg number.
  template <typename detector_t>
  constexpr bool
  operator()(const detector_t&, const Particle& particle) const
  {
    return (particle.pdg() != Pdg);
  }
};

}  // namespace Fatras
