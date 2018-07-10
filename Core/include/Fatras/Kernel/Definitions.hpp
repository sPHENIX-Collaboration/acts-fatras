// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
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

namespace constants {

// @todo multiply with units

/// KOverA factor in Bethe-Bloch equation [MeV*cm2/gram]
constexpr double ka_BetheBloch = 30.7075;

/// Fine structure constexprant
constexpr double alpha = 1. / 137.;

/// Multiple scattering parameters for MIPS
constexpr double main_RutherfordScott = 13.6;
constexpr double log_RutherfordScott = 0.038;

/// Multiple scattering parameters for electrons
constexpr double main_RossiGreisen = 17.5;
constexpr double log_RossiGreisen = 0.125;

} // namespace constants

/// The information to be writtern out per hit surface
/// - this is the bare information on a surface
struct SensitiveHit {
  const Acts::Surface *surface = nullptr;
  Acts::Vector3D position;
  Acts::Vector3D direction;
  double value = 0.;
};

/// A simple selector to select sensitive surfaces
/// this is for charged particles
struct SensitiveSelector {
  /// boolean operator()
  ///
  /// checks if the geometry id of the
  /// surface has a sensitive bit set
  bool operator()(const Acts::Surface &sf) const {
    return (sf.geoID().value(Acts::GeometryID::sensitive_mask));
  }
};

/// A simple struct to switch sensitive surfaces off
struct NoneSelector {
  /// boolean operator()
  ///
  /// ignores what it needs
  bool operator()(const Acts::Surface &sf) const { return false; }
};

} // namespace Fatras
