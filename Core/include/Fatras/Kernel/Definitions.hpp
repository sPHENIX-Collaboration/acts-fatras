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
