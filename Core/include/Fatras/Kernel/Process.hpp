// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Definitions.hpp"
#include "Fatras/Kernel/Definitions.hpp"
#include "Fatras/Kernel/RandomNumberDistributions.hpp"
#include <cmath>

namespace Fatras {

/// @class Process
///
/// This is plugin for physics processes
///  - scattering
///  - energy loss
///  - pair production
///  - hadronic interaction
///  - decay
///
/// To be plugged into the PhysicsList
///
/// The type (and actual trigger) of the particle
/// and interaction is steered via the Selector list
/// for in and out.
template <typename Physics, typename SelectorIn, typename SelectorOut,
          typename SelectorChild>

struct Process {

  /// The actual physics that is happening
  Physics process;

  /// The selector list
  SelectorIn selectorIn;
  SelectorOut selectorOut;
  SelectorChild selectorChild;

  /// This is the scattering call operator
  template <typename generator_t, typename detector_t, typename particle_t>
  bool operator()(generator_t &gen, const detector_t &det, particle_t &in,
                  std::vector<particle_t> &out) const {
    // check if the process applies
    if (selectorIn(det, in)) {
      // apply energy loss and get eventual children
      auto children = process(gen, det, in);
      if (children.size()) {
        // copy the children that comply with the child selector
        std::copy_if(
            children.begin(), children.end(), out.begin(),
            [this, det](const particle_t &p) { return selectorChild(det, p); });
      }
    }
    // check if this killed the partilce,
    // or pushed below threshold
    return (!selectorOut(det, in));
  }
};

} // namespace Fatras
