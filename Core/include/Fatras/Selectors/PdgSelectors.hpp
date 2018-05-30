// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

namespace Fatras {

template <int PDG> struct AbsPdgSelector {

  // absolute Pdg selection
  const int saPDG = PDG;

  /// Return true for all particles with | pdg | matching
  /// the selection criteria
  template <typename detector_t, typename particle_t>
  bool operator()(const detector_t &, const particle_t &particle) const {
    return (particle.pdg * particle.pdg == saPDG * saPDG);
  }
};

template <int PDG> struct PdgSelector {

  // Pdg selection
  const int saPDG = PDG;

  /// Return true for all particles with pdg matching
  /// the selection criteria
  template <typename detector_t, typename particle_t>
  bool operator()(const detector_t &, const particle_t &particle) const {
    return (particle.pdg == saPDG);
  }
};

template <int PDG> struct AbsPdgExcluder {

  // absolute Pdg selection
  const int saPDG = PDG;

  /// Return true for all particles with | pdg | matching
  /// the selection criteria
  template <typename detector_t, typename particle_t>
  bool operator()(const detector_t &, const particle_t &particle) const {
    return !(particle.pdg * particle.pdg == saPDG * saPDG);
  }
};

template <int PDG> struct PdgExcluder {

  // Pdg selection
  const int saPDG = PDG;

  /// Return true for all particles with pdg matching
  /// the selection criteria
  template <typename detector_t, typename particle_t>
  bool operator()(const detector_t &, const particle_t &particle) const {
    return !(particle.pdg == saPDG);
  }
};

} // namespace Fatras
