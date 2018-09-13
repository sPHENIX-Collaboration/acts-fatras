// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// This code is based on the B1 example of Geant4

#pragma once

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;

/// @brief Run action class
class MIRunAction : public G4UserRunAction {
public:
  /// @brief Constructor
  MIRunAction();

  /// @brief Destructor
  virtual ~MIRunAction() = default;

  /// @brief Initializer of a run
  virtual void BeginOfRunAction(const G4Run *);

  /// @brief Finalizer of a run
  ///
  /// @param [in] run Provides the number of events in the run
  virtual void EndOfRunAction(const G4Run *run);
};
