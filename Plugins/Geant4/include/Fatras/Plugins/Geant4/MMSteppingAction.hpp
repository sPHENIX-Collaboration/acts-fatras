// This file is part of the ACTS project.
//
// Copyright (C) 2017 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// MMMaterialStepAction.h
///////////////////////////////////////////////////////////////////

#pragma once

#include <vector>
#include "Acts/Plugins/MaterialPlugins/MaterialStep.hpp"
#include "G4UserSteppingAction.hh"
#include "globals.hh"

namespace FW {
namespace G4 {

  /// @class MMSteppingAction
  ///
  /// @brief Collects the MaterialStep entities
  ///
  /// The MMSteppingAction class is the implementation of the
  /// Geant4 class SteppingAction. It creates extracts the weighted material
  /// of every step and collects all material steps.

  class MMSteppingAction : public G4UserSteppingAction
  {
  public:
    /// Constructor
    MMSteppingAction();

    /// Destructor
    ~MMSteppingAction() override;

    /// Static access method
    static MMSteppingAction*
    Instance();

    /// Interface Method doing the step
    /// @note it creates and collects the MaterialStep entities
    /// @param step is teh Geant4 step of the particle
    void
    UserSteppingAction(const G4Step* step) final override;

    /// Interface reset method
    /// @note it clears the collected step vector
    void
    Reset();

    /// Access to the collected MaterialStep entities
    std::vector<Acts::MaterialStep>
    materialSteps()
    {
      return m_steps;
    }

  private:
    /// Instance of the SteppingAction
    static MMSteppingAction* fgInstance;

    /// The collected MaterialStep entities
    std::vector<Acts::MaterialStep> m_steps;
  };

}  // namespace G4
}  // namespace FW
