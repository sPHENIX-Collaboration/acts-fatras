// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// This code is based on the B1 example of Geant4

#pragma once

#include "G4UserSteppingAction.hh"

class MIEventAction;
struct MIparticle;
class G4LogicalVolume;

/// @brief Stepping action class. This class serves as collector of the final state.
class MISteppingAction : public G4UserSteppingAction
{
  public:
	/// @brief Constructor
	///
	/// @param [in, out] eventAction EventAction handler that stores the final state
	/// @param [in] thickness Thickness of the material
    MISteppingAction(MIEventAction* eventAction, double thickness);
    
    /// @brief Destructor
    virtual ~MISteppingAction() = default;

    /// @brief Method from the base class for boundary checks of particles
    /// 
    /// @param [in] step Step calculated in Geant4
    virtual void 
    UserSteppingAction(const G4Step* step);

  private:
	// Event handler
    MIEventAction*  fEventAction;
    // Checker if collection should be performed in the volume
    G4LogicalVolume* fScoringVolume;
    // Thickness of the material
    double m_thickness;
};