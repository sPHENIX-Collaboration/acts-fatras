// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// This code is based on the B1 example of Geant4

#pragma once

#include "G4VUserDetectorConstruction.hh"
#include "G4Material.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;

/// @brief Detector construction class to define materials and geometry.
class MIDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
	/// @brief Constructor
	///
	/// @param [in] material Material of the detector
	/// @param [in] thickness Thickness of the material
    MIDetectorConstruction(G4Material* material, double thickness);
    
    /// @brief Destructor
    virtual ~MIDetectorConstruction() = default;

	/// @brief Constructs the world
	///
	/// @return Pointer to the world
    virtual G4VPhysicalVolume* 
    Construct();
    
    /// @brief Getter of the scoring volume
    ///
    /// @return Pointer to the scoriung volume
    G4LogicalVolume* 
    GetScoringVolume() const;

  protected:
	// Pointer to the scoring volume
    G4LogicalVolume*  fScoringVolume;
    // Pointer to the detector material
    G4Material* m_material;
    // Thickness of the material
    double m_thickness;
};

