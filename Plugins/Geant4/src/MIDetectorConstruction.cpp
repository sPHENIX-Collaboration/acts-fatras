// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// This code is based on the B1 example of Geant4

#include "Fatras/Plugins/Geant4/MIDetectorConstruction.hpp"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

MIDetectorConstruction::MIDetectorConstruction(G4Material *material,
                                               double thickness)
    : G4VUserDetectorConstruction(), fScoringVolume(0), m_material(material),
      m_thickness(thickness) {}

G4VPhysicalVolume *MIDetectorConstruction::Construct() {
  // Get nist material manager
  G4NistManager *nist = G4NistManager::Instance();

  // Envelope and world parameters
  G4double env_sizeXY = 2. * m, env_sizeZ = 2. * m_thickness * cm / 10.;
  G4Material *env_mat = nist->FindOrBuildMaterial("G4_Galactic");

  // Option to switch on/off checking of volumes overlaps
  G4bool checkOverlaps = true;

  // World
  G4Box *solidWorld = new G4Box("World", // its name
                                0.5 * env_sizeXY, 0.5 * env_sizeXY,
                                0.5 * env_sizeZ); // its size

  G4LogicalVolume *logicWorld = new G4LogicalVolume(solidWorld, // its solid
                                                    env_mat,    // its material
                                                    "World");   // its name

  G4VPhysicalVolume *physWorld =
      new G4PVPlacement(0,               // no rotation
                        G4ThreeVector(), // at (0,0,0)
                        logicWorld,      // its logical volume
                        "World",         // its name
                        0,               // its mother  volume
                        false,           // no boolean operation
                        0,               // copy number
                        checkOverlaps);  // overlaps checking

  // Envelope
  G4Box *solidEnv = new G4Box("Envelope", // its name
                              0.5 * env_sizeXY, 0.5 * env_sizeXY,
                              0.5 * env_sizeZ); // its size

  G4LogicalVolume *logicEnv = new G4LogicalVolume(solidEnv,    // its solid
                                                  env_mat,     // its material
                                                  "Envelope"); // its name

  new G4PVPlacement(0,               // no rotation
                    G4ThreeVector(), // at (0,0,0)
                    logicEnv,        // its logical volume
                    "Envelope",      // its name
                    logicWorld,      // its mother  volume
                    false,           // no boolean operation
                    0,               // copy number
                    checkOverlaps);  // overlaps checking

  // Detector
  G4ThreeVector posDetector = G4ThreeVector(0., 0., 0.25 * env_sizeZ);
  G4Box *solidDetector = new G4Box("Detector", 0.5 * env_sizeXY,
                                   0.5 * env_sizeXY, 0.25 * env_sizeZ);
  G4LogicalVolume *logicDetector =
      new G4LogicalVolume(solidDetector, m_material, "Detector");
  new G4PVPlacement(0, posDetector, logicDetector, "Detector", logicEnv, false,
                    0, checkOverlaps);

  // Set Shape2 as scoring volume
  fScoringVolume = logicDetector;

  // always return the physical World
  return physWorld;
}

G4LogicalVolume *MIDetectorConstruction::GetScoringVolume() const {
  return fScoringVolume;
}
