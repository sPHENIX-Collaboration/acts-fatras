// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///  Boost include(s)
#define BOOST_TEST_MODULE ParametricNuclearInteraction Tests

#include <boost/test/included/unit_test.hpp>
// leave blank line

#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialProperties.hpp"

#include "Particle.hpp"
#include "Fatras/Physics/HadronicInteraction/ParametricNuclearInt.hpp"
#include "Fatras/Kernel/Process.hpp"
#include "Fatras/Kernel/PhysicsList.hpp"

#include "include/B1DetectorConstruction.hh"
#include "include/B1ActionInitialization.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "QBBC.hh"
#include "Randomize.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleGun.hh"
    
#include <fstream>
#include <random>

#include "globals.hh"
#include "G4VModularPhysicsList.hh"

namespace tt = boost::test_tools;

namespace Fatras {
namespace Test {

/// @brief Random number generator
double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

// Material and particle
std::string material = "G4_Be";
std::string gunAmmo = "pi+";

// List of physics processes
G4VModularPhysicsList* physicsList = new QBBC;
// Manager of geant4
G4UImanager* UImanager = G4UImanager::GetUIpointer();
G4RunManager* runManager = new G4RunManager;

// File writer
std::ofstream ofsResetter;

//~ double l0 = 394.133 / 10.; // L0 of Be in cm
//~ std::vector<double> mass = {0.1395701, 0.1349766, 0.1395701, 939.56563 * 1e-3, 938.27231 * 1e-3}; // pi-, pi0, pi+, n, p

/// Test the scattering implementation
BOOST_AUTO_TEST_CASE(step_extension_vacuum_test) {
// Loop over different configurations
for(unsigned int index = 0; index < 2; index++)
{
	// Dice thickness and momentum
	// L0 will be used by selecting certain materials
	double detectorThickness = fRand(0.01 * 394.133 / 10., 0.5 * 394.133 / 10.); // 0.01 - 2.
	double p = fRand(0.5, 4.); // 0.5 - 20.

	// Setup of data writing
	ofsResetter.open("geant4out_" + std::to_string(index) + ".txt");
	ofsResetter << "run: " << index << "\t" << detectorThickness << "\t" << p << "\t" << material << "\t" << gunAmmo << std::endl;
	ofsResetter.close();

	// Initial direction
	double x = 0., y = 0., z = 1.;
	Acts::Vector3D direction = Acts::Vector3D(x, y, z).normalized();

	// Build detector
	B1DetectorConstruction* detConstr = new B1DetectorConstruction(material, detectorThickness);
	// Set action (and therewith data recording/writing)
	B1ActionInitialization* actionInit = new B1ActionInitialization(detectorThickness, gunAmmo, p * direction.x() * GeV, p * direction.y() * GeV, p * direction.z() * GeV, index);

	// Initialize everything once
	if(index == 0)
	{
		G4Random::setTheEngine(new CLHEP::RanecuEngine);
		physicsList->SetVerboseLevel(0);
		runManager->SetVerboseLevel(0);
		runManager->SetUserInitialization(physicsList);
		runManager->SetUserInitialization(detConstr);
		runManager->SetUserInitialization(actionInit);
	}
	else
	{
		// Set world and actions
		runManager->DefineWorldVolume(detConstr->Construct());
		runManager->SetUserInitialization(actionInit);
	}
		// Launch
		runManager->Initialize();
		UImanager->ApplyCommand("/tracking/verbose 0");
		runManager->BeamOn(10);
	
	// Active delete to reduce memory leaks
	delete(detConstr);
	delete(actionInit);
}
}
} // namespace Test

} // namespace Fatras

