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

#include <boost/test/data/test_case.hpp>
// leave blank line

#include <boost/test/output_test_stream.hpp>
// leave blank line

#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialProperties.hpp"

#include "Fatras/Kernel/Particle.hpp"
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
#include <chrono>

namespace bdata = boost::unit_test::data;
namespace tt = boost::test_tools;

namespace Fatras {
namespace Test {
	
std::string material = "G4_Be";
std::string gunAmmo = "kaon0";

G4VModularPhysicsList* physicsList = new QBBC;
G4UImanager* UImanager = G4UImanager::GetUIpointer();
G4RunManager* runManager = new G4RunManager;
  
std::ofstream ofsResetter;
std::ofstream ofsRuntime("runtime.txt");

//~ double l0 = 394.133 / 10.; // Be
//~ std::vector<double> mass = {0.1395701, 0.1349766, 0.1395701, 939.56563 * 1e-3, 938.27231 * 1e-3}; // pi-, pi0, pi+, n, p

/// Test the scattering implementation
//~ BOOST_DATA_TEST_CASE(
    //~ ParamNucularInt_test_,
        //~ bdata::random((bdata::seed = 21,
                       //~ bdata::distribution =
                           //~ std::uniform_real_distribution<>(0.01 * 394.133 / 10., 2. * 394.133 / 10.))) ^
        //~ bdata::random((bdata::seed = 22,
                       //~ bdata::distribution =
                           //~ std::uniform_real_distribution<>(0.5 * 0.1349766, 20. * 0.1349766))) ^
        //~ bdata::xrange(1000),
    //~ detectorThickness, p, index) {
BOOST_DATA_TEST_CASE(
    ParamNucularInt_test_,
        bdata::random((bdata::seed = 21,
                       bdata::distribution =
                           std::uniform_real_distribution<>(0.01 * 394.133 / 10., 0.5 * 394.133 / 10.))) ^ // 0.01 - 2.
        bdata::random((bdata::seed = 22,
                       bdata::distribution =
                           std::uniform_real_distribution<>(0.5, 4.))) ^ // 0.5 - 20.
        bdata::xrange(100),
    detectorThickness, p, index) {

std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();

if(index == 0)
	ofsResetter.open("geant4outNeu.txt");
else
	ofsResetter.open("geant4outNeu.txt", std::ofstream::app);
ofsResetter << "run: " << index << "\t" << detectorThickness << "\t" << p << "\t" << material << "\t" << gunAmmo << std::endl;
ofsResetter.close();

double x = 0., y = 0., z = 1.;
Acts::Vector3D direction = Acts::Vector3D(x, y, z).unit();

B1DetectorConstruction* detConstr = new B1DetectorConstruction(material, detectorThickness);
B1ActionInitialization* actionInit = new B1ActionInitialization(detectorThickness, gunAmmo, p * direction.x() * GeV, p * direction.y() * GeV, p * direction.z() * GeV);

if(index == 0)
{
	G4Random::setTheEngine(new CLHEP::RanecuEngine);
	physicsList->SetVerboseLevel(0);
	runManager->SetVerboseLevel(0);
	runManager->SetUserInitialization(physicsList);
	//~ runManager->SetUserInitialization(new B1DetectorConstruction(material, detectorThickness));
	runManager->SetUserInitialization(detConstr);
	//~ runManager->SetUserInitialization(new B1ActionInitialization(detectorThickness, gunAmmo, p * direction.x() * GeV, p * direction.y() * GeV, p * direction.z() * GeV));
	runManager->SetUserInitialization(actionInit);
}
else
{
	runManager->DefineWorldVolume(detConstr->Construct());
	runManager->SetUserInitialization(actionInit);
}

	runManager->Initialize();
	UImanager->ApplyCommand("/tracking/verbose 0");
	runManager->BeamOn(10000);

	
//~ std::cout << "munni: " << UImanager->GetCurrentValues("/gun/particle") << std::endl;

std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();

ofsRuntime << index << "\t" << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << std::endl;
//thickness, l0, pdg, p
// Thicknes & p can be sampled
// L0 will be used by selecting certain materials
// PDG code will be set for certain particles -> could be performed in a for-loop
delete(detConstr);
delete(actionInit);
}

} // namespace Test

} // namespace Fatras

