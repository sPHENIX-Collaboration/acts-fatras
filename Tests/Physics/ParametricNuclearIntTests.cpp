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

#include "../Common/Particle.hpp"
#include "ParametricNuclearIntStub.hpp"
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
    
#include <fstream>
#include <random>
#include <chrono>
#include <array>

namespace bdata = boost::unit_test::data;
namespace tt = boost::test_tools;

namespace Fatras {

namespace Test {

struct FixValueGenerator {

	FixValueGenerator(double value) : m_value(value) {}

	double m_value;
	
	double operator()()
	{
		return m_value;
	}
};

// the generator
typedef std::mt19937 Generator;

/// Generator in [0,1]
struct MyGenerator {

MyGenerator(int samen)
{
	m_generator.seed(samen);
}

// standard generator
Generator m_generator;

	double operator()()
	{
		return (double) m_generator() / (double) m_generator.max();
	}
};

/// The selector
struct MySelector {

  /// call operator
  template <typename detector_t, typename particle_t>
  bool operator()(const detector_t &, const particle_t &) const {
    return true;
  }
};

std::string material = "G4_Be";
std::string gunAmmo = "pi+";
double detectorThickness = 5.; // in [cm]

std::ofstream ofs("fatrasout.txt");
std::ofstream ofsResetter("geant4out.txt");

  G4RunManager* runManager = new G4RunManager;
  G4VModularPhysicsList* physicsList = new QBBC;
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* parDef;

BOOST_AUTO_TEST_CASE(ParamNucularInt_functions_test_) {

ParametricNuclearIntStub pnis;

{
	BOOST_TEST(pnis.nuclearInteractionProbStub(4. * Acts::units::_GeV, 0.5, 2112) > 0.);
	BOOST_TEST(pnis.nuclearInteractionProbStub(2. * Acts::units::_GeV, 0.4, 2112) > 0.);
	BOOST_TEST(pnis.nuclearInteractionProbStub(2. * Acts::units::_GeV, 0.6, 2112) > 0.);

	BOOST_TEST(pnis.nuclearInteractionProbStub(4. * Acts::units::_GeV, 0.5, 2212) > 0.);
	BOOST_TEST(pnis.nuclearInteractionProbStub(2. * Acts::units::_GeV, 0.4, 2212) > 0.);
	BOOST_TEST(pnis.nuclearInteractionProbStub(2. * Acts::units::_GeV, 0.6, 2212) > 0.);

	BOOST_TEST(pnis.nuclearInteractionProbStub(5. * Acts::units::_GeV, 0.5, -211) > 0.);
	BOOST_TEST(pnis.nuclearInteractionProbStub(3. * Acts::units::_GeV, 0.4, -211) > 0.);
	BOOST_TEST(pnis.nuclearInteractionProbStub(3. * Acts::units::_GeV, 0.6, -211) > 0.);

	BOOST_TEST(pnis.nuclearInteractionProbStub(5. * Acts::units::_GeV, 0.5, 211) > 0.);
	BOOST_TEST(pnis.nuclearInteractionProbStub(3. * Acts::units::_GeV, 0.4, 211) > 0.);
	BOOST_TEST(pnis.nuclearInteractionProbStub(3. * Acts::units::_GeV, 0.6, 211) > 0.);

	BOOST_TEST(pnis.nuclearInteractionProbStub(5. * Acts::units::_GeV, 0.5, 321) > 0.);
	BOOST_TEST(pnis.nuclearInteractionProbStub(3. * Acts::units::_GeV, 0.4, 321) > 0.);
	BOOST_TEST(pnis.nuclearInteractionProbStub(3. * Acts::units::_GeV, 0.6, 321) > 0.);

	BOOST_TEST(pnis.nuclearInteractionProbStub(5. * Acts::units::_GeV, 0.5, -321) > 0.);
	BOOST_TEST(pnis.nuclearInteractionProbStub(3. * Acts::units::_GeV, 0.4, -321) > 0.);
	BOOST_TEST(pnis.nuclearInteractionProbStub(3. * Acts::units::_GeV, 0.6, -321) > 0.);

	BOOST_TEST(pnis.nuclearInteractionProbStub(5. * Acts::units::_GeV, 0.4, 130) > 0.);
	BOOST_TEST(pnis.nuclearInteractionProbStub(5. * Acts::units::_GeV, 0.6, 130) > 0.);
	BOOST_TEST(pnis.nuclearInteractionProbStub(3. * Acts::units::_GeV, 0.4, 130) > 0.);
	BOOST_TEST(pnis.nuclearInteractionProbStub(3. * Acts::units::_GeV, 0.6, 130) > 0.);

	BOOST_TEST(pnis.nuclearInteractionProbStub(5. * Acts::units::_GeV, 0.4, 310) > 0.);
	BOOST_TEST(pnis.nuclearInteractionProbStub(5. * Acts::units::_GeV, 0.6, 310) > 0.);
	BOOST_TEST(pnis.nuclearInteractionProbStub(3. * Acts::units::_GeV, 0.4, 310) > 0.);
	BOOST_TEST(pnis.nuclearInteractionProbStub(3. * Acts::units::_GeV, 0.6, 310) > 0.);

	BOOST_TEST(pnis.nuclearInteractionProbStub(5. * Acts::units::_GeV, 0.4, 311) > 0.);
	BOOST_TEST(pnis.nuclearInteractionProbStub(5. * Acts::units::_GeV, 0.6, 311) > 0.);
	BOOST_TEST(pnis.nuclearInteractionProbStub(3. * Acts::units::_GeV, 0.4, 311) > 0.);
	BOOST_TEST(pnis.nuclearInteractionProbStub(3. * Acts::units::_GeV, 0.6, 311) > 0.);

	BOOST_TEST(pnis.nuclearInteractionProbStub(1. * Acts::units::_GeV, 0.5, 1337) == 0.);
}

FixValueGenerator fvg0(0.);
FixValueGenerator fvg1(1.);

Acts::Vector3D pos{0., 0., 0.};
Acts::Vector3D mom{0., 0., 1. * Acts::units::_GeV};
double mass = 0.1395701 * Acts::units::_GeV;
double charge = -1.;
int pdg = -211;
int barcode = 0.;
double time = 0.;
Particle p(pos, mom, mass, charge, pdg, barcode, time);

pnis.nuclearInteractionStub(fvg0, 0.5, p); 

//~ hadronSurvivesStub(const double momentum, const double thickness, const int pdg)
//~ multiplicityProbStub(const double momentum, const double thickness, const int pdg, const unsigned int mult) 
//~ multiplicityStub(generator_t& generator, const double thickness, particle_t& particle) 
//~ particleCompositionStub(generator_t& generator, const int pdg, const unsigned int nParticles) 
//~ kinematicsStub(generator_t& generator, particle_t& particle, const std::vector<int>& particlesPDGs)
//~ finalStateHadronsStub(generator_t& generator, const double thicknessInL0, particle_t& particle)
}  


/// Test the scattering implementation
//~ BOOST_DATA_TEST_CASE(
    //~ ParamNucularInt_test_,
    //~ bdata::random(
        //~ (bdata::seed = 20,
         //~ bdata::distribution = std::uniform_real_distribution<>(0., 1.))) ^
        //~ bdata::random(
            //~ (bdata::seed = 21,
             //~ bdata::distribution = std::uniform_real_distribution<>(0., 1.))) ^
        //~ bdata::random(
            //~ (bdata::seed = 22,
             //~ bdata::distribution = std::uniform_real_distribution<>(0., 1.))) ^
        //~ bdata::random((bdata::seed = 23,
                       //~ bdata::distribution =
                           //~ std::uniform_real_distribution<>(0.5, 10.5))) ^
        //~ bdata::xrange(10000),
    //~ x, y, z, p, index) {
		

/**
if(index ==0)
{
  G4Random::setTheEngine(new CLHEP::RanecuEngine);
  physicsList->SetVerboseLevel(0);
  runManager->SetVerboseLevel(0);
  runManager->SetUserInitialization(new B1DetectorConstruction(material, detectorThickness));
  runManager->SetUserInitialization(physicsList);
  runManager->SetUserInitialization(new B1ActionInitialization(detectorThickness));
	ofsResetter.close();
	parDef = particleTable->FindParticle(gunAmmo);
}
// TODO: record timings
// TODO: G4NistManager could allow access to material properties

// some material
G4NistManager* nist = G4NistManager::Instance();
G4Material* g4mat = nist->FindOrBuildMaterial(material);
Acts::Material actsMaterial = Acts::Material(g4mat->GetRadlen(), g4mat->GetNuclearInterLength(), g4mat->GetA() * mole / g, g4mat->GetZ(), g4mat->GetDensity() * cm3 / kg);

	double x = 0., y = 0., z = 1., p = 1.;

  // create the particle and set the momentum
  /// position at 0.,0.,0
  Acts::Vector3D position{0., 0., 0.};
  Acts::Vector3D direction = Acts::Vector3D(x, y, z).unit();
  // p of 1 GeV
  Acts::Vector3D momentum =
      p * Acts::units::_GeV * direction;

  UImanager->ApplyCommand("/run/initialize");
  UImanager->ApplyCommand("/gun/particle " + gunAmmo);
  //~ UImanager->ApplyCommand("/gun/energy " + std::to_string(p) + " GeV");
  //~ UImanager->ApplyCommand("/gun/direction " + std::to_string(direction.x())
		//~ + " " + std::to_string(direction.y())
		//~ + " " + std::to_string(direction.z()));
  UImanager->ApplyCommand("/gun/momentum " + std::to_string(p * direction.x())
	+ " " + std::to_string(p * direction.y())
	+ " " + std::to_string(p * direction.z()) + " GeV");
  UImanager->ApplyCommand("/gun/position 0. 0. 0.");
  UImanager->ApplyCommand("/gun/time 0.");
  UImanager->ApplyCommand("/tracking/verbose 1");
  UImanager->ApplyCommand("/run/beamOn 1");



MyGenerator mg(index);

  // a detector
  Acts::MaterialProperties detector(actsMaterial, detectorThickness * Acts::units::_cm);

std::cout << detector << std::endl;

  // create the particle
  Particle particle(position, momentum, parDef->GetPDGMass() * Acts::units::_MeV, parDef->GetPDGCharge(), parDef->GetPDGEncoding(), 1);

ParametricNuclearInt paramNuclInt;

std::vector<Particle> par = paramNuclInt(mg, detector, particle);
for(size_t i = 0; i < par.size(); i++)
	ofs << par[i].pdg() << "\t" << par[i].m() << "\t" << par[i].q() << "\t" << par[i].E() << "\t" 
		<< par[i].position().x() << "\t" << par[i].position().y() << "\t" << par[i].position().z() << "\t"
		<< par[i].momentum().x() << "\t" << par[i].momentum().y() << "\t" << par[i].momentum().z() << std::endl;
ofs << "*" << std::endl;

  typedef MySelector All;
  std::vector<Particle> outgoing;
  typedef Process<ParametricNuclearInt, All, All, All> HadronProcess;
  PhysicsList<HadronProcess> hsPhysicsList;
  hsPhysicsList(mg, detector, particle, outgoing);
  
}
**/
} // namespace Test

} // namespace Fatras

