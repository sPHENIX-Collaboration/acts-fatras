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

// PhysicsList:
#include "globals.hh"
#include "G4VModularPhysicsList.hh"
#include "G4DataQuestionaire.hh"
#include "G4BetheBlochModel.hh"
#include "G4EmParameters.hh"
#include "G4MuIonisation.hh"
#include "G4hIonisation.hh"

namespace bdata = boost::unit_test::data;
namespace tt = boost::test_tools;

namespace Fatras {
namespace Test {

class G4EmParticleList {

public:

  explicit G4EmParticleList()
  {
  pNames = 
    { 
        "gamma",            "e-",           "e+",           "mu+",        "mu-",
          "pi+",           "pi-",        "kaon+",         "kaon-",     "proton",
  "anti_proton",         "alpha",          "He3",    "GenericIon",         "B+",
           "B-",            "D+",           "D-",           "Ds+",        "Ds-",
     "anti_He3",    "anti_alpha","anti_deuteron","anti_lambda_c+","anti_omega-",
"anti_sigma_c+","anti_sigma_c++",  "anti_sigma+",   "anti_sigma-","anti_triton",
   "anti_xi_c+",      "anti_xi-",     "deuteron",     "lambda_c+",     "omega-",
     "sigma_c+",     "sigma_c++",       "sigma+",        "sigma-",       "tau+",
         "tau-",        "triton",        "xi_c+",           "xi-"
    };
}

  ~G4EmParticleList(){}

  const std::vector<G4String>& PartNames() const {return pNames;}

private:
  std::vector<G4String>  pNames; 
  
};

class BetheBlochPhysics : public G4VPhysicsConstructor // From G4EmStandardPhysics
{
public:
	explicit BetheBlochPhysics(G4int ver=0, const G4String& name="") : G4VPhysicsConstructor("FuckYouG4"), verbose(ver)
	{
	  G4EmParameters* param = G4EmParameters::Instance();
	  param->SetDefaults();
	  param->SetVerbose(verbose);
	  SetPhysicsType(2);
	}
	virtual ~BetheBlochPhysics() {}

	virtual void ConstructParticle()
	{
		G4MuonPlus::MuonPlus();
		G4Gamma::Gamma();
		G4Electron::Electron();
		G4Positron::Positron();
		G4PionPlus::PionPlusDefinition();
		G4PionMinus::PionMinusDefinition();
	}
  virtual void ConstructProcess()
{
  if(verbose > 1) {
    G4cout << "### " << GetPhysicsName() << " Construct Processes " << G4endl;
  }
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
  
  // Add standard EM Processes
  G4ParticleTable* table = G4ParticleTable::GetParticleTable();
  G4EmParticleList* partList = new G4EmParticleList();
  for(const auto& particleName : partList->PartNames()) {
    G4ParticleDefinition* particle = table->FindParticle(particleName);
    if (!particle) { continue; }
    if (particleName == "mu+" ||
               particleName == "mu-"    ) {
      ph->RegisterProcess(new G4MuIonisation(), particle);
    } else if (particleName == "pi+" ||
               particleName == "pi-" ) {
      ph->RegisterProcess(new G4hIonisation(), particle);
    }
  }
}

private:
  G4int  verbose;
};

class MyPhysicsList : public G4VModularPhysicsList
{
public:

  MyPhysicsList(G4int ver = 1, const G4String& type = "MyPhysicsList")
  {
	  G4DataQuestionaire it(photon, neutronxs);
	  
	  G4cout << "<<< Reference Physics List MyPhysicsList " <<G4endl;	

	  defaultCutValue = 0.7*mm;  
	  SetVerboseLevel(ver);
	  
     RegisterPhysics( new BetheBlochPhysics(ver) );


	  //~ RegisterPhysics(new G4EmModelActivator());
  }

  virtual ~MyPhysicsList(){}

  virtual void SetCuts(){  SetCutsWithDefault();  }

private:

  // copy constructor and hide assignment operator
  MyPhysicsList(MyPhysicsList &);
  MyPhysicsList & operator=(const MyPhysicsList &right);

//~ QBBC::QBBC( G4int ver, const G4String&)
//~ {
  //~ G4DataQuestionaire it(photon, neutronxs);
  //~ G4cout << "<<< Reference Physics List QBBC "
	 //~ <<G4endl;	

  //~ defaultCutValue = 0.7*mm;  
  //~ SetVerboseLevel(ver);

  //~ // EM Physics
  //~ RegisterPhysics( new G4EmStandardPhysics(ver) );

  //~ // Synchroton Radiation & GN Physics
  //~ RegisterPhysics( new G4EmExtraPhysics(ver) );

  //~ // Decays
  //~ RegisterPhysics( new G4DecayPhysics(ver) );

   //~ // Hadron Physics
  //~ RegisterPhysics( new G4HadronElasticPhysicsXS(ver) );

  //~ RegisterPhysics( new G4StoppingPhysics(ver) );

  //~ RegisterPhysics( new G4IonPhysics(ver) );

  //~ RegisterPhysics( new G4HadronInelasticQBBC(ver));

  //~ // Neutron tracking cut
  //~ RegisterPhysics( new G4NeutronTrackingCut(ver) );
//~ }
};

std::string material = "G4_Be";
std::string gunAmmo = "pi+";

G4VModularPhysicsList* physicsList = new MyPhysicsList;
G4UImanager* UImanager = G4UImanager::GetUIpointer();
G4RunManager* runManager = new G4RunManager;

BOOST_AUTO_TEST_CASE(step_actor_test)
{
std::ofstream ofs("geant4energyloss.txt");
ofs.close();
double detectorThickness = 100.;
double p = 5. * Acts::units::_GeV;

double x = 0., y = 0., z = 1.;
Acts::Vector3D direction = Acts::Vector3D(x, y, z).normalized();

B1DetectorConstruction* detConstr = new B1DetectorConstruction(material, detectorThickness);
B1ActionInitialization* actionInit = new B1ActionInitialization(detectorThickness, gunAmmo, p * direction.x() * GeV, p * direction.y() * GeV, p * direction.z() * GeV);

//~ if(index == 0)
//~ {
	G4Random::setTheEngine(new CLHEP::RanecuEngine);
	physicsList->SetVerboseLevel(0);
	runManager->SetVerboseLevel(0);
	runManager->SetUserInitialization(physicsList);
	runManager->SetUserInitialization(detConstr);
	runManager->SetUserInitialization(actionInit);
//~ }
//~ else
//~ {
	//~ runManager->DefineWorldVolume(detConstr->Construct());
	//~ runManager->SetUserInitialization(actionInit);
//~ }
	runManager->Initialize();
	UImanager->ApplyCommand("/tracking/verbose 0");
	runManager->BeamOn(10000);
}

} // namespace Test

} // namespace Fatras

