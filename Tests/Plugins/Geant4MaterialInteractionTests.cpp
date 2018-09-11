// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///  Boost include(s)
#define BOOST_TEST_MODULE Geant4MaterialInteraction Tests

#include <boost/test/included/unit_test.hpp>
// leave blank line

#include <boost/test/data/test_case.hpp>
// leave blank line

#include <boost/test/output_test_stream.hpp>
// leave blank line

#include "Geant4MaterialInteractionStub.hpp"
#include "Particle.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"

#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"

namespace bdata = boost::unit_test::data;
namespace tt = boost::test_tools;

namespace Fatras {
namespace Test {

BOOST_DATA_TEST_CASE(Geant4MaterialInteraction_test_, bdata::xrange(1), index) {

Acts::Vector3D position(0., 0., 0.);
Acts::Vector3D momentum(0., 0., 1. * Acts::units::_GeV);
double mass = 0.1395701 * Acts::units::_GeV;
double charge = 1.;
int pdg = 211;
Particle particle(position, momentum, mass, charge, pdg);

Geant4MaterialInteractionStub g4mis;

G4ParticleDefinition* parDef = g4mis.convertParticleToG4Stub(particle);

BOOST_TEST(parDef->GetPDGEncoding() == pdg);
BOOST_TEST(parDef->GetPDGMass() / GeV * Acts::units::_GeV == mass);
BOOST_TEST(parDef->GetPDGCharge() == charge);

G4ParticleGun* pGun = g4mis.createParticleGunStub(particle);

BOOST_TEST(pGun->GetParticleDefinition()->GetPDGEncoding() == pdg);
BOOST_TEST(pGun->GetParticleDefinition()->GetPDGMass() / GeV * Acts::units::_GeV == mass);
BOOST_TEST(pGun->GetParticleDefinition()->GetPDGCharge() == charge);
BOOST_TEST(pGun->GetParticleMomentumDirection().x() == 0.);
BOOST_TEST(pGun->GetParticleMomentumDirection().y() == 0.);
BOOST_TEST(pGun->GetParticleMomentumDirection().z() / MeV == 1. * Acts::units::_GeV);

}
} // namespace Test
} // namespace Fatras

