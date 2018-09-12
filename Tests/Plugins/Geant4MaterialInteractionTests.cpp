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
#include "Acts/Material/Material.hpp"

#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4Material.hh"

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
BOOST_TEST(pGun->GetParticleMomentumDirection().z() / MeV * Acts::units::_GeV == 1. * Acts::units::_GeV);

Acts::Material mat(352.8, 407., 9.012, 4., 1.848 * Acts::units::_g / (Acts::units::_cm * Acts::units::_cm * Acts::units::_cm));

G4Material* matG4 = g4mis.convertMaterialToG4Stub(mat);

BOOST_TEST(matG4->GetA() * mole / g == mat.A());
BOOST_TEST(matG4->GetZ() == mat.Z());
BOOST_TEST(matG4->GetDensity() * cm3 / g * Acts::units::_g / (Acts::units::_cm * Acts::units::_cm * Acts::units::_cm) == mat.rho());

std::vector<B1particle> bps;
B1particle bp;
bp.momentum = {0., 0., 2. * Acts::units::_GeV};
bp.mass = 2. * mass;
bp.charge = 2. * charge;
bp.pdg = 2. * pdg;
bps.push_back(std::move(bp));

std::vector<Particle> particles;
g4mis.convertParticlesFromG4Stub(bps, particle, particles);

BOOST_TEST(particles.size() == 1);
BOOST_TEST(particles[0].position() == Acts::Vector3D::Zero(3));
BOOST_TEST(particles[0].momentum() == bp.momentum);
double tmpD = particles[0].m();
BOOST_TEST(tmpD == bp.mass);
tmpD = particles[0].q();
BOOST_TEST(tmpD == bp.charge);
int tmpI = particles[0].pdg();
BOOST_TEST(tmpI == bp.pdg);

particles.clear();
particles = g4mis(particle, mat, 1. * Acts::units::_m);
BOOST_TEST(particles.size() != 0);
}
} // namespace Test
} // namespace Fatras

