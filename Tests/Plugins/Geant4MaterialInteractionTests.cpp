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

#include "Acts/Material/Material.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"
#include "Geant4MaterialInteractionStub.hpp"
#include "Particle.hpp"

#include "G4Material.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"

#include "G4HadronElasticPhysics.hh"
#include "G4VPhysicsConstructor.hh"

namespace utf = boost::unit_test;
namespace bdata = boost::unit_test::data;
namespace tt = boost::test_tools;

namespace Fatras {
namespace Test {

BOOST_AUTO_TEST_CASE(Geant4MaterialInteraction_test_, *utf::tolerance(1e-10)) {

  // Create particle
  Acts::Vector3D position(0., 0., 0.);
  Acts::Vector3D momentum(0., 0., 1. * Acts::units::_GeV);
  double mass = 0.1395701 * Acts::units::_GeV;
  double charge = 1.;
  int pdg = 211;
  Particle particle(position, momentum, mass, charge, pdg);

  Geant4MaterialInteractionStub g4mis;
  std::pair<double, double> angleStub(0., 0.);

  // Test conversion of particle
  G4ParticleDefinition *parDef = g4mis.convertParticleToG4Stub(particle);

  BOOST_TEST(parDef->GetPDGEncoding() == pdg);
  BOOST_TEST(parDef->GetPDGMass() / GeV * Acts::units::_GeV == mass);
  BOOST_TEST(parDef->GetPDGCharge() == charge);

  // Test particle gun creation
  G4ParticleGun *pGun = g4mis.createParticleGunStub(particle, angleStub);

  BOOST_TEST(pGun->GetParticleDefinition()->GetPDGEncoding() == pdg);
  BOOST_TEST(pGun->GetParticleDefinition()->GetPDGMass() / GeV *
                 Acts::units::_GeV ==
             mass);
  BOOST_TEST(pGun->GetParticleDefinition()->GetPDGCharge() == charge);
  BOOST_TEST(pGun->GetParticleMomentumDirection().x() == 0.);
  BOOST_TEST(pGun->GetParticleMomentumDirection().y() == 0.);
  BOOST_TEST(pGun->GetParticleMomentumDirection().z() / MeV *
                 Acts::units::_GeV ==
             1. * Acts::units::_GeV);

  // Create material and test its conversion
  Acts::Material mat(
      352.8, 407., 9.012, 4.,
      1.848 * Acts::units::_g /
          (Acts::units::_cm * Acts::units::_cm * Acts::units::_cm));

  G4Material *matG4 = g4mis.convertMaterialToG4Stub(mat);

  BOOST_TEST(matG4->GetA() * mole / g == mat.A());
  BOOST_TEST(matG4->GetZ() == mat.Z());
  BOOST_TEST(matG4->GetDensity() * cm3 / g * Acts::units::_g /
                 (Acts::units::_cm * Acts::units::_cm * Acts::units::_cm) ==
             mat.rho());

  // Test back conversion
  std::vector<MIparticle> bps;
  MIparticle bp;
  bp.momentum = {0., 0., 2. * GeV};
  bp.mass = 2. * 0.1395701 * GeV;
  bp.charge = 2. * charge;
  bp.pdg = 2. * pdg;
  bps.push_back(std::move(bp));

  std::vector<Particle> particles;
  g4mis.convertParticlesFromG4Stub(bps, particle, angleStub, particles);

  BOOST_TEST(particles.size() == 1);
  BOOST_TEST(particles[0].position() == Acts::Vector3D::Zero(3));
  BOOST_TEST(particles[0].momentum() == bp.momentum * Acts::units::_GeV / GeV);
  double tmpD = particles[0].m();
  BOOST_TEST(tmpD == bp.mass * Acts::units::_GeV / GeV);
  tmpD = particles[0].q();
  BOOST_TEST(tmpD == bp.charge);
  int tmpI = particles[0].pdg();
  BOOST_TEST(tmpI == bp.pdg);

  // Test angle transformation
  Acts::Vector3D normal(1., 1., 1.);
  std::pair<double, double> angles =
      g4mis.angleOfNormalVectorStub(normal, momentum);
  BOOST_TEST(std::isfinite(angles.first));
  BOOST_TEST(std::isfinite(angles.second));

  normal = {1., 0., 1.};
  angles = g4mis.angleOfNormalVectorStub(normal, momentum);
  BOOST_TEST(std::isfinite(angles.first));
  BOOST_TEST(std::isfinite(angles.second));

  normal = {0., 0., 1.};
  angles = g4mis.angleOfNormalVectorStub(normal, momentum);
  BOOST_TEST(angles.first == 0.);
  BOOST_TEST(angles.second == 0.);

  Acts::Vector3D negMomentum(0., 0., -1.);
  std::pair<double, double> angles2 =
      g4mis.angleOfNormalVectorStub(normal, negMomentum);
  BOOST_TEST(angles2.first == M_PI);
  BOOST_TEST(angles2.second == 0.);

  normal = {0., 1., 1.};
  angles = g4mis.angleOfNormalVectorStub(normal, momentum);
  BOOST_TEST(std::isfinite(angles.first));
  BOOST_TEST(std::isfinite(angles.second));

  normal = {1., 1., 0.};
  angles = g4mis.angleOfNormalVectorStub(normal, momentum);
  BOOST_TEST(std::isfinite(angles.first));
  BOOST_TEST(std::isfinite(angles.second));

  negMomentum = {-1., -1., 0.};
  angles2 = g4mis.angleOfNormalVectorStub(normal, negMomentum);
  BOOST_TEST(angles2.first - M_PI == angles.first);
  BOOST_TEST(angles2.second == angles.second);

  // Test of the main call operator
  normal = {1., 1., 1.};
  particles.clear();
  particles = g4mis(particle, mat, 1. * Acts::units::_m, normal);
  BOOST_TEST(particles.size() != 0);

  Acts::Vector3D normalZero(0., 0., 0.);
  particles.clear();
  particles = g4mis(particle, mat, 1. * Acts::units::_m, normalZero);
  BOOST_TEST(particles.size() == 0);

  particles.clear();
  particles = g4mis(particle, mat, 0.001 * Acts::units::_nm, normal);
  BOOST_TEST(particles.size() == 1);
  BOOST_TEST(particles[0].momentum().x() == particle.momentum().x(),
             tt::tolerance(1e-8));
  BOOST_TEST(particles[0].momentum().y() == particle.momentum().y(),
             tt::tolerance(1e-8));
  BOOST_TEST(particles[0].momentum().z() == particle.momentum().z(),
             tt::tolerance(1e-8));

  // Stability test for wrong particle data
  Particle wrongParticle(position, momentum, mass, charge, 0);
  BOOST_TEST(!g4mis.convertParticleToG4Stub(wrongParticle));
  BOOST_TEST(!g4mis.createParticleGunStub(wrongParticle, angleStub));
  particles.clear();
  particles = g4mis(wrongParticle, mat, 1. * Acts::units::_m, normal);
  BOOST_TEST(particles.size() == 0);

  // Stability test for wrong material data
  Acts::Material wrongMat1(
      352.8, 407., -1, 4.,
      1.848 * Acts::units::_g /
          (Acts::units::_cm * Acts::units::_cm * Acts::units::_cm));
  BOOST_TEST(!g4mis.convertMaterialToG4Stub(wrongMat1));
  particles.clear();
  particles = g4mis(particle, wrongMat1, 1. * Acts::units::_m, normal);
  BOOST_TEST(particles.size() == 0);

  Acts::Material wrongMat2(
      352.8, 407., 9.012, -1.,
      1.848 * Acts::units::_g /
          (Acts::units::_cm * Acts::units::_cm * Acts::units::_cm));
  BOOST_TEST(!g4mis.convertMaterialToG4Stub(wrongMat2));
  particles.clear();
  particles = g4mis(particle, wrongMat2, 1. * Acts::units::_m, normal);
  BOOST_TEST(particles.size() == 0);

  Acts::Material wrongMat3(
      352.8, 407., 9.012, 4.,
      -1. * Acts::units::_g /
          (Acts::units::_cm * Acts::units::_cm * Acts::units::_cm));
  BOOST_TEST(!g4mis.convertMaterialToG4Stub(wrongMat3));
  particles.clear();
  particles = g4mis(particle, wrongMat3, 1. * Acts::units::_m, normal);
  BOOST_TEST(particles.size() == 0);

  // Stability test for wrong thickness
  particles.clear();
  particles = g4mis(particle, mat, -1. * Acts::units::_m, normal);
  BOOST_TEST(particles.size() == 0);
}
} // namespace Test
} // namespace Fatras
