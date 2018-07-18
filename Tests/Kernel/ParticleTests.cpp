// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///  Boost include(s)
#define BOOST_TEST_MODULE Fatras Definitions Tests

#include <boost/test/included/unit_test.hpp>
// leave blank line

#include <boost/test/data/test_case.hpp>
// leave blank line

#include <boost/test/output_test_stream.hpp>
// leave blank line

#include "Acts/Material/Material.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"
#include "Fatras/Kernel/Particle.hpp"

namespace bdata = boost::unit_test::data;
namespace tt = boost::test_tools;

namespace Fatras {

namespace Test {

// This tests the implementation of the particle info
BOOST_AUTO_TEST_CASE(Particle_tests) {

  /// position at 0.,0.,0
  Acts::Vector3D position{0., 0., 0.};
  // pT of 1 GeV
  Acts::Vector3D momentum{1. * Acts::units::_GeV, 0., 0.};
  // positively charged
  double q = 1.;
  double m = 105.658367 * Acts::units::_MeV; // muon mass

  // create the particle
  Particle particle(position, momentum, m, q, 13, 1);

  // test the energy conservation
  BOOST_CHECK_CLOSE(particle.E, 1.0055663531150525 * Acts::units::_GeV, 10e-5);
  // test the beta factor
  BOOST_CHECK_CLOSE(particle.beta, 0.9944644596571782, 10e-5);
  // test the mass
  BOOST_CHECK_EQUAL(particle.m, m);
  // test the pdg id
  BOOST_CHECK_EQUAL(particle.pdg, 13);

  // test magnitude of momentum
  BOOST_TEST(particle.pT == 1000. * Acts::units::_MeV);
  BOOST_TEST(particle.pT == particle.p);
}

} // namespace Test
} // namespace Fatras
