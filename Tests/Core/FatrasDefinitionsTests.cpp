// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///  Boost include(s)
#define BOOST_TEST_MODULE AbortList Tests

#include <boost/test/included/unit_test.hpp>
// leave blank line

#include <boost/test/data/test_case.hpp>
// leave blank line

#include <boost/test/output_test_stream.hpp>
// leave blank line

#include "Fatras/Kernel/FatrasDefinitions.hpp"
#include "ACTS/Material/Material.hpp"
#include "ACTS/Utilities/Units.hpp"
#include "ACTS/Utilities/Definitions.hpp"

namespace bdata = boost::unit_test::data;
namespace tt    = boost::test_tools;

namespace Fatras {

namespace Test {
  
  // This tests the implementation of the detector info 
  BOOST_AUTO_TEST_CASE(DetectorInfo_tests)
  {
    
    Acts::Material mat(1.,2.,3.,4.,5.);
    double pathlength = 0.5;
  
    DetectorInfo detector;
    detector.material   = mat;
    detector.pathLength = pathlength;
    
    // This is a simple container 
    BOOST_TEST(detector.material == mat);
    BOOST_TEST(detector.pathLength == pathlength);    
  }

  // This tests the implementation of the particle info 
  BOOST_AUTO_TEST_CASE(ParticleInfo_tests)
  {
    
    /// position at 0.,0.,0
    Acts::Vector3D position{0.,0.,0.};
    // pT of 1 GeV 
    Acts::Vector3D momentum{1.*Acts::units::_GeV,0.,0.};
    // positively charged
    double q = 1.;
    double m = 105.658367 * Acts::units::_MeV;  // muon mass
    
    // create the particle 
    ParticleInfo particle(position,momentum,q,m,13,1);

    // test the energy conservation
    BOOST_CHECK_CLOSE(particle.E,1.0055663531150525*Acts::units::_GeV,10e-5);
    // test the beta factof
    BOOST_CHECK_CLOSE(particle.beta, 0.9944644596571782, 10e-5);
    
    // test magnitude of momentum
    BOOST_TEST(particle.pT == 1000.*Acts::units::_MeV);
    BOOST_TEST(particle.pT == particle.p);  
  }
 
}  // namespace Test
}  // namespace Fatras
