// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///  Boost include(s)
#define BOOST_TEST_MODULE KinemtaicCast Tests

#include <boost/test/included/unit_test.hpp>
// leave blank line

#include <boost/test/data/test_case.hpp>
// leave blank line

#include <boost/test/output_test_stream.hpp>
// leave blank line

#include "ACTS/Utilities/Units.hpp"
#include "Fatras/Kernel/FatrasDefinitions.hpp"
#include "Fatras/Selectors/KinematicCasts.hpp"

namespace bdata = boost::unit_test::data;
namespace tt    = boost::test_tools;

namespace Fatras {

typedef ParticleInfo Particle;

namespace Test {
    
  double m_pion = 134.9766 * Acts::units::_MeV;    // pi0 rest mass

  // This tests the implementation of kinematic cast operators
  BOOST_AUTO_TEST_CASE(Kinematic_cast_tests)
  {

    Acts::Vector3D position(0.,0.,0.);
    Acts::Vector3D momentum_c(1500.,0.,0);

    // e central electron
    Particle pion_c(position, momentum_c, -1., m_pion);
    Acts::Vector3D position_fwd(0.,0.,100.);    
    Acts::Vector3D momentum_fwd(10., 10., 1500.);  
    
    Particle pion_fwd(position_fwd, momentum_fwd, -1., m_pion);
      
    // the list of possible casts  
    casts::eta    eta_c;
    casts::absEta absEta_c;
    casts::pT     pT_c;
    casts::p      p_c;
    casts::E      E_c;
    casts::vR     vR_c;
    casts::vZ     vZ_c;
      
    // test the central     
    BOOST_TEST(eta_c(pion_c),0.);
    BOOST_TEST(absEta_c(pion_c),0.);
    BOOST_TEST(pT_c(pion_c),1500.);
    BOOST_TEST(p_c(pion_c),1500.);
    BOOST_CHECK(E_c(pion_c) > p_c(pion_c));
    
    BOOST_CHECK_CLOSE(vR_c(pion_c), 0., 10e-5);
    BOOST_CHECK_CLOSE(vZ_c(pion_c), 0., 10e-5);
    
    // test the forward
    BOOST_CHECK(eta_c(pion_fwd) > eta_c(pion_c));
    BOOST_TEST(vZ_c(pion_fwd), 100.);
    
  }
  
} // end of namespace Test
} // end of namespace Fatras
