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
#include "Fatras/Selectors/SelectorHelpers.hpp"

namespace bdata = boost::unit_test::data;
namespace tt    = boost::test_tools;

namespace Fatras {

typedef ParticleInfo Particle;
typedef DetectorInfo Detector;


namespace Test {
    
  double m_pion = 134.9766 * Acts::units::_MeV;    // pi0 rest mass

  // This tests the implementation of kinematic cast operators
  BOOST_AUTO_TEST_CASE(SelectorHelper_tests)
  {

    Detector detector;
    
    Acts::Vector3D position(0.,0.,0.);
    Acts::Vector3D momentum_c(1500.,0.,0);

    // e central electron
    Particle pion_c(position, momentum_c, -1., m_pion);
    Acts::Vector3D position_fwd(0.,0.,100.);    
    Acts::Vector3D momentum_fwd(10., 10., 1500.);  
    
    Particle pion_fwd(position_fwd, momentum_fwd, -1., m_pion);


    Acts::Vector3D position_bwd(0.,0.,0.);    
    Acts::Vector3D momentum_bwd(10., 10., -1500.);  
    
    Particle pion_bwd(position_bwd, momentum_bwd, -1., m_pion);
      
    // the list of possible casts  
    casts::eta    eta_c;
    casts::absEta etaAbs_c;
      
    // A minimum of 0.5 Eta is required
    Min<casts::eta> minEta05;
    minEta05.valMin = 0.5;
      
    Min<casts::absEta> minEtaAbs05;
    minEtaAbs05.valMin = 0.5;
    
    // the central will fail both
    BOOST_CHECK(!minEta05(detector,pion_c));
    BOOST_CHECK(!minEtaAbs05(detector,pion_c));
    
    // the forward will satisfy both
    BOOST_CHECK(minEta05(detector,pion_fwd));
    BOOST_CHECK(minEtaAbs05(detector,pion_fwd));
    
    // A maximum of 0.5 Eta is required
    Max<casts::eta> maxEta05;
    maxEta05.valMax = 0.5;
    
    // the central will satisfy both
    BOOST_CHECK(maxEta05(detector,pion_c));
    BOOST_CHECK(maxEta05(detector,pion_c));
    
    // the forward will fail both
    BOOST_CHECK(!maxEta05(detector,pion_fwd));
    BOOST_CHECK(!maxEta05(detector,pion_fwd));
    
    // a range test
    Range<casts::eta> rangeEtaM0;
    rangeEtaM0.valMin = -6.;
    rangeEtaM0.valMax = -0.5;

    Range<casts::absEta> rangeEtaM1;
    rangeEtaM1.valMin = -6.;
    rangeEtaM1.valMax = -0.5;
    
    // the central will fail both
    BOOST_CHECK(!rangeEtaM0(detector,pion_c));
    BOOST_CHECK(!rangeEtaM1(detector,pion_c));
    
    // the forward will fail both
    BOOST_CHECK(!rangeEtaM0(detector,pion_fwd));
    BOOST_CHECK(!rangeEtaM1(detector,pion_fwd));
    
    // the backeard will satsify the eta cast, not the abs(eta)
    BOOST_CHECK(rangeEtaM0(detector,pion_bwd));
    BOOST_CHECK(!rangeEtaM1(detector,pion_bwd));
    
  }
  
} // end of namespace Test
} // end of namespace Fatras
