// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///  Boost include(s)
#define BOOST_TEST_MODULE PdgSelectors Tests

#include <boost/test/included/unit_test.hpp>
// leave blank line

#include <boost/test/data/test_case.hpp>
// leave blank line

#include <boost/test/output_test_stream.hpp>
// leave blank line

#include "Fatras/Kernel/FatrasDefinitions.hpp"
#include "Fatras/Kernel/SelectorList.hpp"
#include "Fatras/Selectors/PdgSelectors.hpp"

namespace bdata = boost::unit_test::data;
namespace tt    = boost::test_tools;

namespace Fatras {

typedef DetectorInfo Detector;
typedef ParticleInfo Particle;

namespace Test {
  
  Detector detector;
  
  double m_muon = 0.51099891 * Acts::units::_MeV;  // electron mass
  double m_e    = 105.658367 * Acts::units::_MeV;  // muon mass
  double m_pion = 134.9766 * Acts::units::_MeV;    // pi0 rest mass

  // This tests the implementation of the pdg selcetors
  BOOST_AUTO_TEST_CASE(PdgSelectors_test)
  {

    Acts::Vector3D position(0.,0.,0.);
    Acts::Vector3D momentum(1500.,0.,0);

    Particle electron(position, momentum, -1., m_e);
    electron.pdg = 11;
    
    Particle positron(position, momentum, -1., m_e);
    positron.pdg = -11;
    
    Particle muon(position, momentum, -1., m_muon);
    muon.pdg     = 13;

    Particle antimuon(position, momentum, 1., m_muon);
    antimuon.pdg = -13;
    
    AbsPdgSelector<11> epSelector;
    
    BOOST_CHECK(epSelector(detector,electron));
    BOOST_CHECK(epSelector(detector,positron));
    BOOST_CHECK(!epSelector(detector,muon));
    BOOST_CHECK(!epSelector(detector,antimuon));
    
    PdgSelector<13> muSelector;

    BOOST_CHECK(!muSelector(detector,electron));
    BOOST_CHECK(!muSelector(detector,positron));
    BOOST_CHECK(muSelector(detector,muon));
    BOOST_CHECK(!muSelector(detector,antimuon));
    
    AbsPdgExcluder<11> epExcluder;
    
    BOOST_CHECK(!epExcluder(detector,electron));
    BOOST_CHECK(!epExcluder(detector,positron));
    BOOST_CHECK(epExcluder(detector,muon));
    BOOST_CHECK(epExcluder(detector,antimuon));
    
    PdgExcluder<13> muExcluder;

    BOOST_CHECK(muExcluder(detector,electron));
    BOOST_CHECK(muExcluder(detector,positron));
    BOOST_CHECK(!muExcluder(detector,muon));
    BOOST_CHECK(muExcluder(detector,antimuon));
    
    SelectorListOR<PdgSelector<13>,PdgSelector<11>> emuSelection;
    
    BOOST_CHECK(emuSelection(detector,electron));
    BOOST_CHECK(emuSelection(detector,muon));
    BOOST_CHECK(!emuSelection(detector,positron));
    BOOST_CHECK(!emuSelection(detector,antimuon));
    
    
  }
  
  
} // end of namespace Test

} // end of namespace Fatras

