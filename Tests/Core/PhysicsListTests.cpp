// This file is part of the ACTS project.
//
// Copyright (C) 2017-2018 ACTS project team
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

#include "Fatras/PhysicsList.hpp"

namespace bdata = boost::unit_test::data;
namespace tt    = boost::test_tools;

namespace Fatras {

namespace Test {

  /// needed are : cache, generator, detector, particle
  struct Cache_type
  {
  };
  
  struct Generator_type 
  {
  };
  
  struct Detector_type
  {
  };
  
  struct Particle_type
  {
  };
  
  
  // This tests the implementation of the AbortList
  // and the standard aborters
  BOOST_AUTO_TEST_CASE(PhysicsLists_test)
  {
    //
    Cache_type                 cache;
    Generator_type             generator;
    Detctor_type               detector;
    Particle_type              incoming;
    std::vector<Particle_tyep> outoing;
    
    /// empty physics_list
    typedef PhysicsList<> ProcessLess;
    ProcessLess sterile;
    BOOST_TEST(sterile(cache,generator,detector,incoming,outgoing));
    
  }

}  // namespace Test
}  // namespace Fatras
