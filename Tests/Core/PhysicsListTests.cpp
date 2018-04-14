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

#include "Fatras/Kernel/PhysicsList.hpp"

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
  
  /// Physics process that does not trigger a break
  struct SterileProcess {

    /// call operator 
    template <typename cache_t,
              typename generator_t,
              typename detector_t,
              typename particle_t>
    bool
    operator()(cache_t&, 
               generator_t&,
               const detector_t&,
               const particle_t&,             
               std::vector<particle_t>&) const 
    { return false; }
  
  };
  
  
  /// Physics process that DOES trigger a break
  struct FatalProcess {

    /// call operator 
    template <typename cache_t,
              typename generator_t,
              typename detector_t,
              typename particle_t>
    bool
    operator()(cache_t&, 
               generator_t&,
               const detector_t&,
               const particle_t&,             
               std::vector<particle_t>&) const 
    { return true; }
  
  };
  
  // This tests the implementation of the AbortList
  // and the standard aborters
  BOOST_AUTO_TEST_CASE(PhysicsLists_test)
  {
    //
    Cache_type                 cache;
    Generator_type             generator;
    Detector_type               detector;
    Particle_type              in;
    std::vector<Particle_type> out;
    
    /// empty physics_list
    typedef PhysicsList<> ProcessLess;
    ProcessLess emptyList;
    
    /// sterile test should never send the abort command
    BOOST_TEST(!emptyList(cache,generator,detector,in,out));
    
    // now create a single sterile process
    typedef PhysicsList<SterileProcess> SterileList;
    SterileList sterileProcess;
    
    /// sterile test should not send the abort command
    BOOST_TEST(!sterileProcess(cache,generator,detector,in,out));
    
    // now create a single fatal process
    typedef PhysicsList<FatalProcess> FatalList;
    FatalList fatalProcess;
    
    /// fatal test should send  abort command
    BOOST_TEST(fatalProcess(cache,generator,detector,in,out));
    
    // now create a list of a sterile and fatal process
    typedef PhysicsList<SterileProcess,FatalProcess> SterileFatalList;
    SterileFatalList stfaProcess;
    
    /// fatal test should send  abort command
    BOOST_TEST(stfaProcess(cache,generator,detector,in,out));
    
  }

}  // namespace Test
}  // namespace Fatras
