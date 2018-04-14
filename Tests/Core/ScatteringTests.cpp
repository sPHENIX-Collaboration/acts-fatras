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

#include <algorithm>
#include "Fatras/Kernel/SelectorList.hpp"
#include "Fatras/Kernel/PhysicsList.hpp"
#include "Fatras/Processes/Scattering.hpp"
#include "ACTS/Utilities/Definitions.hpp"

namespace bdata = boost::unit_test::data;
namespace tt    = boost::test_tools;

namespace Fatras {

namespace Test {

  /// the generator
  typedef std::mt19937 Generator;

  /// The particle definition
  struct Particle
  {
    Acts::Vector3D momentum;
    double         p = 0.;
    
    Particle(const Acts::Vector3D& pvec) :
      momentum(pvec),
      p(pvec.mag())
    {}
    
  };
  
  /// The selector 
  struct Selector {
    
    /// call operator 
    template <typename particle_t>
    bool
    operator()(const particle_t&) const 
    { return true; }
  
  };
  
  /// The detector
  struct Detector {
    
  };
  
  /// The scattering formula
  struct ScatteringAngle {

    double cvalue = 0.15;

    template <typename generator_t,
              typename detector_t,
              typename particle_t>
    double 
    operator()(generator_t&,
               const detector_t&,
               const particle_t&) const 
    {
      return cvalue;
    }
    
  };
  
  /// 
  BOOST_DATA_TEST_CASE(
      Scattering_test_,
      bdata::random((bdata::seed = 20,
                     bdata::distribution
                     = std::uniform_real_distribution<>(0.,1.)))
          ^ bdata::random((bdata::seed = 21,
                           bdata::distribution
                           = std::uniform_real_distribution<>(0.,1.)))
          ^ bdata::random((bdata::seed = 22,
                           bdata::distribution
                           = std::uniform_real_distribution<>(0.,1.)))
          ^ bdata::random((bdata::seed = 23,
                           bdata::distribution
                           = std::uniform_int_distribution<>(0, 1)))
          ^ bdata::xrange(100),
      x,
      y,
      z,
      p,
      index)
  {

    // standard generator
    Generator generator;
    
    // Dummy detctor
    Detector detector;

    // create the particle and set the momentum
    Acts::Vector3D startMomentum = p * Acts::Vector3D(x,y,z).unit();
    Particle particle(startMomentum);
      
    // outgoing particles (always none for scattering)
    std::vector<Particle> outgoing;
    
    // The select all list
    typedef SelectorList<Selector> SelectAll;
    Scattering<ScatteringAngle, SelectAll> cScattering;
    
    // scattering is not allowed to throw abort command
    BOOST_CHECK(!cScattering(generator,detector,particle,outgoing));
    
    // check the the particle momentum magnitude is identical
    BOOST_CHECK_CLOSE(particle.momentum.mag(),startMomentum.mag(),10e-3);
    
    
  }

}  // namespace Test
}  // namespace Fatras
