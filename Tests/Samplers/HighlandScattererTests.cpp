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

#include <random>
#include "Fatras/Kernel/FatrasDefinitions.hpp"
#include "Fatras/Samplers/Scattering/HighlandScatterer.hpp"

namespace bdata = boost::unit_test::data;
namespace tt    = boost::test_tools;

namespace Fatras {

namespace Test {
  
  // let's get some material
  Material berilium = Material(352.8, 407., 9.012, 4., 1.848e-3);
  
  /// Test the scattering implementation
  BOOST_DATA_TEST_CASE(
      HighlandScattering_test_,
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
                           = std::uniform_int_distribution<>(0.5, 1.5)))
          ^ bdata::xrange(100),
      x,
      y,
      z,
      p,
      index)
  {

    /// the generator
    typedef std::mt19937 Generator;
    typedef ParticleInfo Particle;
    typedef DetectorInfo Detector;
    
    // standard generator
    Generator generator;
    
    // a detector with 1 mm Be
    Detector detector;
    detector.material = berilium;
    detector.thickness = 1 * Acts::units::_mm;

    // create the particle and set the momentum
    /// position at 0.,0.,0
    Acts::Vector3D position{0.,0.,0.};
    // pT of 1 GeV 
    Acts::Vector3D momentum = p * Acts::units::_GeV * Acts::Vector3D(x,y,z).unit();
    // positively charged
    double q = -1.;
    double m = 105.658367 * Acts::units::_MeV;  // muon mass
    
    // create the particle 
    Particle particle(position,momentum,q,m,13,1);
    
    // make the highland scatterer
    HighlandScatterer hscat;
    
    double angle = hscat(generator, detector, particle);
    
    
    
  }
  
  
} // end of namespace Test

} // end of namespace Fatras
