// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///  Boost include(s)
#define BOOST_TEST_MODULE ParametricNuclearInteraction Tests

#include <boost/test/included/unit_test.hpp>
// leave blank line

#include <boost/test/data/test_case.hpp>
// leave blank line

#include <boost/test/output_test_stream.hpp>
// leave blank line

#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialProperties.hpp"

#include "Fatras/Kernel/Particle.hpp"
#include "Fatras/Physics/HadronicInteraction/ParametricNuclearInt.hpp"
#include "Fatras/Kernel/Process.hpp"
#include "Fatras/Kernel/PhysicsList.hpp"

#include <fstream>
#include <random>

namespace bdata = boost::unit_test::data;
namespace tt = boost::test_tools;

namespace Fatras {

namespace Test {

// the generator
typedef std::mt19937 Generator; // TODO: range?

/// Generator in [0,1]
struct MyGenerator {

MyGenerator(int samen)
{
	generator.seed(samen);
}

// standard generator
Generator generator;

	double operator()()
	{
		return (double) generator() / (double) generator.max();
	}
};

/// The selector
struct MySelector {

  /// call operator
  template <typename detector_t, typename particle_t>
  bool operator()(const detector_t &, const particle_t &) const {
    return true;
  }
};

// some material
Acts::Material berilium = Acts::Material(352.8, 407., 9.012, 4., 1.848e-3);

bool writeOut = true;

std::ofstream ofs("Nuculars.txt", std::ofstream::out | std::ofstream::app);

/// Test the scattering implementation
//~ BOOST_DATA_TEST_CASE(
    //~ ParamNucularInt_test_,
    //~ bdata::random(
        //~ (bdata::seed = 20,
         //~ bdata::distribution = std::uniform_real_distribution<>(0., 1.))) ^
        //~ bdata::random(
            //~ (bdata::seed = 21,
             //~ bdata::distribution = std::uniform_real_distribution<>(0., 1.))) ^
        //~ bdata::random(
            //~ (bdata::seed = 22,
             //~ bdata::distribution = std::uniform_real_distribution<>(0., 1.))) ^
        //~ bdata::random((bdata::seed = 23,
                       //~ bdata::distribution =
                           //~ std::uniform_real_distribution<>(0.5, 10.5))) ^
        //~ bdata::xrange(10000),
    //~ x, y, z, p, index) {
		
BOOST_DATA_TEST_CASE(
    ParamNucularInt_test_, bdata::xrange(10000), index) {
	
	double x = 0., y = 0., z = 1., p = 10.;
MyGenerator mg(index);

  // a detector with 1 mm Be
  Acts::MaterialProperties detector(berilium, 1. * Acts::units::_mm);

  // create the particle and set the momentum
  /// position at 0.,0.,0
  Acts::Vector3D position{0., 0., 0.};
  // p of 1 GeV
  Acts::Vector3D momentum =
      p * Acts::units::_GeV * Acts::Vector3D(x, y, z).unit();
  // positively charged
  double q = -1.;
  double m = 134.9766 * Acts::units::_MeV; // pion mass

  // create the particle
  Particle particle(position, momentum, m, q, -211, 1);

ParametricNuclearInt paramNuclInt;

std::vector<Particle> par = paramNuclInt(mg, detector, particle);
//~ ofs << index << "\t" << par.size() << std::endl;
//~ for(size_t i = 0; i < par.size(); i++)
	//~ ofs << par[i].pdg << "\t" << par[i].m << "\t" << par[i].q << "\t" << par[i].E << "\t" 
		//~ << par[i].p << "\t" << par[i].momentum.x() << "\t" << par[i].momentum.y() << "\t" 
		//~ << par[i].momentum.z() << std::endl;


  typedef MySelector All;
  std::vector<Particle> outgoing;
  typedef Process<ParametricNuclearInt, All, All, All> HadronProcess;
  PhysicsList<HadronProcess> hsPhysicsList;
  hsPhysicsList(mg, detector, particle, outgoing);
  
}

} // namespace Test

} // namespace Fatras

