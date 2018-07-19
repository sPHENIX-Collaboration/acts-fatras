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
//~ #include "Fatras/Kernel/PhysicsList.hpp"

#include <fstream>
#include <random>

namespace bdata = boost::unit_test::data;
namespace tt = boost::test_tools;

namespace Fatras {

namespace Test {

// the generator
typedef std::mt19937 Generator; // TODO: range?

/// Generator [0,1]
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

// some material
Acts::Material berilium = Acts::Material(352.8, 407., 9.012, 4., 1.848e-3);

//~ bool write_csv = true;

//~ std::ofstream os("ScatteringAngles.csv",
                 //~ std::ofstream::out | std::ofstream::trunc);

/// Test the scattering implementation
BOOST_DATA_TEST_CASE(
    ParamNucularInt_test_,
    bdata::random(
        (bdata::seed = 20,
         bdata::distribution = std::uniform_real_distribution<>(0., 1.))) ^
        bdata::random(
            (bdata::seed = 21,
             bdata::distribution = std::uniform_real_distribution<>(0., 1.))) ^
        bdata::random(
            (bdata::seed = 22,
             bdata::distribution = std::uniform_real_distribution<>(0., 1.))) ^
        bdata::random((bdata::seed = 23,
                       bdata::distribution =
                           std::uniform_real_distribution<>(0.5, 10.5))) ^
        bdata::xrange(100),
    x, y, z, p, index) {
		
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
  Particle particle(position, momentum, q, m, -211, 1);

ParametricNuclearInt::Config cfg;
cfg.m_hadronInteractionFromX0 = false; // TODO: this flag can also be used with true;
cfg.m_hadronInteractionProbabilityScale = 1.;
cfg.MAXHADINTCHILDREN = 100000;
cfg.m_minimumHadOutEnergy = 0.;
ParametricNuclearInt paramNuclInt(cfg);

std::vector<Particle> par = paramNuclInt(mg, detector, particle);
std::cout << "#particles: " << par.size() << std::endl;
for(size_t i = 0; i < par.size(); i++)
	std::cout << par[i].pdg << "\t" << par[i].E << std::endl;
	


// TODO: Process and PhysicsList test
}

} // namespace Test

} // namespace Fatras

