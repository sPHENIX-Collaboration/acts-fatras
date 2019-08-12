// This file is part of the Acts project.
//
// Copyright (C) 2018-2019 Acts project team
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

#include "../Common/Particle.hpp"
#include "Fatras/Kernel/PhysicsList.hpp"
#include "Fatras/Kernel/Process.hpp"
#include "Fatras/Physics/HadronicInteraction/ParametricNuclearInt.hpp"

#include <fstream>
#include <random>

namespace bdata = boost::unit_test::data;
namespace tt = boost::test_tools;

namespace Fatras {

namespace Test {

// the generator
typedef std::mt19937 Generator;

struct {

  double operator()() {
    return (double)m_generator() /
           (double)(m_generator.max() - m_generator.min());
  }

  // standard generator
  Generator m_generator;
} myGenerator;

// some material
const Acts::Material berylium = Acts::Material(
    352.8 * Acts::units::_mm, 421. * Acts::units::_mm, 9.012, 4.,
    1.848 / (Acts::units::_cm * Acts::units::_cm * Acts::units::_cm));

/// The selector
struct Selector {

  /// call operator
  template <typename detector_t, typename particle_t>
  bool operator()(const detector_t, const particle_t &) const {
    return true;
  }
};

/// Test the scattering implementation
BOOST_DATA_TEST_CASE(
    ParametricNuclInt_test_,
    bdata::random((bdata::seed = 23,
                   bdata::distribution = std::uniform_real_distribution<>(
                       0.1 * Acts::units::_GeV, 4. * Acts::units::_GeV))) ^
        bdata::random(
            (bdata::seed = 24,
             bdata::distribution = std::uniform_real_distribution<>(0.5, 2.))) ^
        bdata::xrange(10000),
    p, radLengths, index) {

  const Acts::MaterialProperties detector(berylium, berylium.L0() * radLengths);

  Acts::Vector3D pos{0., 0., 0.};
  Acts::Vector3D mom{0., 0., p};
  double mass = 0.1395701 * Acts::units::_GeV;
  double charge = -1.;
  int pdg = -211;
  int barcode = 0.;
  double time = 0.;
  Particle particle(pos, mom, mass, charge, pdg, barcode, time);

  // Accept everything
  typedef Selector All;
  // Define the processes with selectors
  typedef Process<ParametricNuclearInt, All, All, All> ParamNucIntProcess;

  // now check the EnergyLoss as a PhysicsList
  typedef PhysicsList<ParamNucIntProcess> PhysList;
  PhysList physicsList;

  std::vector<Particle> outgoing;
  physicsList(myGenerator, detector, particle, outgoing);
}
} // namespace Test

} // namespace Fatras
