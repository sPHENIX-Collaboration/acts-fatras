// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///  Boost include(s)
#define BOOST_TEST_MODULE Geant4MaterialInteraction Tests

#include <boost/test/included/unit_test.hpp>
// leave blank line

#include <boost/test/data/test_case.hpp>
// leave blank line

#include <boost/test/output_test_stream.hpp>
// leave blank line

#include "Geant4MaterialInteractionStub.hpp"

namespace bdata = boost::unit_test::data;
namespace tt = boost::test_tools;

namespace Fatras {
namespace Test {

BOOST_DATA_TEST_CASE(
    ParamNucularInt_test_,
        bdata::random((bdata::seed = 21,
                       bdata::distribution =
                           std::uniform_real_distribution<>(0.01 * 394.133 / 10., 2. * 394.133 / 10.))) ^
        bdata::random((bdata::seed = 22,
                       bdata::distribution =
                           std::uniform_real_distribution<>(0.5, 20.))) ^
        bdata::xrange(1000),
    detectorThickness, p, index) {

	Geant4MaterialInteractionStub g4mi;

}
} // namespace Test
} // namespace Fatras

