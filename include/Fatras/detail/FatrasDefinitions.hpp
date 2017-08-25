///////////////////////////////////////////////////////////////////
// FatrasDefinitions.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_FATRAS_FATRASDEFINITIONS_H
#define ACTS_FATRAS_FATRASDEFINITIONS_H

#include <cmath>
#include "ACTS/EventData/ParticleDefinitions.hpp"
#include "ACTS/Utilities/Units.hpp"

namespace Fatras {

namespace constants {

// @todo multiply with units

/// KOverA factor in Bethe-Bloch equation [MeV*cm2/gram]
const double ka_BetheBloch = 30.7075;

/// Fine structure constant
const double alpha = 1. / 137.;

/// Multiple scattering paramters
const double main_RutherfordScott = 13.6;
const double log_RutherfordScott = 0.038;

const double main_RossiGreisen = 17.5;
const double log_RossiGreisen = 0.125;

// statics doubles used for calculations
const double sqrtTwo = sqrt(2.);
const double oneOverThree = 1. / 3.;
}
}

#endif
