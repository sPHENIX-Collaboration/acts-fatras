///////////////////////////////////////////////////////////////////
// FatrasDefinitions.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_FATRAS_FATRASDEFINITIONS_H
#define ACTS_FATRAS_FATRASDEFINITIONS_H

#include <cmath>

#include <ACTS/EventData/ParticleDefinitions.hpp>
#include <ACTS/Utilities/Units.hpp>

namespace Fatras {
namespace constants {

  // @todo multiply with units

  /// KOverA factor in Bethe-Bloch equation [MeV*cm2/gram]
  constexpr double ka_BetheBloch = 30.7075;

  /// Fine structure constexprant
  constexpr double alpha = 1. / 137.;

  /// Multiple scattering paramters
  constexpr double main_RutherfordScott = 13.6;
  constexpr double log_RutherfordScott  = 0.038;

  constexpr double main_RossiGreisen = 17.5;
  constexpr double log_RossiGreisen  = 0.125;

}  // namespace constants
}  // namespace Fatras

#endif
