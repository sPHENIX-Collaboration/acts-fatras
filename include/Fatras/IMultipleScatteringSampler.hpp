///////////////////////////////////////////////////////////////////
// IMultipleScatteringSampler.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_FATRAS_IMULTIPLESCATTERINGSAMPLER_H
#define ACTS_FATRAS_IMULTIPLESCATTERINGSAMPLER_H 1

// ACTS include
#include "ACTS/EventData/ParticleDefinitions.hpp"

namespace Acts {
class MaterialProperties;
}

namespace Fatras {

/// @class IMultipleScatteringSampler
///
/// Interface class IMultipleScatteringSampler
///
/// @author Andreas Salzburger <Andreas.Salzburger@cern.ch>
/// @author Noemi Calace       <Noemi.Calace@cern.ch>

template <class RandomGenerator>
class IMultipleScatteringSampler {
 public:
  /// Virtual destructor
  virtual ~IMultipleScatteringSampler() = default;

  /// Calculate the theta introduced by multiple scattering
  /// @param[in] mat The material
  /// @param[in] momentum The value of the momentum
  /// @param[in] pathcorrection The correction due to the incident angle
  /// @param[in] particle The particle type
  virtual double simTheta(RandomGenerator& randomGenerator,
                          const Acts::MaterialProperties& mat, double momentum,
                          double pathcorrection,
                          Acts::ParticleType particle = Acts::pion) const = 0;
};
}

#endif  // ACTS_FATRAS_IMULTIPLESCATTERINGSAMPLER_H
