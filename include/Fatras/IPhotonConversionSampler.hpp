///////////////////////////////////////////////////////////////////
// IPhotonConversionSampler.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_FATRASINTERFACES_IPHOTONCONVERSIONSAMPLER_H
#define ACTS_FATRASINTERFACES_IPHOTONCONVERSIONSAMPLER_H 1

#include <ACTS/Utilities/Definitions.hpp>

namespace Acts {
class InteractionVertex;
}

namespace Fatras {

/// @class IPhotonConversionSampler
/// Interface definition for the handling of photon conversion,
/// to be used by the MaterialInteractionEngine
///
/// @author Sarka Todorova <Sarka.Todorova@cern.ch>
/// @author Noemi Calace   <Noemi.Calace@cern.ch>

class IPhotonConversionSampler {
 public:
  /// Virtual destructor
  virtual ~IPhotonConversionSampler() = default;

  /// Interface for processing of the presampled conversion on layer
  /// @param[in] time @todo what exacteley
  /// @param[in] position The position
  /// @param[in] momentum The momentum
  /// @return A vector of interaction vertices
  virtual std::vector<Acts::InteractionVertex> doConversion(
      double time, const Acts::Vector3D& position,
      const Acts::Vector3D& momentum) const = 0;
};
}

#endif  // ACTS_FATRASINTERFACES_IPHOTONCONVERSIONSAMPLER_H
