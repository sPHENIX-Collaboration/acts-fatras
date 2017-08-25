///////////////////////////////////////////////////////////////////
// IHadronicInteractionSampler.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_FATRAS_IHADRONICINTERACTIONSAMPLER_H
#define ACTS_FATRAS_IHADRONICINTERACTIONSAMPLER_H 1

// ACTS include
#include "ACTS/EventData/ParticleDefinitions.hpp"
#include "ACTS/Utilities/Definitions.hpp"

namespace Acts {
class InteractionVertex;
}

namespace Fatras {

/// @class IHadronicInteractionSampler
/// Interface definition for the handling of nuclear/hadronic interactions,
/// to be used by the MC based material effects updater
///
/// @author Andreas Salzburger <Andreas.Salzburger@cern.ch>
/// @author Noemi Calace   <Noemi.Calace@cern.ch>

class IHadronicInteractionSampler {
 public:
  /// Virtual destructor
  virtual ~IHadronicInteractionSampler() = default;

  /// Interface for processing of the presampled nuclear interactions on layer
  /// @param[in] time The time @todo what exactly?
  /// @param[in] position The position
  /// @param[in] momentum The momentum
  /// @param[in] particle The particle type
  /// @return A vector of interaction vertices
  virtual std::vector<Acts::InteractionVertex> doHadronicInteraction(
      double time, const Acts::Vector3D& position,
      const Acts::Vector3D& momentum,
      Acts::ParticleType particle = Acts::pion) const = 0;
};
}

#endif  // ACTS_FATRAS_IHADRONICINTERACTIONSAMPLER_H
