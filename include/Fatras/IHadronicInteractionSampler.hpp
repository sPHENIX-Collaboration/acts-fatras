///////////////////////////////////////////////////////////////////
// IHadronicInteractionSampler.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_FATRAS_IHADRONICINTERACTIONSAMPLER_H
#define ACTS_FATRAS_IHADRONICINTERACTIONSAMPLER_H 1

#include <Acts/EventData/ParticleDefinitions.hpp>
#include <Acts/Utilities/Definitions.hpp>

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

template <class RandomGenerator>
class IHadronicInteractionSampler
{
public:
  /// Virtual destructor
  virtual ~IHadronicInteractionSampler() = default;

  /// Interface for processing of the presampled nuclear interactions on layer
  /// @paramT[in] randomGenerator - the provided generator
  /// @param[in] time the time of the interaction
  /// @param[in] position The position
  /// @param[in] momentum The momentum
  /// @param[in] particle The particle type
  /// @return A vector of process vertices vertices ( in case of mulitple )
  virtual std::vector<Acts::ProcessVertex>
  doHadronicInteraction(RandomGenerator&      randomGenerator,
                        double                time,
                        const Acts::Vector3D& position,
                        const Acts::Vector3D& momentum,
                        Acts::ParticleType    particle = Acts::pion) const = 0;
};
}

#endif  // ACTS_FATRAS_IHADRONICINTERACTIONSAMPLER_H
