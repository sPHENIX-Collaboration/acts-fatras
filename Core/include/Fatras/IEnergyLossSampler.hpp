///////////////////////////////////////////////////////////////////
// IEnergyLossSampler.h, Acts project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_FATRAS_IENERGYLOSSSAMPLER_H
#define ACTS_FATRAS_IENERGYLOSSSAMPLER_H 1

#include <ACTS/EventData/ParticleDefinitions.hpp>
#include <ACTS/Utilities/Definitions.hpp>

namespace Acts {
class MaterialProperties;
}

namespace Fatras {

class EnergyLoss;

/// @class IEnergyLossSampler
///
/// @brief Interface class for energy loss samplers
///
/// @author Andreas Salzburger <Andreas.Salzburger@cern.ch>
/// @author Noemi Calace       <Noemi.Calace@cern.ch>

template <class RandomGenerator>
class IEnergyLossSampler
{
public:
  /// Virtual destructor
  virtual ~IEnergyLossSampler() = default;

  /// Calculation of energy loss along a certain path
  /// Using dEdX and integrating along pathlength
  /// - The sign depends on the given propagation direction
  /// - Units: [MeV]
  /// @note assuming constant dEdX during for the path
  /// @param[in] mat The material properties, containing the thickness of the
  /// material
  /// @param[in] momentum The momentum of the paricle
  /// @param[in] pathcorrection Correction to incident angle
  /// @param[in] dir The direction of the propagation
  /// @param[in] particle The particle type
  /// @return The energy loss
  virtual EnergyLoss
  energyLoss(RandomGenerator&                randomGenerator,
             const Acts::MaterialProperties& mat,
             double                          momentum,
             double                          pathcorrection,
             Acts::NavigationDirection             dir      = Acts::forward,
             Acts::ParticleType              particle = Acts::pion) const = 0;

  /// dEdX calculation when providing MaterialProperties,
  /// a momentum, and a ParicleHypothesis
  /// - Units: [Mev/mm]
  /// @param mat The material properties, containing the thickness of the
  /// material
  /// @param momentum The momentum of the paricle
  /// @param particle The particle type
  /// @return The energy loss per length unit
  /*  virtual double dEdX(const Acts::MaterialProperties& mat, double momentum,
                        Acts::ParticleType particle = Acts::pion) const = 0;*/
};

}  // end of namespace

#endif  // ACTS_FATRAS_IENERGYLOSSSAMPLER_H
