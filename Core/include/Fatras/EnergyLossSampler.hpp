///////////////////////////////////////////////////////////////////
// EnergyLossSampler.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_FATRAS_ENERGYLOSSSAMPLER_H
#define ACTS_FATRAS_ENERGYLOSSSAMPLER_H 1

#include <memory>

#include <ACTS/EventData/ParticleDefinitions.hpp>
#include <ACTS/Material/MaterialProperties.hpp>
#include <ACTS/Utilities/Definitions.hpp>
#include <ACTS/Utilities/Logger.hpp>

#include "Fatras/IEnergyLossSampler.hpp"

namespace Fatras {

class EnergyLoss;

/// @class EnergyLossSampler
///
/// Sampler calculating the energy loss
///
/// Sampler for a eloss of a track. It uses the Bethe-Bloch calculation
/// and extends the IEnergyLossSampler interface.
///
/// @author Andreas Salzburger <Andreas.Salzburger@cern.ch>
/// @author Noemi Calace       <Noemi.Calace@cern.ch>

template <class RandomGenerator>
class EnergyLossSampler : virtual public IEnergyLossSampler<RandomGenerator>
{
public:
  /// @struct Config
  /// Configuration config for the EnergyLossSampler
  struct Config
  {
    /// Scalor for the mean energy loss
    double scalorMOP = 1.;  // G4: 0.962
    /// Scalor for the energy loss uncertainty
    double scalorSigma = 1.;  // G4: 0.6755

    Config() {}
  };

  /// Constructor with AlgTool parameters
  /// @param[in] config Configuration
  /// @param[in] logger The logging instance
  EnergyLossSampler(const Config&                       config,
                    std::unique_ptr<const Acts::Logger> logger
                    = Acts::getDefaultLogger("EnergyLossSampler",
                                             Acts::Logging::INFO));

  /// Destructor
  ~EnergyLossSampler() = default;

  // clang-format off
  /// @copydoc IEnergyLossSampler::energyLoss(const Acts::MaterialProperties&,double,double,Acts::NavigationDirection,Acts::ParticleType) const
  // clang-format on
  EnergyLoss
  energyLoss(RandomGenerator&                randomGenerator,
             const Acts::MaterialProperties& mat,
             double                          momentum,
             double                          pathcorrection,
             Acts::NavigationDirection       dir      = Acts::forward,
             Acts::ParticleType              particle = Acts::pion) const final;

  // @todo Merge with functions in ACTS core - only one bethe and one heitler
  // formula!
  // clang-format off
  /// @copydoc IEnergyLossSampler::energyLoss(const Acts::MaterialProperties&,double,Acts::ParticleType particleHypothesis) const
  // clang-format on
  /*  double dEdX(const Acts::MaterialProperties& materialProperties,
                double momentum,
                Acts::ParticleType particleHypothesis = Acts::pion) const
     final;*/

  /// Set configuration method
  /// @param[in] config The configuration object
  void
  setConfiguration(const Config& config);

  /// Get configuration method
  Config
  getConfiguration() const
  {
    return m_config;
  }

  /// Set logging instance
  ///
  /// @param[in] logger the logging instance to be set
  void
  setLogger(std::unique_ptr<const Acts::Logger> logger);

protected:
  /// Configuration object
  Config m_config;
  /// Struct of Particle masses
  Acts::ParticleMasses m_particleMasses;

private:
  /// Apply dEdX according to Bethe-Bloch
  /// @param mat The material properties, containing the thickness of the
  /// material
  /// @param gamma The lorentz factor
  /// @param beta The ratio velocity to speed of light
  /// @param particle The particle type
  /// @return The energy loss per length unit
  /* double dEdXBetheBloch(const Acts::MaterialProperties& mat, double gamma,
                         double beta,
                         Acts::ParticleType particle = Acts::pion) const;

   /// Apply dEdX according to Bethe-Heitler
   /// @param mat The material properties, containing the thickness of the
   /// material
   /// @param initialE The initial particle energy
   /// @param particle The particle type
   /// @return The energy loss per length unit
   double dEdXBetheHeitler(const Acts::MaterialProperties& mat, double initialE,
                           Acts::ParticleType particle = Acts::pion) const;
 */
  const Acts::Logger&
  logger() const
  {
    return *m_logger;
  }
  /// logger instance
  std::unique_ptr<const Acts::Logger> m_logger;
};
}

/// define the templated function
#include "detail/EnergyLossSampler.ipp"

#endif  //  ACTS_FATRAS_ENERGYLOSSSAMPLER_H
