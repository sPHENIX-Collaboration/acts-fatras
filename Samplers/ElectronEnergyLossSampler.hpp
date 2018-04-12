///////////////////////////////////////////////////////////////////
// ElectronEnergyLossSampler.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_FATRSTOOLS_ELECTRONENERGYLOSSSAMPLER_H
#define ACTS_FATRSTOOLS_ELECTRONENERGYLOSSSAMPLER_H 1

#include <memory>

#include "../include/Fatras/RandomNumberDistributions.hpp"
#include "ACTS/EventData/ParticleDefinitions.hpp"
#include "ACTS/Extrapolation/detail/MaterialInteraction.hpp"
#include "ACTS/Material/MaterialProperties.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/Logger.hpp"
#include "Fatras/IEnergyLossSampler.hpp"

namespace Fatras {

class EnergyLoss;

/** @class ElectronEnergyLossSampler
 *
 * Sampler for a eloss of a track
 * It uses the Bethe-Heitler calculation for electrons
 * it extends the IElectronEnergyLossSampler interface
 *
 * @author Andreas Salzburger <Andreas.Salzburger@cern.ch>
 * @author Noemi Calace       <Noemi.Calace@cern.ch>
 */

class ElectronEnergyLossSampler : virtual public IEnergyLossSampler
{
public:
  /** @struct Config
    Configuration config for the EnergyLossSampler
    */
  struct Config
  {
    /** Random Generator service  */
    std::shared_ptr<IRandomNumbers> randomNumbers;
    /** the one free parameter to scale */
    double scaleFactor;

    Config() : randomNumbers(nullptr), scaleFactor(1.) {}
  };

  /** Constructor with AlgTool parameters */
  ElectronEnergyLossSampler(
      const Config&                 elConfig,
      std::unique_ptr<Acts::Logger> logger
      = Acts::getDefaultLogger("ElectronEnergyLossSampler",
                               Acts::Logging::INFO));

  /** Destructor */
  ~ElectronEnergyLossSampler();

  /** IEnergyLossSampler public method to compute dEdX */
  double
  dEdX(const Acts::MaterialProperties& materialProperties,
       double                          momentum,
       Acts::ParticleType particleHypothesis = Acts::pion) const final;

  /** IEnergyLossSampler public method to compute the mean and variance of the
   * energy loss */
  Fatras::EnergyLoss
  energyLoss(const Acts::MaterialProperties& mat,
             double                          momentum,
             double                          pathcorrection,
             Acts::NavigationDirection             dir      = Acts::forward,
             Acts::ParticleType              particle = Acts::pion) const final;

  /** Set configuration method */
  void
  setConfiguration(const Config& eeConfig);

  /** Get configuration method */
  Config
  getConfiguration() const;

  /// Set logging instance
  ///
  /// @param logger the logging instance to be set
  void
  setLogger(std::unique_ptr<Acts::Logger> logger);

protected:
  Config m_config;  //!< configuration object

private:
  const Acts::Logger&
  logger() const
  {
    return *m_logger;
  }
  /// logger instance
  std::unique_ptr<Acts::Logger> m_logger;

private:
  /** Private method to compute the Bethe-Heitler PDF */
  std::vector<double>
  betheHeitlerPDF(double pathLength) const;
};

inline double
ElectronEnergyLossSampler::dEdX(const Acts::MaterialProperties&,
                                double,
                                Acts::ParticleType) const
{
  return 0;
}

/** Return the configuration object */
inline ElectronEnergyLossSampler::Config
ElectronEnergyLossSampler::getConfiguration() const
{
  return m_config;
}
}
#endif  //  ACTS_FATRSTOOLS_ELECTRONENERGYLOSSSAMPLER_H
