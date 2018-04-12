///////////////////////////////////////////////////////////////////
// HadronicInteractionParametricSampler.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_FATRASTOOLS_HADRONICINTERACTIONPARAMETRICSAMPLER_H
#define ACTS_FATRASTOOLS_HADRONICINTERACTIONPARAMETRICSAMPLER_H 1

// FATRAS
#include <memory>

#include "Acts/EventData/ParticleDefinitions.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Fatras/IHadronicInteractionSampler.hpp"

namespace Fatras {

///  @class HadronicInteractionParametricSampler
///
/// Parametric implementation of nuclear interactions to be used
/// in Fatras. The parameterisation is gathered from Geant4.
///
/// @author Andreas Salzburger <Andreas.Salzburger@cern.ch>
/// @author Carsten Magass     <Carsten.Magass@cern.ch>
/// @author Noemi Calace       <Noemi.Calace@cern.ch>
///

template <class RandomGenerator>
class HadronicInteractionParametricSampler
    : virtual public IHadronicInteractionSampler<RandomGenerator>
{
public:
  ///  @struct Config
  ///  Configuration of this Samples
  struct Config
  {
    int    processCode;          //!< process code */
    double minimumHadOutEnergy;  //!< hadronic interaction setting */
    bool   cutChain;             //!<

    Config()
      : processCode(1)
      , minimumHadOutEnergy(100. * Acts::units::_MeV)
      , cutChain(true)
    {
    }
  };

  /// Constructor
  /// @param [in] hiConfig the configuration object
  /// @param [in] logger the logger instance
  HadronicInteractionParametricSampler(
      const Config&                       hiConfig,
      std::unique_ptr<const Acts::Logger> logger
      = Acts::getDefaultLogger("HadronicInteractionSampler",
                               Acts::Logging::INFO));

  /// Destructor
  virtual ~HadronicInteractionParametricSampler();

  /// processing of the presampled nuclear interactions on layer
  /// This method returns the particle's children
  /// @paramT [in] the random generator
  /// @param [in] time is the time of the interaction
  /// @param [in] position of the interaction vertex
  /// @param [in] momentum is the momentum of the incoming particle
  /// @param [in] particle type - nomen est omen
  /// @return of the ProcessVertex
  std::vector<Acts::ProcessVertex>
  doHadronicInteraction(RandomGenerator&      randomGenerator,
                        double                time,
                        const Acts::Vector3D& position,
                        const Acts::Vector3D& momentum,
                        Acts::ParticleType particle = Acts::pion) const final;

  /// Set configuration method
  void
  setConfiguration(const Config& eeConfig);

  /// Set logging instance
  ///
  /// @param logger the logging instance to be set
  void
  setLogger(std::unique_ptr<const Acts::Logger> logger);

protected:
  Config m_config;  //!< configuration object
  /// Struct of Particle masses
  Acts::ParticleMasses m_particleMasses;

private:
  /// Collect secondaries for layer material update
  std::vector<Acts::ProcessVertex>
  getHadState(RandomGenerator&      randomGenerator,
              double                time,
              double                p,
              const Acts::Vector3D& vertex,
              const Acts::Vector3D& particleDir,
              Acts::ParticleType    particle) const;

  const Acts::Logger&
  logger() const
  {
    return *m_logger;
  }

  /// logger instance
  std::unique_ptr<const Acts::Logger> m_logger;
};

}  // end of namespace

/// Define the templated function
#include "detail/HadronicInteractionParametricSampler.ipp"

#endif  // ACTS_FATRASTOOLS_HADRONICINTERACTIONPARAMETRICSAMPLER_H
