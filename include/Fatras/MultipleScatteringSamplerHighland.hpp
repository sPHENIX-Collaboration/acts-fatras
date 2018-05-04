///////////////////////////////////////////////////////////////////
// MultipleScatteringSamplerHighland.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_FATRASTOOLS_MULTIPLESCATTERINGSAMPLERHIGHLAND_H
#define ACTS_FATRASTOOLS_MULTIPLESCATTERINGSAMPLERHIGHLAND_H 1

#include <memory>

#include <Acts/EventData/ParticleDefinitions.hpp>
#include <Acts/Material/MaterialProperties.hpp>
#include <Acts/Utilities/Definitions.hpp>

#include "Fatras/IMultipleScatteringSampler.hpp"

namespace Fatras {

/// @class MultipleScatteringSamplerHighland
///
/// The Formula used is the highland formula for the projected scattering angle:
///
/// @f$ \theta_{ms} = \frac{13.6MeV}{p}\cdot\sqrt{t/X_{0}}[1 +
/// 0.038\ln(t/X_{0})]
/// @f$
///
/// What is returned is the square of the expectation value of the deflection
/// @f$ < (\theta_ms)^2 > = \sigma_ms^2 @f$
///
/// @author Andreas Salzburger <Andreas.Salzburger@cern.ch>
/// @author Noemi Calace       <Noemi.Calace@cern.ch>
///
template <class RandomGenerator>
class MultipleScatteringSamplerHighland
    : virtual public IMultipleScatteringSampler<RandomGenerator>
{
public:
  /// @struct Config
  /// Configuration Struct for the MultipleScattering Sampler
  struct Config
  {
    /// Include the log term
    bool log_include;

    Config() : log_include(true) {}
  };

  /// Constructor
  /// @param[in] config The configuration object
  MultipleScatteringSamplerHighland(const Config& config);

  /// Virtual destructor
  virtual ~MultipleScatteringSamplerHighland() = default;

  /// Implementation according to the RutherFord-Scott Formula
  // clang-format off
  /// @copydoc IMultipleScatteringSampler::simTheta(const Acts::MaterialProperties&,double,double,Acts::ParticleType) const
  // clang-format on
  double
  simTheta(RandomGenerator&                randomGenerator,
           const Acts::MaterialProperties& mat,
           double                          p,
           double                          pathcorrection,
           Acts::ParticleType              particle = Acts::pion) const final;

  /// Set configuration method
  /// @param[in] config The configuration object
  void
  setConfiguration(const Config& config);

  /// Get configuration method
  /// @return The configuration object
  Config
  getConfiguration() const
  {
    return m_config;
  }

protected:
  /// The configuraiton object
  Config m_config;
  /// Struct of Particle masses
  Acts::ParticleMasses m_particleMasses;
};

}  // end of namespace

/// Define the templated function
#include "detail/MultipleScatteringSamplerHighland.ipp"

#endif  // ACTS_FATRASTOOLS_MULTIPLESCATTERINGSAMPLERHIGHLAND_H
