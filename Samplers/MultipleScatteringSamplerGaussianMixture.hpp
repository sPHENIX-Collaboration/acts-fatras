///////////////////////////////////////////////////////////////////
// MultipleScatteringSamplerGaussianMixture.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_FATRASTOOLS_MULTIPLESCATTERINGSAMPLERGAUSSIANMIXTURE_H
#define ACTS_FATRASTOOLS_MULTIPLESCATTERINGSAMPLERGAUSSIANMIXTURE_H 1

#include <memory>
#include "ACTS/EventData/ParticleDefinitions.hpp"
#include "ACTS/Material/MaterialProperties.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "Fatras/IMultipleScatteringSampler.hpp"

namespace Fatras {

/** @class MultipleScatteringSamplerGaussianMixture
 *
 * ========= Gaussian mixture model Fruehwirth, Regler Nucl. Inst. Methods A
 * 456(2001) =========
 * @author Andreas Salzburger <Andreas.Salzburger@cern.ch>
 * @author Noemi Calace       <Noemi.Calace@cern.ch>
 * @author Artem Basalaev     <Artem.Baralaev@cern.ch>
 *
 */

template <class RandomNumbers>
class MultipleScatteringSamplerGaussianMixture
    : virtual public IMultipleScatteringSampler
{
public:
  /** Config
  Configuration object for this MultipleScatteringSampler*/
  struct Config
  {
    std::shared_ptr<RandomNumbers>
         randomNumbers;        /** Random Generator service  */
    bool log_include;          /** boolean switch to include log term  */
    bool optGaussianMixtureG4; /** modifies the Fruehwirth/Regler model to fit
                                  with G4 */
    double
        gausMixSigma1_a0; /** gaussian mixture paramters for sigma and epsion */
    double
        gausMixSigma1_a1; /** gaussian mixture paramters for sigma and epsion */
    double
           gausMixSigma1_a2; /** gaussian mixture paramters for sigma and epsion */
    double gausMixEpsilon_a0; /** gaussian mixture paramters for sigma and
                                 epsion */
    double gausMixEpsilon_a1; /** gaussian mixture paramters for sigma and
                                 epsion */
    double gausMixEpsilon_a2; /** gaussian mixture paramters for sigma and
                                 epsion */
    double gausMixEpsilon_b0; /** gaussian mixture paramters for sigma and
                                 epsion */
    double gausMixEpsilon_b1; /** gaussian mixture paramters for sigma and
                                 epsion */
    double gausMixEpsilon_b2; /** gaussian mixture paramters for sigma and
                                 epsion */

    Config()
      : randomNumbers(nullptr)
      , log_include(true)
      , optGaussianMixtureG4(true)
      , gausMixSigma1_a0(8.471e-1)
      , gausMixSigma1_a1(3.347e-2)
      , gausMixSigma1_a2(-1.843e-3)
      , gausMixEpsilon_a0(4.841e-2)
      , gausMixEpsilon_a1(6.348e-3)
      , gausMixEpsilon_a2(6.096e-4)
      , gausMixEpsilon_b0(-1.908e-2)
      , gausMixEpsilon_b1(1.106e-1)
      , gausMixEpsilon_b2(-5.729e-3)
    {
    }
  };

  /** AlgTool like constructor */
  MultipleScatteringSamplerGaussianMixture(const Config& msConfig);

  /**Virtual destructor*/
  virtual ~MultipleScatteringSamplerGaussianMixture();

  /** Calculate the theta introduced by multiple scattering,
   *          according to the RutherFord-Scott Formula
   */
  double
  simTheta(const Acts::MaterialProperties& mat,
           double                          p,
           double                          pathcorrection,
           Acts::ParticleType              particle = Acts::pion) const final;

  /** Set configuration method */
  void
  setConfiguration(const Config& eeConfig);

  /** Get configuration method */
  Config
  getConfiguration() const;

protected:
  Config m_config;  //!< configuration object
};

/** Return the configuration object */
inline MultipleScatteringSamplerGaussianMixture::Config
MultipleScatteringSamplerGaussianMixture::getConfiguration() const
{
  return m_config;
}

}  // end of namespace

//!< define the templated function
#include "detail/MultipleScatteringSamplerGaussianMixture.ipp"

#endif  // ACTS_FATRASTOOLS_MULTIPLESCATTERINGSAMPLERGAUSSIANMIXTURE_H
