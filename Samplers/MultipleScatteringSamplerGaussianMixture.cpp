///////////////////////////////////////////////////////////////////
// MultipleScatteringSamplerGaussianMixture.cxx, ACTS project
///////////////////////////////////////////////////////////////////

// class header include
#include "Fatras/MultipleScatteringSamplerGaussianMixture.hpp"
#include "Acts/EventData/ParticleDefinitions.hpp"
#include "Fatras/MaterialInteractionEngine.hpp"
#include "Fatras/detail/FatrasDefinitions.hpp"

// ============================= Gaussian mixture model =============

// constructor
template <class RandomNumbers>
Fatras::MultipleScatteringSamplerGaussianMixture<RandomNumbers>::
    MultipleScatteringSamplerGaussianMixture(
        const MultipleScatteringSamplerGaussianMixture::Config& msConfig)
  : m_config()
{
  setConfiguration(msConfig);
}

// destructor
template <class RandomNumbers>
Fatras::MultipleScatteringSamplerGaussianMixture<RandomNumbers>::
    ~MultipleScatteringSamplerGaussianMixture()
{
}

template <class RandomNumbers>
void
Fatras::MultipleScatteringSamplerGaussianMixture<RandomNumbers>::
    setConfiguration(
        const Fatras::MultipleScatteringSamplerGaussianMixture::Config&
            msConfig)
{
  //!< @TODO update to configuration checking
  m_config = msConfig;
}

template <class RandomNumbers>
double
Fatras::MultipleScatteringSamplerGaussianMixture<RandomNumbers>::simTheta(
    const Acts::MaterialProperties& mat,
    double                          p,
    double                          pathcorrection,
    Acts::ParticleType              particle) const
{
  if (mat.thicknessInX0() <= 0. || particle == Acts::geantino) return 0.;

  // make sure the path correction is positive to avoid a floating point
  // exception
  pathcorrection *= pathcorrection < 0. ? (-1.) : (1);

  // scale the path length to the radiation length
  double t = pathcorrection * mat.thicknessInX0();

  // kinematics (relativistic)
  double m    = s_particleMasses.mass[particle];
  double E    = sqrt(p * p + m * m);
  double beta = p / E;

  double sigma2(0.);

  if (particle != Acts::electron) {
    // the highland formula
    sigma2 = s_main_RutherfordScott / (beta * p);

    if (m_config.log_include) sigma2 *= (1. + s_log_RutherfordScott * log(t));

    sigma2 *= (sigma2 * t);
  }

  else {
    // Electron multiple scattering effects - see Highland NIM 129 (1975)
    // p497-499
    // (Highland extension to the Rossi-Greisen formulation)
    // NOTE: The formula can be extended by replacing the momentum^2 term with
    // pi * pf
    sigma2 = s_main_RossiGreisen / (beta * p);
    sigma2 *= (sigma2 * t);

    if (m_config.log_include) {
      double factor = 1. + s_log_RossiGreisen * log10(10. * t);
      factor *= factor;
      sigma2 *= factor;
    }
  }

  // d_0'
  double dprime     = t / (beta * beta);
  double log_dprime = log(dprime);
  // d_0''
  double log_dprimeprime = log(std::pow(mat.averageZ(), 2.0 / 3.0) * dprime);
  // get epsilon
  double epsilon = log_dprimeprime < 0.5
      ? m_config.gausMixEpsilon_a0
          + m_config.gausMixEpsilon_a1 * log_dprimeprime
          + m_config.gausMixEpsilon_a2 * log_dprimeprime * log_dprimeprime
      : m_config.gausMixEpsilon_b0
          + m_config.gausMixEpsilon_b1 * log_dprimeprime
          + m_config.gausMixEpsilon_b2 * log_dprimeprime * log_dprimeprime;
  // the standard sigma
  double sigma1square = m_config.gausMixSigma1_a0
      + m_config.gausMixSigma1_a1 * log_dprime
      + m_config.gausMixSigma1_a2 * log_dprime * log_dprime;
  // G4 optimised / native double Gaussian model
  if (!m_config.optGaussianMixtureG4) sigma2 = 225. * dprime / (p * p);
  // throw the random number core/tail
  if (m_config.randomNumbers->drawUniform() < epsilon) {
    sigma2 *= (1. - (1. - epsilon) * sigma1square) / epsilon;
  }
  // @todo replace with GaussZiggurat later
  return s_sqrtTwo * sqrt(sigma2) * m_config.randomNumbers->drawGauss();
}
