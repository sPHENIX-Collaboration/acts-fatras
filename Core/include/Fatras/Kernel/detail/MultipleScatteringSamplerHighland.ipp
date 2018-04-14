///////////////////////////////////////////////////////////////////
// MultipleScatteringSamplerHighland.cxx, ACTS project
///////////////////////////////////////////////////////////////////

#include <Acts/EventData/ParticleDefinitions.hpp>
#include <Acts/Utilities/MaterialInteraction.hpp>
#include "Fatras/RandomNumberDistributions.hpp"
#include "Fatras/detail/FatrasDefinitions.hpp"

// constructor
template <class RandomGenerator>
Fatras::MultipleScatteringSamplerHighland<RandomGenerator>::
    MultipleScatteringSamplerHighland(
        const MultipleScatteringSamplerHighland::Config& config)
  : Fatras::IMultipleScatteringSampler<RandomGenerator>(), m_config()
{
  setConfiguration(config);
}

template <class RandomGenerator>
void
Fatras::MultipleScatteringSamplerHighland<RandomGenerator>::setConfiguration(
    const Fatras::MultipleScatteringSamplerHighland<RandomGenerator>::Config&
        config)
{
  /// @todo update to configuration checking
  m_config = config;
}

template <class RandomGenerator>
double
Fatras::MultipleScatteringSamplerHighland<RandomGenerator>::simTheta(
    RandomGenerator&                randomGenerator,
    const Acts::MaterialProperties& mat,
    double                          p,
    double                          pathcorrection,
    Acts::ParticleType              particle) const
{
  // Create a random gauss distribution between in the intervall [0,1]
  Fatras::GaussDist gaussDist(0., 1.);

  if (mat.thicknessInX0() <= 0. || particle == Acts::undefined) return 0.;
  // make sure the path correction is positive to avoid a floating point
  // exception
  pathcorrection *= pathcorrection < 0. ? (-1.) : (1);

  // scale the path length to the radiation length
  double t = pathcorrection * mat.thicknessInX0();

  // kinematics (relativistic)
  double m    = m_particleMasses.mass[particle];
  double E    = sqrt(p * p + m * m);
  double beta = p / E;

  double sigma2(0.);

  double sigma = Acts::sigmaMS(t, p, beta);
  sigma2       = sigma * sigma;

  if (particle == Acts::electron) {
    // Electron multiple scattering effects - see Highland NIM 129 (1975)
    // p497-499
    // (Highland extension to the Rossi-Greisen formulation)
    // NOTE: The formula can be extended by replacing the momentum^2 term with
    // pi * pf
    sigma2 = constants::main_RossiGreisen / (beta * p);
    sigma2 *= (sigma2 * t);

    if (m_config.log_include) {
      double factor = 1. + constants::log_RossiGreisen * log10(10. * t);
      factor *= factor;
      sigma2 *= factor;
    }
  }
  // returned scaled by the projection factor
  // @todo before gaussZiggurat was used, currently use only gauss until it is
  // implemented
  double returnV = M_SQRT2 * std::sqrt(sigma2) * gaussDist(randomGenerator);
  return returnV;
}
