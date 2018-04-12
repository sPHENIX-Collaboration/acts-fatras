///////////////////////////////////////////////////////////////////
// EnergyLossSampler.cxx, ACTS project
///////////////////////////////////////////////////////////////////

#include <ACTS/EventData/ParticleDefinitions.hpp>
#include <ACTS/Utilities/MaterialInteraction.hpp>

#include "Fatras/EnergyLoss.hpp"
#include "Fatras/EnergyLossSampler.hpp"
#include "Fatras/RandomNumberDistributions.hpp"
#include "Fatras/detail/FatrasDefinitions.hpp"

template <class RandomGenerator>
Fatras::EnergyLossSampler<RandomGenerator>::EnergyLossSampler(
    const Fatras::EnergyLossSampler<RandomGenerator>::Config& elConfig,
    std::unique_ptr<const Acts::Logger>                       logger)
  : Fatras::IEnergyLossSampler<RandomGenerator>()
  , m_config()
  , m_logger(std::move(logger))
{
  setConfiguration(elConfig);
}

template <class RandomGenerator>
void
Fatras::EnergyLossSampler<RandomGenerator>::setConfiguration(
    const Fatras::EnergyLossSampler<RandomGenerator>::Config& elConfig)
{
  // @todo update to configuration checking
  m_config = elConfig;
}

template <class RandomGenerator>
void
Fatras::EnergyLossSampler<RandomGenerator>::setLogger(
    std::unique_ptr<const Acts::Logger> newLogger)
{
  m_logger = std::move(newLogger);
}

template <class RandomGenerator>
Fatras::EnergyLoss
Fatras::EnergyLossSampler<RandomGenerator>::energyLoss(
    RandomGenerator&                randomGenerator,
    const Acts::MaterialProperties& materialProperties,
    double                          momentum,
    double                          pathCorrection,
    Acts::NavigationDirection             direction,
    Acts::ParticleType              particleHypothesis) const
{
  Fatras::EnergyLoss sampledEloss(0., 0.);
  if (particleHypothesis == Acts::undefined) {
    ACTS_WARNING("undefined ParticleType, energy loss calculation cancelled");
    return sampledEloss;
  }

  // calculate the path length
  double pathLength = pathCorrection * materialProperties.thickness();
  // Evaluate the energy loss and its sigma
  auto eLoss = Acts::ionizationEnergyLoss_mop(momentum,
                                              materialProperties.material(),
                                              particleHypothesis,
                                              m_particleMasses,
                                              pathLength);
  double energyLoss = eLoss.first;
  // the uncertainty of the mean energy loss
  double energyLossSigma = eLoss.second;
  // Create a random landau distribution between in the intervall [0,1]
  Fatras::LandauDist landauDist(0., 1.);
  double             landau = landauDist(randomGenerator);
  // Simulate the energy loss
  double simulatedDeltaE = fabs(energyLoss) * m_config.scalorMOP
      + energyLossSigma * landau * m_config.scalorSigma;

  // giving the right sign to the energy loss
  // the sign depends on the given propagation direction

  sampledEloss.update(
      -1. * direction * simulatedDeltaE, energyLossSigma, 0., 0., false);

  // protection due to straggling - maximum energy loss is E-m
  double m = m_particleMasses.mass[particleHypothesis];
  double E = sqrt(momentum * momentum + m * m);

  if ((sampledEloss.deltaE() + E) < m) {
    // particle stopping - rest energy
    float dRad_rest = m - E - sampledEloss.deltaE();
    sampledEloss.update(0., 0., dRad_rest, 0., false);
  }
  return sampledEloss;
}

// public interface method
/*template <class RandomGenerator>
double Fatras::EnergyLossSampler<RandomGenerator>::dEdX(
    const Acts::MaterialProperties& mat, double p,
    Acts::ParticleType particle) const {
  if (particle == Acts::undefined || particle == Acts::nonInteracting)
    return 0.;

  // preparation of kinetic constants
  double m = m_particleMasses.mass[particle];
  double E = sqrt(p * p + m * m);
  double beta = p / E;
  double gamma = E / m;

  // add ionization and radiation
  double dEdX = dEdXBetheBloch(mat, beta, gamma, particle) +
                dEdXBetheHeitler(mat, E, particle);

  // add e+e- pair production and photonuclear effect for muons at energies
  // above 8 GeV
  if ((particle == Acts::muon) && (E > 8000.)) {
    if (E < 1.e6) {
      dEdX += -0.5345 / mat.averageX0() + 6.803e-5 * E / mat.averageX0() +
              2.278e-11 * E * E / mat.averageX0() -
              9.899e-18 * E * E * E / mat.averageX0();  // E below 1 TeV
    } else {
      dEdX += -2.986 / mat.averageX0() +
              9.253e-5 * E / mat.averageX0();  // E above 1 TeV
    }
  }

  return dEdX;
}

template <class RandomGenerator>
double Fatras::EnergyLossSampler<RandomGenerator>::dEdXBetheBloch(
    const Acts::MaterialProperties& mat, double beta, double gamma,
    Acts::ParticleType particle) const {
  if (particle == Acts::undefined || particle == Acts::nonInteracting)
    return 0.;

  if (mat.averageZ() == 0. || mat.zOverAtimesRho() == 0.) {
    ACTS_ERROR(
        "empty material properties pass to the EnergyLossSampler:Z,zOAtr:"
        << mat.averageZ() << "," << mat.zOverAtimesRho());
    return 0.;
  }

  // 16 eV * Z**0.9 - bring to MeV
  double iPot = 16.e-6 * std::pow(mat.averageZ(), 0.9);
  // and the electron mass in MeV
  double me = m_particleMasses.mass[Acts::electron];
  double m = m_particleMasses.mass[particle];
  double eta2 = beta * gamma;

  // K/A*Z = 0.5 * 30.7075MeV/(g/mm2) * Z/A * rho[g/mm3]  / scale to mm by this
  double kaz = 0.5 * constants::ka_BetheBloch * mat.zOverAtimesRho();

  if (particle != Acts::electron) {
    // density effect, only valid for high energies (gamma > 10 -> p > 1GeV for
    // muons)
    double delta = 0.;

    // ST replace with STEP-like coding
    // @todo check with Sharka
    // // high energy density effect --- will be ramped up linearly
    // double eplasma = 28.816e-6 * sqrt(1000.*0.5);
    // delta = 2.*log(eplasma/iPot) + log(eta2) - 0.5;
    // if (eta2 < 100.){
    // delta *= (eta2-3.)/97.;
    // }

    eta2 *= eta2;

    if (gamma > 10.) {
      double eplasma = 28.816e-6 * sqrt(1000. * mat.zOverAtimesRho());
      delta = 2. * log(eplasma / iPot) + log(eta2) - 1.;
    }

    // mass fraction
    double mfrac = me / m;
    // tmax - cut off energy
    double tMax = 2. * eta2 * me / (1. + 2. * gamma * mfrac + mfrac * mfrac);
    // divide by beta^2 for non-electrons
    kaz /= beta * beta;
    // return
    return kaz * (log(2. * me * eta2 * tMax / (iPot * iPot)) -
                  2. * (beta * beta) - delta);
  }
  // for electrons use slightly different BetheBloch adaption
  // see Stampfer, et al, "Track Fitting With Energy Loss", Comp. Pyhs. Comm. 79
  // (1994), 157-164
  return kaz * (2. * log(2. * me / iPot) + 3. * log(gamma) - 1.95);
}

template <class RandomGenerator>
double Fatras::EnergyLossSampler<RandomGenerator>::dEdXBetheHeitler(
    const Acts::MaterialProperties& mat, double initialE,
    Acts::ParticleType particle) const {
  if (particle == Acts::undefined || particle == Acts::nonInteracting)
    return 0.;

  double mfrac =
      (m_particleMasses.mass[Acts::electron] / m_particleMasses.mass[particle]);
  mfrac *= mfrac;

  return initialE / mat.averageX0() * mfrac;
}*/
