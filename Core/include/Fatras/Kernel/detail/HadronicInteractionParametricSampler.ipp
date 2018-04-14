///////////////////////////////////////////////////////////////////
// HadronicInteractionParametricSampler.cxx, ACTS project
///////////////////////////////////////////////////////////////////

#include "Acts/EventData/ParticleDefinitions.hpp"
#include "Fatras/HadronicInteractionParametricSampler.hpp"
#include "Fatras/RandomNumberDistributions.hpp"
#include "Fatras/detail/FatrasDefinitions.hpp"

// constructor
template <class RandomGenerator>
Fatras::HadronicInteractionParametricSampler<RandomGenerator>::
    HadronicInteractionParametricSampler(
        const HadronicInteractionParametricSampler<RandomGenerator>::Config&
                                            hiConfig,
        std::unique_ptr<const Acts::Logger> logger)
  : IHadronicInteractionSampler<RandomGenerator>()
  , m_config()
  , m_logger(std::move(logger))
{
  setConfiguration(hiConfig);
}

// destructor
template <class RandomGenerator>
Fatras::HadronicInteractionParametricSampler<RandomGenerator>::
    ~HadronicInteractionParametricSampler()
{
}

template <class RandomGenerator>
void
Fatras::HadronicInteractionParametricSampler<RandomGenerator>::setConfiguration(
    const HadronicInteractionParametricSampler<RandomGenerator>::Config&
        hiConfig)
{
  //!< @TODO update to configuration checking
  m_config = hiConfig;
}

template <class RandomGenerator>
void
Fatras::HadronicInteractionParametricSampler<RandomGenerator>::setLogger(
    std::unique_ptr<const Acts::Logger> newLogger)
{
  m_logger = std::move(newLogger);
}

template <class RandomGenerator>
std::vector<Acts::ProcessVertex>
Fatras::HadronicInteractionParametricSampler<RandomGenerator>::
    doHadronicInteraction(RandomGenerator&      randomGenerator,
                          double                time,
                          const Acts::Vector3D& position,
                          const Acts::Vector3D& momentum,
                          Acts::ParticleType    particle) const
{
  return getHadState(randomGenerator,
                     time,
                     momentum.mag(),
                     position,
                     momentum.unit(),
                     particle);
}

template <class RandomGenerator>
std::vector<Acts::ProcessVertex>
Fatras::HadronicInteractionParametricSampler<RandomGenerator>::getHadState(
    RandomGenerator&      randomGenerator,
    double                time,
    double                p,
    const Acts::Vector3D& vertex,
    const Acts::Vector3D& particleDir,
    Acts::ParticleType    particle) const
{

  // the flat distribution
  Fatras::UniformDist flatDist(0., 1.);

  std::vector<Acts::ProcessVertex> children;

  // sampling of hadronic interaction
  double m = m_particleMasses.mass[particle];
  double E = sqrt(p * p + m * m);

  // get the maximum multiplicity
  double multiplicity_max = 0.25 * E / Acts::units::_GeV + 18.;

  // multiplicity distribution
  double randx, randy, arg = 0.;

  double p1 = 0.;
  double p2 = 0.;

  if (E > 15. * Acts::units::_GeV) {
    p1 = 8.69;
    p2 = 2.34;
  } else {
    p1 = 6.77;
    p2 = 2.02;
  }

  for (;;) {
    randx = 30. * flatDist(randomGenerator);
    randy = 1. * flatDist(randomGenerator);
    arg   = exp(-0.5 * ((randx - p1) / p2 + exp(-(randx - p1) / p2)));
    if (randy < arg && randx > 3 && randx < multiplicity_max) break;
  }

  randx *= (1.2 - 0.4 * exp(-p));  // trying to adjust

  int Npart = (int)randx;

  // protection against Npart < 3
  if (Npart < 3) {
    ACTS_VERBOSE(
        "[ had ] Number of particles smaller than 3, parameterisation not "
        "valid."
        << " Doing Nothing");
    return children;
  }

  ACTS_VERBOSE("[ had ] interaction of " << particle << " with " << Npart
                                         << " outgoing particles ");

  // record the interaction

  // ------ now create the new hadrons ------
  ACTS_DEBUG("[ had ] create hadronic shower for particle " << particle);

  // create the genParticles

  ACTS_VERBOSE(
      "[ had ] incoming particle energy | mass | momentum " << E << " | " << m
                                                            << " | "
                                                            << p
                                                            << " | ");

  /*!>@TODO: If isSecondary&&m_config.cutChain do not save the children
   * In ATLAS this is done checking the parent barcode
   * if (m_config.cutChain && ( parent->barcode()>100000 || parent->barcode()==0
   * ) ) */
  bool isSecondary = false;
  if (m_config.cutChain && isSecondary) {
    ACTS_VERBOSE(
        "[ had ] interaction initiated by a secondary particle, no children "
        "saved ");
    return children;
  }

  int gen_part = 0;

  // new sampling: sample particle type and energy in the CMS frame of outgoing
  // particles
  // creation of shower particles
  double                          chargedist = 0.;
  std::vector<double>             charge(Npart);
  std::vector<Acts::ParticleType> childType(Npart);
  std::vector<double>             newm(Npart);
  std::vector<int>                pdgid(Npart);

  // children type sampling  : simplified
  // double pif = 0.19;
  // double nef = 0.20;
  // double prf = 0.20;

  // sample heavy particles (alpha) but don't save
  double pif = 0.10;
  double nef = 0.30;
  double prf = 0.30;

  if (particle == Acts::pion || particle == Acts::kaon || particle == Acts::pi0
      || particle == Acts::k0) {
    pif = 0.15;
    nef = 0.25;
    prf = 0.25;
  }
  if (particle == Acts::proton) {
    pif = 0.06;
    nef = 0.25;
    prf = 0.35;
  }
  if (particle == Acts::neutron) {
    pif = 0.03;
    nef = 0.35;
    prf = 0.17;
  }

  for (int i = 0; i < Npart; i++) {
    chargedist = flatDist(randomGenerator);
    if (chargedist < pif) {
      charge[i]    = 0.;
      childType[i] = Acts::pi0;
      newm[i]      = m_particleMasses.mass[Acts::pi0];  // MeV
      pdgid[i]     = 111;
      continue;
    }
    if (chargedist < 2 * pif) {
      charge[i]    = 1.;
      childType[i] = Acts::pion;
      newm[i]      = m_particleMasses.mass[Acts::pion];  // MeV
      pdgid[i]     = 211;
      continue;
    }
    if (chargedist < 3 * pif) {
      charge[i]    = -1.;
      childType[i] = Acts::pion;
      newm[i]      = m_particleMasses.mass[Acts::pion];  // MeV
      pdgid[i]     = -211;
      continue;
    }
    if (chargedist < 3 * pif + nef) {
      charge[i]    = 0.;
      childType[i] = Acts::neutron;
      newm[i]      = 939.565;  // MeV
      pdgid[i]     = 2112;     // neutron
      continue;
    }
    if (chargedist < 3 * pif + nef + prf) {
      charge[i]    = 1.;
      childType[i] = Acts::proton;
      newm[i]      = m_particleMasses.mass[Acts::proton];  // MeV
      pdgid[i]     = 2212;
      continue;
    }
    charge[i]    = 2.;
    childType[i] = Acts::proton;
    newm[i]      = 4000.;
    pdgid[i]     = 20000;
  }

  // move the incoming particle type forward
  //!>@TODO Do we really need this?
  // If so we need to get the parent charge
  // Asking the eCell??
  //
  //   if ( childType[0] != particle ) {
  //     for (int i=1; i<Npart; i++) {
  //       if (childType[i]==particle) {
  //         childType[i]=childType[0];
  //         childType[0]=particle;
  //         double cho = charge[i];
  //         charge[i]=charge[0];
  //         charge[0]=parent ? parent->charge() : cho;
  // 	newm[i]=m_particleMasses.mass[childType[i]]; // MeV
  // 	newm[0]=m_particleMasses.mass[childType[0]]; // MeV
  //         break;
  //       }
  //     }
  //   }

  std::vector<double> mom(Npart);
  std::vector<double> th(Npart);
  std::vector<double> ph(Npart);

  // sample first particle energy fraction and random momentum direction
  double eps = 2. / Npart;
  double rnd = flatDist(randomGenerator);
  mom[0]     = 0.5 * pow(eps, rnd);
  th[0]      = acos(2 * flatDist(randomGenerator) - 1.);
  ph[0]      = 2 * M_PI * flatDist(randomGenerator);

  // toss particles around in a way which preserves the total momentum
  // (0.,0.,0.) at this point
  //!>@TODO shoot first particle along the impact direction preferentially

  Acts::Vector3D ptemp(mom[0] * sin(th[0]) * cos(ph[0]),
                       mom[0] * sin(th[0]) * sin(ph[0]),
                       mom[0] * cos(th[0]));
  double ptot = mom[0];

  double theta = 0.;
  double phi   = 0.;
  for (int i = 1; i < Npart - 2; i++) {
    eps    = 1. / (Npart - i);
    mom[i] = (eps + flatDist(randomGenerator) * (1 - eps)) * (1 - ptot);
    if (ptemp.mag() < 1 - ptot) {
      while (fabs(ptemp.mag() - mom[i]) > 1 - ptot - mom[i]) {
        mom[i] = (eps + flatDist(randomGenerator) * (1 - eps)) * (1 - ptot);
      }
    }
    // max p remaining
    double p_rem  = 1 - ptot - mom[i];
    double cthmax = fmin(
        1.,
        (-ptemp.mag() * ptemp.mag() - mom[i] * mom[i] + p_rem * p_rem) / 2
            / ptemp.mag()
            / mom[i]);
    double rnd = flatDist(randomGenerator);
    theta      = acos((cthmax + 1.) * rnd - 1.);
    phi        = 2 * M_PI * flatDist(randomGenerator);
    Acts::Vector3D test(
        sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
    // create the rotation
    //!>@TODO Check if this is right
    Acts::RotationMatrix3D rotation;
    rotation = Acts::AngleAxis3D(ptemp.phi(), Acts::Vector3D::UnitZ())
        * Acts::AngleAxis3D(ptemp.theta(), Acts::Vector3D::UnitY());
    // create the transform
    Acts::Transform3D transform(rotation, Acts::Vector3D(0., 0., 0.));
    Acts::Vector3D    dnew = transform * test;
    th[i]                  = dnew.theta();
    ph[i]                  = dnew.phi();
    ptemp += mom[i] * dnew;
    ptot += mom[i];
  }

  eps            = 0.5;
  mom[Npart - 2] = pow(eps, flatDist(randomGenerator)) * (1 - ptot);
  mom[Npart - 1] = 1 - ptot - mom[Npart - 2];

  if (ptemp.mag() < 1 - ptot) {
    while (mom[Npart - 1] + mom[Npart - 2] < ptemp.mag()) {
      mom[Npart - 2] = pow(eps, flatDist(randomGenerator)) * (1 - ptot);
      mom[Npart - 1] = 1 - ptot - mom[Npart - 2];
    }
  }
  if (ptemp.mag() < fabs(mom[Npart - 1] - mom[Npart - 2])) {
    double diff    = ptemp.mag() * flatDist(randomGenerator);
    double sum     = mom[Npart - 1] - mom[Npart - 2];
    mom[Npart - 2] = 0.5 * (sum + diff);
    mom[Npart - 1] = 0.5 * (sum - diff);
  }
  double cth = (-ptemp.mag() * ptemp.mag() - mom[Npart - 2] * mom[Npart - 2]
                + mom[Npart - 1] * mom[Npart - 1])
      / 2 / ptemp.mag() / mom[Npart - 2];
  if (fabs(cth) > 1.) cth = (cth > 0.) ? 1. : -1.;

  theta = acos(cth);
  phi   = 2 * M_PI * flatDist(randomGenerator);
  Acts::Vector3D test(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
  // create the rotation
  //!>@TODO Check if this is right
  Acts::RotationMatrix3D rotation;
  rotation = Acts::AngleAxis3D(ptemp.phi(), Acts::Vector3D::UnitZ())
      * Acts::AngleAxis3D(ptemp.theta(), Acts::Vector3D::UnitY());
  // create the transform
  Acts::Transform3D transform(rotation, Acts::Vector3D(0., 0., 0.));
  Acts::Vector3D    dnew = transform * test;
  th[Npart - 2]          = dnew.theta();
  ph[Npart - 2]          = dnew.phi();
  ptemp += mom[Npart - 2] * dnew;
  Acts::Vector3D dlast = -ptemp;
  th[Npart - 1]        = dlast.theta();
  ph[Npart - 1]        = dlast.phi();

  // particle sampled, rotate, boost and save final state
  double etot = 0.;
  for (int i = 0; i < Npart; i++)
    etot += sqrt(mom[i] * mom[i] + newm[i] * newm[i]);
  double summ = 0.;
  for (int i = 0; i < Npart; i++) summ += newm[i];

  // std::cout <<"hadronic interaction: current energy, expected :"<< etot
  // <<","<< sqrt(summ*summ+2*summ*p+m*m)<< std::endl;
  // rescale (roughly) to the expected energy
  float scale = sqrt(summ * summ + 2 * summ * p + m * m) / etot;
  etot        = 0.;
  for (int i = 0; i < Npart; i++) {
    mom[i] *= scale;
    etot += sqrt(mom[i] * mom[i] + newm[i] * newm[i]);
  }

  /*!>@TODO Implement this with Eigen vector
   * the code is commented for the moment
   * And using std::vector < Acts::Vector3D > */

  typedef Acts::ActsVector<double, 4> ActsLorentzVector;
  std::vector<ActsLorentzVector> childBoost(Npart);

  // the boost vector
  ActsLorentzVector bv(p * particleDir.unit().x(),
                       p * particleDir.unit().y(),
                       p * particleDir.unit().z(),
                       sqrt(etot * etot + p * p));

  Acts::Vector3D in(0., 0., 0.);
  Acts::Vector3D fin(0., 0., 0.);
  for (int i = 0; i < Npart; i++) {
    // direction in the cneter of mass frame
    Acts::Vector3D dirCms(
        sin(th[i]) * cos(ph[i]), sin(th[i]) * sin(ph[i]), cos(th[i]));
    // child momentum in the center of mass frame
    Acts::Vector3D childP = mom[i] * dirCms;
    in += childP;
    // child lorentz vector in cms frame
    ActsLorentzVector childLorenz(childP.x(),
                                  childP.y(),
                                  childP.z(),
                                  sqrt(mom[i] * mom[i] + newm[i] * newm[i]));
    double bx = bv.x() / bv[3];
    double by = bv.y() / bv[3];
    double bz = bv.z() / bv[3];

    // Boost this Lorentz vector
    double b2     = bx * bx + by * by + bz * bz;
    double lgamma = 1.0 / sqrt(1.0 - b2);
    double bp
        = bx * childLorenz.x() + by * childLorenz.y() + bz * childLorenz.z();
    double gamma2 = b2 > 0 ? (lgamma - 1.0) / b2 : 0.0;

    double bChildPx
        = (childLorenz[0] + gamma2 * bp * bx + lgamma * bx * childLorenz[3]);
    double bChildPy
        = (childLorenz[1] + gamma2 * bp * by + lgamma * by * childLorenz[3]);
    double bChildPz
        = (childLorenz[2] + gamma2 * bp * bz + lgamma * bz * childLorenz[3]);
    double bChildT = (lgamma * (childLorenz[3] + bp));

    childBoost[i] = ActsLorentzVector(bChildPx, bChildPy, bChildPz, bChildT);
    fin += Acts::Vector3D(bChildPx, bChildPy, bChildPz);
  }

  // Add children to the vector of children
  std::vector<Acts::ParticleProperties> pIngoing = {};
  std::vector<Acts::ParticleProperties> pOutgoing;
  unsigned short                        numChildren = 0;

  for (int i = 0; i < Npart; i++) {
    if (pdgid[i] < 10000) {
      /*!>@TODO Getting the momentum from the boost */
      Acts::Vector3D childP = Acts::Vector3D(
          childBoost[i].x(), childBoost[i].y(), childBoost[i].z());
      if (childP.mag() > m_config.minimumHadOutEnergy) {
        // create the particle and increase the number of children
        // auto part = convert(childType[i], charge[i], false);
        pOutgoing.push_back(
            Acts::ParticleProperties(childP, m, charge[i], pdgid[i], 100000));
        numChildren++;
      }
      // increase the number of generated particles
      gen_part++;
    }
  }  // particle loop

  children.push_back(Acts::ProcessVertex(
      vertex, time, m_config.processCode, pIngoing, pOutgoing));
  ACTS_VERBOSE("[ had ] it was kinematically possible to create "
               << gen_part
               << " shower particles and "
               << numChildren
               << " particles have been collected ");

  return children;
}
