///////////////////////////////////////////////////////////////////
// ElectronEnergyLossSampler.cxx, ACTS project
///////////////////////////////////////////////////////////////////

// class header include
#include "FATRAS/ElectronEnergyLossSampler.h"
#include "FATRAS/EnergyLoss.h"
// EventData module
#include "ACTS/EventData/ParticleProperties.h"

// static partilce masses
Acts::ParticleMasses Fatras::ElectronEnergyLossSampler::s_particleMasses;

Fatras::ElectronEnergyLossSampler::ElectronEnergyLossSampler( const Fatras::ElectronEnergyLossSampler::Config& elConfig )
  : Fatras::IEnergyLossSampler(),
    m_config()
{
    setConfiguration(elConfig);
}

Fatras::ElectronEnergyLossSampler::~ElectronEnergyLossSampler()
{}


void Fatras::ElectronEnergyLossSampler::setConfiguration(const Fatras::ElectronEnergyLossSampler::Config& elConfig ) 
{
    //!< @TODO update to configuration checking
   m_config = elConfig;   
}

Fatras:EnergyLoss Fatras::ElectronEnergyLossSampler::energyLoss( const Acts::MaterialProperties& materialProperties,
							                                     double pInitial,
							                                     double pathCorrection,
							                                     Acts::PropDirection direction,
							                                     Acts::ParticleHypothesis) const
{  
    // start with a default energy loss
  Acts::EnergyLoss sampledEloss;
  
  double pathLength = pathCorrection*materialProperties.thicknessInX0();
  
  if (pathLength==0.) return sampledEloss;
  
  double p    = pInitial;
  double me   = s_particleMasses.mass[Acts::electron];
  double E    = sqrt(p*p+me*me);
  
  // the following formulas are imported from STEP
  // preparation of kinetic constants
 
  double beta  = p/E;
  double gamma = E/me;
  double eta2 = beta*gamma; eta2 *= eta2;  
  
  //Ionization - Bethe-Bloch
  double I = 16.e-6 * std::pow(materialProperties.averageZ(),0.9); //16 eV * Z**0.9 - bring to MeV
  
  double delta = 0.;
  if (gamma > 10.) {
    double eplasma = 28.816e-6 * sqrt(1000.*materialProperties.zOverAtimesRho());
    delta = 2.*log(eplasma/I) + log(eta2) - 1.;
  }
  
  //K/A*Z = 0.5 * 30.7075MeV/(g/mm2) * Z/A * rho[g/mm3]  / scale to mm by this
  double kaz = 0.5*30.7075*materialProperties.zOverAtimesRho();
  double kazL = kaz*pathCorrection;
  
  // for electrons use slightly different BetheBloch adaption
  // see Stampfer, et al, "Track Fitting With Energy Loss", Comp. Pyhs. Comm. 79 (1994), 157-164
  //
  // the landau sigmaL is path length dependent
  //    PDG formula 32.11 for MOP value from http://http://pdg.lbl.gov/2014/reviews/rpp2014-rev-passage-particles-matter.pdf
  //
  
  double MOP =  -kazL*kaz*(log(2.*me*eta2/I) + log(kazL/I) + 0.2 - (beta*beta) - delta);

  double energyLossSigma = 0.424*4.*kazL; //0.424: scale factor from sigma to FWHM
  
  double simulatedDeltaE = fabs(MOP)+energyLossSigma*m_config.randomNumbers->draw(Fatras::Landau); 

  //Bethe-Heitler for electron brem description as described here:
  // "A Gaussian-mixture approximation of the Bethe–Heitler model of electron energy loss by bremsstrahlung"
  // R. Frühwirth
  
  double u = m_config.randomNumbers->draw(Acts::Gamma, pathLength/log(2.), 1.);
  double z = exp( -1. * u );
  double deltaE(0.);
  if ( direction == Acts::alongMomentum )
    deltaE = m_config.scaleFactor*E * ( z - 1. );
  else
    deltaE = m_config.scaleFactor*E * ( 1. / z - 1. );
  
  simulatedDeltaE+=fabs(deltaE);
  
  // protection due to straggling - maximum energy loss is E
  simulatedDeltaE = simulatedDeltaE*simulatedDeltaE > E*E ? E : simulatedDeltaE;
  
  // giving the right sign to the energy loss
  // The sign depends on the given propagation direction 
  sampledEloss.update(-1.*direction*simulatedDeltaE, energyLossSigma, 0., 0., false);
  
  MSG_VERBOSE( "[electron Eloss] created random deltaP as : " << sampledEloss->deltaE() );
  
  return sampledEloss;
     
}