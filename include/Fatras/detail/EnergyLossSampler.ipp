///////////////////////////////////////////////////////////////////
// EnergyLossSampler.cxx, ACTS project
///////////////////////////////////////////////////////////////////

// FATRAS includes
#include "FATRAS/Simulation/EnergyLossSampler.h"
#include "FATRAS/Simulation/EnergyLoss.h"
#include "FATRAS/Simulation/detail/FatrasDefinitions.h"
// ACTS includes
#include "ACTS/EventData/ParticleDefinitions.h"
#include "ACTS/Utilities/MsgMacros.h"

Fatras::EnergyLossSampler::EnergyLossSampler( const Fatras::EnergyLossSampler::Config& elConfig )
: Fatras::IEnergyLossSampler(),
  m_config()
{
    setConfiguration(elConfig);
}

Fatras::EnergyLossSampler::~EnergyLossSampler()
{}

void Fatras::EnergyLossSampler::setConfiguration(const Fatras::EnergyLossSampler::Config& elConfig ) 
{
    //!< @TODO update to configuration checking
   m_config = elConfig;   
}

Fatras::EnergyLoss Fatras::EnergyLossSampler::energyLoss( const Acts::MaterialProperties& materialProperties,
					                                      double momentum,
					                                      double pathCorrection,
					                                      Acts::PropDirection direction,
					                                      Acts::ParticleType particleHypothesis) const
{
  Fatras::EnergyLoss sampledEloss(0., 0.);
  
  if (particleHypothesis == Acts::undefined ) {
    MSG_WARNING( "undefined ParticleType, energy loss calculation cancelled" );
    return sampledEloss;
  }
  
  double pathLength = pathCorrection * materialProperties.thicknessInX0()*materialProperties.x0();  
  
  double energyLossSigma = 0.;
  double kazL    = 0.;
  
  // Evaluate the energy loss and its sigma
  double energyLoss = s_interactionFormulae.PDG_energyLoss_ionization(momentum, &(materialProperties.material()), particleHypothesis, energyLossSigma, kazL, pathLength);
  double simulatedDeltaE = fabs(energyLoss)+energyLossSigma*m_config.randomNumbers->draw(Fatras::Landau); 
 
  // giving the right sign to the energy loss
  // The sign depends on the given propagation direction 
  sampledEloss.update(-1.*direction*simulatedDeltaE, energyLossSigma, 0., 0., false);
  
  // protection due to straggling - maximum energy loss is E-m
  double m     = s_particleMasses.mass[particleHypothesis];
  double E     = sqrt(momentum*momentum+m*m);
  
  if (sampledEloss.deltaE()+E<m ) {    
    // particle stopping - rest energy
    float dRad_rest = m-E-sampledEloss.deltaE();
    sampledEloss.update(0.,0.,dRad_rest,0.,false);   
  }
   
  MSG_VERBOSE( "[eloss] created random deltaP as : " << sampledEloss.deltaE() );
 
  return sampledEloss;

}

// public interface method
double Fatras::EnergyLossSampler::dEdX(const Acts::MaterialProperties& mat,
                                     double p,
                                     Acts::ParticleType particle) const
{
  if (particle == Acts::undefined || particle == Acts::nonInteracting) return 0.; 
   
  // preparation of kinetic constants
  double m     = s_particleMasses.mass[particle];
  double E     = sqrt(p*p+m*m);
  double beta  = p/E;
  double gamma = E/m;

  // add ionization and radiation
  double dEdX = dEdXBetheBloch(mat, beta, gamma, particle) + dEdXBetheHeitler(mat, E, particle);

  // add e+e- pair production and photonuclear effect for muons at energies above 8 GeV
  if ((particle == Acts::muon) && (E > 8000.)) {
    if (E < 1.e6) {
      dEdX += - 0.5345/mat.x0() + 6.803e-5*E/mat.x0() + 2.278e-11*E*E/mat.x0() - 9.899e-18*E*E*E/mat.x0(); //E below 1 TeV
    } else {
      dEdX += - 2.986/mat.x0() + 9.253e-5*E/mat.x0(); //E above 1 TeV
    }
  }

  return dEdX;
}

double Fatras::EnergyLossSampler::dEdXBetheBloch(const Acts::MaterialProperties& mat,
                                                 double beta,
                                                 double gamma,
                                                 Acts::ParticleType particle) const
{

  if (particle == Acts::undefined || particle == Acts::nonInteracting ) return 0.;

  if (mat.averageZ()==0. || mat.zOverAtimesRho()==0. ) {
    MSG_ERROR("empty material properties pass to the EnergyLossSampler:Z,zOAtr:"<<mat.averageZ()<<","<<mat.zOverAtimesRho());
    return 0.;
  }

  // 16 eV * Z**0.9 - bring to MeV
  double iPot = 16.e-6 * std::pow(mat.averageZ(),0.9);
  // and the electron mass in MeV
  double me    = s_particleMasses.mass[Acts::electron];
  double m     = s_particleMasses.mass[particle];
  double eta2  = beta*gamma;
 
  // K/A*Z = 0.5 * 30.7075MeV/(g/mm2) * Z/A * rho[g/mm3]  / scale to mm by this
  double kaz = 0.5*s_ka_BetheBloch*mat.zOverAtimesRho();
     
  if (particle != Acts::electron){

    // density effect, only valid for high energies (gamma > 10 -> p > 1GeV for muons)
    double delta = 0.;

    /* ST replace with STEP-like coding  
    @TODO check with Sharka 
    // high energy density effect --- will be ramped up linearly
    double eplasma = 28.816e-6 * sqrt(1000.*0.5);
    delta = 2.*log(eplasma/iPot) + log(eta2) - 0.5;
    if (eta2 < 100.){
    delta *= (eta2-3.)/97.; 
    }
    */
      
    eta2 *= eta2; 

    if (gamma>10.) {
      double eplasma = 28.816e-6 * sqrt(1000.*mat.zOverAtimesRho());
      delta = 2.*log(eplasma/iPot) + log(eta2) - 1.;
    }
    
    // mass fraction
    double mfrac = me/m;
    // tmax - cut off energy
    double tMax = 2.*eta2*me/(1.+2.*gamma*mfrac+mfrac*mfrac);
    // divide by beta^2 for non-electrons
    kaz /= beta*beta;
    // return
    return  kaz*(log(2.*me*eta2*tMax/(iPot*iPot)) - 2.*(beta*beta) - delta);
  
  }
  // for electrons use slightly different BetheBloch adaption
  // see Stampfer, et al, "Track Fitting With Energy Loss", Comp. Pyhs. Comm. 79 (1994), 157-164
  return kaz*(2.*log(2.*me/iPot)+3.*log(gamma) - 1.95);
  
}

double Fatras::EnergyLossSampler::dEdXBetheHeitler(const Acts::MaterialProperties& mat,
                                                   double initialE, 
                                                   Acts::ParticleType particle) const
{

  if (particle == Acts::undefined || particle == Acts::nonInteracting ) return 0.;

  double mfrac = (s_particleMasses.mass[Acts::electron]/s_particleMasses.mass[particle]); mfrac *= mfrac;
  
  return initialE/mat.x0()*mfrac;
}