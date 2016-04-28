///////////////////////////////////////////////////////////////////
// MultipleScatteringSamplerGaussianMixture.cxx, ACTS project
///////////////////////////////////////////////////////////////////

// class header include
#include "FATRAS/MultipleScatteringSamplerGaussianMixture.h"
// ACTS include
#include "ACTS/EventData/ParticleProperties.h"

// static particle masses
Acts::ParticleMasses Fatras::MultipleScatteringSamplerGaussianMixture::s_particleMasses;
// static doubles
double Fatras::MultipleScatteringSamplerGaussianMixture::s_main_RutherfordScott = 13.6;
double Fatras::MultipleScatteringSamplerGaussianMixture::s_log_RutherfordScott  =  0.038;

double Fatras::MultipleScatteringSamplerGaussianMixture::s_main_RossiGreisen = 17.5;
double Fatras::MultipleScatteringSamplerGaussianMixture::s_log_RossiGreisen  =  0.125;

// ============================= Gaussian mixture model =============
double Fatras::MultipleScatteringSamplerGaussianMixture::s_gausMixSigma1_a0  =  8.471e-1;
double Fatras::MultipleScatteringSamplerGaussianMixture::s_gausMixSigma1_a1  =  3.347e-2;
double Fatras::MultipleScatteringSamplerGaussianMixture::s_gausMixSigma1_a2  = -1.843e-3;

double Fatras::MultipleScatteringSamplerGaussianMixture::s_gausMixEpsilon_a0 =  4.841e-2;
double Fatras::MultipleScatteringSamplerGaussianMixture::s_gausMixEpsilon_a1 =  6.348e-3;
double Fatras::MultipleScatteringSamplerGaussianMixture::s_gausMixEpsilon_a2 =  6.096e-4;

double Fatras::MultipleScatteringSamplerGaussianMixture::s_gausMixEpsilon_b0 = -1.908e-2;
double Fatras::MultipleScatteringSamplerGaussianMixture::s_gausMixEpsilon_b1 =  1.106e-1;
double Fatras::MultipleScatteringSamplerGaussianMixture::s_gausMixEpsilon_b2 = -5.729e-3;

double Fatras::MultipleScatteringSamplerGaussianMixture::s_projectionFactor  =  sqrt(2.);

// constructor
Fatras::MultipleScatteringSamplerGaussianMixture::MultipleScatteringSamplerGaussianMixture(const MultipleScatteringSamplerGaussianMixture::Config& msConfig) 
: m_config()
{
    setConfiguration(msConfig);
}

// destructor
Fatras::MultipleScatteringSamplerGaussianMixture::~MultipleScatteringSamplerGaussianMixture()
{}

void Fatras::MultipleScatteringSamplerGaussianMixture::seConfiguration(const Fatras::MultipleScatteringSamplerGaussianMixture::Config& msConfig ) 
{
    //!< @TODO update to configuration checking
   m_config = msConfig;   
}

double Fatras::MultipleScatteringSamplerGaussianMixture::simTheta(const Acts::MaterialProperties& mat,
								                                  double p,
								                                  double pathcorrection,
								                                  Acts::ParticleHypothesis particle) const
{
  if (mat.thicknessInX0()<=0. || particle==Acts::geantino) return 0.;
 
  // make sure the path correction is positive to avoid a floating point exception
  pathcorrection *= pathcorrection < 0. ? (-1.) : (1) ;
 
  // scale the path length to the radiation length
  double t = pathcorrection * mat.thicknessInX0();

  // kinematics (relativistic)
  double m    = s_particleMasses.mass[particle];
  double E    = sqrt(p*p + m*m);
  double beta = p/E;
 
  double sigma2(0.);
 
  if (particle != Acts::electron) {

    // the highland formula
    sigma2 = s_main_RutherfordScott/(beta*p);

    if (m_config.log_include)
      sigma2 *= (1.+s_log_RutherfordScott*log(t));
   
    sigma2 *= (sigma2*t);
  }
 
  else {
   
    // Electron multiple scattering effects - see Highland NIM 129 (1975) p497-499
    // (Highland extension to the Rossi-Greisen formulation)
    // NOTE: The formula can be extended by replacing the momentum^2 term with pi * pf
    sigma2 = s_main_RossiGreisen / ( beta * p );
    sigma2 *= (sigma2*t);
   
    if ( m_config.log_include ) {
      double factor = 1. + s_log_RossiGreisen * log10( 10. * t );
      factor *= factor;
      sigma2 *= factor;
    }
  }
 
  // d_0'
  double dprime          = t/(beta*beta);
  double log_dprime      = log(dprime);
  // d_0''
  double log_dprimeprime = log(std::pow(mat.averageZ(),2.0/3.0) * dprime);
  // get epsilon
  double epsilon = log_dprimeprime < 0.5 ?
    s_gausMixEpsilon_a0 + s_gausMixEpsilon_a1*log_dprimeprime + s_gausMixEpsilon_a2*log_dprimeprime*log_dprimeprime :
    s_gausMixEpsilon_b0 + s_gausMixEpsilon_b1*log_dprimeprime + s_gausMixEpsilon_b2*log_dprimeprime*log_dprimeprime;
  // the standard sigma
  double sigma1square = s_gausMixSigma1_a0 + s_gausMixSigma1_a1*log_dprime + s_gausMixSigma1_a2*log_dprime*log_dprime;
  // G4 optimised / native double Gaussian model
  if (!m_config.optGaussianMixtureG4) sigma2 = 225.*dprime/(p*p);
  // throw the random number core/tail
  if ( m_config.randomNumbers->draw(Acts::Flat) < epsilon) {
    sigma2 *= (1.-(1.-epsilon)*sigma1square)/epsilon;
  }
 
  return s_projectionFactor*sqrt(sigma2)*m_config.randomNumbers->draw(Fatras::GaussZiggurat);
 
}