///////////////////////////////////////////////////////////////////
// MultipleScatteringSamplerHighland.cxx, ACTS project
///////////////////////////////////////////////////////////////////

// FATRAS includes
#include "FATRAS/Simulation/MultipleScatteringSamplerHighland.h"
// ACTS includes
#include "ACTS/EventData/ParticleProperties.h"

// static particle masses
Acts::ParticleMasses Fatras::MultipleScatteringSamplerHighland::s_particleMasses;

// static doubles
double Fatras::MultipleScatteringSamplerHighland::s_main_RutherfordScott = 13.6;
double Fatras::MultipleScatteringSamplerHighland::s_log_RutherfordScott  =  0.038;

double Fatras::MultipleScatteringSamplerHighland::s_main_RossiGreisen    = 17.5;
double Fatras::MultipleScatteringSamplerHighland::s_log_RossiGreisen     =  0.125;

double Fatras::MultipleScatteringSamplerHighland::s_projectionFactor  =  sqrt(2.);

// constructor
Fatras::MultipleScatteringSamplerHighland::MultipleScatteringSamplerHighland(const MultipleScatteringSamplerHighland::Config& msConfig)
: m_config()
{
    setConfiguration(msConfig);
}

// destructor
Fatras::MultipleScatteringSamplerHighland::~MultipleScatteringSamplerHighland()
{}


void Fatras::MultipleScatteringSamplerHighland::setConfiguration(const Fatras::MultipleScatteringSamplerHighland::Config& msConfig ) 
{
    //!< @TODO update to configuration checking
   m_config = msConfig;   
}

double Fatras::MultipleScatteringSamplerHighland::simTheta(const Acts::MaterialProperties& mat,
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
 
  double sigma = m_interactionFormulae.sigmaMS(t, p, beta);
  sigma2 = sigma*sigma;
 
  // Code below will not be used if the parameterization of ActsUtils is used
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
 
  return s_projectionFactor*sqrt(sigma2)*m_config.randomNumbers->draw(Fatras::GaussZiggurat);
 
}