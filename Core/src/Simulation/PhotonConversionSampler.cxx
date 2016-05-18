//////////////////////////////////////////////////////////////
// PhotonConversionSampler.cxx, ACTS project
///////////////////////////////////////////////////////////////////

// FATRAS includes
#include "FATRAS/Simulation/PhotonConversionSampler.h"
#include "FATRAS/Simulation/detail/FatrasDefinitions.h"
#include "FATRAS/Simulation/detail/FatrasHelpers.h"
// ACTS includes
#include "ACTS/EventData/ParticleDefinitions.h"
#include "ACTS/Utilities/MsgMacros.h"

// constructor
Fatras::PhotonConversionSampler::PhotonConversionSampler(const PhotonConversionSampler::Config& pcConfig)
: m_config()
{
    setConfiguration(pcConfig);
}

// destructor
Fatras::PhotonConversionSampler::~PhotonConversionSampler()
{}

void Fatras::PhotonConversionSampler::setConfiguration(const Fatras::PhotonConversionSampler::Config& msConfig ) 
{
  //!< @TODO update to configuration checking
  m_config = msConfig;   
}


/** processing the conversion*/
std::vector<Acts::InteractionVertex> Fatras::PhotonConversionSampler::doConversion(double time,
										                                           const Acts::Vector3D& position, 
										                                           const Acts::Vector3D& momentum) const
{
  double p = momentum.mag();

  // get the energy
  double childEnergy = p*childEnergyFraction(p);
  
  // now get the deflection
  Acts::Vector3D childDir(childDirection(momentum, childEnergy));
  
  // verbose output
  MSG_VERBOSE(  "[ conv ] Child energy simulated as : " << childEnergy << " MeV" );
  
  // calculate the second child direction and return
  return getChildren(time, position, momentum, childEnergy, childDir, Acts::electron);

}


double Fatras::PhotonConversionSampler::childEnergyFraction(double gammaMom) const {

  // @TODO write documentation
  // the fraction
  double epsilon0      = s_particleMasses.mass[Acts::electron]/gammaMom;
  // some needed manipolations
  double Z             = 13.;
  double oneOverZpow   = 1./pow(Z,s_oneOverThree);
  double alphaZsquare  = (s_alpha*s_alpha*Z*Z);
  // now f(Z) - from Z and s_alpha
  double fZ            = alphaZsquare*(1./(1.+alphaZsquare)+0.20206-0.0369*alphaZsquare+0.0083*alphaZsquare*alphaZsquare);
  // delta_max
  double deltaMax      = exp((42.24-fZ)*.1195)-0.952;
  // delta_min
  double deltaMin      = 4.*epsilon0*136.*oneOverZpow; 
  // the minimum fraction
  double epsilon1      = 0.5-0.5*sqrt(1.-deltaMin/deltaMax);
  double epsilonMin    = epsilon1 > epsilon0 ? epsilon1 : epsilon0;
  // calculate phi1 / phi2 - calculate from deltaMin
  double Phi1          = phi1(deltaMin);
  double Phi2          = phi2(deltaMin);
  // then calculate F10/F20
  double F10           = 3.*Phi1 - Phi2 - fZ;
  double F20           = 1.5*Phi1 - 0.5*Phi2 - fZ;
  // and finally calucate N1, N2
  double N1            = (0.25-epsilonMin+epsilonMin*epsilonMin)*F10;
  double N2            = 1.5*F20;
  // ------------ decide wich one to take 
  if ( N1/(N1+N2) < m_config.randomNumbers->draw(Fatras::Flat) ) {  
    // sample from f1,g1 distribution
    for ( ; ; ){
      double epsilon = 0.5 - (0.5 - epsilonMin)*pow(m_config.randomNumbers->draw(Fatras::Flat),s_oneOverThree);
      // prepare for the rejection check
      double delta   = 136.*epsilon0*oneOverZpow/(epsilon-epsilon*epsilon);
      double F1 = 3.*phi1(delta)-phi2(delta)-fZ;   
      // reject ? - or redo the exercise 
      if (F1/F10 > m_config.randomNumbers->draw(Fatras::Flat)) return m_config.childEnergyScaleFactor*epsilon;
    }
  } else {
    // sample from f2,g2 distribution
    for ( ; ; ){
      double epsilon = epsilonMin + (0.5-epsilonMin)*m_config.randomNumbers->draw(Fatras::Flat);
      // prepare for the rejection check
      double delta   = 136.*epsilon0*oneOverZpow/(epsilon-epsilon*epsilon);
      double F2 = 1.5*phi1(delta)-0.5*phi2(delta)-fZ;   
     // reject ? - or redo the exercise 
     if (F2/F20 > m_config.randomNumbers->draw(Fatras::Flat)) return m_config.childEnergyScaleFactor*epsilon;  
    }
  }

}

Acts::Vector3D Fatras::PhotonConversionSampler::childDirection(const Acts::Vector3D& gammaMom,
							     double childE) const
{
    // --------------------------------------------------
    // Following the Geant4 approximation from L. Urban
    // the azimutal angle
    double psi    =  2.*M_PI*m_config.randomNumbers->draw(Fatras::Flat);
    
    // the start of the equation
    double theta = s_particleMasses.mass[Acts::electron]/childE;
    // follow 
    double a = 0.625; // 5/8
    //double d = 27.;

    double r1 = m_config.randomNumbers->draw(Fatras::Flat);
    double r2 = m_config.randomNumbers->draw(Fatras::Flat);
    double r3 = m_config.randomNumbers->draw(Fatras::Flat);

    double u =  -log(r2*r3)/a;
    
    theta *= (r1 < 0.25 ) ? u : u*s_oneOverThree; // 9./(9.+27) = 0.25

     MSG_VERBOSE( "[ conv ] Simulated angle to photon    = " << theta << "." );

    // more complex but "more true"
    Acts::Vector3D newDirection(gammaMom.unit());
    double x = -newDirection.y();
    double y = newDirection.x();
    double z = 0.;
    // if it runs along the z axis - no good ==> take the x axis
    if (newDirection.z()*newDirection.z() > 0.999999)       
        x = 1.;
    // deflector direction
    //!>@TODO Check if this is right
    Acts::Vector3D deflector(x,y,z);
    // rotate the new direction for scattering using theta and arbitrarily in psi             
    // create the rotation
    Acts::RotationMatrix3D rotation;
    rotation = Acts::AngleAxis3D(theta, deflector)*Acts::AngleAxis3D(psi, gammaMom);
    // create the transform
    Acts::Transform3D transform(rotation, Acts::Vector3D(0., 0., 0.));
    // get the new direction
    newDirection = transform*newDirection;
    
    return newDirection;
}

std::vector<Acts::InteractionVertex> Fatras::PhotonConversionSampler::getChildren( double time,
                                                                               const Acts::Vector3D& vertex,
                                                                               const Acts::Vector3D& photonMomentum,
                                                                               double child1Energy, const Acts::Vector3D& child1Direction,
                                                                               Acts::ParticleType childType) const
{
    std::vector<Acts::InteractionVertex> children;
 
    // child 1 momentum
    double p1 = sqrt(child1Energy*child1Energy-s_particleMasses.mass[childType]*s_particleMasses.mass[childType]);    

    // now properly : energy-momentum conservation
    // child 2 momentum
    Acts::Vector3D child2Direction= (photonMomentum - p1*child1Direction).unit();
    double p2 = (photonMomentum - p1*child1Direction).mag();
    
    // charge sampling
    double charge1, charge2;
    charge1 = charge2 = 0.;
    if (m_config.randomNumbers->draw(Fatras::Flat)>0.5) {
      charge1 = -1.;
      charge2 =  1.;
    }
    else {
      charge1 =  1.;
      charge2 = -1.;
    }

    int    pdg1  = convertToPdg(childType, charge1, false);
    int    pdg2  = convertToPdg(childType, charge2, false);

    std::vector<Acts::ParticleProperties> pOutgoing;
    
    if (p1>m_config.minChildEnergy)
      pOutgoing.push_back(Acts::ParticleProperties(p1*child1Direction, pdg1));
    
    if (p2>m_config.minChildEnergy)
      pOutgoing.push_back(Acts::ParticleProperties(p2*child2Direction, pdg2));
    
    if (pOutgoing.size()>0)
      children.push_back(Acts::InteractionVertex(vertex, time, m_config.processCode, pOutgoing));

    return children;
}

