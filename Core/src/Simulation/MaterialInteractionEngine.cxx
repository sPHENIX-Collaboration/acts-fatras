///////////////////////////////////////////////////////////////////
// MaterialInteractionEngine.cxx, ACTS project
///////////////////////////////////////////////////////////////////

// FATRAS includes
#include "FATRAS/Simulation/MaterialInteractionEngine.h"
#include "FATRAS/Simulation/EnergyLoss.h"
#include "FATRAS/Simulation/detail/FatrasDefinitions.h"

// ACTS includes
#include "ACTS/Layers/Layer.h"
#include "ACTS/Material/SurfaceMaterial.h"
// STL
#include <sstream>
#include <math.h>

// constructor
Fatras::MaterialInteractionEngine::MaterialInteractionEngine(const MaterialInteractionEngine::Config& miConfig):
  Acts::IMaterialEffectsEngine(),
  m_config(),
  m_projectionFactor(sqrt(2.)/2.)
{
  setConfiguration(miConfig);   
}    

// destructor
Fatras::MaterialInteractionEngine::~MaterialInteractionEngine()
{}

// set configuration
void Fatras::MaterialInteractionEngine::setConfiguration(const MaterialInteractionEngine::Config& miConfig)
{
    //!< @TODO introduce configuration checking
    m_config = miConfig;
}

/** neutral extrapolation */
Acts::ExtrapolationCode Fatras::MaterialInteractionEngine::handleMaterial(Acts::ExCellNeutral& eCell,
                                                                          Acts::PropDirection dir,
                                                                          Acts::MaterialUpdateStage matupstage) const
{
  EX_MSG_DEBUG(++eCell.navigationStep, "handleMaterial", "neut", "handleMaterial for neutral particle called.");
  return handleMaterialT<Acts::NeutralParameters> (eCell, dir, matupstage); 
}


/** charged extrapolation */
Acts::ExtrapolationCode Fatras::MaterialInteractionEngine::handleMaterial(Acts::ExCellCharged& eCell,
                                                                          Acts::PropDirection dir,
                                                                          Acts::MaterialUpdateStage matupstage) const
{
    EX_MSG_DEBUG(++eCell.navigationStep, "handleMaterial", "char", "handleMaterial for charge particle called.");
    return handleMaterialT<Acts::TrackParameters> (eCell, dir, matupstage);    
}


/** neutral extrapolation - returning the same parameters */
const Acts::NeutralParameters* Fatras::MaterialInteractionEngine::electroMagneticInteraction(const Acts::NeutralParameters& parameters,
											                                                 Acts::ExCellNeutral&,
											                                                 Acts::PropDirection,
                                                                                             const Acts::MaterialProperties&,
                                                                                             double, double, double) const
{
  // don't do anything, no EM physics for neutral particles 
  return (&parameters);
}

/** charged extrapolation */
const Acts::TrackParameters* Fatras::MaterialInteractionEngine::electroMagneticInteraction(const Acts::TrackParameters& parameters,
										                                                   Acts::ExCellCharged& eCell,
										                                                   Acts::PropDirection dir,
                                                                                           const Acts::MaterialProperties& mprop,
                                                                                           double dX0,
										                                                   double pathCorrection,
										                                                   double mFraction) const
{
  // for readability
  const Acts::Surface* mSurface = eCell.materialSurface;
  const Acts::Layer*   mLayer   = eCell.leadLayer;
 
  // get the kinematics
  double p    = parameters.momentum().mag();
  double newP = 0.;
  double m    = m_particleMasses.mass[eCell.particleType];
  double E    = sqrt(p*p+m*m);

  // and let's check if there's acutally something to do
  // radiation and ionization preceed the presampled interaction (if any)
  if (m_config.energyLossSampler || m_config.doMultipleScattering) {    
    double thicknessInX0          = mprop.thicknessInX0();
    // a simple cross-check if the parameters are the initial ones
    ActsVectorD<5>  uParameters = parameters.parameters();
    // energy loss emulation        
    if (m_config.energyLossSampler || (eCell.particleType==Acts::electron && m_config.elEnergyLossSampler) ) {
      // smeared/presampled energy loss
      Acts::EnergyLoss eloss = (eCell.particleType==Acts::electron && m_config.elEnergyLossSampler) ? 
                                m_config.elEnergyLossSampler->energyLoss(*materialProperties, p, dX0/thicknessInX0, dir, eCell.particleType) : 
                                m_config.energyLossSampler->energyLoss(*materialProperties, p, dX0/thicknessInX0, dir, eCell.particleType);
      // Electron case 
      if (eCell.particleType==Acts::electron && m_config.recordBremPhoton) {
	    // ionization update first @TODO check with Sharka what about the minimumMomentum ?
        newP = E + eloss.meanIoni() > m ? sqrt((E+eloss.meanIoni())*(E+eloss.meanIoni())-m*m) : 0.;
        // update QOP
        uParameters[Acts::eQOP] = parameters.charge()/newP;
        // radiation
        if (newP > m_config.minimumMomentum) // mFraction used to estimate material thickness for brem photons to be seen
            radiate(uParameters, eCell, dX0, mFraction, pathCorrection*thicknessInX0);   
        // save the actual radiation loss
        float nqOp = uParameters[Acts::eQOP];
        float radLoss = fabs(1./nqOp) - newP;
        eloss.update(0.,0.,radLoss-eloss.meanRad(),eloss.meanRad()-radLoss);       
      } else {
        // MIP case  
        // calculate the new momentum
        newP = E + eloss.deltaE() > m ? sqrt((E+eloss.deltaE())*(E+eloss.deltaE())-m*m) : 0.
        // update QOP
        uParameters[Acts::eQOP] = parameters.charge()/newP;
      }
      
      EX_MSG_VERBOSE(eCell.navigationStep, "layer",  mLayer->geoID().value(), "Energy Loss evaluation : E, deltaE:" << E << "," << eloss.deltaE() );
           
      // return a nullptr if we fall under the minimum cut
      if (newP == 0. || ( newP < m_conf.minParticleMomentum ||   ))

           
//!>@TODO Is this needed here?
//       if ( 1./fabs(uParameters[Acts::eQOP]) < m_config.minimumMomentum ) {
// 	EX_MSG_VERBOSE( eCell.navigationStep, "layer",  mLayer->geoID().value(), "momentum less than minimum momentum. Returning SuccessMaterialLimit");
// 	return Acts::ExtrapolationCode::SuccessMaterialLimit;
//       }
//    }   
   
    if ( m_config.doMultipleScattering && thicknessInX0>0 ) {
      // Evaluate simTheta for multiple scattering and update the track parameters
      double simTheta = m_config.multipleScatteringSampler->simTheta(*materialProperties, p, dX0/thicknessInX0, eCell.particleType);
      //do the update -> You need 2 evaluation of simTheta. The second one is used to calculate deltaphi in multipleScatteringUpdate
      multipleScatteringUpdate(*(eCell.leadParameters), uParameters, simTheta,
                               m_config.multipleScatteringSampler->simTheta(*materialProperties, p, dX0/thicknessInX0, eCell.particleType));
     
    }
    
    // now either create new ones or update - only start parameters can not be updated 
    if (eCell.leadParameters != eCell.startParameters ){
      EX_MSG_VERBOSE(eCell.navigationStep, "layer",  mLayer->geoID().value(), "material update on non-initial parameters.");
      //!>@TODO how to update parameters ?!?
      // parameters.updateParameters(uParameters,uCovariance);
    } else {
      EX_MSG_VERBOSE(eCell.navigationStep, "layer",  mLayer->geoID().value(), "material update on initial parameters, creating new ones.");
      // create new parameters
      const Acts::Surface& tSurface = parameters.associatedSurface();
      const Acts::TrackParameters* tParameters = new Acts::BoundParameters(std::move(uCovariance),uParameters,tSurface);
      // these are newly created
      return tParameters;
    }
  }
  return (&parameters);
  
}


//// interaction for neutral particles //
//std::vector<Acts::InteractionVertex> Fatras::MaterialInteractionEngine::interact(Acts::ExCellNeutral& eCell, 
//										 const Acts::Material&) const
//{
//  // vector of children
//  std::vector<Acts::InteractionVertex> children;
//  
//  // get the process from the cell
//  int process = eCell.materialProcess;
//  if ( process==0 ) return children;
//  
//  // get the ParticleType
//  Acts::ParticleType particle = eCell.particleType;
//  
//  // get the position and the momentum from the parameters
//  const Acts::NeutralParameters* parm=eCell.leadParameters;
//  Acts::Vector3D position=parm->position();
//  Acts::Vector3D momentum=parm->momentum();
//
//  if (m_config.doConversion and process==14) { // photon conversion
//    
//    children = m_config.conversionSampler->doConversion(eCell.time, position, momentum); 
//
//  } else if (m_config.doHadronicInteraction and process==121) { // hadronic interaction
//
//    children =  m_config.hadronicInteractionSampler->doHadronicInteraction(eCell.time, position, momentum, particle);
//
//  } else  if (m_config.doDecay and process == 201 ) { // decay  
//    //!>@TODO Implement decay sampler
//    MSG_DEBUG("Decay not implemented yet");
//      
//  }
//
//  return children; 
//  
//}
//
//// interaction for charged particles //
//std::vector<Acts::InteractionVertex> Fatras::MaterialInteractionEngine::interact(Acts::ExCellCharged& eCell, 
//										 const Acts::Material&) const
//{
//  // vector of children
//  std::vector<Acts::InteractionVertex> children;
//  
//  // get the process from the cell
//  int process = eCell.materialProcess;
//  if ( process==0 ) return children;
//  
//  // get the ParticleType
//  Acts::ParticleType particle = eCell.particleType;
//  
//  // get the position and the momentum from the parameters
//  const Acts::TrackParameters* parm=eCell.leadParameters;
//  Acts::Vector3D position=parm->position();
//  Acts::Vector3D momentum=parm->momentum();
//  
//  if (m_config.doPositronAnnihilation and process == 5 ) {     // positron annihilation
//
//    double mass = m_config.particleMasses.mass[particle];   
//    double fmin = mass/momentum.mag();     
//    double fr = 1.-pow(fmin,randomNumbers->draw(Fatras::Flat));
//    
//    std::vector<ParticleProperties> pOutgoing = { ParticleProperties(fr*momentum, 22),
//                                                  ParticleProperties((1.-fr)*momentum, 22)};
//    
//    children.push_back(Acts::InteractionVertex(position, eCell.time, process, pOutgoing));
//    
//  } else if (m_config.doHadronicInteraction and process==121) {    // hadronic interaction
//
//    children =  m_config.hadronicInteractionSampler->doHadronicInteraction(eCell.time, position, momentum, particle);
//
//  } else  if (m_config.doDecay and process == 201 ) { // decay  
//    //!>@TODO Implement decay sampler
//    MSG_DEBUG("Decay not implemented yet");
//      
//  }
//
//  return children; 
//  
//}
//
//// updating parameters with multiple scattering effects
//void Fatras::MaterialInteractionEngine::multipleScatteringUpdate(const Acts::TrackParameters& pars,
//							         ActsVectorD<5> & parameters ,
//							         double simTheta, double num_config.deltaPhi) const
//{   
//  // parametric scattering - independent in x/y
//  if (m_config.parametricScattering){
//    EX_MSG_VERBOSE("[msupdate]", "MultipleScatteringUpdate", "", "Using parametric scattering." );
//    // the initial values
//    double theta =  parameters[Acts::eTHETA];
//    double phi   =  parameters[Acts::ePHI];
//    double sinTheta   = (sin(theta)*sin(theta) > 10e-10) ? sin(theta) : 1.; 
//   
//    // sample them in an independent way
//    double deltaTheta = m_config.projectionFactor*simTheta;
//    double deltaPhi   = m_config.projectionFactor*num_config.deltaPhi/sinTheta;
//   
//    phi += deltaPhi;
//    if (phi >= M_PI) phi -= M_PI;
//    else if (phi < -M_PI) phi += M_PI;
//    if (theta > M_PI) theta -= M_PI;
//   
//    EX_MSG_VERBOSE("[msupdate]", "MultipleScatteringUpdate", "", "deltaPhi / deltaTheta = " << deltaPhi << " / " << deltaTheta );
//   
//    // assign the new values
//    parameters[Acts::ePHI]   = phi;   
//    parameters[Acts::eTHETA] = fabs(theta + deltaTheta);
//
//  } else {
//    double thetaMs = simTheta;
//    double psi     = 2.*M_PI*randomNumbers->draw(Fatras::Flat);
//    // more complex but "more true"
//    Acts::Vector3D newDirection(pars.momentum().unit());
//    double x = -newDirection.y();
//    double y = newDirection.x();
//    double z = 0.;
//    // if it runs along the z axis - no good ==> take the x axis
//    if (newDirection.z()*newDirection.z() > 0.999999)       
//        x = 1.; y=0.;
//    // deflector direction
//    //!>@TODO Check if this is right
//    Acts::Vector3D deflector(x,y,z);
//    // rotate the new direction for scattering using theta and arbitrarily in psi             
//    // create the rotation
//    Acts::RotationMatrix3D rotation;
//    rotation = Acts::AngleAxis3D(thetaMs, deflector)*Acts::AngleAxis3D(psi, pars.momentum().unit());
//    EX_MSG_VERBOSE("[msupdate]", "MultipleScatteringUpdate", "", "deltaPsi / deltaTheta = " << psi << " / " << thetaMs );
//    // create the transform
//    Acts::Transform3D transform(rotation, Acts::Vector3D(0., 0., 0.));
//    // get the new direction
//    newDirection = transform*newDirection;
//    // assign the new values
//    parameters[Acts::ePHI]   = newDirection.phi();   
//    parameters[Acts::eTHETA] = newDirection.theta();
//  }
//
//}
//


// radiative effects
void Fatras::MaterialInteractionEngine::radiate( ActsVectorD<5> & parm ,
						                         Acts::ExCellCharged& eCell, 
						                         float pathLim, 
						                         float mFr,
						                         float refX) const 
{
  // sample energy loss and free path independently
  double path = 0.;
  double p = 1./ fabs(parm[Acts::eQOP]);
  
  Acts::Vector3D eDir = eCell.leadParameters->momentum().unit();
  Acts::Vector3D ePos = eCell.leadParameters->position();
  
  while ( path < pathLim && p>m_config.minimumMomentum ) {
  
    double rndx = randomNumbers->draw(Fatras::Flat);
  
    // sample visible fraction of the mother momentum taken according to 1/f 
    double eps = fmin(10.,m_config.minimumMomentum)/p;
  
    double z = pow(eps,pow(rndx,exp(1.1)));         
  
    // convert into scaling factor for mother momentum
    z = (1.- z);
  
    // turn into path   
    double dx = -0.7*log(z);     // adjust for mean of exp(-x) 
  
    // resolve the case when there is not enough material left in the layer
    if ( path+dx > pathLim ) {
      double rndp = randomNumbers->draw(Fatras::Flat);
      if (rndp > (pathLim-path)/dx){       
        (parm)[Acts::eQOP] = (parm)[Acts::eQOP]>0 ? 1/p : -1/p;
        mFr += (pathLim-path)/refX;
        path = pathLim;
        break;                   // radiation loop ended         
      }
      path += dx*rndp;
      mFr += dx*rndp/refX;
     
    } else {
      path+=dx;
      mFr += dx/refX;
    }
    if ( p*(1-z) > m_config.minimumBremPhotonMomentum ) {
  
      double deltaP = (1-z)*p;
      collectBremPhoton(eCell,p,deltaP,ePos,eDir);
      p *=z ;
  
      EX_MSG_VERBOSE("", "radiate", "", "brem photon emitted " << deltaP<<":updated e momentum:"<< p   );
    }   
  }
  
  parm[Acts::eQOP]    = (parm)[Acts::eQOP] > 0 ? 1/p : -1./p;
  parm[Acts::eTHETA]  = eDir.theta();
  parm[Acts::ePHI]    = eDir.phi();
  // 
  return;
}


//// generating BremPhoton
//void Fatras::MaterialInteractionEngine::collectBremPhoton(Acts::ExCellCharged& eCell,
//						          double pElectron,
//						          double gammaE,
//						          const Acts::Vector3D& vertex,
//                                                          Acts::Vector3D& particleDir) const
//{
//  // ------------------------------------------------------
//  // simple approach
//  // (a) simulate theta uniform within the opening angle of the relativistic Hertz dipole
//  //      theta_max = 1/gamma
//  // (b)Following the Geant4 approximation from L. Urban -> encapsulate that later
//  //      the azimutal angle
// 
//  double psi    =  2.*M_PI*randomNumbers->draw(Fatras::Flat);
// 
//  // the start of the equation
//  double theta = 0.;
//  double m = m_config.particleMasses.mass[eCell.particleType];
//  double E = sqrt(pElectron*pElectron + m*m);
//  
//  if (m_config.uniformHertzDipoleAngle) {
//    // the simplest simulation
//    theta = m/E * randomNumbers->draw(Fatras::Flat); 
//  } else {
//    // ----->
//    theta = m/E;
//    // follow
//    double a = 0.625; // 5/8
// 
//    double r1 = randomNumbers->draw(Fatras::Flat);
//    double r2 = randomNumbers->draw(Fatras::Flat);
//    double r3 = randomNumbers->draw(Fatras::Flat);
// 
//    double u =  -log(r2*r3)/a;
// 
//    theta *= (r1 < 0.25 ) ? u : u*m_config.oneOverThree; // 9./(9.+27) = 0.25
//  }
//
//  EX_MSG_VERBOSE("[ brem ]", "BremPhoton", "", "Simulated angle to electron    = " << theta << "." );
//
//  double th = particleDir.theta()-theta;
//  double ph = particleDir.phi();
//  if ( th<0.) { th *=-1; ph += M_PI; }
//  
//  Acts::Vector3D newDirection(sin(th)*cos(ph),sin(th)*sin(ph),cos(th));
//  //!>@TODO Check if this is right
//  // rotate the new direction for scattering using theta and arbitrarily in psi             
//  // create the rotation
//  Acts::RotationMatrix3D rotation;
//  rotation = Acts::AngleAxis3D(psi, particleDir);
//  // create the transform
//  Acts::Transform3D transform(rotation, Acts::Vector3D(0., 0., 0.));
//  // get the new direction
//  newDirection = transform*newDirection;
//  
//  // recoil / save input momentum for validation
//  Acts::Vector3D inEl(pElectron*particleDir);   
//  particleDir = (particleDir*pElectron- gammaE*newDirection).unit();
//  
//  std::vector<ParticleProperties> pOutgoing = {ParticleProperties(gammaE*newDirection, 22)};
//  eCell.interactionVertices.push_back(Acts::InteractionVertex(vertex, eCell.time, m_config.bremProcessCode, pOutgoing));
//  
//  //!>@TODO Evaluate the remaining material for children interaction on the same layer
//  // And propagating to the children
//}
