///////////////////////////////////////////////////////////////////
// MaterialInteractionEngine.icc, ACTS project
///////////////////////////////////////////////////////////////////
       
#include "ACTS/Layers/Layer.hpp"
#include "ACTS/Material/SurfaceMaterial.hpp"
#include "Fatras/MaterialInteractionEngine.hpp"
#include <math.h>
#include <sstream>
#include "ACTS/Layers/Layer.hpp"
#include "ACTS/Material/SurfaceMaterial.hpp"
#include "Fatras/EnergyLoss.hpp"
#include "Fatras/detail/FatrasDefinitions.hpp"
#include "Fatras/RandomNumberDistributions.hpp"

// constructor
template <class RandomGenerator>
Fatras::MaterialInteractionEngine<RandomGenerator>::MaterialInteractionEngine(
    const MaterialInteractionEngine::Config& miConfig,
    std::unique_ptr<const Acts::Logger> logger)
    : Acts::IMaterialEffectsEngine(),
      m_config(),
      m_logger(std::move(logger)),
      m_randomGenerator(nullptr),
      m_projectionFactor(sqrt(2.) / 2.),
      m_particleMasses()
       {
  setConfiguration(miConfig);
}

// set configuration
template <class RandomGenerator>
void Fatras::MaterialInteractionEngine<RandomGenerator>::setConfiguration(
    const MaterialInteractionEngine::Config& miConfig) {
  //!< @todo introduce configuration checking
  m_config = miConfig;
}

template <class RandomGenerator>
void Fatras::MaterialInteractionEngine<RandomGenerator>::setLogger(
    std::unique_ptr<const Acts::Logger> newLogger) {
  m_logger = std::move(newLogger);
}

template <class RandomGenerator>
void Fatras::MaterialInteractionEngine<RandomGenerator>::setRandomGenerator(RandomGenerator& randomGenerator){
    m_randomGenerator = &randomGenerator;
}

/** neutral extrapolation */
template <class RandomGenerator>
Acts::ExtrapolationCode
Fatras::MaterialInteractionEngine<RandomGenerator>::handleMaterial(Acts::ExCellNeutral& eCell, const Acts::Surface* msurface, Acts::PropDirection dir,
    Acts::MaterialUpdateStage matupstage) const {
  EX_MSG_DEBUG(++eCell.navigationStep, "handleMaterial", "neut",
               "handleMaterial for neutral particle called.");
  return handleMaterialT<Acts::NeutralParameters>(eCell, msurface, dir, matupstage);
}

/** charged extrapolation */
template <class RandomGenerator>
Acts::ExtrapolationCode
Fatras::MaterialInteractionEngine<RandomGenerator>::handleMaterial(Acts::ExCellCharged& eCell, const Acts::Surface* msurface, Acts::PropDirection dir,
    Acts::MaterialUpdateStage matupstage) const {
  EX_MSG_DEBUG(++eCell.navigationStep, "handleMaterial", "char",
               "handleMaterial for charged particle called.");
  return handleMaterialT<Acts::TrackParameters>(eCell, msurface, dir, matupstage);
}

template <class RandomGenerator>
template <class T> Acts::ExtrapolationCode Fatras::MaterialInteractionEngine<RandomGenerator>::handleMaterialT(Acts::ExtrapolationCell<T>& eCell,
                                                                                                             const Acts::Surface* msurface,
                                                                                                             Acts::PropDirection dir,
                                                                                                             Acts::MaterialUpdateStage matupstage) const
{
    
    // for readability
    // parameters are the lead parameters
    // by definition the material surface is the one the parameters are on
    const Acts::Surface& mSurface = msurface ? (*msurface) : eCell.leadParameters->referenceSurface();
    size_t approachID  = mSurface.geoID().value(Acts::GeometryID::approach_mask);
    size_t sensitiveID = mSurface.geoID().value(Acts::GeometryID::sensitive_mask);
    // approach of sensitive
    std::string surfaceType = sensitiveID ? "sensitive" : "surface";
    size_t      surfaceID   = sensitiveID ? sensitiveID : approachID;
    if (!m_randomGenerator) {
        EX_MSG_FATAL(++eCell.navigationStep,
                                         surfaceType,
                                         surfaceID,
                                         "Random generator not set! Please set it using setRandomGenerator() and reset it for every execution!");
        return Acts::ExtrapolationCode::FailureConfiguration;
    }
    // the Extrapolator made sure that the layer is the lead layer && the parameters are the lead parameters
    if (mSurface.associatedMaterial()) {
        EX_MSG_DEBUG(
                     ++eCell.navigationStep,
                     surfaceType,
                     surfaceID,
                     "handleMaterial called - collect material.");
        // path correction - the length of material seen in the given direction
        double pathCorrection = mSurface.pathCorrection(eCell.leadParameters->position(),dir*(eCell.leadParameters->momentum()));
        // the relative direction wrt to the layer
        Acts::PropDirection rlDir = (pathCorrection > 0. ? Acts::alongMomentum : Acts::oppositeMomentum);
        // multiply by the pre-and post-update factor
        double mFactor = mSurface.associatedMaterial()->factor(rlDir, matupstage);
        if (mFactor == 0.){
            EX_MSG_VERBOSE(eCell.navigationStep, surfaceType,  surfaceID, "material collection with "  << (matupstage > 0. ? "pre " : "post ")  << "factor 0.");
            // return the parameters untouched
            return Acts::ExtrapolationCode::InProgress;
        }
        pathCorrection *= mFactor;
        // screen output
        EX_MSG_VERBOSE(eCell.navigationStep, surfaceType,  surfaceID, "material update with corr factor = " << pathCorrection);
        // get the actual material bin
        const Acts::MaterialProperties* materialProperties = mSurface.associatedMaterial()->material(eCell.leadParameters->position());
        // and let's check if there's acutally something to do
        // if the material is filled perform the material interaction
        if ( materialProperties && materialProperties->thicknessInX0()>0) {
            // the fraction for re-entry of particles
            float mFraction=0.;
            // call the main method
            return processOnSurfaceT(eCell, msurface, dir, *materialProperties, pathCorrection, mFraction);
        }
    }
    return Acts::ExtrapolationCode::InProgress;
    
}

template <class RandomGenerator>
template <class T> Acts::ExtrapolationCode Fatras::MaterialInteractionEngine<RandomGenerator>::processOnSurfaceT(Acts::ExtrapolationCell<T>& eCell,                            const Acts::Surface* msurface,
                                                                                                             Acts::PropDirection dir,
                                                                                                             const Acts::MaterialProperties& mprop,
                                                                                                             double pathCorrection,
                                                                                                             float& mFraction) const
{
    // get the material itself & its parameters
    const Acts::Material& material = mprop.material();
    double thicknessInX0           = mprop.thicknessInX0();
    double thicknessInL0           = mprop.thicknessInL0();
    
    // electromagnetic interaction
  /*  std::cout << "!!Old parameters position: " << eCell.leadParameters->position().x() << "," << eCell.leadParameters->position().y() << "," << eCell.leadParameters->position().z() << std::endl;*/
    // @todo why parameters needed. when eCell handed over?
    auto newParameters = electroMagneticInteraction(*eCell.leadParameters, eCell, msurface,
                                                    dir, mprop, thicknessInX0,
                                                    pathCorrection, mFraction);
   // if (!newParameters) std::cout << "!newParameters" << std::endl;

  /*  std::cout << "!!Old parameters direction: " << eCell.leadParameters->momentum().x() << "," << eCell.leadParameters->momentum().y() << "," << eCell.leadParameters->momentum().z() << std::endl;
   
    std::cout << "!!New parameters direction: " << newParameters->momentum().x() << "," << newParameters->momentum().y() << "," << newParameters->momentum().z() << std::endl;
  */
    const Acts::Vector3D& stepPosition = newParameters->position();
    eCell.stepMaterial(std::move(newParameters),stepPosition,*msurface,pathCorrection,&mprop);
    //@todo could be only needed for output later
    /*
     // for readability
     // parameters are the lead parameters
     // by definition the material surface is the one the parameters are on
    const Acts::Surface& mSurface = msurface ? (*mSurface) : eCell.leadParameters->referenceSurface();
    size_t approachID  = mSurface.geoID().value(GeometryID::approach_mask);
    size_t sensitiveID = mSurface.geoID().value(GeometryID::sensitive_mask);
    // approach of sensitive
    std::string surfaceType = sensitiveID ? "sensitive" : "surface";
    size_t      surfaceID   = sensitiveID ? sensitiveID : approachID;
    */
    
    // figure out if particle stopped in the layer and recalculate path limit
    // - the mFraction determines what's left when re-entering the surface within one pass
   /* bool doInteraction = false;
    float dX0 = (1.-mFraction)*pathCorrection*thicknessInX0;
    float dL0 = (1.-mFraction)*pathCorrection*thicknessInL0;
    
    // electromagnetic interaction @todo check with ST and document
    // @todo: check: before check was eCell.materialProcess < 100 now eCell.interactionProcess < 100 - what does this actually mean?
    if ( eCell.materialLimitX0 > 0. && eCell.interactionProcess < 100 &&
        eCell.materialX0+dX0 >= eCell.materialLimitX0) {
        // @todo : whats this remaining material? Why are there material limits?
        // the remaing path in X0
        float x0rem = eCell.materialLimitX0 - eCell.materialX0;
        // calculate the remaining (to be passed by children) L0
        dX0 *= x0rem > 0. ? x0rem/dX0 : 1.;
        // the remaining material in dX0
        if ( x0rem > 0. ) dX0 = x0rem;
        // interaction to be done as material limit will be passed
        doInteraction = true;
    }
    // hadronic interaction @todo check with ST and document
    else if ( eCell.materialLimitL0 > 0 && eCell.materialProcess > 100 &&
             eCell.materialL0+dL0 >= eCell.materialLimitL0 ) {
        // the remaining potential LO for this layer
        float l0rem = eCell.materialLimitX0 - eCell.materialL0;
        // calculate the remaining (to be passed by children) X0
        dL0 *= l0rem > 0. ? l0rem/dL0 : 1.;
        // the remaining material in L0
        if ( l0rem > 0.) dL0 = l0rem;
        // interaction to be done as material limit will be passed
        doInteraction = true;
    }
    
    // check if material filling was requested - this is mainly for validation
    if (eCell.checkConfigurationMode(Acts::ExtrapolationMode::CollectMaterial)) {
        EX_MSG_VERBOSE(eCell.navigationStep, "layer",  mLayer->geoID().value(), "collecting material of [t/X0] = " << thicknessInX0);
        eCell.stepMaterial(mSurface, mLayer, eCell.leadParameters->position(), (1.-mFraction)*pathCorrection, mprop);
    } else {
        // always just record the material
        EX_MSG_VERBOSE(eCell.navigationStep, "layer",  mLayer->geoID().value(), "adding material of [t/X0] = " << thicknessInX0);
        eCell.addMaterial((1.-mFraction)*pathCorrection, mprop);
    }
    
    eCell.leadParameters;
    
    if (eCell.leadParameters->momentum().mag() < m_config.particleMinMomentum)
        return Acts::ExtrapolationCode::SuccessMaterialLimit;
    
    // Update the material fraction.
    //!>@todo  This should be used when children re-interact on the same layer
    mFraction += dX0/pathCorrection/thicknessInX0;
  */
    /*if ( doInteraction ) {   // interaction with particle stopping
     
     // Interact in the layer and return the vector of InteractionVertex
     std::vector<Acts::InteractionVertex> vertices = interact(eCell, material);
     
     //!>@todo Evaluate the remaining material for children interaction on the same layer
     // And propagating to the children
     
     std::vector<Acts::ParticleProperties> surviving;
     for (auto& vertex : vertices) {
     surviving.clear();
     for (auto& child : vertex.outgoingParticles()) {
     // if the momentum of the child is less than the minimum --> continue
     if (child.momentum().mag()>m_minimumMomentum) surviving.push_back(child);
     }
     if (surviving.size()>0)
     eCell.interactionVertices.push_back(Acts::InteractionVertex(vertex.vertex(), vertex.interactionTime(), vertex.interactionType(), surviving));
     }
     return Acts::ExtrapolationCode::SuccessMaterialLimit;
     }*/
    
    return Acts::ExtrapolationCode::InProgress;   
    
}



/** neutral extrapolation - returning the same parameters */
template <class RandomGenerator>
std::unique_ptr<const Acts::NeutralParameters>
Fatras::MaterialInteractionEngine<RandomGenerator>::electroMagneticInteraction(const Acts::NeutralParameters& parameters, Acts::ExCellNeutral&, const Acts::Surface* msurface,
    Acts::PropDirection, const Acts::MaterialProperties&, double, double,
    double) const {
  // don't do anything, no EM physics for neutral particles
    return (nullptr);//std::make_unique<const Acts::NeutralParameters>(&parameters)
}

/** charged extrapolation */
template <class RandomGenerator>
std::unique_ptr<const Acts::TrackParameters>
Fatras::MaterialInteractionEngine<RandomGenerator>::electroMagneticInteraction(const Acts::TrackParameters& parameters, Acts::ExCellCharged& eCell, const Acts::Surface* msurface,
    Acts::PropDirection dir, const Acts::MaterialProperties& mprop, double dX0,
    double pathCorrection, double mFraction) const {
  // for readability
    const Acts::Surface& mSurface = msurface ? (*msurface) : eCell.leadParameters->referenceSurface();
    size_t approachID  = mSurface.geoID().value(Acts::GeometryID::approach_mask);
    size_t sensitiveID = mSurface.geoID().value(Acts::GeometryID::sensitive_mask);
    // approach of sensitive
    std::string surfaceType = sensitiveID ? "sensitive" : "surface";
    size_t      surfaceID   = sensitiveID ? sensitiveID : approachID;
    
    // get the parameters
    double p = parameters.momentum().mag();
    double newP = 0.;
    double m = m_particleMasses.mass[eCell.particleType];
    double E = sqrt(p * p + m * m);
    // @todo only needed for msc or energyloss
    double thicknessInX0 = mprop.thicknessInX0();
    // the parameters to be returned
    Acts::ActsVectorD<5> uParameters = parameters.parameters();
    // the covariance to be returned
    // @todo use nGlobalPars not 5!
    std::unique_ptr<Acts::ActsSymMatrixD<5>> uCovariance
    = nullptr;
    
    // Multiple scattering
    if (m_config.multipleScatteringSampler && thicknessInX0 > 0){
        double simTheta = m_config.multipleScatteringSampler->simTheta(*m_randomGenerator,
                                                                       mprop, p, dX0 / thicknessInX0,
                                                                       eCell.particleType);
        // do the update -> You need 2 evaluation of simTheta. The second one is used to calculate deltaphi in multipleScatteringUpdate
        
        multipleScatteringUpdate(*(eCell.leadParameters), uParameters, simTheta,
                                 m_config.multipleScatteringSampler->simTheta(*m_randomGenerator,
                                                                              mprop, p, dX0 / thicknessInX0,
                                                                              eCell.particleType));

    }
    
    if (m_config.energyLossSampler || (eCell.particleType == Acts::electron &&
                                       m_config.energyLossSamplerElectrons)) {
        // smeared/presampled energy loss
        Fatras::EnergyLoss eloss =
        (eCell.particleType == Acts::electron &&
         m_config.energyLossSamplerElectrons)
        ? m_config.energyLossSamplerElectrons->energyLoss(*m_randomGenerator,
                                                   mprop, p, dX0 / thicknessInX0, dir,
                                                   eCell.particleType)
        : m_config.energyLossSampler->energyLoss(*m_randomGenerator,mprop, p,
                                                 dX0 / thicknessInX0,
                                                 dir,
                                                 eCell.particleType);
        //@todo add electron case with brem photon and radiation
        // for now make no distinction
        // MIP case
        // calculate the new momentum
       /* std::cout << "MaterialInteractionEngine::m: " << m << "Energy Loss evaluation : E, deltaE:" << E << ","
        << eloss.deltaE() << ", E + eloss.deltaE(): " << E + eloss.deltaE() << std::endl;*/
        newP = (E + eloss.deltaE()) > m
        ? sqrt((E + eloss.deltaE()) * (E + eloss.deltaE()) - m * m)
        : 0.;
     //   std::cout << "MaterialInteractionEngine::newP: " << newP << ", old p: " << p << std::endl;
        
     //   if ((float)newP == (float)p) std::cout << "MaterialInteractionEngine::newP == p!" << std::endl;
    //    if (newP == 0.) std::cout << "MaterialInteractionEngine::newP == 0.!" << std::endl;
        
        // update QOP
        uParameters[Acts::eQOP] = parameters.charge() / newP;
     //   if (fabs(eloss.deltaE()) < 0.001) std::cout << "MaterialInteractionEngine::ELoss smaller 0.001" << std::endl;
        
        EX_MSG_VERBOSE(eCell.navigationStep, surfaceType,  surfaceID, "Energy Loss evaluation : E, deltaE:" << E << ","
                       << eloss.deltaE());
    }

    
    // ------------------------------------------------------------------

    // get the kinematics
  /*  double p = parameters.momentum().mag();
    double newP = 0.;
    double m = m_particleMasses.mass[eCell.particleType];
    double E = sqrt(p * p + m * m);

    // and let's check if there's acutally something to do
    // radiation and ionization preceed the presampled interaction (if any)
    if (m_config.energyLossSampler || m_config.doMultipleScattering) {
      double thicknessInX0 = mprop.thicknessInX0();
      // a simple cross-check if the parameters are the initial ones
      ActsVectorD<5> uParameters = parameters.parameters();
      // energy loss emulation
      if (m_config.energyLossSampler || (eCell.particleType == Acts::electron &&
                                         m_config.elEnergyLossSampler)) {
        // smeared/presampled energy loss
        Acts::EnergyLoss eloss =
            (eCell.particleType == Acts::electron &&
    m_config.elEnergyLossSampler)
                ? m_config.elEnergyLossSampler->energyLoss(
                      *materialProperties, p, dX0 / thicknessInX0, dir,
                      eCell.particleType)
                : m_config.energyLossSampler->energyLoss(*materialProperties, p,
                                                         dX0 / thicknessInX0,
    dir,
                                                         eCell.particleType);
        // Electron case
        if (eCell.particleType == Acts::electron && m_config.recordBremPhoton) {
          // ionization update first @TODO check with Sharka what about the
          // minimumMomentum ?
          newP =
              E + eloss.meanIoni() > m
                  ? sqrt((E + eloss.meanIoni()) * (E + eloss.meanIoni()) - m *
    m)
                  : 0.;
          // update QOP
          uParameters[Acts::eQOP] = parameters.charge() / newP;
          // radiation
          if (newP > m_config.minimumMomentum)  // mFraction used to estimate
                                                // material thickness for brem
                                                // photons to be seen
            radiate(uParameters, eCell, dX0, mFraction,
                    pathCorrection * thicknessInX0);
          // save the actual radiation loss
          float nqOp = uParameters[Acts::eQOP];
          float radLoss = fabs(1. / nqOp) - newP;
          eloss.update(0., 0., radLoss - eloss.meanRad(),
                       eloss.meanRad() - radLoss);
        } else {
          // MIP case
          // calculate the new momentum
          newP = E + eloss.deltaE() > m
                     ? sqrt((E + eloss.deltaE()) * (E + eloss.deltaE()) - m * m)
                     : 0.
                       // update QOP
                       uParameters[Acts::eQOP] = parameters.charge() / newP;
        }

        EX_MSG_VERBOSE(eCell.navigationStep, "layer", mLayer->geoID().value(),
                       "Energy Loss evaluation : E, deltaE:" << E << ","
                                                             << eloss.deltaE());

        // return a nullptr if we fall under the minimum cut
        if (newP == 0. || (newP < m_conf.minParticleMomentum ||))

          //!>@TODO Is this needed here?
          //       if ( 1./fabs(uParameters[Acts::eQOP]) <
          //       m_config.minimumMomentum ) {
          // 	EX_MSG_VERBOSE( eCell.navigationStep, "layer",
          // mLayer->geoID().value(), "momentum less than minimum momentum.
          // Returning SuccessMaterialLimit");
          // 	return Acts::ExtrapolationCode::SuccessMaterialLimit;
          //       }
          //    }

          if (m_config.doMultipleScattering && thicknessInX0 > 0) {
            // Evaluate simTheta for multiple scattering and update the track
            // parameters
            double simTheta = m_config.multipleScatteringSampler->simTheta(
                *materialProperties, p, dX0 / thicknessInX0,
    eCell.particleType);
            // do the update -> You need 2 evaluation of simTheta. The second
    one
            // is used to calculate deltaphi in multipleScatteringUpdate
            multipleScatteringUpdate(
                *(eCell.leadParameters), uParameters, simTheta,
                m_config.multipleScatteringSampler->simTheta(
                    *materialProperties, p, dX0 / thicknessInX0,
                    eCell.particleType));
          }

        // now either create new ones or update - only start parameters can not
    be
        // updated
        if (eCell.leadParameters != eCell.startParameters) {
          EX_MSG_VERBOSE(eCell.navigationStep, "layer", mLayer->geoID().value(),
                         "material update on non-initial parameters.");
          //!>@TODO how to update parameters ?!?
          // parameters.updateParameters(uParameters,uCovariance);
        } else {
          EX_MSG_VERBOSE(
              eCell.navigationStep, "layer", mLayer->geoID().value(),
              "material update on initial parameters, creating new ones.");
          // create new parameters
          const Acts::Surface& tSurface = parameters.associatedSurface();
          const Acts::TrackParameters* tParameters = new Acts::BoundParameters(
              std::move(uCovariance), uParameters, tSurface);
          // these are newly created
          return tParameters;
        }
      }
      return (&parameters);*/
    std::unique_ptr<const Acts::TrackParameters> tParameters = std::make_unique<Acts::BoundParameters>(std::move(uCovariance), uParameters, *msurface);  // only for now
    return std::move(tParameters);
}

//// interaction for neutral particles //
// std::vector<Acts::InteractionVertex>
// Fatras::MaterialInteractionEngine::interact(Acts::ExCellNeutral& eCell,
//										 const
// Acts::Material&)
// const
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
//    children = m_config.conversionSampler->doConversion(eCell.time,
//    position, momentum);
//
//  } else if (m_config.doHadronicInteraction and process==121) { // hadronic
//  interaction
//
//    children =
//    m_config.hadronicInteractionSampler->doHadronicInteraction(eCell.time,
//    position, momentum, particle);
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
// std::vector<Acts::InteractionVertex>
// Fatras::MaterialInteractionEngine::interact(Acts::ExCellCharged& eCell,
//										 const
// Acts::Material&)
// const
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
//  if (m_config.doPositronAnnihilation and process == 5 ) {     // positron
//  annihilation
//
//    double mass = m_config.particleMasses.mass[particle];
//    double fmin = mass/momentum.mag();
//    double fr = 1.-pow(fmin,randomNumbers->draw(Fatras::Flat));
//
//    std::vector<ParticleProperties> pOutgoing = {
//    ParticleProperties(fr*momentum, 22),
//                                                  ParticleProperties((1.-fr)*momentum,
//                                                  22)};
//
//    children.push_back(Acts::InteractionVertex(position, eCell.time,
//    process, pOutgoing));
//
//  } else if (m_config.doHadronicInteraction and process==121) {    //
//  hadronic interaction
//
//    children =
//    m_config.hadronicInteractionSampler->doHadronicInteraction(eCell.time,
//    position, momentum, particle);
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

template <class RandomGenerator>
 void Fatras::MaterialInteractionEngine<RandomGenerator>::multipleScatteringUpdate(const Acts::TrackParameters& pars,
                                                                  Acts::ActsVectorD<5>& parameters,
                                                                  double simTheta,
                                                                  double num_deltaPhi) const
{
    
  // parametric scattering - independent in x/y
  if (m_config.parametricScattering){
    EX_MSG_VERBOSE("[msupdate]", "MultipleScatteringUpdate", "", "Using parametric scattering." );
    // the initial values
    double theta =  parameters[Acts::eTHETA];
    double phi   =  parameters[Acts::ePHI];
    double sinTheta   = (sin(theta)*sin(theta) > 10e-10) ? sin(theta) : 1.;

  //    std::cout << "old theta" << theta << std::endl;
    // @todo whats the projectionfactor?
    // sample them in an independent way
    double deltaTheta = m_projectionFactor*simTheta;
  //    std::cout << "projectionfactor: " << m_projectionFactor << " simTheta: " << simTheta << std::endl;
 //     std::cout << "deltaTheta: " << deltaTheta << std::endl;
    double deltaPhi   =
    m_projectionFactor*num_deltaPhi/sinTheta;

    phi += deltaPhi;
    if (phi >= M_PI) phi -= M_PI;
    else if (phi < -M_PI) phi += M_PI;
    if (theta > M_PI) theta -= M_PI;

    EX_MSG_VERBOSE("[msupdate]", "MultipleScatteringUpdate", "", "deltaPhi / deltaTheta = " << deltaPhi << " / " << deltaTheta );

    // assign the new values
    parameters[Acts::ePHI]   = phi;
    parameters[Acts::eTHETA] = fabs(theta + deltaTheta);
   
  } else {
    // Create a random uniform distribution between in the intervall [0,1]
    Fatras::UniformDist uniformDist(0., 1.);
    //@todo test this non parametric way - not tested yet
    double thetaMs = simTheta;
    double psi     = 2.*M_PI*uniformDist(*m_randomGenerator);
      
  //    std::cout << "MaterialInteractionEngine::thetaMs: " << thetaMs << std::endl;
  //    std::cout << "MaterialInteractionEngine::Psi: " << psi << std::endl;
    // more complex but "more true"
    Acts::Vector3D newDirection(pars.momentum().unit());
  //    std::cout << "MaterialInteractionEngine::Old Direction: " <<newDirection << std::endl;
    double x = -newDirection.y();
    double y = newDirection.x();
    double z = 0.;
    // if it runs along the z axis - no good ==> take the x axis
    if (newDirection.z()*newDirection.z() > 0.999999)
        x = 1.; y=0.;
    // deflector direction
    //!>@todo Check if this is right
    Acts::Vector3D deflector(x,y,z);
    // rotate the new direction for scattering using theta and arbitraril in psi
    // create the rotation
    Acts::RotationMatrix3D rotation;
    rotation = Acts::AngleAxis3D(thetaMs, deflector)*Acts::AngleAxis3D(psi,
    pars.momentum().unit());
    EX_MSG_VERBOSE("[msupdate]", "MultipleScatteringUpdate", "", "deltaPsi / deltaTheta = " << psi << " / " << thetaMs );
    // create the transform
    Acts::Transform3D transform(rotation, Acts::Vector3D(0., 0., 0.));
    // get the new direction
    newDirection = transform*newDirection;
 //     std::cout << "MaterialInteractionEngine::New Direction: " <<newDirection << std::endl;
    // assign the new values
    parameters[Acts::ePHI]   = newDirection.phi();
    parameters[Acts::eTHETA] = newDirection.theta();
  }
}


// radiative effects
/*void Fatras::MaterialInteractionEngine::radiate(ActsVectorD<5>& parm,
                                                Acts::ExCellCharged& eCell,
                                                float pathLim, float mFr,
                                                float refX) const {
  // sample energy loss and free path independently
  double path = 0.;
  double p = 1. / fabs(parm[Acts::eQOP]);

  Acts::Vector3D eDir = eCell.leadParameters->momentum().unit();
  Acts::Vector3D ePos = eCell.leadParameters->position();

  while (path < pathLim && p > m_config.minimumMomentum) {
    double rndx = randomNumbers->draw(Fatras::Flat);

    // sample visible fraction of the mother momentum taken according to 1/f
    double eps = fmin(10., m_config.minimumMomentum) / p;

    double z = pow(eps, pow(rndx, exp(1.1)));

    // convert into scaling factor for mother momentum
    z = (1. - z);

    // turn into path
    double dx = -0.7 * log(z);  // adjust for mean of exp(-x)

    // resolve the case when there is not enough material left in the layer
    if (path + dx > pathLim) {
      double rndp = randomNumbers->draw(Fatras::Flat);
      if (rndp > (pathLim - path) / dx) {
        (parm)[Acts::eQOP] = (parm)[Acts::eQOP] > 0 ? 1 / p : -1 / p;
        mFr += (pathLim - path) / refX;
        path = pathLim;
        break;  // radiation loop ended
      }
      path += dx * rndp;
      mFr += dx * rndp / refX;

    } else {
      path += dx;
      mFr += dx / refX;
    }
    if (p * (1 - z) > m_config.minimumBremPhotonMomentum) {
      double deltaP = (1 - z) * p;
      collectBremPhoton(eCell, p, deltaP, ePos, eDir);
      p *= z;

      EX_MSG_VERBOSE("", "radiate", "", "brem photon emitted "
                                            << deltaP
                                            << ":updated e momentum:" << p);
    }
  }

  parm[Acts::eQOP] = (parm)[Acts::eQOP] > 0 ? 1 / p : -1. / p;
  parm[Acts::eTHETA] = eDir.theta();
  parm[Acts::ePHI] = eDir.phi();
  //
  return;
}
*/
//// generating BremPhoton
// void
// Fatras::MaterialInteractionEngine::collectBremPhoton(Acts::ExCellCharged&
// eCell,
//						          double pElectron,
//						          double gammaE,
//						          const Acts::Vector3D&
// vertex,
//                                                          Acts::Vector3D&
//                                                          particleDir) const
//{
//  // ------------------------------------------------------
//  // simple approach
//  // (a) simulate theta uniform within the opening angle of the relativistic
//  Hertz dipole
//  //      theta_max = 1/gamma
//  // (b)Following the Geant4 approximation from L. Urban -> encapsulate that
//  later
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
//    theta *= (r1 < 0.25 ) ? u : u*m_config.oneOverThree; // 9./(9.+27) =
//    0.25
//  }
//
//  EX_MSG_VERBOSE("[ brem ]", "BremPhoton", "", "Simulated angle to electron
//  = " << theta << "." );
//
//  double th = particleDir.theta()-theta;
//  double ph = particleDir.phi();
//  if ( th<0.) { th *=-1; ph += M_PI; }
//
//  Acts::Vector3D newDirection(sin(th)*cos(ph),sin(th)*sin(ph),cos(th));
//  //!>@TODO Check if this is right
//  // rotate the new direction for scattering using theta and arbitrarily in
//  psi
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
//  std::vector<ParticleProperties> pOutgoing =
//  {ParticleProperties(gammaE*newDirection, 22)};
//  eCell.interactionVertices.push_back(Acts::InteractionVertex(vertex,
//  eCell.time, m_config.bremProcessCode, pOutgoing));
//
//  //!>@TODO Evaluate the remaining material for children interaction on the
//  same layer
//  // And propagating to the children
//}
