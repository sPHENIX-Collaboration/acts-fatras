///////////////////////////////////////////////////////////////////
// MaterialInteractionEngine.icc, ACTS project
///////////////////////////////////////////////////////////////////
       
#include "ACTS/Layers/Layer.h"
#include "ACTS/Material/SurfaceMaterial.h"

template <class T> Acts::ExtrapolationCode Fatras::MaterialInteractionEngine::handleMaterialT( Acts::ExtrapolationCell<T>& eCell,
                                                                                               Acts::PropDirection dir,
                                                                                               Acts::MaterialUpdateStage matupstage) const
{
  // for readability
  const Acts::Surface* mSurface = eCell.materialSurface;
  const Acts::Layer*   mLayer   = eCell.leadLayer;
  // the Extrapolator made sure that the layer is the lead layer && the parameters are the lead parameters
  if (mSurface && mSurface->surfaceMaterial()) {
    EX_MSG_DEBUG( ++eCell.navigationStep, "layer",  mLayer->geoID().value(), "handleMaterial for neutral parameters called - collect material.");
    // path correction             
    double pathCorrection = mSurface->pathCorrection(eCell.leadParameters->position(),dir*(eCell.leadParameters->momentum()));
    // the relative direction wrt with the layer
    Acts::PropDirection rlDir = (pathCorrection > 0. ? Acts::alongMomentum : Acts::oppositeMomentum);
    // multiply by the pre-and post-update factor
    double mFactor = mSurface->surfaceMaterial()->factor(rlDir, matupstage);
    if (mFactor == 0.){
      EX_MSG_VERBOSE(eCell.navigationStep, "layer",  mLayer->geoID().value(), "material collection with "  << (matupstage > 0. ? "pre " : "post ")  << "factor 0.");
      // return the parameters untouched
      return Acts::ExtrapolationCode::InProgress;
    }
    pathCorrection = mFactor*pathCorrection;
    // screen output
    EX_MSG_VERBOSE(eCell.navigationStep, "layer",  mLayer->geoID().value(), "material update with corr factor = " << pathCorrection);
    // get the actual material bin
    const Acts::MaterialProperties* materialProperties = mSurface->surfaceMaterial()->material(eCell.leadParameters->position());
    // and let's check if there's acutally something to do
    // if the material is filled perform the material interaction
    if ( materialProperties && materialProperties->thicknessInX0()>0) {
        // the fraction for re-entry of particles
        float mFraction=0.; 
        // call the main method
        return processOnLayerT(eCell, dir, mFraction);
    }
  }
  return Acts::ExtrapolationCode::InProgress;
}

template <class T> Acts::ExtrapolationCode Acts::MaterialInteractionEngine::processOnLayerT( Acts::ExtrapolationCell<T>& eCell,
												                                             Acts::PropDirection dir,
                                                                                             const Acts::MaterialProperties& mprop,
                                                                                             double pathCorrection,
												                                             float& mFraction) const
{
  // for readability
  const Acts::Surface* mSurface = eCell.materialSurface;
  const Acts::Layer*   mLayer   = eCell.leadLayer;
  
  // get the material itself & its parameters
  const Acts::Material& material = mprop.material();
  double thicknessInX0           = mprop.thicknessInX0();
  double thicknessInL0           = mprop.thicknessInL0();
 
  // figure out if particle stopped in the layer and recalculate path limit
  // - the mFraction determines what's left when re-entering the surface within one pass
  bool doInteraction = false;
  float dX0 = (1.-mFraction)*pathCorrection*thicknessInX0;
  float dL0 = (1.-mFraction)*pathCorrection*thicknessInL0;
  
  // electromagnetic interaction @TODO check with ST and document
  if ( eCell.materialLimitX0 > 0. && eCell.materialProcess < 100 && 
        eCell.materialX0+dX0 >= eCell.materialLimitX0) {
        // the remaing path in X0   
        float x0rem = eCell.materialLimitX0 - eCell.materialX0;
        // calculate the remaining (to be passed by children) L0
        dL0 *= x0rem > 0. ? x0rem/dX0 : 1.;
        // the remaining material in dX0
        if ( x0rem > 0. ) dX0 = x0rem;
        // interaction to be done as material limit will be passed
        doInteraction = true;
  } 
  // hadronic interaction @TODO check with ST and document
  else if ( eCell.materialLimitL0 > 0 && eCell.materialProcess > 100 && 
            eCell.materialL0+dL0 >= eCell.materialLimitL0 ) {
        // the remaining potential LO for this layer        
        float l0rem = eCell.materialLimitX0 - eCell.materialL0;
        // calculate the remaining (to be passed by children) X0
        dX0 *= l0rem > 0. ? l0rem/dL0 : 1.;
        // the remaining material in L0
        if ( l0rem > 0.) dL0 = l0rem;
        // interaction to be done as material limit will be passed
        doInteraction = true;
  }
  
  // check if material filling was requested - this is mainly for validation
  if (eCell.checkConfigurationMode(Acts::ExtrapolationMode::CollectMaterial)) {
    EX_MSG_VERBOSE(eCell.navigationStep, "layer",  mLayer->geoID().value(), "collecting material of [t/X0] = " << thicknessInX0); 
    eCell.stepMaterial(*mSurface, mLayer, eCell.leadParameters->position(), (1.-mFraction)*pathCorrection, materialProperties);
  } else {
    // always just record the material  
    EX_MSG_VERBOSE(eCell.navigationStep, "layer",  mLayer->geoID().value(), "adding material of [t/X0] = " << thicknessInX0);
    eCell.addMaterial((1.-mFraction)*pathCorrection, materialProperties);
  }
  
  // update the leadig track parameters
  const T* = updateTrackParameters(*eCell.leadParameters,eCell,dir,dX0,pathCorrection,mFraction);
  
  
  eCell.leadParameters;
  
  if (eCell.leadParameters->momentum().mag() < m_minimumMomentum ) 
      return Acts::ExtrapolationCode::SuccessMaterialLimit;
  
  // Update the material fraction.
  //!>@TODO  This should be used when children re-interact on the same layer
  mFraction += dX0/pathCorrection/thicknessInX0; 
   
  if ( doInteraction ) {   // interaction with particle stopping
 
    // Interact in the layer and return the vector of InteractionVertex
    std::vector<Acts::InteractionVertex> vertices = interact(eCell, material);
    
    //!>@TODO Evaluate the remaining material for children interaction on the same layer
    // And propagating to the children
    
    std::vector<ParticleProperties> surviving;
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
  }
  
  return Acts::ExtrapolationCode::InProgress;   
  
}