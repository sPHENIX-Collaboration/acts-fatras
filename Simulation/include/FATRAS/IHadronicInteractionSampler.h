///////////////////////////////////////////////////////////////////
// IHadronicInteractionSampler.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_FATRASINTERFACES_IHADRONICINTERACTIONSAMPLER_H
#define ACTS_FATRASINTERFACES_IHADRONICINTERACTIONSAMPLER_H 1

// ACTS include
#include "ACTS/Utilities/Definitions.h"
#include "ACTS/EventData/ParticleHypothesis.h"

namespace Acts {
    class InteractionVertex;
}

namespace Fatras {
  
  /**
   * @class IHadronicInteractionSampler
   * Interface definition for the handling of nuclear/hadronic interactions,
   * to be used by the MC based material effects updater
   * 
   * @author Andreas Salzburger <Andreas.Salzburger@cern.ch>
   * @author Noemi Calace   <Noemi.Calace@cern.ch>
   * 
  */
    class IHadronicInteractionSampler {

  public:
    
    /** Virtual destructor */    
    virtual ~IHadronicInteractionSampler() {}
    
    /** interface for processing of the presampled nuclear interactions on layer*/
    virtual std::vector<Acts::InteractionVertex> doHadronicInteraction(double time, 
								                                       const Acts::Vector3D& position, 
								                                       const Acts::Vector3D& momentum, 
								                                       Acts::ParticleHypothesis particle=Acts::pion) const = 0;

  };

}

#endif // ACTS_FATRASINTERFACES_IHADRONICINTERACTIONSAMPLER_H
