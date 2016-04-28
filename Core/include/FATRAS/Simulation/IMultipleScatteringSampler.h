///////////////////////////////////////////////////////////////////
// IMultipleScatteringSampler.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_FATRASINTERFACES_IMULTIPLESCATTERINGSAMPLER_H
#define ACTS_FATRASINTERFACES_IMULTIPLESCATTERINGSAMPLER_H 1

// ACTS include
#include "ACTS/EventData/ParticleHypothesis.h"

namespace Acts {
    class MaterialProperties;
}

namespace Fatras {
 
  /** @class IMultipleScatteringSampler
   * 
   * Interface class IMultipleScatteringSampler
   * 
   * @author Andreas Salzburger <Andreas.Salzburger@cern.ch>
   * @author Noemi Calace       <Noemi.Calace@cern.ch>
   * 
   */
  
  class IMultipleScatteringSampler {

  public:
    /**Virtual destructor*/
    virtual ~IMultipleScatteringSampler(){}
    
    virtual double simTheta (const Acts::MaterialProperties& mat,
			                 double momentum,
			                 double pathcorrection,
			                 Acts::ParticleHypothesis particle = Acts::pion) const = 0;    
  };

}

#endif // ACTS_FATRASINTERFACES_IMULTIPLESCATTERINGSAMPLER_H
