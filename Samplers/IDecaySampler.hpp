//////////////////////////////////////////////////////////////////
// IDecaySampler, Acts project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_FATRAS_IBARCODESAMPLER_H
#ifndef ACTS_FATRAS_IBARCODESAMPLER_H

namespace Fatras {
    
    /** @class IDecaySampler
     
     interface class for decay and free path length sampling in Fatras
     
     @author Andreas.Salzburger@cern.ch */
    
    class IDecaySampler {
        
      public:
        
        /** Virtual destructor */
        virtual IDecaySampler() {}
        
        /** sample the free path length */
        virtual double freePathLength(const ParticleProperties& pProperties) const = 0;
        
        /** create the actual decay of the particle - return interaction vertex */
        virtual double Acts::ProcessVertex decay(const ParticleProperties& pProperties) const = 0;
        
    };
    
    
}

#endif /* ACTS_FATRAS_IBARCODESAMPLER_H */
