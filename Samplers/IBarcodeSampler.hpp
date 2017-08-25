//////////////////////////////////////////////////////////////////
// IBarcodeSampler.h, Acts project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_FATRAS_IBARCODESAMPLER_H
#ifndef ACTS_FATRAS_IBARCODESAMPLER_H

namespace Fatras {
    
    /** @class IBarcodeSampler 
     
     interface class for Barcode creation in Fatras
     
     @author Andreas.Salzburger@cern.ch */
    
    class IBarcodeSampler {
        
      public:
        
        /** Virtual destructor */
        virtual ~IBarcodeSampler() {}
        
        
        
        /** interface for processing of the presampled conversion on layer*/
        virtual barcode_type generateChildBarcode(barcode_type mother, ) const;
        
    };
     
    
}

#endif /* ACTS_FATRAS_IBARCODESAMPLER_H */
