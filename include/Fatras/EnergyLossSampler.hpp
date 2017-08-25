///////////////////////////////////////////////////////////////////
// EnergyLossSampler.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_FATRSTOOLS_ENERGYLOSSSAMPLER_H
#define ACTS_FATRSTOOLS_ENERGYLOSSSAMPLER_H 1

// ACTS include
#include "ACTS/Utilities/Definitions.h"
#include "ACTS/EventData/ParticleDefinitions.h"
#include "ACTS/Extrapolation/detail/MaterialInteraction.h"
#include "ACTS/Material/MaterialProperties.h"
// Fatras module
#include "FATRAS/Common/IRandomNumbers.h"
#include "FATRAS/Simulation/IEnergyLossSampler.h"
// STD
#include <memory>

namespace Fatras {
  
  class EnergyLoss;

  /** @class EnergyLossSampler
   * 
   * Sampler for a eloss of a track
   * It uses the Bethe-Bloch calculation
   * it extends the IEnergyLossSampler interface 
   * 
   * @author Andreas Salzburger <Andreas.Salzburger@cern.ch>
   * @author Noemi Calace       <Noemi.Calace@cern.ch>
   */

  class EnergyLossSampler : virtual public IEnergyLossSampler {
   
    public:
      /** @struct Config 
        Configuration config for the EnergyLossSampler
        */
      struct Config {
          /** Random Generator service  */
          std::shared_ptr<IRandomNumbers>      randomNumbers;  
      };      
   
      /** Constructor with AlgTool parameters */
      EnergyLossSampler(const Config& elConfig);
      
      /** Destructor */
      ~EnergyLossSampler();
      
      /** IEnergyLossSampler public method to compute dEdX */
      double dEdX( const Acts::MaterialProperties& materialProperties, 
                   double momentum,
	  	           Acts::ParticleType particleHypothesis = Acts::pion ) const final;
      
      /** IEnergyLossSampler public method to compute the mean and variance of the energy loss */
      EnergyLoss energyLoss( const Acts::MaterialProperties& mat,
	  	                   double momentum,
	  	                   double pathcorrection,
	  	                   Acts::PropDirection dir=Acts::alongMomentum,
	  	                   Acts::ParticleType particle=Acts::pion) const final;
    
      /** Set configuration method */
      void setConfiguration(const Config& eeConfig);
    
      /** Get configuration method */
      Config getConfiguration() const;                                 
    
    protected:
      Config            m_config; //!< configuration object                      
         
    private:
      /** apply dEdX according to Bethe-Bloch */
      double dEdXBetheBloch(const Acts::MaterialProperties& mat,
	  		               double gamma,
	  		               double beta,
	  		               Acts::ParticleType particle= Acts::pion) const;
	
      /** apply dEdX according to Bethe-Heitler */		  
      double dEdXBetheHeitler(const  Acts::MaterialProperties& mat,
	  		                  double initialE,
	  		                  Acts::ParticleType particle= Acts::pion) const;
      
       
};

  /** Return the configuration object */    
  inline EnergyLossSampler::Config EnergyLossSampler::getConfiguration() const { return m_config; }

} 
#endif //  ACTS_FATRSTOOLS_ENERGYLOSSSAMPLER_H