///////////////////////////////////////////////////////////////////
// EnergyLossSampler.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_FATRSTOOLS_ENERGYLOSSSAMPLER_H
#define ACTS_FATRSTOOLS_ENERGYLOSSSAMPLER_H 1

// ACTS include
#include "ACTS/Utilities/Definitions.h"
#include "ACTS/EventData/ParticleHypothesis.h"
#include "ACTS/Extrapolation/detail/MaterialInteraction.h"
#include "ACTS/Material/MaterialProperties.h"
// Fatras module
#include "FATRAS/IEnergyLossSampler.h"
#include "FATRAS/IRandomNumbers.h"
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
      EnergyLossSampler(const Congif& elConfig);
      
      /** Destructor */
      ~EnergyLossSampler();
      
      /** IEnergyLossSampler public method to compute dEdX */
      double dEdX( const Acts::MaterialProperties& materialProperties, 
                   double momentum,
	  	           Acts::ParticleHypothesis particleHypothesis = Acts::pion ) const final;
      
      /** IEnergyLossSampler public method to compute the mean and variance of the energy loss */
      Acts::EnergyLoss energyLoss( const Acts::MaterialProperties& mat,
	  			                   double momentum,
	  			                   double pathcorrection,
	  			                   Acts::PropDirection dir=Acts::alongMomentum,
	  			                   Acts::ParticleHypothesis particle=Acts::pion) const final;
    
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
	  		               Acts::ParticleHypothesis particle= Acts::pion) const;
	
      /** apply dEdX according to Bethe-Heitler */		  
      double dEdXBetheHeitler(const  Acts::MaterialProperties& mat,
	  		                  double initialE,
	  		                  Acts::ParticleHypothesis particle= Acts::pion) const;
      
      /** the formulas for energy loss evaluation */
      Acts::MaterialInteraction                  m_interactionFormulae;     
      
      /** struct of Particle masses  */
      static Acts::ParticleMasses                s_particleMasses; 
      
      /** KOverA factor in Bethe-Bloch equation [MeV*cm2/gram] */
      static double                              s_ka_BetheBloch;          
};

  /** Return the configuration object */    
  inline EnergyLossSampler::Config EnergyLossSampler::getConfiguration() const { return m_config; }

} 
#endif //  ACTS_FATRSTOOLS_ENERGYLOSSSAMPLER_H