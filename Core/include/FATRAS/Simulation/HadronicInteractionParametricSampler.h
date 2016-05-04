///////////////////////////////////////////////////////////////////
// HadronicInteractionParametricSampler.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_FATRASTOOLS_HADRONICINTERACTIONPARAMETRICSAMPLER_H
#define ACTS_FATRASTOOLS_HADRONICINTERACTIONPARAMETRICSAMPLER_H 1

// FATRAS
#include "FATRAS/Common/IRandomNumbers.h"
#include "FATRAS/Simulation/IHadronicInteractionSampler.h"
// ACTS
#include "ACTS/Utilities/Definitions.h"
#include "ACTS/EventData/ParticleDefinitions.h"
// STD
#include <memory>
  
namespace Fatras {
  
  class InteractionVertex;

  /** @class HadronicInteractionParametricSampler
   * 
   * Parametric implementation of nuclear interactions to be used
   * in Fatras. The parameterisation is gathered from Geant4.
   * 
   * @author Andreas Salzburger <Andreas.Salzburger@cern.ch>
   * @author Carsten Magass     <Carsten.Magass@cern.ch>
   * @author Noemi Calace       <Noemi.Calace@cern.ch>
   *                 
   */
   
  class HadronicInteractionParametricSampler : virtual public IHadronicInteractionSampler {
  
  public:
    /** @struct Config   
     Configuration of this Samples*/
    struct Config {  
        
        std::shared_ptr<IRandomNumbers>   randomNumbers;          //!< Random Generator service */
        int                               processCode;            //!< process code */
        double                            minimumHadOutEnergy;    //!< hadronic interaction setting */
        bool                              cutChain;               //!< 
        
        Config() :
          randomNumbers(nullptr),
          processCode(1),
          minimumHadOutEnergy(100),
          cutChain(true)
        {}
        
    };       
    
    /** Constructor */
    HadronicInteractionParametricSampler(const Config& hiConfig);
    
    /** Destructor */
    virtual ~HadronicInteractionParametricSampler();

    /** processing of the presampled nuclear interactions on layer
     * This method returns the particle's children
     */
    std::vector<Acts::InteractionVertex> doHadronicInteraction(double time,
							                                   const Acts::Vector3D& position, 
							                                   const Acts::Vector3D& momentum,
							                                   Acts::ParticleType particle=Acts::pion) const final;
    
     /** Set configuration method */
     void setConfiguration(const Config& eeConfig);
   
     /** Get configuration method */
     Config getConfiguration() const;                                 
   
   protected:
     Config            m_config; //!< configuration object                      
         
   private:
     /** collect secondaries for layer material update */
     std::vector<Acts::InteractionVertex> getHadState( double time, double p,
						                               const Acts::Vector3D& vertex,
						                               const Acts::Vector3D& particleDir,
						                               Acts::ParticleType particle ) const;
						       
           
   };
   
   /** Return the configuration object */    
   inline HadronicInteractionParametricSampler::Config HadronicInteractionParametricSampler::getConfiguration() const { return m_config; }
}

#endif // ACTS_FATRASTOOLS_HADRONICINTERACTIONPARAMETRICSAMPLER_H