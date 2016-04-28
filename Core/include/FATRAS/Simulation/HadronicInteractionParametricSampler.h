///////////////////////////////////////////////////////////////////
// HadronicInteractionParametricSampler.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_FATRASTOOLS_HADRONICINTERACTIONPARAMETRICSAMPLER_H
#define ACTS_FATRASTOOLS_HADRONICINTERACTIONPARAMETRICSAMPLER_H 1

// FATRAS
#include "FATRAS/Common/IRandomNumbers.h"
#include "FATRAS/Simulation/IHadronicInteractionSampler.h"
#include "FATRAS/Simulation/detail/PdgToParticleHypothesis..h"
// ACTS
#include "ACTS/Utilities/Definitions.h"
#include "ACTS/EventData/ParticleHypothesis.h"
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
        
        ServiceHandle<Acts::IRandomNumbers>     randomNumbers;          //!< Random Generator service */
        int                                     processCode;            //!< process code */
        double                                  minimumHadOutEnergy;    //!< hadronic interaction setting */
        bool                                    cutChain;               //!< 
        
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
							                                   Acts::ParticleHypothesis particle=Acts::pion) const final;
 
   private:
     /** collect secondaries for layer material update */
     std::vector<Acts::InteractionVertex> getHadState( double time, double p,
						                               const Acts::Vector3D& vertex,
						                               const Acts::Vector3D& particleDir,
						                               Acts::ParticleHypothesis particle ) const;
						       

     /** struct of Particle Masses */
     static ParticleMasses                   s_particleMasses;
     static PdgToParticleHypothesis          s_pdgToHypo;
           
   };
   
   /** Return the configuration object */    
   inline HadronicInteractionParametricSampler::Config HadronicInteractionParametricSampler::getConfiguration() const { return m_config; }
}

#endif // ACTS_FATRASTOOLS_HADRONICINTERACTIONPARAMETRICSAMPLER_H