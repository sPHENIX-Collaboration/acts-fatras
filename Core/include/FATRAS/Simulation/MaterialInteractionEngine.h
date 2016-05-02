///////////////////////////////////////////////////////////////////
// MaterialInteractionEngine.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_FATRASSERVICES_FATRASMATERIALEFFECTSENGINE_H
#define ACTS_FATRASSERVICES_FATRASMATERIALEFFECTSENGINE_H 1

// ACTS include
#include "ACTS/Extrapolation/IMaterialEffectsEngine.h"
#include "ACTS/Extrapolation/ExtrapolationCell.h"
#include "ACTS/Extrapolation/MaterialInteraction.h"
#include "ACTS/Extrapolation/MaterialUpdateMode.h"
#include "ACTS/Extrapolation/detail/ExtrapolationMacros.h"
#include "ACTS/Utilities/PropDirection.h"
#include "ACTS/EventData/TrackParameters.h"
#include "ACTS/EventData/NeutralParameters.h"
// FATRAS include
#include "FATRAS/IEnergyLossSampler.h"
#include "FATRAS/IMultipleScatteringSampler.h"
#include "FATRAS/IPhotonConversionSampler.h"
#include "FATRAS/IHadronicInteractionSampler.h"
 
namespace Fatras {
  
  /** @class MaterialInteractionEngine
   * 
   * Material effects engine interface for charged and neutral (fast track simulation) ,
   * the update is alwyas on the:
   * - eCell.leadParmaeters && eCell.leadLayer
   * - if eCell.leadPameters == eCell.startParamters clone to new parameters
   * else : update the new parameters
   * 
   * @author Andreas Salzburger <Andreas.Salzburger@cern.ch>
   * @author Noemi Calace       <Noemi.Calace@cern.ch>
   * 
   */
  class MaterialInteractionEngine : virtual public Acts::IMaterialEffectsEngine {
    public:
      /** @struct Config 
          Configuration struct for the MaterialInteractionEngine 
        
          If a sampler is not defined, the given process it not performed.
        */    
      struct Config {

          std::shared_ptr<IRandomNumbers>              randomNumbers;                    //!< Random Generator service 
                                                                                         
          double                                       inParticleMinMomentum;            //!< incoming particle: min momentum
          double                                       inParticleMinMomentumT;           //!< incoming particle: min transverse momentum (additional, optional)  
          bool                                         inParticleKillBelowCut;           //!< incoming particle: boolean switch to kill it if it falls below, will throw SuccessUpdate 
                                                                                         
          std::shared_ptr<IEnergyLossSampler>          energyLossSampler;                //!< Sampler: IEnergyLossSampler for ionisation loss, no energy loss if not defined
                                                                                         
          std::shared_ptr<IEnergyLossSampler>          energyLossSamplerElectrons;       //!< Sampler: dedicated electron energy samples, use standard one if not defined
          bool                                         recordBremPhoton;                 //!< boolean switch to record the brem photons
          double                                       outBremPhotonMinMomentum;         //!< outgoing: minimal momentum for a photon to be recorded
          double                                       outBremPhotonMinMomentumT;        //!< outgoing: minimal transverse momentum for a photon (additional, optional)
                                                                                         
          std::shared_ptr<IMultipleScatteringSampler>  multipleScatteringSampler;        //!< multiple scattering sampler
          std::shared_ptr<IPhotonConversionSampler>    conversionSampler;                //!< conversion sampler
          double                                       outConversionProductMinMomentum;  //!< outgoing: minimum momentum for conversion product
          double                                       outConversionProductMinMomentumT; //!< outgoing: minimum transvere momentum for conversion product (additional, optional)
          
          std::shared_ptr<IHadronicInteractionSampler> hadronicInteractionSampler;       //!< hadronic interaction sampler 
          double                                       outHadIntProductMinMomentum;      //!< outgoing: minimum momentum for HI products
          double                                       outHadIntProductMinMomentumT;     //!< outgoint: minim transverse momentum for for HI products (additional, optional)
          
          std::string                                  prefix;                           //!< output prefix
          std::string                                  postfix;                          //!< output postfix
          
          Config() :
            randomNumbers(nullptr),                                                    
            inParticleMinMomentum(50.),          
            inParticleMinMomentumT(-1.),      
            inParticleKillBelowCut(true),                                           
            energyLossSampler(nullptr),                                               
            energyLossSamplerElectrons(nullptr),      
            recordBremPhoton(true),                
            outBremPhotonMinMomentum(50.),        
            outBremPhotonMinMomentumT(-1.),                                       
            multipleScatteringSampler(nullptr),       
            conversionSampler(nullptr),               
            outConversionProductMinMomentum(50.), 
            outConversionProductMinMomentumT(-1.),
            hadronicInteractionSampler(nullptr),
            outHadIntProductMinMomentum(50.),     
            outHadIntProductMinMomentumT(-1.) 
            prefix("[MI] - "),
            postfix(" - ")
          {}   
      };

      /** Constructor */
      MaterialInteractionEngine(const Config&miConfig);

      /** Destructor */
      ~MaterialInteractionEngine();

      /** charged extrapolation */
      Acts::ExtrapolationCode handleMaterial(Acts::ExCellCharged& ecCharged,
                                             Acts::PropDirection dir=Acts::alongMomentum,
                                             Acts::MaterialUpdateStage matupstage=Acts::fullUpdate) const final;

      /** neutral extrapolation - only for Fatras, dummy implementation here */
      Acts::ExtrapolationCode handleMaterial(Acts::ExCellNeutral& ecNeutral,
                                             Acts::PropDirection dir=Acts::alongMomentum,
                                             Acts::MaterialUpdateStage matupstage=Acts::fullUpdate) const final;

      /** Set configuration method */
      void setConfiguration(const Config& eeConfig);
      
      /** Get configuration method */
      Config getConfiguration() const;                                 
                               
    protected:

//      /** main templated handleMaterialT method - to be called by the concrete type ones */
//      template <class T> Acts::ExtrapolationCode handleMaterialT( Acts::ExtrapolationCell<T>& eCell,
//							                                      Acts::PropDirection dir=Acts::alongMomentum,
//							                                      Acts::MaterialUpdateStage matupstage=Acts::fullUpdate) const;
//		
//      /** process the material on the layer, check if in-layer interaction needs to occur */    					    
//      Acts::ExtrapolationCode processMaterialOnLayer(Acts::ExCellCharged& ecCharged,
//					                                 Acts::PropDirection dir,
//					                                 float& mFraction) const;
//
//      /** process the material on the layer, check if in-layer interaction needs to occur */    					    
//      Acts::ExtrapolationCode processMaterialOnLayer(Acts::ExCellNeutral& ecNeutral,
//					                                 Acts::PropDirection dir,
//					                                 float& mFraction) const;
//
//      /** process the material on the layer, check if in-layer interaction needs to occur - templated */    					    
//      template <class T> Acts::ExtrapolationCode processMaterialOnLayerT( Acts::ExtrapolationCell<T>& eCell,
//								                                          Acts::PropDirection dir,
//								                                          float& mFraction) const;
//			
//      /** update the TrackParameters accordingly */       					    
//      const Acts::TrackParameters* updateTrackParameters(const Acts::TrackParameters& parameters,
//							                             Acts::ExCellCharged& eCell,
//							                             Acts::PropDirection dir,
//							                             double dX0,
//							                             double pathCorrection,
//							                             double mFraction) const;
//		
//      /** update the TrackParameters accordingly */       					            					 
//      const Acts::NeutralParameters* updateTrackParameters(const Acts::NeutralParameters& parameters,
//							                               Acts::ExCellNeutral& eCell,
//							                               Acts::PropDirection dir,
//							                               double dX0,
//							                               double pathCorrection,
//							                               double mFraction) const;
//		
//      /** create the interaction for charged */       					    
//      std::vector<Acts::InteractionVertex> interact(Acts::ExCellCharged& eCell, const Acts::Material&) const;
//
//      /** create the interaction for neutral */       					    
//      std::vector<Acts::InteractionVertex> interact(Acts::ExCellNeutral& eCell, const Acts::Material&) const;
//      
//      /** multiple coulomb scattering */       					    
//      void multipleScatteringUpdate(const Acts::TrackParameters& pars,
//                                    ActsVectorD<5>& parameters,
//                                    double simTheta, 
//                                    double num_deltaPhi) const;
//
//      /** Radiate a brem photon */
//      void radiate( ActsVectorD<5> & parm ,
//		            Acts::ExCellCharged& eCell, 
//		            float pathLim, 
//		            float mFr,
//		            float refX) const;
//
//      /** Collect the children */
//      void collectBremPhoton(Acts::ExCellCharged& eCell,
//			                 double pElectron,
//			                 double gammaE,
//			                 const Vector3D& vertex,
//			                 Acts::Vector3D& particleDir) const;
      
       /** ExtrapolationEngine config object */
       Config                m_config;

       Acts::ParticleMasses  particleMasses;                //!< struct of Particle Masses
                             
      /** useful for the angle calculation of the brem photon */ 
      double                 m_oneOverThree;
           
      /** useful for the multiple scattering calculation */ 
      double                 m_projectionFactor;
    
  };
  
  /** Return the configuration object */    
  inline MaterialInteractionEngine::Config MaterialInteractionEngine::getConfiguration() const { return m_config; }
  

} // end of namespace

//!< define the templated function   
// #include "FATRAS/Simulation/detail/MaterialInteractionEngine.icc" 

#endif // ACTS_FATRASSERVICES_FATRASMATERIALEFFECTSENGINE_H
