///////////////////////////////////////////////////////////////////
// MaterialInteractionEngine.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_FATRASSERVICES_FATRASMATERIALEFFECTSENGINE_H
#define ACTS_FATRASSERVICES_FATRASMATERIALEFFECTSENGINE_H 1

// ACTS include
#include "ACTS/Extrapolation/IMaterialEffectsEngine.h"
#include "ACTS/Extrapolation/ExtrapolationCell.h"
#include "ACTS/Extrapolation/MaterialUpdateMode.h"
#include "ACTS/Extrapolation/detail/ExtrapolationMacros.h"
#include "ACTS/Utilities/Definitions.h"
#include "ACTS/EventData/TrackParameters.h"
#include "ACTS/EventData/NeutralParameters.h"
// FATRAS include
#include "FATRAS/Common/IRandomNumbers.h"
#include "FATRAS/Simulation/IEnergyLossSampler.h"
#include "FATRAS/Simulation/IMultipleScatteringSampler.h"
#include "FATRAS/Simulation/IPhotonConversionSampler.h"
#include "FATRAS/Simulation/IHadronicInteractionSampler.h"
 
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
   * @uathor Sharka Todorova    <Sarka.Todorova@cern.ch>
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
                                                                                         
          double                                       particleMinMomentum;            //!< incoming particle: min momentum (if cut > 0.)
          double                                       particleMinMomentumT;           //!< incoming particle: min transverse momentum (if cut > 0.)  
          bool                                         particleKillBelowCut;           //!< incoming particle: kill it if it falls below, will throw SuccessUpdateKill 
                                                                                         
          std::shared_ptr<IEnergyLossSampler>          energyLossSampler;                //!< Sampler: IEnergyLossSampler for ionisation loss, no energy loss if not defined
                                                                                         
          std::shared_ptr<IEnergyLossSampler>          energyLossSamplerElectrons;       //!< Sampler: dedicated electron energy samples, use standard one if not defined
          bool                                         recordBremPhoton;                 //!< boolean switch to record the brem photons
          double                                       outBremPhotonMinMomentum;         //!< outgoing: minimal momentum for a photon to be recorded (if cut > 0.)
          double                                       outBremPhotonMinMomentumT;        //!< outgoing: minimal transverse momentum for a photon (if cut > 0.)
                                                                                         
          std::shared_ptr<IMultipleScatteringSampler>  multipleScatteringSampler;        //!< multiple scattering sampler
          std::shared_ptr<IPhotonConversionSampler>    conversionSampler;                //!< conversion sampler
          double                                       outConversionProductMinMomentum;  //!< outgoing: minimum momentum for conversion product (if cut > 0.)
          double                                       outConversionProductMinMomentumT; //!< outgoing: minimum transvere momentum for conversion product (if cut > 0.)
          
          std::shared_ptr<IHadronicInteractionSampler> hadronicInteractionSampler;       //!< hadronic interaction sampler 
          double                                       outHadIntProductMinMomentum;      //!< outgoing: minimum momentum for HI products (if cut > 0.)
          double                                       outHadIntProductMinMomentumT;     //!< outgoint: minim transverse momentum for for HI products (if cut > 0.)
          
          std::string                                  prefix;                           //!< output prefix
          std::string                                  postfix;                          //!< output postfix
          
          Config() :
            randomNumbers(nullptr),                                                    
            particleMinMomentum(-1.),          
            particleMinMomentumT(50.),      
            particleKillBelowCut(true),                                           
            energyLossSampler(nullptr),                                               
            energyLossSamplerElectrons(nullptr),      
            recordBremPhoton(true),                
            outBremPhotonMinMomentum(-1.),        
            outBremPhotonMinMomentumT(50.),                                       
            multipleScatteringSampler(nullptr),       
            conversionSampler(nullptr),               
            outConversionProductMinMomentum(-1.), 
            outConversionProductMinMomentumT(50.),
            hadronicInteractionSampler(nullptr),
            outHadIntProductMinMomentum(-1.),     
            outHadIntProductMinMomentumT(50.), 
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
      void setConfiguration(const Config& miConfig);
      
      /** Get configuration method */
      Config getConfiguration() const;                                 
                               
    protected:

      /** main templated handleMaterialT method - to be called by the concrete type ones */
      template <class T> Acts::ExtrapolationCode handleMaterialT( Acts::ExtrapolationCell<T>& eCell,
							                                      Acts::PropDirection dir=Acts::alongMomentum,
							                                      Acts::MaterialUpdateStage matupstage=Acts::fullUpdate) const;
		
      /** process the material on the layer, check if in-layer interaction needs to occur - templated */    					    
      template <class T> Acts::ExtrapolationCode processOnLayerT( Acts::ExtrapolationCell<T>& eCell,
								                                  Acts::PropDirection dir,
                                                                  const Acts::MaterialProperties& mprop,
                                                                  double pathCorrection,
                                                                  float& mFraction) const;
      		
      /** update the TrackParameters accordingly - charged parameters */       					    
      const Acts::TrackParameters* electroMagneticInteraction(const Acts::TrackParameters& parameters,
							                                  Acts::ExCellCharged& eCell,
							                                  Acts::PropDirection dir,
                                                              const Acts::MaterialProperties& mprop,
							                                  double dX0,
							                                  double pathCorrection,
							                                  double mFraction) const;
		
      /** update the TrackParameters accordingly - neutral parameters */       					            					 
      const Acts::NeutralParameters* electroMagneticInteraction(const Acts::NeutralParameters& parameters,
							                                    Acts::ExCellNeutral& eCell,
							                                    Acts::PropDirection dir,
                                                                const Acts::MaterialProperties& mprop,
							                                    double dX0,
							                                    double pathCorrection,
							                                    double mFraction) const;
		
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
                                                                  
                                                                  
       /** Check for momentum cuts - checks on momentum and transverse momentum cuts on Acts::Vector3D& */
       bool passesMomentumCuts(const Acts::Vector3D& momentum, double cutM, double cutT) const; 
       
       
       /** ExtrapolationEngine config object */
       Config                m_config;
           
      /** useful for the multiple scattering calculation */ 
      double                 m_projectionFactor;
    
  };
  
  /** Return the configuration object */    
  inline MaterialInteractionEngine::Config MaterialInteractionEngine::getConfiguration() const { return m_config; }
                                                                                                                  
 
  

} // end of namespace

//!< define the templated function   
#include "FATRAS/Simulation/detail/MaterialInteractionEngine.icc" 

#endif // ACTS_FATRASSERVICES_FATRASMATERIALEFFECTSENGINE_H
