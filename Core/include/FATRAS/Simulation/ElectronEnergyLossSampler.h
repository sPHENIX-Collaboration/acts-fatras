///////////////////////////////////////////////////////////////////
// ElectronEnergyLossSampler.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_FATRSTOOLS_ELECTRONENERGYLOSSSAMPLER_H
#define ACTS_FATRSTOOLS_ELECTRONENERGYLOSSSAMPLER_H 1

// FATRAS include
#include "FATRAS/Common/IRandomNumbers.h"
#include "FATRAS/Simulation/IEnergyLossSampler.h"
// ACTS include
#include "ACTS/Utilities/Definitions.h"
#include "ACTS/EventData/ParticleHypothesis.h"
#include "ACTS/Extrapolation/detail/MaterialInteraction.h"
#include "ACTS/Material/MaterialProperties.h"
// STD
#include <memory>

namespace Fatras {
  
  class EnergyLoss;

  /** @class ElectronEnergyLossSampler
   * 
   * Sampler for a eloss of a track
   * It uses the Bethe-Heitler calculation for electrons
   * it extends the IElectronEnergyLossSampler interface 
   * 
   * @author Andreas Salzburger <Andreas.Salzburger@cern.ch>
   * @author Noemi Calace       <Noemi.Calace@cern.ch>
   */

  class ElectronEnergyLossSampler : virtual public IEnergyLossSampler {
   
  public:
    /** @struct Config 
      Configuration config for the EnergyLossSampler
      */
    struct Config {
        /** Random Generator service  */
        std::shared_ptr<IRandomNumbers>      randomNumbers;  
        /** the one free parameter to scale */
        double                               scaleFactor;
        
        Config() :
          randomNumbers(nullptr),
          scaleFactor(1.)
        {}
        
    };      
 
    /** Constructor with AlgTool parameters */
    ElectronEnergyLossSampler(const Config& elConfig);
   
    /** Destructor */
    ~ElectronEnergyLossSampler();
    
    /** IEnergyLossSampler public method to compute dEdX */
    double dEdX( const Acts::MaterialProperties& materialProperties,
		         double momentum,
		         Acts::ParticleHypothesis particleHypothesis = Acts::pion ) const final;
   
    /** IEnergyLossSampler public method to compute the mean and variance of the energy loss */
    Fatras::EnergyLoss energyLoss( const Acts::MaterialProperties& mat,
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
      
    /** Private method to compute the Bethe-Heitler PDF */
    std::vector<double> betheHeitlerPDF( double pathLength ) const;

    /** the formulas for energy loss evaluation */
    Acts::MaterialInteraction                  m_interactionFormulae;     
    
    /** struct of Particle masses  */
    static Acts::ParticleMasses                s_particleMasses; 
    
    /** KOverA factor in Bethe-Bloch equation [MeV*cm2/gram] */
    static double                              s_ka_BetheBloch;          
};

inline double ElectronEnergyLossSampler::dEdX(const Acts::MaterialProperties&, double, Acts::ParticleHypothesis) const
{ return 0; }

/** Return the configuration object */    
inline ElectronEnergyLossSampler::Config ElectronEnergyLossSampler::getConfiguration() const { return m_config; }

} 
#endif //  ACTS_FATRSTOOLS_ELECTRONENERGYLOSSSAMPLER_H