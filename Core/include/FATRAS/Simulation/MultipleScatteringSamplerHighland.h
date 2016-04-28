///////////////////////////////////////////////////////////////////
// MultipleScatteringSamplerHighland.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_FATRASTOOLS_MULTIPLESCATTERINGSAMPLERHIGHLAND_H
#define ACTS_FATRASTOOLS_MULTIPLESCATTERINGSAMPLERHIGHLAND_H 1
 
// ACTS include
#include "ACTS/Utilities/Definitions.h"
#include "ACTS/Material/MaterialProperties.h"
#include "ACTS/EventData/ParticleHypothesis.h"
#include "ACTS/Extrapolation/detail/MaterialInteraction.h"
// FATRAS include
#include "FATRAS/Common/IRandomNumbers.h"
#include "FATRAS/Simulation/IMultipleScatteringSampler.h"
// STD
#include <memory>

namespace Fatras {
  
  /** @class MultipleScatteringSamplerHighland
   * 
   * The Formula used is the highland formula for the projected scattering angle :
   * 
   * @f$ \theta_{ms} = \frac{13.6MeV}{p}\cdot\sqrt{t/X_{0}}[1 + 0.038\ln(t/X_{0})] @f$
   * 
   * What is returned is the square of the expectation value of the deflection
   * @f$ < (\theta_ms)^2 > = \sigma_ms^2 @f$
   * 
   * @author Andreas Salzburger <Andreas.Salzburger@cern.ch>
   * @author Noemi Calace       <Noemi.Calace@cern.ch>
   * 
   */
 
  class MultipleScatteringSamplerHighland : virtual public IMultipleScatteringSampler {
     
    public:
      /** @struct Config 
          Configuration Struct for the MultipleScattering Sampler */
      struct Config {
          
        std::shared_ptr<IRandomNumbers>    randomNumbers;   //!< the Random number service
        bool                                 log_include; //!< include the log term

        Config() :
          randomNumbers(nullptr),
          log_include(true)
        {}
      };
  
      /** AlgTool like constructor */
      MultipleScatteringSamplerHighland(const Config& msConfig);
     
      /**Virtual destructor*/
      virtual ~MultipleScatteringSamplerHighland();
     
      /** Calculate the theta introduced by multiple scattering,
       *          according to the RutherFord-Scott Formula           
       */
      double simTheta(const Acts::MaterialProperties& mat,
                      double p,
                      double pathcorrection,
                      Acts::ParticleHypothesis particle=Acts::pion) const final;
     
      /** Set configuration method */
      void setConfiguration(const Config& eeConfig);
  
      /** Get configuration method */
      Config getConfiguration() const;                                 
  
    protected:
      Config                            m_config; // the configuraiton object

    private:
      /** the formulas for multiple scattering evaluation */
      Acts::MaterialInteraction         m_interactionFormulae;     
      
      /** struct of Particle Masses */
      static Acts::ParticleMasses       s_particleMasses;
     
      /** main factor of Rutherford-Scott formula */
      static double                     s_main_RutherfordScott;  
      /** log factor of Rutherford-Scott formula */
      static double                     s_log_RutherfordScott;   
                                  
      /** main factor for Rossi-Greisen formula */
      static double                     s_main_RossiGreisen;     
      /** main factor for Rossi-Greisen formula */
      static double                     s_log_RossiGreisen;      
                                  
      /** projection factor to scale the projected angle out of the plane */
      static double                     s_projectionFactor;      
     
    };

    /** Return the configuration object */    
    inline MultipleScatteringSamplerHighland::Config MultipleScatteringSamplerHighland::getConfiguration() const { return m_config; }
 
} // end of namespace


#endif // ACTS_FATRASTOOLS_MULTIPLESCATTERINGSAMPLERHIGHLAND_H