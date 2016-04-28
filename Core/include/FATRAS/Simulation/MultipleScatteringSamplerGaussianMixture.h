///////////////////////////////////////////////////////////////////
// MultipleScatteringSamplerGaussianMixture.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_FATRASTOOLS_MULTIPLESCATTERINGSAMPLERGAUSSIANMIXTURE_H
#define ACTS_FATRASTOOLS_MULTIPLESCATTERINGSAMPLERGAUSSIANMIXTURE_H 1

// FATRAS 
#include "FATRAS/IRandomNumbers.h"
#include "FATRAS/IMultipleScatteringSampler.h"
// ACTS
#include "ACTS/Utilities/Definitions.h"
#include "ACTS/EventData/ParticleHypothesis.h"
#include "ACTS/Material/MaterialProperties.h"

namespace Fatras {
 
  /** @class MultipleScatteringSamplerGaussianMixture
   * 
   * @author Andreas Salzburger <Andreas.Salzburger@cern.ch>
   * @author Noemi Calace       <Noemi.Calace@cern.ch>
   * @author Artem Basalaev     <Artem.Baralaev@cern.ch>     
   * 
   */
 
  class MultipleScatteringSamplerGaussianMixture : virtual public IMultipleScatteringSampler {
     
  public:
      /** Config 
      Configuration object for this MultipleScatteringSampler*/
      struct Config {
          std::shared_ptr<Acts::IRandomNumbers>  randomNumbers;   /** Random Generator service  */
          bool                                     log_include;   /** boolean switch to include log term  */
          bool                            optGaussianMixtureG4;   /** modifies the Fruehwirth/Regler model to fit with G4 */

      };
      
      /** AlgTool like constructor */
      MultipleScatteringSamplerGaussianMixture(const Config& msConfig);
     
      /**Virtual destructor*/
      virtual ~MultipleScatteringSamplerGaussianMixture();
     
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
      Config            m_config; //!< configuration object     
     
  private:
      
      /** struct of Particle Masses */
      static Acts::ParticleMasses s_particleMasses;
     
      /** main factor of Rutherford-Scott formula */
      static double               s_main_RutherfordScott;  
      /** log factor of Rutherford-Scott formula */
      static double               s_log_RutherfordScott;   
                                  
      /** main factor for Rossi-Greisen formula */
      static double               s_main_RossiGreisen;     
      /** main factor for Rossi-Greisen formula */
      static double               s_log_RossiGreisen;      
      
      /** ========= Gaussian mixture model Fruehwirth, Regler Nucl. Inst. Methods A 456(2001) ========= */
      /** Gaussian mixture model: Sigma parameter a0 */
      static double         s_gausMixSigma1_a0;     
      /** Gaussian mixture model: Sigma parameter a1 */
      static double         s_gausMixSigma1_a1;     
      /** Gaussian mixture model: Sigma parameter a2 */
      static double         s_gausMixSigma1_a2;     
      
      /** Gaussian mixture model: Epsilon parameter a0 */
      static double         s_gausMixEpsilon_a0;     
      /** Gaussian mixture model: Epsilon parameter a1 */
      static double         s_gausMixEpsilon_a1;     
      /** Gaussian mixture model: Epsilon parameter a2 */
      static double         s_gausMixEpsilon_a2;     
     
      /** Gaussian mixture model: Epsilon parameter b0 */
      static double         s_gausMixEpsilon_b0;     
      /** Gaussian mixture model: Epsilon parameter b1 */
      static double         s_gausMixEpsilon_b1;     
      /** Gaussian mixture model: Epsilon parameter b2 */
      static double         s_gausMixEpsilon_b2;    
       
      /** projection factor to scale the projected angle out of the plane */
      static double         s_projectionFactor;      

     
    };
 
 
} // end of namespace


#endif // ACTS_FATRASTOOLS_MULTIPLESCATTERINGSAMPLERGAUSSIANMIXTURE_H