///////////////////////////////////////////////////////////////////
// MultipleScatteringSamplerGeneralMixture.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_FATRASTOOLS_MULTIPLESCATTERINGSAMPLERGENERALMIXTURE_H
#define ACTS_FATRASTOOLS_MULTIPLESCATTERINGSAMPLERGENERALMIXTURE_H 1
 
// FATRAS
#include "FATRAS/Common/IRandomNumbers.h"
#include "FATRAS/Simulation/IMultipleScatteringSampler.h"
// ACTS
#include "ACTS/Utilities/Definitions.h"
#include "ACTS/EventData/ParticleDefinitions.h"
#include "ACTS/Material/MaterialProperties.h"
// STD
#include <memory>

namespace Fatras {
  
   /** @class MultipleScatteringSamplerGeneralMixture
   *
   * ========= General mixture model Fruehwirth, M. Liendl. Comp. Phys. Comm. 141 (2001) 230-246 ========= 
   *
   * @TODO write documentation, move away from double* interface
   *     
   * @author Artem Basalaev     <Artem.Baralaev@cern.ch>     
   * @author Andreas Salzburger <Andreas.Salzburger@cern.ch>
   * @author Noemi Calace       <Noemi.Calace@cern.ch>
   * 
   */
 
  class MultipleScatteringSamplerGeneralMixture : virtual public IMultipleScatteringSampler {
     
    public:
      /** @struct Config 
          Configuration Struct for the MultipleScattering Sampler */
      struct Config {
          
        std::shared_ptr<IRandomNumbers>    randomNumbers;   //!< the Random number service
        bool                               log_include; //!< include the log term
        
        double                             genMixtureScale; //!< numberically derived factor on mixture scale
      
        Config() :
          randomNumbers(nullptr),
          log_include(true),
          genMixtureScale(0.608236)
        {}
      };
      
      /** AlgTool like constructor */
      MultipleScatteringSamplerGeneralMixture(const Config& msConfig);
     
      /**Virtual destructor*/
      virtual ~MultipleScatteringSamplerGeneralMixture();
     
      /** Calculate the theta introduced by multiple scattering,
       *          according to the RutherFord-Scott Formula           
       */
      double simTheta(const Acts::MaterialProperties& mat,
                      double p,
                      double pathcorrection,
                      Acts::ParticleType particle=Acts::pion) const final;

      /** Set configuration method */
      void setConfiguration(const Config& msConfig);
  
      /** Get configuration method */
      Config getConfiguration() const;                                 
  
    protected:
      Config     m_config; // the configuraiton object
     
    private:
      /** General mixture model: get parameters for single gaussian simulation */
      double*    getGaussian(double beta, double p,double dOverX0, double scale) const;
      /** General mixture model: get parameters for gaussian mixture */
      double*    getGaussmix(double beta, double p,double dOverX0,double Z, double scale) const;
      /** General mixture model: get parameters for semi-gaussian mixture */
      double*    getSemigauss(double beta,double p,double dOverX0,double Z, double scale) const;
      /** General mixture model: simulate semi-gaussian mixture */
      double    simGaussmix(double* scattering_params) const;
      /** General mixture model: simulate gaussian mixture */
      double    simSemigauss(double* scattering_params) const;
                
    };
 
    /** Return the configuration object */    
    inline MultipleScatteringSamplerGeneralMixture::Config MultipleScatteringSamplerGeneralMixture::getConfiguration() const { return m_config; }
 
} // end of namespace


#endif // ACTS_FATRASTOOLS_MULTIPLESCATTERINGSAMPLERGENERALMIXTURE_H