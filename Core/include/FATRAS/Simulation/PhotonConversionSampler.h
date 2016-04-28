///////////////////////////////////////////////////////////////////
// PhotonConversionSampler.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_FATRASTOOLS_PHOTONCONVERSIONSAMPLER_H
#define ACTS_FATRASTOOLS_PHOTONCONVERSIONSAMPLER_H 1

// ACTS includes
#include "ACTS/Utilities/Definitions.h"
#include "ACTS/EventData/ParticleHypothesis.h"

// FATRAS includes
#include "FATRAS/Common/IRandomNumbers.h"
#include "FATRAS/Simulation/IPhotonConversionSampler.h"
#include "FATRAS/Simulation/detail/PdgToParticleHypothesis.h"

// STL
#include <algorithm>
#include <memory>


namespace Acts {
    class InteractionVertex;
}

namespace Fatras {
    
  /**
   * @class PhotonConversionSampler
   * 
   * The PhotonConversionSampler called by the FatrasMaterialEffecsEngine
   * 
   * @TODO write documentation 
   * 
   * @author Sarka Todorova <Sarka.Todorova@cern.ch>
   * @author Noemi Calace   <Noemi.Calace@cern.ch>
   * @author Andreas Salzburger <Andreas.Salzburger@cern.ch>
   * 
  */
  class PhotonConversionSampler : virtual public IPhotonConversionSampler {

    public:
      /** @struct Config 
        Configuration struct for this PhotonConversionSampler */
      struct Config {
          /** Random Generator service  */
          std::shared_ptr<IRandomNumbers>    randomNumbers;  
          int                                processCode;
          double                             minChildEnergy;
          double                             childEnergyScaleFactor;
          
          Config() :
            randomNumbers(nullptr), 
            processCode(1),
            minChildEnergy(50.),
            childEnergyScaleFactor(1.)
          {}
          
      };  
    
      /** Constructor */
      PhotonConversionSampler(const Config& pcConfig);
          
      /** Destructor */    
      ~PhotonConversionSampler();
      
      /** processing of the presampled conversion on layer 
       * This methods return the photon's children
       */
      std::vector<Acts::InteractionVertex> doConversion(double time,
	  					                                const Acts::Vector3D& position , 
	  					                                const Acts::Vector3D& momentum) const final;
    
      /** Set configuration method */
      void setConfiguration(const Config& pcConfig);
    
      /** Get configuration method */
      Config getConfiguration() const;                                 
    
    protected:
      
        /** Method to calculate the fracton of energy for the child */
      double childEnergyFraction(double gammaMom) const;
      
      Acts::Vector3D childDirection(const Acts::Vector3D& gammaMom, double childE) const;
      
      std::vector<Acts::InteractionVertex> getChildren( double time, 
                                                        const Acts::Vector3D& vertex,
                                                        const Acts::Vector3D& photonMomentum,
                                                        double childEnergy, 
                                                        const Acts::Vector3D& childDirection,
                                                        Acts::ParticleHypothesis childType) const;
      /** helper functions for the Phi1/phi2 */
      double phi1(double delta) const;
      
      /** helper functions for the Phi1/phi2 */
      double phi2(double delta) const;
      
      /** the congifuration */
      Config                                     m_config;
      
      /** statics for interaction */
      static Acts::ParticleMasses                s_particleMasses;
      static PdgToParticleHypothesis             s_pdgToHypo;
      static double                              s_alpha;
      static double                              s_oneOverThree;
      
  };
  
  /** Return the configuration object */    
  inline PhotonConversionSampler::Config PhotonConversionSampler::getConfiguration() const { return m_config; }
  
  /** phi1 value of the first electron - @TODO write documentation */
  inline double PhotonConversionSampler::phi1(double delta) const {
    if (delta <= 1.)
      return 20.867 - 3.242 * delta  + 0.625*delta*delta;
    else
      return 21.12 - 4.184*log(delta+0.952);
  }
  /** phi1 value of the second poton - @TODO write documentation */
  inline double PhotonConversionSampler::phi2(double delta) const {
    if (delta <= 1.)
      return 20.209 - 1.930 * delta  - 0.086*delta*delta;
    return 21.12 - 4.184*log(delta+0.952);
  }
  
  
}

#endif // ACTS_FATRASTOOLS_PHOTONCONVERSIONSAMPLER_H