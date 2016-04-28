///////////////////////////////////////////////////////////////////
// IEnergyLossSampler.h, Acts project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_FATRASINTERFACES_IENERGYLOSSSAMPLER_H
#define ACTS_FATRASINTERFACES_IENERGYLOSSSAMPLER_H 1

// ACTS includes
#include "ACTS/EventData/ParticleHypothesis.h"
#include "ACTS/Utilities/Definitions.h"

namespace Acts {
    class MaterialProperties;
}

namespace Fatras {
  
  class EnergyLoss;

  /** @class IEnergyLossSampler
   *
   * Interface class IEnergyLossSampler
   * 
   * @author Andreas Salzburger <Andreas.Salzburger@cern.ch>
   * @author Noemi Calace       <Noemi.Calace@cern.ch>
   */

  class IEnergyLossSampler {
    
  public:
    /**Virtual destructor*/
    virtual ~IEnergyLossSampler(){}
    
    /** deltaE calculation
     * using dEdX and integrating along pathlength,
     * assuming constant dEdX during for the path.
     * - The sign depends on the given propagation direction 
     * - Units: [MeV]
    */
    virtual EnergyLoss energyLoss(const Acts::MaterialProperties& mat,
				                  double momentum,
				                  double pathcorrection,
				                  Acts::PropDirection dir=Acts::alongMomentum,
				 	              Acts::ParticleHypothesis particle=Acts::pion) const = 0;  
					 
    /** dEdX calculation when providing MaterialProperties,
     * a momentum, and a ParicleHypothesis. 
     * - Units: [Mev/mm]
     */
    virtual double dEdX(const Acts::MaterialProperties& mat,
			            double momentum, Acts::ParticleHypothesis particle=Acts::pion) const = 0;
   
  };

} // end of namespace

#endif // ACTS_FATRASINTERFACES_IENERGYLOSSSAMPLER_H