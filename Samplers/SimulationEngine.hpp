//////////////////////////////////////////////////////////////////
// SimulationEngine.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_FATRAS_SIMULATIONENGINE_H
#define ACTS_FATRAS_SIMULATIONENGINE_H 1

#include <memory>
#include "Acts/EventData/ParticleDefinitions.hpp"
#include "Acts/Extrapolation/IExtrapolationEngine.hpp"

namespace Acts {
class IExtrapolationEngine;
}

namespace Fatras {

class IProcessSampler;

class ParticleBranch
{
};

/** @class SimulationEngine

 This is the Fatras simulation engine,
 it takes a particle and returns a particle tree

    @author Andreas.Salzburger@cern.ch
 */

class SimulationEngine
{
public:
  /** @struct Config

   Configuraiton object for the Simulation Engine,
   it takes the configured Extrapolator */

  struct Config
  {
    std::shared_ptr<IExtrapolationEngine>
        transportEngine;  //!< the transport engine

    Config() : transportEngine(nullptr) {}
  };

  /** Constructor */
  SimulationEngine(const Config& config);

  /** Destructor */
  ~SimulationEngine() {}

  /** public method: simulate and return what you have*/
  std::pair < Acts::ParticleProperties,
      std::vector<Acts::InteractionVertex>
      simulate(const Acts::ParticleProperties& pProperties) const;
};
}

#endif
