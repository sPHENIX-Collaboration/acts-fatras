///////////////////////////////////////////////////////////////////
// MaterialInteractionEngine.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_FATRAS_FATRASMATERIALEFFECTSENGINE_H
#define ACTS_FATRAS_FATRASMATERIALEFFECTSENGINE_H 1

#include <ACTS/EventData/NeutralParameters.hpp>
#include <ACTS/EventData/ParticleDefinitions.hpp>
#include <ACTS/EventData/TrackParameters.hpp>
#include <ACTS/Extrapolation/ExtrapolationCell.hpp>
#include <ACTS/Extrapolation/IMaterialEffectsEngine.hpp>
#include <ACTS/Extrapolation/MaterialUpdateMode.hpp>
#include <ACTS/Extrapolation/detail/ExtrapolationMacros.hpp>
#include <ACTS/Utilities/Definitions.hpp>
#include <ACTS/Utilities/Logger.hpp>

#include "Fatras/IEnergyLossSampler.hpp"
#include "Fatras/IHadronicInteractionSampler.hpp"
#include "Fatras/IMultipleScatteringSampler.hpp"
#include "Fatras/IPhotonConversionSampler.hpp"

namespace Fatras {

/// @class MaterialInteractionEngine
///
/// Material effects engine interface for charged and neutral (fast track
/// simulation) ,
/// the update is alwyas on the:
/// - eCell.leadParmaeters && eCell.leadLayer
/// - if eCell.leadPameters == eCell.startParamters clone to new parameters
/// else : update the new parameters
///
/// @author Andreas Salzburger <Andreas.Salzburger@cern.ch>
/// @author Noemi Calace       <Noemi.Calace@cern.ch>
/// @author Sharka Todorova    <Sarka.Todorova@cern.ch>

template <class RandomGenerator>
class MaterialInteractionEngine : virtual public Acts::IMaterialEffectsEngine {
 public:
  /// @struct Config
  /// Configuration struct for the MaterialInteractionEngine
  ///
  /// If a sampler is not defined, the given process it not performed
  struct Config {
    /// Incoming particle: min momentum (if cut > 0.)
    double particleMinMomentum;
    /// Incoming particle: min transverse momentum (if cut > 0.)
    double particleMinMomentumT;
    /// Incoming particle: kill it if it falls below, will throw
    /// SuccessUpdateKill
    bool particleKillBelowCut;
    /// Sampler: IEnergyLossSampler for ionisation loss, no energy loss if not
    /// defined
    std::shared_ptr<IEnergyLossSampler<RandomGenerator>> energyLossSampler;
    /// Sampler: dedicated electron energy samples, use standard one if not
    /// defined
    std::shared_ptr<IEnergyLossSampler<RandomGenerator>>
        energyLossSamplerElectrons;
    /// Boolean switch to record the brem photons
    bool recordBremPhoton;
    /// Outgoing: minimal momentum for a photon to be recorded (if cut > 0.)
    double outBremPhotonMinMomentum;
    /// Outgoing: minimal transverse momentum for a photon (if cut > 0.)
    double outBremPhotonMinMomentumT;
    /// Multiple scattering sampler
    std::shared_ptr<IMultipleScatteringSampler<RandomGenerator>>
        multipleScatteringSampler;
    /// Conversion sampler
    std::shared_ptr<IPhotonConversionSampler> conversionSampler;
    /// Outgoing: minimum momentum for
    double outConversionProductMinMomentum;
    /// conversion product (if cut > 0.)
    /// Outgoing: minimum transvere
    double outConversionProductMinMomentumT;
    /// momentum for conversion product
    /// (if cut > 0.)

    /// Hadronic interaction sampler
    std::shared_ptr<IHadronicInteractionSampler> hadronicInteractionSampler;
    /// Outgoing: minimum momentum for HIproducts (if cut > 0.)
    double outHadIntProductMinMomentum;
    /// Outgoint: minim transverse momentum for for HI products (ifcut > 0.)
    double outHadIntProductMinMomentumT;
    /// Output prefix
    std::string prefix;
    /// Output postfix
    std::string postfix;
    /// @todo
    bool parametricScattering = false;

    Config()
        : particleMinMomentum(-1.),
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
          postfix(" - ") {}
  };

  /// Constructor
  /// @param[in] The configuration object
  /// @param[in] Optional logger object
  MaterialInteractionEngine(const Config& miConfig,
                            std::unique_ptr<const Acts::Logger> logger =
                                Acts::getDefaultLogger("MaterialEffectsEngine",
                                                       Acts::Logging::INFO));

  /// Destructor
  ~MaterialInteractionEngine() = default;

  /// Handle material for charged extrapolation
  /// @param[in,out] ecCharged The extrapolation cell for charged parameters.
  /// The exrapolation cell and its parameters will be updated.
  /// @param[in] msurface The material surface
  /// @param[in] dir The propagation direction (along or oppoaite momentum)
  /// @param[in] matupstage Indicating if the material update should be done
  /// fully or split before and after the surface.
  Acts::ExtrapolationCode handleMaterial(
      Acts::ExCellCharged& ecCharged, const Acts::Surface* msurface = nullptr,
      Acts::PropDirection dir = Acts::alongMomentum,
      Acts::MaterialUpdateStage matupstage = Acts::fullUpdate) const final;

  /// Handle material for neutral extrapolation - only for Fatras, dummy
  /// implementation here
  /// @param[in,out] ecNeutral The extrapolation cell for neutral parameters.
  /// The exrapolation cell and its parameters will be updated.
  /// @param[in] msurface The material surface
  /// @param[in] dir The propagation direction (along or oppoaite momentum)
  /// @param[in] matupstage Indicating if the material update should be done
  /// fully or split before and after the surface.
  Acts::ExtrapolationCode handleMaterial(
      Acts::ExCellNeutral& ecNeutral, const Acts::Surface* msurface = nullptr,
      Acts::PropDirection dir = Acts::alongMomentum,
      Acts::MaterialUpdateStage matupstage = Acts::fullUpdate) const final;

  /// Set configuration method
  /// @param[in] config The configuration object
  void setConfiguration(const Config& config);

  /// Hand back the configuartion
  /// @return The configuartion object
  Config getConfiguration() const { return m_config; }

  /// Set logging instance
  /// @param[in] logger The logging instance to be set
  void setLogger(std::unique_ptr<const Acts::Logger> logger);

  /// (Re)set the random number generator
  /// @note This generator should be an algorithm local generator in order to
  /// assure thread safety. Therefore in every execution of the underlying
  /// algorithm a new generator should be generated and set using
  /// this function
  void setRandomGenerator(RandomGenerator& randomGenerator);

 private:
  /// main templated handleMaterialT method - to be called by the concrete type
  /// ones
  /// @param[in,out] eCll The extrapolation cell for track parameters.
  /// The exrapolation cell and its parameters will be updated.
  /// @param[in] msurface The material surface
  /// @param[in] dir The propagation direction (along or oppoaite momentum)
  /// @param[in] matupstage Indicating if the material update should be done
  /// fully or split before and after the surface.
  template <class T>
  Acts::ExtrapolationCode handleMaterialT(
      Acts::ExtrapolationCell<T>& eCell,
      const Acts::Surface* msurface = nullptr,
      Acts::PropDirection dir = Acts::alongMomentum,
      Acts::MaterialUpdateStage matupstage = Acts::fullUpdate) const;

  /// Process the material on the surface
  /// @param[in,out] eCll The extrapolation cell for track parameters.
  /// The exrapolation cell and its parameters will be updated.
  /// @param[in] msurface The material surface
  /// @param[in] dir The propagation direction (along or oppoaite momentum)
  /// @param[in] mprop The material properties
  /// @param[in] pathCorrection The thickness correction to the incident angle
  /// @param[in] mFraction @todo
  template <class T>
  Acts::ExtrapolationCode processOnSurfaceT(
      Acts::ExtrapolationCell<T>& eCell, const Acts::Surface* msurface,
      Acts::PropDirection dir, const Acts::MaterialProperties& mprop,
      double pathCorrection, float& mFraction) const;

  /// Update the TrackParameters accordingly - charged parameters
  /// @param[in] The input track parameters
  /// @param[in,out] eCll The extrapolation cell for track parameters.
  /// The exrapolation cell and its parameters will be updated.
  /// @param[in] msurface The material surface
  /// @param[in] dir The propagation direction (along or oppoaite momentum)
  /// @param[in] mprop The material properties
  /// @param[in] dX0 The thickness in X0 @todo why mat prop and this handed
  /// over?
  /// @param[in] pathCorrection The thickness correction to the incident angle
  /// @param[in] mFraction @todo
  /// @return The output track parameters after electromagnetic interaction
  std::unique_ptr<const Acts::TrackParameters> electroMagneticInteraction(
      const Acts::TrackParameters& parameters, Acts::ExCellCharged& eCell,
      const Acts::Surface* msurface, Acts::PropDirection dir,
      const Acts::MaterialProperties& mprop, double dX0, double pathCorrection,
      double mFraction) const;

  /// update the TrackParameters accordingly - neutral parameters
  /// dummy implementation
  /// don't do anything, no EM physics for neutral particles! @todo remove?
  std::unique_ptr<const Acts::NeutralParameters> electroMagneticInteraction(
      const Acts::NeutralParameters& parameters, Acts::ExCellNeutral& eCell,
      const Acts::Surface* msurface, Acts::PropDirection dir,
      const Acts::MaterialProperties& mprop, double dX0, double pathCorrection,
      double mFraction) const;

  //      /** create the interaction for charged */
  //      std::vector<Acts::InteractionVertex> interact(Acts::ExCellCharged&
  //      eCell, const Acts::Material&) const;
  //
  //      /** create the interaction for neutral */
  //      std::vector<Acts::InteractionVertex> interact(Acts::ExCellNeutral&
  //      eCell, const Acts::Material&) const;
  //

  /// Multiple coulomb scattering
  /// @param[in] pars The input charged track parameters
  /// @param[out] parameters The output track parameters
  /// @todo make it consistent either as return or function parameter for all
  /// functions
  /// @param[in] simTheta The simulated deviation in theta due to multiple
  /// scattering
  /// @param[in] num_deltaPhi @todo rename
  void multipleScatteringUpdate(const Acts::TrackParameters& pars,
                                Acts::ActsVectorD<5>& parameters,
                                double simTheta, double num_deltaPhi) const;

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

  /** Check for momentum cuts - checks on momentum and transverse momentum cuts
   * on Acts::Vector3D& */
  /// @todo is that used and needed?
  bool passesMomentumCuts(const Acts::Vector3D& momentum, double cutM,
                          double cutT) const;

  const Acts::Logger& logger() const { return *m_logger; }

  ///  charged extrapolation
  ///  depending on the MaterialUpdateStage:
  ///    - postUpdate : creates a new unique_ptr and stores them as step
  /// parameters
  ///    - preUpdate | fullUpdate : manipulates the parameters and returns a
  /// nullptr
  ///   nothing to do (e.g. no material) : return nullptr */
  /*void updateTrackParameters(const Acts::TrackParameters& parameters,
                             Acts::ExCellCharged& eCell,
                             Acts::PropDirection dir,
                             Acts::MaterialUpdateStage matupstage) const;
*/

  /// configuartion object
  Config m_config;
  /// logger instance
  std::unique_ptr<const Acts::Logger> m_logger;
  /// The random number generator
  /// @note This generator should be an algorithm local generator in order to
  /// assure thread safety. Therefore in every execution of the underlying
  /// algorithm a new generator should be generated and set using
  /// setRandomGenerator().
  RandomGenerator* m_randomGenerator;
  /// useful for the multiple scattering calculation
  double m_projectionFactor;
  /// Struct of Particle masses
  Acts::ParticleMasses m_particleMasses;
};
}  // end of namespace

//!< define the templated function
#include "detail/MaterialInteractionEngine.ipp"

#endif  // ACTS_FATRAS_FATRASMATERIALEFFECTSENGINE_H
