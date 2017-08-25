//////////////////////////////////////////////////////////////////
// EnergyLoss.h, Acts project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_FATRAS_ENERGYLOSS_H
#define ACTS_FATRAS_ENERGYLOSS_H 1

#include <math.h>
#include <cassert>
#include <iostream>

namespace Fatras {

/// @class EnergyLoss
///
/// @brief The energy loss description
///
/// Energy loss through ionisation and/or radiation leads to a change
/// (reduction) of the momentum. It uncertainty can be asymmetric in this
/// class. The quantity is energy since the calculation from energy to
/// momentum can be done better inside the FatrasMaterialEffectsEngine
/// (which knows the particle hypothesis).
///
/// @author Common tracking software group
/// @author Noemi Calace <Noemi.Calace@cern.ch>

class EnergyLoss {
 public:
  /// (Default) Constructor
  /// @param[in] deltaE The energy loss
  /// @param[in] sigmaDeltaE The uncertainty of the energy loss
  /// @param[in] sigmaMinusDeltaE The negative side of sigmaDeltaE
  /// @param[in] sigmaPlusDeltaE The positive side of sigmaDeltaE
  EnergyLoss(double deltaE = 0.0, double sigmaDeltaE = 0.0,
             double sigmaMinusDeltaE = 0.0, double sigmaPlusDeltaE = 0.0);

  /// Constructor
  /// @param[in] deltaE The energy loss
  /// @param[in] sigmaDeltaE The uncertainty of the energy loss
  /// @param[in] meanIonization The mean ionization
  /// @param[in] sigmaIonization The uncertainty of the mean ionization
  /// @param[in] meanRadiation The mean radiation
  /// @param[in] sigmaRadiation The uncertainty of the eman radiation
  EnergyLoss(double deltaE, double sigmaDeltaE, double meanIonization,
             double sigmaIonization, double meanRadiation,
             double sigmaRadiation);

  /// Constructor
  /// @param[in] deltaE The energy loss
  /// @param[in] sigmaDeltaE The uncertainty of the energy loss
  /// @param[in] sigmaMinusDeltaE The negative side of sigmaDeltaE
  /// @param[in] sigmaPlusDeltaE The positive side of sigmaDeltaE
  /// @param[in] meanIonization The mean ionization
  /// @param[in] sigmaIonization The uncertainty of the mean ionization
  /// @param[in] meanRadiation The mean radiation
  /// @param[in] sigmaRadiation The uncertainty of the eman radiation
  /// @param[in] length The length along which the energy loss happened
  EnergyLoss(double deltaE, double sigmaDeltaE, double sigmaMinusDeltaE,
             double sigmaPlusDeltaE, double meanIonization,
             double sigmaIonization, double meanRadiation,
             double sigmaRadiation, double length);

  /// Destructor
  virtual ~EnergyLoss() = default;

  /// Implicit constructor
  virtual EnergyLoss* clone() const;

  /// @return The energy loss @f$ \Delta E @f$
  double deltaE() const;

  /// @return The symmatric error @f$ \sigma(\Delta E) @f$
  double sigmaDeltaE() const;

  /// @return The negative side @f$ \sigma(\Delta E) @f$
  double sigmaMinusDeltaE() const;

  /// @return The positive side @f$ \sigma(\Delta E) @f$
  double sigmaPlusDeltaE() const;

  /// @return The mean ionization
  double meanIonization() const;
  /// @return The uncertainty of the mean ionization
  double sigmaIonization() const;
  /// @return The mean radiation
  double meanRadiation() const;
  /// @return The uncertainty of the eman radiation
  double sigmaRadiation() const;
  /// @return The length along which the energy loss happened
  double length() const;

  /// Update from mean values
  /// @brief Adds the new parameters to the old parameters
  /// @param[in] meanIonization The mean ionization
  /// @param[in] sigmaIonization The uncertainty of the mean ionization
  /// @param[in] meanRadiation The mean radiation
  /// @param[in] sigmaRadiation The uncertainty of the mean radiation
  /// @param[in] mpv If set to true the most probable value (which differs from
  /// the
  /// mean for the landau distributed energy loss) instead of the mean value
  /// will be taken to calculate the energy loss
  void update(double meanIonization, double sigmaIonization,
              double meanRadiation, double sigmaRadiation, bool mpv = false);

  /// Update energy loss
  /// @brief Adds the new parameters to the old parameters
  /// @param[in] eLoss The energy loss
  /// @param[in] mpv If set to true the most probable value (which differs from
  /// the
  /// mean for the landau distributed energy loss) instead of the mean value
  /// will be taken to calculate the energy loss
  void update(EnergyLoss& eLoss, bool mpv = false);

 private:
  /// @f$ \Delta E @f$ - the estimated or measured energy loss
  double m_deltaE;
  /// < @f$ \sigma(\Delta E) @f$ - error on the energy loss
  double m_sigmaDeltaE;
  /// < @f$ \sigma(\Delta E) @f$ - negative error on the energy loss
  double m_sigmaMinusDeltaE;
  /// < @f$ \sigma(\Delta E) @f$ - positive error on the energy loss
  double m_sigmaPlusDeltaE;
  /// Mean value for ionization
  double m_meanIonization;
  /// Sigma for ionization
  double m_sigmaIonization;
  /// Mean value for radiation
  double m_meanRadiation;
  /// Sigma for radiation
  double m_sigmaRadiation;
  /// 3D length of material
  double m_length;
};

/// Overload of << operator for both, MsgStream and std::ostream for debug
/// output
std::ostream& operator<<(std::ostream& sl, const EnergyLoss& eloss);

inline EnergyLoss* EnergyLoss::clone() const { return new EnergyLoss(*this); }

inline double EnergyLoss::deltaE() const { return m_deltaE; }

inline double EnergyLoss::sigmaDeltaE() const { return m_sigmaDeltaE; }

inline double EnergyLoss::sigmaMinusDeltaE() const {
  return m_sigmaMinusDeltaE;
}

inline double EnergyLoss::sigmaPlusDeltaE() const { return m_sigmaPlusDeltaE; }

inline double EnergyLoss::meanIonization() const { return m_meanIonization; }

inline double EnergyLoss::sigmaIonization() const { return m_sigmaIonization; }

inline double EnergyLoss::meanRadiation() const { return m_meanRadiation; }

inline double EnergyLoss::sigmaRadiation() const { return m_sigmaRadiation; }

inline double EnergyLoss::length() const {
  return m_length;
}  // length can be positive and negative like Eloss depending on (back)tracking

inline void EnergyLoss::update(double meanIonization, double sigmaIonization,
                               double meanRadiation, double sigmaRadiation,
                               bool mpv) {
  m_meanIonization += meanIonization;
  m_meanRadiation += meanRadiation;
  m_sigmaIonization += sigmaIonization;
  m_sigmaRadiation += sigmaRadiation;
  m_deltaE += mpv ? 0.9 * meanIonization + 0.15 * meanRadiation
                  : meanIonization + meanRadiation;
  m_sigmaDeltaE = sqrt(m_sigmaIonization * m_sigmaIonization +
                       m_sigmaRadiation * m_sigmaRadiation);
}

inline void EnergyLoss::update(EnergyLoss& eloss, bool mpv) {
  m_meanIonization += eloss.meanIonization();
  m_meanRadiation += eloss.meanRadiation();
  m_sigmaIonization += eloss.sigmaIonization();
  m_sigmaRadiation += eloss.sigmaRadiation();
  m_deltaE += mpv ? 0.9 * eloss.meanIonization() + 0.15 * eloss.meanRadiation()
                  : eloss.meanIonization() + eloss.meanRadiation();
  m_sigmaDeltaE = sqrt(m_sigmaIonization * m_sigmaIonization +
                       m_sigmaRadiation * m_sigmaRadiation);
}
}
#endif  // ACTS_FATRAS_ENERGYLOSS_H
