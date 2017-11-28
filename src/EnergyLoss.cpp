///////////////////////////////////////////////////////////////////
// EnergyLoss.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

#include "Fatras/EnergyLoss.hpp"

Fatras::EnergyLoss::EnergyLoss(double deltaE, double sigmaDeltaE,
                               double sMinusDeltaE, double sPlusDeltaE)
    : m_deltaE(deltaE),
      m_sigmaDeltaE(sigmaDeltaE),
      m_sigmaMinusDeltaE(sMinusDeltaE > 0.0 ? sMinusDeltaE : sigmaDeltaE),
      m_sigmaPlusDeltaE(sPlusDeltaE > 0.0 ? sPlusDeltaE : sigmaDeltaE),
      m_meanIonization(0.),
      m_sigmaIonization(0.),
      m_meanRadiation(0.),
      m_sigmaRadiation(0.),
      m_length(-std::numeric_limits<double>::min()) {}

Fatras::EnergyLoss::EnergyLoss(double deltaE, double sigmaDeltaE,
                               double sMinusDeltaE, double sPlusDeltaE,
                               double meanIonization, double sigmaIonization,
                               double meanRadiation, double sigmaRadiation,
                               double length)
    : m_deltaE(deltaE),
      m_sigmaDeltaE(sigmaDeltaE),
      m_sigmaMinusDeltaE(sMinusDeltaE),
      m_sigmaPlusDeltaE(sPlusDeltaE),
      m_meanIonization(meanIonization),
      m_sigmaIonization(sigmaIonization),
      m_meanRadiation(meanRadiation),
      m_sigmaRadiation(sigmaRadiation),
      m_length(length) {}

Fatras::EnergyLoss::EnergyLoss(double deltaE, double sigmaDeltaE,
                               double meanIonization, double sigmaIonization,
                               double meanRadiation, double sigmaRadiation)
    : m_deltaE(deltaE),
      m_sigmaDeltaE(sigmaDeltaE),
      m_sigmaMinusDeltaE(0.),
      m_sigmaPlusDeltaE(0.),
      m_meanIonization(meanIonization),
      m_sigmaIonization(sigmaIonization),
      m_meanRadiation(meanRadiation),
      m_sigmaRadiation(sigmaRadiation),
      m_length(-std::numeric_limits<double>::min()) {}

std::ostream& Fatras::operator<<(std::ostream& sl,
                                 const Fatras::EnergyLoss& eloss) {
  sl << "EnergyLoss :   ( delta(E), sigma(dE) ) = \t"
     << "(" << eloss.deltaE() << ", \t" << eloss.sigmaDeltaE() << ")";
  return sl;
}
