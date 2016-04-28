///////////////////////////////////////////////////////////////////
// EnergyLoss.cxx, ACTS project
///////////////////////////////////////////////////////////////////

// FATRAS includes
#include "FATRAS/Simulation/EnergyLoss.h"
// output
#include <iostream>
#include <iomanip>
#include <string>
#include <limits>

Fatras::EnergyLoss::EnergyLoss(double deltaE, double sigmaDeltaE,
                               double sMinusDeltaE, double sPlusDeltaE)  :
  m_deltaE(deltaE),
  m_sigmaDeltaE(sigmaDeltaE),
  m_sigmaMinusDeltaE(sMinusDeltaE>0.0?sMinusDeltaE:sigmaDeltaE),
  m_sigmaPlusDeltaE(sPlusDeltaE>0.0?sPlusDeltaE:sigmaDeltaE),
  m_mean_ioni(0.),
  m_sig_ioni(0.),
  m_mean_rad(0.),
  m_sig_rad(0.),
  m_length(-std::numeric_limits<double>::min())
{ }

Fatras::EnergyLoss::EnergyLoss(double deltaE, double sigmaDeltaE,
                               double sMinusDeltaE, double sPlusDeltaE, 
                               double mean_ioni, double sig_ioni,
                               double mean_rad, double sig_rad, double length)  :
  m_deltaE(deltaE),
  m_sigmaDeltaE(sigmaDeltaE),
  m_sigmaMinusDeltaE(sMinusDeltaE),
  m_sigmaPlusDeltaE(sPlusDeltaE),
  m_mean_ioni(mean_ioni),
  m_sig_ioni(sig_ioni),
  m_mean_rad(mean_rad),
  m_sig_rad(sig_rad),
  m_length(length)
{ }

Fatras::EnergyLoss::EnergyLoss(double deltaE, double sigmaDeltaE,
                               double mean_ioni, double sig_ioni,
                               double mean_rad, double sig_rad)  :
  m_deltaE(deltaE),
  m_sigmaDeltaE(sigmaDeltaE),
  m_sigmaMinusDeltaE(0.),
  m_sigmaPlusDeltaE(0.),
  m_mean_ioni(mean_ioni),
  m_sig_ioni(sig_ioni),
  m_mean_rad(mean_rad),
  m_sig_rad(sig_rad),
  m_length(-std::numeric_limits<double>::min())
{ }

std::ostream& Fatras::operator << ( std::ostream& sl, const Fatras::EnergyLoss& eloss)
{ 
  sl << "EnergyLoss :   ( delta(E), sigma(dE) ) = \t"
     << "("<< eloss.deltaE()<<", \t"<< eloss.sigmaDeltaE() << ")";
  return sl; 
}
