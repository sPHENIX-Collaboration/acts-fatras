// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Fatras/Kernel/Definitions.hpp"
#include "Fatras/Kernel/RandomNumberDistributions.hpp"
#include "Acts/Utilities/Units.hpp"
#include <math.h>

namespace Fatras {

/// The struct for the physics list
///
///
struct ParametricNuclearInt {

struct Config()
{
	bool m_hadronInteractionFromX0;
	double m_hadronInteractionProbabilityScale;
	unsigned int MAXHADINTCHILDREN;
	double m_minimumHadOutEnergy;
}

/// @brief Constructor with given configuration
/// @param [in] cfg Configuration file
ParametricNuclearInt(Config& cfg);

/// @brief Calculates the hadronic with a given probability given by the properties of the material and the incoming particle
///
/// @tparam generator_t data type of the random number generator
/// @tparam material_t data type of the material
/// @tparam particle_t data type of the particle
/// @param generator random number generator
/// @param material penetrated material
/// @param pathCorrection scaling factor for the passed material
/// @param particle penetrating particle
template <typename generator_t, typename material_t, typename particle_t>
std::vector<particle_t> 
hadronicInteraction(generator_t& generator, const material_t& material, double pathCorrection, particle_t& particle) const;

  /// Call operator
  ///
  /// @tparam generator_t is a random number generator type
  /// @tparam detector_t is the detector information type
  /// @tparam particle_t is the particle information type
  ///
  /// @param[in] generator is the random number generator
  /// @param[in] detector the detector information
  /// @param[in] particle the particle which is being scattered
  ///
  /// @return eventually produced photons
  template <typename generator_t, typename detector_t, typename particle_t>
  std::vector<particle_t> operator()(generator_t &generator,
                                     const detector_t &detector,
                                     particle_t &particle) const
	{
		return hadronicInteraction(generator, detector, particle);
	}
                                     
private:

/// @brief Calculates the absorption length for various hadrons
///
/// @tparam material_t data type of the material
/// @tparam particle_t data type of the particle
/// @param [in] material material that is penetrated
/// @param [in] particle particle that penetrates the material
///
/// @return absorption length
template <typename material_t, typename particle_t>
double 
absorptionLength(const material_t* matertial, particle_t& particle) const;

/// @brief Dices the number of particle candidates that leave the detector
///
/// @tparam generator_t data type of the random number generator
/// @tparam particle_t data type of the particle
/// @param [in] generator random number generator
/// @param [in] particle ingoing particle
/// @return number of outgoing candidates
template<typename generator_t, typename particle_t>
int
diceNumberOfParticles(generator_t& generator, particle_t& particle) const;

/// @brief Creates the particle candidates that leave the detector
///
/// @tparam generator_t data type of the random number generator
/// @tparam particle_t data type of the particle
/// @param [in] generator random number generator
/// @param [in] particle ingoing particle
/// @param [in] particles list of created particles
template<typename generator_t, typename particle_t>
void
createMultiplicity(generator_t& generator, particle_t& particle, std::vector<particle_t>& particles) const;

/// @brief Creates the kinematics of a list of particles
///
/// @tparam generator_t data type of the random number generator
/// @tparam particle_t data type of the particle
/// @param [in] generator random number generator
/// @param [in] particles list of created particles
template<typename generator_t, typename particle_t>
void
kinematics(generator_t& generator, std::vector<particle_t>& particles) const;

/// @brief Selects created particles if they fulfill the criteria
///
/// @tparam particle_t data type of the particle
/// @param [in] particles list of created particles
template<typename particle_t>
void
selectionOfCollection(std::vector<particle_t>& particles) const;

/// @brief Calculates the hadron interactions of a particle
///
/// @tparam generator_t data type of the random number generator
/// @tparam particle_t data type of the particle
/// @param [in] generator random number generator
/// @param [in] time starting time of the evolution
/// @param [in] particle particle that interacts
///
/// @return vector of outgoing particles
template<typename generator_t, typename particle_t>
std::vector<particle_t> 
getHadronState(generator_t& generator, double time, particle_t& particle) const;

// TODO: funktion waere geil, die eine vertex list ausgibt
// TODO: renaming von templates
// Configuration storage
Config m_cfg;
};

ParametricNuclearInt::ParametricNuclearInt(ParametricNuclearInt::Config& cfg) : m_cfg(cfg)
{
}

template <typename material_t, typename particle_t>
double 
ParametricNuclearInt::absorptionLength(const material_t* matertial, particle_t& particle) const 
{
  double al = material->l0();

  if(particle.pdg == 211 || particle.pdg == -211 || particle.pdg == 321 || particle.pdg == -321 || particle.pdg == 111 || particle.pdg == 311)
    al *= 1. / (1. + exp(-(particle.p - 270.) * (particle.p - 270.) / 7200.)); // TODO: da kann man sicherlich noch etwas optimieren

  if(particle.pdg == 2212 || particle.pdg == 2112) al *= 0.7;
  if(particle.pdg == 211 || particle.pdg == -211 || particle.pdg == 111) al *= 0.9;

  return al;
}

template<typename generator_t, typename particle_t>
int
ParametricNuclearInt::diceNumberOfParticles(generator_t& generator, particle_t& particle) const
{
  // sampling of hadronic interaction
  double E = sqrt(particle.p * particle.p + particle.m * particle.m);
  // get the maximum multiplicity    
  double multiplicity_max = 0.25 * E / 1000. + 18.;
  // multiplicity distribution
  double randx, randy, arg = 0.;
  double p1 = 0.;
  double p2 = 0.;
  if (E > 15000.) 
  {
    p1 = 8.69;
    p2 = 2.34;
  } 
  else 
  {
    p1 = 6.77;
    p2 = 2.02;
  }
  
  for (;;) {
    randx = 30. * generator();
    randy = 1.  * generator();
    arg = exp(-0.5 * ((randx - p1) / p2 + exp(-(randx - p1) / p2)));
    if (randy < arg && randx > 3 && randx < multiplicity_max) break;
  }
  
  randx *= (1.2 - 0.4 * exp(-0.001 * p));     // trying to adjust
  
  return (int)randx;
}

template<typename generator_t, typename particle_t>
void
ParametricNuclearInt::createMultiplicity(generator_t& generator, particle_t& particle, std::vector<particle_t>& particles) const
{    
  // new sampling: sample particle type and energy in the CMS frame of outgoing particles
  // creation of shower particles
  double chargedist = 0.;
  
  // sample heavy particles (alpha) but don't save  
  double pif = 0.10; 
  double nef = 0.30;
  double prf = 0.30;
  
  if(particle.pdg == 211 || particle.pdg == -211 || particle.pdg == 321 || particle.pdg == -321 || particle.pdg == 111 || particle.pdg == 311) 
  { // TODO: hier muss wohl ein else fuer andere particles her
      pif = 0.15;
      nef = 0.25;
      prf = 0.25;
  }
  if(particle.pdg == 2212 ) 
  {
    pif = 0.06;
    nef = 0.25;
    prf = 0.35;
  }
  if(particle.pdg == 2112 ) 
  {
    pif = 0.03;
    nef = 0.35;
    prf = 0.17;
  }
  
  // Source of masses: Geant4
  for(int i=0; i<Npart; i++) {
    chargedist  = generator();
    if(chargedist < pif) 
    {
		particles[i].q = 0.
		particles[i].pdg = 111;
		particles[i].m = 0.1349766 * Acts::units::_GeV;
		continue;
    }
    if(chargedist < 2 * pif) 
    {
		particles[i].q = Acts::units::_e;
		particles[i].pdg = 211;
		particles[i].m = 0.1395701 * Acts::units::_GeV;
		continue;
    }
    if(chargedist < 3 * pif) 
    {
		particles[i].q = -Acts::units::_e;
		particles[i].pdg = -211;
		particles[i].m = 0.1395701 * Acts::units::_GeV;
		continue;
    }
    if(chargedist < 3 * pif + nef) 
    {
		particles[i].q = 0.;
		particles[i].pdg = 2112;
		particles[i].m = 939.56563 * Acts::units::_MeV;
		continue;
    }
    if(chargedist < 3 * pif + nef + prf) 
    {
		particles[i].q = Acts::units::_e;
		particles[i].pdg = 2212;
		particles[i].m = 938.27231 * Acts::units::_MeV;
		continue;
    }
    particles[i].q = 2.;
    particles[i].pdg = 20000;
    particles[i].m = 4. * Acts::units::_GeV;
  }

  // move the incoming particle type forward
  if(particles[0].pdg != particle.pdg) 
    for(int i = 1; i < particles.size(); i++)
      if(particles[i].pdg == particle.pdg)
      {
        particles[i] = particles[0];
        particles[0] = particle;
        break;
      }
}

template<typename generator_t, typename particle_t>
void
ParametricNuclearInt::kinematics(generator_t& generator, std::vector<particle_t>& particles, particle_t& particle) const
{
	unsigned int Npart = particles.size();
  std::vector<double> mom(Npart);
  std::vector<double> th(Npart);
  std::vector<double> ph(Npart);

  // sample first particle energy fraction and random momentum direction
  double eps = 2. / Npart;
  double rnd  = generator();
  mom[0] = 0.5 * pow(eps, rnd);          
  th[0]  = acos(2 * generator() - 1.);
  ph[0]  = 2 * M_PI * generator();
  
  // toss particles around in a way which preserves the total momentum (0.,0.,0.) at this point
  // TODO shoot first particle along the impact direction preferentially

  Acts::Vector3D ptemp(mom[0]*sin(th[0])*cos(ph[0]),mom[0]*sin(th[0])*sin(ph[0]),mom[0]*cos(th[0]));
  double ptot = mom[0];
  
  double theta = 0.; double phi = 0.; 
  for (int i = 1; i < Npart - 2; i++) 
  {
    eps = 1. / (Npart - i); 
    mom[i] = (eps + generator() * (1 - eps)) * (1 - ptot); 
    if(ptemp.mag() < 1 - ptot) 
      while(fabs(ptemp.mag() - mom[i]) > 1 - ptot - mom[i])
		mom[i] = (eps + generator() * (1 - eps)) * (1 - ptot);
    
    // max p remaining
    double p_rem = 1 - ptot-mom[i];
    double cthmax = fmin(1.,(-ptemp.mag()*ptemp.mag()-mom[i]*mom[i]+p_rem*p_rem)/2/ptemp.mag()/mom[i]);
    double rnd  = generator();
    theta = acos( (cthmax+1.)*rnd-1.);          
    phi = 2 * M_PI * generator();
    Acts::Vector3D test(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta));
    Acts::RotationMatrix3D rotY, rotZ;
    rotY << cos(ptemp.theta()) << 0. << sin(ptemp.theta())
		 << 0. << 1. << 0.
		 << -sin(ptemp.theta()) << 0. << cos(ptemp.theta());
	rotZ << cos(ptemp.phi()) << -sin(ptemp.phi()) << 0.
		 << sin(ptemp.phi() << cos(ptemp.phi()) << 0.
		 << 0. << 0. << 1.;
    Acts::Vector3D dnewHep = rotZ * rotY * test;
    Acts::Vector3D dnew(dnewHep.x(), dnewHep.y(), dnewHep.z());
    th[i] = dnew.theta();    
    ph[i] = dnew.phi();          
    ptemp += mom[i] * dnew;
    ptot += mom[i];
  }
  
  eps = 0.5; 
  mom[Npart-2] = pow(eps, generator()) * (1 - ptot);
  mom[Npart-1] = 1 - ptot - mom[Npart - 2];
  
  if(ptemp.mag() < 1 - ptot) 
    while(mom[Npart-1]+mom[Npart-2]<ptemp.mag()) 
    { 
      mom[Npart-2] = pow(eps, generator()) * (1 - ptot);
      mom[Npart-1] = 1 - ptot - mom[Npart - 2];
    }
    
  if (ptemp.mag()<fabs(mom[Npart-1]-mom[Npart-2]) ) {
    double diff = ptemp.mag() * generator();
    double sum = mom[Npart - 1] - mom[Npart - 2];
    mom[Npart - 2] = 0.5 * (sum + diff);  
    mom[Npart - 1] = 0.5 * (sum - diff);  
  }
  double cth =(-ptemp.mag()*ptemp.mag()-mom[Npart-2]*mom[Npart-2]+mom[Npart-1]*mom[Npart-1])/2/ptemp.mag()/mom[Npart-2];
  if (fabs(cth)>1.) cth = (cth>0.) ? 1. : -1.;
  
  theta = acos(cth);
  phi = 2 * M_PI * generator();
  Acts::Vector3D test(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta));
    Acts::RotationMatrix3D rotY, rotZ;
    rotY << cos(ptemp.theta()) << 0. << sin(ptemp.theta())
		 << 0. << 1. << 0.
		 << -sin(ptemp.theta()) << 0. << cos(ptemp.theta());
	rotZ << cos(ptemp.phi()) << -sin(ptemp.phi()) << 0.
		 << sin(ptemp.phi() << cos(ptemp.phi()) << 0.
		 << 0. << 0. << 1.;
  Acts::Vector3D dnewHep = rotZ * rotY * test;
  Acts::Vector3D dnew(dnewHep.x(), dnewHep.y(), dnewHep.z());
  
  th[Npart - 2]=dnew.theta();    
  ph[Npart - 2]=dnew.phi();    
  ptemp += mom[Npart - 2] * dnew;
  Acts::Vector3D dlast = -ptemp;
  th[Npart - 1] = dlast.theta(); 
  ph[Npart - 1] = dlast.phi();
  
  
  // particle sampled, rotate, boost and save final state
  double etot = 0.;
  for (int i=0;i<Npart; i++) 
	etot += sqrt(mom[i]*mom[i]+newm[i]*newm[i]);
  double summ = 0.;
  for (int i=0;i<Npart; i++) 
	summ += newm[i];

  // rescale (roughly) to the expected energy
  float scale = sqrt(summ*summ+2*summ*p+m*m)/etot;
  etot = 0.;
  for (int i=0;i<Npart; i++) {
    mom[i] *= scale;
    etot += sqrt(mom[i]*mom[i]+newm[i]*newm[i]);
  }
  
  // Source: http://www.apc.univ-paris7.fr/~franco/g4doxy4.10/html/_lorentz_vector_8cc_source.html - boostvector()
  Acts::Vector3D bv = particle.momentum / sqrt(etot*etot+p*p); // TODO: Why such an energy term?

  std::vector<CLHEP::HepLorentzVector> childBoost(Npart);
  
  for (int i = 0; i < Npart; i++) 
  {
    Acts::Vector3D dirCms(sin(th[i])*cos(ph[i]),sin(th[i])*sin(ph[i]),cos(th[i])); 
    particles[i].momentum = mom[i] * dirCms;
    particles[i].E = sqrt(mom[i]*mom[i]+newm[i]*newm[i]);
	particles[i].boost(bv);
  }   
}

template<typename particle_t>
void
ParametricNuclearInt::selectionOfCollection(std::vector<particle_t> particles) const
{
  // child particle vector for TruthIncident
  //  Reserve space for as many paricles as created due to hadr. int.
  //  However, the number of child particles for TruthIncident might be
  //  smaller due to (momentum) cuts
  unsigned int m_hadIntChildren = 0;
  
  for (int i=0; i<Npart; i++) {
    if (particles[i].pdg > 10000 || childP.mag() < m_cfg.m_minimumHadOutEnergy)
		particles.erase(particle.begin() + i);
  } // particle loop
}

template<typename generator_t, typename particle_t>
std::vector<particle_t> 
ParametricNuclearInt::getHadronState(generator_t& generator, particle_t& particle) const
{  
	std::vector<particle_t> chDef; 
  
  // Calculate multiplicity
  int Npart = diceNumberOfParticles(generator, particle);
  
  // protection against Npart < 3
  if (Npart < 3)
    return chDef;
  else
	chDef.resize(Npart);
	
	createMultiplicity(generator, particle, chDef);
	
	kinematics(generator, chDef, particle);
	selectionOfCollection(chDef);
  
  return chDef;
}

template <typename generator_t, typename material_t, typename particle_t>
std::vector<particle_t> 
ParametricNuclearInt::hadronicInteraction(generator_t& generator, const material_t& material, double pathCorrection, particle_t& particle) const
{
	const material_t* extMprop = dynamic_cast<const material_t*>(&material);
	double prob = 0.;

	// m_hadIntProbScale is used later, not here
	if (extMprop && !m_cfg.m_hadronInteractionFromX0) 
	{	  
		double al = absorptionLength(extMprop, particle);  // in mm
	
	    if (al > 0.) 
			prob = exp(-pathCorrection * extMprop->thickness() / al);
	    else
			prob = exp(-pathCorrection * extMprop->thicknessInL0());
	} 
	else 
		// using approximation lambda = 0.37 * Z * X0 instead -- giving a warning message
		prob = exp(-pathCorrection * mprop.thicknessInX0() / (0.37 * mprop.averageZ()));
  
	// apply a global scalor of the probability
	// (1. - prob) is generally O(0.01), so this is the right way to scale it
	if (generator() < (1. - prob) * m_hadronInteractionProbabilityScale) 
		return getHadronState(generator, particle);
  
	// no hadronic interactions were computed
	return {};  
}

} // namespace Fatras
