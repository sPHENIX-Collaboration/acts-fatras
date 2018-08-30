// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"
#include <math.h>
// TODO: explicit calculations may become an interface
namespace Fatras {

/// The struct for the physics list
///
struct ParametricNuclearInt {

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
                                     particle_t &particle) const;
                                
protected:

/// @brief Calculates and states if a nuclear interaction occurs
///
/// @tparam generator_t data type of the random number generator
/// @tparam material_t data type of the material
/// @tparam particle_t data type of the particle
/// @param [in] generator is the random number generator
/// @param [in] material material that is penetrated
/// @param [in] particle particle that penetrates the material
///
/// @return boolean result if a nuclear interaction occurs
template <typename generator_t, typename material_t, typename particle_t>
bool 
nuclearInteraction(generator_t& generator, const material_t& matertial, particle_t& particle) const;

/// @brief Calculates the probability of a nuclear interaction
///
/// @tparam generator_t data type of the random number generator
/// @param [in] thickness Thickness of the material in terms of L0
/// @param [in] momentum Particles momentum
/// @param [in] pdg PDG code of the particle
///
/// @return boolean result if a nuclear interaction occurs
double 
nuclearInteractionProb(const double thickness, const double momentum, const int pdg) const;

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
kinematics(generator_t& generator, std::vector<particle_t>& particles, particle_t& particle) const;

/// @brief Calculates the hadron interactions of a particle
///
/// @tparam generator_t data type of the random number generator
/// @tparam particle_t data type of the particle
/// @param [in] generator random number generator
/// @param [in] particle particle that interacts
///
/// @return vector of outgoing particles
template<typename generator_t, typename particle_t>
std::vector<particle_t> 
getHadronState(generator_t& generator, particle_t& particle) const;
};

template <typename generator_t, typename material_t, typename particle_t>
bool 
ParametricNuclearInt::nuclearInteraction(generator_t& generator, const material_t& material, particle_t& particle) const 
{
	double prob = nuclearInteractionProb(extMprop->thickness() / material.averageL0(), 
  double al = extMprop->thickness() / material.averageL0();

	if(particle.pdg() == 211 || particle.pdg() == -211 || particle.pdg() == 321 || particle.pdg() == -321 || particle.pdg() == 111 || particle.pdg() == 311)
		al *= 1. / (1. + exp(-(particle.p() - 270.) * (particle.p() - 270.) / 7200.));

	if(particle.pdg() == 2212 || particle.pdg() == 2112) al *= 0.7;
	if(particle.pdg() == 211 || particle.pdg() == -211 || particle.pdg() == 111) al *= 0.9;

	const material_t* extMprop = &material;
	double prob = 0.;

    if (al > 0.) 
		prob = exp(-extMprop->thickness() / al);
    else
		prob = exp(-extMprop->averageL0());

  return generator() < (1. - prob) * 0.5;
}

// TODO: move this function to .cpp
double
ParametricNuclearInt::nuclearInteractionProb(const double thickness, const double momentum, const int pdg) const
{

	//p/n
	//~ //return (1. - exp(-val[1] * 1.1804 * (1.6309e-07 + exp(-(val[0] - 12.833) * (val[0] - 12.833) / 221.47)))) * 0.96211; // volle range
	//~ if(val[0] > 3.)
		//~ return (1. - exp(-val[1] * 0.98835 * (0.07491 + exp(-(val[0] - 13.744) * (val[0] - 13.744) / 805.25)))) * 0.97985; // >3
	//~ else
		//~ if(val[1] > 0.5)
			//~ return 0.05703 + (1. - exp(-val[1] * 1.1889 * (1.1714e-06 + exp(-(val[0] - 2.497) * (val[0] - 2.497) / 2.0915)))) * 0.75229; // < 3
		//~ else
			//~ return 0.009525 + (1. - exp(-val[1] * 1.4277 * (0.011844 + exp(-(val[0] - 2.3332) * (val[0] - 2.3332) / 3.9307)))) * 0.71596;

	//pi+/pi-
	//~ return (1. - exp(-val[1] * 0.87514 * (4.8592e-06 + exp(-(val[0] - 12.05) * (val[0] - 12.05) / 483.61)))) * 0.935931; // volle range
	//~ if(val[0] > 4.)
		//~ return (1. - exp(-val[1] * 0.55001 * (0.46996 + exp(-(val[0] - 8.5732) * (val[0] - 8.5732) / 1800.)))) * 0.96459; // > 4
	//~ else 
		//~ if(val[1] > 0.5)
			//~ return 0.03192 + (1. - exp(-val[1] * 1.0276 * (0.028738 + exp(-(val[0] - 3.1018) * (val[0] - 3.1018) / 8.866)))) * 0.76954; // < 4
		//~ else
			//~ return 0.003446 + (1. - exp(-val[1] * 1.0482 * (0.074884 + exp(-(val[0] - 2.277) * (val[0] - 2.277) / 22.708)))) * 0.74133;

	//kaon+
	//~ return (1. - exp(-val[1] * 0.86505 * (1.5586e-07 + exp(-(val[0] - 11.243) * (val[0] - 11.243) / 577.7)))) * 0.84876; // volle range
	//~ if(val[0] > 4.)
		//~ return (1. - exp(-val[1] * 0.13597 * (4.8514 + exp(-(val[0] - 8.1613e-05) * (val[0] - 8.1613e-05) / 771.03)))) * 0.88989; // > 4
	//~ else 
		//~ if(val[1] > 0.5)
			//~ return 0.009418 + (1. - exp(-val[1] * 0.83191 * (7.8499e-10 + exp(-(val[0] - 2.9685) * (val[0] - 2.9685) / 4.0707)))) * 0.74356; // < 4
		//~ else
			//~ return -0.004841 + (1. - exp(-val[1] * 1.5404 * (0.00046885 + exp(-(val[0] - 2.8789) * (val[0] - 2.8789) / 8.8744)))) * 0.44758;
			
	//kaon-
	//~ if(val[0] > 4.)
		//~ return (1. - exp(-val[1] * 0.13597 * (4.8514 + exp(-(val[0] - 8.1613e-05) * (val[0] - 8.1613e-05) / 771.03)))) * 0.88989; // > 4
	//~ else 
		//~ if(val[1] > 0.5)
			//~ return 0.0004233 + (1. - exp(-val[1] * 0.70435 * (0.40021 + exp(-(val[0] - 5.8872e-05) * (val[0] - 5.8872e-05) / 107.32)))) * 0.82217; // < 4
		//~ else
			//~ return 0.003083 + (1. - exp(-val[1] * 3.6515 * (0.14931 + exp(-(val[0] - 88.407) * (val[0] - 88.407) / 15.98)))) * 1.3231;
			
	//kaon0
	if(val[0] > 4.)
		if(val[1] < 0.5 && val[0] < 12.)
			return (1. - exp(-val[1] * 0.62966 * (1.8493 + exp(-(val[0] - 8.9774e-08) * (val[0] - 8.9774e-08) / 44.93)))) * 0.72643; // > 4
		else
			return (1. - exp(-val[1] * 0.13597 * (4.8514 + exp(-(val[0] - 8.1613e-05) * (val[0] - 8.1613e-05) / 771.03)))) * 0.88989;
	else 
		if(val[1] > 0.5)
			return 0.2789 + (1. - exp(-val[1] * 0.98166 * (6.3897e-08 + exp(-(val[0] - 3.1151) * (val[0] - 3.1151) / 4.08)))) * 0.46121; // < 4
		else
			return 0.01415 + (1. - exp(-val[1] * 3.8071 * (0.59535 + exp(-(val[0] - 0.0004938) * (val[0] - 0.0004938) / 10.21)))) * 0.42625;
			
}

template<typename generator_t, typename particle_t>
int
ParametricNuclearInt::diceNumberOfParticles(generator_t& generator, particle_t& particle) const
{
  // sampling of hadronic interaction
  double E = sqrt(particle.p() * particle.p() + particle.m() * particle.m());
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
  
  randx *= (1.2 - 0.4 * exp(-0.001 * particle.p()));     // trying to adjust

  return (int)randx;
}

template<typename generator_t, typename particle_t>
void
ParametricNuclearInt::createMultiplicity(generator_t& generator, particle_t& particle, std::vector<particle_t>& particles) const
{    
  //~ // new sampling: sample particle type and energy in the CMS frame of outgoing particles
  //~ // creation of shower particles
  //~ double chargedist = 0.;
  
  //~ // sample heavy particles (alpha) but don't save  
  //~ double pif = 0.10; 
  //~ double nef = 0.30;
  //~ double prf = 0.30;
  
  //~ if(particle.pdg() == 211 || particle.pdg() == -211 || particle.pdg() == 321 || particle.pdg() == -321 || particle.pdg() == 111 || particle.pdg() == 311) 
  //~ {
      //~ pif = 0.15;
      //~ nef = 0.25;
      //~ prf = 0.25;
  //~ }
  //~ if(particle.pdg() == 2212) 
  //~ {
    //~ pif = 0.06;
    //~ nef = 0.25;
    //~ prf = 0.35;
  //~ }
  //~ if(particle.pdg() == 2112) 
  //~ {
    //~ pif = 0.03;
    //~ nef = 0.35;
    //~ prf = 0.17;
  //~ }
  
  //~ // Source of masses: Geant4
  //~ for(unsigned int i = 0; i < particles.size(); i++) {
    //~ chargedist  = generator();
    //~ if(chargedist < pif) 
    //~ {
		//~ particles[i].q() = 0.;
		//~ particles[i].pdg() = 111;
		//~ particles[i].m() = 0.1349766 * Acts::units::_GeV;
		//~ continue;
    //~ }
    //~ if(chargedist < 2 * pif) 
    //~ {
		//~ particles[i].q() = Acts::units::_e;
		//~ particles[i].pdg() = 211;
		//~ particles[i].m() = 0.1395701 * Acts::units::_GeV;
		//~ continue;
    //~ }
    //~ if(chargedist < 3 * pif) 
    //~ {
		//~ particles[i].q() = -Acts::units::_e;
		//~ particles[i].pdg() = -211;
		//~ particles[i].m() = 0.1395701 * Acts::units::_GeV;
		//~ continue;
    //~ }
    //~ if(chargedist < 3 * pif + nef) 
    //~ {
		//~ particles[i].q() = 0.;
		//~ particles[i].pdg() = 2112;
		//~ particles[i].m() = 939.56563 * Acts::units::_MeV;
		//~ continue;
    //~ }
    //~ if(chargedist < 3 * pif + nef + prf) 
    //~ {
		//~ particles[i].q() = Acts::units::_e;
		//~ particles[i].pdg() = 2212;
		//~ particles[i].m() = 938.27231 * Acts::units::_MeV;
		//~ continue;
    //~ }
    //~ particles[i].q() = 2.;
    //~ particles[i].pdg() = 20000;
    //~ particles[i].m() = 4. * Acts::units::_GeV;
  //~ }

  //~ // move the incoming particle type forward
  //~ if(particles[0].pdg() != particle.pdg()) 
    //~ for(unsigned int i = 1; i < particles.size(); i++)
      //~ if(particles[i].pdg() == particle.pdg())
      //~ {
        //~ particles[i] = particles[0];
        //~ particles[0] = particle;
        //~ break;
      //~ }
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
  for (unsigned int i = 1; i < Npart - 2; i++) 
  {
    eps = 1. / (Npart - i); 
    mom[i] = (eps + generator() * (1 - eps)) * (1 - ptot); 
    if(ptemp.mag() < 1 - ptot) 
      while(fabs(ptemp.mag() - mom[i]) > 1 - ptot - mom[i])
		mom[i] = (eps + generator() * (1 - eps)) * (1 - ptot);
    
    // max p remaining
    double p_rem = 1 - ptot-mom[i];
    double cthmax = fmin(1.,(-ptemp.mag()*ptemp.mag()-mom[i]*mom[i]+p_rem*p_rem)/2/ptemp.mag()/mom[i]);
    rnd  = generator();
    theta = acos( (cthmax+1.)*rnd-1.);          
    phi = 2 * M_PI * generator();
    Acts::Vector3D test(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta));
    Acts::RotationMatrix3D rotY, rotZ;
    rotY << cos(ptemp.theta()), 0., sin(ptemp.theta()),
		0., 1., 0.,
		-sin(ptemp.theta()), 0., cos(ptemp.theta());
	rotZ << cos(ptemp.phi()), -sin(ptemp.phi()), 0.,
		sin(ptemp.phi()), cos(ptemp.phi()), 0.,
		0., 0., 1.;
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
    rotY << cos(ptemp.theta()), 0., sin(ptemp.theta()),
		 0., 1., 0.,
		 -sin(ptemp.theta()), 0., cos(ptemp.theta());
	rotZ << cos(ptemp.phi()), -sin(ptemp.phi()), 0.,
		 sin(ptemp.phi()), cos(ptemp.phi()), 0.,
		 0., 0., 1.;
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
  for (unsigned int i = 0; i < Npart; i++) 
	etot += sqrt(mom[i] * mom[i] + particles[i].m() * particles[i].m());
  double summ = 0.;
  for (unsigned int i = 0; i < Npart; i++) 
	summ += particles[i].m();

  // rescale (roughly) to the expected energy
  float scale = sqrt(summ*summ+2*summ*particle.p()+particle.m() * particle.m())/etot;
  etot = 0.;
  for (unsigned int i = 0; i < Npart; i++) {
    mom[i] *= scale;
    etot += sqrt(mom[i] * mom[i] + particles[i].m() * particles[i].m());
  }
  
  // Source: http://www.apc.univ-paris7.fr/~franco/g4doxy4.10/html/_lorentz_vector_8cc_source.html - boostvector()
  Acts::Vector3D bv = particle.momentum() / sqrt(etot * etot + particle.p() * particle.p()); // TODO: Why such an energy term?
  
  /// New particle structure does not allow setting values in the old way
  //~ for (unsigned int i = 0; i < Npart; i++) 
  //~ {
    //~ Acts::Vector3D dirCms(sin(th[i])*cos(ph[i]),sin(th[i])*sin(ph[i]),cos(th[i])); 
    //~ particles[i].momentum = mom[i] * dirCms;
    //~ particles[i].p = particles[i].momentum.mag();
    //~ particles[i].E = sqrt(mom[i] * mom[i] + particles[i].m * particles[i].m);
	//~ particles[i].boost(bv);
  //~ }
}

template<typename generator_t, typename particle_t>
std::vector<particle_t> 
ParametricNuclearInt::getHadronState(generator_t& generator, particle_t& particle) const
{  
	std::vector<particle_t> chDef; 
  
  // Calculate multiplicity
  int Npart = diceNumberOfParticles(generator, particle);
  
  // protection against Npart < 3
  if (Npart >= 3)
	{
	chDef.resize(Npart);

	createMultiplicity(generator, particle, chDef);
	
	kinematics(generator, chDef, particle);
  }
  else
	chDef.push_back(particle); // Return at least the leading particle
  return chDef;
}

template <typename generator_t, typename detector_t, typename particle_t>
std::vector<particle_t> ParametricNuclearInt::operator()(generator_t &generator,
                                     const detector_t &detector,
                                     particle_t &particle) const
{
	if (nuclearInteraction(generator, detector, particle))
		return getHadronState(generator, particle);
 
	// no hadronic interactions were computed
	return {particle};  // Return the incoming particle
}
	
} // namespace Fatras
