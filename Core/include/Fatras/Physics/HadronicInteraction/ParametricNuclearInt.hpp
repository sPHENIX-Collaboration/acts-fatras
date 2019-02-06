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
#include <vector>
#include <array>

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

/// @brief Calculates the probability of a nuclear interaction
///
/// @param [in] momentum Particles momentum in GeV
/// @param [in] thickness Thickness of the material in terms of L0
/// @param [in] pdg PDG code of the particle
///
/// @return boolean result if a nuclear interaction occurs
double 
nuclearInteractionProb(const double momentum, const double thickness, const int pdg) const;

/// @brief Calculates the probability for a certain multiplicity
///
/// @param [in] thickness Thickness of the material in terms of L0
/// @param [in] momentum Particles momentum in GeV
/// @param [in] pdg PDG code of the particle
/// @param [in] mult Multiplicity of the final state
///
/// @return Probability of the given configuration
double
multiplicityProb(const double momentum, const double thickness, const int pdg, const unsigned int mult) const;

/// @brief Dices the number of particle candidates that leave the detector
///
/// @tparam generator_t data type of the random number generator
/// @tparam material_t data type of the material
/// @tparam particle_t data type of the particle
/// @param [in] generator random number generator
/// @param [in] thickness Thickness of the material in terms of L0
/// @param [in] particle ingoing particle
///
/// @return number of outgoing candidates
template<typename generator_t, typename particle_t>
unsigned int
multiplicity(generator_t& generator, const double thickness, particle_t& particle) const;

/// @brief Creates the particle types that leave the detector
///
/// @tparam generator_t data type of the random number generator
/// @tparam particle_t data type of the particle
/// @param [in] generator random number generator
/// @param [in] pdg PDG code of the ingoing particle
/// @param [in] nParticles Number of outgoing particles
///
/// @return list of outgoing particle PDGs
template<typename generator_t>
std::vector<int>
particleComposition(generator_t& generator, const int pdg, const unsigned int nParticles) const;

/// @brief Creates the kinematics of a list of particles
///
/// @tparam generator_t data type of the random number generator
/// @tparam particle_t data type of the particle
/// @param [in] generator random number generator
/// @param [in] particle incoming particle
/// @param [in] particlePDGs list of created particle PDGs
///
/// @return vector of outgoing particles
template<typename generator_t, typename particle_t>
std::vector<particle_t>
kinematics(generator_t& generator, particle_t& particle, const std::vector<int>& particlesPDGs) const;

/// @brief Calculates the hadron interactions of a particle
///
/// @tparam generator_t data type of the random number generator
/// @tparam particle_t data type of the particle
/// @param [in] generator random number generator
/// @param [in] thickness Thickness of the material in terms of L0
/// @param [in] particle particle that interacts
///
/// @return vector of outgoing particles
template<typename generator_t, typename particle_t>
std::vector<particle_t> 
finalStateHadrons(generator_t& generator, const double thicknessInL0, particle_t& particle) const;
};

template<typename generator_t, typename particle_t>
unsigned int
ParametricNuclearInt::multiplicity(generator_t& generator, const double thickness, particle_t& particle) const
{
	const double dice = generator();
	size_t mult = 2;
	double cumulativeProb = 0.;
	
	// Adding probabilites of the multiplicities
	while(dice > cumulativeProb)
	{
		cumulativeProb += multiplicityProb(particle.p(), thickness, particle.pdg(), mult);
		mult++;
	}
	
	return --mult;
}

template<typename generator_t>
std::vector<int>
ParametricNuclearInt::particleComposition(generator_t& generator, const int pdg, const unsigned int nParticles) const
{    	
	//~ std::vector<int> result;
	//~ result.reserve(nParticles);
	//~ double dice;
	//~ unsigned int index;
	
	//~ switch(pdg)
	//~ {
		//~ //k-
		//~ case -321:
		//~ {
			//~ while(result.size() < nParticles)
			//~ {
				//~ dice = generator();
				//~ for(index = 0; index < 8; index++) // 8 particles can be produced
				//~ {
					//~ if(cProbsKm[index] > dice)
					//~ {
						//~ result.push_back(pdgLookUp[index]);
						//~ break;
					//~ }	
				//~ }
			//~ }
			//~ break;
		//~ }
		
		//~ //pi-
		//~ case -211:
		//~ {
			//~ while(result.size() < nParticles)
			//~ {
				//~ dice = generator();
				//~ for(index = 0; index < 8; index++) // 8 particles can be produced
				//~ {
					//~ if(cProbsPim[index] > dice)
					//~ {
						//~ result.push_back(pdgLookUp[index]);
						//~ break;
					//~ }	
				//~ }
			//~ }
			//~ break;
		//~ }
		
		//~ //pi+
		//~ case 211:
		//~ {
			//~ while(result.size() < nParticles)
			//~ {
				//~ dice = generator();
				//~ for(index = 0; index < 8; index++) // 8 particles can be produced
				//~ {
					//~ if(cProbsPip[index] > dice)
					//~ {
						//~ result.push_back(pdgLookUp[index]);
						//~ break;
					//~ }	
				//~ }
			//~ }
			//~ break;
		//~ }
		
		//~ //k0
		//~ case 130:
		//~ case 310:
		//~ case 311:
		//~ {
			//~ while(result.size() < nParticles)
			//~ {
				//~ dice = generator();
				//~ for(index = 0; index < 8; index++) // 8 particles can be produced
				//~ {
					//~ if(cProbsK0[index] > dice)
					//~ {
						//~ result.push_back(pdgLookUp[index]);
						//~ break;
					//~ }	
				//~ }
			//~ }
			//~ break;
		//~ }
		
		//~ //k+
		//~ case 321:
		//~ {
			//~ while(result.size() < nParticles)
			//~ {
				//~ dice = generator();
				//~ for(index = 0; index < 8; index++) // 8 particles can be produced
				//~ {
					//~ if(cProbsKp[index] > dice)
					//~ {
						//~ result.push_back(pdgLookUp[index]);
						//~ break;
					//~ }	
				//~ }
			//~ }
			//~ break;
		//~ }
		
		//~ //n
		//~ case 2112:
		//~ {
			//~ while(result.size() < nParticles)
			//~ {
				//~ dice = generator();
				//~ for(index = 0; index < 8; index++) // 8 particles can be produced
				//~ {
					//~ if(cProbsN[index] > dice)
					//~ {
						//~ result.push_back(pdgLookUp[index]);
						//~ break;
					//~ }	
				//~ }
			//~ }
			//~ break;
		//~ }
		
		//~ //p
		//~ case 2212:
		//~ {
			//~ while(result.size() < nParticles)
			//~ {
				//~ dice = generator();
				//~ for(index = 0; index < 8; index++) // 8 particles can be produced
				//~ {
					//~ if(cProbsP[index] > dice)
					//~ {
						//~ result.push_back(pdgLookUp[index]);
						//~ break;
					//~ }	
				//~ }
			//~ }
		//~ }
	//~ }

	//~ return result;


	return {};
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
std::vector<particle_t>
ParametricNuclearInt::kinematics(generator_t& generator, particle_t& particle, const std::vector<int>& particlesPDGs) const
{
	// TODO: whole function
	// TODO: boost should be a part of this class
	
	//~ unsigned int Npart = particles.size();
  //~ std::vector<double> mom(Npart);
  //~ std::vector<double> th(Npart);
  //~ std::vector<double> ph(Npart);

  //~ // sample first particle energy fraction and random momentum direction
  //~ double eps = 2. / Npart;
  //~ double rnd  = generator();
  //~ mom[0] = 0.5 * pow(eps, rnd);          
  //~ th[0]  = acos(2 * generator() - 1.);
  //~ ph[0]  = 2 * M_PI * generator();
  
  //~ // toss particles around in a way which preserves the total momentum (0.,0.,0.) at this point
  //~ // TODO shoot first particle along the impact direction preferentially

  //~ Acts::Vector3D ptemp(mom[0]*sin(th[0])*cos(ph[0]),mom[0]*sin(th[0])*sin(ph[0]),mom[0]*cos(th[0]));
  //~ double ptot = mom[0];
  
  //~ double theta = 0.; double phi = 0.; 
  //~ for (unsigned int i = 1; i < Npart - 2; i++) 
  //~ {
    //~ eps = 1. / (Npart - i); 
    //~ mom[i] = (eps + generator() * (1 - eps)) * (1 - ptot); 
    //~ if(ptemp.mag() < 1 - ptot) 
      //~ while(fabs(ptemp.mag() - mom[i]) > 1 - ptot - mom[i])
		//~ mom[i] = (eps + generator() * (1 - eps)) * (1 - ptot);
    
    //~ // max p remaining
    //~ double p_rem = 1 - ptot-mom[i];
    //~ double cthmax = fmin(1.,(-ptemp.mag()*ptemp.mag()-mom[i]*mom[i]+p_rem*p_rem)/2/ptemp.mag()/mom[i]);
    //~ rnd  = generator();
    //~ theta = acos( (cthmax+1.)*rnd-1.);          
    //~ phi = 2 * M_PI * generator();
    //~ Acts::Vector3D test(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta));
    //~ Acts::RotationMatrix3D rotY, rotZ;
    //~ rotY << cos(ptemp.theta()), 0., sin(ptemp.theta()),
		//~ 0., 1., 0.,
		//~ -sin(ptemp.theta()), 0., cos(ptemp.theta());
	//~ rotZ << cos(ptemp.phi()), -sin(ptemp.phi()), 0.,
		//~ sin(ptemp.phi()), cos(ptemp.phi()), 0.,
		//~ 0., 0., 1.;
    //~ Acts::Vector3D dnewHep = rotZ * rotY * test;
    //~ Acts::Vector3D dnew(dnewHep.x(), dnewHep.y(), dnewHep.z());
    //~ th[i] = dnew.theta();    
    //~ ph[i] = dnew.phi();          
    //~ ptemp += mom[i] * dnew;
    //~ ptot += mom[i];
  //~ }
  
  //~ eps = 0.5; 
  //~ mom[Npart-2] = pow(eps, generator()) * (1 - ptot);
  //~ mom[Npart-1] = 1 - ptot - mom[Npart - 2];
  
  //~ if(ptemp.mag() < 1 - ptot) 
    //~ while(mom[Npart-1]+mom[Npart-2]<ptemp.mag()) 
    //~ { 
      //~ mom[Npart-2] = pow(eps, generator()) * (1 - ptot);
      //~ mom[Npart-1] = 1 - ptot - mom[Npart - 2];
    //~ }
    
  //~ if (ptemp.mag()<fabs(mom[Npart-1]-mom[Npart-2]) ) {
    //~ double diff = ptemp.mag() * generator();
    //~ double sum = mom[Npart - 1] - mom[Npart - 2];
    //~ mom[Npart - 2] = 0.5 * (sum + diff);  
    //~ mom[Npart - 1] = 0.5 * (sum - diff);  
  //~ }
  //~ double cth =(-ptemp.mag()*ptemp.mag()-mom[Npart-2]*mom[Npart-2]+mom[Npart-1]*mom[Npart-1])/2/ptemp.mag()/mom[Npart-2];
  //~ if (fabs(cth)>1.) cth = (cth>0.) ? 1. : -1.;
  
  //~ theta = acos(cth);
  //~ phi = 2 * M_PI * generator();
  //~ Acts::Vector3D test(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta));
    //~ Acts::RotationMatrix3D rotY, rotZ;
    //~ rotY << cos(ptemp.theta()), 0., sin(ptemp.theta()),
		 //~ 0., 1., 0.,
		 //~ -sin(ptemp.theta()), 0., cos(ptemp.theta());
	//~ rotZ << cos(ptemp.phi()), -sin(ptemp.phi()), 0.,
		 //~ sin(ptemp.phi()), cos(ptemp.phi()), 0.,
		 //~ 0., 0., 1.;
  //~ Acts::Vector3D dnewHep = rotZ * rotY * test;
  //~ Acts::Vector3D dnew(dnewHep.x(), dnewHep.y(), dnewHep.z());
  
  //~ th[Npart - 2]=dnew.theta();    
  //~ ph[Npart - 2]=dnew.phi();    
  //~ ptemp += mom[Npart - 2] * dnew;
  //~ Acts::Vector3D dlast = -ptemp;
  //~ th[Npart - 1] = dlast.theta(); 
  //~ ph[Npart - 1] = dlast.phi();
  
  
  //~ // particle sampled, rotate, boost and save final state
  //~ double etot = 0.;
  //~ for (unsigned int i = 0; i < Npart; i++) 
	//~ etot += sqrt(mom[i] * mom[i] + particles[i].m() * particles[i].m());
  //~ double summ = 0.;
  //~ for (unsigned int i = 0; i < Npart; i++) 
	//~ summ += particles[i].m();

  //~ // rescale (roughly) to the expected energy
  //~ float scale = sqrt(summ*summ+2*summ*particle.p()+particle.m() * particle.m())/etot;
  //~ etot = 0.;
  //~ for (unsigned int i = 0; i < Npart; i++) {
    //~ mom[i] *= scale;
    //~ etot += sqrt(mom[i] * mom[i] + particles[i].m() * particles[i].m());
  //~ }
  
  //~ // Source: http://www.apc.univ-paris7.fr/~franco/g4doxy4.10/html/_lorentz_vector_8cc_source.html - boostvector()
  //~ Acts::Vector3D bv = particle.momentum() / sqrt(etot * etot + particle.p() * particle.p()); // TODO: Why such an energy term?
  
  /// New particle structure does not allow setting values in the old way
  //~ for (unsigned int i = 0; i < Npart; i++) 
  //~ {
    //~ Acts::Vector3D dirCms(sin(th[i])*cos(ph[i]),sin(th[i])*sin(ph[i]),cos(th[i])); 
    //~ particles[i].momentum = mom[i] * dirCms;
    //~ particles[i].p = particles[i].momentum.mag();
    //~ particles[i].E = sqrt(mom[i] * mom[i] + particles[i].m * particles[i].m);
	//~ particles[i].boost(bv);
  //~ }
  return {};
}

template<typename generator_t, typename particle_t>
std::vector<particle_t> 
ParametricNuclearInt::finalStateHadrons(generator_t& generator, const double thickness, particle_t& particle) const
{
	// Calculate multiplicity
	const unsigned int Npart = multiplicity(generator, thickness, particle);
	
	// Calculate particle types
	const std::vector<int> particlePDGs = particleComposition(generator, particle, Npart);

	// Calculate the kinematics
	return kinematics(generator, particle, particlePDGs);
}

template <typename generator_t, typename detector_t, typename particle_t>
std::vector<particle_t> ParametricNuclearInt::operator()(generator_t& generator,
                                     const detector_t& detector,
                                     particle_t &particle) const
{
	// Test applicable PDG codes
	for(const auto& pc : pdgCodes)
		if(particle.pdg() == pc)
		{
			const double thicknessInL0 = detector.thickness() / detector.averageL0();
			
			// If a nuclear interaction occurs ...
			if (generator() < nuclearInteractionProb(particle.p(), thicknessInL0, particle.pdg()))
				// ... calculate the final state hadrons
				return finalStateHadrons(generator, thicknessInL0, particle);
		 
			break;
		}
	// No hadronic interactions occured
	return {particle};
}
	
} // namespace Fatras
