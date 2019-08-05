// This file is part of the Acts project.
//
// Copyright (C) 2018-2019 Acts project team
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
#include <list>

#include <iostream>
#include <fstream>

std::ofstream ofs;

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
  /// @param [in] generator is the random number generator
  /// @param [in] detector the detector information
  /// @param [in] particle the particle which is being scattered
  ///
  /// @return eventually produced hadrons
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

/// @brief Evaluates the fraction E_{out} / E_{in} of a single outgoing particle
/// @note This function is based on the inverse sampling method and therefore requires a uniform distributed random number in [0,1].
///
/// @param [in] cProb Uniform random number from cumulative probability distribution
/// @param [in] scaling Fit parameter to scale the dependency on the multiplicity
/// @param [in] n Demands the sampling from the distribution of the outgoing particle with the nth highest energy
///
/// @return Energy fraction of the outgoing particle
double
energyFraction(const double cProb, const double scaling, const unsigned int n) const;

/// @brief Samples the energies of the ougoing particles
///
/// @tparam generator_t data type of the random number generator
/// @tparam particle_t data type of the particle
/// @param [in] generator random number generator
/// @param [in] pdg PDG code of the ingoing particle
/// @param [in] nParticles Number of outgoing particles
///
/// @return list of outgoing particle energy fractions E_{out} / E_{in}
template<typename generator_t>
std::vector<double>
energyFractions(generator_t& generator, const int pdg, const unsigned int nParticles) const;

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
	size_t mult = 0;
	double cumulativeProb = 0.;
	double currentProb = 0.;

	// Adding probabilites of the multiplicities
	while(dice > cumulativeProb)
	{
		currentProb = multiplicityProb(particle.p(), thickness, particle.pdg(), mult);
		if(currentProb < 1e-4 || mult > 18) // TODO: are there better cuts?
			return mult;
		
		cumulativeProb += currentProb;
		
		mult++;
	}
	
	return (mult > 0) ? --mult : 0;
}

template<typename generator_t, typename particle_t>
std::vector<particle_t> 
ParametricNuclearInt::finalStateHadrons(generator_t& generator, const double thickness, particle_t& particle) const
{

	// Calculate multiplicity
	const unsigned int Npart = multiplicity(generator, thickness, particle);

	// Easy exit if nothing gets out
	if(Npart == 0)
		return {};

	// Calculate particle types
	const std::vector<int> particlePDGs = particleComposition(generator, particle.pdg(), Npart);

	// Calculate the kinematics
	return kinematics(generator, particle, particlePDGs);
}	
} // namespace Fatras
#include "Fatras/Physics/HadronicInteraction/detail/ParametricNuclearInt.ipp"