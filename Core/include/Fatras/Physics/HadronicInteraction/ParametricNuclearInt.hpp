// This file is part of the Acts project.
//
// Copyright (C) 2018-2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Definitions.hpp"
#include <math.h>
#include <vector>
#include <array>
#include <list>

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
   
private:

/// @brief Calculates the probability of a nuclear interaction
///
/// @param [in] momentum Particles momentum in GeV
/// @param [in] thickness Thickness of the material in terms of L0
/// @param [in] pars Parametrisation 
///
/// @return boolean result if a nuclear interaction occurs
double 
nuclearInteractionProb(const double momentum, const double thickness, const std::array<double, 6>& pars) const;

/// @brief Calculates the probability for a certain multiplicity
///
/// @param [in] thickness Thickness of the material in terms of L0
/// @param [in] momentum Particles momentum in GeV
/// @param [in] pars Parametrisation
/// @param [in] mult Multiplicity of the final state
///
/// @return Probability of the given configuration
double
multiplicityProb(const double momentum, const double thickness, const std::array<double, 7>& pars, const unsigned int mult) const;

/// @brief Dices the number of particle candidates that leave the detector
///
/// @tparam generator_t data type of the random number generator
/// @tparam material_t data type of the material
/// @tparam particle_t data type of the particle
/// @param [in] generator random number generator
/// @param [in] thickness Thickness of the material in terms of L0
/// @param [in] particle ingoing particle
/// @param [in] pars Parametrisation
///
/// @return number of outgoing candidates
template<typename generator_t, typename particle_t>
unsigned int
multiplicity(generator_t& generator, const double thickness, particle_t& particle, const std::array<double, 7>& pars) const;

/// @brief Creates the particle types that leave the detector
///
/// @tparam generator_t data type of the random number generator
/// @tparam particle_t data type of the particle
/// @param [in] generator random number generator
/// @param [in] particleLookUp Look-up table for PID production probabilities to create
/// @param [in] nParticles Number of outgoing particles
///
/// @return list of outgoing particle PDGs
template<typename generator_t>
std::vector<int>
particleComposition(generator_t& generator, const std::list<std::pair<double, int>>& particleLookUp, const unsigned int nParticles) const;

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
/// @tparam parameters_t Type of the container of parametrizations
/// @param [in] generator random number generator
/// @param [in] pdg PDG code of the ingoing particle
/// @param [in] nParticles Number of outgoing particles
/// @param [in] params Parametrisations
///
/// @return list of outgoing particle energy fractions E_{out} / E_{in}
template<typename generator_t, typename parameters_t>
std::vector<double>
energyFractions(generator_t& generator, const parameters_t& params, const unsigned int nParticles) const;

/// @brief Evaluate the probability the receive a certain value for cos(theta)
///
/// @param [in] cosTheta The value of cos(theta)
/// @param [in] fitParameters Fit parameters to evaluate the probability
///
/// @return The probability for a certain value of cos(theta)
double
cosThetaProbability(double cosTheta, const std::array<double, 6>& fitParameters) const;

/// @brief Sample a theta angle given a set of fit parameters for a certain model
///
/// @tparam generator_t data type of the random number generator
/// @param [in] generator random number generator
/// @param [in] fitParameters Fit parameters of the underlying model
///
/// @return Theta angle
template<typename generator_t>
double
sampleTheta(generator_t& generator, const std::array<double, 6>& fitParameters) const;

/// @brief Creates the kinematics of a list of particles
///
/// @tparam generator_t data type of the random number generator
/// @tparam particle_t data type of the particle
/// @tparam parameters_t Type of the container of parametrizations
/// @param [in] generator random number generator
/// @param [in] particle incoming particle
/// @param [in] particlePDGs list of created particle PDGs
/// @param [in] params Parametrisations
///
/// @return vector of outgoing particles
template<typename generator_t, typename particle_t, typename parameters_t>
std::vector<particle_t>
kinematics(generator_t& generator, particle_t& particle, const std::vector<int>& particlesPDGs, const parameters_t& params) const;

/// @brief Calculates the hadron interactions of a particle
///
/// @tparam generator_t data type of the random number generator
/// @tparam particle_t data type of the particle
/// @tparam parameters_t Type of the container of parametrizations
/// @param [in] generator random number generator
/// @param [in] thickness Thickness of the material in terms of L0
/// @param [in] particle particle that interacts
/// @param [in] params Parametrisations of the interaction
///
/// @return vector of outgoing particles
template<typename generator_t, typename particle_t, typename parameters_t>
std::vector<particle_t> 
finalStateHadrons(generator_t& generator, const double thicknessInL0, particle_t& particle, const parameters_t& params) const;
};

template<typename generator_t, typename particle_t>
unsigned int
ParametricNuclearInt::multiplicity(generator_t& generator, const double thickness, particle_t& particle, const std::array<double, 7>& pars) const
{
	const double dice = generator();
	size_t mult = 0;
	double cumulativeProb = 0.;
	double currentProb = 0.;

	// Adding probabilites of the multiplicities
	while(dice > cumulativeProb)
	{
		currentProb = multiplicityProb(particle.p(), thickness, pars, mult);
		if(currentProb < 1e-4 || mult > 18) // TODO: are there better cuts?
			return mult;
		
		cumulativeProb += currentProb;
		
		mult++;
	}
	
	return (mult > 0) ? --mult : 0;
}

template<typename generator_t>
std::vector<int>
ParametricNuclearInt::particleComposition(generator_t& generator, const std::list<std::pair<double, int>>& particleLookUp, const unsigned int nParticles) const
{
	// Setup of result container
	std::vector<int> result;
	result.reserve(nParticles);
	double dice;

	// Find the list of probabilities
	std::list<std::pair<double, int>>::const_iterator cit;

	// Loop and insert PDG codes
	while(result.size() < nParticles)
	{
		dice = generator();
		// Search for fitting PDG code
		for(cit = particleLookUp.begin(); cit != particleLookUp.end(); cit++)
		{
			// Insert PGD code
			if(dice < cit->first)
			{
				result.push_back(cit->second);
				break;
			}
		}
	}
	
	return result;
}

template<typename generator_t, typename parameters_t>
std::vector<double>
ParametricNuclearInt::energyFractions(generator_t& generator, const parameters_t& params, const unsigned int nParticles) const
{
	// Storage of the resulting energies
	std::vector<double> result;
	result.resize(nParticles);

	double eFraction = 0.;
		
	// Extract the fit parameters
	const std::array<double, 10>& scalingFactors = params.energyScaling;

	// Repeat sampling as long as the energies are too big
	while(true)
	{
		double sumFractions = 0.;

		// Sample the energies from distribution and store them
		for(unsigned int n = 0; n < nParticles; n++)
		{
			if(nParticles <= 10)
			{
				eFraction = energyFraction(generator(), scalingFactors[n], n + 1);
			}
			else
			{
				// Extract the fit parameters to estimate the fit parameters for sampling
				const std::pair<double, double>& fitParameters = params.energyScalingExtrapolation;
				eFraction = energyFraction(generator(), fitParameters.first + (n + 1) * fitParameters.second, n + 1);
			}
			result[n] = eFraction;
			sumFractions += eFraction;
		}
		
		// Test if energies are <= the initial energy
		if(sumFractions <= 1.)
			break;
	}
	
	return result;
}

template<typename generator_t>
double
ParametricNuclearInt::sampleTheta(generator_t& generator, const std::array<double, 6>& fitParameters) const
{
	double cosTheta = generator();
	// Imprtance sampling until condition is fulfilled
	while(generator() > cosThetaProbability(cosTheta, fitParameters))
	{
		cosTheta = generator();
	}
	// Return theta
	return std::acos(cosTheta);
}

template<typename generator_t, typename particle_t, typename parameters_t>
std::vector<particle_t> 
ParametricNuclearInt::finalStateHadrons(generator_t& generator, const double thickness, particle_t& particle, const parameters_t& params) const
{
	// Calculate multiplicity
	const unsigned int Npart = multiplicity(generator, thickness, particle, params.multiplicity);

	// Easy exit if nothing gets out
	if(Npart == 0)
		return {};

	// Calculate particle types
	const std::vector<int> particlePDGs = particleComposition(generator, params.particleTypes, Npart);

	// Calculate the kinematics
	return kinematics(generator, particle, particlePDGs, params);
}
} // namespace Fatras
#include "Fatras/Physics/HadronicInteraction/detail/ParametricNuclearInt.ipp"