// This file is part of the Acts project.
//
// Copyright (C) 2018-2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <map>

namespace Fatras {
namespace detail{

const double sqrt2pi = std::sqrt(2. * M_PI);

constexpr std::array<int, 5> pdgCodes = {-211, 111, 211, 2112, 2212};

/// Parameters used to estimate the probability for a nuclear interaction
const std::map<int, std::array<double, 6>> probability =
{

// pi-
	{-211, {-0.85589, 1.0763, -0.028606, 0.01827, 1.3097, 0.081749}},
// pi0
	{111, {0,0,0,0,0,0}},
// pi+
	{211, {-0.904434, 0.985991, -0.015039, 0.036966, 1.31977, 0.12179}},
// neutron
	{2112, {-1.04484, 0.67312, 0.079093, 0.42005, 1.8368, 0.92685}},
// proton
	{2212, {-1.01212, 0.717381, 0.075032, 0.35375, 1.89725, 0.83433}}
};

/// Parameters used to estimate the multiplicity in a nuclear interaction
const std::map<int, std::array<double, 7>> multiplicity =
{
// pi-
	{-211, {1.9677, -0.399229, -0.0405634, 0.915227, 1.39859, 0.130268, 0.0292009}},
// pi0
	{111, {0,0,0,0,0,0}},
// pi+
	{211, {1.22082, -0.661119, 0., 0.880236, 1.28554, 0.18008, 0.}},
// neutron
	{2112, {1.8136, -0.453892, 0., 0.900732, 1.187129, 0.125797, 0.}},
// proton
	{2212, {0.679744, -1.18508, 0.157405, 1.07033, 1.09336, -0.119505, 0.0505715}}
};

/// Cumulative probabilities for the production of resulting particle in a nuclear interaction
const std::map<int, std::list<std::pair<double, int>>> particleTypes = 
{
// pi-
	{-211, {
	std::make_pair(0.58345, -211),
	std::make_pair(0.585245, 130),
	std::make_pair(0.612815, 211),
	std::make_pair(0.614008, 321),
	std::make_pair(0.949433, 2112),
	std::make_pair(0.996385, 2212)
	}},
// pi0
	{111, {
	std::make_pair(0.0745813, -211),
	std::make_pair(0.129518, 211),
	std::make_pair(0.86305, 2112),
	std::make_pair(0.997921, 2212)
	}},
// pi+
	{211, {
	std::make_pair(0.037348, -211),
	std::make_pair(0.0384697, 130),
	std::make_pair(0.634316, 211),
	std::make_pair(0.636931, 321),
	std::make_pair(0.926136, 2112),
	std::make_pair(0.996833, 2212)
	}},
// neutron
	{2112, {
	std::make_pair(0.0381688, -211),
	std::make_pair(0.0516587, 211),
	std::make_pair(0.91314, 2112),
	std::make_pair(0.99883, 2212)
	
	}},
// proton
	{2212, {
	std::make_pair(0.0170427, -211),
	std::make_pair(0.0457174, 211),
	std::make_pair(0.378015, 2112),
	std::make_pair(0.998838, 2212)
	}}
};

/// Scaling factors of the energy distribution function
const std::map<int, std::array<double, 10>> energyScaling =
{
// pi-
{-211,
{1.43911, 3.03515, 6.24957, 13.4978, 35.7948, 53.0301, 63.4815, 72.3156, 80.5419, 88.7695}},
{111,
{0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}},
// pi+
{211,
{1.48089, 3.11388, 6.53058, 14.2392, 38.2195, 54.059, 63.3495, 71.2761, 78.8044, 86.3353}},
// neutron
{2112,
{0.984621, 2.5168, 5.44376, 12.6065, 41.0249, 58.18, 69.3694, 79.4628, 88.9836, 98.8031}},
// proton
{2212, 
{1.06923, 2.75259, 5.86034, 13.6034, 42.9559, 58.9314, 69.3068, 78.6077, 87.4014, 95.5143}}
};

} // namespace detail

template<typename generator_t>
std::vector<int>
ParametricNuclearInt::particleComposition(generator_t& generator, const int pdg, const unsigned int nParticles) const
{
	// Setup of result container
	std::vector<int> result;
	result.reserve(nParticles);
	double dice;

	// Find the list of probabilities
	const std::list<std::pair<double, int>>& particleLookUp = detail::particleTypes.at(pdg);
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
	
	// TODO: Don't know about the following lines
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

template<typename generator_t>
std::vector<double>
ParametricNuclearInt::energyFractions(generator_t& generator, const int pdg, const unsigned int nParticles) const
{
	if(nParticles < 10)
	{
	// Extract the fit parameters
	const std::array<double, 10>& scalingFactors = detail::energyScaling.at(pdg);
	
	// Storage of the resulting energies
	std::vector<double> result;
	result.resize(nParticles);

	// Repeat sampling as long as the energies are too big
	while(true)
	{
		double sumFractions = 0.;

		// Sample the energies from distribution and store them
		for(unsigned int n = 0; n < nParticles; n++)
		{
			double eFraction = energyFraction(generator(), scalingFactors[n], n + 1);
			result[n] = eFraction;
			sumFractions += eFraction;
		}
		
		// Test if energies are <= the initial energy
		if(sumFractions <= 1.)
			break;
	}
	
	return result;
}
else // TODO: implement linear functions
return {};
}

template <typename generator_t, typename detector_t, typename particle_t>
std::vector<particle_t> 
ParametricNuclearInt::operator()(generator_t& generator,
                                     const detector_t& detector,
                                     particle_t& particle) const
{
	// Test applicable PDG codes
	for(const auto& pc : detail::pdgCodes)
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
