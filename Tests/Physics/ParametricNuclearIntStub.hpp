// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Fatras/Physics/HadronicInteraction/ParametricNuclearInt.hpp"

namespace Fatras {
	
namespace Test {
	
struct ParametricNuclearIntStub : public ParametricNuclearInt {

double 
nuclearInteractionProbStub(const double momentum, const double thickness, const int pdg) const
{
	return nuclearInteractionProb(momentum, thickness, pdg);
}

double
multiplicityProbStub(const double momentum, const double thickness, const int pdg, const unsigned int mult) const
{
	return multiplicityProb(momentum, thickness, pdg, mult);
}

template<typename generator_t, typename particle_t>
unsigned int
multiplicityStub(generator_t& generator, const double thickness, particle_t& particle) const
{
	return multiplicity(generator, thickness, particle);
}

template<typename generator_t>
std::vector<int>
particleCompositionStub(generator_t& generator, const int pdg, const unsigned int nParticles) const
{
	return particleComposition(generator, pdg, nParticles);
}

template<typename generator_t, typename particle_t>
std::vector<particle_t>
kinematicsStub(generator_t& generator, particle_t& particle, const std::vector<int>& particlesPDGs) const
{
	return kinematics(generator, particle, particlesPDGs);
}

template<typename generator_t, typename particle_t>
std::vector<particle_t> 
finalStateHadronsStub(generator_t& generator, const double thicknessInL0, particle_t& particle) const
{
	return finalStateHadrons(generator, thicknessInL0, particle);
}
};
	
} // namespace Fatras
}
