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

namespace Fatras {

/// The struct for the physics list
///
///
struct ParametricNuclearInt {

struct Config()
{
	bool m_hadronInteractionFromX0;
	double m_hadronInteractionProbabilityScale;
}

template <typename material_t, typename particle_t>
double absorptionLength(const material_t* matertial, double p, particle_t& particle) const 
{
  double al = material->l0();

  if(particle.pdg == 211 || particle.pdg == -211 || particle.pdg == 321 || particle.pdg == 111 || particle.pdg == 311) // TODO: nur K^+?
    al *= 1. / (1. + exp(-0.5 * (p - 270.) * (p - 270.) / 3600.)); // TODO: da kann man sicherlich noch etwas optimieren

  if(particle.pdg == 2212 || particle.pdg == 2212) al *= 0.7;
  if(particle.pdg == 211 || particle.pdg == -211 || particle.pdg == 111) al *= 0.9;

  return al;
}

template<typename particle_t>
std::vector<particle_t> recordHadronState(double time, particle_t& particle ) const {
	
	std::vector<particle_t> particleVector = getHadronState(time, particle);
	//~ ISF::ISFParticleVector ispVec=getHadState(parent, time, p, vertex, particleDir, particle);
	  
	// having no secondaries does not necessarily mean the interaction did not take place : TODO : add flag into ::getHadState
	//  if (!ispVec.size()) return false;
	  
	// push onto ParticleStack
	if (particleVector.size()) {
		for (unsigned int ic = 0; ic < particleVector.size(); ic++) 
		{
			if(!particleVector[ic]->getTruthBinding()) { // TODO: wtf is das? - binding von fatras zu mc truth
				particleVector[ic]->setTruthBinding(new ISF::TruthBinding(*parent->getTruthBinding()));
		    }
	    }  
	}
	return particleVector;
}

template <typename generator_t, typename material_t, typename particle_t>
std::vector<particle_t> hadronicInteraction(generator_t& generator, const material_t& material, double pathCorrection, particle_t& particle) const
{
	const material_t* extMprop = dynamic_cast<const material_t*>(&material);
	double prob = 0.;

	// m_hadIntProbScale is used later, not here
	if (extMprop && !m_cfg.m_hadronInteractionFromX0) {
		  
	double al = absorptionLength(extMprop, p, particle);  // in mm

    if (al > 0.) prob = exp(-pathCorrection * extMprop->thickness() / al);
    else       prob = exp(-pathCorrection * extMprop->thicknessInL0()); 

	} else {
		// using approximation lambda = 0.37 * Z * X0 instead -- giving a warning message
		prob = exp(-pathCorrection * mprop.thicknessInX0() / (0.37 * mprop.averageZ()));
	}
  
	// apply a global scalor of the probability
	// (1. - prob) is generally O(0.01), so this is the right way to scale it
	// TODO fix time info (if needed)
	if (generator() < (1. - prob) * m_hadronInteractionProbabilityScale) 
		return recordHadronState(0., particle);
  
	// no hadronic interactions were computed
	return {};  
}

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
                                     particle_t &particle) const {

    // todo return photons
    return {};
  }
};

} // namespace Fatras
