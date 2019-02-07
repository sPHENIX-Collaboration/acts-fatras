// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Fatras/Physics/HadronicInteraction/ParametricNuclearInt.hpp"
#include "Fatras/Physics/HadronicInteraction/detail/ParametricNuclearInt.ipp"

double
Fatras::ParametricNuclearInt::nuclearInteractionProb(const double momentum, const double thickness, const int pdg) const
{
	const std::array<double, 6>& pars = detail::probability.at(pdg);

	const double shapeThickness = exp(thickness * pars[0]);
	
	const double shapeMomentum = pars[1] + pars[2] * momentum + pars[3] / (detail::sqrt2pi * pars[5]) * exp(-(momentum - pars[4]) * (momentum - pars[4]) / (2. * pars[5]));
	
	return (1. - shapeThickness) * shapeMomentum;
	return 0.;
}

double
Fatras::ParametricNuclearInt::multiplicityProb(const double momentum, const double thickness, const int pdg, const unsigned int mult) const
{
	const std::array<double, 7>& pars = detail::multiplicity.at(pdg);
	
	const double diff = (i - pars[0] + pars[1] * momentum + pars[2] * momentum * momentum) / (pars[3] + pars[5] * momentum + pars[6] * momentum * momentum);
	const double landau = exp(-0.5 * (diff + exp(-diff)));
		
	return landau / (pars[4] * detail::sqrt2pi);
}



