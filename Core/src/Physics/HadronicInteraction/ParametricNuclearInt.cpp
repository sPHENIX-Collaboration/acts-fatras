// This file is part of the Acts project.
//
// Copyright (C) 2018-2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Fatras/Physics/HadronicInteraction/ParametricNuclearInt.hpp"
#include "Fatras/Physics/HadronicInteraction/detail/ParametricNuclearInt.ipp"

double
Fatras::ParametricNuclearInt::nuclearInteractionProb(const double momentum, const double thickness, const std::array<double, 6>& pars) const
{
	const double shapeThickness = exp(thickness * pars[0]);
	
	const double shapeMomentum = pars[1] + pars[2] * momentum + pars[3] / (detail::sqrt2pi * pars[5]) * exp(-(momentum - pars[4]) * (momentum - pars[4]) / (2. * pars[5]));
	
	return (1. - shapeThickness) * shapeMomentum;
}

double
Fatras::ParametricNuclearInt::multiplicityProb(const double momentum, const double thickness, const std::array<double, 7>& pars, const unsigned int mult) const
{	
	const double diff = (mult - pars[0] + pars[1] * momentum + pars[2] * momentum * momentum) / (pars[3] + pars[5] * momentum + pars[6] * momentum * momentum);
	
	const double landau = exp(-0.5 * (diff + exp(-diff)));
		
	return landau / (pars[4] * detail::sqrt2pi);
}

double 
Fatras::ParametricNuclearInt::energyFraction(const double cProb, const double scaling, const unsigned int n) const
{
	const double ln = std::log(1 - cProb);
	
	return -ln / (scaling * n - ln);
}

double
Fatras::ParametricNuclearInt::cosThetaProbability(double cosTheta, const std::array<double, 6>& fitParameters) const
{
	return fitParameters[0] + cosTheta * fitParameters[1] + cosTheta * cosTheta * fitParameters[2] + cosTheta * cosTheta * cosTheta * fitParameters[3] + 
			cosTheta * cosTheta * cosTheta * cosTheta * fitParameters[4] + cosTheta * cosTheta * cosTheta * cosTheta * cosTheta * fitParameters[5];
}