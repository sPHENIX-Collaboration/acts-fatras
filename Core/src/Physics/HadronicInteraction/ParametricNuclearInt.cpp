// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Fatras/Physics/HadronicInteraction/ParametricNuclearInt.hpp"

double
Fatras::ParametricNuclearInt::nuclearInteractionProb(const double thickness, const double momentum, const int pdg) const
{
	//p/n
	if(pdg == 2112 || pdg == 2212)
	{	
		if(momentum > 3.)
			return (1. - exp(-thickness * 0.98835 * (0.07491 + exp(-(momentum - 13.744) * (momentum - 13.744) / 805.25)))) * 0.97985;
		else
			if(thickness > 0.5)
				return 0.05703 + (1. - exp(-thickness * 1.1889 * (1.1714e-06 + exp(-(momentum - 2.497) * (momentum - 2.497) / 2.0915)))) * 0.75229;
			else
				return 0.009525 + (1. - exp(-thickness * 1.4277 * (0.011844 + exp(-(momentum - 2.3332) * (momentum - 2.3332) / 3.9307)))) * 0.71596;
	}

	//pi+/pi-
	if(pdg == 211 || pdg == -211)
	{
		if(momentum > 4.)
			return (1. - exp(-thickness * 0.55001 * (0.46996 + exp(-(momentum - 8.5732) * (momentum - 8.5732) / 1800.)))) * 0.96459;
		else 
			if(thickness > 0.5)
				return 0.03192 + (1. - exp(-thickness * 1.0276 * (0.028738 + exp(-(momentum - 3.1018) * (momentum - 3.1018) / 8.866)))) * 0.76954;
			else
				return 0.003446 + (1. - exp(-thickness * 1.0482 * (0.074884 + exp(-(momentum - 2.277) * (momentum - 2.277) / 22.708)))) * 0.74133;
	}
	
	//kaon+
	if(pdg == 321)
	{
		if(momentum > 4.)
			return (1. - exp(-thickness * 0.13597 * (4.8514 + exp(-(momentum - 8.1613e-05) * (momentum - 8.1613e-05) / 771.03)))) * 0.88989;
		else 
			if(thickness > 0.5)
				return 0.009418 + (1. - exp(-thickness * 0.83191 * (7.8499e-10 + exp(-(momentum - 2.9685) * (momentum - 2.9685) / 4.0707)))) * 0.74356;
			else
				return -0.004841 + (1. - exp(-thickness * 1.5404 * (0.00046885 + exp(-(momentum - 2.8789) * (momentum - 2.8789) / 8.8744)))) * 0.44758;
	}
	
	//kaon-
	if(pdg == -321)
	{
		if(momentum > 4.)
			return (1. - exp(-thickness * 0.13597 * (4.8514 + exp(-(momentum - 8.1613e-05) * (momentum - 8.1613e-05) / 771.03)))) * 0.88989;
		else 
			if(thickness > 0.5)
				return 0.0004233 + (1. - exp(-thickness * 0.70435 * (0.40021 + exp(-(momentum - 5.8872e-05) * (momentum - 5.8872e-05) / 107.32)))) * 0.82217;
			else
				return 0.003083 + (1. - exp(-thickness * 3.6515 * (0.14931 + exp(-(momentum - 88.407) * (momentum - 88.407) / 15.98)))) * 1.3231;
	}
			
	//kaon0
	if(pdg == 311)
	{
		if(momentum > 4.)
			if(thickness < 0.5 && momentum < 12.)
				return (1. - exp(-thickness * 0.62966 * (1.8493 + exp(-(momentum - 8.9774e-08) * (momentum - 8.9774e-08) / 44.93)))) * 0.72643;
			else
				return (1. - exp(-thickness * 0.13597 * (4.8514 + exp(-(momentum - 8.1613e-05) * (momentum - 8.1613e-05) / 771.03)))) * 0.88989;
		else 
			if(thickness > 0.5)
				return 0.2789 + (1. - exp(-thickness * 0.98166 * (6.3897e-08 + exp(-(momentum - 3.1151) * (momentum - 3.1151) / 4.08)))) * 0.46121;
			else
				return 0.01415 + (1. - exp(-thickness * 3.8071 * (0.59535 + exp(-(momentum - 0.0004938) * (momentum - 0.0004938) / 10.21)))) * 0.42625;
	}

	return 0.;
}
