// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Fatras/Physics/HadronicInteraction/ParametricNuclearInt.hpp"

double
Fatras::ParametricNuclearInt::nuclearInteractionProb(const double momentum, const double thickness, const int pdg) const
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

double
Fatras::ParametricNuclearInt::multiplicityProb(const double momentum, const double thickness, const int pdg, const unsigned int mult) const
{
	double scaleP, scaleD, scaling, exponent;
	// TODO: should be a switch
	//k0
	if(pdg == 311)
	{
		if(thickness < 0.2)
		{
			scaleP = -2.115 * exp(0.5228 * momentum);
			scaleD = 8.719 * exp(-1.623 * thickness);
			scaling = 43.01 + scaleD + scaleP;
			exponent = -0.5 * ((mult + 0.9625) / 0.3452 + exp(-(mult + 0.9625) / 0.3452));
		}
		else	
			if(thickness < 0.3)
			{
				scaleD = -6.897 * thickness;
				scaleP = 8.345 * momentum - 6.957 * momentum * momentum * 0.5;
				scaling = 104.4 + scaleP + scaleD;
				exponent = -0.5 * ((mult + 3.15) / 0.4815 + exp(-(mult + 3.15) / 0.4815));
			}
			else
			{   
				scaleP = 11.16 * exp(-4.377 * momentum);
				scaleD = -6.143 * exp(-4.123 * thickness);
				scaling = 14.76 + scaleD + scaleP;
				exponent = -0.5 * ((mult + 6.29) / 1.103 + exp(-(mult + 6.29) / 1.103));
			}
	}
	
	//k-
	if(pdg == -321)
	{
		if(momentum < 2.5 && thickness < 0.2)
		{
			scaleD = 0.1642 * thickness;
			scaleP = 13.34 * momentum - 8.865 * momentum * momentum * 0.5;
			scaling = 91.39 + scaleP + scaleD;
			exponent = -0.5 * ((mult + 21.48) / 1.909 + exp(-(mult + 21.48) / 1.909));
		}
		else
		{
			scaleP = -0.04855 * exp(-2.079  * momentum);
			scaleD = 0.1731 * exp(0.2436 * thickness);
			scaling = 0.06844 + scaleD + scaleP;
			exponent = -0.5 * ((mult - 3.357) / 1.598 + exp(-(mult - 3.357) / 1.598));
		}
	}
	
	//k+
	if(pdg == 321)
	{
		if(momentum < 1.5)
		{
			scaleD = -0.4838 * thickness;
			scaleP = 24.71 * momentum - 27.69 * momentum * momentum * 0.5;
			scaling = 35.8 + scaleP + scaleD;
			exponent = -0.5 * ((mult + 6.571) / 0.9084 + exp(-(mult + 6.571) / 0.9084));
		}
		else
		{
			scaleP = -0.6534 * exp(-2.729 * momentum);
			scaleD = 0.001245 * exp(2.788 * thickness);
			scaling = 0.3315 + scaleD + scaleP;
			exponent = -0.5 * ((mult - 2.752) / 1.182 + exp(-(mult - 2.752) / 1.182));
		}
	}
	
	//n
	if(pdg == 2112)
	{
		if(momentum < 1.5)
		{
			scaleD = -1.002 * thickness;
			scaleP = 12.5 * momentum - 14.79 * momentum * momentum * 0.5;
			scaling = 32.15 + scaleP + scaleD;
			exponent = -0.5 * ((mult + 5.3) / 0.826 + exp(-(mult + 5.3) / 0.826));
		}
		else
		{
			scaleP = -0.7808 * exp(-3.069 * momentum);
			scaleD = 0.009081 * exp(-19.4 * thickness);
			scaling = 0.3876 + scaleD + scaleP;
			exponent = -0.5 * ((mult - 2.466) / 1.04 + exp(-(mult - 2.466) / 1.04));
		}
	}

	//p 
	if(pdg == 2212)
	{
		if(momentum < 1.5)
		{
			scaleD = 1.317 * thickness;
			scaleP = 14.43 * momentum - 15.41 * momentum * momentum * 0.5;
			scaling = 15.29 + scaleP + scaleD;
			exponent = -0.5 * ((mult + 2.931) / 0.6564 + exp(-(mult + 2.931) / 0.6564));
		}
		else
		{
			scaleP = -2.223 * exp(-3.626 * momentum);
			scaleD = 0.0614 * exp(-57.61 * thickness);
			scaling = 0.5872 + scaleD + scaleP;
			exponent = -0.5 * ((mult - 1.113) / 1.032 + exp(-(mult - 1.113) / 1.032));
		}
	}
	
	//pi-
	if(pdg == -211)
	{
		if(momentum < 1.5)
		{
			scaleD = 0.207 * thickness;
			scaleP = 4.581 * momentum -3.31 * momentum * momentum * 0.5;
			scaling = 7.191 + scaleP + scaleD;
			exponent = -0.5 * ((mult + 5.983) / 1.1683 + exp(-(mult + 5.983) / 1.1683));
		}
		else
		{     

			scaleP = -1.094 * exp(-4.398 * momentum);
			scaleD = -0.02462 * exp(-4.425 * thickness);
			scaling = 0.3094 + scaleD + scaleP;
			exponent = -0.5 * ((mult - 2.316) / 1.473 + exp(-(mult - 2.316) / 1.473));
		}
	}

	//pi+
	if(pdg == 211)
	{
		if(momentum < 2.5)
		{
			scaleD = 0.207 * thickness;
			scaleP = 4.581 * momentum -3.31 * momentum * momentum * 0.5;
			scaling = 7.191 + scaleP + scaleD;
			exponent = -0.5 * ((mult + 5.983) / 1.1683 + exp(-(mult + 5.983) / 1.1683));
		}
		else
		{     
			scaleP = -0.786 * exp(-0.2199 * momentum);
			scaleD = -0.03931 * exp(-11.16 * thickness);
			scaling = 0.7287 + scaleD + scaleP;
			exponent = -0.5 * ((mult - 3.485) / 1.292 + exp(-(mult - 3.485) / 1.292));
		}
	}
	return exp(exponent) * scaling;
}
