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
	// TODO: move sqrt(2*PI) out

	const std::array<double, 6>& pars = detail::probability.at(pdg);
	//~ switch(pdg)
	//~ {
		//~ case 211:
		//~ {
			//~ pars = &detail::probPip;
			//~ break;
		//~ }
		//~ case -211:
		//~ {
			//~ pars = &detail::probPim;
			//~ break;
		//~ }
		//~ case 111:
		//~ {
			//~ pars = &detail::probPi0;
			//~ break;
		//~ }
		//~ case 2212:
		//~ {
			//~ pars = &detail::probP;
			//~ break;
		//~ }
		//~ case 2112:
		//~ {
			//~ pars = &detail::probN;
			//~ break;
		//~ }
		//~ default:
			//~ return 0.;
	//~ }
	
	const double shapeThickness = exp(thickness * pars[0]);
	
	const double shapeMomentum = pars[1] + pars[2] * momentum + pars[3] / (sqrt(2. * M_PI) * pars[5]) * exp(-(momentum - pars[4]) * (momentum - pars[4]) / (2. * pars[5]));
	
	return (1. - shapeThickness) * shapeMomentum;
	return 0.;
}

double
Fatras::ParametricNuclearInt::multiplicityProb(const double momentum, const double thickness, const int pdg, const unsigned int mult) const
{
	double scaleP, scaleD, scaling, exponent;
	
	// Multiplicity depends on the particle type
	switch(pdg)
	{
		//k0
		case 130:
		case 310:
		case 311:
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
				break;
		}
		
		//k-
		case -321:
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
			break;
		}
		
		//k+
		case 321:
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
			break;
		}
		
		//n
		case 2112:
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
			break;
		}

		//p
		case 2212:
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
			break;
		}
		
		//pi-
		case -211:
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
			break;
		}

		//pi+
		case 211:
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
			break;
		}
	}
	return exp(exponent) * scaling;
}



