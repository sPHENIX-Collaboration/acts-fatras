// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <map>
#include <array>
#include <list>

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

const std::map<int, std::list<std::pair<double, int>> particleTypes = 
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
	{111, {
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

	//~ // Look up table for produceable hadrons 
	//~ const std::array<double, 8> pdgLookUp = {-321, -211, 111, 211, 310, 321, 2112, 2212};
	
	//~ // Probabilities for the production of different particles
	//~ // Order in list: k-, pi-, pi0, pi+, k0, k+, n, p (as the PDG code)
	//~ const std::array<double, 8> probsKm = {0.131426, 0.108972, 2.15513e-05, 0.0735287, 0.00333459 + 0.0122508, 0.00437491, 0.586811, 0.072719};
	//~ const std::array<double, 8> probsKp = {0.00193334, 0.0857213, 2.62148e-05, 0.0949544, 0.00502494 + 0.0176528, 0.168827, 0.54008, 0.0818682};
	//~ const std::array<double, 8> probsK0 = {0.00880667, 0.111725, 3.91104e-05, 0.103377, 0.0368607 + 0.0782085, 0.0132087, 0.56815, 0.0746813};
	//~ const std::array<double, 8> probsPim = {0.00331695, 0.228171, 3.58508e-06, 0.0750272, 0.00127534 + 0.00512345, 0.00588292, 0.609138, 0.0678355};
	//~ const std::array<double, 8> probsPip = {0.00317051, 0.0880611, 1.59362e-06, 0.229114, 0.00128486 + 0.00513405, 0.00696432, 0.575764, 0.0862595};
	//~ const std::array<double, 8> probsP = {0.00102691, 0.0770944, 1.72117e-06, 0.0848894, 0.00051197 + 0.00272571, 0.00363871, 0.630106, 0.195689};
	//~ const std::array<double, 8> probsN = {0.00104396, 0.0940003, 1.66574e-06, 0.065843, 0.000550299 + 0.00277906, 0.00333905, 0.719939, 0.10825};

	//~ // Cumulative probabilities
	//~ const std::array<double, 8> cProbsKm = {0.131426, 0.240398, 0.24042, 0.313948, 0.329534, 0.333909, 0.92072, 0.993439};
	//~ const std::array<double, 8> cProbsKp = {0.00193334, 0.0876546, 0.0876809, 0.182635, 0.205313, 0.37414, 0.91422, 0.996088};
	//~ const std::array<double, 8> cProbsK0 = {0.00880667, 0.120532, 0.120571, 0.223948, 0.339017, 0.352226, 0.920376, 0.995057};
	//~ const std::array<double, 8> cProbsPim = {0.00331695, 0.231488, 0.231492, 0.306519, 0.312918, 0.3188, 0.927938, 0.995774};
	//~ const std::array<double, 8> cProbsPip = {0.00317051, 0.0912316, 0.0912332, 0.320347, 0.326766, 0.33373, 0.909494, 0.995754};
	//~ const std::array<double, 8> cProbsP = {0.00102691, 0.0781213, 0.078123, 0.163012, 0.16625, 0.169889, 0.799995, 0.995684};
	//~ const std::array<double, 8> cProbsN = {0.00104396, 0.0950443, 0.0950459, 0.160889, 0.164218, 0.167557, 0.887496, 0.995746};


} // namespace detail	
} // namespace Fatras
