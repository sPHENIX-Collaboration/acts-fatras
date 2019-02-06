// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <array>

namespace Fatras {
namespace detail{

/// Parameters used to estimate the probability for a nuclear interaction
// pi+
	constexpr std::array<double, 6> probPip = {-0.904434, 0.985991, -0.015039, 0.036966, 1.31977, 0.12179};
// pi-
	constexpr std::array<double, 6> probPim = {-0.85589, 1.0763, -0.028606, 0.01827, 1.3097, 0.081749};
// pi0
	constexpr std::array<double, 6> probPi0 = {0,0,0,0,0,0};
// proton
	constexpr std::array<double, 6> probP = {-1.01212, 0.717381, 0.075032, 0.35375, 1.89725, 0.83433};
// neutron
	constexpr std::array<double, 6> probN = {-1.04484, 0.67312, 0.079093, 0.42005, 1.8368, 0.92685};
	
} // namespace detail	
} // namespace Fatras
