// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <vector>

namespace Fatras {

class Geant4MaterialInteraction
{
public:

	Geant4MaterialInteraction() = default;
	
	template<typename particle_t, typename material_t>
	std::vector<particle_t> 
	operator()(const particle_t& particle, const material_t& material) const;
	
private:

	

};

template<typename particle_t, typename material_t>
std::vector<particle_t>
Geant4MaterialInteraction::operator()(const particle_t&, const material_t& material) const
{
	return {};
}

}
