// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Fatras/Plugins/Geant4/Geant4MaterialInteraction.hpp"

#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4Material.hh"

namespace Fatras {
namespace Test{

class Geant4MaterialInteractionStub : public Geant4MaterialInteraction
{
public:
	template<typename particle_t>
	G4ParticleDefinition*
	convertParticleToG4Stub(const particle_t& particle) const
	{
		return convertParticleToG4(particle);
	}
	
	template<typename particle_t>
	G4ParticleGun*
	createParticleGunStub(const particle_t& particle) const
	{
		return createParticleGun(particle);
	}
	
	template<typename material_t>
	G4Material*
	convertMaterialToG4Stub(const material_t& material) const
	{
		return convertMaterialToG4(material);
	}
};
}
}
