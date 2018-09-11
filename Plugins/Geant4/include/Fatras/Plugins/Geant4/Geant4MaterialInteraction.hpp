// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"

#include <vector>
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"

#include "G4Material.hh"

namespace Fatras {

class Geant4MaterialInteraction
{
public:

	Geant4MaterialInteraction() = default;
	
	template<typename particle_t, typename material_t>
	std::vector<particle_t> 
	operator()(const particle_t& particle, const material_t& material) const;
	
private:

	template<typename particle_t>
	G4ParticleDefinition*
	convertParticleToG4(const particle_t& particle) const;
	
	template<particle_t>
	G4ParticleGun*
	createParticleGun(G4ParticleDefinition* particleG4, const particle_t& particle) const;
	
	template<typename particle_t>
	particle_t
	convertParticleFromG4(const G4ParticleDefinition* particleG4) const;
	
	template<typename material_t>
	G4Material*
	convertMaterialToG4(const material& material) const;
	
	template<typename material_t>
	material_t
	convertMaterialFromG4(const G4Material* materialG4) const;
};

template<typename particle_t>
G4ParticleDefintion*
convertParticleToG4(const particle_t* particle) const
{
	if(particle->pdg() != 0)
	{
		G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable(); // TODO: will be reconstructed all the time, can be stored probably
		return particleTable->FindParticle(particle->pdg());
	}
	return nullptr;
}

template<typename particle_t>
G4ParticleGun*
createParticleGun(const particle_t& particle) const
{
	G4ParticleGun* pGun = new G4ParticleGun(1);
	
	G4ParticleDefinition* parDef = convertParticleToG4(particle);
	pGun->SetParticleDefinition(parDef);
	
	Acts::Vector3D momentum = particle.momentum();
	double scaleActsToG4 = MeV / Acts::units::_MeV 
	
	pGun->SetParticleMomentum({momentum.x() * scaleActsToG4, momentum.y() * scaleActsToG4, momentum.z() * scaleActsToG4});
	pGun->SetParticlePosition({0., 0., 0.,});
	pGun->SetParticleTime(0.); // TODO: passed path in L0,X0 and time missing
	return pGun;
}

  
template<typename particle_t, typename material_t>
std::vector<particle_t>
Geant4MaterialInteraction::operator()(const particle_t&, const material_t& material) const
{
	G4ParticleGun* pGun = createParticleGun(particle);
	
	return {};
}

}
