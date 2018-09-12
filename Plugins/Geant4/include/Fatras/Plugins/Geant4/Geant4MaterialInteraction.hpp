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
#include "Fatras/Plugins/Geant4/B1ActionInitialization.hpp"
#include "Fatras/Plugins/Geant4/B1DetectorConstruction.hpp"

#include <vector>
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
#include "QBBC.hh"
#include "G4Material.hh"

namespace Fatras {

class Geant4MaterialInteraction
{
public:

	Geant4MaterialInteraction();
	
	~Geant4MaterialInteraction();
	
	template<typename particle_t, typename material_t>
	std::vector<particle_t> 
	operator()(particle_t& particle, const material_t& material, const double thickness) const;
	
protected:

	template<typename particle_t>
	G4ParticleDefinition*
	convertParticleToG4(const particle_t& particle) const;
	
	template<typename particle_t>
	G4ParticleGun*
	createParticleGun(const particle_t& particle) const;
	
	template<typename material_t>
	G4Material*
	convertMaterialToG4(const material_t& material) const;
	
	template<typename particle_t>
	void
	convertParticlesFromG4(const std::vector<B1particle>& particlesG4, particle_t& particleIn, std::vector<particle_t>& particles) const;
	
private:

	G4RunManager* runManager;
	G4VModularPhysicsList* physicsList;
	G4ParticleTable* particleTable;
};

Geant4MaterialInteraction::Geant4MaterialInteraction()
{
	runManager = new G4RunManager;
	physicsList = new QBBC;
	runManager->SetUserInitialization(physicsList);
	particleTable = G4ParticleTable::GetParticleTable();
}

Geant4MaterialInteraction::~Geant4MaterialInteraction()
{
	delete(physicsList);	
	//~ delete(runManager);
}

template<typename particle_t>
G4ParticleDefinition*
Geant4MaterialInteraction::convertParticleToG4(const particle_t& particle) const
{
	if(particle.pdg() != 0)
		return particleTable->FindParticle(particle.pdg());
	return nullptr;
}

template<typename particle_t>
G4ParticleGun*
Geant4MaterialInteraction::createParticleGun(const particle_t& particle) const
{
	G4ParticleGun* pGun = new G4ParticleGun(1);
	
	G4ParticleDefinition* parDef = convertParticleToG4(particle);
	pGun->SetParticleDefinition(parDef);
	
	Acts::Vector3D momentum = particle.momentum();
	double scaleActsToG4 = MeV / Acts::units::_MeV;

	//TODO: Angles of incoming particles
	pGun->SetParticleMomentum({momentum.x() * scaleActsToG4, momentum.y() * scaleActsToG4, momentum.z() * scaleActsToG4});
	pGun->SetParticlePosition({0., 0., 0.,});
	pGun->SetParticleTime(0.); // TODO: passed path in L0,X0 and time missing
	return pGun;
}

template<typename material_t>
G4Material*
Geant4MaterialInteraction::convertMaterialToG4(const material_t& material) const
{
	return new G4Material("Material", material.Z(), material.A() * g / mole, material.rho() * Acts::units::_cm * Acts::units::_cm * Acts::units::_cm / Acts::units::_g * g / cm3);
}

template<typename particle_t>
void
Geant4MaterialInteraction::convertParticlesFromG4(const std::vector<B1particle>& particlesG4, particle_t& particleIn, std::vector<particle_t>& particles) const
{
	Acts::Vector3D momentum;
	// TODO: Angles of ougoing particles
	for(const B1particle& bp : particlesG4)
	{
		momentum = {bp.momentum[0], bp.momentum[1], bp.momentum[2]};
		particle_t p(particleIn.position(), momentum, bp.mass, bp.charge, bp.pdg);
		particles.push_back(std::move(p));
	}
}

template<typename particle_t, typename material_t>
std::vector<particle_t>
Geant4MaterialInteraction::operator()(particle_t& particle, const material_t& material, const double thickness) const
{	
	double materialThickness = thickness * mm / Acts::units::_mm;
	
	G4ParticleGun* pGun = createParticleGun(particle);
	B1ActionInitialization* actionInit = new B1ActionInitialization(materialThickness, pGun);
	
	G4Material* materialG4 = convertMaterialToG4(material);
	B1DetectorConstruction* detConstr = new B1DetectorConstruction(materialG4, materialThickness);
	
	runManager->SetUserInitialization(detConstr);
	runManager->SetUserInitialization(actionInit);
	runManager->Initialize();
	runManager->BeamOn(1);
	
	std::vector<particle_t> particles;
	convertParticlesFromG4(actionInit->particles(), particle, particles);
	
	delete(actionInit);
	delete(pGun);
	delete(materialG4);
	delete(detConstr);	
	
	return particles;
}

}
