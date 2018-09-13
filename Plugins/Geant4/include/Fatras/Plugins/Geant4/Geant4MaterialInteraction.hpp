// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"
#include "Fatras/Plugins/Geant4/MIActionInitialization.hpp"
#include "Fatras/Plugins/Geant4/MIDetectorConstruction.hpp"

#include <vector>
#include <cmath>
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
#include "QBBC.hh"
#include "G4Material.hh"

namespace Fatras {

/// @brief This class is a converter for material interaction of particles.
/// Given an Acts described particle and material, this class converts the
/// environment into terms of Geant4 and calculates the interaction.
/// The result of the interaction is transformed back into terms of Acts
class Geant4MaterialInteraction
{
public:

	/// @brief Constructor
	Geant4MaterialInteraction();
	
	/// @brief Destructor
	~Geant4MaterialInteraction();
	
	/// @brief Call parameter for the material interaction
	///
	/// @tparam particle_t Type of the particle
	/// @tparam material_t Type of the material
	/// @param [in] particle Ingoing particle
	/// @param [in] material Penentrated material
	/// @param [in] thickness Thickness of the material
	/// @param [in] normalVector Normal vector of the penetrated material
	/// @return Vector containing all outgoing particles
	template<typename particle_t, typename material_t>
	std::vector<particle_t> 
	operator()(particle_t& particle, const material_t& material, const double thickness, const Acts::Vector3D normalVector) const;
	
protected:

	/// @brief Calculates the angle of the normal vector for coordinate transformation
	///
	/// @param [in] normalVector Normal vector of the penetrated material
	/// @param [in] momentum Momentum vector of the particle
	/// @return Pair consisting of Theta and Phi in spherical coordinates for the transformation
	std::pair<double, double>
	angleOfNormalVector(const Acts::Vector3D normalVector, const Acts::Vector3D momentum) const;

	/// @brief Converts a particle into a Geant4 particle
	///
	/// @tparam particle_t Type of the particle
	/// @param [in] particle Particle that is converted
	/// @return Pointer to the corresponding Geant4 particle
	template<typename particle_t>
	G4ParticleDefinition*
	convertParticleToG4(const particle_t& particle) const;
	
	/// @brief Constructs a Geant4 particle gun
	///
	/// @tparam particle_t Type of the particle
	/// @param [in] particle Ammo of the gun
	/// @param [in] angles Angles of the penetrated material for transformation
	/// @return Pointer to the particle gun
	template<typename particle_t>
	G4ParticleGun*
	createParticleGun(const particle_t& particle, const std::pair<double, double>& angles) const;
	
	/// @brief Converts material into Geant4 material
	///
	/// @tparam material_t Type of the material
	/// @param [in] material Material that is penetrated
	/// @return Pointer to the corresponding Geant4 material
	template<typename material_t>
	G4Material*
	convertMaterialToG4(const material_t& material) const;
	
	/// @brief Converts Geant4 particles back
	///
	/// @tparam particle_t Type of the particle
	/// @param [in] particlesG4 Vector of recorded outgoing particles
	/// @param [in] particleIn Ingoing particle
	/// @param [in] angles Angles of the penentrated material for transformation
	/// @param [out] particles Vector of outgoing particles
	template<typename particle_t>
	void
	convertParticlesFromG4(const std::vector<MIparticle>& particlesG4, const particle_t& particleIn, const std::pair<double, double>& angles, std::vector<particle_t>& particles) const;
	
private:

	// Geant4 run manager
	G4RunManager* m_runManager;
	// List of physics processes in Geant4
	G4VModularPhysicsList* m_physicsList;
	// Table of particles in Geant4
	G4ParticleTable* m_particleTable;
};

template<typename particle_t>
G4ParticleDefinition*
Geant4MaterialInteraction::convertParticleToG4(const particle_t& particle) const
{
	// Check if pdg code is provided, return related particle or nothing
	if(particle.pdg() != 0 && std::isfinite(particle.pdg()))
	{
		G4ParticleDefinition* parDef = m_particleTable->FindParticle(particle.pdg());
		parDef->SetPDGLifeTime(parDef->GetPDGLifeTime() - particle.t());
		return parDef;
	}
	return nullptr;
}

template<typename particle_t>
G4ParticleGun*
Geant4MaterialInteraction::createParticleGun(const particle_t& particle, const std::pair<double, double>& angles) const
{	
	// Create particle
	G4ParticleDefinition* parDef = convertParticleToG4(particle);

	if(parDef)
	{
		// Build gun
		G4ParticleGun* pGun = new G4ParticleGun(1);
		pGun->SetParticleDefinition(parDef);
		
		// Set initial kinematics
		Acts::Vector3D momentum = particle.momentum();
		
		// Correct angles
		double theta = std::acos(momentum.z() / momentum.norm()) - angles.first;;
			
		double phi;
		if(momentum.x() == 0. && momentum.y() == 0.)
			phi = -angles.second;
		else
			phi = std::atan(momentum.y() / momentum.x()) - angles.second;
			
		double scaleActsToG4 = MeV / Acts::units::_MeV;

		pGun->SetParticleMomentum({particle.p() * std::sin(theta) * std::cos(phi) * scaleActsToG4, 
								   particle.p() * std::sin(theta) * std::sin(phi) * scaleActsToG4, 
								   particle.p() * std::cos(theta) * scaleActsToG4});
								   
		pGun->SetParticlePosition({0., 0., 0.,});
		pGun->SetParticleTime(0.);
		return pGun;
	}
	return nullptr;
}

template<typename material_t>
G4Material*
Geant4MaterialInteraction::convertMaterialToG4(const material_t& material) const
{
	// Translate material if valid
	if(material.Z() < 0. || material.A() < 0. || material.rho() < 0.
		|| !std::isfinite(material.Z()) || !std::isfinite(material.A()) || !std::isfinite(material.rho()))
		return nullptr;
	return new G4Material("", material.Z(), material.A() * g / mole, material.rho() * Acts::units::_cm * Acts::units::_cm * Acts::units::_cm / Acts::units::_g * g / cm3);
}

template<typename particle_t>
void
Geant4MaterialInteraction::convertParticlesFromG4(const std::vector<MIparticle>& particlesG4, const particle_t& particleIn, const std::pair<double, double>& angles, std::vector<particle_t>& particles) const
{
	Acts::Vector3D momentum;
	const double scaleG4ToActs = Acts::units::_GeV / GeV;
	double p, theta, phi;

	// Translate every particle from Geant4
	for(const MIparticle& bp : particlesG4)
	{
		// Correct angles
		p = bp.momentum.norm();
		theta = std::acos(bp.momentum.z() / p) - angles.first; // Modified due to symmetry of cos
		if(bp.momentum.x() == 0. && bp.momentum.y() == 0.)
			phi = angles.second;
		else
			phi = std::atan(bp.momentum.y() / bp.momentum.x()) + angles.second;

		momentum = {p * std::sin(theta) * std::cos(phi) * scaleG4ToActs, 
				    p * std::sin(theta) * std::sin(phi) * scaleG4ToActs,
				    p * std::cos(theta) * scaleG4ToActs};
				    
		particle_t part(particleIn.position(), momentum, bp.mass * Acts::units::_GeV / GeV, bp.charge, bp.pdg);
		particles.push_back(std::move(part));
	}
}

template<typename particle_t, typename material_t>
std::vector<particle_t>
Geant4MaterialInteraction::operator()(particle_t& particle, const material_t& material, const double thickness, const Acts::Vector3D normalVector) const
{	
	// Protection against invalid normal vectors
	if(normalVector == Acts::Vector3D::Zero(3))
		return {};
		
	double materialThickness = thickness * mm / Acts::units::_mm;
	std::pair<double, double> angles = angleOfNormalVector(normalVector, particle.momentum());

	// Load the gun and build material
	G4ParticleGun* pGun = createParticleGun(particle, angles);
	G4Material* materialG4 = convertMaterialToG4(material);
	
	if(pGun && materialG4 && materialThickness > 0.)
	{
		// Configure the Process
		MIActionInitialization* actionInit = new MIActionInitialization(materialThickness, pGun);			
		MIDetectorConstruction* detConstr = new MIDetectorConstruction(materialG4, materialThickness);

		m_runManager->SetUserInitialization(detConstr);
		m_runManager->DefineWorldVolume(detConstr->Construct());

		m_runManager->SetUserInitialization(actionInit);
		m_runManager->Initialize();
		m_runManager->BeamOn(1);
	
		// Collect the result
		std::vector<particle_t> particles;
		convertParticlesFromG4(actionInit->particles(), particle, angles, particles);

		// Free memory
		delete(detConstr);
		delete(actionInit);
		
		return particles;
	}
	return {};
}

}
