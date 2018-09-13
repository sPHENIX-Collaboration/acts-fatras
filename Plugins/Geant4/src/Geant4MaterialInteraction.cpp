// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Fatras/Plugins/Geant4/Geant4MaterialInteraction.hpp"

namespace Fatras {

Geant4MaterialInteraction::Geant4MaterialInteraction()
{
	// Initialize Geant4 managers
	m_runManager = new G4RunManager;
	m_physicsList = new QBBC;
	m_runManager->SetUserInitialization(m_physicsList);
	m_particleTable = G4ParticleTable::GetParticleTable();
}

Geant4MaterialInteraction::~Geant4MaterialInteraction()
{
	// Free heap memory
	delete(m_physicsList);
}

std::pair<double, double>
Geant4MaterialInteraction::angleOfNormalVector(const Acts::Vector3D normalVector, const Acts::Vector3D momentum) const
{
	double tanphi = normalVector.y() / normalVector.x();
	if(std::isnan(tanphi))
		tanphi = 0.;
	else
		tanphi = std::atan(tanphi);
		
	// Check direction, invert direction if both vectors point in opposite direction
	if(normalVector.dot(momentum) >= 0.)
		return std::make_pair(std::acos(normalVector.z() / normalVector.norm()), tanphi);
	return std::make_pair(std::acos(normalVector.z() / normalVector.norm()) + M_PI, tanphi);		
}
}
