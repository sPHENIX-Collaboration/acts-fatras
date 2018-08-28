//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: B1SteppingAction.cc 74483 2013-10-09 13:37:06Z gcosmo $
//
/// \file B1SteppingAction.cc
/// \brief Implementation of the B1SteppingAction class

#include "../include/B1SteppingAction.hh"
#include "../include/B1EventAction.hh"
#include "../include/B1DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"

#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1SteppingAction::B1SteppingAction(B1EventAction* eventAction, double thickness)
: G4UserSteppingAction(),
  fEventAction(eventAction),
  fScoringVolume(0),
  m_thickness(thickness)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1SteppingAction::~B1SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1SteppingAction::UserSteppingAction(const G4Step* step)
{
//~ std::cout << "dicki hoppenstett: " << m_thickness * cm << std::endl;
	if(fabs((step->GetPostStepPoint()->GetPosition().z() - m_thickness * cm) / (m_thickness * cm)) < 1e-15
		|| step->GetPostStepPoint()->GetPosition().z() < 0. * cm
		|| fabs((step->GetPostStepPoint()->GetPosition().x() - 10. * cm) / (10. * cm)) < 1e-15
		|| fabs((step->GetPostStepPoint()->GetPosition().y() - 10. * cm) / (10. * cm)) < 1e-15)
	{
		B1particle p;
		p.position[0] = step->GetPostStepPoint()->GetPosition().x();
		p.position[1] = step->GetPostStepPoint()->GetPosition().y();
		p.position[2] = step->GetPostStepPoint()->GetPosition().z();
		p.momentum[0] = step->GetPostStepPoint()->GetMomentum().x();
		p.momentum[1] = step->GetPostStepPoint()->GetMomentum().y();
		p.momentum[2] = step->GetPostStepPoint()->GetMomentum().z();
		p.pdg = step->GetTrack()->GetDynamicParticle()->GetPDGcode();
		p.energy = step->GetPostStepPoint()->GetTotalEnergy();
		p.mass = step->GetPostStepPoint()->GetMass();
		p.charge = step->GetPostStepPoint()->GetCharge();
		p.trackid = step->GetTrack()->GetTrackID();
		p.parentid = step->GetTrack()->GetParentID();
		fEventAction->AddParticle(p);
	}
	
  if (!fScoringVolume) { 
    const B1DetectorConstruction* detectorConstruction
      = static_cast<const B1DetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    fScoringVolume = detectorConstruction->GetScoringVolume();   
  }

  // get volume of the current step
  G4LogicalVolume* volume 
    = step->GetPreStepPoint()->GetTouchableHandle()
      ->GetVolume()->GetLogicalVolume();

  // check if we are in scoring volume
  if (volume != fScoringVolume) return;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

