/// \file  SteppingAction.cc
/// \brief Implementation of the SteppingAction class

#include "SteppingAction.hh"
#include "EventAction.hh"

#include "G4Track.hh"
#include "G4Step.hh"
#include "G4Event.hh"
#include "G4Box.hh"
#include "G4OpticalPhoton.hh"
#include "G4TouchableHistory.hh"
#include "G4GeometryTolerance.hh"

SteppingAction::SteppingAction(EventAction* ea):fEvent(ea){}

SteppingAction::~SteppingAction() {}

void SteppingAction::UserSteppingAction(const G4Step * aStep){
	if(aStep->GetTrack()->GetParticleDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()){
		G4TouchableHistory* theTouchable = 
			(G4TouchableHistory*) (aStep->GetPreStepPoint()->GetTouchable());
		if(theTouchable->GetVolume()->GetName() != "World"){
			if(theTouchable->GetVolume(1)->GetName() == "Fiber"){
				
				if(aStep->GetPostStepPoint()->GetStepStatus() == fGeomBoundary){
					G4double kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
					G4double dimensionx = ((G4Box*) theTouchable->GetVolume(1)->GetLogicalVolume()->GetSolid())->GetXHalfLength();
					G4double dimensiony = ((G4Box*) theTouchable->GetVolume(1)->GetLogicalVolume()->GetSolid())->GetYHalfLength();
					G4double dimensionz = ((G4Box*) theTouchable->GetVolume(1)->GetLogicalVolume()->GetSolid())->GetZHalfLength();
					G4ThreeVector worldPos = aStep->GetPostStepPoint()->GetPosition();
					G4ThreeVector localPos = theTouchable->GetHistory()->GetTransform(1).TransformPoint(worldPos);

					if(std::fabs(localPos.x() - dimensionx) < kCarTolerance &&
					   aStep->GetTrack()->GetMomentumDirection().x() > 0){
						fEvent->AddOutW();
					}
					else if(std::fabs(localPos.x() + dimensionx) < kCarTolerance &&
					   aStep->GetTrack()->GetMomentumDirection().x() < 0){
						fEvent->AddOutW();
					}
					else if(std::fabs(localPos.y() - dimensiony) < kCarTolerance &&
					   aStep->GetTrack()->GetMomentumDirection().y() > 0){
						fEvent->AddOut();
					}
					else if(std::fabs(localPos.y() + dimensiony) < kCarTolerance &&
					   aStep->GetTrack()->GetMomentumDirection().y() < 0){
						fEvent->AddOut();
					}
					else if(std::fabs(localPos.z() - dimensionz) < kCarTolerance &&
					   aStep->GetTrack()->GetMomentumDirection().z() > 0){
						fEvent->AddOutW();
					}
					else if(std::fabs(localPos.z() + dimensionz) < kCarTolerance &&
					   aStep->GetTrack()->GetMomentumDirection().z() < 0){
						fEvent->AddOutW();
					}
				}
				if(aStep->GetTrack()->GetTrackStatus() != fAlive){
						fEvent->AddAbsW();
				}

			}
		}
	}
	else if(aStep->GetTrack()->GetVolume()->GetName() != "Core") {
		const std::vector<const G4Track*>* secondaries = aStep->GetSecondaryInCurrentStep();
		if(secondaries->size() > 0){
			for(unsigned int i = 0; i < secondaries->size(); i++){
				if(secondaries->at(i)->GetParentID() > 0){
					if(secondaries->at(i)->GetDynamicParticle()->GetParticleDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()){
						fEvent->AddGamma();
					}
				}
			}
		}
	}
}
