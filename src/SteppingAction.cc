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
					G4AffineTransform momentumTransform = theTouchable->GetHistory()->GetTopTransform();
					momentumTransform.SetNetTranslation(G4ThreeVector(0,0,0));
					G4ThreeVector momentumDir = momentumTransform.TransformPoint(aStep->GetTrack()->GetMomentumDirection());

					if(std::fabs(localPos.x() - dimensionx) < kCarTolerance &&
					   momentumDir.x() > 0){
						fEvent->AddOutW();
					}
					else if(std::fabs(localPos.x() + dimensionx) < kCarTolerance &&
					   momentumDir.x() < 0){
						fEvent->AddOutW();
					}
					else if(std::fabs(localPos.y() - dimensiony) < kCarTolerance &&
					   momentumDir.y() > 0){
						fEvent->AddOut();
					}
					else if(std::fabs(localPos.y() + dimensiony) < kCarTolerance &&
					   momentumDir.y() < 0){
						fEvent->AddOut();
					}
					else if(std::fabs(localPos.z() - dimensionz) < kCarTolerance &&
					   momentumDir.z() > 0){
						fEvent->AddOutW();
					}
					else if(std::fabs(localPos.z() + dimensionz) < kCarTolerance &&
					   momentumDir.z() < 0){
						fEvent->AddOutW();
					}
				}
				if(aStep->GetTrack()->GetTrackStatus() != fAlive){
						fEvent->AddAbsW();
				}

			}
		}
	}
	else if(aStep->GetTrack()->GetTrackID() == 1) {
		if(aStep->GetTrack()->GetVolume()->GetName() != "Core"){
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

		G4StepPoint* preStep = aStep->GetPreStepPoint();
		G4TouchableHistory* thePreTouchable = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
		if(thePreTouchable->GetVolume()->GetName() != "World"){
			if(thePreTouchable->GetVolume(1)->GetName() == "Fiber"){
				fEvent->AddTrackLength(aStep->GetStepLength());
				if(aStep->IsFirstStepInVolume() && preStep->GetStepStatus() == fGeomBoundary){
					

					G4double kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
					G4double dimensionx = ((G4Box*) thePreTouchable->GetVolume(1)->GetLogicalVolume()->GetSolid())->GetXHalfLength();
					G4double dimensiony = ((G4Box*) thePreTouchable->GetVolume(1)->GetLogicalVolume()->GetSolid())->GetYHalfLength();
					G4double dimensionz = ((G4Box*) thePreTouchable->GetVolume(1)->GetLogicalVolume()->GetSolid())->GetZHalfLength();
					G4ThreeVector worldPos = aStep->GetPreStepPoint()->GetPosition();
					G4ThreeVector localPos = thePreTouchable->GetHistory()->GetTopTransform().TransformPoint(worldPos);
					G4AffineTransform momentumTransform = thePreTouchable->GetHistory()->GetTopTransform();
					momentumTransform.SetNetTranslation(G4ThreeVector(0,0,0));
					G4ThreeVector momentumDir = momentumTransform.TransformPoint(aStep->GetPreStepPoint()->GetMomentumDirection());


					/// Surfaces:
					///          - x>0: DS right. Surf = 0
					///          - x<0: DS left.  Surf = 1
					///          - y>0: DS up.    Surf = 2
					///          - y<0: DS down.  Surf = 3
					///          - z>0: DS front. Surf = 4
					///          - z<0: DS back.  Surf = 5
					if(std::fabs(localPos.x() - dimensionx) < kCarTolerance &&
					   momentumDir.x() < 0){
						fEvent->SetThetaIn(acos(-momentumDir.x()));
					}
					else if(std::fabs(localPos.x() + dimensionx) < kCarTolerance &&
					   momentumDir.x() > 0){
						fEvent->SetThetaIn(acos(momentumDir.x()));
					}
					else if(std::fabs(localPos.y() - dimensiony) < kCarTolerance &&
					   momentumDir.y() < 0){
						fEvent->SetThetaIn(acos(-momentumDir.y()));
					}
					else if(std::fabs(localPos.y() + dimensiony) < kCarTolerance &&
					   momentumDir.y() > 0){
						fEvent->SetThetaIn(acos(momentumDir.y()));
					}
					else if(std::fabs(localPos.z() - dimensionz) < kCarTolerance &&
					   momentumDir.z() < 0){
						fEvent->SetThetaIn(acos(-momentumDir.z()));
					}
					else if(std::fabs(localPos.z() + dimensionz) < kCarTolerance &&
					   momentumDir.z() > 0){
						fEvent->SetThetaIn(acos(momentumDir.z()));
					}
					fEvent->SetPrimaryChannel(thePreTouchable->GetReplicaNumber(2));
				}
			}
		}
	}
}
