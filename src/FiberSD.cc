/// \file  ScinSD.cc
/// \brief Implementation of the FiberSD class

#include "FiberSD.hh"
#include "FiberHit.hh"
#include "RunAction.hh"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4VTouchable.hh"
#include "G4ios.hh"
#include "G4VProcess.hh"
#include "G4RunManager.hh"

#include "G4Box.hh"
#include "G4GeometryTolerance.hh"
#include "G4VSolid.hh"


#include "TVector3.h"

FiberSD::FiberSD(G4String name) : G4VSensitiveDetector(name), fEin(0), fEdep(0), fDelta(0), fNgamma(0), fNgammaOut(0){
	fFiberCollection = nullptr;
	collectionName.insert("fiberCollection");
}

FiberSD::~FiberSD(){}

void FiberSD::Initialize(G4HCofThisEvent* hitsCE){
	fFiberCollection = new FiberHitsCollection(SensitiveDetectorName, collectionName[0]);

	// Putting all the hits in the same place
	static G4int hitsCID = -1;
	if(hitsCID<0){
		hitsCID = GetCollectionID(0);
	}
	hitsCE->AddHitsCollection(hitsCID, fFiberCollection);
}

///aStep->GetTrack()->GetParticleDefinition()->GetParticleName() == "e+" && 
G4bool FiberSD::ProcessHits(G4Step *aStep, G4TouchableHistory*){
	if(aStep->GetTrack()->GetTrackID() == 1 && aStep->GetTrack()->GetDynamicParticle()->GetParticleDefinition()->GetParticleName() != "gamma"){
		G4double edep = aStep->GetTotalEnergyDeposit();
		G4double delta = aStep->GetPostStepPoint()->GetKineticEnergy() - aStep->GetPreStepPoint()->GetKineticEnergy() + edep; 

		fEdep += edep;
		fDelta -= delta;
		
		G4double ein = 0, eout = 0;
		G4StepPoint* preStep = aStep->GetPreStepPoint();
		G4StepPoint* postStep = aStep->GetPostStepPoint();
		
		// Counting photons
		const std::vector<const G4Track*>* secondaries = aStep->GetSecondaryInCurrentStep();
		if(secondaries->size() > 0){
			for(unsigned int i = 0; i < secondaries->size(); i++){
				if(secondaries->at(i)->GetParentID() > 0){
					if(secondaries->at(i)->GetDynamicParticle()->GetParticleDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()){
						fNgamma += 1;
						
					}
				}
			}
		}

		// Saving e+ characteristics				
		if(aStep->IsFirstStepInVolume() && fEin == 0){
			ein = preStep->GetKineticEnergy();
			fEin = ein;
		}
		
		eout = postStep->GetKineticEnergy();

		if(postStep->GetStepStatus() == fGeomBoundary){
//			aStep->GetTrack()->SetTrackStatus(fStopAndKill);
			return true;
		}
		else if (eout == 0.){
			aStep->GetTrack()->SetTrackStatus(fStopAndKill);
			return true;
		}
		
		return false;
	}
	else if(aStep->GetTrack()->GetParticleDefinition() != G4OpticalPhoton::OpticalPhotonDefinition()) {
		const std::vector<const G4Track*>* secondaries = aStep->GetSecondaryInCurrentStep();
		if(secondaries->size() > 0){
			for(unsigned int i = 0; i < secondaries->size(); i++){
				if(secondaries->at(i)->GetParentID() > 0){
					if(secondaries->at(i)->GetDynamicParticle()->GetParticleDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()){
						fNgamma += 1;
					}
				}
			}
		}
		return false;
	}
	return false;
}

void FiberSD::EndOfEvent(G4HCofThisEvent*){
	FiberHit* Hit = new FiberHit();
	Hit->SetEin(fEin);
	Hit->SetEdep(fEdep);
	Hit->SetEdelta(fDelta);
	Hit->SetNgamma(fNgamma);
//	Hit->SetNgammaOut(fNgammaOut);
	fFiberCollection->insert(Hit);
	fEdep = 0;
	fEin = 0;
	fDelta = 0;
	fNgamma = 0;
	fNgammaOut = 0;
}

void FiberSD::clear(){}

void FiberSD::DrawAll(){}

void FiberSD::PrintAll(){}


	
		
