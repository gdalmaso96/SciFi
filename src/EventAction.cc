/// \file  EventAction.cc
/// \brief Implementation of the EventAction class


#include "EventAction.hh"
#include "FiberHit.hh"
#include "PixelHit.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"

#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4UnitsTable.hh"
#include "G4THitsMap.hh"
#include "G4SystemOfUnits.hh"

EventAction::EventAction(RunAction* runAction) : G4UserEventAction(), fRunAction(runAction), fCollIDFiber(-1), fCollIDSiPM(-1), fOut(false), fAbs(false), fPhID(0), fEventID(0), fEvID(-1){
}

EventAction::~EventAction(){}

void EventAction::BeginOfEventAction(const G4Event* event){
	fNgamma    = 0;
	fNgammaOut = 0;
	fOut = 0;
	fAbs = 0;
	fEventID = event->GetEventID();
}

void EventAction::EndOfEventAction(const G4Event* event){
	// Hits collections
	G4HCofThisEvent*HCE = event->GetHCofThisEvent();
	if(!HCE) return;
	// Get hits collections IDs
	if(fCollIDFiber < 0){
		G4SDManager* SDMan = G4SDManager::GetSDMpointer();
		fCollIDFiber = SDMan->GetCollectionID("fiberCollection");
	}
	
	FiberHitsCollection* FiberHitCollection = (FiberHitsCollection*) (HCE->GetHC(fCollIDFiber));

	if(fCollIDSiPM < 0){
		G4SDManager* SDMan = G4SDManager::GetSDMpointer();
		fCollIDSiPM = SDMan->GetCollectionID("pixelCollection");
	}
	
	PixelHitsCollection* PixelHitCollection = (PixelHitsCollection*) (HCE->GetHC(fCollIDSiPM));

	FiberHit* fiberHit;
	G4int N = FiberHitCollection->entries();
	PixelHit* pixelHit;
	G4int N1 = PixelHitCollection->entries();
	assert(N == N1);
	
	for(int i = 0; i < N; i++){
		fiberHit = (*FiberHitCollection)[i];
		pixelHit = (*PixelHitCollection)[i];
		if(fEvID < 0){
			fEvID = fEventID;
			fRunAction->SetEin(fiberHit->GetEin());
			fRunAction->SetEdep(fiberHit->GetEdep());
			fRunAction->SetEdelta(fiberHit->GetEdelta());
			fRunAction->SetID(fEvID);
			fRunAction->SetNgamma(fiberHit->GetNgamma() + fNgamma);
			fRunAction->SetNgammaOut(fNgammaOut);
			fRunAction->SetNAbs(fAbs);
			fRunAction->SetNOut(fOut);
			//fRunAction->SetNTransition(fTransition);
			fRunAction->SetChannel(pixelHit->GetChannel());
			fRunAction->SetCells(pixelHit->GetCells());
			fRunAction->SetCellTime(pixelHit->GetCellTime());
			fRunAction->SetOCTFlag(pixelHit->GetOCTFlag());
			fRunAction->SetDNFlag(pixelHit->GetDNFlag());
			
			(fRunAction->GetTreePtr())->Fill();
		}
		fiberHit->Clear();
		pixelHit->Clear();
	}
	if(fEvID%100 == 0) std::cout << "Event n.: " << fEvID << std::endl;
	fRunAction->AdvanceGunTime();
	fEvID = -1;
}
