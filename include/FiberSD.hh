/// \file  FiberSD.hh
/// \brief Definition of the FiberSD class

#ifndef FiberSD_h
#define FiberSD_h 1

#include "FiberHit.hh"

#include "G4ThreeVector.hh"
#include "G4VSensitiveDetector.hh"
#include "G4OpticalPhoton.hh"

#include "TVector3.h"

class G4Step;
class G4HCofThisEvent;
class G4VLogicalVolume;

class FiberSD : public G4VSensitiveDetector{
	public:
		FiberSD(G4String name);
		virtual ~FiberSD();

		virtual void Initialize(G4HCofThisEvent*);
		virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory*);
		virtual void EndOfEvent(G4HCofThisEvent*);
		virtual void clear();
		virtual void DrawAll();
		virtual void PrintAll();
		
	private:
		FiberHitsCollection* fFiberCollection;
		G4double fEin, fEdep, fDelta;
		G4int fNgamma, fNgammaOut;
};

#endif


