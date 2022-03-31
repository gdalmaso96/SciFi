/// \file  EventAction.hh
/// \brief Definition of the EventAction class

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

#include "RunAction.hh"

#include "TTree.h"
#include "TFile.h"

class EventAction : public G4UserEventAction{
	public:
		EventAction(RunAction* runAction);
		virtual ~EventAction();

		virtual void BeginOfEventAction(const G4Event*);
		virtual void   EndOfEventAction(const G4Event*);

		inline void SetThetaIn(G4double val){fThetaIn = val;}
		inline void AddTrackLength(G4double val){fTrackLength += val;}
		inline void SetPrimaryChannel(G4int val){fPrimaryChannel = val;}
		inline void SetSecondaryChannel(G4int val){fSecondaryChannel = val;}
		inline void SetSurfIn(G4int val){fSurfIn = val;}

		inline G4double GetThetaIn(){return fThetaIn;}
		inline G4double GetTrackLength(){return fTrackLength;}
		inline G4int GetPrimaryChannel(){return fPrimaryChannel;}
		inline G4int GetSecondaryChannel(){return fSecondaryChannel;}
		inline G4int GetSurfIn(){return fSurfIn;}

		inline void AddGamma(){fNgamma++;}
		inline void AddOut(){fNgammaOut++;}
		inline void AddOutW(){fOut ++;}
		inline void AddAbsW(){fAbs ++;}

	private:
		RunAction* fRunAction;

		G4int fCollIDFiber;
		G4int fCollIDSiPM;

		G4double fThetaIn, fTrackLength;
		G4int fPrimaryChannel;
		G4int fSecondaryChannel;
		G4int fSurfIn;

		G4int fOut, fAbs;
		G4int fNgamma, fNgammaOut;
		G4int fPhID, fEventID;
		G4int fEvID; // to register each event just once

		TTree* fTree;
};

#endif


