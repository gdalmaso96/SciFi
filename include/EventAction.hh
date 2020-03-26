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

		inline void AddGamma(){fNgamma++;}
		inline void AddOut(){fNgammaOut++;}
		inline void AddOutW(){fOut ++;}
		inline void AddAbsW(){fAbs ++;}

	private:
		RunAction* fRunAction;

		G4int fCollIDFiber;
		G4int fCollIDSiPM;

		G4int fOut, fAbs;
		G4int fNgamma, fNgammaOut;
		G4int fPhID, fEventID;
		G4int fEvID; // to register each event just once

		TTree* fTree;
};

#endif


