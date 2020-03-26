/// \file  SteppingAction.hh
/// \brief Definition of the SteppingAction class

#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"
#include <vector>

class DetectorConstruction;
class SteppingActionMessenger;

class G4Track;
class G4StepPoint;
class EventAction;

class G4OpBoundaryProcess;

class SteppingAction : public G4UserSteppingAction{
	public:

		SteppingAction(EventAction*);
		virtual ~SteppingAction();

		virtual void UserSteppingAction(const G4Step*);

	private:
		EventAction* fEvent;
};

#endif
