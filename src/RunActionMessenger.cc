/// \file  RunActionMessenger.cc
/// \brief Implemenatation of the RunActionMessenger class

#include "RunActionMessenger.hh"
#include "RunAction.hh"

#include "globals.hh"
#include "G4SystemOfUnits.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

RunActionMessenger::RunActionMessenger(RunAction* action) : G4UImessenger(), fRunAction(action){
	fAnalysisDirectory = new G4UIdirectory("/Analysis/");
	fAnalysisDirectory->SetGuidance("UI commands to specify analysis options.");

	fCmdFileName = new G4UIcmdWithAString("/Analysis/SetFileName", this);
	fCmdFileName->SetGuidance("Choose analysis file name.");
	fCmdFileName->SetParameterName("Name", false);
	fCmdFileName->AvailableForStates(G4State_Idle);
}

RunActionMessenger::~RunActionMessenger(){
	delete fCmdFileName;
	delete fAnalysisDirectory;
}

void RunActionMessenger::SetNewValue(G4UIcommand* command, G4String newValue){
	if(command == fCmdFileName){
		fRunAction->SetFileName(newValue);
	}
}
