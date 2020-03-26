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

	fPrimaryDirectory = new G4UIdirectory("/Primary/");
	fPrimaryDirectory->SetGuidance("UI commands to specify gun options.");

	fCmdFileName = new G4UIcmdWithAString("/Analysis/SetFileName", this);
	fCmdFileName->SetGuidance("Choose analysis file name.");
	fCmdFileName->SetParameterName("Name", false);
	fCmdFileName->AvailableForStates(G4State_Idle);

	fCmdOCT = new G4UIcmdWithABool("/Element/det/OCT", this);
	fCmdOCT->SetGuidance("Activate optical cross talk among pixels in SiPMs");
	fCmdOCT->SetParameterName("OCT", false);
	fCmdOCT->AvailableForStates(G4State_Idle);

	fCmdDN = new G4UIcmdWithABool("/Element/det/DN", this);
	fCmdDN->SetGuidance("Activate dark noise in SiPMs");
	fCmdDN->SetParameterName("DN", false);
	fCmdDN->AvailableForStates(G4State_Idle);

	fCmdGunTime = new G4UIcmdWithADoubleAndUnit("/Primary/Rate", this);
	fCmdGunTime->SetGuidance("Set beam rate.");
	fCmdGunTime->SetParameterName("rate", false);
	fCmdGunTime->SetUnitCategory("Rate");
	fCmdGunTime->SetDefaultValue(1.9e9*hertz);
	fCmdGunTime->AvailableForStates(G4State_Idle);

}

RunActionMessenger::~RunActionMessenger(){
	delete fCmdGunTime;
	delete fCmdFileName;
	delete fCmdOCT;
	delete fCmdDN;
	delete fPrimaryDirectory;
	delete fAnalysisDirectory;
}

void RunActionMessenger::SetNewValue(G4UIcommand* command, G4String newValue){
	if(command == fCmdFileName){
		fRunAction->SetFileName(newValue);
	}
	else if (command == fCmdOCT){
		fRunAction->SetCmdOCT(fCmdOCT->GetNewBoolValue(newValue));
	}
	else if (command == fCmdDN){
		fRunAction->SetCmdDN(fCmdDN->GetNewBoolValue(newValue));
	}
	else if (command == fCmdGunTime){
		fRunAction->SetGunTimeMean(1 / fCmdGunTime->GetNewDoubleValue(newValue));
	}
}
