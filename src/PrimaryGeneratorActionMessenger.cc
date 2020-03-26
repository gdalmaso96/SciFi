/// \file PrimaryGeneratorActionMessenger.cc
/// \brief Implementation of the PrimaryGeneratorActionMessenger class

#include "PrimaryGeneratorActionMessenger.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithABool.hh"


PrimaryGeneratorActionMessenger::PrimaryGeneratorActionMessenger(PrimaryGeneratorAction* Det) : G4UImessenger(), fPrimaryGeneratorAction(Det){
	fBeamX = new G4UIcmdWithADoubleAndUnit("/Matrix/BeamX", this);
	fBeamX->SetGuidance("Set beam sigma x");
	fBeamX->SetUnitCategory("Length");
	fBeamX->SetParameterName("size_x", false);
	fBeamX->AvailableForStates(G4State_Idle);

	fBeamY = new G4UIcmdWithADoubleAndUnit("/Matrix/BeamY", this);
	fBeamY->SetGuidance("Set beam sigma y");
	fBeamY->SetUnitCategory("Length");
	fBeamY->SetParameterName("size_y", false);
	fBeamY->AvailableForStates(G4State_Idle);

	fBeam = new G4UIcmdWithADoubleAndUnit("/Matrix/Beam", this);
	fBeam->SetGuidance("Set beam sigma. It will set same sigma for x and y");
	fBeam->SetUnitCategory("Length");
	fBeam->SetParameterName("size", false);
	fBeam->AvailableForStates(G4State_Idle);

	fTheta = new G4UIcmdWithADouble("/Matrix/Theta", this);
	fTheta->SetGuidance("Set beam theta. It rotates the beam around z direction");
	fTheta->SetParameterName("theta", false);
	fTheta->AvailableForStates(G4State_Idle);

	fPositionX = new G4UIcmdWithADoubleAndUnit("/Matrix/PositionX", this);
	fPositionX->SetGuidance("Set beam position x");
	fPositionX->SetUnitCategory("Length");
	fPositionX->SetParameterName("position_x", false);
	fPositionX->AvailableForStates(G4State_Idle);

	fPositionY = new G4UIcmdWithADoubleAndUnit("/Matrix/Positiony", this);
	fPositionY->SetGuidance("Set beam position y");
	fPositionY->SetUnitCategory("Length");
	fPositionY->SetParameterName("position_y", false);
	fPositionY->AvailableForStates(G4State_Idle);

	fReal = new G4UIcmdWithABool("/Matrix/Real", this);
	fReal->SetGuidance("Active mixed muon and positron beam");
	fReal->SetParameterName("real", false);
	fReal->AvailableForStates(G4State_Idle);
}

PrimaryGeneratorActionMessenger::~PrimaryGeneratorActionMessenger(){
	delete fBeamX;
	delete fBeamY;
	delete fBeam;
	delete fTheta;
	delete fPositionX;
	delete fPositionY;
	delete fReal;
}

void PrimaryGeneratorActionMessenger::SetNewValue(G4UIcommand* command, G4String newValue){
	if(command == fBeamX){
		fPrimaryGeneratorAction->SetBeamX(fBeamX->GetNewDoubleValue(newValue));
	}
	else if(command == fBeamY){
		fPrimaryGeneratorAction->SetBeamY(fBeamY->GetNewDoubleValue(newValue));
	}
	else if(command == fBeam){
		fPrimaryGeneratorAction->SetBeamX(fBeamX->GetNewDoubleValue(newValue));
		fPrimaryGeneratorAction->SetBeamY(fBeamY->GetNewDoubleValue(newValue));
	}
	else if(command == fTheta){
		fPrimaryGeneratorAction->SetTheta(fTheta->GetNewDoubleValue(newValue));
	}
	else if(command == fPositionX){
		fPrimaryGeneratorAction->SetPositionX(fPositionX->GetNewDoubleValue(newValue));
	}
	else if(command == fPositionY){
		fPrimaryGeneratorAction->SetPositionY(fPositionY->GetNewDoubleValue(newValue));
	}
	else if(command == fReal){
		fPrimaryGeneratorAction->SetReal(fReal->GetNewBoolValue(newValue));
	}
}



