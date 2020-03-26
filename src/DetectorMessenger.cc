/// \file DetectorMessenger.cc
/// \brief Implementation of the DetectorMessenger class

#include "DetectorMessenger.hh"
#include "DetectorConstruction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"


DetectorMessenger::DetectorMessenger(DetectorConstruction* Det) : G4UImessenger(), fDetectorConstruction(Det){
	fFiberDirectory = new G4UIdirectory("/SciFi/");
	fFiberDirectory->SetGuidance("UI commands specific to this simulation.");
	
	fDetDirectory = new G4UIdirectory("SciFi/det/");
	fDetDirectory->SetGuidance("Detector construction control");
	
	fFiberWidthCmd = new G4UIcmdWithADoubleAndUnit("/SciFi/det/FiberWidth", this);
	fFiberWidthCmd->SetGuidance("Define Fiber width");
	fFiberWidthCmd->SetParameterName("FiberWidth", false);
	fFiberWidthCmd->SetUnitCategory("Length");
	fFiberWidthCmd->AvailableForStates(G4State_Idle);

	fSiPMmodelCmd = new G4UIcmdWithAString("/SciFi/det/SiPMmodel", this);
	fSiPMmodelCmd->SetGuidance("Set SiPM model");
	fSiPMmodelCmd->SetParameterName("model", false);
	fSiPMmodelCmd->SetCandidates("75PE || 50PE || 25PE || 75CS || 50CS || 25CS");
	fSiPMmodelCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
	
	fFiberLengthCmd = new G4UIcmdWithADoubleAndUnit("/SciFi/det/FiberLength", this);
	fFiberLengthCmd->SetGuidance("Define Fiber length");
	fFiberLengthCmd->SetParameterName("FiberLength", false);
	fFiberLengthCmd->SetUnitCategory("Length");
	fFiberLengthCmd->AvailableForStates(G4State_Idle);
	
	fMaterialCmd = new G4UIcmdWithAString("/SciFi/det/Material", this);
	fMaterialCmd->SetGuidance("Set scintillating material");
	fMaterialCmd->SetParameterName("material", false);
	fMaterialCmd->SetCandidates("BCF10 || BCF12 || BCF20");
	fMaterialCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

	fSurfaceCmd = new G4UIcmdWithABool("/SciFi/det/Surface", this);
	fSurfaceCmd->SetGuidance("Activate ground surface");
	fSurfaceCmd->SetParameterName("surface", false);
	fSurfaceCmd->AvailableForStates(G4State_Idle);
}

DetectorMessenger::~DetectorMessenger(){
	delete fFiberWidthCmd;
	delete fFiberLengthCmd;
	delete fSiPMmodelCmd;
	delete fMaterialCmd;
	delete fSurfaceCmd;
	delete fDetDirectory;
	delete fFiberDirectory;
}

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue){
	if(command == fFiberWidthCmd){
		fDetectorConstruction->SetWidth(fFiberWidthCmd->GetNewDoubleValue(newValue));
	}
	else if (command == fFiberLengthCmd){
		fDetectorConstruction->SetLength(fFiberLengthCmd->GetNewDoubleValue(newValue));
	}
	else if (command == fMaterialCmd){
		fDetectorConstruction->SetMaterial(newValue);
	}
	else if (command == fSurfaceCmd){
		fDetectorConstruction->SetSurface(fSurfaceCmd->GetNewBoolValue(newValue));
	}
	else if(command == fSiPMmodelCmd){
		fDetectorConstruction->SetSiPMmodel(newValue);
	}
}

