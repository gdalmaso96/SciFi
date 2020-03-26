/// \file element.cc
/// \brief Main of the single element simulation

#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"

#include "G4RunManager.hh"

#include <TTree.h>
#include "G4UImanager.hh"
#include "FTFP_BERT.hh"
#include "G4OpticalPhysics.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "Randomize.hh"
#include "G4TScoreNtupleWriter.hh"

#include "G4Types.hh"

int main(int argc, char** argv){
	// Detect interactive mode (if no macro) and define UI session
	G4UIExecutive* ui = 0;
	if(argc == 1){
	ui = new G4UIExecutive(argc, argv);
	}

	G4RunManager*   runManager = new G4RunManager;


	// Set mandatory initialization classes
	runManager->SetUserInitialization(new DetectorConstruction);

	G4VModularPhysicsList* physicsList = new FTFP_BERT;
	physicsList->ReplacePhysics(new G4EmStandardPhysics_option4());
	G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics();
	opticalPhysics->SetTrackSecondariesFirst(kCerenkov, true);
	opticalPhysics->SetTrackSecondariesFirst(kScintillation, true);

	physicsList->RegisterPhysics(opticalPhysics);	
	runManager->SetUserInitialization(physicsList);

	// Set user action initialization
	runManager->SetUserInitialization(new ActionInitialization());

	//Initialize visualization
	G4VisManager* visManager = new G4VisExecutive;
	visManager->Initialize();

	// Get the pointer to the User Interface manager
	G4UImanager* UImanager = G4UImanager::GetUIpointer();


	// Process macro or start UI session
	if(!ui){
		// batch mode
		G4String command  = "/control/execute ";
				G4String fileName = argv[1];
		UImanager->ApplyCommand(command + fileName);
	}
	else{
		// interactive mode
		UImanager->ApplyCommand("/control/execute init_vis.mac");
		/*
		if(ui->IsGUI()){
		UImanager->ApplyCommand("/control/execute gui.mac");
		}
		*/
		ui->SessionStart();
		delete ui;
	}

	// Job termination
	delete visManager;
	delete runManager;
}


