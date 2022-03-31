/// \file  RunAction.cc
/// \brief Implementation of the RunAction class

#include "RunAction.hh"
#include "RunActionMessenger.hh"
#include "PrimaryGeneratorAction.hh"

#include "TFile.h"
#include "TTree.h"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"


RunAction::RunAction() : G4UserRunAction(), fData(nullptr), fTree(nullptr), fCmdOCT(false), fCmdDN(false), fGunTime(0), fDNTime(0), fGunTimeMean(1/(1.9e9*CLHEP::hertz)), fDNTimeMean(1/(90*CLHEP::kilohertz)), fName("./data.root"){
	fMessenger = new RunActionMessenger(this);
}

RunAction::~RunAction(){
	delete fMessenger;
}

void RunAction::BeginOfRunAction(const G4Run*){
	fGunTime = 0;
	fDNTime.resize(84);
	for(int i = 0; i < 84; i++){
		fDNTime.at(i) = 0;
		this->AdvanceDNTime(i);
	}

	fData = TFile::Open(fName, "RECREATE");
	// New Tree
	fTree = new TTree();

	// Defining tree branches
	fTree->Branch( "ein",  &fEin);
	fTree->Branch("edep", &fEdep);
	fTree->Branch("delta", &fDelta);
	fTree->Branch("ThetaIn", &fThetaIn);
	fTree->Branch("TrackLength", &fTrackLength);
	fTree->Branch("eventID", &fID);
	fTree->Branch("Ngamma", &fNgamma);
	fTree->Branch("NgammaOut", &fNgammaOut);
	fTree->Branch("PrimaryChannel", &fPrimaryChannel);
	fTree->Branch("SecondaryChannel", &fSecondaryChannel);
	fTree->Branch("SecondaryID", &fSecondaryID);
	fTree->Branch("SurfIn", &fSurfIn);
	fTree->Branch("Out", &fOut);
	fTree->Branch("Abs", &fAbs);
//	fTree->Branch("Transition", &fTransition);
	fTree->Branch("Channel", &fChannel);
	fTree->Branch("Cells", &fCells);
	fTree->Branch("CellTime", &fCellTime);
	fTree->Branch("OCTflag", &fOCTflag);
	fTree->Branch("DNflag", &fDNflag);
	fTree->Branch("GunTime", &fGunTime);
}

void RunAction::EndOfRunAction(const G4Run*){
	fData->cd();
	//fTree->Print();
	fTree->Write("T");
	fData->Close();
}


