/// \file  PixelSD.cc
/// \brief Implementation of the PixelSD class

#include "PixelSD.hh"
#include "PixelHit.hh"
#include "RunAction.hh"

#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4VTouchable.hh"
#include "G4ios.hh"
#include "G4VProcess.hh"
#include "G4RunManager.hh"
#include "G4GeometryTolerance.hh"
#include "G4Box.hh"

#include "Randomize.hh"


class DetectorConstruction;


PixelSD::PixelSD(G4String name, G4String model) : G4VSensitiveDetector(name), fChannel(0), fVoltage(56), fDetEffGain(0), fPhotonGain(0), fOCT(0){
	fPixelCollection = nullptr;
	collectionName.insert("pixelCollection");
	
	filename[0] = "";
	filename[1] = "";
	filename[2] = "";
	filename[3] = "";
	fNbOfPixels = 0;
	fModel = model;
}

PixelSD::~PixelSD(){}

G4double PixelSD::GetAbsProbability(G4double val){
	G4double temp = 0;
	G4int N =fDetEffx.size();
	for(int i = 0; i < N; i++){
		if(i == 0){
			if((fDetEffx.at(1) - fDetEffx.at(0)) / 2 > val - fDetEffx.at(0)){
				if(val > fDetEffx.at(0)) temp = fDetEff.at(0);
				break;
			}
		}

		else if(i == N - 1){
			if((fDetEffx.at(N - 1) - fDetEffx.at(N - 2)) / 2 > fDetEffx.at(N - 1) - val){
				if(val < fDetEffx.at(N - 1)) temp = fDetEff.at(N - 1);
				break;
			}
		}

		else{
			if((fDetEffx.at(i) - fDetEffx.at(i - 1)) / 2 > fDetEffx.at(i) - val && 
			   (fDetEffx.at(i + 1) - fDetEffx.at(i)) / 2 > val - fDetEffx.at(i)){
				temp = fDetEff.at(i);
				break;
			}
		}
	}
	return temp * fDetEffGain;
}

void PixelSD::SetDetEffGain(){
	std::vector<G4double> DetEffGain;
	std::vector<G4double> DetEffGainx;
	G4double x, y;

	std::ifstream myfile;
	myfile.open("../tables/" + this->GetFileName(1));
	while(true){
		myfile >> x >> y;
		DetEffGainx.push_back(x); // x is a voltage (over the breakdown threashold)
		DetEffGain.push_back(y);
		if(myfile.eof()) break;
	}
	myfile.close();

	G4double temp = 0;
	G4int N =DetEffGainx.size();
	G4double val = fVoltage - 53;
	for(int i = 0; i < N; i++){
		if(i == 0){
			if((DetEffGainx.at(1) - DetEffGainx.at(0)) / 2 > val - DetEffGainx.at(0)){
				temp = DetEffGain.at(0);
				break;
			}
		}

		else if(i == N - 1){
			if((DetEffGainx.at(N - 1) - DetEffGainx.at(N - 2)) / 2 > DetEffGainx.at(N - 1) - val){
				temp = DetEffGain.at(N - 1);
				break;
			}
		}
		else{
			if((DetEffGainx.at(i) - DetEffGainx.at(i - 1)) / 2 > DetEffGainx.at(i) - val && 
			   (DetEffGainx.at(i + 1) - DetEffGainx.at(i)) / 2 > val - DetEffGainx.at(i)){
				temp = DetEffGain.at(i);
				break;
			}
		}
	}

	fDetEffGain = temp / 0.5;
	DetEffGain.clear();
	DetEffGainx.clear();	
}


void PixelSD::SetPhotonGain(){
	std::vector<G4double> PhotonGain;
	std::vector<G4double> PhotonGainx;
	G4double x, y;
	std::ifstream myfile;
	myfile.open("../tables/" + this->GetFileName(3));
	while(true){
		myfile >> x >> y;
		PhotonGainx.push_back(x); // x is a voltage (over the breakdown threashold)
		PhotonGain.push_back(y);
		if(myfile.eof()) break;
	}
	myfile.close();

	G4double temp = 0;
	G4int N =PhotonGainx.size();
	G4double val = fVoltage - 53;
	for(int i = 0; i < N; i++){
		if(i == 0){
			if((PhotonGainx.at(1) - PhotonGainx.at(0)) / 2 > val - PhotonGainx.at(0)){
				temp = PhotonGain.at(0);
				break;
			}
		}

		else if(i == N - 1){
			if((PhotonGainx.at(N - 1) - PhotonGainx.at(N - 2)) / 2 > PhotonGainx.at(N - 1) - val){
				temp = PhotonGain.at(N - 1);
				break;
			}
		}

		else{
			if((PhotonGainx.at(i) - PhotonGainx.at(i - 1)) / 2 > PhotonGainx.at(i) -val                    &&  (PhotonGainx.at(i + 1) - PhotonGainx.at(i)) / 2 > val - PhotonGainx.at(i)){
				temp = PhotonGain.at(i);
				break;
			}
		}
	}

	fPhotonGain = temp;
	PhotonGain.clear();
	PhotonGainx.clear();
}

void PixelSD::SetOCT(){
	std::vector<G4double> OCT;
	std::vector<G4double> OCTx;
	G4double x, y;
	std::ifstream myfile;
	myfile.open("../tables/" + this->GetFileName(2));
	while(true){
		myfile >> x >> y;
		OCTx.push_back(x); // x is a voltage (over the breakdown threashold)
		OCT.push_back(y);
		if(myfile.eof()) break;
	}
	myfile.close();

	G4double temp = 0;
	G4int N =OCTx.size();
	G4double val = fVoltage - 53;
	for(int i = 0; i < N; i++){
		if(i == 0){
			if((OCTx.at(1) - OCTx.at(0)) / 2 > val - OCTx.at(0)){
				temp = OCT.at(0);
				break;
			}
		}

		else if(i == N - 1){
			if((OCTx.at(N - 1) - OCTx.at(N - 2)) / 2 > OCTx.at(N - 1) - val){
				temp = OCT.at(N - 1);
				break;
			}
		}

		else{
			if((OCTx.at(i) - OCTx.at(i - 1)) / 2 > OCTx.at(i) -val                    &&  (OCTx.at(i + 1) - OCTx.at(i)) / 2 > val - OCTx.at(i)){
				temp = OCT.at(i);
				break;
			}
		}
	}

	fOCT = temp;
	OCT.clear();
	OCTx.clear();
}

void PixelSD::DefineProperties(){
	G4int i = 0, j = 0; //i is for the pitch, j is for the window material;
	if (fModel == "75PE") {i = 0; j = 1;}
	else if (fModel == "50PE") {i = 1; j = 1;}
	else if (fModel == "25PE") {i = 2; j = 1;}
	else if (fModel == "75CS") {i = 0; j = 0;}
	else if (fModel == "50CS") {i = 1; j = 0;}
	else if (fModel == "25CS") {i = 2; j = 0;};

	SetNbOfPixels(Model::NbPixelsX[i] * Model::NbPixelsY[i], Model::NbPixelsX[i], Model::NbPixelsY[i]);

	filename[0] = Model::eff_name[i + j * 3];
	filename[1] = Model::eff_gain_name[i];
	filename[2] = Model::OCT_gain_name[i];
	filename[3] = Model::photon_gain_name[i];
		
	std::ifstream myfile;
	myfile.open("../tables/" + this->GetFileName(0));
	G4double x, y;
	while(true){
		myfile >> x >> y;
		fDetEffx.push_back(1239.84197/x); // x is the light wave length
		fDetEff.push_back(y);
		if(myfile.eof()) break;
	}
	myfile.close();
	std::reverse(fDetEffx.begin(), fDetEffx.end());
	std::reverse(fDetEff.begin(), fDetEff.end());

	fFillFactor = Model::FillFactor[i];
	
	SetVoltage(53 + Model::OVoltage[i]);
	SetDetEffGain();
	SetPhotonGain();
	SetOCT();
	std::cout << "OCT = " << fOCT << std::endl;
	fOCT = fOCT * Model::OCT_factor[i];
	std::cout << "OCT times factor = " << fOCT << std::endl;
}

void PixelSD::Initialize(G4HCofThisEvent* hitsCE){
	fPixelCollection = new PixelHitsCollection(SensitiveDetectorName, collectionName[0]);

	static G4int hitsCID =-1;
	if(hitsCID < 0) hitsCID = GetCollectionID(0);
	hitsCE->AddHitsCollection(hitsCID, fPixelCollection);

	fCmdOCT = ((RunAction*) G4RunManager::GetRunManager()->GetUserRunAction())->GetCmdOCT();
	fCmdDN  = ((RunAction*) G4RunManager::GetRunManager()->GetUserRunAction())->GetCmdDN();
}


G4bool PixelSD::ProcessHits(G4Step *aStep, G4TouchableHistory*){
	if(aStep->GetTrack()->GetParticleDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()){

		G4TouchableHandle theTouchable = aStep->GetPreStepPoint()->GetTouchableHandle();
		G4double kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
		G4ThreeVector stppos1= aStep->GetPreStepPoint()->GetPosition();
    		G4ThreeVector localpos1 = 
      			theTouchable->GetHistory()->GetTopTransform().TransformPoint(stppos1);
		G4VPhysicalVolume* physVol = aStep->GetPreStepPoint()->GetTouchable()->GetVolume(0);
		G4double dimensions = - ((G4Box*) physVol->GetLogicalVolume()->GetSolid())->GetZHalfLength();
		G4double dimensionsy = - ((G4Box*) physVol->GetLogicalVolume()->GetSolid())->GetYHalfLength();
		G4double dimensionsx = - ((G4Box*) physVol->GetLogicalVolume()->GetSolid())->GetXHalfLength();

		G4ThreeVector momentum = aStep->GetPostStepPoint()->GetMomentumDirection();
		G4AffineTransform momentumTransform = theTouchable->GetHistory()->GetTopTransform();
		momentumTransform.SetNetTranslation(G4ThreeVector(0,0,0));
		G4ThreeVector localMomentum = momentumTransform.TransformPoint(momentum);

		if(aStep->GetPreStepPoint()->GetStepStatus() == fGeomBoundary && 
		   std::fabs(localpos1.z() + dimensions) < kCarTolerance){

//		if(aStep->GetPreStepPoint()->GetStepStatus() == fGeomBoundary){
			G4double energy = aStep->GetPreStepPoint()->GetKineticEnergy();
			G4double a = G4UniformRand();
			G4double p = GetAbsProbability(energy/CLHEP::eV);
/*
			std::cout << "pos w " << stppos1 << std::endl;
			std::cout << "pos l " << localpos1 << std::endl;
			std::cout << "mom w " << momentum << std::endl;
			std::cout << "mom l " << localMomentum << std::endl;
			std::cout << "fiber " << aStep->GetPreStepPoint()->GetTouchable()->GetReplicaNumber(2) << std::endl;
			std::cout << "sipm " << aStep->GetPreStepPoint()->GetTouchable()->GetReplicaNumber(1) << std::endl << std::endl;
*/
			fOCTflag = 0;
			if(a < p){ // SiPM Fill Factor
				
				fChannel = aStep->GetPreStepPoint()->GetTouchable()->GetReplicaNumber(2) * 2 + aStep->GetPreStepPoint()->GetTouchable()->GetReplicaNumber(1);

				RunAction* action = (RunAction*) G4RunManager::GetRunManager()->GetUserRunAction();
				G4double ptime = aStep->GetPreStepPoint()->GetGlobalTime() + action->GetGunTime();
				G4int replica = int((localpos1.x() + 1.3/2)/1.3 * fNbOfPixelsX) + 
				      fNbOfPixelsY * int((localpos1.y() + 1.3/2)/1.3 * fNbOfPixelsY);

			
				if(fCmdDN){
					while(action->GetDNTime(fChannel) < ptime){
						G4int DNcell = int(G4UniformRand() * this->GetNbOfPixels());
						fChannelvec.push_back(fChannel);
						fCells.push_back(DNcell);
						fCellTime.push_back(action->GetDNTime(fChannel));
						if(fCmdOCT && G4UniformRand() < fOCT){
							G4double cost = G4UniformRand() * 2 - 1;
							G4double sint = sqrt(1 - cost * cost);
							G4double phi  = G4UniformRand() * 2 * CLHEP::pi;
							G4ThreeVector dir = G4ThreeVector(sint * cos(phi), sint * sin(phi), cost);
							G4DynamicParticle* dynOCT = new G4DynamicParticle(G4OpticalPhoton::OpticalPhotonDefinition(), dir, 2*CLHEP::eV);
							G4Track* OCT = new G4Track(dynOCT, aStep->GetPreStepPoint()->GetGlobalTime(), aStep->GetPreStepPoint()->GetPosition());
							OCT->SetParentID(-replica - 1);
							G4TrackVector* newTrack = aStep->GetfSecondary();
							newTrack->push_back(OCT);
							fOCTflagvec.push_back(1);
						}
						else fOCTflagvec.push_back(0);
						fDNflagvec.push_back(1);
						
						action->AdvanceDNTime(fChannel);
					}
				}
				if(fCmdOCT && G4UniformRand() < fOCT){
					G4double cost = G4UniformRand() * 2 - 1;
					G4double sint = sqrt(1 - cost * cost);
					G4double phi  = G4UniformRand() * 2 * CLHEP::pi;
					G4ThreeVector dir = G4ThreeVector(sint * cos(phi), sint * sin(phi), cost);
					G4DynamicParticle* dynOCT = new G4DynamicParticle(G4OpticalPhoton::OpticalPhotonDefinition(), dir, 2*CLHEP::eV);
					G4Track* OCT = new G4Track(dynOCT, aStep->GetPreStepPoint()->GetGlobalTime(), aStep->GetPreStepPoint()->GetPosition());
					OCT->SetParentID(-replica - 1);
					G4TrackVector* newTrack = aStep->GetfSecondary();
					newTrack->push_back(OCT);
				}
				fChannelvec.push_back(fChannel);
				fCells.push_back(replica);
				fCellTime.push_back(ptime);
				fOCTflagvec.push_back(fOCTflag);
				fDNflagvec.push_back(0);
			}
			aStep->GetTrack()->SetTrackStatus(fStopAndKill);
		}

		else if(aStep->GetPreStepPoint()->GetStepStatus() == fGeomBoundary && 
		        std::fabs(localpos1.z() - dimensions) < kCarTolerance &&
		        localMomentum.getZ() < 0){
			aStep->GetTrack()->SetTrackStatus(fStopAndKill);
		}
		else if(aStep->GetPreStepPoint()->GetStepStatus() == fGeomBoundary && 
		        std::fabs(localpos1.y() + dimensionsy) < kCarTolerance &&
		        localMomentum.getY() > 0){
			aStep->GetTrack()->SetTrackStatus(fStopAndKill);
		}
		else if(aStep->GetPreStepPoint()->GetStepStatus() == fGeomBoundary && 
		        std::fabs(localpos1.y() - dimensionsy) < kCarTolerance &&
		        localMomentum.getY() < 0){
			aStep->GetTrack()->SetTrackStatus(fStopAndKill);
		}		
		else if(aStep->GetPreStepPoint()->GetStepStatus() == fGeomBoundary && 
		        std::fabs(localpos1.x() + dimensionsx) < kCarTolerance &&
		        localMomentum.getX() > 0){
			aStep->GetTrack()->SetTrackStatus(fStopAndKill);
		}
		else if(aStep->GetPreStepPoint()->GetStepStatus() == fGeomBoundary && 
		        std::fabs(localpos1.x() - dimensionsx) < kCarTolerance &&
		        localMomentum.getX() < 0){
			aStep->GetTrack()->SetTrackStatus(fStopAndKill);
		}

	}
	
	return false;
}

void PixelSD::EndOfEvent(G4HCofThisEvent*){
	PixelHit* Hit = new PixelHit();
	Hit->SetChannel(fChannelvec);
	Hit->SetCells(fCells);
	Hit->SetCellTime(fCellTime);
	Hit->SetOCTFlag(fOCTflagvec);
	Hit->SetDNFlag(fDNflagvec);
	fPixelCollection->insert(Hit);
	fChannelvec.clear();
	fCells.clear();
	fCellTime.clear();
	fOCTflagvec.clear();
	fDNflagvec.clear();
}

void PixelSD::clear(){}

void PixelSD::DrawAll(){}

void PixelSD::PrintAll(){}


