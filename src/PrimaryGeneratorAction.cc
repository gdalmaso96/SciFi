/// \file  PrimaryGeneraotrAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class

#include "PrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "TMath.h"

#include "poly34.hh"


PrimaryGeneratorAction::PrimaryGeneratorAction() : G4VUserPrimaryGeneratorAction(), fParticleGun(0), fBeamX(2*CLHEP::cm), fBeamY(2*CLHEP::cm), fTheta(0), fPositionX(0), fPositionY(0), fReal(false){
	fProbability = 1.9 / (1.9 + 0.08);
	G4int n_particle = 1;
	fParticleGun = new G4ParticleGun(n_particle);
	fMessenger = new PrimaryGeneratorActionMessenger(this);

	//default particle kinematic
	auto particle = G4ParticleTable::GetParticleTable()->FindParticle("e+");
	fParticleGun->SetParticleDefinition(particle);
	fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., -1.));
	fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,0.));
	fParticleGun->SetParticleEnergy(2.2*MeV);
	std::ofstream myfile;
	myfile.open("out.txt");
	myfile.close();

}


PrimaryGeneratorAction::~PrimaryGeneratorAction(){
	delete fParticleGun;
	delete fMessenger;
}

G4double PrimaryGeneratorAction::Michel(){
	G4double xmin = 2.711/53.3;
	G4double A = 0.5 - pow(xmin, 3)*(1 - xmin/2);

	G4double csi = 0;
	G4int res = 0;
	G4double e = 0.5 - A*csi;
	G4double x[4] = {0., 0., 0., 0.};
	std::ofstream myfile;
	myfile.open("out.txt", std::ios::out | std::ios::app);

	while(true){
		csi = G4UniformRand();
		e = 0.5 - A*csi;
		res = SolveP4(x, 2, 0, 0, e * 2);
		if(res > 0){
			for(int i = 0; i < res; i++){
				if(-x[i] > xmin && -x[i] < 1){
					myfile << -x[i] << std::endl;
					myfile.close();
					return -x[i] * CLHEP::MeV*53.3 - 0.511*CLHEP::MeV;
				}
			}
			std::cout << std::endl;
		}
	}
}
		


void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent){
	G4ParticleDefinition* particle = fParticleGun->GetParticleDefinition();
	
	// In order to avoid dependence of PrimaryGeneratorAction
	// on DetectorConstruction class we get world volume from G4LogicalVolumeStore.

	G4double world_size = 0;
	auto worldLV = G4LogicalVolumeStore::GetInstance()->GetVolume("World");

	G4Box* worldBox = dynamic_cast<G4Box*>(worldLV->GetSolid());

	world_size = worldBox->GetXHalfLength();

	if(fReal){
		if(G4UniformRand() < fProbability){
			particle = G4ParticleTable::GetParticleTable()->FindParticle("e+");
			fParticleGun->SetParticleEnergy(Michel());
		}
		else{
			particle = G4ParticleTable::GetParticleTable()->FindParticle("mu+");
			fParticleGun->SetParticleEnergy(sqrt(105.66*105.66 + 28*28)*CLHEP::MeV - 105.66*CLHEP::MeV);
		}
	}
	// Set gun position
	double x = 0, y = 0, tempx = 0, tempy = 0;
	tempx = G4RandGauss::shoot(fPositionX, fBeamX);
	tempy = G4RandGauss::shoot(fPositionY, fBeamY);

	double rho = fTheta;
	x =  tempx*cos(rho) + tempy*sin(rho);
	y = -tempx*sin(rho) + tempy*cos(rho);
	fParticleGun->SetParticleDefinition(particle);
	fParticleGun->SetParticlePosition(G4ThreeVector(x, y, world_size));
	fParticleGun->GeneratePrimaryVertex(anEvent);
}


