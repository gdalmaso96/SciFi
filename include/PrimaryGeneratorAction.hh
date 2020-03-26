/// \file  PrimaryGeneratorAction.hh
/// \brief Definition of PrimaryGeneratorAction class

#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "RunAction.hh"
#include "PrimaryGeneratorActionMessenger.hh"

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;

///Primary generator action class with particle gun
///
/// The default energy is 2 MeV, monochromatic positron

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction{
    public:
        PrimaryGeneratorAction();
        virtual ~PrimaryGeneratorAction();

        //method of the base class
        virtual void GeneratePrimaries(G4Event*);
	const G4ParticleGun* GetParticleGun() const {return fParticleGun;}

	inline void SetBeamX(G4double val){fBeamX = val;}
	inline void SetBeamY(G4double val){fBeamY = val;}
	inline void SetTheta(G4double val){fTheta = val;}
	inline void SetPositionX(G4double val){fPositionX = val;}
	inline void SetPositionY(G4double val){fPositionY = val;}

	inline void SetReal(G4bool val){
		fReal = val;
		if(val){
			((RunAction*) G4RunManager::GetRunManager()->GetUserRunAction())->SetGunTimeMean(1 / (1.98));
		}
		else{
			((RunAction*) G4RunManager::GetRunManager()->GetUserRunAction())->SetGunTimeMean(1 / (1.9));
		}
	}

	G4double Michel();
    
    private:
        G4ParticleGun* fParticleGun; // pointer to G4 gun class
	PrimaryGeneratorActionMessenger* fMessenger;
	G4double fBeamX, fBeamY, fTheta, fPositionX, fPositionY, fProbability;

	G4bool fReal;
};

#endif


