/// \file PrimaryGeneratorActionMessenger.hh
/// \brief Definition of the PrimaryGeneratorActionMessenger class

#ifndef PrimaryGeneratorActionMessenger_h
#define PrimaryGeneratorActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class PrimaryGeneratorAction;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithADouble;
class G4UIcmdWithABool;

/// it implements command:
/// - /Element/det/SetCrysSize value unit
/// - /Element/det/SiPMmodel string

class PrimaryGeneratorActionMessenger : public G4UImessenger{
	public:
		PrimaryGeneratorActionMessenger(PrimaryGeneratorAction*);
		virtual ~PrimaryGeneratorActionMessenger();
		
		virtual void SetNewValue(G4UIcommand*, G4String);
	
	private:
		PrimaryGeneratorAction* fPrimaryGeneratorAction;
		
		G4UIdirectory* fMatrixDirectory;
		G4UIdirectory* fDetDirectory;
		
		G4UIcmdWithADoubleAndUnit* fBeamX;
		G4UIcmdWithADoubleAndUnit* fBeamY;
		G4UIcmdWithADoubleAndUnit* fBeam;
		G4UIcmdWithADouble*        fTheta;
		G4UIcmdWithADoubleAndUnit* fPositionX;
		G4UIcmdWithADoubleAndUnit* fPositionY;
		G4UIcmdWithABool         * fReal;
};

#endif


