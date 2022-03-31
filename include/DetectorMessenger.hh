/// \file DetectorMessenger.hh
/// \brief Definition of the DetectorMessenger class

#ifndef DetectorMessenger_h
#define DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithABool;
class G4UIcmdWithADoubleAndUnit;

/// it implements command:
/// - /Element/det/SetCrysSize value unit
/// - /Element/det/SiPMmodel string

class DetectorMessenger : public G4UImessenger{
	public:
		DetectorMessenger(DetectorConstruction*);
		virtual ~DetectorMessenger();
		
		virtual void SetNewValue(G4UIcommand*, G4String);
	
	private:
		DetectorConstruction* fDetectorConstruction;
		
		G4UIdirectory* fFiberDirectory;
		G4UIdirectory* fDetDirectory;
		
		G4UIcmdWithADoubleAndUnit* fFiberWidthCmd;
		G4UIcmdWithADoubleAndUnit* fFiberLengthCmd;
		G4UIcmdWithADoubleAndUnit* fLayerDinstanceCmd;
		G4UIcmdWithAString* fSiPMmodelCmd;
		G4UIcmdWithAString* fMaterialCmd;
		G4UIcmdWithABool* fSurfaceCmd;
};

#endif


