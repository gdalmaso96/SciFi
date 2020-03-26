/// \file  DetectorConstruction.hh
/// \brief Definition or DetectorConstruction class

#ifndef DetectorConstruction_h
#define DecectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4Cache.hh"
#include "G4RunManager.hh"

#include "FiberSD.hh"
#include "PixelSD.hh"

class G4VPhysicalVolume;
class G4VLogicallVolume;
class G4Box;
class G4GenericMessenger;
class G4Material;
class G4Element;
class G4MaterialPropertiesTable;
class DetectorMessenger;

/// Detector construction class to define geometry

class DetectorConstruction : public G4VUserDetectorConstruction{
    
    public:
        DetectorConstruction();
        virtual ~DetectorConstruction();

    public:
        virtual G4VPhysicalVolume* Construct();
        virtual void ConstructSDandField();

    public:
	void SetMaterial(G4String);
	void SetSurface(G4bool);
	void SetWidth(G4double);
	void SetLength(G4double);
	void SetSiPMmodel(G4String);

    private:
        // methods	
	
	DetectorMessenger* fDetectorMessenger;
	void DefineMaterials();

	G4double fFiberWidth, fFiberLength;

	G4String fmodel;
	G4int fNbOfPixelsX, fNbOfPixelsY, fNbOfPixels;
	G4double fSiPM_sizeXY, fSiPM_sizeZ, fSiPM_windowZ;

	G4LogicalVolume* fLogicFiber;
	G4LogicalVolume* fLogicPixel;

	G4Element* fH;
	G4Element* fC;
	G4Element* fN;
	G4Element* fO;
	G4Element* fSie;
	G4Material* fSi;
	G4Material* fSiResin;
	G4Material* fEpResin;
	G4Material* fMaterial;
	G4Material* fMaterialWindow;
	G4Material* fAir;
	G4Material* fVacuum;
	G4Material* fBCF10;
	G4Material* fBCF12;
	G4Material* fBCF20;
	G4Material* fFClad;
	G4Material* fSClad;

	G4VPhysicalVolume* DefineVolumes();
	G4MaterialPropertiesTable* fBCF10_mt;
	G4MaterialPropertiesTable* fBCF12_mt;
	G4MaterialPropertiesTable* fBCF20_mt;

	G4bool fCheckOverlaps, fCmdSurface;

	G4Cache<FiberSD*> fFiber_SD;
	G4Cache<PixelSD*> fPixel_SD;

};

#endif
