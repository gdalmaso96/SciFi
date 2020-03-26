/// \file  DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "RunAction.hh"
#include "SiPMModel.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4AutoDelete.hh"

#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"

#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4GenericMessenger.hh"

/// Constructor
DetectorConstruction::DetectorConstruction() : G4VUserDetectorConstruction(), fFiberWidth(0.25*mm), fFiberLength(20*cm),fSiPM_sizeXY(1.3*mm), fCheckOverlaps(true), fCmdSurface(false){
	fDetectorMessenger = new DetectorMessenger(this);
	DefineMaterials();
	SetSiPMmodel("50CS");
	SetMaterial("BCF10");
}

/// Deconstructor
DetectorConstruction :: ~DetectorConstruction(){
	delete fDetectorMessenger;
}


/// Implementation

G4VPhysicalVolume* DetectorConstruction::Construct(){
	// Define volumes
	return DefineVolumes();
}

void DetectorConstruction::DefineMaterials(){
	G4NistManager* nist = G4NistManager::Instance();
	G4double a; // atomic mass
	G4double z; // atomic number
	G4double density;
	
	/// Elements
	fH  = new G4Element( "H",  "H", z =  1., a =   1.01*g/mole);
	fC  = new G4Element( "C",  "C", z =  6., a =  12.01*g/mole);
	fN  = new G4Element( "N",  "N", z =  7., a =  14.01*g/mole);
	fO  = new G4Element( "O",  "O", z =  8., a =  16.00*g/mole);
	fSie= new G4Element("Si", "Si", z = 14., a = 28.0855*g/mole);
	
	/// Materials
	// BCF10
	fBCF10 = new G4Material("BCF10", density = 1.05*g/cm3, 2);
	fBCF10->AddElement(fC, 485);
	fBCF10->AddElement(fH, 482);

	// BCF12
	fBCF12 = new G4Material("BCF12", density = 1.05*g/cm3, 2);
	fBCF12->AddElement(fC, 485);
	fBCF12->AddElement(fH, 482);

	// BCF20
	fBCF20 = new G4Material("BCF20", density = 1.05*g/cm3, 2);
	fBCF20->AddElement(fC, 485);
	fBCF20->AddElement(fH, 482);

	// First Cladding: PMMA
	fFClad = new G4Material("FClad", density = 1.2*g/cm3, 3);
	fFClad->AddElement(fC, 5);
	fFClad->AddElement(fH, 8);
	fFClad->AddElement(fO, 2);

	// Second Cladding: PMMA EMA
	fSClad = new G4Material("SClad", density = 1.2*g/cm3, 3);
	fSClad->AddElement(fC, 5);
	fSClad->AddElement(fH, 8);
	fSClad->AddElement(fO, 2);

	// Silicon resin
	fSiResin = new G4Material("SiResin",z=1.,a=1.01*g/mole, 
		     		 density = universe_mean_density, kStateGas,
				 0.1 * kelvin, 1.e-19 * pascal);

	// Epoxy resin
	fEpResin = new G4Material("EpoxyResin",z=1.,a=1.01*g/mole, 
		     		 density = universe_mean_density, kStateGas,
				 0.1 * kelvin, 1.e-19 * pascal);

	// Air
	fAir = new G4Material("Air", density = 0.0010*g/cm3, 2);
	fAir->AddElement(fN, 70 * perCent);
	fAir->AddElement(fO, 30 * perCent);

	// Vacuuum
	fVacuum = new G4Material("Vacuum",z=1.,a=1.01*g/mole, 
		     		 density = universe_mean_density, kStateGas,
				 0.1 * kelvin, 1.e-19 * pascal);

	/// Material properties tables
	//  BCF10 optics
	std::vector<G4double> energy, scint;
	G4double tempe, tempscint;
	std::ifstream myfile;
	myfile.open("../tables/BCF10_light_out.txt");
	while(true){
		myfile >> tempe >> tempscint;
		energy.push_back(1239.84197/tempe);
		scint.push_back(tempscint);
		
		if(myfile.eof()) break;
	}
	myfile.close();

	assert(energy.size() == scint.size());
	const G4int bcf10 = int(energy.size());

	G4double* BCF10_Energy = new G4double[bcf10];
	G4double* BCF10_SCINT = new G4double[bcf10];

	G4double* BCF10_RIND = new G4double[bcf10];
	G4double* BCF10_ABSL = new G4double[bcf10];
	
	for(int i = 0; i < bcf10; i++){
		BCF10_Energy[i] = energy.at(i)*eV;
		BCF10_SCINT[i] = scint.at(i);
		BCF10_RIND[i] = 1.6;
		BCF10_ABSL[i] = 220*cm;
	}
	
	energy.clear();
	scint.clear();
	
	
	assert(sizeof(BCF10_SCINT) == sizeof(BCF10_Energy));
	
	assert(sizeof(BCF10_RIND) == sizeof(BCF10_Energy));
	
	assert(sizeof(BCF10_ABSL) == sizeof(BCF10_Energy));

	fBCF10_mt = new G4MaterialPropertiesTable();
	fBCF10_mt->AddProperty(       "RINDEX", BCF10_Energy,  BCF10_RIND, bcf10);
	fBCF10_mt->AddProperty(    "ABSLENGTH", BCF10_Energy,  BCF10_ABSL, bcf10);
	fBCF10_mt->AddProperty("FASTCOMPONENT", BCF10_Energy, BCF10_SCINT, bcf10);
	
	fBCF10_mt->AddConstProperty("SCINTILLATIONYIELD",        8000./MeV);
	fBCF10_mt->AddConstProperty(   "RESOLUTIONSCALE",                 1);
	fBCF10_mt->AddConstProperty(  "SLOWTIMECONSTANT",            2.7*ns);
	fBCF10_mt->AddConstProperty(  "FASTTIMECONSTANT",            2.7*ns);
	fBCF10_mt->AddConstProperty(        "YIELDRATIO",                 0.);
	
	fBCF10->SetMaterialPropertiesTable(fBCF10_mt);

	//  Set Birks Constant
	fBCF10->GetIonisation()->SetBirksConstant(0.117645*mm/MeV);

	//  BCF12 optics
	myfile.open("../tables/BCF12_light_out.txt");
	while(true){
		myfile >> tempe >> tempscint;
		energy.push_back(1239.84197/tempe);
		scint.push_back(tempscint);
		
		if(myfile.eof()) break;
	}
	myfile.close();

	assert(energy.size() == scint.size());
	const G4int bcf12 = int(energy.size());

	G4double* BCF12_Energy = new G4double[bcf12];
	G4double* BCF12_SCINT = new G4double[bcf12];

	G4double* BCF12_RIND = new G4double[bcf12];
	G4double* BCF12_ABSL = new G4double[bcf12];
	
	for(int i = 0; i < bcf12; i++){
		BCF12_Energy[i] = energy.at(i)*eV;
		BCF12_SCINT[i] = scint.at(i);
		BCF12_RIND[i] = 1.6;
		BCF12_ABSL[i] = 270*cm;
	}
	
	energy.clear();
	scint.clear();
	
	
	assert(sizeof(BCF12_SCINT) == sizeof(BCF12_Energy));
	
	assert(sizeof(BCF12_RIND) == sizeof(BCF12_Energy));
	
	assert(sizeof(BCF12_ABSL) == sizeof(BCF12_Energy));

	fBCF12_mt = new G4MaterialPropertiesTable();
	fBCF12_mt->AddProperty(       "RINDEX", BCF12_Energy,  BCF12_RIND, bcf12);
	fBCF12_mt->AddProperty(    "ABSLENGTH", BCF12_Energy,  BCF12_ABSL, bcf12);
	fBCF12_mt->AddProperty("FASTCOMPONENT", BCF12_Energy, BCF12_SCINT, bcf12);
	
	fBCF12_mt->AddConstProperty("SCINTILLATIONYIELD",        8000./MeV);
	fBCF12_mt->AddConstProperty(   "RESOLUTIONSCALE",                 1);
	fBCF12_mt->AddConstProperty(  "SLOWTIMECONSTANT",            3.2*ns);
	fBCF12_mt->AddConstProperty(  "FASTTIMECONSTANT",            3.2*ns);
	fBCF12_mt->AddConstProperty(        "YIELDRATIO",                 0.);
	
	fBCF12->SetMaterialPropertiesTable(fBCF12_mt);

	//  Set Birks Constant
	fBCF12->GetIonisation()->SetBirksConstant(0.117645*mm/MeV);

	//  BCF20 optics
	myfile.open("../tables/BCF20_light_out.txt");
	while(true){
		myfile >> tempe >> tempscint;
		energy.push_back(1239.84197/tempe);
		scint.push_back(tempscint);
		
		if(myfile.eof()) break;
	}
	myfile.close();

	assert(energy.size() == scint.size());
	const G4int bcf20 = int(energy.size());

	G4double* BCF20_Energy = new G4double[bcf20];
	G4double* BCF20_SCINT = new G4double[bcf20];

	G4double* BCF20_RIND = new G4double[bcf20];
	G4double* BCF20_ABSL = new G4double[bcf20];
	
	for(int i = 0; i < bcf20; i++){
		BCF20_Energy[i] = energy.at(i)*eV;
		BCF20_SCINT[i] = scint.at(i);
		BCF20_RIND[i] = 1.6;
		BCF20_ABSL[i] = 350*cm;
	}
	
	energy.clear();
	scint.clear();
	
	
	assert(sizeof(BCF20_SCINT) == sizeof(BCF20_Energy));
	
	assert(sizeof(BCF20_RIND) == sizeof(BCF20_Energy));

	assert(sizeof(BCF20_ABSL) == sizeof(BCF20_Energy));

	fBCF20_mt = new G4MaterialPropertiesTable();
	fBCF20_mt->AddProperty(       "RINDEX", BCF20_Energy,  BCF20_RIND, bcf20);
	fBCF20_mt->AddProperty(    "ABSLENGTH", BCF20_Energy,  BCF20_ABSL, bcf20);
	fBCF20_mt->AddProperty("FASTCOMPONENT", BCF20_Energy, BCF20_SCINT, bcf20);
	
	fBCF20_mt->AddConstProperty("SCINTILLATIONYIELD",        8000./MeV);
	fBCF20_mt->AddConstProperty(   "RESOLUTIONSCALE",                 1);
	fBCF20_mt->AddConstProperty(  "SLOWTIMECONSTANT",            2.7*ns);
	fBCF20_mt->AddConstProperty(  "FASTTIMECONSTANT",            2.7*ns);
	fBCF20_mt->AddConstProperty(        "YIELDRATIO",                 0.);
	
	fBCF20->SetMaterialPropertiesTable(fBCF20_mt);

	//  Set Birks Constant
	fBCF20->GetIonisation()->SetBirksConstant(0.117645*mm/MeV);

	// Vacuum and air
	G4double vacuum_Energy[] = {1.5*eV, 4.*eV};
	const G4int vacnum = sizeof(vacuum_Energy) / sizeof(G4double);
	
	G4double vRIND = 1.;
	G4double vacuum_RIND[] = {vRIND, vRIND};
	assert(sizeof(vacuum_RIND) == sizeof(vacuum_Energy));
	
	G4MaterialPropertiesTable* vacuum_mt = new G4MaterialPropertiesTable();
	vacuum_mt->AddProperty("RINDEX", vacuum_Energy, vacuum_RIND, vacnum);
	fVacuum->SetMaterialPropertiesTable(vacuum_mt);
	fAir   ->SetMaterialPropertiesTable(vacuum_mt);

	// Cladding
	G4double fCladRIND = 1.49;
	G4double fClad_RIND[] = {fCladRIND, fCladRIND};
	assert(sizeof(fClad_RIND) = sizeof(vacuum_Energy));

	G4MaterialPropertiesTable* fClad_mt = new G4MaterialPropertiesTable();
	fClad_mt->AddProperty("RINDEX", vacuum_Energy, fClad_RIND, vacnum);
	fFClad->SetMaterialPropertiesTable(fClad_mt);


	G4double sCladRIND = 1.42;
	G4double sClad_RIND[] = {sCladRIND, sCladRIND};
	assert(sizeof(sClad_RIND) = sizeof(vacuum_Energy));

	G4MaterialPropertiesTable* sClad_mt = new G4MaterialPropertiesTable();
	sClad_mt->AddProperty("RINDEX", vacuum_Energy, sClad_RIND, vacnum);
	fSClad->SetMaterialPropertiesTable(sClad_mt);

	// Silicium
	fSi = nist->FindOrBuildMaterial("G4_Si");

	G4double Si_Energy[] = {.5*eV, 9.*eV};
	const G4int Sinum = sizeof(vacuum_Energy) / sizeof(G4double);

	G4double Si_RIND[] = {3.4, 3.4};
	assert(sizeof(Si_RIND) == sizeof(Si_Energy));

	G4MaterialPropertiesTable* Si_mt = new G4MaterialPropertiesTable();
	Si_mt->AddProperty("RINDEX", Si_Energy, Si_RIND, Sinum);
	fSi->SetMaterialPropertiesTable(Si_mt);

	// Silicon resin
	G4double SiRes_RIND[] = {1.41, 1.41};
	assert(sizeof(SiRes_RIND) == sizeof(Si_Energy));
	
	G4MaterialPropertiesTable* SiRes_mt = new G4MaterialPropertiesTable();
	SiRes_mt->AddProperty("RINDEX", Si_Energy, SiRes_RIND, Sinum);
	fSiResin->SetMaterialPropertiesTable(SiRes_mt);

	// Silicium
	G4double Ep_RIND[] = {1.55, 1.55};
	assert(sizeof(Ep_RIND) == sizeof(Si_Energy));
	
	G4MaterialPropertiesTable* Ep_mt = new G4MaterialPropertiesTable();
	Ep_mt->AddProperty("RINDEX", Si_Energy, Ep_RIND, Sinum);
	fEpResin->SetMaterialPropertiesTable(Ep_mt);

}

void DetectorConstruction::SetMaterial(G4String name){
	if(name == "BCF10") fMaterial = fBCF10;
	else if(name == "BCF12") fMaterial = fBCF12;
	else if(name == "BCF20") fMaterial = fBCF20;
	G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void DetectorConstruction::SetSiPMmodel(G4String name){
	G4int i = 0, j = 0; //i is for the pitch, j is for the window material;
	if (name == "75PE") {i = 0; j = 1;}
	else if (name == "50PE") {i = 1; j = 1;}
	else if (name == "25PE") {i = 2; j = 1;}
	else if (name == "75CS") {i = 0; j = 0;}
	else if (name == "50CS") {i = 1; j = 0;}
	else if (name == "25CS") {i = 2; j = 0;};
	
	fmodel = name;
	if(j == 1) fMaterialWindow = fEpResin;
	else if(j == 0) fMaterialWindow = fSiResin;

	fNbOfPixelsX = Model::NbPixelsX[i];
	fNbOfPixelsY = Model::NbPixelsY[i];
	
	fSiPM_sizeZ = Model::SiPM_size_Z[j];
	fSiPM_windowZ = Model::window_size_Z[j];

	if(G4RunManager::GetRunManager()->GetUserRunAction()){
		((RunAction*) G4RunManager::GetRunManager()->GetUserRunAction())->SetDNTimeMean(1 / Model::dark_noise_rate[i]);
	}
	G4RunManager::GetRunManager()->ReinitializeGeometry();
}

G4VPhysicalVolume* DetectorConstruction::DefineVolumes(){

	/// MATERIALS AND PARAMETERS
	
	// Crystal parameters

	// World
	G4double world_sizeX = fFiberLength + 2*fSiPM_sizeZ + fSiPM_windowZ + 2*mm;
	G4double world_sizeY = fFiberLength + 2*fSiPM_sizeZ + fSiPM_windowZ + 2*mm;
	G4double world_sizeZ = 2 * fSiPM_sizeXY + 7*mm;

	// World
	G4Box* SolidWorld = new G4Box("World", 0.5*world_sizeX, 0.5*world_sizeY, 0.5*world_sizeZ);
	G4LogicalVolume* logicWorld = new G4LogicalVolume(SolidWorld, fVacuum, "World");
	G4PVPlacement* physWorld = new G4PVPlacement(0, G4ThreeVector(), logicWorld, "World", 0, false, 0, fCheckOverlaps);

	// Fiber
	G4Box* solidFiber = new G4Box("Fiber", 0.5 * fFiberWidth,
				               0.5 * fFiberLength,
					       0.5 * fFiberWidth);
	G4LogicalVolume* logicFiber = new G4LogicalVolume(solidFiber, fVacuum, "Fiber");
	
	// Core
	G4Box* solidCore = new G4Box("Core", 0.5 * fFiberWidth * 0.94,
				              0.5 * fFiberLength * 0.94,
					      0.5 * fFiberWidth * 0.94);
	G4LogicalVolume* logicCore = new G4LogicalVolume(solidCore, fMaterial, "Core");
	fLogicFiber = logicCore;

	// First Cladding
	G4Box* temp = new G4Box("Temp",  0.5 * fFiberWidth * 0.98,
			                 0.5 * fFiberLength * 0.98,
				         0.5 * fFiberWidth * 0.98);
	G4SubtractionSolid* solidClad = new G4SubtractionSolid("fClad", temp, solidCore, 0, G4ThreeVector());
	G4LogicalVolume* logicFClad = new G4LogicalVolume(solidClad, fFClad, "fClad");

	// Second Cladding
	G4SubtractionSolid* solidSClad = new G4SubtractionSolid("sClad", solidFiber, temp, 0, G4ThreeVector());
	G4LogicalVolume* logicSClad = new G4LogicalVolume(solidSClad, fSClad, "sClad");

	// SiPM
	G4Box* solidSiPM = new G4Box("Pixel", 0.5 * fSiPM_sizeXY, 0.5 * fSiPM_sizeXY, 0.5 * (fSiPM_sizeZ + fSiPM_windowZ));
	G4LogicalVolume* logicSiPM = new G4LogicalVolume(solidSiPM, fVacuum, "SiPM");

	// SiPM window
	G4Box* solidWindow = new G4Box("Window", 0.5 * fSiPM_sizeXY, 0.5 * fSiPM_sizeXY, 0.5 * fSiPM_windowZ);
	G4LogicalVolume* logicWindow = new G4LogicalVolume(solidWindow, fMaterialWindow, "Window");

	// SiPM Pixel
	G4Box* solidPixel = new G4Box("Pixel", 0.5*fSiPM_sizeXY, 0.5*fSiPM_sizeXY, 0.5*fSiPM_sizeZ);
	G4LogicalVolume* logicPixel = new G4LogicalVolume(solidPixel, fSi, "PixelLV");
	fLogicPixel = logicPixel;

	// Element
	G4Box* solidElement = new G4Box("Element", 0.5 * fSiPM_sizeXY, 0.5 * (fFiberLength + fSiPM_sizeZ + fSiPM_windowZ), 0.5 * fSiPM_sizeXY);
	G4LogicalVolume* logicElement = new G4LogicalVolume(solidElement, fVacuum, "Element");
	// Placing volumes

	G4ThreeVector pixel_pos = G4ThreeVector(0, 0, -0.5 * (fSiPM_windowZ));
	new G4PVPlacement(0, pixel_pos, logicPixel, "Pixel", logicSiPM, false, 0, fCheckOverlaps);

	new G4PVPlacement(0, G4ThreeVector(0, 0, 0.5 * fSiPM_sizeZ), logicWindow, "Window", logicSiPM, false, 0, fCheckOverlaps);
	G4VPhysicalVolume* physCore = new G4PVPlacement(0,G4ThreeVector(), logicCore, "Core", logicFiber, false, 0, fCheckOverlaps);
	G4VPhysicalVolume* physFClad = new G4PVPlacement(0,G4ThreeVector(), logicFClad, "fClad", logicFiber, false, 0, fCheckOverlaps);
	G4VPhysicalVolume* physSClad = new G4PVPlacement(0,G4ThreeVector(), logicSClad, "sClad", logicFiber, false, 0, fCheckOverlaps);

	new G4PVPlacement(0,G4ThreeVector(), logicFiber, "Fiber", logicElement, false, 0, fCheckOverlaps);

	G4RotationMatrix* Rotation = new G4RotationMatrix();
	Rotation->rotateX(-90*deg);
	Rotation->rotateY(0.*deg);
	Rotation->rotateZ(0.*deg);
	new G4PVPlacement(Rotation, G4ThreeVector(0,0.5 * (fFiberLength + fSiPM_sizeZ + fSiPM_windowZ), 0), logicSiPM, "SiPM", logicElement, false, 0, fCheckOverlaps);

	Rotation = new G4RotationMatrix();
	Rotation->rotateX(90*deg);
	Rotation->rotateY(0.*deg);
	Rotation->rotateZ(0.*deg);
	new G4PVPlacement(Rotation, G4ThreeVector(0,-0.5 * (fFiberLength + fSiPM_sizeZ + fSiPM_windowZ),0), logicSiPM, "SiPM", logicElement, false, 1, fCheckOverlaps);

	G4double fibersDinstance = 5*mm;
	// Vertical Fibers
	for(int i = 0; i < 21; i++){
		G4ThreeVector posFiber = G4ThreeVector(-5*cm + fibersDinstance * i, 0, -2.5*mm);
		new G4PVPlacement(0, posFiber, logicElement, "Element", logicWorld, false, i, fCheckOverlaps);
	}

	// Orizontal Fibers
	Rotation = new G4RotationMatrix();
	Rotation->rotateX(0.*deg);
	Rotation->rotateY(0.*deg);
	Rotation->rotateZ(90*deg);
	for(int i = 0; i < 21; i++){
		G4ThreeVector posFiber = G4ThreeVector(0, 5*cm - fibersDinstance * i, 2.5*mm);
		new G4PVPlacement(Rotation, posFiber, logicElement, "Element", logicWorld, false, i + 21, fCheckOverlaps);
	}		

	logicCore->SetVisAttributes(G4Color(0,0,1,0.8));
	logicFClad->SetVisAttributes(G4Color(0.5,0.5,0.5,0.4));
	logicSClad->SetVisAttributes(G4Color(0.,0.,0.,0.3));
	logicWorld->SetVisAttributes(G4Color(1, 1, 1, 0.2));
	logicSiPM->SetVisAttributes(G4VisAttributes::GetInvisible());
	logicElement->SetVisAttributes(G4VisAttributes::GetInvisible());
	fLogicFiber->SetVisAttributes(G4VisAttributes::GetInvisible());
	logicWindow->SetVisAttributes(G4Colour(0.,0., 0.5, 0.3));
	logicPixel->SetVisAttributes(G4Colour(0.,0.,1., 0.8));
    

	// Surface
	// Core Surface
	if(fCmdSurface){
		G4OpticalSurface* OpCoreSurface = new G4OpticalSurface("CoreSurface");
		OpCoreSurface->SetModel(glisur);
		OpCoreSurface->SetType(dielectric_dielectric);
		OpCoreSurface->SetFinish(ground);
		OpCoreSurface->SetPolish(0.985);

		G4LogicalBorderSurface* CoreSurface = new G4LogicalBorderSurface("CoreSurface", physCore, physFClad, OpCoreSurface);

		// First Cladding Surface
		G4OpticalSurface* OpFCladSurface = new G4OpticalSurface("FCladSurface");
		OpFCladSurface->SetModel(glisur);
		OpFCladSurface->SetType(dielectric_dielectric);
		OpFCladSurface->SetFinish(ground);
		OpFCladSurface->SetPolish(0.98);

		G4LogicalBorderSurface* FCladSurface = new G4LogicalBorderSurface("FCladSurface", physFClad, physSClad, OpFCladSurface);

		// Second Cladding Surface
		G4OpticalSurface* OpSCladSurface = new G4OpticalSurface("SCladSurface");
		OpSCladSurface->SetModel(glisur);
		OpSCladSurface->SetType(dielectric_dielectric);
		OpSCladSurface->SetFinish(ground);
		OpSCladSurface->SetPolish(0.5);

		G4LogicalBorderSurface* SCladSurface = new G4LogicalBorderSurface("SCladSurface", physSClad, physWorld, OpSCladSurface);
	}

	return physWorld;
}

void DetectorConstruction::SetSurface(G4bool val){
	fCmdSurface = val;
	G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void DetectorConstruction::SetWidth(G4double size){
	fFiberWidth = size;

	G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void DetectorConstruction::SetLength(G4double size){
	fFiberLength = size;

	G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void DetectorConstruction::ConstructSDandField(){
	if(!fFiber_SD.Get()){
		G4cout << "Construction /Det/FiberSD" << G4endl;
		FiberSD* fiber_SD = new FiberSD("/Det/FiberSD");
		fFiber_SD.Put(fiber_SD);
	};
	G4SDManager::GetSDMpointer()->AddNewDetector(fFiber_SD.Get());
	SetSensitiveDetector(fLogicFiber, fFiber_SD.Get());

	if(!fPixel_SD.Get()){
		G4cout << "Contruction /Det/PixelSD" << G4endl;
		PixelSD* pixel_SD = new PixelSD("Det/PixelSD", fmodel);
		pixel_SD->SetModel(fmodel);
		pixel_SD->DefineProperties();
		fPixel_SD.Put(pixel_SD);
	}

	else{
		PixelSD* pixel_SD = fPixel_SD.Get();
		pixel_SD->SetModel(fmodel);
		pixel_SD->DefineProperties();
	}

	G4SDManager::GetSDMpointer()->AddNewDetector(fPixel_SD.Get());
	SetSensitiveDetector(fLogicPixel, fPixel_SD.Get());

}

