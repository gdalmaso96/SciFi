/// \file  FiberHit.cc
/// \brief Implementation of the FiberHit class

#include "FiberHit.hh"
#include "G4ios.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"

G4ThreadLocal G4Allocator<FiberHit>* FiberHitAllocator = nullptr;

FiberHit::FiberHit() : fEin(0.), fEdep(0.), fDelta(0), fNgamma(0), fNgammaOut(0), fPhysVol(nullptr){}

FiberHit::FiberHit(G4VPhysicalVolume* pVol) : fPhysVol(pVol){}

FiberHit::~FiberHit(){}

FiberHit::FiberHit(const FiberHit &right) : G4VHit(){
	fEin  = right.fEin;
	fEdep = right.fEdep;
	fDelta = right.fDelta;
	fNgamma = right.fNgamma;
	fNgammaOut = right.fNgammaOut;
	fPhysVol = right.fPhysVol;
}

const FiberHit& FiberHit::operator=(const FiberHit &right){
	fEin  = right.fEin;
	fEdep = right.fEdep;
	fDelta = right.fDelta;
	fNgamma = right.fNgamma;
	fNgammaOut = right.fNgammaOut;
	return* this;
}

G4bool FiberHit::operator==(const FiberHit &) const{
	return false;
}

void FiberHit::Draw(){}

void FiberHit::Print(){}


