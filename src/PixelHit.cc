/// \file  PixelHit.cc
/// \brief Implementation of the PixelHit class

#include "PixelHit.hh"
#include "G4ios.hh"
#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"

G4ThreadLocal G4Allocator<PixelHit>* PixelHitAllocator = nullptr;

PixelHit::PixelHit() : fChannel(0){}

PixelHit::~PixelHit(){}

PixelHit::PixelHit(const PixelHit &right) : G4VHit(){
	fChannel = right.fChannel; // active channel
	fCells = right.fCells; // vector of active cells
	fCellTime = right.fCellTime; // vector of arrive time in cell
	fOCTflag = right.fOCTflag;
	fDNflag = right.fDNflag;
}

const PixelHit& PixelHit::operator=(const PixelHit &right){
	fChannel = right.fChannel; // active channel
	fCells = right.fCells; // vector of active cells
	fCellTime = right.fCellTime; // vector of arrive time in cell
	fOCTflag = right.fOCTflag;
	fDNflag = right.fDNflag;
	return* this;
}


G4bool PixelHit::operator==(const PixelHit &) const{
	return false;
}

void PixelHit::Draw(){}

void PixelHit::Print(){}
