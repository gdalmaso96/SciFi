/// \file  FiberHit.hh
/// \brief Definition of the FiberHit class

#ifndef FiberHit_h
#define FiberHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"

#include "TVector3.h"

#include "tls.hh"

#include <vector>

class FiberHit : public G4VHit{
	public:
		FiberHit();
		FiberHit(G4VPhysicalVolume* pVol);
		virtual ~FiberHit();

		FiberHit(const FiberHit &right);
		const FiberHit& operator=(const FiberHit &right);
		G4bool operator==(const FiberHit &right) const;

		inline void *operator new(size_t);
		inline void operator delete(void *aHit);

		virtual void Draw();
		virtual void Print();

		inline void SetEin (G4double val){fEin  = val;}
		inline void SetEdep(G4double val){fEdep = val;}
		inline void SetEdelta(G4double val){fDelta = val;}
		inline void SetSecondaryID(G4int val){fSecondaryID = val;}
		
		inline G4double GetEin (){return fEin;}
		inline G4double GetEdep(){return fEdep;}
		inline G4double GetEdelta(){return fDelta;}
		inline G4int GetSecondaryID(){return fSecondaryID;}

		inline void SetNgamma(G4int val){fNgamma = val;}
		inline G4int GetNgamma(){return fNgamma;}

		inline void SetNgammaOut(G4int val){fNgammaOut = val;}
		inline G4int GetNgammaOut(){return fNgammaOut;}

		inline void Clear(){fEin = 0; fEdep = 0; fDelta = 0; fThetaIn = 0; fNgamma = 0; fNgammaOut = 0; fPrimaryChannel = 0; fSecondaryID = -1;}
		inline const G4VPhysicalVolume* GetPhysV(){return fPhysVol;}

	private:
		G4double fEin, fEdep, fDelta, fThetaIn;
		G4int fNgamma, fNgammaOut, fPrimaryChannel, fSecondaryID;
		const G4VPhysicalVolume* fPhysVol;
};

typedef G4THitsCollection<FiberHit> FiberHitsCollection;

extern G4ThreadLocal G4Allocator<FiberHit>* FiberHitAllocator;

inline void* FiberHit::operator new(size_t){
	if(!FiberHitAllocator) FiberHitAllocator = new G4Allocator<FiberHit>;
	return (void*) FiberHitAllocator->MallocSingle();
}

inline void FiberHit::operator delete(void* aHit){
	FiberHitAllocator->FreeSingle((FiberHit*) aHit);
}

#endif


