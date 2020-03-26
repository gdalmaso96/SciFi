/// \file  ScintHit.hh
/// \brief Definition of the ScintHit class

#ifndef ScintHit_h
#define ScintHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"

#include "TVector3.h"

#include "tls.hh"

#include <vector>

class ScintHit : public G4VHit{
	public:
		ScintHit();
		ScintHit(G4VPhysicalVolume* pVol);
		virtual ~ScintHit();

		ScintHit(const ScintHit &right);
		const ScintHit& operator=(const ScintHit &right);
		G4bool operator==(const ScintHit &right) const;

		inline void *operator new(size_t);
		inline void operator delete(void *aHit);

		virtual void Draw();
		virtual void Print();

		inline void SetEin (G4double val){fEin  = val;}
		inline void SetEdep(G4double val){fEdep = val;}
		inline void SetEdelta(G4double val){fDelta = val;}
		
		inline G4double GetEin (){return fEin;}
		inline G4double GetEdep(){return fEdep;}
		inline G4double GetEdelta(){return fDelta;}

		inline void SetNgammaSec(G4int ngammasec){fNgammaSec = ngammasec;}
		inline G4int GetNgammaSec(){return fNgammaSec;}


		inline void Clear(){fEin = 0; fEdep = 0; fDelta = 0; fNgammaSec = 0;}
		inline const G4VPhysicalVolume* GetPhysV(){return fPhysVol;}

	private:
		G4double fEin, fEdep, fDelta;
		G4int fNgammaSec;
		const G4VPhysicalVolume* fPhysVol;
};

typedef G4THitsCollection<ScintHit> ScintHitsCollection;

extern G4ThreadLocal G4Allocator<ScintHit>* ScintHitAllocator;

inline void* ScintHit::operator new(size_t){
	if(!ScintHitAllocator) ScintHitAllocator = new G4Allocator<ScintHit>;
	return (void*) ScintHitAllocator->MallocSingle();
}

inline void ScintHit::operator delete(void* aHit){
	ScintHitAllocator->FreeSingle((ScintHit*) aHit);
}

#endif


