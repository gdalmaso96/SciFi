/// \file  PixelHit.hh
/// \brief Definition of the PixelHit class

#ifndef PixelHit_h
#define PixelHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"

#include "tls.hh"

#include<vector>

class PixelHit : public G4VHit{
	public:
		PixelHit();
		virtual ~PixelHit();

		PixelHit(const PixelHit &right);
		const PixelHit& operator=(const PixelHit &right);
		G4bool operator==(const PixelHit &right) const;

		inline void *operator new(size_t);
		inline void operator delete(void *aHit);

		virtual void Draw();
		virtual void Print();

		inline void SetChannel(std::vector<G4int> val){fChannel = val;}
		inline void SetCells(std::vector<G4int> val){fCells = val;}
		inline void SetCellTime(std::vector<G4double> val){fCellTime = val;}
		inline void SetOCTFlag(std::vector<G4int> val){fOCTflag = val;}
		inline void SetDNFlag(std::vector<G4int> val){fDNflag = val;}

		inline std::vector<G4int> GetChannel(){return fChannel;}
		inline std::vector<G4int> GetCells(){return fCells;}
		inline std::vector<G4double> GetCellTime(){return fCellTime;}
		inline std::vector<G4int> GetOCTFlag(){return fOCTflag;}
		inline std::vector<G4int> GetDNFlag(){return fDNflag;}

		inline void Clear(){fChannel.clear(); fCells.clear(); fCellTime.clear(); fOCTflag.clear(); fDNflag.clear();}

	private:
		std::vector<G4int> fChannel;
		std::vector<G4int> fCells;
		std::vector<G4double> fCellTime;
		std::vector<G4int> fOCTflag;
		std::vector<G4int> fDNflag;
};

typedef G4THitsCollection<PixelHit> PixelHitsCollection;

extern G4ThreadLocal G4Allocator<PixelHit>* PixelHitAllocator;

inline void* PixelHit::operator new(size_t){
	if(!PixelHitAllocator) PixelHitAllocator = new G4Allocator<PixelHit>;
	return (void*) PixelHitAllocator->MallocSingle();
}

inline void PixelHit::operator delete(void* aHit){
	PixelHitAllocator->FreeSingle((PixelHit*) aHit);
}

#endif


