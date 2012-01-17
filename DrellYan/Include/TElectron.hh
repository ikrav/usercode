#ifndef EWKANA_NTUPLER_TELECTRON_HH
#define EWKANA_NTUPLER_TELECTRON_HH

#include <TObject.h>

namespace mithep
{
  class TElectron : public TObject
  {
    public:
      TElectron(){}
      ~TElectron(){}
    
      Float_t pt, eta, phi;        // kinematics
      Float_t trkIso03;            // track isolation
      Float_t emIso03;             // ECAL-based isolation
      Float_t hadIso03;            // HCAL-based isolation
      Float_t pfIso03, pfIso04;    // Particle Flow isolation
      Float_t pfPx, pfPy;          // Matching Particle Flow candidate (px,py)
      Float_t d0, dz;              // impact parameter
      Float_t scEt, scEta, scPhi;  // supercluster
      Float_t HoverE;              // H / E
      Float_t EoverP;              // E / p
      Float_t fBrem;               // brem fraction  
      Float_t deltaEtaIn;          // eta difference between track (at vertex) and SC
      Float_t deltaPhiIn;          // phi difference between track (at vertex) and SC
      Float_t sigiEtaiEta;         // eta-width of shower in number of crystals
      Float_t partnerDeltaCot;     // cot(theta) difference with conversion partner track	
      Float_t partnerDist;         // distance in x-y plane to nearest conversion partner track
      Float_t mva;                 // MVA electron ID
      Int_t   q;                   // charge
      UInt_t  nExpHitsInner;       // number of hits expected before first hit      	             
      UInt_t  scID;                // supercluster ID (for matching to photon superclusters)
      UInt_t  trkID;               // tracker track ID (for matching to muons)
      ULong_t hltMatchBits;        // bits for matching with HLT primitives
      Bool_t  isEcalDriven;        // is ECAL seeded electron?
      Bool_t  isConv;              // is conversion? (vertexing method)
          
    ClassDef(TElectron,2)
  };
}
#endif
