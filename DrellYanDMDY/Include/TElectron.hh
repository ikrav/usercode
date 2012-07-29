#ifndef EWKANA_NTUPLER_TELECTRON_HH
#define EWKANA_NTUPLER_TELECTRON_HH

#include <TObject.h>

namespace mithep
{
  class TElectron : public TObject
  {
    public:
      TElectron():
      pt(0), eta(0), phi(0), trkIso03(0), emIso03(0), hadIso03(0),	     
      chIso_00_01(0), chIso_01_02(0), chIso_02_03(0), chIso_03_04(0), chIso_04_05(0),
      gammaIso_00_01(0), gammaIso_01_02(0), gammaIso_02_03(0), gammaIso_03_04(0), gammaIso_04_05(0),
      neuHadIso_00_01(0), neuHadIso_01_02(0), neuHadIso_02_03(0), neuHadIso_03_04(0), neuHadIso_04_05(0),
      pfPt(0), pfEta(0), pfPhi(0), d0(0), dz(0), scEt(0), scEta(0), scPhi(0),
      ecalE(0), HoverE(0), EoverP(0), fBrem(0), deltaEtaIn(0), deltaPhiIn(0), sigiEtaiEta(0),
      partnerDeltaCot(0), partnerDist(0), mva(0), q(0), nExpHitsInner(0), scID(0), trkID(0),
      typeBits(0), hltMatchBits(0), isConv(0)
      {}
      ~TElectron(){}
    
      Float_t pt, eta, phi;        // kinematics
      Float_t trkIso03;            // track isolation
      Float_t emIso03;             // ECAL-based isolation
      Float_t hadIso03;            // HCAL-based isolation
      Float_t chIso_00_01;         // Particle Flow charged isolation
      Float_t chIso_01_02;
      Float_t chIso_02_03;
      Float_t chIso_03_04;
      Float_t chIso_04_05;
      Float_t gammaIso_00_01;      // Particle Flow gamma isolation
      Float_t gammaIso_01_02;
      Float_t gammaIso_02_03;
      Float_t gammaIso_03_04;
      Float_t gammaIso_04_05;
      Float_t neuHadIso_00_01;     // Particle Flow neutral hadron isolation
      Float_t neuHadIso_01_02;
      Float_t neuHadIso_02_03;
      Float_t neuHadIso_03_04;
      Float_t neuHadIso_04_05;
      Float_t pfPt, pfEta, pfPhi;  // Matching Particle Flow candidate kinematics
      Float_t d0, dz;              // impact parameter
      Float_t scEt, scEta, scPhi;  // supercluster
      Float_t ecalE;               // ECAL energy
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
      UInt_t  typeBits;            // bits for electron type
      ULong_t hltMatchBits;        // bits for matching with HLT primitives
      Bool_t  isConv;              // is conversion? (vertexing method)
          
    ClassDef(TElectron,2)
  };
}
#endif
