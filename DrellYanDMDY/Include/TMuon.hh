#ifndef EWKANA_NTUPLER_TMUON_HH
#define EWKANA_NTUPLER_TMUON_HH

#include <TObject.h>

namespace mithep 
{
  class TMuon : public TObject
  {
    public:
      TMuon():
      pt(0), ptErr(0), eta(0), phi(0), staPt(0), staEta(0), staPhi(0), trkIso03(0), emIso03(0), hadIso03(0),
      chIso_00_01(0), chIso_01_02(0), chIso_02_03(0), chIso_03_04(0), chIso_04_05(0),
      gammaIso_00_01(0), gammaIso_01_02(0), gammaIso_02_03(0), gammaIso_03_04(0), gammaIso_04_05(0),
      neuHadIso_00_01(0), neuHadIso_01_02(0), neuHadIso_02_03(0), neuHadIso_03_04(0), neuHadIso_04_05(0),
      pfPt(0), pfEta(0), pfPhi(0), d0(0), dz(0), tkNchi2(0), muNchi2(0), trkKink(0), gblKink(0), q(0),
      nValidHits(0), qualityBits(0), typeBits(0), nTkHits(0), nPixHits(0), nSeg(0), nMatch(0), trkID(0), hltMatchBits(0)	     
      {}
      ~TMuon(){} 
  
      Float_t pt, ptErr, eta, phi;    // kinematics
      Float_t staPt, staEta, staPhi;  // standalone muon measurements
      Float_t trkIso03;	              // track isolation
      Float_t emIso03;	              // ECAL-based isolation
      Float_t hadIso03;	              // HCAL-based isolation
      Float_t chIso_00_01;            // Particle Flow charged isolation
      Float_t chIso_01_02;
      Float_t chIso_02_03;
      Float_t chIso_03_04;
      Float_t chIso_04_05;
      Float_t gammaIso_00_01;         // Particle Flow gamma isolation
      Float_t gammaIso_01_02;
      Float_t gammaIso_02_03;
      Float_t gammaIso_03_04;
      Float_t gammaIso_04_05;
      Float_t neuHadIso_00_01;        // Particle Flow neutral hadron isolation
      Float_t neuHadIso_01_02;
      Float_t neuHadIso_02_03;
      Float_t neuHadIso_03_04;
      Float_t neuHadIso_04_05;
      Float_t pfPt, pfEta, pfPhi;     // Matching Particle Flow candidate kinematics
      Float_t d0, dz;                 // impact parameter
      Float_t tkNchi2;	              // track chi^2/ndf 
      Float_t muNchi2;	              // global muon chi^2/ndf
      Float_t trkKink;                // kink of tracker track
      Float_t gblKink;                // kink of global track
      Int_t   q;		      // charge
      Int_t   nValidHits;	      // number of valid hits in muon system
      UInt_t  qualityBits;            // bits for various muon quality criteria
      UInt_t  typeBits;	              // global muon, tracker muon, or standalone muon
      UInt_t  nTkHits;	              // number of inner tracker hits
      UInt_t  nPixHits;	              // number of pixel hits
      UInt_t  nSeg;  	              // number of muon segments
      UInt_t  nMatch;                 // number of muon chambers matched to segments      
      UInt_t  trkID;                  // tracker track ID
      ULong_t hltMatchBits;           // bits for matching with HLT primitives 

    ClassDef(TMuon,2)
  };  
}
#endif
