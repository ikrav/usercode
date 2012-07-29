#ifndef EWKANA_NTUPLER_TDIELECTRON_HH
#define EWKANA_NTUPLER_TDIELECTRON_HH

#include <TObject.h>

namespace mithep
{
  class TDielectron : public TObject
  {
    public:
      TDielectron():
      mass(0), pt(0), y(0), phi(0),
      pt_1(0), eta_1(0), phi_1(0), trkIso03_1(0), emIso03_1(0), hadIso03_1(0),
      chIso_00_01_1(0), chIso_01_02_1(0), chIso_02_03_1(0), chIso_03_04_1(0), chIso_04_05_1(0),
      gammaIso_00_01_1(0), gammaIso_01_02_1(0), gammaIso_02_03_1(0), gammaIso_03_04_1(0), gammaIso_04_05_1(0),
      neuHadIso_00_01_1(0), neuHadIso_01_02_1(0), neuHadIso_02_03_1(0), neuHadIso_03_04_1(0), neuHadIso_04_05_1(0),
      pfPt_1(0), pfEta_1(0), pfPhi_1(0), d0_1(0), dz_1(0), scEt_1(0), scEta_1(0), scPhi_1(0),
      ecalE_1(0), HoverE_1(0), EoverP_1(0), fBrem_1(0), deltaEtaIn_1(0), deltaPhiIn_1(0), sigiEtaiEta_1(0),
      partnerDeltaCot_1(0), partnerDist_1(0), mva_1(0), q_1(0), nExpHitsInner_1(0), scID_1(0), trkID_1(0),
      typeBits_1(0), hltMatchBits_1(0), isConv_1(0),
      pt_2(0), eta_2(0), phi_2(0), trkIso03_2(0), emIso03_2(0), hadIso03_2(0),
      chIso_00_01_2(0), chIso_01_02_2(0), chIso_02_03_2(0), chIso_03_04_2(0), chIso_04_05_2(0),
      gammaIso_00_01_2(0), gammaIso_01_02_2(0), gammaIso_02_03_2(0), gammaIso_03_04_2(0), gammaIso_04_05_2(0),
      neuHadIso_00_01_2(0), neuHadIso_01_02_2(0), neuHadIso_02_03_2(0), neuHadIso_03_04_2(0), neuHadIso_04_05_2(0),
      pfPt_2(0), pfEta_2(0), pfPhi_2(0), d0_2(0), dz_2(0), scEt_2(0), scEta_2(0), scPhi_2(0),
      ecalE_2(0), HoverE_2(0), EoverP_2(0), fBrem_2(0), deltaEtaIn_2(0), deltaPhiIn_2(0), sigiEtaiEta_2(0),
      partnerDeltaCot_2(0), partnerDist_2(0), mva_2(0), q_2(0), nExpHitsInner_2(0), scID_2(0), trkID_2(0),
      typeBits_2(0), hltMatchBits_2(0), isConv_2(0)
      {}
      ~TDielectron(){} 
   
      Float_t mass, pt, y, phi;  // dielectron kinematics
      
      // leading electron
      Float_t pt_1, eta_1, phi_1;      
      Float_t trkIso03_1;          
      Float_t emIso03_1;           
      Float_t hadIso03_1;          
      Float_t chIso_00_01_1;
      Float_t chIso_01_02_1;
      Float_t chIso_02_03_1;
      Float_t chIso_03_04_1;
      Float_t chIso_04_05_1;
      Float_t gammaIso_00_01_1;
      Float_t gammaIso_01_02_1;
      Float_t gammaIso_02_03_1;
      Float_t gammaIso_03_04_1;
      Float_t gammaIso_04_05_1;
      Float_t neuHadIso_00_01_1;
      Float_t neuHadIso_01_02_1;
      Float_t neuHadIso_02_03_1;
      Float_t neuHadIso_03_04_1;
      Float_t neuHadIso_04_05_1;
      Float_t pfPt_1, pfEta_1, pfPhi_1;        
      Float_t d0_1, dz_1;            
      Float_t scEt_1, scEta_1, scPhi_1;
      Float_t ecalE_1;
      Float_t HoverE_1;            
      Float_t EoverP_1;            
      Float_t fBrem_1;             
      Float_t deltaEtaIn_1;        
      Float_t deltaPhiIn_1;        
      Float_t sigiEtaiEta_1;       
      Float_t partnerDeltaCot_1;   
      Float_t partnerDist_1;       
      Float_t mva_1;            
      Int_t   q_1;                 
      UInt_t  nExpHitsInner_1;     
      UInt_t  scID_1;              
      UInt_t  trkID_1;        
      UInt_t  typeBits_1;     
      ULong_t hltMatchBits_1;      
      Bool_t  isConv_1;
            
      // lagging electron
      Float_t pt_2, eta_2, phi_2;      
      Float_t trkIso03_2;          
      Float_t emIso03_2;           
      Float_t hadIso03_2;          
      Float_t chIso_00_01_2;
      Float_t chIso_01_02_2;
      Float_t chIso_02_03_2;
      Float_t chIso_03_04_2;
      Float_t chIso_04_05_2;
      Float_t gammaIso_00_01_2;
      Float_t gammaIso_01_02_2;
      Float_t gammaIso_02_03_2;
      Float_t gammaIso_03_04_2;
      Float_t gammaIso_04_05_2;
      Float_t neuHadIso_00_01_2;
      Float_t neuHadIso_01_02_2;
      Float_t neuHadIso_02_03_2;
      Float_t neuHadIso_03_04_2;
      Float_t neuHadIso_04_05_2;
      Float_t pfPt_2, pfEta_2, pfPhi_2;        
      Float_t d0_2, dz_2;            
      Float_t scEt_2, scEta_2, scPhi_2;
      Float_t ecalE_2;
      Float_t HoverE_2;            
      Float_t EoverP_2;            
      Float_t fBrem_2;             
      Float_t deltaEtaIn_2;        
      Float_t deltaPhiIn_2;        
      Float_t sigiEtaiEta_2;       
      Float_t partnerDeltaCot_2;   
      Float_t partnerDist_2;       
      Float_t mva_2;            
      Int_t   q_2;                 
      UInt_t  nExpHitsInner_2;     
      UInt_t  scID_2;              
      UInt_t  trkID_2;       
      UInt_t  typeBits_2;      
      ULong_t hltMatchBits_2;    
      Bool_t  isConv_2;

    ClassDef(TDielectron,2)
  };
}
#endif
