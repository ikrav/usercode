#ifndef EWKANA_NTUPLER_TDIELECTRON_HH
#define EWKANA_NTUPLER_TDIELECTRON_HH

#include <TObject.h>

namespace mithep
{
  class TDielectron : public TObject
  {
    public:
      TDielectron(){}
      ~TDielectron(){} 
   
      Float_t mass, pt, y, phi;  // dielectron kinematics
      
      // leading electron
      Float_t pt_1, eta_1, phi_1;      
      Float_t trkIso03_1;          
      Float_t emIso03_1;           
      Float_t hadIso03_1;          
      Float_t pfIso03_1, pfIso04_1;  
      Float_t pfPx_1, pfPy_1;        
      Float_t d0_1, dz_1;            
      Float_t scEt_1, scEta_1, scPhi_1;
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
      ULong_t hltMatchBits_1;      
      Bool_t  isEcalDriven_1;      
      Bool_t  isConv_1;
            
      // lagging electron
      Float_t pt_2, eta_2, phi_2;      
      Float_t trkIso03_2;          
      Float_t emIso03_2;           
      Float_t hadIso03_2;          
      Float_t pfIso03_2, pfIso04_2;  
      Float_t pfPx_2, pfPy_2;        
      Float_t d0_2, dz_2;            
      Float_t scEt_2, scEta_2, scPhi_2;
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
      ULong_t hltMatchBits_2;      
      Bool_t  isEcalDriven_2;      
      Bool_t  isConv_2;

    ClassDef(TDielectron,1)
  };
}
#endif
