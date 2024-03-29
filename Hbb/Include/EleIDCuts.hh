#ifndef ELEIDCUTS_HH
#define ELEIDCUTS_HH

#include <TMath.h>
#include "TElectron.hh"
#include "TDielectron.hh"
#include <cassert>

/*
 * [0] WP95
 * [1] WP90
 * [2] WP85
 * [3] WP80
 * [4] WP70
 * [5] WP60
 *
 */
Double_t _MissingHits[6]     = {     1,     1,     1,     0,     0,     0 };
Double_t _Dist[6]            = {     0,  0.02,  0.02,  0.02,  0.02,  0.02 };
Double_t _DCot[6]            = {     0,  0.02,  0.02,  0.02,  0.02,  0.02 };

Double_t _CombIsoEB[6]       = {  0.15,  0.10,  0.09,  0.07,  0.04,  0.03 };
Double_t _TrkIsoEB[6]        = {  0.15,  0.12,  0.09,  0.09,  0.05,  0.04 };
Double_t _EcalIsoEB[6]       = {  2.00,  0.09,  0.08,  0.07,  0.06,  0.04 };
Double_t _HcalIsoEB[6]       = {  0.12,  0.10,  0.10,  0.10,  0.03,  0.03 };
Double_t _SigmaiEtaiEtaEB[6] = {  0.01,  0.01,  0.01,  0.01,  0.01,  0.01 };
Double_t _DPhiEB[6]          = {  0.80,  0.80,  0.06,  0.06,  0.03, 0.025 };
Double_t _DEtaEB[6]          = { 0.007, 0.007, 0.006, 0.004, 0.004, 0.004 };
Double_t _HoEEB[6]           = {  0.15,  0.12,  0.04,  0.04, 0.025, 0.025 };

Double_t _CombIsoEE[6]       = {   0.1,  0.07,  0.06,  0.06,  0.03,  0.02 };
Double_t _TrkIsoEE[6]        = {  0.08,  0.05,  0.05,  0.04, 0.025, 0.025 };
Double_t _EcalIsoEE[6]       = {  0.06,  0.06,  0.05,  0.05, 0.025,  0.02 };
Double_t _HcalIsoEE[6]       = {  0.05,  0.03, 0.025, 0.025,  0.02,  0.02 };
Double_t _SigmaiEtaiEtaEE[6] = {  0.03,  0.03,  0.03,  0.03,  0.03,  0.03 };
Double_t _DPhiEE[6]          = {   0.7,   0.7,  0.04,  0.03,  0.02,  0.02 };
Double_t _DEtaEE[6]          = {  0.01, 0.009, 0.007, 0.007, 0.005, 0.005 };
Double_t _HoEEE[6]           = {  0.07,  0.05, 0.025, 0.025, 0.025, 0.025 };


Bool_t passWP(const mithep::TDielectron *dielectron, const Bool_t useCombIso, const Int_t iWP);
Bool_t passWP(const mithep::TElectron *electron, const Bool_t useCombIso, const Int_t iWP);

Bool_t passWP95(const mithep::TDielectron *dielectron, const Bool_t useCombIso=kFALSE) { return passWP(dielectron, useCombIso, 0); }
Bool_t passWP95(const mithep::TElectron   *electron,   const Bool_t useCombIso=kFALSE) { return passWP(electron,   useCombIso, 0); }
Bool_t passWP90(const mithep::TDielectron *dielectron, const Bool_t useCombIso=kFALSE) { return passWP(dielectron, useCombIso, 1); }
Bool_t passWP90(const mithep::TElectron   *electron,   const Bool_t useCombIso=kFALSE) { return passWP(electron,   useCombIso, 1); }
Bool_t passWP85(const mithep::TDielectron *dielectron, const Bool_t useCombIso=kFALSE) { return passWP(dielectron, useCombIso, 2); }
Bool_t passWP85(const mithep::TElectron   *electron,   const Bool_t useCombIso=kFALSE) { return passWP(electron,   useCombIso, 2); }
Bool_t passWP80(const mithep::TDielectron *dielectron, const Bool_t useCombIso=kFALSE) { return passWP(dielectron, useCombIso, 3); }
Bool_t passWP80(const mithep::TElectron   *electron,   const Bool_t useCombIso=kFALSE) { return passWP(electron,   useCombIso, 3); }
Bool_t passWP70(const mithep::TDielectron *dielectron, const Bool_t useCombIso=kFALSE) { return passWP(dielectron, useCombIso, 4); }
Bool_t passWP70(const mithep::TElectron   *electron,   const Bool_t useCombIso=kFALSE) { return passWP(electron,   useCombIso, 4); }
Bool_t passWP60(const mithep::TDielectron *dielectron, const Bool_t useCombIso=kFALSE) { return passWP(dielectron, useCombIso, 5); }
Bool_t passWP60(const mithep::TElectron   *electron,   const Bool_t useCombIso=kFALSE) { return passWP(electron,   useCombIso, 5); }

// 2011 cuts
Bool_t passSmurf(const mithep::TDielectron *dielectron);
Bool_t passSmurf(const mithep::TElectron   *electron);

Bool_t passWP95ID2011(const mithep::TDielectron *dielectron);
Bool_t passWP95ID2011(const mithep::TElectron   *electron);

// Implementation
Bool_t passWP(const mithep::TDielectron *dielectron, const Bool_t useCombIso, const Int_t iWP)
{
  assert(dielectron);
  
  const Double_t kGAP_LOW  = 1.4442;
  const Double_t kGAP_HIGH = 1.566;

  if((fabs(dielectron->scEta_1)>kGAP_LOW) && (fabs(dielectron->scEta_1)<kGAP_HIGH)) return kFALSE;
  if((fabs(dielectron->scEta_2)>kGAP_LOW) && (fabs(dielectron->scEta_2)<kGAP_HIGH)) return kFALSE;

  const Bool_t isB1 = (fabs(dielectron->scEta_1)<kGAP_LOW);
  const Bool_t isB2 = (fabs(dielectron->scEta_2)<kGAP_LOW);
  
  // conversion rejection
  if(dielectron->nExpHitsInner_1 > _MissingHits[iWP]) return kFALSE;
  if(dielectron->nExpHitsInner_2 > _MissingHits[iWP]) return kFALSE;
  if((fabs(dielectron->partnerDist_1) < _Dist[iWP]) && (fabs(dielectron->partnerDeltaCot_1) < _DCot[iWP])) return kFALSE;
  if((fabs(dielectron->partnerDist_2) < _Dist[iWP]) && (fabs(dielectron->partnerDeltaCot_2) < _DCot[iWP])) return kFALSE;	  
  	  
  // barrel/endcap dependent requirments      
  if(isB1) {  // barrel
    if(useCombIso) {
      Double_t iso = (dielectron->trkIso03_1 + TMath::Max(dielectron->emIso03_1-1,Float_t(0)) + dielectron->hadIso03_1)/dielectron->pt_1;
      if(iso > _CombIsoEB[iWP]) return kFALSE;
    } else {
      if(dielectron->trkIso03_1	> _TrkIsoEB[iWP]*(dielectron->pt_1))  return kFALSE;
      if(dielectron->emIso03_1	> _EcalIsoEB[iWP]*(dielectron->pt_1)) return kFALSE;
      if(dielectron->hadIso03_1	> _HcalIsoEB[iWP]*(dielectron->pt_1)) return kFALSE;
    }
    if(dielectron->sigiEtaiEta_1      > _SigmaiEtaiEtaEB[iWP]) return kFALSE;
    if(fabs(dielectron->deltaPhiIn_1) > _DPhiEB[iWP])	       return kFALSE;
    if(fabs(dielectron->deltaEtaIn_1) > _DEtaEB[iWP])  	       return kFALSE;
    if(dielectron->HoverE_1	      > _HoEEB[iWP])	       return kFALSE;
  
  } else {  // endcap
    if(useCombIso) {
      Double_t iso = (dielectron->trkIso03_1 + dielectron->emIso03_1 + dielectron->hadIso03_1)/dielectron->pt_1;
      if(iso > _CombIsoEE[iWP]) return kFALSE;
    } else {
      if(dielectron->trkIso03_1	> _TrkIsoEE[iWP]*(dielectron->pt_1))  return kFALSE;
      if(dielectron->emIso03_1	> _EcalIsoEE[iWP]*(dielectron->pt_1)) return kFALSE;
      if(dielectron->hadIso03_1	> _HcalIsoEE[iWP]*(dielectron->pt_1)) return kFALSE;
    }
    if(dielectron->sigiEtaiEta_1      > _SigmaiEtaiEtaEE[iWP]) return kFALSE;
    if(fabs(dielectron->deltaPhiIn_1) > _DPhiEE[iWP])	       return kFALSE;
    if(fabs(dielectron->deltaEtaIn_1) > _DEtaEE[iWP])  	       return kFALSE;
    if(dielectron->HoverE_1	      > _HoEEE[iWP])  	       return kFALSE;
  }

  if(isB2) {  // barrel
    if(useCombIso) {
      Double_t iso = (dielectron->trkIso03_2 + TMath::Max(dielectron->emIso03_2-1,Float_t(0)) + dielectron->hadIso03_2)/dielectron->pt_2;
      if(iso > _CombIsoEB[iWP]) return kFALSE;
    } else {
      if(dielectron->trkIso03_2	> _TrkIsoEB[iWP]*(dielectron->pt_2))  return kFALSE;
      if(dielectron->emIso03_2	> _EcalIsoEB[iWP]*(dielectron->pt_2)) return kFALSE;
      if(dielectron->hadIso03_2	> _HcalIsoEB[iWP]*(dielectron->pt_2)) return kFALSE;
    }
    if(dielectron->sigiEtaiEta_2      > _SigmaiEtaiEtaEB[iWP]) return kFALSE;
    if(fabs(dielectron->deltaPhiIn_2) > _DPhiEB[iWP])	       return kFALSE;
    if(fabs(dielectron->deltaEtaIn_2) > _DEtaEB[iWP])  	       return kFALSE;
    if(dielectron->HoverE_2	      > _HoEEB[iWP])	       return kFALSE;
  
  } else {  // endcap
    if(useCombIso) {
      Double_t iso = (dielectron->trkIso03_2 + dielectron->emIso03_2 + dielectron->hadIso03_2)/dielectron->pt_2;
      if(iso > _CombIsoEE[iWP]) return kFALSE;
    } else {
      if(dielectron->trkIso03_2	> _TrkIsoEE[iWP]*(dielectron->pt_2))  return kFALSE;
      if(dielectron->emIso03_2	> _EcalIsoEE[iWP]*(dielectron->pt_2)) return kFALSE;
      if(dielectron->hadIso03_2	> _HcalIsoEE[iWP]*(dielectron->pt_2)) return kFALSE;
    }
    if(dielectron->sigiEtaiEta_2      > _SigmaiEtaiEtaEE[iWP]) return kFALSE;
    if(fabs(dielectron->deltaPhiIn_2) > _DPhiEE[iWP])	       return kFALSE;
    if(fabs(dielectron->deltaEtaIn_2) > _DEtaEE[iWP])  	       return kFALSE;
    if(dielectron->HoverE_2	      > _HoEEE[iWP])  	       return kFALSE;
  }
  
  return kTRUE;
}

Bool_t passWP(const mithep::TElectron *electron, const Bool_t useCombIso, const Int_t iWP)
{
  assert(electron);
  
  const Double_t kGAP_LOW  = 1.4442;
  const Double_t kGAP_HIGH = 1.566;

  if((fabs(electron->scEta)>kGAP_LOW) && (fabs(electron->scEta)<kGAP_HIGH)) return kFALSE;
  
  // conversion rejection
  if(electron->nExpHitsInner > _MissingHits[iWP]) return kFALSE;
  if((fabs(electron->partnerDist) < _Dist[iWP]) && (fabs(electron->partnerDeltaCot) < _DCot[iWP])) return kFALSE;
  	  
  // barrel/endcap dependent requirments      
  if(fabs(electron->scEta)<kGAP_LOW) {  // barrel
    if(useCombIso) {
      Double_t iso = (electron->trkIso03 + TMath::Max(electron->emIso03-1,Float_t(0)) + electron->hadIso03)/electron->pt;
      if(iso > _CombIsoEB[iWP]) return kFALSE;
    } else {
      if(electron->trkIso03 > _TrkIsoEB[iWP]*(electron->pt))  return kFALSE;
      if(electron->emIso03  > _EcalIsoEB[iWP]*(electron->pt)) return kFALSE;
      if(electron->hadIso03 > _HcalIsoEB[iWP]*(electron->pt)) return kFALSE;
    }
    if(electron->sigiEtaiEta      > _SigmaiEtaiEtaEB[iWP]) return kFALSE;
    if(fabs(electron->deltaPhiIn) > _DPhiEB[iWP])	   return kFALSE;
    if(fabs(electron->deltaEtaIn) > _DEtaEB[iWP])  	   return kFALSE;
    if(electron->HoverE	          > _HoEEB[iWP])	   return kFALSE;
  
  } else {  // endcap
    if(useCombIso) {
      Double_t iso = (electron->trkIso03 + electron->emIso03 + electron->hadIso03)/electron->pt;
      if(iso > _CombIsoEE[iWP]) return kFALSE;
    } else {
      if(electron->trkIso03 > _TrkIsoEE[iWP]*(electron->pt))  return kFALSE;
      if(electron->emIso03  > _EcalIsoEE[iWP]*(electron->pt)) return kFALSE;
      if(electron->hadIso03 > _HcalIsoEE[iWP]*(electron->pt)) return kFALSE;
    }
    if(electron->sigiEtaiEta      > _SigmaiEtaiEtaEE[iWP]) return kFALSE;
    if(fabs(electron->deltaPhiIn) > _DPhiEE[iWP])	   return kFALSE;
    if(fabs(electron->deltaEtaIn) > _DEtaEE[iWP])  	   return kFALSE;
    if(electron->HoverE	          > _HoEEE[iWP])  	   return kFALSE;
  }
  
  return kTRUE;
}

Bool_t passSmurf(const mithep::TDielectron *dielectron)
{
  if(fabs(dielectron->d0_1) > 0.02) return kFALSE;
  if(fabs(dielectron->d0_2) > 0.02) return kFALSE;
  if(fabs(dielectron->dz_1) > 0.1)  return kFALSE;  
  if(fabs(dielectron->dz_2) > 0.1)  return kFALSE;
  
  // conversion rejection
  if(dielectron->nExpHitsInner_1 > 0) return kFALSE;
  if(dielectron->nExpHitsInner_2 > 0) return kFALSE;
  if(dielectron->isConv_1)            return kFALSE;  
  if(dielectron->isConv_2)            return kFALSE;
  
  // barrel/endcap dependent requirments      
  if(fabs(dielectron->scEta_1)<1.479) {
    // barrel
    if(dielectron->pfIso04_1 > 0.13*(dielectron->pt_1)) return kFALSE;
    
    if(dielectron->pt_1>20) {
      if(dielectron->sigiEtaiEta_1	> 0.01)  return kFALSE;
      if(fabs(dielectron->deltaPhiIn_1) > 0.06)  return kFALSE;
      if(fabs(dielectron->deltaEtaIn_1) > 0.004) return kFALSE;
      if(dielectron->HoverE_1	        > 0.04)  return kFALSE;
      
    } else {
      if(dielectron->sigiEtaiEta_1	> 0.01)  return kFALSE;
      if(fabs(dielectron->deltaPhiIn_1) > 0.03)  return kFALSE;
      if(fabs(dielectron->deltaEtaIn_1) > 0.004) return kFALSE;
      if(dielectron->HoverE_1	        > 0.025) return kFALSE;    
    }
    
  } else {
    // endcap
    if(dielectron->pfIso04_1 > 0.09*(dielectron->pt_1)) return kFALSE;
    
    if(dielectron->pt_1>20) {
      if(dielectron->sigiEtaiEta_1	> 0.03)  return kFALSE;
      if(fabs(dielectron->deltaPhiIn_1) > 0.03)  return kFALSE;
      if(fabs(dielectron->deltaEtaIn_1) > 0.007) return kFALSE;
      if(dielectron->HoverE_1	        > 0.10)  return kFALSE;
      
    } else {
      if(dielectron->sigiEtaiEta_1	> 0.03)  return kFALSE;
      if(fabs(dielectron->deltaPhiIn_1) > 0.02)  return kFALSE;
      if(fabs(dielectron->deltaEtaIn_1) > 0.005) return kFALSE;
      if(dielectron->HoverE_1	        > 0.10)  return kFALSE;      
    }
  }
  
  if(dielectron->pt_1 < 20 &&
     !((dielectron->fBrem_1>0.15) || (fabs(dielectron->scEta_1)<1 && dielectron->EoverP_1>0.95)))
    return kFALSE;
  
  if(fabs(dielectron->scEta_2)<1.479) {
    // barrel
    if(dielectron->pfIso04_2 > 0.13*(dielectron->pt_2)) return kFALSE;
    
    if(dielectron->pt_2>20) {
      if(dielectron->sigiEtaiEta_2	> 0.01)  return kFALSE;
      if(fabs(dielectron->deltaPhiIn_2) > 0.06)  return kFALSE;
      if(fabs(dielectron->deltaEtaIn_2) > 0.004) return kFALSE;
      if(dielectron->HoverE_2	        > 0.04)  return kFALSE;
      
    } else {
      if(dielectron->sigiEtaiEta_2	> 0.01)  return kFALSE;
      if(fabs(dielectron->deltaPhiIn_2) > 0.03)  return kFALSE;
      if(fabs(dielectron->deltaEtaIn_2) > 0.004) return kFALSE;
      if(dielectron->HoverE_2	        > 0.025) return kFALSE;    
    }
    
  } else {
    // endcap
    if(dielectron->pfIso04_2 > 0.09*(dielectron->pt_2)) return kFALSE;
    
    if(dielectron->pt_2>20) {
      if(dielectron->sigiEtaiEta_2	> 0.03)  return kFALSE;
      if(fabs(dielectron->deltaPhiIn_2) > 0.03)  return kFALSE;
      if(fabs(dielectron->deltaEtaIn_2) > 0.007) return kFALSE;
      if(dielectron->HoverE_2	        > 0.10)  return kFALSE;
      
    } else {
      if(dielectron->sigiEtaiEta_2	> 0.03)  return kFALSE;
      if(fabs(dielectron->deltaPhiIn_2) > 0.02)  return kFALSE;
      if(fabs(dielectron->deltaEtaIn_2) > 0.005) return kFALSE;
      if(dielectron->HoverE_2	        > 0.10)  return kFALSE;      
    }
  }
  
  if(dielectron->pt_2 < 20 && 
     !((dielectron->fBrem_2>0.15) || (fabs(dielectron->scEta_2)<1 && dielectron->EoverP_2>0.95)))
    return kFALSE;
  
  return kTRUE;
}

Bool_t passSmurf(const mithep::TElectron *electron)
{
  if(fabs(electron->d0) > 0.02) return kFALSE;
  if(fabs(electron->dz) > 0.1)  return kFALSE;
  
  // conversion rejection
  if(electron->nExpHitsInner > 0) return kFALSE;
  if(electron->isConv)            return kFALSE;
     
  // barrel/endcap dependent requirments      
  if(fabs(electron->scEta)<1.479) {
    // barrel
    if(electron->pfIso04 > 0.13*(electron->pt)) return kFALSE;
     
    if(electron->pt>20) {
      if(electron->sigiEtaiEta	    > 0.01)  return kFALSE;
      if(fabs(electron->deltaPhiIn) > 0.06)  return kFALSE;
      if(fabs(electron->deltaEtaIn) > 0.004) return kFALSE;
      if(electron->HoverE	    > 0.04)  return kFALSE;
    
    } else {
      if(electron->sigiEtaiEta	    > 0.01)  return kFALSE;
      if(fabs(electron->deltaPhiIn) > 0.03)  return kFALSE;
      if(fabs(electron->deltaEtaIn) > 0.004) return kFALSE;
      if(electron->HoverE	    > 0.025) return kFALSE;    
    }
  
  } else {
    // endcap
    if(electron->pfIso04 > 0.09*(electron->pt)) return kFALSE;
     
    if(electron->pt>20) {
      if(electron->sigiEtaiEta	    > 0.03)  return kFALSE;
      if(fabs(electron->deltaPhiIn) > 0.03)  return kFALSE;
      if(fabs(electron->deltaEtaIn) > 0.007) return kFALSE;
      if(electron->HoverE	    > 0.10)  return kFALSE;
    
    } else {
      if(electron->sigiEtaiEta	    > 0.03)  return kFALSE;
      if(fabs(electron->deltaPhiIn) > 0.02)  return kFALSE;
      if(fabs(electron->deltaEtaIn) > 0.005) return kFALSE;
      if(electron->HoverE	    > 0.10)  return kFALSE;      
    }
  }
  
  if(electron->pt < 20)
    return ((electron->fBrem>0.15) || (fabs(electron->scEta)<1 && electron->EoverP>0.95));
  
  return kTRUE;
}

Bool_t passWP95ID2011(const mithep::TElectron   *electron){

  // barrel/endcap dependent requirments      
  if(fabs(electron->scEta)<1.479) {
    // barrel
    if(electron->sigiEtaiEta	  > 0.01)  return kFALSE;
    if(fabs(electron->deltaPhiIn) > 0.8 )  return kFALSE;
    if(fabs(electron->deltaEtaIn) > 0.007) return kFALSE;
    if(electron->HoverE	          > 0.15) return kFALSE;    
  } else {
    // endcap
    if(electron->sigiEtaiEta	  > 0.03)  return kFALSE;
    if(fabs(electron->deltaPhiIn) > 0.7)  return kFALSE;
    if(fabs(electron->deltaEtaIn) > 0.01) return kFALSE;
    if(electron->HoverE	          > 0.15)  return kFALSE;
    
  }
  return kTRUE;  
}

Bool_t passWP95ID2011(const mithep::TDielectron *dielectron){

  // barrel/endcap dependent requirments      
  if(fabs(dielectron->scEta_1)<1.479) {
    // barrel
    if(dielectron->sigiEtaiEta_1	> 0.01)  return kFALSE;
    if(fabs(dielectron->deltaPhiIn_1)   > 0.8 )  return kFALSE;
    if(fabs(dielectron->deltaEtaIn_1)   > 0.007) return kFALSE;
    if(dielectron->HoverE_1	        > 0.15)  return kFALSE;
  } else {
    // endcap
    if(dielectron->sigiEtaiEta_1	> 0.03)  return kFALSE;
    if(fabs(dielectron->deltaPhiIn_1)   > 0.7 )  return kFALSE;
    if(fabs(dielectron->deltaEtaIn_1)   > 0.01) return kFALSE;
    if(dielectron->HoverE_1	        > 0.15)  return kFALSE;      
  }
  
  if(fabs(dielectron->scEta_2)<1.479) {
    // barrel
    if(dielectron->sigiEtaiEta_2	> 0.01)  return kFALSE;
    if(fabs(dielectron->deltaPhiIn_2)   > 0.8 )  return kFALSE;
    if(fabs(dielectron->deltaEtaIn_2)   > 0.007) return kFALSE;
    if(dielectron->HoverE_2	        > 0.15 ) return kFALSE;        
  } else {
    // endcap
    if(dielectron->sigiEtaiEta_2	> 0.03)  return kFALSE;
    if(fabs(dielectron->deltaPhiIn_2)   > 0.7 )  return kFALSE;
    if(fabs(dielectron->deltaEtaIn_2)   > 0.01) return kFALSE;
    if(dielectron->HoverE_2	        > 0.15)  return kFALSE;      
  }
  
  return kTRUE;
}

#endif
