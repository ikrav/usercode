#ifndef ELEIDCUTS_HH
#define ELEIDCUTS_HH

#include <TMath.h>
#include "../Include/TElectron.hh"
#include "../Include/TDielectron.hh"
#include "../Include/TEventInfo.hh"
#include <cassert>
#include <iostream>

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

Bool_t passSmurf(const mithep::TDielectron *dielectron);
Bool_t passSmurf(const mithep::TElectron   *electron);

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

  printf("smurfID is deprecated in this code version\n");
  assert(0);

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
    // pfIso04 is not in the ntuple now!
//     if(dielectron->pfIso04_1 > 0.13*(dielectron->pt_1)) return kFALSE;
    
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
    // pfIso04 is not in the ntuple now!
//     if(dielectron->pfIso04_1 > 0.09*(dielectron->pt_1)) return kFALSE;
    
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
    // pfIso04 is not in the ntuple now!
//     if(dielectron->pfIso04_2 > 0.13*(dielectron->pt_2)) return kFALSE;
    
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
    // pfIso04 is not in the ntuple now!
//     if(dielectron->pfIso04_2 > 0.09*(dielectron->pt_2)) return kFALSE;
    
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

  printf("smurfID is deprecated in this code version\n");
  assert(0);
  if(fabs(electron->d0) > 0.02) return kFALSE;
  if(fabs(electron->dz) > 0.1)  return kFALSE;
  
  // conversion rejection
  if(electron->nExpHitsInner > 0) return kFALSE;
  if(electron->isConv)            return kFALSE;
     
  // barrel/endcap dependent requirments      
  if(fabs(electron->scEta)<1.479) {
    // barrel
    // pfIso04 is not in the ntuple now!
//     if(electron->pfIso04 > 0.13*(electron->pt)) return kFALSE;
     
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
    // pfIso04 is not in the ntuple now!
//     if(electron->pfIso04 > 0.09*(electron->pt)) return kFALSE;
     
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

//
// The cuts defined below are recommended by EGM in spring 2012 for
// analyses based on 2011 data. Taken from:
//    https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaCutBasedIdentification

typedef enum {
  WP_VETO = 0,  // similar to VBTF 95
  WP_LOOSE = 1, // similar to VBTF 90
  WP_MEDIUM = 2,// similar to VBTF 80
  WP_TIGHT = 3  // similar to VBTF 70
} WorkingPointType;

Double_t _EGM2011_DEtaEB[4]          = {0.007,  0.007,  0.004,  0.004};
Double_t _EGM2011_DPhiEB[4]          = {0.8,    0.15,   0.06,   0.03};
Double_t _EGM2011_SigmaIetaIetaEB[4] = {0.01,   0.01,   0.01,   0.01};
Double_t _EGM2011_HoEEB[4]           = {0.15,   0.12,   0.12,   0.12};
Double_t _EGM2011_RelPFIsoEB[4]      = {0.15,   0.15,   0.15,   0.10};

Double_t _EGM2011_DEtaEE[4]          = {0.01,   0.009,   0.007,   0.005};
Double_t _EGM2011_DPhiEE[4]          = {0.7,    0.10,    0.03,    0.02};
Double_t _EGM2011_SigmaIetaIetaEE[4] = {0.03,   0.03,    0.03,    0.03};
Double_t _EGM2011_HoEEE[4]           = {-1,     0.10,    0.10,    0.10}; // -1 means cut not applied
Double_t _EGM2011_RelPFIsoEEHighPt[4]= {0.15,   0.15,    0.15,    0.10}; // For endcaps the cut is different for pt>20 and pt<20
Double_t _EGM2011_RelPFIsoEELowPt[4] = {0.15,   0.10,    0.10,    0.07};

Double_t _EGM2011_InvEMinusInvP[4]   = {-1, 0.05, 0.05, 0.05}; // This is cut value for fabs(1/E - 1/p). -1 means it is not applied

Double_t _EGM2011_D0Vtx[4]           = {0.04,   0.02,    0.02,    0.02};
Double_t _EGM2011_DZVtx[4]           = {0.2,    0.2,     0.1,     0.1};

Double_t _EGM2011_MissingHits[4]     = {    -1,     1,     1,     0}; // -1 means cut is not applied
Double_t _EGM2011_PassConvVFitCut[4] = {    -1,     1,     1,     1}; // -1 means N/A, 1 means to apply the cut

// The effective area constants for rho correction for pile-up.
// The constants are recommended by EGM for use on 2011 and 2012 
// data:
//     https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaEARhoCorrection
// This is a single effective area that covers photons and neutral hadrons.
// estimated on 2011 data for dR<0.3.
const int _nEtaBinsForAeff = 7;
Double_t _etaLimitsForAeff[_nEtaBinsForAeff+1] = {0.0, 1.0, 1.479, 2.0, 2.2, 2.3, 2.4, 100.0};
Double_t _AeffDR03[_nEtaBinsForAeff] = {0.10, 0.12, 0.085, 0.11, 0.12, 0.12, 0.13};


Bool_t passEGM2011(const mithep::TDielectron *dielectron, WorkingPointType wp, double rho);
Bool_t passEGM2011(const mithep::TElectron   *electron, WorkingPointType wp, double rho);

Bool_t passEGM2011(const mithep::TElectron *electron, WorkingPointType wp, double rho)
{
  if(fabs(electron->d0) > _EGM2011_D0Vtx[wp] ) return kFALSE;
  if(fabs(electron->dz) > _EGM2011_DZVtx[wp] )  return kFALSE;
  
  // conversion rejection
  if(_EGM2011_MissingHits[wp] != -1      // this means that the cut IS to be applied
     && electron->nExpHitsInner > _EGM2011_MissingHits[wp]) return kFALSE;
  if(_EGM2011_PassConvVFitCut[wp] != -1  // this means that the cut IS to be applied
     && electron->isConv)            return kFALSE;

  // Cut on fabs(1/E - 1/p)
//    double theta = 2*atan(exp( - electron->scEta ) );
//    double scE = fabs( electron->scEt / sin( theta ) ); //fabs is probably not really needed, but just in case
//    // NOTE: THIS IS NOT TRULY PROPER IMPLEMENTATION OF THE CUT,
//    // the energy should be not SC, but of GsfElectron
//    double invEMinusInvP = (1/scE)*fabs( 1 - electron->EoverP );
   double invEMinusInvP = (1/electron->ecalE)*fabs( 1 - electron->EoverP );
//    printf("  old %f      new %f   ecalE %f \n", invEMinusInvP, invEMinusInvP1, electron->ecalE);
  if( _EGM2011_InvEMinusInvP[wp] != -1   // this means that the cut IS to be applied
      && invEMinusInvP > _EGM2011_InvEMinusInvP[wp] ) return kFALSE;

  // Compute PF isolation:
  // Add energy in the rings with dR 0.0-0.1, 0.1-0.2, 0.2-0.3 
  // to get the total deposits in the cone 0.3.

  double chIso = electron->chIso_00_01 + electron->chIso_01_02 
    + electron->chIso_02_03;
  double gammaIso = electron->gammaIso_00_01 + electron->gammaIso_01_02 
    + electron->gammaIso_02_03;
  double neuHadIso = electron->neuHadIso_00_01 + electron->neuHadIso_01_02 
    + electron->neuHadIso_02_03;
  double Aeff = 0;
  for(int i=0; i<_nEtaBinsForAeff; i++){
    if ( fabs(electron->scEta) >= _etaLimitsForAeff[i] 
	 && fabs(electron->scEta) < _etaLimitsForAeff[i+1]){
      Aeff = _AeffDR03[i];
      break;
    }
  }
  
  double relPFIso03 = ( chIso + TMath::Max( gammaIso + neuHadIso - rho*Aeff, 0.0) )
    / electron->pt;
  
  // barrel/endcap dependent requirments      
  if(fabs(electron->scEta)<1.479) {
    // barrel
    if( relPFIso03 > _EGM2011_RelPFIsoEB[wp] ) return kFALSE;
    
    if(fabs(electron->deltaEtaIn) > _EGM2011_DEtaEB[wp]          )  return kFALSE;
    if(fabs(electron->deltaPhiIn) > _EGM2011_DPhiEB[wp]          )  return kFALSE;
    if(electron->sigiEtaiEta	  > _EGM2011_SigmaIetaIetaEB[wp] )  return kFALSE;
    if(electron->HoverE	          > _EGM2011_HoEEB[wp]           )  return kFALSE;
      
  } else {
    // endcap

    if(electron->pt>20) {
      if( relPFIso03 > _EGM2011_RelPFIsoEEHighPt[wp] ) return kFALSE;
    } else {
      if( relPFIso03 > _EGM2011_RelPFIsoEELowPt[wp] ) return kFALSE;
    }
    
    if(fabs(electron->deltaEtaIn) > _EGM2011_DEtaEE[wp]          )  return kFALSE;
    if(fabs(electron->deltaPhiIn) > _EGM2011_DPhiEE[wp]          )  return kFALSE;
    if(electron->sigiEtaiEta	  > _EGM2011_SigmaIetaIetaEE[wp] )  return kFALSE;
    if(_EGM2011_HoEEE[wp] != -1
       && electron->HoverE	          > _EGM2011_HoEEE[wp]   )  return kFALSE;

  }
  
  return kTRUE;
}


Bool_t passEGM2011(const mithep::TDielectron *dielectron, WorkingPointType wp, double rho)
{
  if(fabs(dielectron->d0_1) > _EGM2011_D0Vtx[wp] )  return kFALSE;
  if(fabs(dielectron->dz_1) > _EGM2011_DZVtx[wp] )  return kFALSE;
  if(fabs(dielectron->d0_2) > _EGM2011_D0Vtx[wp] )  return kFALSE;
  if(fabs(dielectron->dz_2) > _EGM2011_DZVtx[wp] )  return kFALSE;
  
  // conversion rejection
  if(_EGM2011_MissingHits[wp] != -1      // this means that the cut IS to be applied
     && dielectron->nExpHitsInner_1 > _EGM2011_MissingHits[wp]) return kFALSE;
  if(_EGM2011_PassConvVFitCut[wp] != -1  // this means that the cut IS to be applied
     && dielectron->isConv_1)            return kFALSE;
  if(_EGM2011_MissingHits[wp] != -1      // this means that the cut IS to be applied
     && dielectron->nExpHitsInner_2 > _EGM2011_MissingHits[wp]) return kFALSE;
  if(_EGM2011_PassConvVFitCut[wp] != -1  // this means that the cut IS to be applied
     && dielectron->isConv_2)            return kFALSE;

  // Cut on fabs(1/E - 1/p)
//   double theta_1 = 2*atan(exp( - dielectron->scEta_1 ) );
//   // Below, fabs is probably not really needed, but just in case
//   double scE_1 = fabs( dielectron->scEt_1 / sin( theta_1 ) ); 
//   // NOTE: THIS IS NOT TRULY PROPER IMPLEMENTATION OF THE CUT,
//   // the energy should be not SC, but of GsfDielectron
//   double invEMinusInvP_1 = (1/scE_1)*fabs( 1 - dielectron->EoverP_1 );
  double invEMinusInvP_1 = (1/dielectron->ecalE_1)*fabs( 1 - dielectron->EoverP_1 );
  if( _EGM2011_InvEMinusInvP[wp] != -1   // this means that the cut IS to be applied
      && invEMinusInvP_1 > _EGM2011_InvEMinusInvP[wp] ) return kFALSE;
//   double theta_2 = 2*atan(exp( - dielectron->scEta_2 ) );
//   // Below, fabs is probably not really needed, but just in case
//   double scE_2 = fabs( dielectron->scEt_2 / sin( theta_2 ) ); 
//   // NOTE: THIS IS NOT TRULY PROPER IMPLEMENTATION OF THE CUT,
//   // the energy should be not SC, but of GsfDielectron
//   double invEMinusInvP_2 = (1/scE_2)*fabs( 1 - dielectron->EoverP_2 );
  double invEMinusInvP_2 = (1/dielectron->ecalE_2)*fabs( 1 - dielectron->EoverP_2 );
  if( _EGM2011_InvEMinusInvP[wp] != -1   // this means that the cut IS to be applied
      && invEMinusInvP_2 > _EGM2011_InvEMinusInvP[wp] ) return kFALSE;

  // Compute PF isolation:
  // Add energy in the rings with dR 0.0-0.1, 0.1-0.2, 0.2-0.3 
  // to get the total deposits in the cone 0.3.
  double chIso_1 = dielectron->chIso_00_01_1 + dielectron->chIso_01_02_1
    + dielectron->chIso_02_03_1;
  double gammaIso_1 = dielectron->gammaIso_00_01_1 + dielectron->gammaIso_01_02_1
    + dielectron->gammaIso_02_03_1;
  double neuHadIso_1 = dielectron->neuHadIso_00_01_1 + dielectron->neuHadIso_01_02_1
    + dielectron->neuHadIso_02_03_1;
  double Aeff_1 = 0;
  for(int i=0; i<_nEtaBinsForAeff; i++){
    if ( fabs(dielectron->scEta_1) >= _etaLimitsForAeff[i] 
	 && fabs(dielectron->scEta_1) < _etaLimitsForAeff[i+1]){
      Aeff_1 = _AeffDR03[i];
      break;
    }
  }
  double relPFIso03_1 = ( chIso_1 + 
			  TMath::Max( gammaIso_1 + neuHadIso_1 - rho*Aeff_1, 0.0) )
    / dielectron->pt_1;
  
  double chIso_2 = dielectron->chIso_00_01_2 + dielectron->chIso_01_02_2
    + dielectron->chIso_02_03_2;
  double gammaIso_2 = dielectron->gammaIso_00_01_2 + dielectron->gammaIso_01_02_2
    + dielectron->gammaIso_02_03_2;
  double neuHadIso_2 = dielectron->neuHadIso_00_01_2 + dielectron->neuHadIso_01_02_2
    + dielectron->neuHadIso_02_03_2;
  double Aeff_2 = 0;
  for(int i=0; i<_nEtaBinsForAeff; i++){
    if ( fabs(dielectron->scEta_2) >= _etaLimitsForAeff[i] 
	 && fabs(dielectron->scEta_2) < _etaLimitsForAeff[i+1]){
      Aeff_2 = _AeffDR03[i];
      break;
    }
  }
  double relPFIso03_2 = ( chIso_2 + TMath::Max( gammaIso_2 + neuHadIso_2 - rho*Aeff_2, 0.0) )
    / dielectron->pt_2;
  
  // barrel/endcap dependent requirments      
  if(fabs(dielectron->scEta_1)<1.479) {
    // barrel
    if( relPFIso03_1 > _EGM2011_RelPFIsoEB[wp] ) return kFALSE;
    
    if(fabs(dielectron->deltaEtaIn_1) > _EGM2011_DEtaEB[wp]          )  return kFALSE;
    if(fabs(dielectron->deltaPhiIn_1) > _EGM2011_DPhiEB[wp]          )  return kFALSE;
    if(dielectron->sigiEtaiEta_1	  > _EGM2011_SigmaIetaIetaEB[wp] )  return kFALSE;
    if(dielectron->HoverE_1	          > _EGM2011_HoEEB[wp]           )  return kFALSE;
      
  } else {
    // endcap

    if(dielectron->pt_1>20) {
      if( relPFIso03_1 > _EGM2011_RelPFIsoEEHighPt[wp] ) return kFALSE;
    } else {
      if( relPFIso03_1 > _EGM2011_RelPFIsoEELowPt[wp] ) return kFALSE;
    }
    
    if(fabs(dielectron->deltaEtaIn_1) > _EGM2011_DEtaEE[wp]          )  return kFALSE;
    if(fabs(dielectron->deltaPhiIn_1) > _EGM2011_DPhiEE[wp]          )  return kFALSE;
    if(dielectron->sigiEtaiEta_1	  > _EGM2011_SigmaIetaIetaEE[wp] )  return kFALSE;
    if(_EGM2011_HoEEE[wp] != -1 &&
       dielectron->HoverE_1	          > _EGM2011_HoEEE[wp]           )  return kFALSE;

  }
  
  // barrel/endcap dependent requirments for the second electron
  if(fabs(dielectron->scEta_2)<1.479) {
    // barrel
    if( relPFIso03_2 > _EGM2011_RelPFIsoEB[wp] ) return kFALSE;
    
    if(fabs(dielectron->deltaEtaIn_2) > _EGM2011_DEtaEB[wp]          )  return kFALSE;
    if(fabs(dielectron->deltaPhiIn_2) > _EGM2011_DPhiEB[wp]          )  return kFALSE;
    if(dielectron->sigiEtaiEta_2	  > _EGM2011_SigmaIetaIetaEB[wp] )  return kFALSE;
    if(dielectron->HoverE_2	          > _EGM2011_HoEEB[wp]           )  return kFALSE;
      
  } else {
    // endcap

    if(dielectron->pt_2>20) {
      if( relPFIso03_2 > _EGM2011_RelPFIsoEEHighPt[wp] ) return kFALSE;
    } else {
      if( relPFIso03_2 > _EGM2011_RelPFIsoEELowPt[wp] ) return kFALSE;
    }
    
    if(fabs(dielectron->deltaEtaIn_2) > _EGM2011_DEtaEE[wp]          )  return kFALSE;
    if(fabs(dielectron->deltaPhiIn_2) > _EGM2011_DPhiEE[wp]          )  return kFALSE;
    if(dielectron->sigiEtaiEta_2	  > _EGM2011_SigmaIetaIetaEE[wp] )  return kFALSE;
    if(_EGM2011_HoEEE[wp] != -1 &&
       dielectron->HoverE_2	          > _EGM2011_HoEEE[wp]           )  return kFALSE;

  }
  
  return kTRUE;
}


//
// Define electron IDs
//

typedef enum { _EleID_Smurf2011,
	       _EleID_EGM2011_Medium } TEleID_t;

//const TEleID_t _electronID=_EleID_Smurf2011;
// Started working on this but not finished.
const TEleID_t _electronID=_EleID_EGM2011_Medium;

//
// Generic passEleID function
//

template<class EleObj_t>
Bool_t passEleID(const EleObj_t *electron, const mithep::TEventInfo *info=NULL) {
  Bool_t pass=kFALSE;
  switch(_electronID) {
  case _EleID_Smurf2011: pass=passSmurf(electron); break;
  case _EleID_EGM2011_Medium: 
    assert(info);
    pass=passEGM2011(electron,_EleID_EGM2011_Medium,info->rhoLowEta);
  default:
    std::cout << "passEleID: ElectronID is not prepared\n";
  }
  return pass;
}

#endif
