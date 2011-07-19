#ifndef DYTools
#define DYTools

#include <iostream>

#include "EWKAnaDefs.hh"
#include "TElectron.hh"
#include "TDielectron.hh"

namespace DYTools {

  // Tag and probe fitting constants
  enum {COUNTnCOUNT, COUNTnFIT, FITnFIT};
  enum {GSF, ID, HLT};
  //
  // Define mass binning
  //
  const int nMassBins = 13;
  const double massBinLimits[nMassBins+1] = 
    // {20,60,120,600}; // 3 bins
    {15,20,30,40,50,60,76,86,96,106,120,150,200,600}; // 13 bins

  const double etMinLead  = 20;
  const double etMinTrail = 10;

  int findMassBin(double mass){
    
    int result =-1;
    for(int ibin=0; ibin < nMassBins; ibin++){
      if( mass >= massBinLimits[ibin] && mass < massBinLimits[ibin+1]) {
	result = ibin;
	break;
      }
    }
    
    return result;
  };
  
  //
  // Define single electron Pt binning
  //
  const int nPtBins = 5;
  const double ptBinLimits[nPtBins+1] = 
    {10, 20, 30, 40, 50, 500};
  
  int findPtBin(double pt){
    
    int result =-1;
    for(int ibin=0; ibin < nPtBins; ibin++){
      if( pt >= ptBinLimits[ibin] && pt < ptBinLimits[ibin+1]) {
	result = ibin;
	break;
      }
    }
    
    return result;
  };

  //
  // Define Et and Eta binning
  //
  enum {ETBINS1, ETBINS5};
  const int nEtBins1 = 1;
  const double etBinLimits1[nEtBins1 + 1] = 
    {10, 500};
  const int nEtBins5 = 5;
  const double etBinLimits5[nEtBins5 + 1] = 
    {10, 20, 30, 40, 50, 500};

  int getNEtBins(int binning){
    int n=0;
    if( binning == ETBINS1 ){
      n = nEtBins1;
    }else if( binning == ETBINS5 ){
      n = nEtBins5;
    }else{
      printf("ERROR: unknown binning requested\n");
      n=0;
    }
    return n;
  }

  double *getEtBinLimits(int binning){
    int n = getNEtBins(binning);
    double *limitsOut = new double[n+1];
    const double *limits = NULL;
    if( binning == ETBINS1 ){
      limits = etBinLimits1;
    }else if( binning == ETBINS5 ){
      limits = etBinLimits5;
    }else{
      printf("ERROR: unknown binning requested\n");
    }
    for(int i=0; i<=n; i++)
      limitsOut[i] = limits[i];
    
    return limitsOut;
  }

  int findEtBin(double et, int binning){
    
    int result =-1;
    int n = getNEtBins(binning);
    double *limits = getEtBinLimits(binning);
    for(int ibin=0; ibin < n; ibin++){
      if( et >= limits[ibin] && et < limits[ibin+1]) {
	result = ibin;
	break;
      }
    }
    delete limits;
    return result;
  };

  enum {ETABINS1, ETABINS2};
  const int nEtaBins1 = 1;
  const double etaBinLimits1[nEtBins1 + 1] = 
    {0, 2.5000001};
  const int nEtaBins2 = 2;
  const double etaBinLimits2[nEtaBins2 + 1] = 
    {0, 1.479, 2.5000001};

  int getNEtaBins(int binning){
    int n=0;
    if( binning == ETABINS1 ){
      n = nEtaBins1;
    }else if( binning == ETABINS2 ){
      n = nEtaBins2;
    }else{
      printf("ERROR: unknown binning requested\n");
      n=0;
    }
    return n;
  }

  double *getEtaBinLimits(int binning){
    int n = getNEtaBins(binning);
    double *limitsOut = new double[n+1];
    const double *limits = NULL;
    if( binning == ETABINS1 ){
      limits = etaBinLimits1;
    }else if( binning == ETABINS2 ){
      limits = etaBinLimits2;
    }else{
      printf("ERROR: unknown binning requested\n");
    }
    for(int i=0; i<=n; i++)
      limitsOut[i] = limits[i];
    
    return limitsOut;
  }

  int findEtaBin(double eta, int binning){
    
    int result =-1;
    int n = getNEtaBins(binning);
    const double *limits = getEtaBinLimits(binning);
    for(int ibin=0; ibin < n; ibin++){
      if( fabs(eta) >= limits[ibin] && fabs(eta) < limits[ibin+1]) {
	result = ibin;
	break;
      }
    }
    delete limits;
    return result;
  };


  // 
  // Triggers vs run numbers
  //
  enum { UNDEF, REL38X, REL39X};
  enum { DATA, F10MC, W11MC};

  UInt_t triggerBits(int sample, int runNum){
    
    // Just "a" trigger
    UInt_t trigger = kHLT_Ele17_SW_L1R; 

    if( sample == DATA) {
      // Triggers as in WW analysis
      // Actually there are no runs below 136K that we presently use
      if((runNum >= 132440) && (runNum <= 137028)) trigger = kHLT_Photon10_L1R;
      //
      if((runNum >= 136033) && (runNum <= 139980)) trigger = kHLT_Ele10_LW_L1R;
      if((runNum >= 140058) && (runNum <= 141882)) trigger = kHLT_Ele15_SW_L1R;
      if((runNum >= 141956) && (runNum <= 144114)) trigger = kHLT_Ele15_SW_CaloEleId_L1R; 
      if((runNum >= 146428) && (runNum <= 147116)) trigger = kHLT_Ele17_SW_CaloEleId_L1R;
      if((runNum >= 147196) && (runNum <= 148058)) trigger = kHLT_Ele17_SW_TightEleId_L1R;
      if((runNum >= 148819) && (runNum <= 149442)) trigger = kHLT_Ele17_SW_TighterEleIdIsol_L1R;
    } else if( sample == F10MC ) {
      trigger = kHLT_Ele17_SW_CaloEleId_L1R;
    } else if( sample == W11MC ) {
      trigger = kHLT_Ele17_SW_CaloEleId_L1R;
    }else
      std::cout << "DYTools:: Unknown sample" << std::endl;
   
    return trigger;
  }

  //
  // Repackage TDielectron->TElectron
  //
  mithep::TElectron *extractElectron(const mithep::TDielectron *dielectron, int index){
    
    mithep::TElectron *ele = new mithep::TElectron;
    
    if(index == 1){
      ele-> pt                  = dielectron-> pt_1                 ;
      ele-> eta                 = dielectron-> eta_1                ;
      ele-> phi                 = dielectron-> phi_1                ;        
      // Ignore transverse mass of (electron+MET) object
      ele-> caloMt              = 0             ;
      ele-> tcMt                = 0             ;
      ele-> pfMt                = 0             ;  
      
      ele-> trkIso03            = dielectron-> trkIso03_1           ;            
      ele-> emIso03             = dielectron-> emIso03_1            ;             
      ele-> hadIso03            = dielectron-> hadIso03_1           ;            
      ele-> d0                  = dielectron-> d0_1                 ;
      ele-> d0Err               = dielectron-> d0Err_1              ;
      ele-> dz                  = dielectron-> dz_1                 ;       
      ele-> scEt                = dielectron-> scEt_1               ;
      ele-> scEta               = dielectron-> scEta_1              ;
      ele-> scPhi               = dielectron-> scPhi_1              ;  
      ele-> HoverE              = dielectron-> HoverE_1             ;              
      ele-> deltaEtaIn          = dielectron-> deltaEtaIn_1         ;          
      ele-> deltaPhiIn          = dielectron-> deltaPhiIn_1         ;          
      ele-> sigiEtaiEta         = dielectron-> sigiEtaiEta_1        ;         
      ele-> nExpHitsInner       = dielectron-> nExpHitsInner_1      ;       
      ele-> partnerDeltaCot     = dielectron-> partnerDeltaCot_1    ;     
      ele-> partnerDist         = dielectron-> partnerDist_1        ;         
      ele-> partnerRadius       = dielectron-> partnerRadius_1      ;       
      ele-> idLL                = dielectron-> idLL_1               ;                
      ele-> q                   = dielectron-> q_1                  ;                 
      ele-> idCat               = dielectron-> idCat_1              ;              
      ele-> hltMatchBits        = dielectron-> hltMatchBits_1       ;       
      ele-> isConv              = dielectron-> isConv_1             ;             
      ele-> isEcalDriven        = dielectron-> isEcalDriven_1       ;       
      ele->  scID               = -1                                ;               
      ele->  isUsed             = true;       
    }else{
      ele-> pt                  = dielectron-> pt_2                 ;
      ele-> eta                 = dielectron-> eta_2                ;
      ele-> phi                 = dielectron-> phi_2                ;        
      // Ignore transverse mass of (electron+MET) object
      ele-> caloMt              = 0             ;
      ele-> tcMt                = 0             ;
      ele-> pfMt                = 0             ;  
      
      ele-> trkIso03            = dielectron-> trkIso03_2           ;            
      ele-> emIso03             = dielectron-> emIso03_2            ;             
      ele-> hadIso03            = dielectron-> hadIso03_2           ;            
      ele-> d0                  = dielectron-> d0_2                 ;
      ele-> d0Err               = dielectron-> d0Err_2              ;
      ele-> dz                  = dielectron-> dz_2                 ;       
      ele-> scEt                = dielectron-> scEt_2               ;
      ele-> scEta               = dielectron-> scEta_2              ;
      ele-> scPhi               = dielectron-> scPhi_2              ;  
      ele-> HoverE              = dielectron-> HoverE_2             ;              
      ele-> deltaEtaIn          = dielectron-> deltaEtaIn_2         ;          
      ele-> deltaPhiIn          = dielectron-> deltaPhiIn_2         ;          
      ele-> sigiEtaiEta         = dielectron-> sigiEtaiEta_2        ;         
      ele-> nExpHitsInner       = dielectron-> nExpHitsInner_2      ;       
      ele-> partnerDeltaCot     = dielectron-> partnerDeltaCot_2    ;     
      ele-> partnerDist         = dielectron-> partnerDist_2        ;         
      ele-> partnerRadius       = dielectron-> partnerRadius_2      ;       
      ele-> idLL                = dielectron-> idLL_2               ;                
      ele-> q                   = dielectron-> q_2                  ;                 
      ele-> idCat               = dielectron-> idCat_2              ;              
      ele-> hltMatchBits        = dielectron-> hltMatchBits_2       ;       
      ele-> isConv              = dielectron-> isConv_2             ;             
      ele-> isEcalDriven        = dielectron-> isEcalDriven_2       ;       
      ele->  scID               = -1                                ;               
      ele->  isUsed             = true;       
    }      
    
    return ele;
  }

  //
  // Barrel or endcap
  //
  const Double_t ECAL_GAP_LOW  = 1.4442;
  const Double_t ECAL_GAP_HIGH = 1.566;

  bool isBarrel(double eta){
    double result = false;
    if(fabs(eta) <= ECAL_GAP_LOW) 
      result = true;
    return result;
  }

  bool isEndcap(double eta){
    double result = false;
    if(fabs(eta) >= ECAL_GAP_HIGH && fabs(eta)<2.5 ) 
      result = true;
    return result;
  }
  
}

#endif
