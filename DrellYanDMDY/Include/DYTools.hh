#ifndef DYTools
#define DYTools

#include <iostream>
#include <math.h>
#include <assert.h>
#include <TEfficiency.h>

#include "EWKAnaDefs.hh"
#include "TElectron.hh"
#include "TDielectron.hh"


namespace DYTools {

  // Constants that define binning in mass and rapidity
  // Note: bin zero is underflow, overflow is neglected
  const int nMassBins2D = 7;
  const double massBinLimits2D[nMassBins2D+1] = 
    {0, // first bin is underflow
     20, 30, 45, 60, 120, 200, 1500
    }; // overflow is very unlikely, do not account for it
  // Rapidity binning is different for different mass bins
  // Note: this implementation neglects underflow and overflow
  // in rapidity.
  const double yRangeMin =  0.0;
  const double yRangeMax =  2.5;
  const int nYBins[nMassBins2D] = 
    { 25,// underflow, binned like first mass bin 
      25, 25, 25, 25, 25, 10,
    }; // overflow is neglected
  // 
  // 
  // Unfolding matrix binning
  int getNumberOf2DBins(){
    int result = 0;
    for(int i=0; i<nMassBins2D; i++)
      result += nYBins[i];
    return result;
  }
  // Largest number of Ybins
  int getYBinsMax(){
    int nYBMax=nYBins[0];
    for (int i=1; i<nMassBins2D; i++)
      if (nYBins[i]>nYBMax) nYBMax=nYBins[i];
    return nYBMax;
  }
  // Bin limits in rapidity for particular mass slice
  double *getYBinLimits2D(int massBin){
    double *result = 0;
    if( massBin < 0 || massBin >= nMassBins2D ) {
      return result;
    }
    int nYBinsThisSlice = nYBins[massBin];
    result = new double[nYBinsThisSlice+1];
    double delta = (yRangeMax - yRangeMin)/double(nYBinsThisSlice);
    for(int i=0; i<nYBinsThisSlice; i++){
      result[i] = yRangeMin + i * delta;
    }
    result[nYBinsThisSlice] = yRangeMax;
    return result;
  }
//   // Note: this HAS TO BE the total number of all 2D bins, that is
//   // the sum of the contents of the nYBins array above
//   const int nUnfoldingBins = 160;
  
  // Functions that return mass and rapidity bin index
  int findMassBin2D(double m){
    
    int result =-1;
    for(int ibin=0; ibin < nMassBins2D; ibin++){
      if( m >= massBinLimits2D[ibin] && m < massBinLimits2D[ibin+1]) {
	result = ibin;
	break;
      }
    }
    
    return result;
    
  }
  
  int findAbsYBin2D(int massBin, double y){
    
    int result = -1;
    if(massBin < 0 || massBin > nMassBins2D) return result;

    double absY = fabs(y);
    
    int nYBinsThisMassRange = nYBins[massBin];
    // Make the array large enough to accommodate any Y binning
    double yBinLimits[5000];
    double width = (yRangeMax - yRangeMin)/double(nYBinsThisMassRange);
    for(int i=0; i<nYBinsThisMassRange; i++){
      yBinLimits[i] = yRangeMin + i * width;
    }
    // The line below could have been done in the loop above,
    // but it is here done separately to avoid rounding uncertainty.
    yBinLimits[nYBinsThisMassRange] = yRangeMax;
    
    for(int i=0; i < nYBinsThisMassRange; i++){
      if( absY >= yBinLimits[i] && absY < yBinLimits[i+1]) {
	result = i;
	break;
      }
    }
    
    return result;
  }
  
  // This function finds a unique 1D index for (index_m, index_y) pair
  int findIndexFlat(int massBin, int yBin){
    
    int result = -1;
    if( massBin < 0 || massBin > nMassBins2D || yBin < 0 || yBin > nYBins[massBin] )
      return result;
    
    result = 0;
    for(int i=0; i< massBin; i++)
      result += nYBins[i];
    
    result += yBin;
    
    return result;
  }

  // This function finds what is the maximum number of Y bins
  // among the all mass slices.
  int findMaxYBins(){
    int result = 0;
    for(int i=0; i<nMassBins2D; i++)
      if( result < nYBins[i])
	result = nYBins[i];
    return result;
  }
  
  // ------------- 2011 and 2010 content is below ---------------

  // Systematics modes for unfolding and acceptance 
  typedef enum {NORMAL, RESOLUTION_STUDY, FSR_STUDY, ESCALE_RESIDUAL, ESCALE_STUDY }
    TSystematicsStudy_t;

  // Tag and probe fitting constants
  typedef enum {COUNTnCOUNT, COUNTnFIT, FITnFIT} TTnPMethod_t;
  typedef enum {GSF, ID, HLT} TEfficiencyKind_t;
 

  //
  // Define mass binning
  //
  // 2010 mass binning
  const int nMassBins13 = 13;
  const double massBinLimits13[nMassBins13+1] = 
    {15,20,30,40,50,60,76,86,96,106,120,150,200,600}; // 13 bins

  // 2011 mass binning
  const int nMassBins2011 = 40;
  const double massBinLimits2011[nMassBins2011+1] = 
    {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 
     81, 86, 91, 96, 101, 106, 110, 115, 120, 126, 133, 141, 
     150, 160, 171, 185, 200, 220, 243, 273, 320, 380, 440, 
     510, 600, 1000, 1500}; // 40 bins

  const int nMassBinsLumi = 9;
  const double massBinLimitsLumi[nMassBinsLumi+1] = 
    {15,20,30,40,50,60,120,150,200,600}; // 9 bins with Z-peak region singled-out

  const double etMinLead  = 20;
  const double etMinTrail = 10;

  int findMassBin(double mass, int nMassBinsLoc, const double *massBinLimitsLoc){  

    int result =-1;
    for(int ibin=0; ibin < nMassBinsLoc; ibin++){
      if( mass >= massBinLimitsLoc[ibin] && mass < massBinLimitsLoc[ibin+1]) {
	result = ibin;
	break;
      }
    }
    return result;

  };

  int findMassBin2011(double mass) { return findMassBin(mass,nMassBins2011,massBinLimits2011); }
  int findMassBin13(double mass) { return findMassBin(mass,nMassBins13,massBinLimits13); }
  int findMassBinLumi(double mass) { return findMassBin(mass,nMassBinsLumi,massBinLimitsLumi); }

// Declare mass binnings
  typedef enum { _MassBins_Undefined, _MassBins_2010, _MassBins_2011, _MassBins_2011_Lumi, _MassBins_2011_2D } TMassBinning_t;


  // Declare mass binning
  const int study2D=1;
  //DYTools::TMassBinning_t massBinningSet=(study2D) ? _MassBins_2011_2D : _MassBins_2011;
  const int nMassBins=(study2D) ? nMassBins2D : nMassBins2011;
  const double *massBinLimits=(study2D) ? massBinLimits2D : massBinLimits2011;


  // find mass bin idx
  int findMassBin(double mass) { return findMassBin(mass,nMassBins,massBinLimits); }


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
  typedef enum {ETBINS1, ETBINS5} TEtBinSet_t;
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

  typedef enum {ETABINS1, ETABINS2} TEtaBinSet_t;
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


  // Primary vertices 

  //const int nPVBinCount=7;
  //const double nPVLimits[nPVBinCount+1] = { 0., 5., 10., 15., 20., 25., 30., 100. };
  const int nPVBinCount=11;
  const double nPVLimits[nPVBinCount+1] = { 0.5, 2.5, 4.5, 6.5, 8.5, 10.5, 12.5, 14.5, 16.5, 20.5, 24.5, 40.5 };
  //  1., 3., 5., 7., 9., 11., 13., 15., 17., 21., 25., 40. };

  int findPUBin(int nPV) { return findMassBin(double(nPV),nPVBinCount,nPVLimits); }

  // 
  // Triggers vs run numbers
  //
  enum { UNDEF, REL38X, REL39X};
  typedef enum { DATA, MC} TDataKind_t;


  //
  // Repackage TDielectron->TElectron
  //
  mithep::TElectron *extractElectron(const mithep::TDielectron *dielectron, int index){
    
    mithep::TElectron *ele = new mithep::TElectron;
    
    if(index == 1){
      ele-> pt                  = dielectron-> pt_1                 ;
      ele-> eta                 = dielectron-> eta_1                ;
      ele-> phi                 = dielectron-> phi_1                ;        

      ele-> trkIso03            = dielectron-> trkIso03_1           ;            
      ele-> emIso03             = dielectron-> emIso03_1            ;             
      ele-> hadIso03            = dielectron-> hadIso03_1           ;            

      ele->pfIso03              = dielectron-> pfIso03_1 ;          
      ele->pfIso04              = dielectron-> pfIso04_1 ;          
      ele->pfPx                 = dielectron->pfPx_1 ;
      ele->pfPy                 = dielectron->pfPy_1 ;

      ele-> d0                  = dielectron-> d0_1                 ;
      ele-> dz                  = dielectron-> dz_1                 ;       

      ele-> scEt                = dielectron-> scEt_1               ;
      ele-> scEta               = dielectron-> scEta_1              ;
      ele-> scPhi               = dielectron-> scPhi_1              ;  
      ele-> HoverE              = dielectron-> HoverE_1             ;              
      ele-> EoverP              = dielectron-> EoverP_1             ;              
      ele-> fBrem               = dielectron-> fBrem_1             ;              

      ele-> deltaEtaIn          = dielectron-> deltaEtaIn_1         ;          
      ele-> deltaPhiIn          = dielectron-> deltaPhiIn_1         ;          
      ele-> sigiEtaiEta         = dielectron-> sigiEtaiEta_1        ;         
      ele-> nExpHitsInner       = dielectron-> nExpHitsInner_1      ;       
      ele-> partnerDeltaCot     = dielectron-> partnerDeltaCot_1    ;     
      ele-> partnerDist         = dielectron-> partnerDist_1        ;         
//       ele->ellhID               = dielectron->ellhID_1 ;
      ele->mva                  = dielectron->mva_1                 ;

      ele-> q                   = dielectron-> q_1                  ;                 
      ele-> hltMatchBits        = dielectron-> hltMatchBits_1       ;       
      ele-> isConv              = dielectron-> isConv_1             ;             
      ele-> isEcalDriven        = dielectron-> isEcalDriven_1       ;       
      ele->  scID               = -1                                ;               
      ele->  trkID              = -1 ;
    }else{
      ele-> pt                  = dielectron-> pt_2                 ;
      ele-> eta                 = dielectron-> eta_2                ;
      ele-> phi                 = dielectron-> phi_2                ;        
      
      ele-> trkIso03            = dielectron-> trkIso03_2           ;            
      ele-> emIso03             = dielectron-> emIso03_2            ;             
      ele-> hadIso03            = dielectron-> hadIso03_2           ;            

      ele->pfIso03              = dielectron-> pfIso03_2 ;          
      ele->pfIso04              = dielectron-> pfIso04_2 ;          
      ele->pfPx                 = dielectron->pfPx_2 ;
      ele->pfPy                 = dielectron->pfPy_2 ;

      ele-> d0                  = dielectron-> d0_2                 ;
      ele-> dz                  = dielectron-> dz_2                 ;       

      ele-> scEt                = dielectron-> scEt_2               ;
      ele-> scEta               = dielectron-> scEta_2              ;
      ele-> scPhi               = dielectron-> scPhi_2              ;  
      ele-> HoverE              = dielectron-> HoverE_2             ;              
      ele-> EoverP              = dielectron-> EoverP_2             ;              
      ele-> fBrem               = dielectron-> fBrem_2             ;              

      ele-> deltaEtaIn          = dielectron-> deltaEtaIn_2         ;          
      ele-> deltaPhiIn          = dielectron-> deltaPhiIn_2         ;          
      ele-> sigiEtaiEta         = dielectron-> sigiEtaiEta_2        ;         
      ele-> nExpHitsInner       = dielectron-> nExpHitsInner_2      ;       
      ele-> partnerDeltaCot     = dielectron-> partnerDeltaCot_2    ;     
      ele-> partnerDist         = dielectron-> partnerDist_2        ;         
//       ele->ellhID               = dielectron->ellhID_2 ;
      ele->mva                  = dielectron->mva_2                 ;

      ele-> q                   = dielectron-> q_2                  ;                 
      ele-> hltMatchBits        = dielectron-> hltMatchBits_2       ;       
      ele-> isConv              = dielectron-> isConv_2             ;             
      ele-> isEcalDriven        = dielectron-> isEcalDriven_2       ;       
      ele->  scID               = -1                                ;               
      ele->  trkID              = -1 ;
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

  //
  // Several functions to calculate efficiency or ratio
  //
  enum {EFF_POISSON, EFF_BINOMIAL, EFF_CLOPPER_PEARSON};
  void calcEfficiency(double nPass, double nTotal, int method, 
		      double &eff, double &effErrLow, double &effErrHigh)
  {

    eff        = 0;
    effErrLow  = 0;
    effErrHigh = 0;

    if(nTotal == 0) return;
      
    // The efficiency is always just the ratio, it is in the error
    // calculation that the methods differ
    eff = nPass/nTotal;

    if( nPass > nTotal && method != EFF_POISSON){
      printf("calcEfficiency::WARNING requested efficiency calculation for nPass > nTotal,\n");
      printf("     switching to EFF_SIMPLE_DIVIDE method\n");
      method = EFF_POISSON;
    }

    if(method == EFF_POISSON){
      // Numerator and denominator are assumed uncorrelated
      if( nPass != 0 )
	effErrLow =  eff*sqrt( 1/nPass + 1/nTotal );
      else
	effErrLow = 0;
      effErrHigh = effErrLow;
    } else if (method == EFF_BINOMIAL){
      // Numerator is a subset of denominator
      effErrLow = sqrt(eff*(1-eff)/nTotal);
      effErrHigh = effErrLow;
    } else if (method == EFF_CLOPPER_PEARSON) {
      // Use Clopper-Pearson method implemented in ROOT to get +-1 sigma
      // asymmertic errors.
      // a) Clopper-Pearson method requires integer pass/total.
      // b) C++ does not have a built-in rounding function
      int nPassInt = int((nPass-floor(nPass)<0.5)?floor(nPass):ceil(nPass)+1e-3);
      int nTotalInt = int((nTotal-floor(nTotal)<0.5)?floor(nTotal):ceil(nTotal)+1e-3);
      effErrLow  = eff - TEfficiency::ClopperPearson(nTotalInt, nPassInt, 0.68269, kFALSE);
      effErrHigh = TEfficiency::ClopperPearson(nTotalInt, nPassInt, 0.68269, kTRUE) - eff;
    } else{
      printf("CalcEfficiency::ERROR: unknown method requested\n");
    }

    return;
  }

}

// ------------------------------------------------------------------


#endif


