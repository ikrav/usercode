#ifndef DYTools_HH
#define DYTools_HH
#define DYTools

#include <iostream>
#include <math.h>
#include <assert.h>
#include <TEfficiency.h>

#include "EWKAnaDefs.hh"
#include "TElectron.hh"
#include "TDielectron.hh"

namespace DYTools {

  // ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  // Global variable: analysis is either 1D or 2D
  // 1D analysis collapses the rapidity binning
  // The default choice of binnings for 1D and 2D
  //   assign correct values to 
  //             nMassBins, massBinLimits, nYBins
  // ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

  const int study2D=1;
  const int extendYRangeFor1D=1; // whether |ymax|=9 for 1D study
  const TString analysisTag_USER=""; //(!study2D && extendYRangeFor1D) ? "ymax9" : "";  // extra name to differentiate the analysis files

  // ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  // ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !


  // Global parameters
  const Double_t kECAL_GAP_LOW  = 1.4442;
  const Double_t kECAL_GAP_HIGH = 1.566;
  const double etMinLead  = 20;
  const double etMinTrail = 10;

  const TString strLumiAtECMS="4.7 fb^{-1} at #sqrt{s} = 7 TeV";


  // Constants that define binning in mass and rapidity
  // Note: bin zero is underflow, overflow is neglected
  const int _nMassBins2D = 7;
  const double _massBinLimits2D[_nMassBins2D+1] = 
    {0, // first bin is underflow
     20, 30, 45, 60, 120, 200, 1500
    }; // overflow is very unlikely, do not account for it
  // Rapidity binning is different for different mass bins
  // Note: this implementation neglects underflow and overflow
  // in rapidity.
  const double yRangeMin =  0.0;
  const double yRangeMax =  2.5 + ((!study2D && extendYRangeFor1D) ? 6.5 : 0);
  const int _nYBinsMax2D=25; // the largest division into Y bins
  const int _nYBins2D[_nMassBins2D] = 
    { 25,// underflow, binned like first mass bin 
      25, 25, 25, 25, 25, 10,
    }; // overflow is neglected

  
  // ------------- 2011 and 2010 content is below ---------------

  // Systematics modes for unfolding and acceptance 
  typedef enum {NORMAL, RESOLUTION_STUDY, FSR_STUDY, ESCALE_RESIDUAL, ESCALE_STUDY, ESCALE_STUDY_RND } TSystematicsStudy_t;

  // Tag and probe fitting constants
  typedef enum {COUNTnCOUNT, COUNTnFIT, FITnFIT} TTnPMethod_t;
  typedef enum {RECO=0, ID=1, HLT=2} TEfficiencyKind_t;
 

  //
  // Define mass binning
  //
  // 2010 mass binning
  const int _nMassBins13 = 13;
  const double _massBinLimits13[_nMassBins13+1] = 
    {15,20,30,40,50,60,76,86,96,106,120,150,200,600}; // 13 bins
  const int _nYBins13[_nMassBins13] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };

  // 2011 mass binning
  const int _nMassBins2011 = 40;
  const double _massBinLimits2011[_nMassBins2011+1] = 
    {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 
     81, 86, 91, 96, 101, 106, 110, 115, 120, 126, 133, 141, 
     150, 160, 171, 185, 200, 220, 243, 273, 320, 380, 440, 
     510, 600, 1000, 1500}; // 40 bins
  const int _nYBinsMax2011=1; // the largest division into Y bins
  const int _nYBins2011[_nMassBins2011] = { 
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1
  };

  const int _nMassBinsLumi = 9;
  const double _massBinLimitsLumi[_nMassBinsLumi+1] = 
    {15,20,30,40,50,60,120,150,200,600}; // 9 bins with Z-peak region singled-out
  const int _nYBinsMaxLumi=1; // the largest division into Y bins
  const int _nYBinsLumi[_nMassBinsLumi] = { 1, 1, 1, 1, 1, 1, 1, 1, 1 };


  // generic bin idx finder
  int _findMassBin(double mass, int nMassBinsLoc, const double *massBinLimitsLoc){  

    int result =-1;
    for(int ibin=0; ibin < nMassBinsLoc; ibin++){
      if( mass >= massBinLimitsLoc[ibin] && mass < massBinLimitsLoc[ibin+1]) {
	result = ibin;
	break;
      }
    }
    return result;

  };

  // some derived bin idx finders -- for debug
  int _findMassBin2011(double mass) { return _findMassBin(mass,_nMassBins2011,_massBinLimits2011); }
  int _findMassBin13(double mass) { return _findMassBin(mass,_nMassBins13,_massBinLimits13); }
  int _findMassBinLumi(double mass) { return _findMassBin(mass,_nMassBinsLumi,_massBinLimitsLumi); }

// Declare mass binnings
  typedef enum { _MassBins_Undefined, _MassBins_2010, _MassBins_2011, _MassBins_2011_Lumi, _MassBins_2011_2D } TMassBinning_t;


  // ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  // Declare mass binning for the analysis (1D and 2D case)
  // ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

  const DYTools::TMassBinning_t massBinningSet=(study2D) ? _MassBins_2011_2D : _MassBins_2011;
  const int nMassBins=(study2D) ? _nMassBins2D : _nMassBins2011;
  const double *massBinLimits=(study2D) ? _massBinLimits2D : _massBinLimits2011;
  const int *nYBins=(study2D) ? _nYBins2D : _nYBins2011;
  const int nYBinsMax=(study2D) ? _nYBinsMax2D : _nYBinsMax2011;
  


  const TString study2Dstr=(study2D) ? "2D" : "1D";
  const TString analysisTag=study2Dstr + analysisTag_USER;
  const int nUnfoldingBinsMax= nMassBins * nYBinsMax;

  // some analyses use 1D always
  //const int nMassBins1D=nMassBins2011;
  //const double *massBinLimits1D=massBinLimits2011;

  // 
  // Define functions which should be used in the code
  //

  // find mass bin idx
  int findMassBin(double mass) { return _findMassBin(mass,nMassBins,massBinLimits); }


  int findYBin(int massBin, double y){
  
    int result = -1;
    if( massBin < 0 || massBin > nMassBins) return result;
    if ( y < yRangeMin  ||  y > yRangeMax ) return result;

    int nYBinsThisMassRange = nYBins[massBin];

    int newCode=1;
    if (newCode) {
      result = int( ( y - yRangeMin )* 1.000001*nYBinsThisMassRange/( yRangeMax - yRangeMin ) );
      double dy=((yRangeMax-yRangeMin)/double(nYBinsThisMassRange));
      //std::cout << "y=" << y << ", dy=" << dy  << ", binIdx=" << result << ", yRange=" << (result*dy) << ".." << (result*dy+dy) << "\n";
      if (result < 0 || result >= nYBinsThisMassRange) {
	std::cout << "y=" << y << ", dy=" << dy  << ", binIdx=" << result << ", yRange=" << (result*dy) << ".." << (result*dy+dy) << "\n";
	result=-1;
      }
    }
    else {
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
	if( y >= yBinLimits[i] && y < yBinLimits[i+1]) {
	  result = i;
	  break;
	}
      }
    }
    
    return result;
  }
  
  int findAbsYBin(int massBin, double y){ return findYBin(massBin,fabs(y)); }
  int findYBin(double mass, double y) { return findYBin(findMassBin(mass),y); }
  int findAbsYBin(double mass, double y) { return findYBin(findMassBin(mass),fabs(y)); }

  // return rapidity value at the Ybin center
  double findAbsYValue(int massBin, int yBin) {
    if (massBin>=nMassBins) {
      std::cout << "ERROR in findAbsYValue: massBin=" << massBin << ", nMassBins=" << nMassBins << "\n";
      return 0;
    }
    int ybinCount=nYBins[massBin];
    if (yBin>ybinCount) {
      std::cout << "ERROR in findAbsYValue: massBin=" << massBin << ", yBin=" << yBin << ", nYBins[massBin]=" << ybinCount << "\n";
      return 0;
    }
    double yh=(yRangeMax-yRangeMin)/double(ybinCount);
    return (yBin+0.5)*yh;
  }
  
  // return rapidity range of Ybin
  void findAbsYValueRange(int massBin, int yBin, double &absYMin, double &absYMax) {
    absYMin=0.; absYMax=0.;
    if (massBin>=nMassBins) {
      std::cout << "ERROR in findAbsYValueRange: massBin=" << massBin << ", nMassBins=" << nMassBins << "\n";
      return;
    }
    int ybinCount=nYBins[massBin];
    if (yBin>ybinCount) {
      std::cout << "ERROR in findAbsYValueRange: massBin=" << massBin << ", yBin=" << yBin << ", nYBins[massBin]=" << ybinCount << "\n";
      return;
    }
    double yh=(yRangeMax-yRangeMin)/double(ybinCount);
    absYMin=yBin*yh;
    absYMax=yBin*yh+yh;
    return;
  }
  

  // This function finds a unique 1D index for (index_m, index_y) pair
  int findIndexFlat(int massBin, int yBin){
    
    int result = -1;
    if( massBin < 0 || massBin > nMassBins || yBin < 0 || yBin > nYBins[massBin] )
      return result;
    
    result = 0;
    for(int i=0; i< massBin; i++)
      result += nYBins[i];
    
    result += yBin;
    
    return result;
  }

  // This function finds a unique 1D index for (index_m, index_y) pair
  int findIndexFlat(double mass, double y){
    int massBin=findMassBin(mass);
    int yBin=findAbsYBin(massBin,y);
    return findIndexFlat(massBin,yBin);
  }

  // 
  // 
  // Unfolding matrix binning
  // info: was getNumberOf2DBins
  int getTotalNumberOfBins(){
    int result = 0;
    for(int i=0; i<nMassBins; i++)
      result += nYBins[i];
    return result;
  }

  // Largest number of Ybins
  int findMaxYBins(){
    int nYBMax=nYBins[0];
    for (int i=1; i<nMassBins; i++)
      if (nYBins[i]>nYBMax) nYBMax=nYBins[i];
    return nYBMax;
  }

  // Bin limits in rapidity for a particular mass slice
  double *getYBinLimits(int massBin){
    double *result = 0;
    if( massBin < 0 || massBin >= nMassBins ) {
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

  // NEEDS clean up
//   // Note: this HAS TO BE the total number of all 2D bins, that is
//   // the sum of the contents of the nYBins array above
//   const int nUnfoldingBins = 160;
  


  //
  // Define single electron Pt binning
  //
  /*
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
  */

  //
  // Define Et and Eta binning
  //
  typedef enum {ETBINS1=1, ETBINS5, ETBINS6} TEtBinSet_t;
  const int nEtBins1 = 1;
  const double etBinLimits1[nEtBins1 + 1] = 
    {10, 500};
  const int nEtBins5 = 5;
  const double etBinLimits5[nEtBins5 + 1] = 
    {10, 20, 30, 40, 50, 500};
  const int nEtBins6 = 6;
  const double etBinLimits6[nEtBins6 + 1] = 
    {10, 15, 20, 30, 40, 50, 500};

  int getNEtBins(int binning){
    int n=0;
    if( binning == ETBINS1 ){
      n = nEtBins1;
    }else if( binning == ETBINS5 ){
      n = nEtBins5;
    }else if( binning == ETBINS6 ){
      n = nEtBins6;
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
    }else if( binning == ETBINS6 ){
      limits = etBinLimits6;
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

  typedef enum {ETABINS1, ETABINS2, ETABINS5} TEtaBinSet_t;
  const int nEtaBins1 = 1;
  const double etaBinLimits1[nEtBins1 + 1] = 
    {0, 2.5000001};
  const int nEtaBins2 = 2;
  const double etaBinLimits2[nEtaBins2 + 1] = 
    {0, 1.479, 2.5000001};
  const int nEtaBins5 = 5;
  const double etaBinLimits5[nEtaBins5 + 1] =
    {0, 0.8, 1.4442, 1.566, 2.0, 2.500001 };

  int getNEtaBins(int binning){
    int n=0;
    if( binning == ETABINS1 ){
      n = nEtaBins1;
    }else if( binning == ETABINS2 ){
      n = nEtaBins2;
    }else if( binning == ETABINS5 ){
      n = nEtaBins5;
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
    }else if( binning == ETABINS5 ){
      limits = etaBinLimits5;
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

  int findPUBin(int nPV) { return _findMassBin(double(nPV),nPVBinCount,nPVLimits); }

  // 
  // Triggers vs run numbers
  //
  enum { UNDEF, REL38X, REL39X};
  typedef enum { DATA, MC} TDataKind_t;


  // The commented out code below is not valid on 2011 data since
  // bit constants have changed in EWKAnaDefs.hh
//   UInt_t triggerBits(int sample, int runNum){
    
//     // Just "a" trigger
//     UInt_t trigger = kHLT_Ele17_SW_L1R; 

//     if( sample == DATA) {
//       // Triggers as in WW analysis
//       // Actually there are no runs below 136K that we presently use
//       if((runNum >= 132440) && (runNum <= 137028)) trigger = kHLT_Photon10_L1R;
//       //
//       if((runNum >= 136033) && (runNum <= 139980)) trigger = kHLT_Ele10_LW_L1R;
//       if((runNum >= 140058) && (runNum <= 141882)) trigger = kHLT_Ele15_SW_L1R;
//       if((runNum >= 141956) && (runNum <= 144114)) trigger = kHLT_Ele15_SW_CaloEleId_L1R; 
//       if((runNum >= 146428) && (runNum <= 147116)) trigger = kHLT_Ele17_SW_CaloEleId_L1R;
//       if((runNum >= 147196) && (runNum <= 148058)) trigger = kHLT_Ele17_SW_TightEleId_L1R;
//       if((runNum >= 148819) && (runNum <= 149442)) trigger = kHLT_Ele17_SW_TighterEleIdIsol_L1R;
//     } else if( sample == F10MC ) {
//       trigger = kHLT_Ele17_SW_CaloEleId_L1R;
//     } else if( sample == W11MC ) {
//       trigger = kHLT_Ele17_SW_CaloEleId_L1R;
//     }else
//       std::cout << "DYTools:: Unknown sample" << std::endl;
   
//     return trigger;
//   }

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

  bool isBarrel(double eta){
    double result = false;
    if(fabs(eta) <= kECAL_GAP_LOW) 
      result = true;
    return result;
  }

  bool isEndcap(double eta){
    double result = false;
    if(fabs(eta) >= kECAL_GAP_HIGH && fabs(eta)<2.5 ) 
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
