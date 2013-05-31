#ifndef DYTools_HH
#define DYTools_HH


#include <iostream>
#include <math.h>
#include <assert.h>
#include <TEfficiency.h>
#include <TString.h>

#include "EWKAnaDefs.hh"
#include "TElectron.hh"
#include "TDielectron.hh"

//#define _check_Zpeak

namespace DYTools {

  // ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  // Global variable: analysis is either 1D or 2D
  // 1D analysis collapses the rapidity binning
  // The default choice of binnings for 1D and 2D
  //   assign correct values to 
  //             nMassBins, massBinLimits, nYBins
  // ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

  const int energy8TeV = 1; // Set this to 0 for 7 TeV settings
  const int study2D=1;
  const int extendYRangeFor1D=1; // whether |ymax|=9 for 1D study
  const TString analysisTag_USER=""; //(!study2D && extendYRangeFor1D) ? "ymax9" : "";  // extra name to differentiate the analysis files

  // ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  // ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !


  // Global parameters, kinematics, etc
  const bool excludeEcalGap     = false; // For 8 TeV we are not excluding the transition region
  // ECAL constants, analysis-invariant
  const Double_t kECAL_MAX_ETA  = 2.5;
  const Double_t kECAL_GAP_LOW  = 1.4442;
  const Double_t kECAL_GAP_HIGH = 1.566;
  const Double_t kECAL_GAP_MIDDLE = 1.479; // This is where the split barrel/endcap happens
  // Max electron eta
  const double electronEtaMax8TeV = 2.4;
  const double electronEtaMax7TeV = 2.5;
  // For 8 TeV, we do unconventional 2.4 to match the muon channel
  const double electronEtaMax = ( energy8TeV == 1) ? electronEtaMax8TeV : electronEtaMax7TeV;
  const double etMinLead  = 20;
  const double etMinTrail = 10;

  const int maxTnPCanvasDivisions=10; // maximum number of canvas divisions in EventScaleFactors macros

  //const TString strLumiAtECMS="4.7 fb^{-1} at #sqrt{s} = 7 TeV";
//   const TString strLumiAtECMS="4.8 fb^{-1} at #sqrt{s} = 7 TeV";
//   const double lumiAtECMS=4839.;
//   const double lumiAtECMS_accuracy=1.;
  const TString strLumiAtECMS="19.8 fb^{-1} at #sqrt{s} = 8 TeV";
  const double lumiAtECMS=19789.;
  const double lumiAtECMS_accuracy=1.;


  // asymmetric et cuts
  inline
  bool goodEtPair(double et1, double et2) {
    return ( (et1>etMinTrail) && (et2>etMinTrail) &&
	     ((et1>etMinLead) || (et2>etMinLead)) ) ? true : false;
  }


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
  // On 7 TeV data, we measured 2D cross section to Y=2.5 
  // On 8 TeV, it is up to Y=2.4 to be consistent with muons
  const double yRangeMax2D = ( energy8TeV == 1) ? 2.4 : 2.5;
  // For 1D, extend max of Y to ~9, for 2D keep ~2.4-2.5:
  const double yRangeMax =  yRangeMax2D + ((!study2D && extendYRangeFor1D) ? 6.5 : 0);
  const int _nYBinsMax2D=25; // the largest division into Y bins
  // For 7 TeV, for m<200 GeV, the range in Y is 0-2.5, and we have 25 bins,
  //   while for 8 TeV, the range in Y is 0-24, and we have 24 bins.
  // For 7 TeV last mass range has dY = 0.25 wide bins, for 8 TeV it is dY = 0.2
  const int nBinsYLowMass  = (energy8TeV == 1) ? 24 : 25;
  const int nBinsYHighMass = (energy8TeV == 1) ? 12 : 10;
  const int _nYBins2D[_nMassBins2D] = 
    { nBinsYLowMass,// underflow, binned like first mass bin 
      nBinsYLowMass, nBinsYLowMass, nBinsYLowMass, nBinsYLowMass, nBinsYLowMass, 
      nBinsYHighMass
    }; // overflow is neglected

  
  // ------------- 2011 and 2010 content is below ---------------

  // Systematics modes for unfolding and acceptance 
  typedef enum {NORMAL, RESOLUTION_STUDY, FSR_STUDY, ESCALE_RESIDUAL, ESCALE_STUDY, ESCALE_STUDY_RND } TSystematicsStudy_t;

  // Tag and probe fitting constants
  typedef enum {COUNTnCOUNT, COUNTnFIT, FITnFIT} TTnPMethod_t;
  typedef enum {RECO=0, ID=1, HLT=2, HLT_leg1, HLT_leg2, HLT_rndTag} TEfficiencyKind_t;
 
  inline 
  bool efficiencyIsHLT(TEfficiencyKind_t eff) {
    bool yes=false;
    switch(eff) {
    case HLT: case HLT_leg1: case HLT_leg2: case HLT_rndTag: yes=true; break;
    default: yes=false;
    }
    return yes;
  }


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


  const int _nMassBinsTest4 = 4;
  const double _massBinLimitsTest4[_nMassBinsTest4+1] = 
    {15,45,60,120,1500}; // 4 bins with Z-peak region singled-out
  const int _nYBinsMaxTest4=1; // the largest division into Y bins
  const int _nYBinsTest4[_nMassBinsTest4] = { 1, 1, 1, 1 };
  const int _nYBinsTest4_2D[_nMassBinsTest4] = { 2, 2, 2, 1 };

  // Z-peak region 1GeV bins
  const int _nMassBinsZpeak = 60;
  const double _massBinLimitsZpeak[_nMassBinsZpeak+1] = 
    { 60, 61, 62, 63, 64, 65, 66, 67, 68, 69,
      70, 71, 72, 73, 74, 75, 76, 77, 78, 79,
      80, 81, 82, 83, 84, 85, 86, 87, 88, 89,
      90, 91, 92, 93, 94, 95, 96, 97, 98, 99,
      100, 101, 102, 103, 104, 105, 106, 107, 108, 109,
      110, 111, 112, 113, 114, 115, 116, 117, 118, 119,
      120 };
  const int _nYBinsMaxZpeak=1; // the largest division into Y bins
  const int _nYBinsZpeak[_nMassBinsZpeak] = { 
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1
  };

  // generic bin idx finder
  inline
  int _findMassBin(double mass, int nMassBinsLoc, const double *massBinLimitsLoc){  

    int result =-1;
    for(int ibin=0; ibin < nMassBinsLoc; ibin++){
      if( (mass >= massBinLimitsLoc[ibin]) && (mass < massBinLimitsLoc[ibin+1])) {
	result = ibin;
	break;
      }
    }
    return result;

  };

  // some derived bin idx finders -- for debug
  inline int _findMassBin2011(double mass) { return _findMassBin(mass,_nMassBins2011,_massBinLimits2011); }
  inline int _findMassBin13(double mass) { return _findMassBin(mass,_nMassBins13,_massBinLimits13); }
  inline int _findMassBinLumi(double mass) { return _findMassBin(mass,_nMassBinsLumi,_massBinLimitsLumi); }
  inline int _findMassBinZpeak(double mass) { return _findMassBin(mass,_nMassBinsZpeak,_massBinLimitsZpeak); }

// Declare mass binnings
  typedef enum { _MassBins_Undefined, _MassBins_2010, _MassBins_2011, _MassBins_2011_Lumi, _MassBins_2011_2D, _MassBins_test4, _MassBins_Zpeak } TMassBinning_t;


  // ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  // Declare mass binning for the analysis (1D and 2D case)
  // ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  

#ifndef _check_Zpeak
  const DYTools::TMassBinning_t massBinningSet=(study2D) ? _MassBins_2011_2D : _MassBins_2011;
  const int nMassBins=(study2D) ? _nMassBins2D : _nMassBins2011;
  const double *massBinLimits=(study2D) ? _massBinLimits2D : _massBinLimits2011;
  const int *nYBins=(study2D) ? _nYBins2D : _nYBins2011;
  const int nYBinsMax=(study2D) ? _nYBinsMax2D : _nYBinsMax2011;
#else
  const DYTools::TMassBinning_t massBinningSet=(study2D) ? _MassBins_2011_2D : _MassBins_Zpeak;
  const int nMassBins=(study2D) ? _nMassBins2D : _nMassBinsZpeak;
  const double *massBinLimits=(study2D) ? _massBinLimits2D : _massBinLimitsZpeak;
  const int *nYBins=(study2D) ? _nYBins2D : _nYBinsZpeak;
  const int nYBinsMax=(study2D) ? _nYBinsMax2D : _nYBinsMaxZpeak;
#endif

  /*
  const DYTools::TMassBinning_t massBinningSet= _MassBins_test4;
  const int nMassBins= _nMassBinsTest4;
  const double *massBinLimits=_massBinLimitsTest4;
  const int *nYBins=(study2D) ? _nYBinsTest4_2D : _nYBinsTest4;
  const int nYBinsMax= (study2D) ? 2 : 1;
  */
  


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
  inline int findMassBin(double mass) { return _findMassBin(mass,nMassBins,massBinLimits); }


  template<class Idx_t>
  inline double findMassBinCenter(Idx_t idx, int warn_on_error=1) { 
    if ((idx<Idx_t(0)) || (idx>=Idx_t(nMassBins))) {
      if (warn_on_error) std::cout << "\n\tDYTools::findMassBinCenter(" << idx << ") - index error\n" << std::endl;
      return 0;
    }
    return 0.5*(massBinLimits[idx] + massBinLimits[idx+1]);
  }

  inline 
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
  
  inline int findAbsYBin(int massBin, double y){ return findYBin(massBin,fabs(y)); }
  inline int findYBin(double mass, double y) { return findYBin(findMassBin(mass),y); }
  inline int findAbsYBin(double mass, double y) { return findYBin(findMassBin(mass),fabs(y)); }

  // return rapidity value at the Ybin center
  inline 
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
  inline 
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
  inline 
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
  inline 
  int findIndexFlat(double mass, double y){
    int massBin=findMassBin(mass);
    int yBin=findAbsYBin(massBin,y);
    return findIndexFlat(massBin,yBin);
  }

  // 
  // 
  // Unfolding matrix binning
  // info: was getNumberOf2DBins
  inline 
  int getTotalNumberOfBins(){
    int result = 0;
    for(int i=0; i<nMassBins; i++)
      result += nYBins[i];
    return result;
  }

  // Largest number of Ybins
  inline 
  int findMaxYBins(){
    int nYBMax=nYBins[0];
    for (int i=1; i<nMassBins; i++)
      if (nYBins[i]>nYBMax) nYBMax=nYBins[i];
    return nYBMax;
  }

  // Bin limits in rapidity for a particular mass slice
  inline 
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
  typedef enum {ETBINS_UNDEFINED=-1, ETBINS1=1, ETBINS5, ETBINS6, ETBINS6alt, ETBINS6altB, ETBINS7, ETBINS7alt, ETBINS7altB, ETBINS8, ETBINS9} TEtBinSet_t;
  const int nEtBins1 = 1;
  const double etBinLimits1[nEtBins1 + 1] = 
    {10, 500};
  const int nEtBins5 = 5;
  const double etBinLimits5[nEtBins5 + 1] = 
    {10, 20, 30, 40, 50, 500};
  const int nEtBins6 = 6;
  const double etBinLimits6[nEtBins6 + 1] = 
    {10, 15, 20, 30, 40, 50, 500};
  const int nEtBins6alt = 6;
  const double etBinLimits6alt[nEtBins6alt + 1] = 
    {10, 20, 30, 40, 50, 75, 500};
  const int nEtBins6altB = 6;
  const double etBinLimits6altB[nEtBins6altB + 1] = 
    {10, 20, 30, 40, 50, 60, 500};
  const int nEtBins7 = 7;
  const double etBinLimits7[nEtBins7 + 1] = 
    {10, 15, 20, 30, 40, 50, 100, 500};
  const int nEtBins7alt = 7;
  const double etBinLimits7alt[nEtBins7alt + 1] = 
    {10, 15, 20, 30, 40, 50, 75, 500};
  const int nEtBins7altB = 7;
  const double etBinLimits7altB[nEtBins7altB + 1] = 
    {10, 20, 30, 40, 50, 75, 100, 500};
  const int nEtBins8 = 8;
  const double etBinLimits8[nEtBins8 + 1] = 
    {10, 15, 20, 30, 40, 50, 75, 100, 500};
  const int nEtBins9 = 9;
  const double etBinLimits9[nEtBins9 + 1] = 
    {10, 15, 20, 30, 40, 50, 75, 100, 125, 500};

  const int nEtBinsMax = nEtBins9;


  inline 
  int getNEtBins(int binning){
    int n=0;
    switch(binning) {
    case ETBINS_UNDEFINED: n = 0; break;
    case ETBINS1: n = nEtBins1; break;
    case ETBINS5: n = nEtBins5; break;
    case ETBINS6: n = nEtBins6; break;
    case ETBINS6alt: n = nEtBins6alt; break;
    case ETBINS6altB: n = nEtBins6altB; break;
    case ETBINS7: n = nEtBins7; break;
    case ETBINS7alt: n = nEtBins7alt; break;
    case ETBINS7altB: n = nEtBins7altB; break;
    case ETBINS8: n = nEtBins8; break;
    case ETBINS9: n = nEtBins9; break;
    default:
      printf("ERROR: unknown binning requested\n");
      n=0;
    }
    return n;
  }

  inline 
  double *getEtBinLimits(int binning){
    int n = getNEtBins(binning);
    double *limitsOut = new double[n+1];
    const double *limits = NULL;
    switch(binning) {
    case ETBINS1: limits=etBinLimits1; break;
    case ETBINS5: limits=etBinLimits5; break;
    case ETBINS6: limits=etBinLimits6; break;
    case ETBINS6alt: limits=etBinLimits6alt; break;
    case ETBINS6altB: limits=etBinLimits6altB; break;
    case ETBINS7: limits=etBinLimits7; break;
    case ETBINS7alt: limits=etBinLimits7alt; break;
    case ETBINS7altB: limits=etBinLimits7altB; break;
    case ETBINS8: limits=etBinLimits8; break;
    case ETBINS9: limits=etBinLimits9; break;
    default:
      printf("ERROR: unknown/undefined binning requested\n");
      assert(0);
    }
    for(int i=0; i<=n; i++)
      limitsOut[i] = limits[i];
    
    return limitsOut;
  }

  inline 
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

  typedef enum {ETABINS_UNDEFINED=-1, ETABINS1=1, ETABINS2, ETABINS2Negs, ETABINS3, ETABINS3Negs, ETABINS5, ETABINS5Negs, ETABINS4test, ETABINS4testNegs,  ETABINS4alt, ETABINS4altNegs, ETABINS5alt, ETABINS5altNegs, ETABINS8alt, ETABINS8altNegs} TEtaBinSet_t;
  const int nEtaBins1 = 1;
  const double etaBinLimits1[nEtBins1 + 1] = 
    {0, 2.5000001};
  const int nEtaBins2 = 2;
  const double etaBinLimits2[nEtaBins2 + 1] = 
    {0, 1.479, 2.5000001};
  const int nEtaBins3 = 3;
  const double etaBinLimits3[nEtaBins3 + 1] = 
    {0, 1.479, 2.0, 2.5000001};
  const int nEtaBins3Negs = 6;
  const double etaBinLimits3Negs[nEtaBins3Negs + 1] = 
    {-2.5000001, -2.0, -1.479, 0., 1.479, 2.0, 2.5000001};
  const int nEtaBins2Negs = 4;
  const double etaBinLimits2Negs[nEtaBins2Negs + 1] = 
    {-2.500001, -1.479, 0, 1.479, 2.5000001};
  const int nEtaBins5 = 5;
  const double etaBinLimits5[nEtaBins5 + 1] =
    {0, 0.8, 1.4442, 1.566, 2.0, 2.500001 };
  const int nEtaBins4test = 4;
  const double etaBinLimits4test[nEtaBins4test + 1] =
    {0, 0.8, 1.479, 2.0, 2.500001 };
  const int nEtaBins4testNegs = 8;
  const double etaBinLimits4testNegs[nEtaBins4testNegs + 1] =
    {-2.50001, -2.0, -1.479, -0.8, 0, 0.8, 1.479, 2.0, 2.500001 };
  const int nEtaBins4alt = 4;
  const double etaBinLimits4alt[nEtaBins4alt + 1 ] = 
    {0, 0.8, 1.479, 2.2, 2.500001 };
  const int nEtaBins4altNegs = 8;
  const double etaBinLimits4altNegs[nEtaBins4altNegs + 1 ] =
    {-2.500001, -2.2, -1.479, -0.8, 0, 0.8, 1.479, 2.2, 2.500001 };
  const int nEtaBins5alt = 5;
  const double etaBinLimits5alt[nEtaBins5alt + 1 ] = 
    {0, 0.8, 1.4442, 1.556, 2.2, 2.500001 };
  const int nEtaBins5altNegs = 10;
  const double etaBinLimits5altNegs[nEtaBins5altNegs + 1 ] =
    {-2.500001, -2.2, -1.556, -1.4442, -0.8, 0, 0.8, 1.4442, 1.556, 2.2, 2.500001 };
  const int nEtaBins8alt = 8;
  const double etaBinLimits8alt[nEtaBins8alt + 1 ] =
    {0, 0.4, 0.8, 1.0, 1.479, 1.8, 2.0, 2.2, 2.50001 };
  const int nEtaBins8altNegs = 16;
  const double etaBinLimits8altNegs[nEtaBins8altNegs + 1 ] =
    {-2.50001, -2.2, -2.0, -1.8, -1.479, -1.0, -0.8, -0.4, 0, 0.4, 0.8, 1.0, 1.479, 1.8, 2.0, 2.2, 2.5 };

  const int nEtaBinsMax= nEtaBins5altNegs;


  inline 
  int getNEtaBins(int binning){
    int n=0;
    switch(binning) {
    case ETABINS_UNDEFINED: n = 0; break;
    case ETABINS1: n = nEtaBins1; break;
    case ETABINS2: n = nEtaBins2; break;
    case ETABINS2Negs: n = nEtaBins2Negs; break;
    case ETABINS3: n = nEtaBins3; break;
    case ETABINS3Negs: n = nEtaBins3Negs; break;
    case ETABINS5: n = nEtaBins5; break;
    case ETABINS4test: n = nEtaBins4test; break;
    case ETABINS4testNegs: n = nEtaBins4testNegs; break;
    case ETABINS4alt: n = nEtaBins4alt; break;
    case ETABINS4altNegs: n = nEtaBins4altNegs; break;
    case ETABINS5alt: n = nEtaBins5alt; break;
    case ETABINS5altNegs: n = nEtaBins5altNegs; break;
    case ETABINS8alt: n = nEtaBins8alt; break;
    case ETABINS8altNegs: n = nEtaBins8altNegs; break;
    default:
      printf("ERROR: unknown binning requested\n");
      assert(0);
      n=0;
    }
    return n;
  }

  inline 
  double *getEtaBinLimits(int binning){
    int n = getNEtaBins(binning);
    double *limitsOut = new double[n+1];
    const double *limits = NULL;
    switch(binning) {
    case ETABINS1: limits = etaBinLimits1; break;
    case ETABINS2: limits = etaBinLimits2; break;
    case ETABINS2Negs: limits = etaBinLimits2Negs; break;
    case ETABINS3: limits = etaBinLimits3; break;
    case ETABINS3Negs: limits = etaBinLimits3Negs; break;
    case ETABINS5: limits = etaBinLimits5; break;
    case ETABINS4test: limits = etaBinLimits4test; break;
    case ETABINS4testNegs: limits = etaBinLimits4testNegs; break;
    case ETABINS4alt: limits = etaBinLimits4alt; break;
    case ETABINS4altNegs: limits = etaBinLimits4altNegs; break;
    case ETABINS5alt: limits = etaBinLimits5alt; break;
    case ETABINS5altNegs: limits = etaBinLimits5altNegs; break;
    case ETABINS8alt: limits = etaBinLimits8alt; break;
    case ETABINS8altNegs: limits = etaBinLimits8altNegs; break;
    default:
      printf("ERROR: unknown/undefined binning requested\n");
      assert(0);
      n=0;
    }
    for (int i=0; i<=n; ++i) {
      limitsOut[i] = limits[i];
    }
    return limitsOut;
  }

  inline
  int signedEtaBinning(int binning) {
    int yes=0;
    switch(binning) {
    case ETABINS1: 
    case ETABINS2: 
    case ETABINS3: 
    case ETABINS5:
    case ETABINS4test:
    case ETABINS4alt:
    case ETABINS5alt:
    case ETABINS8alt:
      yes=0;
      break;
    case ETABINS2Negs: 
    case ETABINS3Negs:
    case ETABINS4testNegs:
    case ETABINS4altNegs:
    case ETABINS5altNegs:
    case ETABINS8altNegs:
      yes=1;
      break;
    default:
      printf("ERROR: unknown/undefined binning requested\n");
      assert(0);
    }
    return yes;
  }

  inline 
  int findEtaBin(double eta, int binning){
    int result =-1;
    int n = getNEtaBins(binning);
    const double *limits = getEtaBinLimits(binning);
    if (!signedEtaBinning(binning) && (eta<0)) eta=-eta;
    for(int ibin=0; ibin < n; ibin++){
      if( eta >= limits[ibin] && eta < limits[ibin+1]) {
	result = ibin;
	break;
      }
    }
    delete limits;
    return result;
  }


  // Primary vertices 

  //const int nPVBinCount=7;
  //const double nPVLimits[nPVBinCount+1] = { 0., 5., 10., 15., 20., 25., 30., 100. };
  const int nPVBinCount=11;
  const double nPVLimits[nPVBinCount+1] = { 0.5, 2.5, 4.5, 6.5, 8.5, 10.5, 12.5, 14.5, 16.5, 20.5, 24.5, 40.5 };
  //  1., 3., 5., 7., 9., 11., 13., 15., 17., 21., 25., 40. };

  inline int findPUBin(int nPV) { return _findMassBin(double(nPV),nPVBinCount,nPVLimits); }
  
  // 
  // Cross section types
  //
  typedef enum { _cs_None=0, _cs_preFsr, _cs_preFsrNorm, 
		 _cs_preFsrDet, _cs_preFsrDetNorm,
		 _cs_preFsrDetErr, _cs_preFsrDetNormErr,
		 _cs_preFsrDetSystErr, _cs_preFsrDetNormSystErr,
		 _cs_postFsr, _cs_postFsrNorm, 
		 _cs_postFsrDet, _cs_postFsrDetNorm } TCrossSectionKind_t;

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
  inline 
  mithep::TElectron *extractElectron(const mithep::TDielectron *dielectron, int index){
    
    mithep::TElectron *ele = new mithep::TElectron;
    
    if(index == 1){
      ele-> pt                  = dielectron-> pt_1                 ;
      ele-> eta                 = dielectron-> eta_1                ;
      ele-> phi                 = dielectron-> phi_1                ;        

      ele-> trkIso03            = dielectron-> trkIso03_1           ;            
      ele-> emIso03             = dielectron-> emIso03_1            ;             
      ele-> hadIso03            = dielectron-> hadIso03_1           ;            

      ele-> chIso_00_01         = dielectron-> chIso_00_01_1        ;
      ele-> chIso_01_02         = dielectron-> chIso_01_02_1        ;
      ele-> chIso_02_03         = dielectron-> chIso_02_03_1        ;
      ele-> chIso_03_04         = dielectron-> chIso_03_04_1        ;
      ele-> chIso_04_05         = dielectron-> chIso_04_05_1        ;

      ele-> gammaIso_00_01      = dielectron-> gammaIso_00_01_1        ;
      ele-> gammaIso_01_02      = dielectron-> gammaIso_01_02_1        ;
      ele-> gammaIso_02_03      = dielectron-> gammaIso_02_03_1        ;
      ele-> gammaIso_03_04      = dielectron-> gammaIso_03_04_1        ;
      ele-> gammaIso_04_05      = dielectron-> gammaIso_04_05_1        ;

      ele-> neuHadIso_00_01     = dielectron-> neuHadIso_00_01_1        ;
      ele-> neuHadIso_01_02     = dielectron-> neuHadIso_01_02_1        ;
      ele-> neuHadIso_02_03     = dielectron-> neuHadIso_02_03_1        ;
      ele-> neuHadIso_03_04     = dielectron-> neuHadIso_03_04_1        ;
      ele-> neuHadIso_04_05     = dielectron-> neuHadIso_04_05_1        ;

      ele->pfPt                 = dielectron->pfPt_1 ;
      ele->pfEta                = dielectron->pfEta_1 ;
      ele->pfPhi                = dielectron->pfPhi_1 ;

      ele-> d0                  = dielectron-> d0_1                 ;
      ele-> dz                  = dielectron-> dz_1                 ;       

      ele-> scEt                = dielectron-> scEt_1               ;
      ele-> scEta               = dielectron-> scEta_1              ;
      ele-> scPhi               = dielectron-> scPhi_1              ;  
      ele-> ecalE               = dielectron-> ecalE_1             ;              
      ele-> HoverE              = dielectron-> HoverE_1             ;              
      ele-> EoverP              = dielectron-> EoverP_1             ;              
      ele-> fBrem               = dielectron-> fBrem_1             ;              

      ele-> deltaEtaIn          = dielectron-> deltaEtaIn_1         ;          
      ele-> deltaPhiIn          = dielectron-> deltaPhiIn_1         ;          
      ele-> sigiEtaiEta         = dielectron-> sigiEtaiEta_1        ;         
      ele-> partnerDeltaCot     = dielectron-> partnerDeltaCot_1    ;     
      ele-> partnerDist         = dielectron-> partnerDist_1        ;         
      ele->mva                  = dielectron->mva_1                 ;

      ele-> q                   = dielectron-> q_1                  ;                 
      ele-> nExpHitsInner       = dielectron-> nExpHitsInner_1      ;       
      ele->  scID               = dielectron-> scID_1               ;               
      ele->  trkID              = dielectron-> trkID_1              ;
      ele->  typeBits           = dielectron-> typeBits_1           ;

      ele-> hltMatchBits        = dielectron-> hltMatchBits_1       ;       
      ele-> isConv              = dielectron-> isConv_1             ;             
    }else{
      ele-> pt                  = dielectron-> pt_2                 ;
      ele-> eta                 = dielectron-> eta_2                ;
      ele-> phi                 = dielectron-> phi_2                ;        

      ele-> trkIso03            = dielectron-> trkIso03_2           ;            
      ele-> emIso03             = dielectron-> emIso03_2            ;             
      ele-> hadIso03            = dielectron-> hadIso03_2           ;            

      ele-> chIso_00_01         = dielectron-> chIso_00_01_2        ;
      ele-> chIso_01_02         = dielectron-> chIso_01_02_2        ;
      ele-> chIso_02_03         = dielectron-> chIso_02_03_2        ;
      ele-> chIso_03_04         = dielectron-> chIso_03_04_2        ;
      ele-> chIso_04_05         = dielectron-> chIso_04_05_2        ;

      ele-> gammaIso_00_01      = dielectron-> gammaIso_00_01_2        ;
      ele-> gammaIso_01_02      = dielectron-> gammaIso_01_02_2        ;
      ele-> gammaIso_02_03      = dielectron-> gammaIso_02_03_2        ;
      ele-> gammaIso_03_04      = dielectron-> gammaIso_03_04_2        ;
      ele-> gammaIso_04_05      = dielectron-> gammaIso_04_05_2        ;

      ele-> neuHadIso_00_01     = dielectron-> neuHadIso_00_01_2        ;
      ele-> neuHadIso_01_02     = dielectron-> neuHadIso_01_02_2        ;
      ele-> neuHadIso_02_03     = dielectron-> neuHadIso_02_03_2        ;
      ele-> neuHadIso_03_04     = dielectron-> neuHadIso_03_04_2        ;
      ele-> neuHadIso_04_05     = dielectron-> neuHadIso_04_05_2        ;

      ele->pfPt                 = dielectron->pfPt_2 ;
      ele->pfEta                = dielectron->pfEta_2 ;
      ele->pfPhi                = dielectron->pfPhi_2 ;

      ele-> d0                  = dielectron-> d0_2                 ;
      ele-> dz                  = dielectron-> dz_2                 ;       

      ele-> scEt                = dielectron-> scEt_2               ;
      ele-> scEta               = dielectron-> scEta_2              ;
      ele-> scPhi               = dielectron-> scPhi_2              ;  
      ele-> ecalE               = dielectron-> ecalE_2             ;              
      ele-> HoverE              = dielectron-> HoverE_2             ;              
      ele-> EoverP              = dielectron-> EoverP_2             ;              
      ele-> fBrem               = dielectron-> fBrem_2             ;              

      ele-> deltaEtaIn          = dielectron-> deltaEtaIn_2         ;          
      ele-> deltaPhiIn          = dielectron-> deltaPhiIn_2         ;          
      ele-> sigiEtaiEta         = dielectron-> sigiEtaiEta_2        ;         
      ele-> partnerDeltaCot     = dielectron-> partnerDeltaCot_2    ;     
      ele-> partnerDist         = dielectron-> partnerDist_2        ;         
      ele->mva                  = dielectron->mva_2                 ;

      ele-> q                   = dielectron-> q_2                  ;                 
      ele-> nExpHitsInner       = dielectron-> nExpHitsInner_2      ;       
      ele->  scID               = dielectron-> scID_2               ;               
      ele->  trkID              = dielectron-> trkID_2              ;
      ele->  typeBits           = dielectron-> typeBits_2           ;

      ele-> hltMatchBits        = dielectron-> hltMatchBits_2       ;       
      ele-> isConv              = dielectron-> isConv_2             ;             
    }      
    
    return ele;
  }

  //
  // Barrel or endcap
  //

  inline 
  bool isBarrel(double eta){
    bool result = false;
    if(fabs(eta) <= kECAL_GAP_MIDDLE) 
      result = true;
    return result;
  }

  inline 
  bool isEndcap(double eta){
    bool result = false;
    if(fabs(eta) >= kECAL_GAP_MIDDLE && fabs(eta)<electronEtaMax ) 
      result = true;
    return result;
  }

  inline
  bool isEcalGap(double eta) {
    return ( fabs(eta) > kECAL_GAP_LOW && fabs(eta) < kECAL_GAP_HIGH );
  }

  inline 
  bool goodEta(double eta) {
    bool result = true;
    // General requirement to be in range
    if( fabs(eta) >= electronEtaMax ) 
      result = false;

    // Exclude ECAL gap if required
    if( isEcalGap(eta) && excludeEcalGap )
      result = false;

    return result;
  }


  inline 
  bool goodEtaPair(double eta1, double eta2) {
    return (goodEta(eta1) && goodEta(eta2));
  }

  inline int goodEtEtaPair(double et1, double eta1, double et2, double eta2) {
    return (goodEtPair(et1,et2) && goodEtaPair(eta1,eta2));
  }

  //
  // Several functions to calculate efficiency or ratio
  //
  enum {EFF_POISSON, EFF_BINOMIAL, EFF_CLOPPER_PEARSON};

  inline 
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

// --------------------------------------------------------

namespace DYTools {

  inline int checkTotalLumi(double lumi) {
    if (fabs(lumi-DYTools::lumiAtECMS)>lumiAtECMS_accuracy) {
      std::cout << "checkTotalLumi detected mismatch: lumi=" << lumi 
		<< ", DYTools::lumiAtECMS=" << DYTools::lumiAtECMS
		<< " +/- " << DYTools::lumiAtECMS_accuracy << std::endl;
      return 0;
    }
    return 1;
  }
 
};

// ------------------------------------------------------------------


#endif
