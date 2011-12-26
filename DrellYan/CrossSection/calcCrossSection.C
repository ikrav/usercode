#include "TROOT.h"
#include "TSystem.h"
#include "TVectorD.h"
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TAxis.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>                  // functions to format standard I/O

#include "../Include/DYTools.hh"
#include "../Include/UnfoldingTools.hh"
#include "../Include/MyTools.hh"        // miscellaneous helper functions
#include "../Include/DYTools.hh"
// This global constants will be filled from 
// the configuration file. This is not great C++...
TString tagDirYields = "";
TString tagDirConstants = "";
Double_t lumi = 0;
const int nMassBinTh=518;

const TString fileDataYields        ("yields_bg-subtracted.root");

const TString fileMcReferenceYields ("yields_MC_unfolding_reference.root");

// This file contains unfolding matrix 
const TString fileUnfoldingConstants("unfolding_constants.root");

// Contains relative unfolding systematic errors
const TString fileUnfoldingErrors("unfolding_systematics.root");

// Contains relative escale systematic errors
const TString fileEscaleErrors("escale_systematics.root");

const TString fileEfficiencyConstants("event_efficiency_constants.root");

const TString fileScaleFactorConstants("scale_factors.root");

const TString fileAcceptanceConstants("acceptance_constants.root");

const TString fileAcceptanceSystematics("theoretical_uncertainties.root");

const TString fileAcceptanceFSRSystematics("acceptance_FSR_systematics.root");

const TString fileFsrCorrectionConstants("fsr_constants.root");

const TString fileFsrCorrectionSansAccConstants("fsr_constants_sans_acc.root");

// Forward declarations
void readData(TVectorD &v, TVectorD &vErr1, TVectorD &vErr2);

void  applyUnfolding(TVectorD &vin, TVectorD &vinStatErr, TVectorD &vinSystErr,
             TVectorD &vout, TVectorD &voutStatErr, TVectorD &voutSystErr);
void  applyUnfoldingToMc(TString fullUnfoldingConstFileName, 
			 TString fullMcRefYieldsFileName);
void  efficiencyCorrection(TVectorD &vin, TVectorD &vinStatErr, TVectorD &vinSystErr,
             TVectorD &vout, TVectorD &voutStatErr, TVectorD &voutSystErr);
void  acceptanceCorrection(TVectorD &vin, TVectorD &vinStatErr, TVectorD &vinSystErr,
             TVectorD &vout, TVectorD &voutStatErr, TVectorD &voutSystErr);
void  fsrCorrection(TVectorD &vin, TVectorD &vinStatErr, TVectorD &vinSystErr,
             TVectorD &vout, TVectorD &voutStatErr, TVectorD &voutSystErr);
void  fsrCorrectionSansAcceptance(TVectorD &vin, TVectorD &vinStatErr, TVectorD &vinSystErr,
				  TVectorD &vout, TVectorD &voutStatErr, TVectorD &voutSystErr);

void  crossSections(TVectorD &vin, TVectorD &vinStatErr, TVectorD &vinSystErr,
		    TVectorD &vout, TVectorD &voutStatErr, TVectorD &voutSystErr,
		    TVectorD &voutNorm, TVectorD &voutNormStatErr, TVectorD &voutNormSystErr);

void  crossSectionsDET(TVectorD &vin, TVectorD &vinStatErr, TVectorD &vinSystErr,
		       TVectorD &vout, TVectorD &voutStatErr, TVectorD &voutSystErr,
		       TVectorD &voutNorm, TVectorD &voutNormStatErr, TVectorD &voutNormSystErr);

void  postFsrCrossSections(TVectorD &vin, TVectorD &vinStatErr, TVectorD &vinSystErr,
			   TVectorD &vout, TVectorD &voutStatErr, TVectorD &voutSystErr,
			   TVectorD &voutNorm, TVectorD &voutNormStatErr, TVectorD &voutNormSystErr);

void  postFsrCrossSectionsDET(TVectorD &vin, TVectorD &vinStatErr, TVectorD &vinSystErr,
			      TVectorD &vout, TVectorD &voutStatErr, TVectorD &voutSystErr,
			      TVectorD &voutNorm, TVectorD &voutNormStatErr, TVectorD &voutNormSystErr);

void printTableForNotes(TVectorD obs, TVectorD obsErr, 
			TVectorD unf, TVectorD unfErr,
			TVectorD ecor, TVectorD ecorErr,
			TVectorD acor, TVectorD acorErr,
			TVectorD fcor, TVectorD fcorErr);

void printAllCorrections();
void printRelativeSystErrors();

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
//Four plots of R-shape at the same picture

void RShapePlot
(TVectorD relCrossSection, TVectorD relCrossSectionStatErr, 
 TVectorD relCrossSectionDET, TVectorD relCrossSectionDETStatErr, 
 TVectorD relPostFsrCrossSection, TVectorD relPostFsrCrossSectionStatErr, 
 TVectorD relPostFsrCrossSectionDET, TVectorD relPostFsrCrossSectionDETStatErr);

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////


const double lowZMass = 60.0;
const double highZMass = 120.0;
void getNormBinRange(int &firstNormBin, int &lastNormBin);

// Some global arrays for convenience.
// These will contain errors on R
// (after unfolding, etc).

// Absolute error values at the point just before unfolding
TVectorD systBackgrBeforeUnfolding(DYTools::nMassBins);

// Relative error values. These are meant to be AFTER unfolding.
TVectorD systBackgrRelative(DYTools::nMassBins);
TVectorD systEscaleRelative(DYTools::nMassBins);
TVectorD systUnfoldRelative(DYTools::nMassBins);
TVectorD systAccTheoryRelative(DYTools::nMassBins);
TVectorD systAcceptanceRelative(DYTools::nMassBins);

TVectorD systEfficiency(DYTools::nMassBins);
TVectorD systOthers(DYTools::nMassBins);

// ---------------------------------------------------------------
// Main function
// ---------------------------------------------------------------

void calcCrossSection(const TString conf){

  // Read from configuration file only the location of the root files
  ifstream ifs;
  ifs.open(conf.Data());
  assert(ifs.is_open());
  string line;
  int state = 0;
  while(getline(ifs,line)) {
    if(line[0]=='#') continue;
    if(state==0){
      stringstream ss1(line); ss1 >> lumi;
      state++;
    }else if(state==1){
      tagDirYields = TString(line);
      state++;
    }else if(state==2){
      tagDirConstants = TString(line);
      break;
    }
  }
  ifs.close();

  TVectorD signalYields(nMassBins);
  TVectorD signalYieldsStatErr(nMassBins);
  TVectorD signalYieldsSystErr(nMassBins);
  
  TVectorD unfoldedYields(nMassBins);
  TVectorD unfoldedYieldsStatErr(nMassBins);
  TVectorD unfoldedYieldsSystErr(nMassBins);
  
  TVectorD effCorrectedYields(nMassBins);
  TVectorD effCorrectedYieldsStatErr(nMassBins);
  TVectorD effCorrectedYieldsSystErr(nMassBins);
  
  TVectorD accCorrectedYields(nMassBins);
  TVectorD accCorrectedYieldsStatErr(nMassBins);
  TVectorD accCorrectedYieldsSystErr(nMassBins);
  
  TVectorD preFsrYields(nMassBins);
  TVectorD preFsrYieldsStatErr(nMassBins);
  TVectorD preFsrYieldsSystErr(nMassBins);
  
  TVectorD preFsrSansAccYields(nMassBins);
  TVectorD preFsrSansAccYieldsStatErr(nMassBins);
  TVectorD preFsrSansAccYieldsSystErr(nMassBins);
  
  TVectorD absCrossSection(nMassBins);
  TVectorD absCrossSectionStatErr(nMassBins);
  TVectorD absCrossSectionSystErr(nMassBins);
  
  TVectorD absCrossSectionDET(nMassBins);
  TVectorD absCrossSectionStatErrDET(nMassBins);
  TVectorD absCrossSectionSystErrDET(nMassBins);
  
  TVectorD relCrossSection(nMassBins);
  TVectorD relCrossSectionStatErr(nMassBins);
  TVectorD relCrossSectionSystErr(nMassBins);

  TVectorD relCrossSectionDET(nMassBins);
  TVectorD relCrossSectionStatErrDET(nMassBins);
  TVectorD relCrossSectionSystErrDET(nMassBins);

  TVectorD absPostFsrCrossSection(nMassBins);
  TVectorD absPostFsrCrossSectionStatErr(nMassBins);
  TVectorD absPostFsrCrossSectionSystErr(nMassBins);
  
  TVectorD absPostFsrCrossSectionDET(nMassBins);
  TVectorD absPostFsrCrossSectionStatErrDET(nMassBins);
  TVectorD absPostFsrCrossSectionSystErrDET(nMassBins);
  
  TVectorD relPostFsrCrossSection(nMassBins);
  TVectorD relPostFsrCrossSectionStatErr(nMassBins);
  TVectorD relPostFsrCrossSectionSystErr(nMassBins);

  TVectorD relPostFsrCrossSectionDET(nMassBins);
  TVectorD relPostFsrCrossSectionStatErrDET(nMassBins);
  TVectorD relPostFsrCrossSectionSystErrDET(nMassBins);

  systBackgrBeforeUnfolding = 0;
  systEfficiency = 0;
  systOthers = 0;

  // Read data yields from file (background subtraction is already done)
  readData(signalYields, signalYieldsStatErr, signalYieldsSystErr);

  // Apply unfolding
  applyUnfolding(signalYields, signalYieldsStatErr, signalYieldsSystErr,
		 unfoldedYields, unfoldedYieldsStatErr, unfoldedYieldsSystErr);

  // Apply efficiency correction
  efficiencyCorrection(unfoldedYields, unfoldedYieldsStatErr, unfoldedYieldsSystErr,
		       effCorrectedYields, effCorrectedYieldsStatErr, effCorrectedYieldsSystErr);

  // Acceptance corrections
  acceptanceCorrection(effCorrectedYields, effCorrectedYieldsStatErr, effCorrectedYieldsSystErr,
		       accCorrectedYields, accCorrectedYieldsStatErr, accCorrectedYieldsSystErr);
  
  // FSR corrections
  fsrCorrection(accCorrectedYields, accCorrectedYieldsStatErr, accCorrectedYieldsSystErr,
		preFsrYields, preFsrYieldsStatErr, preFsrYieldsSystErr);

  // Alternative (for DET shapes): all corrections except for acceptance correction 
  fsrCorrectionSansAcceptance(effCorrectedYields, effCorrectedYieldsStatErr, effCorrectedYieldsSystErr,
			      preFsrSansAccYields, preFsrSansAccYieldsStatErr, preFsrSansAccYieldsSystErr);
 
  // Calculate absolute and relative cross-sections
  crossSections(preFsrYields, preFsrYieldsStatErr, preFsrYieldsSystErr,
		absCrossSection, absCrossSectionStatErr, absCrossSectionSystErr,
		relCrossSection, relCrossSectionStatErr, relCrossSectionSystErr);

  // Calculate absolute and relative cross-sections DET (shapes, no Acc, but with FSR)
  crossSectionsDET(preFsrSansAccYields, preFsrSansAccYieldsStatErr, preFsrSansAccYieldsSystErr,
		   absCrossSectionDET, absCrossSectionStatErrDET, absCrossSectionSystErrDET,
		   relCrossSectionDET, relCrossSectionStatErrDET, relCrossSectionSystErrDET);

  // Also, calculate absolute and relative cross-sections for post-FSR stage
  postFsrCrossSections(accCorrectedYields, accCorrectedYieldsStatErr, accCorrectedYieldsSystErr,
		absPostFsrCrossSection, absPostFsrCrossSectionStatErr, absPostFsrCrossSectionSystErr,
		relPostFsrCrossSection, relPostFsrCrossSectionStatErr, relPostFsrCrossSectionSystErr);

  // calculate absolute and relative cross-sections DET (shapes, no Acc, no FSR corrections)
  postFsrCrossSectionsDET(effCorrectedYields, effCorrectedYieldsStatErr, effCorrectedYieldsSystErr,
			  absPostFsrCrossSectionDET, absPostFsrCrossSectionStatErrDET, absPostFsrCrossSectionSystErrDET,
			  relPostFsrCrossSectionDET, relPostFsrCrossSectionStatErrDET, relPostFsrCrossSectionSystErrDET);

  // Output
  printAllCorrections();
  printRelativeSystErrors();

  printTableForNotes(signalYields      , signalYieldsStatErr, 
		     unfoldedYields    , unfoldedYieldsStatErr,
		     effCorrectedYields, effCorrectedYieldsStatErr,
		     accCorrectedYields, accCorrectedYieldsStatErr,
		     preFsrYields      , preFsrYieldsStatErr);

  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////
  //Four plots of R-shape at the same picture

  RShapePlot
    (relCrossSection, relCrossSectionStatErr, 
     relCrossSectionDET, relCrossSectionStatErrDET, 
     relPostFsrCrossSection, relPostFsrCrossSectionStatErr, 
     relPostFsrCrossSectionDET, relPostFsrCrossSectionStatErrDET);
  
  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////



}

//-----------------------------------------------------------------
// Read data
//-----------------------------------------------------------------
void readData(TVectorD &v, TVectorD &vErr1, TVectorD &vErr2){

  printf("Load data yields\n"); fflush(stdout);
  TFile fileYields   (TString("../root_files/yields/")+tagDirYields+TString("/")+fileDataYields);
  TVectorD YieldsSignal       = *(TVectorD *)fileYields.FindObjectAny("YieldsSignal");
  TVectorD YieldsSignalErr    = *(TVectorD *)fileYields.FindObjectAny("YieldsSignalErr");
  TVectorD YieldsSignalSystErr= *(TVectorD *)fileYields.FindObjectAny("YieldsSignalSystErr");
  TVectorD BinLimitsForYields = *(TVectorD *)fileYields.FindObjectAny("BinLimitsForYields");

  // Check that the binning is consistent
  bool checkResult = true;
  if( v.GetNoElements() != nMassBins ) checkResult = false;
  if( YieldsSignal.GetNoElements() != nMassBins ) checkResult = false;
  for(int i=0; i<nMassBins+1; i++){
    if( massBinLimits[i] != BinLimitsForYields[i] )
      checkResult = false;
  }
  if( !checkResult ){
    printf("ERROR: inconsistent binning in the inputs\n");
    assert(0);
  }else
    printf("readData: Binning in the inputs is consistent\n");

  // Prepare output yields and errors
  for(int i=0; i<nMassBins; i++){
    v[i] = YieldsSignal[i];
    vErr1[i] = YieldsSignalErr[i];
    // Background systematics should be already in, add 
    // energy scale systematics
    vErr2[i] = YieldsSignalSystErr[i];
    systBackgrBeforeUnfolding[i] = YieldsSignalSystErr[i];
  } 

  fileYields.Close();
  return;
}

//-----------------------------------------------------------------
// Unfold
//-----------------------------------------------------------------
void  applyUnfolding(TVectorD &vin, TVectorD &vinStatErr, TVectorD &vinSystErr,
             TVectorD &vout, TVectorD &voutStatErr, TVectorD &voutSystErr)
{

  // Construct file names
  TString fullUnfoldingConstFileName = TString("../root_files/constants/")+tagDirConstants+TString("/")+fileUnfoldingConstants;
  TString fullUnfoldingErrFileName   = TString("../root_files/systematics/")+tagDirConstants+TString("/")+fileUnfoldingErrors;
  TString fullMcRefYieldsFileName    = TString("../root_files/constants/")+tagDirConstants+TString("/")+fileMcReferenceYields;
  
  // First, propagate through unfolding the signal yields with stat and syst errors
  unfolding::unfold(vin, vout, fullUnfoldingConstFileName);
  unfolding::propagateErrorThroughUnfolding(vinStatErr,voutStatErr, fullUnfoldingConstFileName);
  unfolding::propagateErrorThroughUnfolding(vinSystErr,voutSystErr, fullUnfoldingConstFileName);

  // Second, propagate separately systematic error components that need it.
  // These are already included in the total systematic error above in vinSystErr,
  // however we do it separately so that we can quote the breakdown in the
  // table of systematic errors
  TVectorD systBackgr(vin.GetNoElements());
  unfolding::propagateErrorThroughUnfolding(systBackgrBeforeUnfolding, systBackgr, fullUnfoldingConstFileName);

  // The electron energy scale systematics that is loaded here
  // is estimated on the unfolded yields. So we read it in at this time
  // and will add to the outgoing total systematic below in this function.
  TVectorD systEscale(vin.GetNoElements());
  TString fullEscaleSystErrorsFileName = TString("../root_files/systematics/")
    +tagDirYields+TString("/")
    +fileEscaleErrors;
  TFile fileEscaleSystematics(fullEscaleSystErrorsFileName);
  if( ! fileEscaleSystematics.IsOpen()){
    printf("ERROR: required file with escale errors %s is not found!\n",
	   fullEscaleSystErrorsFileName.Data());
    assert(0);
  }
  TVectorD escaleSystematicsPercent 
    = *(TVectorD *)fileEscaleSystematics.FindObjectAny("escaleSystPercent");
  if( escaleSystematicsPercent.GetNoElements() != DYTools::nMassBins){
    printf("ERROR: Wrong binning of the escale systematics array!\n");
    assert(0);
  }
  for(int i=0; i<nMassBins; i++){
    systEscale[i] = (escaleSystematicsPercent[i]/100.0) * vout[i];
  }

  // Pool together the unfolding systematics and add it to the total systematics
  TVectorD systUnfolding(vin.GetNoElements());
  unfolding::calculateTotalUnfoldingSystError(vin, systUnfolding, 
					      fullUnfoldingConstFileName, 
					      fullUnfoldingErrFileName);
  // Add unfolding and escale systematics to the total systematic error
  for(int i=0; i<DYTools::nMassBins; i++){
    voutSystErr[i] = sqrt( voutSystErr[i]*voutSystErr[i] 
			   + systUnfolding[i]*systUnfolding[i]
			   + systEscale[i]*systEscale[i]);
  }

  // After propagating through unfolding all errors that we had on yields before 
  // unfolding we can compute the relative errors of each kind. While unfolding
  // changes relative errors, all subsequent manipulations do not, so we can 
  // save the errors here.
  for(int i=0; i<DYTools::nMassBins; i++){
    systBackgrRelative[i] = systBackgr[i]/vout[i];
    systEscaleRelative[i] = systEscale[i]/vout[i];
    systUnfoldRelative[i] = systUnfolding[i]/vout[i];
  }  

  // Print the result
  printf("\nUNFOLD: Results for the data, yields:\n");
  printf("                   yields observed        after unfolding            \n");
  for(int i=0; i<DYTools::nMassBins; i++){
    printf("%4.0f-%4.0f   %8.1f +- %6.1f +- %5.1f       %8.1f +- %7.1f +- %6.1f\n",
	   DYTools::massBinLimits[i], DYTools::massBinLimits[i+1],
	   vin[i], vinStatErr[i], vinSystErr[i],
	   vout[i], voutStatErr[i], voutSystErr[i]);
  }

  // Last, do a closure test on MC
  applyUnfoldingToMc(fullUnfoldingConstFileName, fullMcRefYieldsFileName);

  return;
}

void applyUnfoldingToMc(TString fullUnfoldingConstFileName, TString fullMcRefYieldsFileName){

  printf("Load MC reference yields\n"); fflush(stdout);
  TFile fileMcRef(fullMcRefYieldsFileName);
  TVectorD yieldsMcFsrOfRec        = *(TVectorD *)fileMcRef.FindObjectAny("yieldsMcFsrOfRec");
  TVectorD yieldsMcRec             = *(TVectorD *)fileMcRef.FindObjectAny("yieldsMcRec");

  int nBins = yieldsMcRec.GetNoElements();
  TVectorD dNdMmcCheck(nBins);
  dNdMmcCheck = 0;
  unfolding::unfold(yieldsMcRec, dNdMmcCheck, fullUnfoldingConstFileName);

  // Report results
  printf("\nUNFOLD: Check on the MC, yields:\n");
  for(int i=0; i<nBins; i++){
    printf("%4.0f-%4.0f   %10.0f    %10.0f\n",
 	   DYTools::massBinLimits[i],DYTools::massBinLimits[i+1],
 	   yieldsMcFsrOfRec[i],
 	   dNdMmcCheck[i]);
  }

  fileMcRef.Close();
}


void  efficiencyCorrection(TVectorD &vin, TVectorD &vinStatErr, TVectorD &vinSystErr,
             TVectorD &vout, TVectorD &voutStatErr, TVectorD &voutSystErr)
{

  // Read efficiency constants
  printf("Efficiency: Load constants\n"); fflush(stdout);
    
  TFile fileConstants(TString("../root_files/constants/")+tagDirConstants+TString("/")+fileEfficiencyConstants);
  TVectorD efficiencyArray    = *(TVectorD *)fileConstants.FindObjectAny("efficiencyArray");
  TVectorD efficiencyErrArray = *(TVectorD *)fileConstants.FindObjectAny("efficiencyErrArray");

  TFile fileScaleConstants(TString("../root_files/constants/")+tagDirConstants+TString("/")+fileScaleFactorConstants);
  TVectorD rhoDataMc    = *(TVectorD *)fileScaleConstants.FindObjectAny("scaleFactorArray");
  TVectorD rhoDataMcErr = *(TVectorD *)fileScaleConstants.FindObjectAny("scaleFactorErrArray");

  // Check that the binning is consistent
  bool checkResult = true;
  if( efficiencyArray.GetNoElements() != nMassBins) checkResult = false;
  if( !checkResult ){
    printf("Efficiency: ERROR: inconsistent binning in the inputs\n");
    assert(0);
  }else
    printf("Efficiency: Binning in the inputs is consistent\n");

  // Apply the correction
  TVectorD systErrorPropagated(nMassBins);
  TVectorD systErrorAdditional(nMassBins);
  systErrorAdditional = 0;
  for(int i=0; i<nMassBins; i++){
    double effFactor = efficiencyArray[i] * rhoDataMc[i];
    double effErr = effFactor
      * sqrt( efficiencyErrArray[i]*efficiencyErrArray[i]/efficiencyArray[i]/efficiencyArray[i]
	      + rhoDataMcErr[i]*rhoDataMcErr[i]/rhoDataMc[i]/rhoDataMc[i]);
    vout[i]        = vin[i] / effFactor;
    // Statistical error propagated
    voutStatErr[i] = vinStatErr[i] / effFactor;
    // Old systematic error, propagated
    systErrorPropagated[i] = vinSystErr[i] / effFactor;
    // Extra systematic error due to the errors on the efficiency and scale factor
    systErrorAdditional[i] = (vin[i]/effFactor)*(effErr/effFactor);
    systEfficiency[i] = systErrorAdditional[i]/vout[i];
    voutSystErr[i] = sqrt(systErrorPropagated[i]*systErrorPropagated[i]
			  + systErrorAdditional[i]*systErrorAdditional[i]);
  }

  printf("\nEfficiency: Results for the data, yields:\n");
  printf("                after unfolding          eff. factors,%%    rho(data/mc)        efficiency-corrected       syst-err-eff, %%\n");
  for(int i=0; i<nMassBins; i++){
    printf("%4.0f-%4.0f   %8.1f +- %6.1f +- %5.1f   %5.1f +- %5.1f  %5.3f +- %5.3f   %8.1f +- %7.1f +- %6.1f            %6.1f\n",
	   massBinLimits[i],massBinLimits[i+1],
	   vin[i], vinStatErr[i], vinSystErr[i],
	   efficiencyArray[i]*100, efficiencyErrArray[i]*100, rhoDataMc[i], rhoDataMcErr[i],
	   vout[i], voutStatErr[i], voutSystErr[i],
	   systErrorAdditional[i]*100.0/vout[i]);
  }

  fileConstants.Close();
  fileScaleConstants.Close();
  return;

}

void  acceptanceCorrection(TVectorD &vin, TVectorD &vinStatErr, TVectorD &vinSystErr,
			   TVectorD &vout, TVectorD &voutStatErr, TVectorD &voutSystErr){
  // Read efficiency constants
  printf("Acceptance: Load constants\n"); fflush(stdout);
    
  TFile fileConstants(TString("../root_files/constants/")+tagDirConstants+TString("/")+fileAcceptanceConstants);
  TVectorD acceptanceArray    = *(TVectorD *)fileConstants.FindObjectAny("acceptanceArray");
  TVectorD acceptanceErrArray = *(TVectorD *)fileConstants.FindObjectAny("acceptanceErrArray");

  TFile fileSystematics(TString("../root_files/systematics/")+tagDirConstants+TString("/")+fileAcceptanceSystematics);
  TVectorD acceptanceTheoryErrArray = *(TVectorD *)fileSystematics.FindObjectAny("acceptanceTheoryErrArray");

  TFile fileAccFSRSyst(TString("../root_files/systematics/")+tagDirConstants+TString("/")+fileAcceptanceFSRSystematics);
  TVectorD acceptanceFSRErrArray = *(TVectorD *)fileAccFSRSyst.FindObjectAny("accSystPercent");

  // Check that the binning is consistent
  bool checkResult = true;
  if( acceptanceArray.GetNoElements() != nMassBins) checkResult = false;
  //    if( acceptanceTheoryErrArray.GetNoElements() != nMassBins) checkResult = false;
  if( !checkResult ){
    printf("Acceptance: ERROR: inconsistent binning in the inputs\n");
    assert(0);
  }else
    printf("Acceptance: Binning in the inputs is consistent\n");

  // Apply the correction
  TVectorD systErrorPropagated(nMassBins);
  TVectorD systErrorAdditional(nMassBins);
  TVectorD systErrorTheory(nMassBins);
  systErrorAdditional = 0;
  for(int i=0; i<nMassBins; i++){
    double accFactor = acceptanceArray[i];
    double accErr    = acceptanceErrArray[i];
    double accThErr  = acceptanceTheoryErrArray[i];
    double accFSRErr = acceptanceFSRErrArray[i]/100;
    systAccTheoryRelative[i]=acceptanceTheoryErrArray[i];
    vout[i]        = vin[i] / accFactor;
    voutStatErr[i] = vinStatErr[i] / accFactor;
    systErrorPropagated[i] = vinSystErr[i]/accFactor;
    systErrorAdditional[i] = (vin[i]/accFactor)*(accErr/accFactor);
    //voutSystErr[i] = sqrt(systErrorPropagated[i]*systErrorPropagated[i] + systErrorAdditional[i]*systErrorAdditional[i]);
    voutSystErr[i] = sqrt(systErrorPropagated[i]*systErrorPropagated[i] + systErrorAdditional[i]*systErrorAdditional[i] 
            + vout[i]*accThErr*vout[i]*accThErr + vout[i]*accFSRErr*vout[i]*accFSRErr);
    systAcceptanceRelative[i]=sqrt(systErrorAdditional[i]*systErrorAdditional[i]/(vout[i]*vout[i]) 
            + accThErr*accThErr + accFSRErr*accFSRErr);
  }

  printf("\nAcceptance: Results for the data, yields:\n");
  printf("                eff-corrected                acc. factors,%%          acceptance-corrected      syst-err-acc%%,   Theory_Err%% \n");
  for(int i=0; i<nMassBins; i++){
    printf("%4.0f-%4.0f   %8.1f +- %7.1f +- %6.1f   %8.3f +- %7.3f  %9.1f +- %8.1f +- %6.1f  %6.2f   %10.2f\n",
	   massBinLimits[i],massBinLimits[i+1],
	   vin[i], vinStatErr[i], vinSystErr[i],
	   acceptanceArray[i]*100, acceptanceErrArray[i]*100, 
	   vout[i], voutStatErr[i], voutSystErr[i],
	   systErrorAdditional[i]*100/vout[i],100*acceptanceTheoryErrArray[i]);
  }

  fileConstants.Close();
  fileSystematics.Close();
  return;
}

void  fsrCorrection(TVectorD &vin, TVectorD &vinStatErr, TVectorD &vinSystErr,
		    TVectorD &vout, TVectorD &voutStatErr, TVectorD &voutSystErr)
{
  // Read efficiency constants
  printf("FsrCorrection: Load constants\n"); fflush(stdout);
    
  TFile fileConstants(TString("../root_files/constants/")+tagDirConstants+TString("/")+fileFsrCorrectionConstants);
  TVectorD fsrCorrectionArray    = *(TVectorD *)fileConstants.FindObjectAny("fsrCorrectionArray");
  TVectorD fsrCorrectionErrArray = *(TVectorD *)fileConstants.FindObjectAny("fsrCorrectionErrArray");

  // Check that the binning is consistent
  bool checkResult = true;
  if( fsrCorrectionArray.GetNoElements() != nMassBins) checkResult = false;
  if( !checkResult ){
    printf("FsrCorrection: ERROR: inconsistent binning in the inputs\n");
    assert(0);
  }else
    printf("FsrCorrection: Binning in the inputs is consistent\n");

  // Apply the correction
  TVectorD systErrorPropagated(nMassBins);
  for(int i=0; i<nMassBins; i++){
    double fsrFactor = fsrCorrectionArray[i];
    double fsrErr    = fsrCorrectionErrArray[i];
    vout[i]        = vin[i] / fsrFactor;
    voutStatErr[i] = vinStatErr[i] / fsrFactor;
    systErrorPropagated[i] = vinSystErr[i] / fsrFactor;
    double systNew = (vin[i]/fsrFactor)*(fsrErr/fsrFactor);
    voutSystErr[i] = sqrt(systErrorPropagated[i]*systErrorPropagated[i] + systNew*systNew);
  }

  printf("\nFsrCorrection: Results for the data, yields:\n");
  printf("                acc-corrected             fsr. factors            fsr-corrected                  sys-err-fsr\n");
  for(int i=0; i<nMassBins; i++){
    printf("%4.0f-%4.0f   %9.1f +- %8.1f +- %6.1f   %5.3f +- %5.3f  %8.1f +- %7.1f +- %5.1f      %6.2f\n",
	   massBinLimits[i],massBinLimits[i+1],
	   vin[i], vinStatErr[i], vinSystErr[i],
	   fsrCorrectionArray[i], fsrCorrectionErrArray[i], 
	   vout[i], voutStatErr[i], voutSystErr[i],
	   fsrCorrectionErrArray[i]*100.0/fsrCorrectionArray[i]);
  }

  fileConstants.Close();
  return;
}

void  fsrCorrectionSansAcceptance(TVectorD &vin, TVectorD &vinStatErr, TVectorD &vinSystErr,
				  TVectorD &vout, TVectorD &voutStatErr, TVectorD &voutSystErr)
{
  // Read efficiency constants
  printf("FsrCorrection: Load constants\n"); fflush(stdout);
  
  TFile fileConstants(TString("../root_files/constants/")+tagDirConstants+TString("/")+fileFsrCorrectionSansAccConstants);
  TVectorD fsrCorrectionArray    = *(TVectorD *)fileConstants.FindObjectAny("fsrCorrectionArray");
  TVectorD fsrCorrectionErrArray = *(TVectorD *)fileConstants.FindObjectAny("fsrCorrectionErrArray");
  
  // Check that the binning is consistent
  bool checkResult = true;
  if( fsrCorrectionArray.GetNoElements() != nMassBins) checkResult = false;
  if( !checkResult ){
    printf("FsrCorrectionSansAcceptance: ERROR: inconsistent binning in the inputs\n");
    assert(0);
  }else
    printf("FsrCorrectionSansAcceptance: Binning in the inputs is consistent\n");
  
  // Apply the correction
  TVectorD systErrorPropagated(nMassBins);
  for(int i=0; i<nMassBins; i++){
    double fsrFactor = fsrCorrectionArray[i];
    double fsrErr    = fsrCorrectionErrArray[i];
    vout[i]        = vin[i] / fsrFactor;
    voutStatErr[i] = vinStatErr[i] / fsrFactor;
    systErrorPropagated[i] = vinSystErr[i] / fsrFactor;
    double systNew = (vin[i]/fsrFactor)*(fsrErr/fsrFactor);
    voutSystErr[i] = sqrt(systErrorPropagated[i]*systErrorPropagated[i] + systNew*systNew);
  }

  printf("\nThis printout below is FSR correction being applied to data without acceptance correction.\n");
  printf("\nFsrCorrectionSansAcceptance: Results for the data, yields:\n");
  printf("                acc-corrected               fsr. factors             fsr-corrected                  sys-err-fsr\n");
  for(int i=0; i<nMassBins; i++){
    printf("%4.0f-%4.0f   %9.1f +- %8.1f +- %6.1f   %5.3f +- %5.3f  %9.1f +- %8.1f +- %6.1f      %6.2f\n",
	   massBinLimits[i],massBinLimits[i+1],
	   vin[i], vinStatErr[i], vinSystErr[i],
	   fsrCorrectionArray[i], fsrCorrectionErrArray[i], 
	   vout[i], voutStatErr[i], voutSystErr[i],
	   fsrCorrectionErrArray[i]*100.0/fsrCorrectionArray[i]);
  }

  fileConstants.Close();
  return;
}

void  crossSections(TVectorD &vin, TVectorD &vinStatErr, TVectorD &vinSystErr,
		    TVectorD &vout, TVectorD &voutStatErr, TVectorD &voutSystErr,
		    TVectorD &voutNorm, TVectorD &voutNormStatErr, TVectorD &voutNormSystErr)
{

  // Find absolute cross-section
  for(int i=0; i<nMassBins; i++) {
    vout[i] = vin[i] / lumi;
    voutStatErr[i] = vinStatErr[i] / lumi;
    voutSystErr[i] = vinSystErr[i] / lumi;
  }

  // Find normalized cross-section
  int low, high;
  getNormBinRange(low, high);
  double xsecReference = 0;
  double xsecReferenceStatErr = 0;
  double xsecReferenceSystErr = 0;
  for( int i=low; i<=high; i++){
    xsecReference += vout[i];
    xsecReferenceStatErr += voutStatErr[i] * voutStatErr[i];
    xsecReferenceSystErr += voutSystErr[i] * voutSystErr[i];
  }
  xsecReferenceStatErr = sqrt(xsecReferenceStatErr);
  xsecReferenceSystErr = sqrt(xsecReferenceSystErr);

  // Find normalized cross-section
  for(int i=0; i<nMassBins; i++) {
    voutNorm[i] = vout[i] / xsecReference;
    voutNormStatErr[i] = voutNorm[i] 
      * sqrt( xsecReferenceStatErr*xsecReferenceStatErr / xsecReference/xsecReference
	      + voutStatErr[i]*voutStatErr[i] / vout[i]/vout[i]);
    voutNormSystErr[i] = voutNorm[i] 
      * sqrt( xsecReferenceSystErr*xsecReferenceSystErr / xsecReference/xsecReference
	      + voutSystErr[i]*voutSystErr[i] / vout[i]/vout[i]);
  }

  TVectorD normXSec  (DYTools::nMassBins);
  TVectorD normXSecErr   (DYTools::nMassBins);
  TVectorD normXSecErrStat   (DYTools::nMassBins);
  TVectorD normXSecErrSyst   (DYTools::nMassBins);
  TVectorD normXSecByBin  (DYTools::nMassBins);
  TVectorD normXSecErrByBin   (DYTools::nMassBins);
  TVectorD normXSecErrByBinStat   (DYTools::nMassBins);
  TVectorD normXSecErrByBinSyst   (DYTools::nMassBins);

  printf("\nPre FSR cross sections: :\n");
  printf("                    absolute                       normalized +- stat +- sys (total)           (1/sigma)(1/dM)norm +-stat +-syst (total) \n");
  for(int i=0; i<nMassBins; i++){
    double binw = massBinLimits[i+1] - massBinLimits[i];
    normXSec[i]=voutNorm[i];
    normXSecErrStat[i]=voutNormStatErr[i];
    normXSecErrSyst[i]=voutNormSystErr[i];
    normXSecErr[i]=sqrt(voutNormStatErr[i] * voutNormStatErr[i] + voutNormSystErr[i] * voutNormSystErr[i]);

    normXSecByBin[i]=voutNorm[i]/binw;
    normXSecErrByBinStat[i]=voutNormStatErr[i]/binw;
    normXSecErrByBinSyst[i]=voutNormSystErr[i]/binw;
    normXSecErrByBin[i]=sqrt(voutNormStatErr[i] * voutNormStatErr[i] + voutNormSystErr[i] * voutNormSystErr[i])/binw;

    printf("%4.0f-%4.0f   %9.1f +- %8.1f +- %6.1f   %9.6f +- %8.6f +- %8.6f  ( %8.6f )  %10.8f +- %10.8f +- %10.8f  ( %10.8f )    \n",
	   massBinLimits[i],massBinLimits[i+1],
	   vout[i], voutStatErr[i], voutSystErr[i],
	   voutNorm[i], voutNormStatErr[i], voutNormSystErr[i],
	   normXSecErr[i],
	   voutNorm[i]/binw, voutNormStatErr[i]/binw, voutNormSystErr[i]/binw,
	   normXSecErrByBin[i]
	   );
  }
  TString outputDir(TString("../root_files/"));
  gSystem->mkdir(outputDir,kTRUE);
  TString xSecResultFileName(outputDir+TString("/xSec_results.root"));
  TFile fa(xSecResultFileName,"recreate");
  normXSec.Write("normXSec");
  normXSecErr.Write("normXSecErr");
  normXSecErrStat.Write("normXSecErrStat");
  normXSecErrSyst.Write("normXSecErrSyst");
  normXSecByBin.Write("normXSecByBin");
  normXSecErrByBin.Write("normXSecErrByBin");
  normXSecErrByBinStat.Write("normXSecErrByBinStat");
  normXSecErrByBinSyst.Write("normXSecErrByBinSyst");
  fa.Close();
  TVectorD normXSecTh      (nMassBinTh);
  TVectorD normXSecThErr   (nMassBinTh);
 const double xsecTh[nMassBinTh] =
 {234.3088990000, 189.5636730000, 154.2339370000, 128.1240230000, 103.3138130000, 89.6849390000, 74.0037640000, 64.8208330000, 55.5298190000, 47.9486180000,
162.2335730000, 87.9219184000, 52.5499102000, 33.4059051000, 22.4922309000, 16.2577520000, 12.4447444000, 2.1452405100, 2.0837222400, 2.0259685600,
2.0091127000, 1.9741387600, 1.9477418100, 1.9413425600, 1.9129359800, 1.9360725600, 1.9730000900, 2.0119539600, 2.0606883600, 2.1536326500,
2.1963802000, 2.3770811800, 2.5336544900, 2.7480582700, 2.9232081300, 3.2941152500, 3.7071633700, 4.1991465600, 5.0589665800, 5.9611963800,
7.2927069000, 9.5088389300, 12.6242665000, 17.8459438000, 26.8549804000, 47.1862349000, 94.3420578000, 193.4938620000, 226.7245430000, 122.1672890000,
57.7738634000, 31.8359832000, 19.5812869000, 12.7883517000, 9.3285398100, 7.0034570000, 5.5074561800, 4.3571058800, 3.5955949200, 3.0099245800,
2.5135813100, 2.1418748000, 1.8890871700, 1.5939531600, 1.4408045800, 1.2911120300, 1.1634632300, 1.0372046300, 0.9186229000, 0.8326196640,
0.7662111410, 0.7150209660, 0.6697220790, 0.6143397590, 0.5682783650, 0.5029516750, 0.4905461680, 0.4526302470, 0.4309129540, 0.4021197160,
0.3768587860, 0.3559064570, 0.3340920250, 0.3154680840, 0.3016650690, 0.2747160610, 0.2685927580, 0.2574719720, 0.2431129250, 0.2298553340,
0.2180590240, 0.2057435990, 0.1967518120, 0.1872268660, 0.1771664660, 0.1740030350, 0.1663401100, 0.1590678370, 0.1459039850, 0.1445823540,
0.1376292810, 0.1311545160, 0.1311336850, 0.1239912920, 0.1180362660, 0.1149302870, 0.1062668660, 0.1081088070, 0.1028612090, 0.1003053510,
0.0937989850, 0.0909142836, 0.0894981253, 0.0863358608, 0.0820310802, 0.0794862585, 0.0773121340, 0.1463999990, 0.1397183200, 0.1299944750,
0.1225224550, 0.1127808020, 0.1100348610, 0.1031091600, 0.0952096627, 0.0918648228, 0.0865855459, 0.0802503113, 0.0749536282, 0.0734436948,
0.0698027879, 0.0659270202, 0.0619270665, 0.0588428705, 0.0553485067, 0.0538149770, 0.0515021235, 0.0473412892, 0.0460651403, 0.0432861343,
0.0427219944, 0.0401107922, 0.0380066631, 0.0371459167, 0.0349539990, 0.0340959551, 0.0323544500, 0.0303585959, 0.0297261195, 0.0287006924,
0.0275281730, 0.0265415395, 0.0249334594, 0.0241343357, 0.0231597500, 0.0222593915, 0.0215118833, 0.0206525641, 0.0191111094, 0.0188281370,
0.0182268621, 0.0178089411, 0.0165238628, 0.0160227432, 0.0157933277, 0.0150331317, 0.0147353842, 0.0142635123, 0.0132977551, 0.0131845078,
0.0124463529, 0.0122147706, 0.0118779377, 0.0110834379, 0.0110813374, 0.0106187710, 0.0100822538, 0.0092940741, 0.0097757750, 0.0092113146,
0.0089035132, 0.0086283755, 0.0084230576, 0.0079958726, 0.0076765890, 0.0074824478, 0.0074846926, 0.0071223561, 0.0069766440, 0.0066123305,
0.0064785037, 0.0062761122, 0.0060616267, 0.0060007252, 0.0057915107, 0.0053961678, 0.0053031536, 0.0050411714, 0.0049715322, 0.0048132170,
0.0047570639, 0.0045781055, 0.0045775359, 0.0042749321, 0.0042580415, 0.0040166830, 0.0040432225, 0.0037669859, 0.0036947583, 0.0035726540,
0.0035402293, 0.0033201286, 0.0034121564, 0.0032878160, 0.0031719171, 0.0028518319, 0.0028657685, 0.0028601061, 0.0027266151, 0.0026423776,
0.0026009907, 0.0025655296, 0.0024587890, 0.0024477675, 0.0023150569, 0.0023981356, 0.0022020057, 0.0022714171, 0.0021538659, 0.0020516644,
0.0019682733, 0.0019201083, 0.0019002827, 0.0018670390, 0.0018276900, 0.0017533989, 0.0017712593, 0.0016882858, 0.0016889946, 0.0016810962,
0.0016416049, 0.0016075080, 0.0015758766, 0.0015407148, 0.0014922012, 0.0014633704, 0.0014350100, 0.0014068149, 0.0013628119, 0.0013336696,
0.0013050581, 0.0012751613, 0.0012443165, 0.0012210761, 0.0012196392, 0.0011942468, 0.0011582827, 0.0011347788, 0.0011125656, 0.0010873757,
0.0010663820, 0.0010453561, 0.0010346534, 0.0010132240, 0.0009943284, 0.0009802604, 0.0009596552, 0.0009380386, 0.0009191294, 0.0009001258,
0.0008705887, 0.0008539690, 0.0008429864, 0.0008285010, 0.0008117547, 0.0007886230, 0.0007744755, 0.0007586493, 0.0007490780, 0.0007362113,
0.0007208215, 0.0007084377, 0.0006949704, 0.0006781920, 0.0006615548, 0.0006503760, 0.0006374327, 0.0006312936, 0.0006181979, 0.0006075874,
0.0005965148, 0.0005859946, 0.0005740594, 0.0005635871, 0.0005526997, 0.0005423626, 0.0005318028, 0.0005213402, 0.0005107406, 0.0005014435,
0.0004927476, 0.0004830730, 0.0004715156, 0.0004632044, 0.0004547493, 0.0004468069, 0.0004357383, 0.0004281298, 0.0004205780, 0.0004121316,
0.0004067348, 0.0003997200, 0.0003912900, 0.0003838492, 0.0003797172, 0.0003733795, 0.0003672665, 0.0003609797, 0.0003545587, 0.0003477687,
0.0003411518, 0.0003346265, 0.0003291462, 0.0003250723, 0.0003197193, 0.0003141892, 0.0003096618, 0.0003048901, 0.0002995848, 0.0002955056,
0.0002902333, 0.0002855474, 0.0002806956, 0.0002755066, 0.0002712229, 0.0002668692, 0.0002625018, 0.0002581071, 0.0002536747, 0.0005858947,
0.0005615481, 0.0005383547, 0.0005162510, 0.0004951796, 0.0004750881, 0.0004559244, 0.0004376402, 0.0004201896, 0.0004035277, 0.0003876141,
0.0003724117, 0.0003578854, 0.0003440010, 0.0003307254, 0.0003180307, 0.0003058867, 0.0002942657, 0.0002831426, 0.0002724940, 0.0002622967,
0.0002525293, 0.0002431718, 0.0002342050, 0.0002256098, 0.0002173694, 0.0002094673, 0.0002018875, 0.0001946154, 0.0001876374, 0.0001809401,
0.0001745104, 0.0001683367, 0.0001624076, 0.0001567121, 0.0001512399, 0.0001459812, 0.0001409269, 0.0001360678, 0.0001313960, 0.0001269031,
0.0001225818, 0.0001184244, 0.0001144242, 0.0001105742, 0.0001068684, 0.0001033008, 0.0000998655, 0.0000965572, 0.0000933706, 0.0000903009,
0.0000873433, 0.0000844929, 0.0000817456, 0.0000790973, 0.0000765442, 0.0000740820, 0.0000717077, 0.0000694177, 0.0000672088, 0.0000650775,
0.0000630207, 0.0000610358, 0.0000591199, 0.0000572703, 0.0000554846, 0.0000537603, 0.0000520951, 0.0000504868, 0.0000489333, 0.0000474324,
0.0000459820, 0.0000445804, 0.0000432257, 0.0000419163, 0.0000406505, 0.0000394266, 0.0000382432, 0.0000370988, 0.0000374046, 0.0000362920,
0.0000352158, 0.0000341745, 0.0000331669, 0.0000321918, 0.0000312481, 0.0000303346, 0.0000294503, 0.0000285942, 0.0000277654, 0.0000269628,
0.0000261856, 0.0000254327, 0.0000247035, 0.0000239971, 0.0000233127, 0.0000226496, 0.0000220071, 0.0000213844, 0.0000207808, 0.0000201959,
0.0000196288, 0.0000190790, 0.0000185459, 0.0000180291, 0.0000175279, 0.0000170419, 0.0000165704, 0.0000161132, 0.0000156696, 0.0000152393,
0.0000148217, 0.0000144165, 0.0000140234, 0.0000136418, 0.0000132715, 0.0000129121, 0.0000125632, 0.0000122246, 0.0000118958, 0.0000115765,
0.0000112665, 0.0000109654, 0.0000106730, 0.0000103891, 0.0000101133, 0.0000098454, 0.0000095851, 0.0000093322, 0.0000090866, 0.0000088478,
0.0000086159, 0.0000083905, 0.0000081714, 0.0000079585, 0.0000077516, 0.0000075504, 0.0000073549, 0.0000071647, 0.0000069799, 0.0000068001,
0.0000066253, 0.0000064554, 0.0000062901, 0.0000061294, 0.0000059731, 0.0000058210, 0.0000056730, 0.0000055291, 0.0000053891, 0.0000052529,
0.0000051203, 0.0000049914, 0.0000048659, 0.0000047438, 0.0000046250, 0.0000045092, 0.0000043967, 0.0000042871, 0.0000041804, 0.0000040765,
0.0000039754, 0.0000038770, 0.0000037811, 0.0000036879, 0.0000035970, 0.0000035085, 0.0000034223, 0.0000033385, 0.0000032567, 0.0000031771,
0.0000030995, 0.0000030240, 0.0000029504, 0.0000028788, 0.0000028089, 0.0000027409, 0.0000026746, 0.0000026100
};
  const double xsecThErr[nMassBinTh] =
  {0.9736090000, 0.4996990000, 0.2661780000, 0.1677890000, 0.1036220000, 0.0759270000, 0.0569500000, 0.0460280000, 0.0377250000, 0.0306410000,
1.5495162000, 0.8373877370, 0.4946190430, 0.3010118260, 0.2176293490, 0.1435376760, 0.1230714480, 0.0195382123, 0.0182393815, 0.0200309467,
0.0189222377, 0.0195810728, 0.0180861909, 0.0166341069, 0.0190250334, 0.0183965323, 0.0177732666, 0.0182856753, 0.0181248944, 0.0206759002,
0.0213877140, 0.0215881383, 0.0230418106, 0.0269436610, 0.0271169456, 0.0280812730, 0.0367120136, 0.0377327575, 0.0488954821, 0.0578837172,
0.0719417462, 0.0814380390, 0.1143605730, 0.1717725120, 0.2417869490, 0.4271109370, 0.9092901730, 1.8468313700, 2.0597027700, 1.0611175100,
0.5680344220, 0.2691415780, 0.1645799110, 0.1239570990, 0.0866360229, 0.0589891264, 0.0470606406, 0.0417405079, 0.0345002135, 0.0275440181,
0.0246678800, 0.0176256792, 0.0187987952, 0.0148526163, 0.0132142304, 0.0114134489, 0.0112532393, 0.0101071024, 0.0085290581, 0.0078123146,
0.0074906271, 0.0071274465, 0.0055320368, 0.0055132184, 0.0054622208, 0.0040592279, 0.0043408925, 0.0043711776, 0.0036933054, 0.0039284209,
0.0035089432, 0.0035193337, 0.0030272953, 0.0028193867, 0.0023443132, 0.0024718899, 0.0026095623, 0.0020955245, 0.0021521617, 0.0021137040,
0.0021792655, 0.0019000435, 0.0017870031, 0.0016710276, 0.0016570840, 0.0016996471, 0.0015449880, 0.0014959439, 0.0014019533, 0.0012713070,
0.0012514640, 0.0012511919, 0.0011529509, 0.0011688082, 0.0011251825, 0.0008981573, 0.0009830698, 0.0009466416, 0.0008643938, 0.0008853900,
0.0007716729, 0.0007037237, 0.0006819547, 0.0007322556, 0.0008026278, 0.0006745350, 0.0007227847, 0.0014271794, 0.0011857603, 0.0011386009,
0.0010919105, 0.0009128826, 0.0010339643, 0.0009021531, 0.0009446595, 0.0008535231, 0.0008452784, 0.0006427124, 0.0006643324, 0.0007114944,
0.0005739513, 0.0005677666, 0.0005824061, 0.0005067671, 0.0005310312, 0.0004989573, 0.0004133245, 0.0004200988, 0.0004534538, 0.0003904281,
0.0003719468, 0.0003798347, 0.0003098710, 0.0003004792, 0.0003305219, 0.0002629889, 0.0002506170, 0.0002360995, 0.0002548150, 0.0002583708,
0.0002260974, 0.0002281204, 0.0002331338, 0.0002077624, 0.0002022034, 0.0002004780, 0.0002068415, 0.0001828855, 0.0001743010, 0.0001735717,
0.0001818101, 0.0001346604, 0.0001470208, 0.0001513490, 0.0001525853, 0.0001488797, 0.0001430393, 0.0001287955, 0.0001215427, 0.0001316728,
0.0001101619, 0.0001214073, 9.1227661200, 0.0001107811, 0.0001098868, 0.0001035738, 9.0098602700, 8.0348374100, 8.5702800600, 9.2658574400,
7.4200980300, 8.3693179600, 9.7224824900, 7.3724100000, 8.9324717500, 7.8299211100, 8.3049467200, 9.7656172800, 6.8202628600, 8.1631544700,
9.5253691000, 7.5094330000, 7.6849936800, 8.0986801400, 7.2567579600, 8.3052302900, 9.0575720100, 7.7101178600, 9.3635354900, 9.7392976800,
9.9909468200, 8.8696947400, 8.6534233100, 9.4992321300, 9.7747226600, 7.5178510000, 8.4587077400, 9.0872197700, 9.2935240900, 8.5787731000,
9.2418234400, 9.1648751400, 8.7193270900, 6.7877812900, 8.9393610800, 8.8264980200, 9.6207098000, 7.4504300900, 8.6558434200, 9.6637993100,
8.1808641700, 8.9855340200, 9.2062324500, 7.3022094700, 9.1940243600, 9.1081396900, 5.7657804000, 7.9295844200, 7.7076240200, 8.2989131800,
7.1685387600, 9.6173827400, 7.6962590800, 9.8827318100, 7.8648121700, 9.9899018800, 9.2578911200, 7.4050496600, 9.1015376000, 4.2796447127,
4.1758018768, 4.0728657066, 3.9899846826, 3.8989264824, 3.9504566276, 3.8450565524, 3.6688903582, 3.5816089612, 3.8971311556, 3.8069781360,
3.7180257425, 3.6444401182, 3.5455576542, 3.4662067814, 3.1220304455, 3.0089892661, 3.1959791791, 3.1251713611, 3.0556935853, 2.9498333026,
2.8846725232, 2.8215416117, 2.5137007556, 2.4571229710, 2.4033326261, 2.4108357043, 2.3563936002, 2.3024482749, 2.2517393427, 2.2021662470,
2.2349093279, 2.1853253042, 2.0430961712, 1.9898534231, 1.9477176359, 1.9282384005, 1.8871746800, 1.8443689570, 1.7637429001, 1.7240642536,
1.6872163559, 1.6513250056, 1.6181721973, 1.5845146314, 1.5879336413, 1.5530176395, 1.5203943205, 1.4594352829, 1.4239527881, 1.3948933913,
1.3644580054, 1.3358040364, 1.3097509718, 1.2845393340, 1.2544832199, 1.2297690313, 1.2047571846, 1.1786230846, 1.1505963419, 1.1272553942,
1.1039830462, 1.0842471792, 1.0482012374, 1.0275303308, 1.0072703274, 9.8746578241, 9.6268851706, 9.4389994279, 9.2585437170, 9.1201352762,
9.0235341616, 8.8502467749, 8.6023235057, 8.4903632481, 8.2575189264, 8.1000252706, 7.9400217982, 7.7887450408, 7.6103827261, 7.4967478274,
7.3398248868, 7.2622237963, 7.1268176873, 6.9626404327, 6.8348130101, 6.7097418763, 6.5724208399, 6.4482890224, 6.3316725980, 6.2739748753,
6.1475755365, 6.0359546120, 5.9263098679, 5.8248703260, 5.7185408812, 5.6155245549, 5.5256664932, 5.4236304724, 5.3247124065, 0.0000137491,
0.0000131316, 0.0000125453, 0.0000119886, 0.0000114596, 0.0000109570, 0.0000104794, 0.0000100254, 0.0000095935, 0.0000091826, 0.0000087913,
0.0000084187, 0.0000080639, 0.0000077260, 0.0000074040, 0.0000070971, 0.0000068044, 0.0000065251, 0.0000062587, 0.0000060045, 0.0000057618,
0.0000055301, 0.0000053088, 0.0000050974, 0.0000048954, 0.0000047023, 0.0000045177, 0.0000043412, 0.0000041723, 0.0000040108, 0.0000038563,
0.0000037084, 0.0000035667, 0.0000034311, 0.0000033013, 0.0000031768, 0.0000030576, 0.0000029433, 0.0000028338, 0.0000027288, 0.0000026281,
0.0000025315, 0.0000024389, 0.0000023500, 0.0000022647, 0.0000021828, 0.0000021042, 0.0000020287, 0.0000019562, 0.0000018865, 0.0000018196,
0.0000017553, 0.0000016936, 0.0000016342, 0.0000015771, 0.0000015222, 0.0000014694, 0.0000014186, 0.0000013697, 0.0000013227, 0.0000012775,
0.0000012340, 0.0000011921, 0.0000011517, 0.0000011129, 0.0000010755, 0.0000010394, 0.0000010047, 0.0000009713, 0.0000009391, 0.0000009081,
0.0000008781, 0.0000008493, 0.0000008215, 0.0000007947, 0.0000007688, 0.0000007439, 0.0000007198, 0.0000006966, 0.0000007007, 0.0000006783,
0.0000006566, 0.0000006357, 0.0000006155, 0.0000005960, 0.0000005772, 0.0000005591, 0.0000005415, 0.0000005246, 0.0000005082, 0.0000004924,
0.0000004772, 0.0000004624, 0.0000004482, 0.0000004344, 0.0000004211, 0.0000004082, 0.0000003958, 0.0000003837, 0.0000003721, 0.0000003608,
0.0000003500, 0.0000003394, 0.0000003293, 0.0000003194, 0.0000003099, 0.0000003007, 0.0000002917, 0.0000002831, 0.0000002747, 0.0000002667,
0.0000002588, 0.0000002512, 0.0000002439, 0.0000002368, 0.0000002299, 0.0000002232, 0.0000002168, 0.0000002105, 0.0000002045, 0.0000001986,
0.0000001929, 0.0000001874, 0.0000001821, 0.0000001769, 0.0000001719, 0.0000001670, 0.0000001623, 0.0000001577, 0.0000001533, 0.0000001490,
0.0000001448, 0.0000001408, 0.0000001369, 0.0000001331, 0.0000001294, 0.0000001258, 0.0000001224, 0.0000001190, 0.0000001157, 0.0000001126,
0.0000001095, 0.0000001065, 0.0000001036, 0.0000001008, 0.0000000981, 0.0000000954, 0.0000000928, 0.0000000903, 0.0000000879, 0.0000000856,
0.0000000833, 0.0000000810, 0.0000000789, 0.0000000768, 0.0000000748, 0.0000000728, 0.0000000709, 0.0000000690, 0.0000000672, 0.0000000654,
0.0000000637, 0.0000000620, 0.0000000604, 0.0000000588, 0.0000000573, 0.0000000558, 0.0000000544, 0.0000000530, 0.0000000516, 0.0000000503,
0.0000000490, 0.0000000477, 0.0000000465, 0.0000000453, 0.0000000442, 0.0000000431, 0.0000000420, 0.0000000409
  };
  for(int iTh=0; iTh<nMassBinTh; iTh++){
      normXSecTh[iTh]=xsecTh[iTh];
      normXSecThErr[iTh]=xsecThErr[iTh];
  }
  TString outputDir1(TString("../root_files/"));
  gSystem->mkdir(outputDir1,kTRUE);
  TString xSecResultFileName1(outputDir1+TString("/xSecTh_results.root"));
  TFile fb(xSecResultFileName1,"recreate");
  normXSecTh.Write("XSecTh");
  normXSecThErr.Write("XSecThErr");
  fb.Close();

  printf("\nPre FSR cross-section in the Z peak from %3.0f to %3.0f:\n",
	 massBinLimits[low], massBinLimits[high+1]);
  printf("           %9.1f +- %8.1f +- %6.1f \n",
	 xsecReference, xsecReferenceStatErr, xsecReferenceSystErr);
//     printf("check %f %f\n", xsecReferenceStatErr, xsecReferenceSystErr);

  return;
}

void  crossSectionsDET(TVectorD &vin, TVectorD &vinStatErr, TVectorD &vinSystErr,
		       TVectorD &vout, TVectorD &voutStatErr, TVectorD &voutSystErr,
		       TVectorD &voutNorm, TVectorD &voutNormStatErr, TVectorD &voutNormSystErr)
{

  // Find absolute cross-section
  for(int i=0; i<nMassBins; i++) {
    vout[i] = vin[i] / lumi;
    voutStatErr[i] = vinStatErr[i] / lumi;
    voutSystErr[i] = vinSystErr[i] / lumi;
  }

  // Find normalized cross-section
  int low, high;
  getNormBinRange(low, high);
  double xsecReference = 0;
  double xsecReferenceStatErr = 0;
  double xsecReferenceSystErr = 0;
  for( int i=low; i<=high; i++){
    xsecReference += vout[i];
    xsecReferenceStatErr += voutStatErr[i] * voutStatErr[i];
    xsecReferenceSystErr += voutSystErr[i] * voutSystErr[i];
  }
  xsecReferenceStatErr = sqrt(xsecReferenceStatErr);
  xsecReferenceSystErr = sqrt(xsecReferenceSystErr);

  // Find normalized cross-section
  for(int i=0; i<nMassBins; i++) {
    voutNorm[i] = vout[i] / xsecReference;
    voutNormStatErr[i] = voutNorm[i] 
      * sqrt( xsecReferenceStatErr*xsecReferenceStatErr / xsecReference/xsecReference
	      + voutStatErr[i]*voutStatErr[i] / vout[i]/vout[i]);
    voutNormSystErr[i] = voutNorm[i] 
      * sqrt( xsecReferenceSystErr*xsecReferenceSystErr / xsecReference/xsecReference
	      + voutSystErr[i]*voutSystErr[i] / vout[i]/vout[i]);
  }

  printf("\nPre FSR DET shape: :\n");
  printf("                    absolute                      normalized +-stat +-sys (total)\n");
  for(int i=0; i<nMassBins; i++){
    printf("%4.0f-%4.0f   %9.1f +- %8.1f +- %6.1f   %11.6f +- %10.6f +- %10.6f  ( %10.6f )\n",
	   massBinLimits[i],massBinLimits[i+1],
	   vout[i], voutStatErr[i], voutSystErr[i],
	   voutNorm[i], voutNormStatErr[i], voutNormSystErr[i],
	   sqrt(voutNormStatErr[i] * voutNormStatErr[i] + voutNormSystErr[i] * voutNormSystErr[i]) );
  }

  return;
}


void  postFsrCrossSections(TVectorD &vin, TVectorD &vinStatErr, TVectorD &vinSystErr,
		    TVectorD &vout, TVectorD &voutStatErr, TVectorD &voutSystErr,
		    TVectorD &voutNorm, TVectorD &voutNormStatErr, TVectorD &voutNormSystErr)
{

  // Find absolute cross-section
  for(int i=0; i<nMassBins; i++) {
    vout[i] = vin[i] / lumi;
    voutStatErr[i] = vinStatErr[i] / lumi;
    voutSystErr[i] = vinSystErr[i] / lumi;
  }

  // Find normalized cross-section
  int low, high;
  getNormBinRange(low, high);
  double xsecReference = 0;
  double xsecReferenceStatErr = 0;
  double xsecReferenceSystErr = 0;
  for( int i=low; i<=high; i++){
    xsecReference += vout[i];
    xsecReferenceStatErr += voutStatErr[i] * voutStatErr[i];
    xsecReferenceSystErr += voutSystErr[i] * voutSystErr[i];
  }
  xsecReferenceStatErr = sqrt(xsecReferenceStatErr);
  xsecReferenceSystErr = sqrt(xsecReferenceSystErr);

  // Find normalized cross-section
  for(int i=0; i<nMassBins; i++) {
    voutNorm[i] = vout[i] / xsecReference;
    voutNormStatErr[i] = voutNorm[i] 
      * sqrt( xsecReferenceStatErr*xsecReferenceStatErr / xsecReference/xsecReference
	      + voutStatErr[i]*voutStatErr[i] / vout[i]/vout[i]);
    voutNormSystErr[i] = voutNorm[i] 
      * sqrt( xsecReferenceSystErr*xsecReferenceSystErr / xsecReference/xsecReference
	      + voutSystErr[i]*voutSystErr[i] / vout[i]/vout[i]);
  }

  printf("\nPost FSR cross sections: :\n");
  printf("                    absolute                     normalized +-stat +-sys (total)\n");
  for(int i=0; i<nMassBins; i++){
    printf("%4.0f-%4.0f   %9.1f +- %8.1f +- %6.1f   %11.6f +- %10.6f +- %10.6f  ( %10.6f )\n",
	   massBinLimits[i],massBinLimits[i+1],
	   vout[i], voutStatErr[i], voutSystErr[i],
	   voutNorm[i], voutNormStatErr[i], voutNormSystErr[i],
	   sqrt(voutNormStatErr[i] * voutNormStatErr[i] + voutNormSystErr[i] * voutNormSystErr[i]) );
  }
  printf("\nPostFsr cross-section in the Z peak from %3.0f to %3.0f:\n",
	 massBinLimits[low], massBinLimits[high+1]);
  printf("           %9.1f +- %8.1f +- %6.1f \n",
	 xsecReference, xsecReferenceStatErr, xsecReferenceSystErr);
//     printf("check %f %f\n", xsecReferenceStatErr, xsecReferenceSystErr);

  return;
}

void  postFsrCrossSectionsDET(TVectorD &vin, TVectorD &vinStatErr, TVectorD &vinSystErr,
			      TVectorD &vout, TVectorD &voutStatErr, TVectorD &voutSystErr,
			      TVectorD &voutNorm, TVectorD &voutNormStatErr, TVectorD &voutNormSystErr)
{
  
  // Find absolute cross-section
  for(int i=0; i<nMassBins; i++) {
    vout[i] = vin[i] / lumi;
    voutStatErr[i] = vinStatErr[i] / lumi;
    voutSystErr[i] = vinSystErr[i] / lumi;
  }
  
  // Find normalized cross-section
  int low, high;
  getNormBinRange(low, high);
  double xsecReference = 0;
  double xsecReferenceStatErr = 0;
  double xsecReferenceSystErr = 0;
  for( int i=low; i<=high; i++){
    xsecReference += vout[i];
    xsecReferenceStatErr += voutStatErr[i] * voutStatErr[i];
    xsecReferenceSystErr += voutSystErr[i] * voutSystErr[i];
  }
  xsecReferenceStatErr = sqrt(xsecReferenceStatErr);
  xsecReferenceSystErr = sqrt(xsecReferenceSystErr);

  // Find normalized cross-section
  for(int i=0; i<nMassBins; i++) {
    voutNorm[i] = vout[i] / xsecReference;
    voutNormStatErr[i] = voutNorm[i] 
      * sqrt( xsecReferenceStatErr*xsecReferenceStatErr / xsecReference/xsecReference
	      + voutStatErr[i]*voutStatErr[i] / vout[i]/vout[i]);
    voutNormSystErr[i] = voutNorm[i] 
      * sqrt( xsecReferenceSystErr*xsecReferenceSystErr / xsecReference/xsecReference
	      + voutSystErr[i]*voutSystErr[i] / vout[i]/vout[i]);
  }

  printf("\nPost FSR DET shape: :\n");
  printf("                    absolute                     normalized +-stat +-sys (total)\n");
  for(int i=0; i<nMassBins; i++){
    printf("%4.0f-%4.0f   %9.1f +- %8.1f +- %8.1f   %11.6f +- %10.6f +- %10.6f  ( %10.6f )\n",
	   massBinLimits[i],massBinLimits[i+1],
	   vout[i], voutStatErr[i], voutSystErr[i],
	   voutNorm[i], voutNormStatErr[i], voutNormSystErr[i],
	   sqrt(voutNormStatErr[i] * voutNormStatErr[i] + voutNormSystErr[i] * voutNormSystErr[i]) );
  }

  return;
}

void printTableForNotes(TVectorD obs, TVectorD obsErr, 
			TVectorD unf, TVectorD unfErr,
			TVectorD ecor, TVectorD ecorErr,
			TVectorD acor, TVectorD acorErr,
			TVectorD fcor, TVectorD fcorErr)
{

  printf("\n\nLatex table for notes\n");
  printf("               obs-bg                  unfolded                 eff-corrected                acc-corrected              fsr-corrected\n");
  for(int i=0; i<nMassBins; i++){
    printf("%4.0f-%4.0f &", massBinLimits[i],massBinLimits[i+1]);
    printf("  $%8.1f \\pm %6.1f$ &", obs[i] , obsErr[i]);
    printf("  $%8.1f \\pm %6.1f$ &", unf[i] , unfErr[i]);
    printf("  $%8.1f \\pm %6.1f$ &", ecor[i], ecorErr[i]);
    printf("  $%9.1f \\pm %8.1f$ &", acor[i], acorErr[i]);
    printf("  $%9.1f \\pm %8.1f$ \\\\", fcor[i], fcorErr[i]);
    printf("\n");
  }

  return;
}


void printAllCorrections(){

  TFile fileConstantsEff(TString("../root_files/constants/")+tagDirConstants+TString("/")+fileEfficiencyConstants);
  TVectorD efficiencyArray    = *(TVectorD *)fileConstantsEff.FindObjectAny("efficiencyArray");
  TVectorD efficiencyErrArray = *(TVectorD *)fileConstantsEff.FindObjectAny("efficiencyErrArray");

  TFile fileScaleConstants(TString("../root_files/constants/")+tagDirConstants+TString("/")+fileScaleFactorConstants);
  TVectorD rhoDataMc    = *(TVectorD *)fileScaleConstants.FindObjectAny("scaleFactorArray");
  TVectorD rhoDataMcErr = *(TVectorD *)fileScaleConstants.FindObjectAny("scaleFactorErrArray");

  TFile fileConstantsAcc(TString("../root_files/constants/")+tagDirConstants+TString("/")+fileAcceptanceConstants);
  TVectorD acceptanceArray    = *(TVectorD *)fileConstantsAcc.FindObjectAny("acceptanceArray");
  TVectorD acceptanceErrArray = *(TVectorD *)fileConstantsAcc.FindObjectAny("acceptanceErrArray");

  TFile fileConstantsFsr(TString("../root_files/constants/")+tagDirConstants+TString("/")+fileFsrCorrectionConstants);
  TVectorD fsrCorrectionArray    = *(TVectorD *)fileConstantsFsr.FindObjectAny("fsrCorrectionArray");
  TVectorD fsrCorrectionErrArray = *(TVectorD *)fileConstantsFsr.FindObjectAny("fsrCorrectionErrArray");

  TFile fileConstantsFsrSansAcc(TString("../root_files/constants/")+tagDirConstants+TString("/")+fileFsrCorrectionSansAccConstants);
  TVectorD fsrCorrectionSansAccArray    = *(TVectorD *)fileConstantsFsrSansAcc.FindObjectAny("fsrCorrectionArray");
  TVectorD fsrCorrectionSansAccErrArray = *(TVectorD *)fileConstantsFsrSansAcc.FindObjectAny("fsrCorrectionErrArray");
  


  printf("\n\nLatex table of various corrections for PAS/paper\n");
  printf("               Acceptance, %%       Acc*Eff, %%          FSR corr, %%         FSR corr in acc, %%\n");
  for(int i=0; i<nMassBins; i++){

    double effFactor = efficiencyArray[i] * rhoDataMc[i];
//     double effErr = effFactor
//       * sqrt( efficiencyErrArray[i]*efficiencyErrArray[i]/efficiencyArray[i]/efficiencyArray[i]
// 	      + rhoDataMcErr[i]*rhoDataMcErr[i]/rhoDataMc[i]/rhoDataMc[i]);

    double accFactor = acceptanceArray[i];
    double accErr    = acceptanceErrArray[i];

    double accEff = accFactor * effFactor;
//     double accEffErr = accEff*sqrt( (effErr/effFactor)*(effErr/effFactor) 
// 				    + (accErr/accFactor)*(accErr/accFactor));

    double fsrFactor = fsrCorrectionArray[i];
    double fsrErr    = fsrCorrectionErrArray[i];

    double fsrFactorSansAcc = fsrCorrectionSansAccArray[i];
    double fsrErrSansAcc    = fsrCorrectionSansAccErrArray[i];

    printf("%4.0f-%4.0f &", massBinLimits[i],massBinLimits[i+1]);
    printf("  $%5.2f \\pm %4.2f$ &", accFactor*100, accErr*100);
    // For acceptance * efficiency, only acceptance contribution
    // to the error is printed. This is because we are reflecting
    // here MC-related errors. The error on the efficiency
    // is dominated by scale factor errors, related to low 
    // statistics of tag and probe in the data. Those errors
    // are explicitly stated in systematics.
    printf("  $%5.2f \\pm %4.2f$ &", accEff*100, accErr*100);
    printf("  $%6.2f \\pm %4.2f$ &", fsrFactor*100, fsrErr*100);
    printf("  $%6.2f \\pm %4.2f$ &", fsrFactorSansAcc*100, fsrErrSansAcc*100);
    printf("\n");
  }
}

void printRelativeSystErrors(){

  printf("\n\nLatex table of relative systematic errors  in percent for PAS/paper\n");
  printf("              Escale  &   Eff.   &   Bkg    &    Unfol &    Acc   &    sum   & AccTheory%%\n");
  for(int i=0; i<nMassBins; i++){

    // Factor out theory error from the total acceptance error
    double systAcceptanceExpRelative 
      = sqrt( systAcceptanceRelative[i]*systAcceptanceRelative[i]
	      - systAccTheoryRelative[i]*systAccTheoryRelative[i]);
    // The "sum" contains only the experimental systematic errors
    double sum = sqrt(systEscaleRelative[i]*systEscaleRelative[i]
		      + systEfficiency[i]*systEfficiency[i]
		      + systBackgrRelative[i]*systBackgrRelative[i]
		      + systUnfoldRelative[i]*systUnfoldRelative[i]
		      + systAcceptanceExpRelative*systAcceptanceExpRelative);
    printf("%4.0f-%4.0f &", massBinLimits[i],massBinLimits[i+1]);
    printf("   $%5.1f$ &  $%5.1f$ &  $%5.1f$ &  $%5.1f$ &  $%5.1f$  &  $%5.1f$ & $%5.1f",
	   100*systEscaleRelative[i], 
	   100*systEfficiency[i], 
	   100*systBackgrRelative[i], 
	   100*systUnfoldRelative[i],
           100*systAcceptanceExpRelative,
	   100*sum,
           100*systAccTheoryRelative[i]);
    printf("\n");
  }
}

  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////
  //Four plots of R-shape at the same picture

  void RShapePlot (TVectorD relCrossSection, TVectorD relCrossSectionStatErr, 
TVectorD relCrossSectionDET, TVectorD relCrossSectionStatErrDET, 
TVectorD relPostFsrCrossSection, TVectorD relPostFsrCrossSectionStatErr, 
TVectorD relPostFsrCrossSectionDET, TVectorD relPostFsrCrossSectionStatErrDET)
{
    
   TCanvas *c1 = new TCanvas("c1","c1",10,10,600,600);
   c1->SetGrid();
   c1->SetLogx(1);
   c1->SetLogy(1);
   c1->SetFillColor(0);
   // draw a frame to define the range
   TMultiGraph *mg = new TMultiGraph();

      // create first graph
   //Pre-FSR -All phase space
   double x[nMassBins];
   double ex[nMassBins];
   /*
   double y1[nMassBins];
   double ey1[nMassBins];
   double y2[nMassBins];
   double ey2[nMassBins];
   double y3[nMassBins];
   double ey3[nMassBins];
   double y4[nMassBins];
   double ey4[nMassBins];
   */
   double* y1;
   double* ey1;
   double* y2;
   double* ey2;
   double* y3;
   double* ey3;
   double* y4;
   double* ey4;
   for (int i=0; i<nMassBins; i++)
     {
       x[i]=(massBinLimits[i+1]+massBinLimits[i])/2;
       ex[i]=(massBinLimits[i+1]-massBinLimits[i])/2;
     }

   y1=relCrossSection.GetMatrixArray();
   ey1=relCrossSectionStatErr.GetMatrixArray();
   y2=relCrossSectionDET.GetMatrixArray();
   ey2=relCrossSectionStatErrDET.GetMatrixArray();
   y3=relPostFsrCrossSection.GetMatrixArray();
   ey3=relPostFsrCrossSectionStatErr.GetMatrixArray();
   y4=relPostFsrCrossSectionDET.GetMatrixArray();
   ey4=relPostFsrCrossSectionStatErrDET.GetMatrixArray();

   const Int_t n = nMassBins;
   
   TGraphErrors *gr1 = new TGraphErrors(n,x,y1,ex,ey1);
   gr1->SetName("PreFSR (d#sigma /dM)/ (d#sigma /dM)_{z}");
   gr1->SetTitle("Pre-FSR (d#sigma /dM)/ (d#sigma /dM)_{z}");
   gr1->SetFillColor(0);
   gr1->SetMarkerColor(kBlack);
   gr1->SetMarkerStyle(20);
   gr1->SetMarkerSize(1.0);

   TGraphErrors *gr2 = new TGraphErrors(n,x,y2,ex,ey2);
   gr2->SetName("PreFSR DET (d#sigma /dM)/ (d#sigma /dM)_{z}");
   gr2->SetTitle("Pre-FSR DET (d#sigma /dM)/ (d#sigma /dM)_{z}");
   gr2->SetFillColor(0);
   gr2->SetMarkerColor(kBlue);
   gr2->SetMarkerStyle(20);
   gr2->SetMarkerSize(1.0);
   gr2->SetLineColor(kBlue);

   TGraphErrors *gr3 = new TGraphErrors(n,x,y3,ex,ey3);
   gr3->SetName("PostFSR (d#sigma /dM)/ (d#sigma /dM)_{z}");
   gr3->SetTitle("Post-FSR (d#sigma /dM)/ (d#sigma /dM)_{z}");
   gr3->SetFillColor(0);
   gr3->SetMarkerColor(kRed);
   gr3->SetMarkerStyle(20);
   gr3->SetMarkerSize(1.0);
   gr3->SetLineColor(kRed);

   TGraphErrors *gr4 = new TGraphErrors(n,x,y4,ex,ey4);
   gr4->SetName("PostFSR DET(d#sigma /dM)/ (d#sigma /dM)_{z}");
   gr4->SetTitle("Post-FSR DEt(d#sigma /dM)/ (d#sigma /dM)_{z}");
   gr4->SetFillColor(0);
   gr4->SetMarkerColor(kGreen);
   gr4->SetMarkerStyle(20);
   gr4->SetMarkerSize(1.0);
   gr4->SetLineColor(kGreen);

   mg->Add(gr4);
   mg->Add(gr3);
   mg->Add(gr2);
   mg->Add(gr1);
   mg->Draw("ap");
  TAxis* yax=mg->GetYaxis();
  yax->SetRangeUser(5e-6,2);
   //mg->GetXaxis()->SetTitle("M_{ee}");
   //mg->GetYaxis()->SetTitle("R-shape");
   //mg->SetName("(d#sigma /dM)/ (d#sigma /dM)_{z}");

   TLegend *leg = new TLegend(.50,.45,.90,.85);
   leg->AddEntry(gr1,"Pre FSR All Phase Space");
   leg->AddEntry(gr2,"Pre FSR Detector Phase space");
   leg->AddEntry(gr3,"Post FSR All Phase Space");
   leg->AddEntry(gr4,"Post FSR Detector Phase space");
   leg->Draw();

}

  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////
void getNormBinRange(int &firstNormBin, int &lastNormBin){

  firstNormBin = -1;
  lastNormBin = -1;
 
  for(int i=0; i<=DYTools::nMassBins; i++){
    if(DYTools::massBinLimits[i] == lowZMass)
      firstNormBin = i;
    if(DYTools::massBinLimits[i] == highZMass)
      lastNormBin = i-1;
  }
  
  if(firstNormBin == -1 || lastNormBin == -1){
    printf("\nERROR: normalization limits are misaligned with mass binning!\n\n");
    assert(0);
  }
  printf("\nCross section normalization is to the bins %d - %d from %5.1f to %5.1f GeV\n", 
	 firstNormBin, lastNormBin, 
	 DYTools::massBinLimits[firstNormBin], DYTools::massBinLimits[lastNormBin+1]);

  return;
}

