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

#include "../Include/DYTools.hh"
#include "../Include/UnfoldingTools.hh"

// This global constants will be filled from 
// the configuration file. This is not great C++...
TString tagDirYields = "";
TString tagDirConstants = "";
Double_t lumi = 0;

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
TVectorD systEscaleBeforeUnfolding(DYTools::nMassBins);

// Relative error values. These are meant to be AFTER unfolding.
TVectorD systBackgrRelative(DYTools::nMassBins);
TVectorD systEscaleRelative(DYTools::nMassBins);
TVectorD systUnfoldRelative(DYTools::nMassBins);

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
  systEscaleBeforeUnfolding = 0;
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

  // We also load systematic errors relevant at this point
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

  // Prepare output yields and errors
  for(int i=0; i<nMassBins; i++){
    v[i] = YieldsSignal[i];
    vErr1[i] = YieldsSignalErr[i];
    // Background systematics should be already in, add 
    // energy scale systematics
    double escaleError = (escaleSystematicsPercent[i]/100.0)*v[i];
    vErr2[i] = sqrt(YieldsSignalSystErr[i] * YieldsSignalSystErr[i]
		    + escaleError*escaleError);
    systBackgrBeforeUnfolding[i] = YieldsSignalSystErr[i];
    systEscaleBeforeUnfolding[i] = escaleError;
  } 

  fileYields.Close();
  fileEscaleSystematics.Close();
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

  // Second, propagate separately several systematic error components.
  // These are already included in the total systematic error above in vinSystErr,
  // however we do it separately so that we can quote the breakdown in the
  // table of systematic errors
  TVectorD systBackgr(vin.GetNoElements());
  TVectorD systEscale(vin.GetNoElements());
  unfolding::propagateErrorThroughUnfolding(systBackgrBeforeUnfolding, systBackgr, fullUnfoldingConstFileName);
  unfolding::propagateErrorThroughUnfolding(systEscaleBeforeUnfolding, systEscale, fullUnfoldingConstFileName);

  // Pool together the unfolding systematics and add it to the total systematics
  TVectorD systUnfolding(vin.GetNoElements());
  unfolding::calculateTotalUnfoldingSystError(vin, systUnfolding, 
					      fullUnfoldingConstFileName, 
					      fullUnfoldingErrFileName);
  // Add unfolding systematics to the total systematic error
  for(int i=0; i<DYTools::nMassBins; i++){
    voutSystErr[i] = sqrt( voutSystErr[i]*voutSystErr[i] + systUnfolding[i]*systUnfolding[i]);
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

  // Check that the binning is consistent
  bool checkResult = true;
  if( acceptanceArray.GetNoElements() != nMassBins) checkResult = false;
  if( !checkResult ){
    printf("Acceptance: ERROR: inconsistent binning in the inputs\n");
    assert(0);
  }else
    printf("Acceptance: Binning in the inputs is consistent\n");

  // Apply the correction
  TVectorD systErrorPropagated(nMassBins);
  TVectorD systErrorAdditional(nMassBins);
  systErrorAdditional = 0;
  for(int i=0; i<nMassBins; i++){
    double accFactor = acceptanceArray[i];
    double accErr    = acceptanceErrArray[i];
    vout[i]        = vin[i] / accFactor;
    voutStatErr[i] = vinStatErr[i] / accFactor;
    systErrorPropagated[i] = vinSystErr[i]/accFactor;
    systErrorAdditional[i] = (vin[i]/accFactor)*(accErr/accFactor);
    voutSystErr[i] = sqrt(systErrorPropagated[i]*systErrorPropagated[i] + systErrorAdditional[i]*systErrorAdditional[i]);
  }

  printf("\nAcceptance: Results for the data, yields:\n");
  printf("                eff-corrected             acc. factors,%%      acceptance-corrected      syst-err-acc, %%\n");
  for(int i=0; i<nMassBins; i++){
    printf("%4.0f-%4.0f   %8.1f +- %7.1f +- %6.1f   %8.3f +- %7.3f  %9.1f +- %8.1f +- %6.1f      %6.2f\n",
	   massBinLimits[i],massBinLimits[i+1],
	   vin[i], vinStatErr[i], vinSystErr[i],
	   acceptanceArray[i]*100, acceptanceErrArray[i]*100, 
	   vout[i], voutStatErr[i], voutSystErr[i],
	   systErrorAdditional[i]*100/vout[i]);
  }

  fileConstants.Close();
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
  printf("                acc-corrected             fsr. factors         fsr-corrected                  sys-err-fsr\n");
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
  printf("                acc-corrected             fsr. factors         fsr-corrected                  sys-err-fsr\n");
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

  printf("\nPre FSR cross sections: :\n");
  printf("                    absolute                   normalized +- stat +- sys (total)\n");
  for(int i=0; i<nMassBins; i++){
    printf("%4.0f-%4.0f   %9.1f +- %8.1f +- %6.1f   %11.6f +- %10.6f +- %10.6f  ( %10.6f )\n",
	   massBinLimits[i],massBinLimits[i+1],
	   vout[i], voutStatErr[i], voutSystErr[i],
	   voutNorm[i], voutNormStatErr[i], voutNormSystErr[i],
	   sqrt(voutNormStatErr[i] * voutNormStatErr[i] + voutNormSystErr[i] * voutNormSystErr[i]) );
  }
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
  printf("                    absolute                   normalized +-stat +-sys (total)\n");
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
  printf("                    absolute                   normalized +-stat +-sys (total)\n");
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
  printf("                    absolute                   normalized +-stat +-sys (total)\n");
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
  printf("               obs-bg    unfolded    eff-corrected    acc-corrected   fsr-corrected\n");
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
  printf("               Acceptance, %%    Acc*Eff, %%     FSR corr, %%       FSR corr in acc, %%\n");
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
  printf("              escale & efficiency & backgrounds & unfolding & sum\n");
  for(int i=0; i<nMassBins; i++){
    double sum = sqrt(systEscaleRelative[i]*systEscaleRelative[i]
		      + systEfficiency[i]*systEfficiency[i]
		      + systBackgrRelative[i]*systBackgrRelative[i]
		      + systUnfoldRelative[i]*systUnfoldRelative[i]);
    printf("%4.0f-%4.0f &", massBinLimits[i],massBinLimits[i+1]);
    printf("   $%5.1f$ &  $%5.1f$ &  $%5.1f$ &  $%5.1f$ &  $%5.1f$", 
	   100*systEscaleRelative[i], 
	   100*systEfficiency[i], 
	   100*systBackgrRelative[i], 
	   100*systUnfoldRelative[i],
	   100*sum);
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

