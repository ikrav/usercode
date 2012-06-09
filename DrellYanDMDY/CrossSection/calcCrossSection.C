#include "TROOT.h"
#include "TSystem.h"
#include "TVectorD.h"
#include "TMatrixD.h"
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
#include "../Include/TriggerSelection.hh"
#include "../Include/CPlot.hh"
#include "../Include/plotFunctions.hh"

using std::string;
using std::stringstream;

template<class T> T SQR(const T &x) { return x*x; }

// This global constants will be filled from 
// the configuration file. This is not great C++...
TString tagDirYields = "";
TString tagDirConstants = "";
Double_t lumi = 0;
const int nMassBinTh=518;

const int includeEScaleSystematics=0;
const int includeUnfoldingSystematics=0;


const int printFSRcorrectionTable=0;
const int printEfficiencyTable=0;
const int printAcceptanceTable=0;
const int printPreFSRCrossSectionTable=0;
const int printPreFSR_DET_ShapeTable=0;
const int printPostFSRCrossSectionTable=0;
const int printPostFSR_DET_ShapeTable=0;

const int printAllCorrectionTable=0;
const int printRelativeSystErrTable=0;
const int callPrintTableForNotes=0;

const TString fileEnd(analysisTag + TString(".root"));

const TString fileDataYields        (TString("yields_bg-subtracted") + fileEnd);
const TString fileMcReferenceYields (TString("yields_MC_unfolding_reference_") + fileEnd);
// This file contains unfolding matrix 
const TString fileUnfoldingConstants(TString("unfolding_constants") + fileEnd);
// Contains relative unfolding systematic errors
const TString fileUnfoldingErrors(TString("unfolding_systematics") + fileEnd);
// Contains relative escale systematic errors
const TString fileEscaleErrors(TString("escale_systematics") + fileEnd);
const TString fileEfficiencyConstants(TString("event_efficiency_constants") + fileEnd);
TString fileScaleFactorConstants(TString("scale_factors_") + fileEnd);
const TString fileAcceptanceConstants(TString("acceptance_constants") + fileEnd);
const TString fileAcceptanceSystematics(TString("theoretical_uncertainties.root"));
const TString fileAcceptanceFSRSystematics(TString("acceptance_FSR_systematics") + fileEnd);
const TString fileFsrCorrectionConstants(TString("fsr_constants_") + fileEnd);
const TString fileFsrCorrectionSansAccConstants(TString("fsr_constants_") + analysisTag + TString("_sans_acc.root"));

// Forward declarations
//-----------------------------------------------------------------

void readData(TMatrixD &v, TMatrixD &vErr1, TMatrixD &vErr2);

void  applyUnfolding(TMatrixD &vin, TMatrixD &vinStatErr, TMatrixD &vinSystErr,
             TMatrixD &vout, TMatrixD &voutStatErr, TMatrixD &voutSystErr);
void  applyUnfoldingToMc(); //TString fullUnfoldingConstFileName, TString fullMcRefYieldsFileName);
void  efficiencyCorrection(TMatrixD &vin, TMatrixD &vinStatErr, TMatrixD &vinSystErr,
             TMatrixD &vout, TMatrixD &voutStatErr, TMatrixD &voutSystErr);
void  acceptanceCorrection(TMatrixD &vin, TMatrixD &vinStatErr, TMatrixD &vinSystErr,
             TMatrixD &vout, TMatrixD &voutStatErr, TMatrixD &voutSystErr);
void  fsrCorrectionBase(TMatrixD &vin, TMatrixD &vinStatErr, TMatrixD &vinSystErr,
			TMatrixD &vout, TMatrixD &voutStatErr, TMatrixD &voutSystErr, 
			const TString &fname, const TString &title);
void  fsrCorrection(TMatrixD &vin, TMatrixD &vinStatErr, TMatrixD &vinSystErr,
             TMatrixD &vout, TMatrixD &voutStatErr, TMatrixD &voutSystErr);
void  fsrCorrectionSansAcceptance(TMatrixD &vin, TMatrixD &vinStatErr, TMatrixD &vinSystErr,
				  TMatrixD &vout, TMatrixD &voutStatErr, TMatrixD &voutSystErr);

void  crossSections(TMatrixD &vin, TMatrixD &vinStatErr, TMatrixD &vinSystErr,
		    TMatrixD &vout, TMatrixD &voutStatErr, TMatrixD &voutSystErr,
		    TMatrixD &voutNorm, TMatrixD &voutNormStatErr, TMatrixD &voutNormSystErr,
		    const TriggerSelection &triggers);

void  crossSectionsDET(TMatrixD &vin, TMatrixD &vinStatErr, TMatrixD &vinSystErr,
		       TMatrixD &vout, TMatrixD &voutStatErr, TMatrixD &voutSystErr,
		       TMatrixD &voutNorm, TMatrixD &voutNormStatErr, TMatrixD &voutNormSystErr);

void  postFsrCrossSections(TMatrixD &vin, TMatrixD &vinStatErr, TMatrixD &vinSystErr,
			   TMatrixD &vout, TMatrixD &voutStatErr, TMatrixD &voutSystErr,
			   TMatrixD &voutNorm, TMatrixD &voutNormStatErr, TMatrixD &voutNormSystErr);

void  postFsrCrossSectionsDET(TMatrixD &vin, TMatrixD &vinStatErr, TMatrixD &vinSystErr,
			      TMatrixD &vout, TMatrixD &voutStatErr, TMatrixD &voutSystErr,
			      TMatrixD &voutNorm, TMatrixD &voutNormStatErr, TMatrixD &voutNormSystErr);

void printTableForNotes(const TMatrixD &obs, const TMatrixD &obsErr, 
			const TMatrixD &unf, const TMatrixD &unfErr,
			const TMatrixD &ecor, const TMatrixD &ecorErr,
			const TMatrixD &acor, const TMatrixD &acorErr,
			const TMatrixD &fcor, const TMatrixD &fcorErr);

void printAllCorrections();
void printRelativeSystErrors();

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
//Four plots of R-shape at the same picture

//void RShapePlot
//(TMatrixD relCrossSection, TMatrixD relCrossSectionStatErr, 
 //TMatrixD relCrossSectionDET, TMatrixD relCrossSectionDETStatErr, 
 //TMatrixD relPostFsrCrossSection, TMatrixD relPostFsrCrossSectionStatErr, 
 //TMatrixD relPostFsrCrossSectionDET, TMatrixD relPostFsrCrossSectionDETStatErr);

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////


const double lowZMass = 60.0;
const double highZMass = 120.0;
void getNormBinRange(int &firstNormBin, int &lastNormBin);

// Some global arrays for convenience.
// These will contain errors on R
// (after unfolding, etc).

int nMaxYBins=DYTools::findMaxYBins();

// Absolute error values at the point just before unfolding
TMatrixD systBackgrBeforeUnfolding(DYTools::nMassBins,nMaxYBins);

// Relative error values. These are meant to be AFTER unfolding.
TMatrixD systBackgrRelative(DYTools::nMassBins,nMaxYBins);
TMatrixD systEscaleRelative(DYTools::nMassBins,nMaxYBins);
TMatrixD systUnfoldRelative(DYTools::nMassBins,nMaxYBins);
TMatrixD systAccTheoryRelative(DYTools::nMassBins,nMaxYBins);
TMatrixD systAcceptanceRelative(DYTools::nMassBins,nMaxYBins);

TMatrixD systEfficiency(DYTools::nMassBins,nMaxYBins);
TMatrixD systOthers(DYTools::nMassBins,nMaxYBins);

// ---------------------------------------------------------------
// Main function
// ---------------------------------------------------------------

void calcCrossSection(const TString conf="../config_files/xsecCalc.conf"){

  // Read from configuration file only the location of the root files
  ifstream ifs;
  ifs.open(conf.Data());
  assert(ifs.is_open());
  string line;
  int state = 0;
  TString triggerSetString;
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
      state++;
    } else if (state==3) {
      triggerSetString = TString(line);
      break;
    }
  }
  ifs.close();

  // Construct the trigger object
  TriggerSelection triggers(triggerSetString, true, 0);
  assert ( triggers.isDefined() );
  // update the name of the file with per-event scale factors
  fileScaleFactorConstants.Insert(fileScaleFactorConstants.Index(".root"),
				  TString("_") + triggers.triggerConditionsName());
  std::cout << "using fileScaleFactorConstants=" << fileScaleFactorConstants << "\n";

  // Last, do a closure test on MC
  applyUnfoldingToMc(); //fullUnfoldingConstFileName, fullMcRefYieldsFileName);


  TMatrixD signalYields(nMassBins,nMaxYBins);
  TMatrixD signalYieldsStatErr(nMassBins,nMaxYBins);
  TMatrixD signalYieldsSystErr(nMassBins,nMaxYBins);
  
  TMatrixD unfoldedYields(nMassBins,nMaxYBins);
  TMatrixD unfoldedYieldsStatErr(nMassBins,nMaxYBins);
  TMatrixD unfoldedYieldsSystErr(nMassBins,nMaxYBins);
  
  TMatrixD effCorrectedYields(nMassBins,nMaxYBins);
  TMatrixD effCorrectedYieldsStatErr(nMassBins,nMaxYBins);
  TMatrixD effCorrectedYieldsSystErr(nMassBins,nMaxYBins);
  
  TMatrixD accCorrectedYields(nMassBins,nMaxYBins);
  TMatrixD accCorrectedYieldsStatErr(nMassBins,nMaxYBins);
  TMatrixD accCorrectedYieldsSystErr(nMassBins,nMaxYBins);
  
  TMatrixD preFsrYields(nMassBins,nMaxYBins);
  TMatrixD preFsrYieldsStatErr(nMassBins,nMaxYBins);
  TMatrixD preFsrYieldsSystErr(nMassBins,nMaxYBins);
  
  TMatrixD preFsrSansAccYields(nMassBins,nMaxYBins);
  TMatrixD preFsrSansAccYieldsStatErr(nMassBins,nMaxYBins);
  TMatrixD preFsrSansAccYieldsSystErr(nMassBins,nMaxYBins);
  
  TMatrixD absCrossSection(nMassBins,nMaxYBins);
  TMatrixD absCrossSectionStatErr(nMassBins,nMaxYBins);
  TMatrixD absCrossSectionSystErr(nMassBins,nMaxYBins);
  
  TMatrixD absCrossSectionDET(nMassBins,nMaxYBins);
  TMatrixD absCrossSectionStatErrDET(nMassBins,nMaxYBins);
  TMatrixD absCrossSectionSystErrDET(nMassBins,nMaxYBins);
  
  TMatrixD relCrossSection(nMassBins,nMaxYBins);
  TMatrixD relCrossSectionStatErr(nMassBins,nMaxYBins);
  TMatrixD relCrossSectionSystErr(nMassBins,nMaxYBins);

  TMatrixD relCrossSectionDET(nMassBins,nMaxYBins);
  TMatrixD relCrossSectionStatErrDET(nMassBins,nMaxYBins);
  TMatrixD relCrossSectionSystErrDET(nMassBins,nMaxYBins);

  TMatrixD absPostFsrCrossSection(nMassBins,nMaxYBins);
  TMatrixD absPostFsrCrossSectionStatErr(nMassBins,nMaxYBins);
  TMatrixD absPostFsrCrossSectionSystErr(nMassBins,nMaxYBins);
  
  TMatrixD absPostFsrCrossSectionDET(nMassBins,nMaxYBins);
  TMatrixD absPostFsrCrossSectionStatErrDET(nMassBins,nMaxYBins);
  TMatrixD absPostFsrCrossSectionSystErrDET(nMassBins,nMaxYBins);
  
  TMatrixD relPostFsrCrossSection(nMassBins,nMaxYBins);
  TMatrixD relPostFsrCrossSectionStatErr(nMassBins,nMaxYBins);
  TMatrixD relPostFsrCrossSectionSystErr(nMassBins,nMaxYBins);

  TMatrixD relPostFsrCrossSectionDET(nMassBins,nMaxYBins);
  TMatrixD relPostFsrCrossSectionStatErrDET(nMassBins,nMaxYBins);
  TMatrixD relPostFsrCrossSectionSystErrDET(nMassBins,nMaxYBins);

  systBackgrBeforeUnfolding = 0;
  systEfficiency = 0;
  systOthers = 0;

  // Read data yields from file (background subtraction is already done)
  std::cout << "read bg-subtracted data yields\n";
  readData(signalYields, signalYieldsStatErr, signalYieldsSystErr);

  // Apply unfolding
  std::cout << "apply unfolding\n";
  applyUnfolding(signalYields, signalYieldsStatErr, signalYieldsSystErr,
		 unfoldedYields, unfoldedYieldsStatErr, unfoldedYieldsSystErr);

  // Apply efficiency correction
  std::cout << "efficiency correction\n";
  efficiencyCorrection(unfoldedYields, unfoldedYieldsStatErr, unfoldedYieldsSystErr,
		       effCorrectedYields, effCorrectedYieldsStatErr, effCorrectedYieldsSystErr);

  // Acceptance corrections
  std::cout << "acceptance correction\n";
  acceptanceCorrection(effCorrectedYields, effCorrectedYieldsStatErr, effCorrectedYieldsSystErr,
		       accCorrectedYields, accCorrectedYieldsStatErr, accCorrectedYieldsSystErr);
  
  // FSR corrections
  std::cout << "fsr correction\n";
  fsrCorrection(accCorrectedYields, accCorrectedYieldsStatErr, accCorrectedYieldsSystErr,
		preFsrYields, preFsrYieldsStatErr, preFsrYieldsSystErr);

  // Alternative (for DET shapes): all corrections except for acceptance correction 
  fsrCorrectionSansAcceptance(effCorrectedYields, effCorrectedYieldsStatErr, effCorrectedYieldsSystErr,
			      preFsrSansAccYields, preFsrSansAccYieldsStatErr, preFsrSansAccYieldsSystErr);
 
  // Calculate absolute and relative cross-sections
  std::cout << "absolute and relative cross sections" << std::endl;
  std::cout << "1. crossSections" << std::endl;
  crossSections(preFsrYields, preFsrYieldsStatErr, preFsrYieldsSystErr,
		absCrossSection, absCrossSectionStatErr, absCrossSectionSystErr,
		relCrossSection, relCrossSectionStatErr, relCrossSectionSystErr,
		triggers);

  // Calculate absolute and relative cross-sections DET (shapes, no Acc, but with FSR)
  std::cout << "2. crossSectionsDET" << std::endl;
  crossSectionsDET(preFsrSansAccYields, preFsrSansAccYieldsStatErr, preFsrSansAccYieldsSystErr,
		   absCrossSectionDET, absCrossSectionStatErrDET, absCrossSectionSystErrDET,
		   relCrossSectionDET, relCrossSectionStatErrDET, relCrossSectionSystErrDET);

  // Also, calculate absolute and relative cross-sections for post-FSR stage
  std::cout << "3. postFsrCrossSection" << std::endl;
  postFsrCrossSections(accCorrectedYields, accCorrectedYieldsStatErr, accCorrectedYieldsSystErr,
		absPostFsrCrossSection, absPostFsrCrossSectionStatErr, absPostFsrCrossSectionSystErr,
		relPostFsrCrossSection, relPostFsrCrossSectionStatErr, relPostFsrCrossSectionSystErr);

  // calculate absolute and relative cross-sections DET (shapes, no Acc, no FSR corrections)
  std::cout << "4. postFsrCrossSectionsDET" << std::endl;
  postFsrCrossSectionsDET(effCorrectedYields, effCorrectedYieldsStatErr, effCorrectedYieldsSystErr,
			  absPostFsrCrossSectionDET, absPostFsrCrossSectionStatErrDET, absPostFsrCrossSectionSystErrDET,
			  relPostFsrCrossSectionDET, relPostFsrCrossSectionStatErrDET, relPostFsrCrossSectionSystErrDET);

  // Output
  std::cout << "\nprint output\n" << std::endl;
  if (printAllCorrectionTable) printAllCorrections();
  if (printRelativeSystErrTable) printRelativeSystErrors();

  if (callPrintTableForNotes) {
  std::cout << "\nprintTablesForNotes\n" << std::endl;
  printTableForNotes(signalYields      , signalYieldsStatErr, 
		     unfoldedYields    , unfoldedYieldsStatErr,
		     effCorrectedYields, effCorrectedYieldsStatErr,
		     accCorrectedYields, accCorrectedYieldsStatErr,
		     preFsrYields      , preFsrYieldsStatErr);
 
   //draw and save plots

/*
   for(int i=0; i<DYTools::nMassBins; i++)
    for (int j=0; j<nYBins[i]; j++)
      {
        relCrossSection(i,j)=i+j; 
        relCrossSectionDET(i,j)=i+j;
        relPostFsrCrossSection(i,j)=i+j;
        relPostFsrCrossSectionDET(i,j)=i+j;
      }
*/
   CPlot::sOutDir="plots" + DYTools::analysisTag;

 
   PlotMatrixVariousBinning(relCrossSection, "relative_CS", "LEGO2",0, CPlot::sOutDir);
   PlotMatrixVariousBinning(relCrossSectionDET, "relative_CS_DET", "LEGO2",0, CPlot::sOutDir);
   PlotMatrixVariousBinning(relPostFsrCrossSection, "relative_postFSR_CS", "LEGO2",0, CPlot::sOutDir);
   PlotMatrixVariousBinning(relPostFsrCrossSectionDET, "relative_postFSR_CS_DET", "LEGO2",0, CPlot::sOutDir);




  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////
  //Four plots of R-shape at the same picture

  std::cout << "\nRShapePlot" << std::endl;
  RShapePlot
    (relCrossSection, relCrossSectionStatErr, 
     relCrossSectionDET, relCrossSectionStatErrDET, 
     relPostFsrCrossSection, relPostFsrCrossSectionStatErr, 
     relPostFsrCrossSectionDET, relPostFsrCrossSectionStatErrDET);
  
  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////

  std::cout << "Leaving" << std::endl;
  return;
}


//-----------------------------------------------------------------
// Read data
//-----------------------------------------------------------------
void readData(TMatrixD &v, TMatrixD &vErr1, TMatrixD &vErr2){

  printf("Load data yields\n"); fflush(stdout);
  TString fileYieldsFullName=TString("../root_files/yields/")+tagDirYields+TString("/")+fileDataYields;
  TFile fileYields   (fileYieldsFullName);
  TMatrixD *YieldsSignalPtr       = (TMatrixD *)fileYields.FindObjectAny("YieldsSignal");
  TMatrixD *YieldsSignalErrPtr    = (TMatrixD *)fileYields.FindObjectAny("YieldsSignalErr");
  TMatrixD *YieldsSignalSystErrPtr= (TMatrixD *)fileYields.FindObjectAny("YieldsSignalSystErr");
  TVectorD *MassBinLimitsForYieldsPtr = (TVectorD *)fileYields.FindObjectAny("massBinning");
  TVectorD *YBinCountsForYieldsPtr = (TVectorD *)fileYields.FindObjectAny("rapidityCounts");

  if (!YieldsSignalPtr || !YieldsSignalErrPtr || !YieldsSignalSystErrPtr ||
      !MassBinLimitsForYieldsPtr || !YBinCountsForYieldsPtr) {
    std::cout << "failed to get at least one of the required objects from file <" << fileDataYields << ">\n";
    assert(0);
  }

  TMatrixD YieldsSignal = *YieldsSignalPtr;
  TMatrixD YieldsSignalErr = *YieldsSignalErrPtr;
  TMatrixD YieldsSignalSystErr = *YieldsSignalSystErrPtr;
  TVectorD MassBinLimitsForYields = *MassBinLimitsForYieldsPtr;
  TVectorD YBinCountsForYields = *YBinCountsForYieldsPtr;

  // Check that the binning is consistent
  bool checkResult = true;
  checkResult=unfolding::checkRangesMY(v,"v in readData");
  if (checkResult) checkResult=unfolding::checkBinningRanges(MassBinLimitsForYields,YBinCountsForYields,fileYieldsFullName);
  if (checkResult) checkResult=unfolding::checkRangesMY(YieldsSignal,"YieldsSignal");
  if( !checkResult ){
    printf("ERROR: inconsistent binning in the inputs\n");
    assert(0);
  }else
    printf("readData: Binning in the inputs is consistent\n");

  // Prepare output yields and errors
  for(int i=0; i<nMassBins; i++){
    for (int yi=0; yi<nYBins[i]; yi++) {
      v[i][yi] = YieldsSignal[i][yi];
      vErr1[i][yi] = YieldsSignalErr[i][yi];
      // Background systematics should be already in, add 
      // energy scale systematics
      vErr2[i][yi] = YieldsSignalSystErr[i][yi];
      systBackgrBeforeUnfolding[i][yi] = YieldsSignalSystErr[i][yi];
    }
  }

  fileYields.Close();
  return;
}

//-----------------------------------------------------------------
// Unfold
//-----------------------------------------------------------------
void  applyUnfolding(TMatrixD &vinM, TMatrixD &vinStatErrM, TMatrixD &vinSystErrM,
             TMatrixD &voutM, TMatrixD &voutStatErrM, TMatrixD &voutSystErrM)
{

  // Construct file names
  TString fullUnfoldingConstFileName = TString("../root_files/constants/")+tagDirConstants+TString("/")+fileUnfoldingConstants;
  TString fullUnfoldingErrFileName   = TString("../root_files/systematics/")+tagDirConstants+TString("/")+fileUnfoldingErrors;

  int nUnfoldingBins = DYTools::getTotalNumberOfBins();
  TVectorD vin(nUnfoldingBins),vinStatErr(nUnfoldingBins),vinSystErr(nUnfoldingBins);
  TVectorD vout(nUnfoldingBins),voutStatErr(nUnfoldingBins),voutSystErr(nUnfoldingBins);

  // First, propagate through unfolding the signal yields with stat and syst errors
  assert(unfolding::unfold(vinM, voutM, fullUnfoldingConstFileName, vin, vout)==1);
  assert(unfolding::propagateErrorThroughUnfolding(vinStatErrM,voutStatErrM, fullUnfoldingConstFileName, vinStatErr, voutStatErr)==1);
  assert(unfolding::propagateErrorThroughUnfolding(vinSystErrM,voutSystErrM, fullUnfoldingConstFileName, vinSystErr, voutSystErr)==1);

  // Second, propagate separately systematic error components that need it.
  // These are already included in the total systematic error above in vinSystErr,
  // however we do it separately so that we can quote the breakdown in the
  // table of systematic errors
  TMatrixD systBackgrM(DYTools::nMassBins,nMaxYBins);
  TVectorD systBackgrBeforeUnfoldingV(nUnfoldingBins), systBackgrV(nUnfoldingBins);
  unfolding::propagateErrorThroughUnfolding(systBackgrBeforeUnfolding, systBackgrM, fullUnfoldingConstFileName, systBackgrBeforeUnfoldingV,systBackgrV);

  // The electron energy scale systematics that is loaded here
  // is estimated on the unfolded yields. So we read it in at this time
  // and will add to the outgoing total systematic below in this function.
  TVectorD systEscaleV(nUnfoldingBins);
  systEscaleV=0;
  if (includeEScaleSystematics) {
  TString fullEscaleSystErrorsFileName = TString("../root_files/systematics/")
    +tagDirYields+TString("/")
    +fileEscaleErrors;
  TFile fileEscaleSystematics(fullEscaleSystErrorsFileName);
  if( ! fileEscaleSystematics.IsOpen()){
    printf("ERROR: required file with escale errors %s is not found!\n",
	   fullEscaleSystErrorsFileName.Data());
    assert(0);
  }
  assert(unfolding::checkBinningArrays(fileEscaleSystematics));
  TMatrixD *escaleSystematicsPercentPtr
    = (TMatrixD *)fileEscaleSystematics.FindObjectAny("escaleSystPercent");
  assert(escaleSystematicsPercentPtr);
  assert(unfolding::checkRangesMY(*escaleSystematicsPercentPtr,"escaleSystPercent"));
  TMatrixD escaleSystematicsPercent = *escaleSystematicsPercentPtr;
  
  systEscaleV=0;
  for(int i=0; i<nMassBins; i++){
    for (int yi=0; yi<nYBins[i]; ++yi) {
      int idx=DYTools::findIndexFlat(i,yi);
      systEscaleV[idx] = (escaleSystematicsPercent[i][yi]/100.0) * vout[idx];
    }
  }
  }

  // Pool together the unfolding systematics and add it to the total systematics
  TVectorD systUnfoldingV(nUnfoldingBins);
  systUnfoldingV=0;
  if (includeUnfoldingSystematics) {
    unfolding::calculateTotalUnfoldingSystErrorFlat(vin, systUnfoldingV, 
						    fullUnfoldingConstFileName, 
						    fullUnfoldingErrFileName);
  }

  // Add unfolding and escale systematics to the total systematic error
  for(int i=0; i<nUnfoldingBins; i++){
    voutSystErr[i] = sqrt( voutSystErr[i]*voutSystErr[i] 
			   + systUnfoldingV[i]*systUnfoldingV[i]
			   + systEscaleV[i]*systEscaleV[i]);
  }

  // After propagating through unfolding all errors that we had on yields before 
  // unfolding we can compute the relative errors of each kind. While unfolding
  // changes relative errors, all subsequent manipulations do not, so we can 
  // save the errors here.
  for(int i=0; i<DYTools::nMassBins; i++){
    for (int yi=0; yi<nYBins[i]; yi++) {
      int idx=findIndexFlat(i,yi);
      systBackgrRelative[i][yi] = systBackgrV[idx]/vout[idx];
      systEscaleRelative[i][yi] = systEscaleV[idx]/vout[idx];
      systUnfoldRelative[i][yi] = systUnfoldingV[idx]/vout[idx];
    }
  }

  // Print the result
  if (0) {
    printf("\nUNFOLD: Results for the data, yields:\n");
    printf(" mass        |y| range      yields observed        after unfolding            \n");
    for(int i=0; i<DYTools::nMassBins; i++){
      double *rapidityBinLimits=DYTools::getYBinLimits(i);
      for (int yi=0; yi<nYBins[i]; ++yi) {
      int idx=findIndexFlat(i,yi);
      printf("%4.0f-%4.0f  %4.2lf-%4.2lf  %8.1f +- %6.1f +- %5.1f       %8.1f +- %7.1f +- %6.1f\n",
	     DYTools::massBinLimits[i], DYTools::massBinLimits[i+1],
	     rapidityBinLimits[yi], rapidityBinLimits[yi+1],
	     vin[idx], vinStatErr[idx], vinSystErr[idx],
	     vout[idx], voutStatErr[idx], voutSystErr[idx]);
      }
      delete rapidityBinLimits;
    }
  }


  // added to the main sequence
  // Last, do a closure test on MC
  //applyUnfoldingToMc(); //fullUnfoldingConstFileName, fullMcRefYieldsFileName);

  return;
}

void applyUnfoldingToMc() { //TString fullUnfoldingConstFileName, TString fullMcRefYieldsFileName)
  TString fullUnfoldingConstFileName = TString("../root_files/constants/")+tagDirConstants+TString("/")+fileUnfoldingConstants;
  TString fullMcRefYieldsFileName    = TString("../root_files/constants/")+tagDirConstants+TString("/")+fileMcReferenceYields;

  printf("Load MC reference yields\n"); fflush(stdout);

  TFile fileMcRef(fullMcRefYieldsFileName);
  if (!fileMcRef.IsOpen()) {
    std::cout << "failed to open a file <" << fullMcRefYieldsFileName << ">\n";
    assert(fileMcRef.IsOpen());
  }
  TVectorD *yieldsMcFsrOfRecPtr        = (TVectorD *)fileMcRef.FindObjectAny("yieldsMcPostFsrRecFIArray");
  TVectorD *yieldsMcRecPtr             = (TVectorD *)fileMcRef.FindObjectAny("yieldsMcPostFsrGenFIArray");
  if (!yieldsMcFsrOfRecPtr || !yieldsMcRecPtr) {
    std::cout << "null pointers from <" << fullMcRefYieldsFileName << ">\n";
    assert(0);
  }
  TVectorD yieldsMcFsrOfRecV = *yieldsMcFsrOfRecPtr;
  TVectorD yieldsMcRecV = *yieldsMcRecPtr;

  int nUnfoldingBins = DYTools::getTotalNumberOfBins();
  TVectorD dNdMmcCheckV(nUnfoldingBins);
  dNdMmcCheckV = 0;

  unfolding::unfold(yieldsMcRecV, dNdMmcCheckV, fullUnfoldingConstFileName);

  // Report results
  if (1) {
    printf("\nUNFOLD: Check on the MC, yields:\n");
    for(int i=0; i<DYTools::nMassBins; i++){
      double *rapidityBinLimits=DYTools::getYBinLimits(i);
      for (int yi=0; yi<nYBins[i]; ++yi) {
	int idx=findIndexFlat(i,yi);
	printf("%4.0f-%4.0f  %4.2lf-%4.2lf  %10.0f    %10.0f\n",
	       DYTools::massBinLimits[i],DYTools::massBinLimits[i+1],
	       rapidityBinLimits[yi],rapidityBinLimits[yi+1],
	       yieldsMcFsrOfRecV[idx],
	       dNdMmcCheckV[idx]);
      }
      delete rapidityBinLimits;
    }
  }

  fileMcRef.Close();
  if (yieldsMcFsrOfRecPtr) delete yieldsMcFsrOfRecPtr;
  if (yieldsMcRecPtr) delete yieldsMcRecPtr;
}


void  efficiencyCorrection(TMatrixD &vin, TMatrixD &vinStatErr, TMatrixD &vinSystErr,
             TMatrixD &vout, TMatrixD &voutStatErr, TMatrixD &voutSystErr)
{
  const int nUnfoldingBins= DYTools::getTotalNumberOfBins();

  // Read efficiency constants
  printf("Efficiency: Load constants\n"); fflush(stdout);
    
  TString fileConstantsFullName=TString("../root_files/constants/")+tagDirConstants+TString("/")+fileEfficiencyConstants;
  TFile fileConstants(fileConstantsFullName);
  TMatrixD* efficiencyArrayPtr    = (TMatrixD *)fileConstants.FindObjectAny("efficiencyArray");
  TMatrixD* efficiencyErrArrayPtr = (TMatrixD *)fileConstants.FindObjectAny("efficiencyErrArray");

  if (!efficiencyArrayPtr || !efficiencyErrArrayPtr) {
    std::cout << "at least one needed object is not present in <" << fileConstantsFullName << ">\n";
    assert(0);
  }

  TString fileScaleConstantsFullName=TString("../root_files/constants/")+tagDirConstants+TString("/")+fileScaleFactorConstants;
  TFile fileScaleConstants(fileScaleConstantsFullName);
  TVectorD* rhoDataMcPtr    = (TVectorD *)fileScaleConstants.FindObjectAny("scaleFactorFlatIdxArray");
  TVectorD* rhoDataMcErrPtr = (TVectorD *)fileScaleConstants.FindObjectAny("scaleFactorErrFlatIdxArray");

  if (!rhoDataMcPtr || !rhoDataMcErrPtr) {
    std::cout << "a least one needed object is not present in <" << fileScaleConstantsFullName << ">\n";
    assert(0);
  }

  TMatrixD efficiencyArray= *efficiencyArrayPtr;
  TMatrixD efficiencyErrArray= *efficiencyErrArrayPtr;
  TVectorD rhoDataMc= *rhoDataMcPtr;
  TVectorD rhoDataMcErr= *rhoDataMcErrPtr;

  // Check that the binning is consistent
  bool checkResult = true;
  if (checkResult) checkResult=unfolding::checkRangesMY(efficiencyArray,"efficiencyArray");
  if (checkResult && (rhoDataMc.GetNoElements()!=nUnfoldingBins)) {
    std::cout << "rhoDataMc has incorrect number of entries (" <<
      rhoDataMc.GetNoElements() << " instead of " << nUnfoldingBins << "\n";
    checkResult=false;
  }
  if( !checkResult ){
    printf("Efficiency: ERROR: inconsistent binning in the inputs\n");
    assert(0);
  }else
    printf("Efficiency: Binning in the inputs is consistent\n");

  // Apply the correction
  TMatrixD systErrorPropagated(nMassBins,nMaxYBins);
  TMatrixD systErrorAdditional(nMassBins,nMaxYBins);
  systErrorAdditional = 0;

  for(int i=0; i<nMassBins; i++){
    for (int yi=0; yi<nYBins[i]; ++yi) {
      double effFactor = efficiencyArray[i][yi] * rhoDataMc[i];
      double effErr = effFactor
	* sqrt( efficiencyErrArray[i][yi]*efficiencyErrArray[i][yi]/efficiencyArray[i][yi]/efficiencyArray[i][yi]
		+ rhoDataMcErr[i]*rhoDataMcErr[i]/rhoDataMc[i]/rhoDataMc[i]);
      vout[i][yi]        = vin[i][yi] / effFactor;
      // Statistical error propagated
      voutStatErr[i][yi] = vinStatErr[i][yi] / effFactor;
      // Old systematic error, propagated
      systErrorPropagated[i][yi] = vinSystErr[i][yi] / effFactor;
      // Extra systematic error due to the errors on the efficiency and scale factor
      systErrorAdditional[i][yi] = (vin[i][yi]/effFactor)*(effErr/effFactor);
      systEfficiency[i][yi] = systErrorAdditional[i][yi]/vout[i][yi];
      voutSystErr[i][yi] = 
	sqrt(systErrorPropagated[i][yi]*systErrorPropagated[i][yi]
	     + systErrorAdditional[i][yi]*systErrorAdditional[i][yi]);
    }
  }

  if (printEfficiencyTable) {
  printf("\nEfficiency: Results for the data, yields:\n");
  printf("                after unfolding            eff. factors,%%       rho(data/mc)         efficiency-corrected        syst-err-eff, %%\n");
  for(int i=0; i<nMassBins; i++){
    double *rapidityBinLimits=DYTools::getYBinLimits(i);
    for (int yi=0; yi<nYBins[i]; ++yi) {
      printf("%4.0f-%4.0f  %4.2lf-%4.2lf  %8.1f +- %7.1f +- %7.1f   %5.2f +- %5.2f      %4.3f +- %4.3f      %8.1f +- %7.1f +- %7.1f        %4.2f\n",
	     massBinLimits[i],massBinLimits[i+1],
	     rapidityBinLimits[yi],rapidityBinLimits[yi+1],
	     vin[i][yi], vinStatErr[i][yi], vinSystErr[i][yi],
	     efficiencyArray[i][yi]*100, efficiencyErrArray[i][yi]*100, 
	     rhoDataMc[i], rhoDataMcErr[i],
	     vout[i][yi], voutStatErr[i][yi], voutSystErr[i][yi],
	     systErrorAdditional[i][yi]*100.0/vout[i][yi]);
    }
    delete rapidityBinLimits;
  }
  }

  fileConstants.Close();
  fileScaleConstants.Close();

  return;

}

void  acceptanceCorrection(TMatrixD &vin, TMatrixD &vinStatErr, TMatrixD &vinSystErr,
			   TMatrixD &vout, TMatrixD &voutStatErr, TMatrixD &voutSystErr){
  // Read efficiency constants
  printf("Acceptance: Load constants\n"); fflush(stdout);
    
  TString fileConstantsFullName=TString("../root_files/constants/")+tagDirConstants+TString("/")+fileAcceptanceConstants;
  TFile fileConstants(fileConstantsFullName);
  TMatrixD *acceptanceMatrixPtr    = (TMatrixD *)fileConstants.FindObjectAny("acceptanceMatrix");
  TMatrixD *acceptanceErrMatrixPtr = (TMatrixD *)fileConstants.FindObjectAny("acceptanceErrMatrix");
  if (!acceptanceMatrixPtr || !acceptanceErrMatrixPtr) {
    std::cout << "at least one object from file <" << fileConstantsFullName << "> is null\n";
    assert(0);
  }

  TMatrixD acceptanceMatrix = *acceptanceMatrixPtr;
  TMatrixD acceptanceErrMatrix = *acceptanceErrMatrixPtr;

  TString fileSystematicsFullName=TString("../root_files/systematics/")+tagDirConstants+TString("/")+fileAcceptanceSystematics;
  TFile fileSystematics(fileSystematicsFullName);
			
  TVectorD *acceptanceTheoryErrArrayPtr = (TVectorD *)fileSystematics.FindObjectAny("acceptanceTheoryErrArray");
  if (!acceptanceTheoryErrArrayPtr) {
    std::cout << "failed to get object from <" << fileSystematicsFullName << ">\n";
    assert(0);
  }
  TVectorD acceptanceTheoryErrArray = *acceptanceTheoryErrArrayPtr;

  TString fileAccFSRSystFullName=TString("../root_files/systematics/")+tagDirConstants+TString("/")+fileAcceptanceFSRSystematics;
  TFile fileAccFSRSyst(fileAccFSRSystFullName);
  TMatrixD *acceptanceFSRErrMatrixPtr = (TMatrixD *)fileAccFSRSyst.FindObjectAny("accSystPercent");
  if (!acceptanceFSRErrMatrixPtr) {
    std::cout << "failed to get object from <" << fileAccFSRSystFullName << ">\n";
    assert(0);
  }
  TMatrixD acceptanceFSRErrMatrix=*acceptanceFSRErrMatrixPtr;

  // Check that the binning is consistent
  bool checkResult = true;
  if (checkResult) checkResult=unfolding::checkRangesMY(acceptanceMatrix,"acceptanceMatrix");
  //if( acceptanceMatrix.GetNoElements() != nMassBins) checkResult = false;
  if( acceptanceTheoryErrArray.GetNoElements() != nMassBins) checkResult = false;
  if( !checkResult ){
    printf("Acceptance: ERROR: inconsistent binning in the inputs\n");
    assert(0);
  }else
    printf("Acceptance: Binning in the inputs is consistent\n");

  // Apply the correction
  TMatrixD systErrorPropagated(nMassBins,nMaxYBins);
  TMatrixD systErrorAdditional(nMassBins,nMaxYBins);
  TMatrixD systErrorTheory(nMassBins,nMaxYBins);
  systErrorAdditional = 0;

  for(int i=0; i<nMassBins; i++){
    for (int yi=0; yi<nYBins[i]; ++yi) {
      double accFactor = acceptanceMatrix[i][yi];
      double accErr    = acceptanceErrMatrix[i][yi];
      double accThErr  = acceptanceTheoryErrArray[i];
      double accFSRErr = acceptanceFSRErrMatrix[i][yi]/100;
      systAccTheoryRelative[i][yi]=acceptanceTheoryErrArray[i];
      vout[i][yi]        = vin[i][yi] / accFactor;
      voutStatErr[i][yi] = vinStatErr[i][yi] / accFactor;
      systErrorPropagated[i][yi] = vinSystErr[i][yi]/accFactor;
      systErrorAdditional[i][yi] = (vin[i][yi]/accFactor)*(accErr/accFactor);
      //voutSystErr[i] = sqrt(systErrorPropagated[i]*systErrorPropagated[i] + systErrorAdditional[i]*systErrorAdditional[i]);
      voutSystErr[i][yi] = sqrt(SQR(systErrorPropagated[i][yi]) + SQR(systErrorAdditional[i][yi])
				+ SQR(vout[i][yi]*accThErr) + SQR(vout[i][yi]*accFSRErr));
      systAcceptanceRelative[i][yi]=sqrt(SQR(systErrorAdditional[i][yi])/SQR(vout[i][yi]) 
					 + SQR(accThErr) + SQR(accFSRErr));
    }
  }

  if (printAcceptanceTable) {
    printf("\nAcceptance: Results for the data, yields:\n");
    printf("                          eff-corrected                 acc. factors,%%         acceptance-corrected        syst-err-acc%%,   Theory_Err%% \n");
    for(int i=0; i<nMassBins; i++){
      double *rapidityBinLimits=DYTools::getYBinLimits(i);
      for (int yi=0; yi<nYBins[i]; ++yi) {
	printf("%4.0f-%4.0f  %4.2lf-%4.2lf  %8.1f +- %7.1f +- %7.1f   %5.2f +- %5.2f      %9.1f +- %8.1f +- %8.1f        %4.2f     %4.2f\n",
	       massBinLimits[i],massBinLimits[i+1],
	       rapidityBinLimits[yi],rapidityBinLimits[yi+1],
	       vin[i][yi], vinStatErr[i][yi], vinSystErr[i][yi],
	       acceptanceMatrix[i][yi]*100, acceptanceErrMatrix[i][yi]*100, 
	       vout[i][yi], voutStatErr[i][yi], voutSystErr[i][yi],
	       systErrorAdditional[i][yi]*100/vout[i][yi],100*acceptanceTheoryErrArray[i]);
      }
    }
  }

  fileConstants.Close();
  fileSystematics.Close();

  return;
}

void  fsrCorrectionBase(TMatrixD &vin, TMatrixD &vinStatErr, TMatrixD &vinSystErr,
			TMatrixD &vout, TMatrixD &voutStatErr, TMatrixD &voutSystErr,
			const TString &fname, const TString &title)
{
  // Read efficiency constants
  std::cout << "FsrCorrectionBase(" << title << "): Load constants" << std::endl;
    
  TFile fileConstants(TString("../root_files/constants/")+tagDirConstants+TString("/")+fname);
  assert(fileConstants.IsOpen());
  TMatrixD *fsrCorrectionMatrixPtr  = (TMatrixD *)fileConstants.FindObjectAny("fsrCorrectionMatrix");
  TMatrixD *fsrCorrectionErrMatrixPtr = (TMatrixD *)fileConstants.FindObjectAny("fsrCorrectionErrMatrix");
  assert(fsrCorrectionMatrixPtr && fsrCorrectionErrMatrixPtr);
  TMatrixD fsrCorrectionMatrix= *fsrCorrectionMatrixPtr;
  TMatrixD fsrCorrectionErrMatrix= *fsrCorrectionErrMatrixPtr;

  // Check that the binning is consistent
  bool checkResult = true;
  if (checkResult) checkResult=unfolding::checkRangesMY(fsrCorrectionMatrix,"fsrCorrectionMatrix");
  //if( fsrCorrectionArray.GetNoElements() != nMassBins) checkResult = false;
  if( !checkResult ){
    printf("FsrCorrectionBase: ERROR: inconsistent binning in the inputs\n");
    assert(0);
  }else
    printf("FsrCorrectionBase: Binning in the inputs is consistent\n");

  // Apply the correction

  TMatrixD systErrorPropagated(nMassBins,nMaxYBins);
  for(int i=0; i<nMassBins; i++){
    for (int yi=0; yi<nYBins[i]; ++yi) {
      double fsrFactor = fsrCorrectionMatrix[i][yi];
      double fsrErr    = fsrCorrectionErrMatrix[i][yi];
      vout[i][yi]        = vin[i][yi] / fsrFactor;
      voutStatErr[i][yi] = vinStatErr[i][yi] / fsrFactor;
      systErrorPropagated[i][yi] = vinSystErr[i][yi] / fsrFactor;
      double systNew = (vin[i][yi]/fsrFactor)*(fsrErr/fsrFactor);
      voutSystErr[i][yi] = sqrt(SQR(systErrorPropagated[i][yi]) + SQR(systNew));
    }
  }

  if (printFSRcorrectionTable) {
  printf("\n"); // FsrCorrection
  std::cout << title <<" : Results for the data, yields:\n";
  printf("                          acc-corrected               fsr. factors            fsr-corrected                  sys-err-fsr\n");
  for(int i=0; i<nMassBins; i++){
    double *rapidityBinLimits=DYTools::getYBinLimits(i);
    for (int yi=0; yi<nYBins[i]; ++yi) {
      printf("%4.0f-%4.0f  %4.2f-%4.2f  %9.1f +- %8.1f +- %8.1f   %4.3f +- %4.3f      %9.1f +- %8.1f +- %8.1f        %5.2f\n",
	     massBinLimits[i],massBinLimits[i+1],
	     rapidityBinLimits[yi],rapidityBinLimits[yi+1],
	     vin[i][yi], vinStatErr[i][yi], vinSystErr[i][yi],
	     fsrCorrectionMatrix[i][yi], fsrCorrectionErrMatrix[i][yi], 
	     vout[i][yi], voutStatErr[i][yi], voutSystErr[i][yi],
	     fsrCorrectionErrMatrix[i][yi]*100.0/fsrCorrectionMatrix[i][yi]);
    }
    delete rapidityBinLimits;
  }
  }

  fileConstants.Close();
  return;
}

void  fsrCorrection(TMatrixD &vin, TMatrixD &vinStatErr, TMatrixD &vinSystErr,
		    TMatrixD &vout, TMatrixD &voutStatErr, TMatrixD &voutSystErr)
{
  if (1) {
    fsrCorrectionBase(vin,vinStatErr,vinSystErr, vout,voutStatErr,voutSystErr,
		      fileFsrCorrectionConstants,"FsrCorrection");
  }
  else {
  // Read efficiency constants
  printf("FsrCorrection: Load constants\n"); fflush(stdout);
    
  TFile fileConstants(TString("../root_files/constants/")+tagDirConstants+TString("/")+fileFsrCorrectionConstants);
  assert(fileConstants.IsOpen());
  TMatrixD *fsrCorrectionMatrixPtr  = (TMatrixD *)fileConstants.FindObjectAny("fsrCorrectionMatrix");
  TMatrixD *fsrCorrectionErrMatrixPtr = (TMatrixD *)fileConstants.FindObjectAny("fsrCorrectionErrMatrix");
  assert(fsrCorrectionMatrixPtr && fsrCorrectionErrMatrixPtr);
  TMatrixD fsrCorrectionMatrix= *fsrCorrectionMatrixPtr;
  TMatrixD fsrCorrectionErrMatrix= *fsrCorrectionErrMatrixPtr;

  // Check that the binning is consistent
  bool checkResult = true;
  if (checkResult) checkResult=unfolding::checkRangesMY(fsrCorrectionMatrix,"fsrCorrectionMatrix");
  //if( fsrCorrectionArray.GetNoElements() != nMassBins) checkResult = false;
  if( !checkResult ){
    printf("FsrCorrection: ERROR: inconsistent binning in the inputs\n");
    assert(0);
  }else
    printf("FsrCorrection: Binning in the inputs is consistent\n");

  // Apply the correction

  TMatrixD systErrorPropagated(nMassBins,nMaxYBins);
  for(int i=0; i<nMassBins; i++){
    for (int yi=0; yi<nYBins[i]; ++yi) {
      double fsrFactor = fsrCorrectionMatrix[i][yi];
      double fsrErr    = fsrCorrectionErrMatrix[i][yi];
      vout[i][yi]        = vin[i][yi] / fsrFactor;
      voutStatErr[i][yi] = vinStatErr[i][yi] / fsrFactor;
      systErrorPropagated[i][yi] = vinSystErr[i][yi] / fsrFactor;
      double systNew = (vin[i][yi]/fsrFactor)*(fsrErr/fsrFactor);
      voutSystErr[i][yi] = sqrt(SQR(systErrorPropagated[i][yi]) + SQR(systNew));
    }
  }

  printf("\nFsrCorrection: Results for the data, yields:\n");
  printf("                          acc-corrected               fsr. factors            fsr-corrected                  sys-err-fsr\n");
  for(int i=0; i<nMassBins; i++){
    double *rapidityBinLimits=DYTools::getYBinLimits(i);
    for (int yi=0; yi<nYBins[i]; ++yi) {
      printf("%4.0f-%4.0f  %4.2f-%4.2f  %9.1f +- %8.1f +- %8.1f   %4.3f +- %4.3f      %9.1f +- %8.1f +- %8.1f        %5.2f\n",
	     massBinLimits[i],massBinLimits[i+1],
	     rapidityBinLimits[yi],rapidityBinLimits[yi+1],
	     vin[i][yi], vinStatErr[i][yi], vinSystErr[i][yi],
	     fsrCorrectionMatrix[i][yi], fsrCorrectionErrMatrix[i][yi], 
	     vout[i][yi], voutStatErr[i][yi], voutSystErr[i][yi],
	     fsrCorrectionErrMatrix[i][yi]*100.0/fsrCorrectionMatrix[i][yi]);
    }
    delete rapidityBinLimits;
  }

  fileConstants.Close();
  }
  return;
}

void  fsrCorrectionSansAcceptance(TMatrixD &vin, TMatrixD &vinStatErr, TMatrixD &vinSystErr,
				  TMatrixD &vout, TMatrixD &voutStatErr, TMatrixD &voutSystErr)
{
  if (1) {
    fsrCorrectionBase(vin,vinStatErr,vinSystErr, vout,voutStatErr,voutSystErr,
		      fileFsrCorrectionSansAccConstants, "FsrCorrectionSansAcc");
  }
  else {
  // Read efficiency constants
  printf("FsrCorrection: Load constants\n"); fflush(stdout);
  
  TFile fileConstants(TString("../root_files/constants/")+tagDirConstants+TString("/")+fileFsrCorrectionSansAccConstants);
  TMatrixD fsrCorrectionArray    = *(TMatrixD *)fileConstants.FindObjectAny("fsrCorrectionArray");
  TMatrixD fsrCorrectionErrArray = *(TMatrixD *)fileConstants.FindObjectAny("fsrCorrectionErrArray");
  
  // Check that the binning is consistent
  bool checkResult = true;
  //if( fsrCorrectionArray.GetNoElements() != nMassBins) checkResult = false;
  if (checkResult) checkResult=unfolding::checkRangesMY(fsrCorrectionArray,"fsrCorrectionArray");
  if( !checkResult ){
    printf("FsrCorrectionSansAcceptance: ERROR: inconsistent binning in the inputs\n");
    assert(0);
  }else
    printf("FsrCorrectionSansAcceptance: Binning in the inputs is consistent\n");

  std::cout << "Error this branch is not ready\n";
  assert(0);
  /*
  // Apply the correction
  TMatrixD systErrorPropagated(nMassBins,nMaxYBins);
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
  printf("                acc-corrected                 fsr. factors             fsr-corrected                  sys-err-fsr\n");
  for(int i=0; i<nMassBins; i++){
    printf("%4.0f-%4.0f   %9.1f +- %8.1f +- %8.1f   %4.3f +- %4.3f      %9.1f +- %8.1f +- %8.1f        %5.2f\n",
	   massBinLimits[i],massBinLimits[i+1],
	   vin[i], vinStatErr[i], vinSystErr[i],
	   fsrCorrectionArray[i], fsrCorrectionErrArray[i], 
	   vout[i], voutStatErr[i], voutSystErr[i],
	   fsrCorrectionErrArray[i]*100.0/fsrCorrectionArray[i]);
  }

  fileConstants.Close();
  */
  }
  return;
}

void  crossSections(TMatrixD &vin, TMatrixD &vinStatErr, TMatrixD &vinSystErr,
		    TMatrixD &vout, TMatrixD &voutStatErr, TMatrixD &voutSystErr,
		    TMatrixD &voutNorm, TMatrixD &voutNormStatErr, TMatrixD &voutNormSystErr,
		    const TriggerSelection &triggers)
{

  // Find absolute cross-section
  for(int i=0; i<nMassBins; i++) {
    for (int yi=0; yi<nYBins[i]; ++yi) {
      vout[i][yi] = vin[i][yi] / lumi;
      voutStatErr[i][yi] = vinStatErr[i][yi] / lumi;
      voutSystErr[i][yi] = vinSystErr[i][yi] / lumi;
    }
  }

  // Find normalized cross-section
  int low, high;
  getNormBinRange(low, high);
  double xsecReference = 0;
  double xsecReferenceStatErr = 0;
  double xsecReferenceSystErr = 0;

  TVectorD xsecRefV(nMaxYBins),xsecRefStatErrV(nMaxYBins), xsecRefSystErrV(nMaxYBins);
  xsecRefV=0; xsecRefStatErrV=0; xsecRefSystErrV=0;
  int lowYBins=nYBins[low];

  for( int i=low; i<=high; i++){
    if (lowYBins!=nYBins[i]) {
      std::cout << "y binning error in crossSections. Correct the code: wrong assumption\n";
      assert(0);
    }
    for (int yi=0; yi<nYBins[i]; ++yi) {
      xsecReference += vout[i][yi];
      xsecReferenceStatErr += voutStatErr[i][yi] * voutStatErr[i][yi];
      xsecReferenceSystErr += voutSystErr[i][yi] * voutSystErr[i][yi];
      xsecRefV[yi] += vout[i][yi];
      xsecRefStatErrV[yi] += voutStatErr[i][yi] * voutStatErr[i][yi];
      xsecRefSystErrV[yi] += voutSystErr[i][yi] * voutSystErr[i][yi];
    }
  }

  xsecReferenceStatErr = sqrt(xsecReferenceStatErr);
  xsecReferenceSystErr = sqrt(xsecReferenceSystErr);
  for (int yi=0; yi<nMaxYBins; ++yi) {
    xsecRefStatErrV[yi]= sqrt(xsecRefStatErrV[yi]);
    xsecRefSystErrV[yi]= sqrt(xsecRefSystErrV[yi]);
  }

  // Find normalized cross-section
  for(int i=0; i<nMassBins; i++) {
    for (int yi=0; yi<nYBins[i]; ++yi) {
      voutNorm[i][yi] = vout[i][yi] / xsecReference;
      voutNormStatErr[i][yi] = voutNorm[i][yi]
	* sqrt( SQR(xsecReferenceStatErr) / SQR(xsecReference)
		+ SQR(voutStatErr[i][yi]) / SQR(vout[i][yi]) );
      voutNormSystErr[i][yi] = voutNorm[i][yi]
	* sqrt( SQR(xsecReferenceSystErr) / SQR(xsecReference)
		+ SQR(voutSystErr[i][yi]) / SQR(vout[i][yi]) );
    }
  }

  TMatrixD normXSec  (DYTools::nMassBins,nMaxYBins);
  TMatrixD normXSecErr   (DYTools::nMassBins,nMaxYBins);
  TMatrixD normXSecErrStat   (DYTools::nMassBins,nMaxYBins);
  TMatrixD normXSecErrSyst   (DYTools::nMassBins,nMaxYBins);
  TMatrixD normXSecByBin  (DYTools::nMassBins,nMaxYBins);
  TMatrixD normXSecErrByBin   (DYTools::nMassBins,nMaxYBins);
  TMatrixD normXSecErrByBinStat   (DYTools::nMassBins,nMaxYBins);
  TMatrixD normXSecErrByBinSyst   (DYTools::nMassBins,nMaxYBins);

  if (printPreFSRCrossSectionTable) {
    printf("\nPre FSR cross sections: :\n");
    printf("                    absolute                       normalized +- stat +- sys (total)           (1/sigma)(1/dM)norm +-stat +-syst (total) \n");
    for(int i=0; i<nMassBins; i++){
      double *rapidityBinLimits=DYTools::getYBinLimits(i);
      for (int yi=0; yi<nYBins[i]; ++yi) {
	double binw = massBinLimits[i+1] - massBinLimits[i];
	normXSec[i][yi]=voutNorm[i][yi];
	normXSecErrStat[i][yi]=voutNormStatErr[i][yi];
	normXSecErrSyst[i][yi]=voutNormSystErr[i][yi];
	normXSecErr[i][yi]=sqrt( SQR(voutNormStatErr[i][yi]) + SQR(voutNormSystErr[i][yi]) );

	normXSecByBin[i][yi]=voutNorm[i][yi]/binw;
	normXSecErrByBinStat[i][yi]=voutNormStatErr[i][yi]/binw;
	normXSecErrByBinSyst[i][yi]=voutNormSystErr[i][yi]/binw;
	normXSecErrByBin[i][yi]=sqrt( SQR(voutNormStatErr[i][yi]) + SQR(voutNormSystErr[i][yi]) )/binw;

	printf("%4.0f-%4.0f  %4.2f-%4.2f    %6.1f +- %4.1f +- %4.1f      %1.6f +- %1.6f +- %1.6f  ( %1.6f )     %1.8f +- %1.8f +- %1.8f  ( %1.8f )    \n",
	       massBinLimits[i],massBinLimits[i+1],
	       rapidityBinLimits[yi],rapidityBinLimits[yi+1],
	       vout[i][yi], voutStatErr[i][yi], voutSystErr[i][yi],
	       voutNorm[i][yi], voutNormStatErr[i][yi], voutNormSystErr[i][yi],
	       normXSecErr[i][yi],
	       voutNorm[i][yi]/binw, voutNormStatErr[i][yi]/binw, voutNormSystErr[i][yi]/binw,
	       normXSecErrByBin[i][yi]
	       );
      }
      delete rapidityBinLimits;
    }
  }

  TString outputDir(TString("../root_files/"));
  gSystem->mkdir(outputDir,kTRUE);
  TString xSecResultFileName(outputDir+TString("/xSec_results_") + analysisTag + TString("_") +
			     triggers.triggerConditionsName() + TString(".root"));
  std::cout << "xSecResultFileName= " << xSecResultFileName << "\n";
  TFile fa(xSecResultFileName,"recreate");
  normXSec.Write("normXSec");
  normXSecErr.Write("normXSecErr");
  normXSecErrStat.Write("normXSecErrStat");
  normXSecErrSyst.Write("normXSecErrSyst");
  normXSecByBin.Write("normXSecByBin");
  normXSecErrByBin.Write("normXSecErrByBin");
  normXSecErrByBinStat.Write("normXSecErrByBinStat");
  normXSecErrByBinSyst.Write("normXSecErrByBinSyst");
  unfolding::writeBinningArrays(fa);
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
  TString xSecResultFileName1(outputDir1+TString("/xSecTh_results_") + analysisTag + TString("_") +
			      triggers.triggerConditionsName() + TString(".root"));
  std::cout << "xSecResultFileName1=" << xSecResultFileName1 << "\n";
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

void  crossSectionsDET(TMatrixD &vin, TMatrixD &vinStatErr, TMatrixD &vinSystErr,
		       TMatrixD &vout, TMatrixD &voutStatErr, TMatrixD &voutSystErr,
		       TMatrixD &voutNorm, TMatrixD &voutNormStatErr, TMatrixD &voutNormSystErr)
{

  // Find absolute cross-section
  for(int i=0; i<nMassBins; i++) {
    for (int yi=0; yi<nYBins[i]; ++yi) {
      vout[i][yi] = vin[i][yi] / lumi;
      voutStatErr[i][yi] = vinStatErr[i][yi] / lumi;
      voutSystErr[i][yi] = vinSystErr[i][yi] / lumi;
    }
  }

  // Find normalized cross-section
  int low, high;
  getNormBinRange(low, high);
  double xsecReference = 0;
  double xsecReferenceStatErr = 0;
  double xsecReferenceSystErr = 0;

  for( int i=low; i<=high; i++){
    for (int yi=0; yi<nYBins[i]; ++yi) {
      xsecReference += vout[i][yi];
      xsecReferenceStatErr += voutStatErr[i][yi] * voutStatErr[i][yi];
      xsecReferenceSystErr += voutSystErr[i][yi] * voutSystErr[i][yi];
    }
  }

  xsecReferenceStatErr = sqrt(xsecReferenceStatErr);
  xsecReferenceSystErr = sqrt(xsecReferenceSystErr);

  // Find normalized cross-section
  for(int i=0; i<nMassBins; i++) {
    for (int yi=0; yi<nYBins[i]; ++yi) {
      voutNorm[i][yi] = vout[i][yi] / xsecReference;
      voutNormStatErr[i][yi] = voutNorm[i][yi]
	* sqrt( SQR(xsecReferenceStatErr) / SQR(xsecReference)
		+ SQR(voutStatErr[i][yi]) / SQR(vout[i][yi]) );
      voutNormSystErr[i][yi] = voutNorm[i][yi] 
	* sqrt( SQR(xsecReferenceSystErr) / SQR(xsecReference)
		+ SQR(voutSystErr[i][yi]) / SQR(vout[i][yi]) );
    }
  }

  if (printPreFSR_DET_ShapeTable) {
    printf("\nPre FSR DET shape: :\n");
    printf("                    absolute                      normalized +-stat +-sys (total)\n");
    for(int i=0; i<nMassBins; i++){
      double *rapidityBinLimits=DYTools::getYBinLimits(i);
      for (int yi=0; yi<nYBins[i]; ++yi) {
	printf("%4.0f-%4.0f  %4.2f-%4.2f    %6.1f +- %4.1f +- %4.1f      %1.6f +- %1.6f +- %1.6f  ( %1.6f )\n",
	       //    printf("%4.0f-%4.0f   %9.1f +- %8.1f +- %8.1f   %1.6f +- %1.6f +- %1.6f  ( %1.6f )\n",
	       massBinLimits[i],massBinLimits[i+1],
	       rapidityBinLimits[yi], rapidityBinLimits[yi+1],
	       vout[i][yi], voutStatErr[i][yi], voutSystErr[i][yi],
	       voutNorm[i][yi], voutNormStatErr[i][yi], voutNormSystErr[i][yi],
	       sqrt(SQR(voutNormStatErr[i][yi]) + SQR(voutNormSystErr[i][yi])) );
      }
      delete rapidityBinLimits;
    }
  }
  printf("\nPreFsr cross-section in the Z peak from %3.0f to %3.0f:\n",
	 massBinLimits[low], massBinLimits[high+1]);
  printf("           %9.1f +- %8.1f +- %6.1f \n",
	 xsecReference, xsecReferenceStatErr, xsecReferenceSystErr);
  return;
}


void  postFsrCrossSections(TMatrixD &vin, TMatrixD &vinStatErr, TMatrixD &vinSystErr,
		    TMatrixD &vout, TMatrixD &voutStatErr, TMatrixD &voutSystErr,
		    TMatrixD &voutNorm, TMatrixD &voutNormStatErr, TMatrixD &voutNormSystErr)
{

  // Find absolute cross-section
  for(int i=0; i<nMassBins; i++) {
    for (int yi=0; yi<nYBins[i]; ++yi) {
      vout[i][yi] = vin[i][yi] / lumi;
      voutStatErr[i][yi] = vinStatErr[i][yi] / lumi;
      voutSystErr[i][yi] = vinSystErr[i][yi] / lumi;
    }
  }

  // Find normalized cross-section
  int low, high;
  getNormBinRange(low, high);
  double xsecReference = 0;
  double xsecReferenceStatErr = 0;
  double xsecReferenceSystErr = 0;
  for( int i=low; i<=high; i++){
    for (int yi=0; yi<nYBins[i]; ++yi) {
      xsecReference += vout[i][yi];
      xsecReferenceStatErr += voutStatErr[i][yi] * voutStatErr[i][yi];
      xsecReferenceSystErr += voutSystErr[i][yi] * voutSystErr[i][yi];
    }
  }

  xsecReferenceStatErr = sqrt(xsecReferenceStatErr);
  xsecReferenceSystErr = sqrt(xsecReferenceSystErr);

  // Find normalized cross-section
  for(int i=0; i<nMassBins; i++) {
    for (int yi=0; yi<nYBins[i]; ++yi) {
      voutNorm[i][yi] = vout[i][yi] / xsecReference;
      voutNormStatErr[i][yi] = voutNorm[i][yi]
	* sqrt( SQR(xsecReferenceStatErr) / SQR(xsecReference)
		+ SQR(voutStatErr[i][yi]) / SQR(vout[i][yi]) );
      voutNormSystErr[i][yi] = voutNorm[i][yi] 
	* sqrt( SQR(xsecReferenceSystErr) / SQR(xsecReference)
		+ SQR(voutSystErr[i][yi]) / SQR(vout[i][yi]) );
    }
  }

  if (printPostFSRCrossSectionTable) {
    printf("\nPost FSR cross section: :\n");
    printf("                    absolute                      normalized +-stat +-sys (total)\n");
    for(int i=0; i<nMassBins; i++){
      double *rapidityBinLimits=DYTools::getYBinLimits(i);
      for (int yi=0; yi<nYBins[i]; ++yi) {
	printf("%4.0f-%4.0f  %4.2f-%4.2f    %6.1f +- %4.1f +- %4.1f      %1.6f +- %1.6f +- %1.6f  ( %1.6f )\n",
	       //    printf("%4.0f-%4.0f   %9.1f +- %8.1f +- %8.1f   %1.6f +- %1.6f +- %1.6f  ( %1.6f )\n",
	       massBinLimits[i],massBinLimits[i+1],
	       rapidityBinLimits[yi], rapidityBinLimits[yi+1],
	       vout[i][yi], voutStatErr[i][yi], voutSystErr[i][yi],
	       voutNorm[i][yi], voutNormStatErr[i][yi], voutNormSystErr[i][yi],
	       sqrt(SQR(voutNormStatErr[i][yi]) + SQR(voutNormSystErr[i][yi])) );
      }
      delete rapidityBinLimits;
    }
  }
  printf("\nPostFsr cross-section in the Z peak from %3.0f to %3.0f:\n",
	 massBinLimits[low], massBinLimits[high+1]);
  printf("           %9.1f +- %8.1f +- %6.1f \n",
	 xsecReference, xsecReferenceStatErr, xsecReferenceSystErr);
//     printf("check %f %f\n", xsecReferenceStatErr, xsecReferenceSystErr);

  return;
}

void  postFsrCrossSectionsDET(TMatrixD &vin, TMatrixD &vinStatErr, TMatrixD &vinSystErr,
			      TMatrixD &vout, TMatrixD &voutStatErr, TMatrixD &voutSystErr,
			      TMatrixD &voutNorm, TMatrixD &voutNormStatErr, TMatrixD &voutNormSystErr)
{
  
  // Find absolute cross-section
  for(int i=0; i<nMassBins; i++) {
    for (int yi=0; yi<nYBins[i]; ++yi) {
      vout[i][yi] = vin[i][yi] / lumi;
      voutStatErr[i][yi] = vinStatErr[i][yi] / lumi;
      voutSystErr[i][yi] = vinSystErr[i][yi] / lumi;
    }
  }
  
  // Find normalized cross-section
  int low, high;
  getNormBinRange(low, high);
  double xsecReference = 0;
  double xsecReferenceStatErr = 0;
  double xsecReferenceSystErr = 0;
  for( int i=low; i<=high; i++){
    for ( int yi=0; yi<nYBins[i]; ++yi) {
      xsecReference += vout[i][yi];
      xsecReferenceStatErr += voutStatErr[i][yi] * voutStatErr[i][yi];
      xsecReferenceSystErr += voutSystErr[i][yi] * voutSystErr[i][yi];
    }
  }
  xsecReferenceStatErr = sqrt(xsecReferenceStatErr);
  xsecReferenceSystErr = sqrt(xsecReferenceSystErr);

  // Find normalized cross-section
  for(int i=0; i<nMassBins; i++) {
    for (int yi=0; yi<nYBins[i]; ++yi) {
      voutNorm[i][yi] = vout[i][yi] / xsecReference;
      voutNormStatErr[i][yi] = voutNorm[i][yi]
	* sqrt( SQR(xsecReferenceStatErr) / SQR(xsecReference)
		+ SQR(voutStatErr[i][yi]) / SQR(vout[i][yi]) );
      voutNormSystErr[i][yi] = voutNorm[i][yi] 
	* sqrt( SQR(xsecReferenceSystErr) / SQR(xsecReference)
		+ SQR(voutSystErr[i][yi]) / SQR(vout[i][yi]) );
    }
  }

  if (printPostFSR_DET_ShapeTable) {
    printf("\nPost FSR DET shape: :\n");
    printf("                    absolute                     normalized +-stat +-sys (total)\n");
    for(int i=0; i<nMassBins; i++){
      double *rapidityBinLimits=DYTools::getYBinLimits(i);
      for (int yi=0; yi<nYBins[i]; ++yi) {
	printf("%4.0f-%4.0f  %4.2lf-%4.2lf    %6.1f +- %4.1f +- %4.1f      %1.6f +- %1.6f +- %1.6f  ( %1.6f )\n",
	       massBinLimits[i],massBinLimits[i+1],
	       rapidityBinLimits[yi], rapidityBinLimits[yi+1],
	       vout[i][yi], voutStatErr[i][yi], voutSystErr[i][yi],
	       voutNorm[i][yi], voutNormStatErr[i][yi], voutNormSystErr[i][yi],
	       sqrt(SQR(voutNormStatErr[i][yi]) + SQR(voutNormSystErr[i][yi])) );
      }
      delete rapidityBinLimits;
    }
  }
  return;
}

void printTableForNotes(const TMatrixD& obs, const TMatrixD& obsErr, 
			const TMatrixD& unf, const TMatrixD& unfErr,
			const TMatrixD& ecor, const TMatrixD& ecorErr,
			const TMatrixD& acor, const TMatrixD& acorErr,
			const TMatrixD& fcor, const TMatrixD& fcorErr)
{

  printf("\n\nLatex table for notes\n");
  printf("               obs-bg                  unfolded                 eff-corrected                acc-corrected              fsr-corrected\n");
  for(int i=0; i<nMassBins; i++){
    for (int yi=0; yi<nYBins[i]; ++yi) {
      printf("$%4.0f-%4.0f$ &", massBinLimits[i],massBinLimits[i+1]);
      printf(" yiRange=%2d &$",yi+1);
      printf(" %8.1f \\pm %6.1f $&$", obs[i][yi] , obsErr[i][yi]);
      printf(" %8.1f \\pm %6.1f $&$", unf[i][yi] , unfErr[i][yi]);
      printf(" %8.1f \\pm %6.1f $&$", ecor[i][yi], ecorErr[i][yi]);
      printf(" %9.1f \\pm %8.1f $&$", acor[i][yi], acorErr[i][yi]);
      printf(" %9.1f \\pm %8.1f $ \\\\", fcor[i][yi], fcorErr[i][yi]);
      printf("\n");
    }
  }

  return;
}


void printAllCorrections(){

  TFile fileConstantsEff(TString("../root_files/constants/")+tagDirConstants+TString("/")+fileEfficiencyConstants);
  TMatrixD *efficiencyArrayPtr    = (TMatrixD *)fileConstantsEff.FindObjectAny("efficiencyArray");
  TMatrixD *efficiencyErrArrayPtr = (TMatrixD *)fileConstantsEff.FindObjectAny("efficiencyErrArray");
  assert(efficiencyArrayPtr); 
  assert(efficiencyErrArrayPtr);

  TFile fileScaleConstants(TString("../root_files/constants/")+tagDirConstants+TString("/")+fileScaleFactorConstants);
  TVectorD *rhoDataMcPtr    = (TVectorD *)fileScaleConstants.FindObjectAny("scaleFactorArray");
  TVectorD *rhoDataMcErrPtr = (TVectorD *)fileScaleConstants.FindObjectAny("scaleFactorErrArray");
  assert(rhoDataMcPtr);
  assert(rhoDataMcErrPtr);

  TFile fileConstantsAcc(TString("../root_files/constants/")+tagDirConstants+TString("/")+fileAcceptanceConstants);
  TMatrixD *acceptanceMatrixPtr    = (TMatrixD *)fileConstantsAcc.FindObjectAny("acceptanceMatrix");
  TMatrixD *acceptanceErrMatrixPtr = (TMatrixD *)fileConstantsAcc.FindObjectAny("acceptanceErrMatrix");
  assert(acceptanceMatrixPtr);
  assert(acceptanceErrMatrixPtr);

  TFile fileConstantsFsr(TString("../root_files/constants/")+tagDirConstants+TString("/")+fileFsrCorrectionConstants);
  TMatrixD *fsrCorrectionArrayPtr    = (TMatrixD *)fileConstantsFsr.FindObjectAny("fsrCorrectionArray");
  TMatrixD *fsrCorrectionErrArrayPtr = (TMatrixD *)fileConstantsFsr.FindObjectAny("fsrCorrectionErrArray");
  assert(fsrCorrectionArrayPtr);
  assert(fsrCorrectionErrArrayPtr);

  TFile fileConstantsFsrSansAcc(TString("../root_files/constants/")+tagDirConstants+TString("/")+fileFsrCorrectionSansAccConstants);
  TMatrixD *fsrCorrectionSansAccArrayPtr    = (TMatrixD *)fileConstantsFsrSansAcc.FindObjectAny("fsrCorrectionArray");
  TMatrixD *fsrCorrectionSansAccErrArrayPtr = (TMatrixD *)fileConstantsFsrSansAcc.FindObjectAny("fsrCorrectionErrArray");
  assert(fsrCorrectionSansAccArrayPtr);
  assert(fsrCorrectionSansAccErrArrayPtr);

  TMatrixD efficiencyArray= *efficiencyArrayPtr;
  TMatrixD efficiencyErrArray= *efficiencyErrArrayPtr;
  TVectorD rhoDataMc= *rhoDataMcPtr;
  TVectorD rhoDataMcErr= *rhoDataMcErrPtr;
  TMatrixD acceptanceMatrix= *acceptanceMatrixPtr;
  TMatrixD acceptanceErrMatrix= *acceptanceErrMatrixPtr;
  TMatrixD fsrCorrectionArray= *fsrCorrectionArrayPtr;
  TMatrixD fsrCorrectionErrArray= *fsrCorrectionErrArrayPtr;
  TMatrixD fsrCorrectionSansAccArray= *fsrCorrectionSansAccArrayPtr;
  TMatrixD fsrCorrectionSansAccErrArray= *fsrCorrectionSansAccErrArrayPtr;
  
  std::string fileName="tbl_allCorrections.tex";
  std::ofstream fout;
  fout.open(fileName.c_str());
  fout << "\n\nLatex table of various corrections for PAS/paper\n";
  fout << "\\begin{table}[tbh]\n";
  fout << "\\begin{center}\\begin{tabular}{ccccccc}\n";
  fout << "\\hline\n";
  fout << " Mass range & rapidity range   &            Acceptance, %%   &  Efficiency, %% &  Acc*Eff, %%    &      FSR corr, %%    &     FSR corr in acc, %%\n";
  fout << "\\\\ \\hline\n";

  char buf[200];
  for(int i=0; i<nMassBins; i++){
    double *rapidityBinLimits=DYTools::getYBinLimits(i);
    for (int yi=0; yi<nYBins[i]; ++yi) {
      double effFactor = efficiencyArray[i][yi] * rhoDataMc[i];
//     double effErr = effFactor
//       * sqrt( efficiencyErrArray[i]*efficiencyErrArray[i]/efficiencyArray[i]/efficiencyArray[i]
// 	      + rhoDataMcErr[i]*rhoDataMcErr[i]/rhoDataMc[i]/rhoDataMc[i]);

      double accFactor = acceptanceMatrix[i][yi];
      double accErr    = acceptanceErrMatrix[i][yi];

      double accEff = accFactor * effFactor;
//     double accEffErr = accEff*sqrt( (effErr/effFactor)*(effErr/effFactor) 
// 				    + (accErr/accFactor)*(accErr/accFactor));

      double fsrFactor = fsrCorrectionArray[i][yi];
      double fsrErr    = fsrCorrectionErrArray[i][yi];

      double fsrFactorSansAcc = fsrCorrectionSansAccArray[i][yi];
      double fsrErrSansAcc    = fsrCorrectionSansAccErrArray[i][yi];

      sprintf(buf,"$ %4.0f-%4.0f $&", massBinLimits[i],massBinLimits[i+1]);
      fout << buf;
      sprintf(buf,"$ %4.2f-%4.2f $&", rapidityBinLimits[yi],rapidityBinLimits[yi+1]);
      fout << buf;
      sprintf(buf,"$ %5.2f \\pm %4.2f $&", accFactor*100, accErr*100);
      fout << buf;
      sprintf(buf,"$ %5.2f $&",effFactor);
      fout << buf;
    // For acceptance * efficiency, only acceptance contribution
    // to the error is printed. This is because we are reflecting
    // here MC-related errors. The error on the efficiency
    // is dominated by scale factor errors, related to low 
    // statistics of tag and probe in the data. Those errors
    // are explicitly stated in systematics.
      sprintf(buf,"$ %5.2f \\pm %4.2f $&", accEff*100, accErr*100);
      fout << buf;
      sprintf(buf,"$ %6.2f \\pm %4.2f $&", fsrFactor*100, fsrErr*100);
      fout << buf;
      sprintf(buf,"$ %6.2f \\pm %4.2f $ \\\\", fsrFactorSansAcc*100, fsrErrSansAcc*100);
      fout << buf << "\n";
    }
    delete rapidityBinLimits;
    fout << "\\hline\n";
  }
  fout << "\\end{tabular}\\end{center}\n\\end{table}\n";
  fout.close();
  std::cout << "file <" << fileName << "> created\n";
  return;
}

void printRelativeSystErrors(){

  std::string fileName="tbl_relativeSystErrors.tex";
  std::ofstream fout;
  fout.open(fileName.c_str());
  fout << "\n\nLatex table of relative systematic errors  in percent for PAS/paper\n";
  fout << "\\begin{table}[tbh]\n";
  fout << "\\begin{center}\\begin{tabular}{ccccccccc}\n";
  fout << "\\hline\n";
  fout << " Mass range & rapidity range   &   Escale  &   Eff.   &   Bkg    &    Unfol &    Acc   &    sum   & AccTheory%%\n";   
  fout << "\\\\ \\hline\n";

  char buf[200];
  for(int i=0; i<nMassBins; i++){
    double *rapidityBinLimits=DYTools::getYBinLimits(i);
    for (int yi=0; yi<nYBins[i]; ++yi) {
      // Factor out theory error from the total acceptance error
      double systAcceptanceExpRelative 
	= sqrt( systAcceptanceRelative[i][yi]*systAcceptanceRelative[i][yi]
		- systAccTheoryRelative[i][yi]*systAccTheoryRelative[i][yi]);
      // The "sum" contains only the experimental systematic errors
      double sum = sqrt(systEscaleRelative[i][yi]*systEscaleRelative[i][yi]
			+ systEfficiency[i][yi]*systEfficiency[i][yi]
			+ systBackgrRelative[i][yi]*systBackgrRelative[i][yi]
			+ systUnfoldRelative[i][yi]*systUnfoldRelative[i][yi]
			+ systAcceptanceExpRelative*systAcceptanceExpRelative);
      sprintf(buf,"$%4.0f-%4.0f $&$ %4.2f-%4.2f $&", massBinLimits[i],massBinLimits[i+1],rapidityBinLimits[yi],rapidityBinLimits[yi+1]);
      fout << buf;
      sprintf(buf,"$  %5.1f  $&$ %5.1f $&$ %5.1f $&$ %5.1f $&$ %5.1f $&$ %5.1f $&$ %5.1f $ \\\\" ,
	      100*systEscaleRelative[i][yi], 
	      100*systEfficiency[i][yi], 
	      100*systBackgrRelative[i][yi], 
	      100*systUnfoldRelative[i][yi],
	      100*systAcceptanceExpRelative,
	      100*sum,
	      100*systAccTheoryRelative[i][yi]);
      fout << buf << "\n";
    }
    delete rapidityBinLimits;
    fout << "\\hline\n";
  }
  fout << "\\end{tabular}\\end{center}\n\\end{table}\n";
  fout.close();
  std::cout << "file <" << fileName << "> created\n";
  return;
}

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

