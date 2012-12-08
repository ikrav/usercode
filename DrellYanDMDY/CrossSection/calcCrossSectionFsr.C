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
#include "../Include/latexPrintouts.hh"
#include "../Include/MitStyleRemix.hh"

using std::string;
using std::stringstream;

template<class T> T SQR(const T &x) { return x*x; }

// This global constants will be filled from 
// the configuration file. This is not great C++...
TString gTagDirYields = "";
TString gTagDirConstants = "";
TString gTagDirScaleFactorConstants = ""; // event scale factors
TString gTagDirXSect = ""; // will be set to tagDirScaleFactorConstants
// but will be placed to ../root_files/${tagDirXSect}/"
Double_t lumi = 0;


typedef enum { _fsrCorr_binByBin=1, _fsrCorr_unfPure, _fsrCorr_unf, _fsrCorr_unfMdf, _fsrCorr_unfGood } TFsrCorrectionType_t;



TString fsrCorrectionName(TFsrCorrectionType_t kind) {
  TString name;
  switch(kind) {
  case _fsrCorr_binByBin: name="_fsrBbB"; break;
  case _fsrCorr_unfPure: name="_fsrUnfPure"; break;
  case _fsrCorr_unf: name="_fsrUnf"; break;
  case _fsrCorr_unfMdf: name="_fsrUnfMdf"; break;
  case _fsrCorr_unfGood: name="_fsrUnfGood"; break;
  default:
    name="_fsrUNKNOWN";
  }
  return name;
}


const int includeEScaleSystematics=1;
const int includeUnfoldingSystematics=1;
const int doClosureTest=1;
const TFsrCorrectionType_t fsrCorrection_BinByBin=_fsrCorr_unfGood;
const int useExactVectorsForMcClosureTest=0;
// Use externally provided errors on norm x-sec, or internally calcylated 13-bin sum from 60 to 120 GeV
const int use1binErrorsForNorm=1;   


const int printFSRcorrectionTable=0;
const int printEfficiencyTable=0;
const int printAcceptanceTable=0;
const int printPreFSRCrossSectionTable=0;
const int printPreFSR_DET_ShapeTable=0;
const int printPostFSRCrossSectionTable=0;
const int printPostFSR_DET_ShapeTable=0;

const int printAllCorrectionTable=0;
const int printRelativeSystErrTable=1;
const int printRelativeSystErrTableForCovarianceMatrix=1;
const int callPrintTableForNotes=0;


TString pathYields;
TString pathConstants;
TString pathSystematics;
TString pathScaleFactors;
TString pathXSect;

TString fnameDataYields;  // bg-subtracted
TString fnameMcReferenceYields; // yields for MC closure test of unfolding
TString fnameMcReferenceYieldsFsr; // yields for MC closure test of FSR unfolding
TString fnameMcReferenceYieldsFsrDET; // yields for MC closure test of FSR unfolding in the DET acceptance

TString fnameUnfoldingConstants; // unfolding matrix for DET response
TString fnameUnfoldingSystErrors; // unfolding (DET response) systematic errors (relative)
TString fnameEscaleSystErrors; // escale systematic errors (relative)

TString fnameEfficiencyConstants; // event_efficiency_constants
TString fnameScaleFactorConstants; // scale_factors
TString fnameAcceptanceConstants; // "acceptance_constants"
TString fnameAcceptanceSystematics; // "theoretical_uncertainties";
TString fnameAcceptanceFSRSystematics; // "acceptance_FSR_systematics"

TString fnameFsrCorrectionConstantsBbB; // FSR correction bin-by-bin
TString fnameFsrCorrectionSansAccConstantsBbB; // FSR correction bin-by-bin in the DET acceptance
TString fnameFsrCorrectionConstantsUnf; // FSR correction by unfolding
TString fnameFsrCorrectionDETConstantsUnf; // FSR correction by unfolding in the DET acceptance
TString fnameFsrDETcorrFactors; // FSR correction factors for DET space


// Forward declarations
//-----------------------------------------------------------------

void initGlobalFileNames(const TriggerSelection &triggers, const TString &tagDirYields, const TString &tagDirConstants, const TString &tagDirScaleFactorConstants, const TString &tagDirXSect, TFsrCorrectionType_t fsrBinByBin);


void readData(TMatrixD &v, TMatrixD &vErr1, TMatrixD &vErr2);
void saveYields(const TMatrixD &cs, const TMatrixD &csErr, const TMatrixD &csErrSyst, const TriggerSelection &triggers, const TString specTag, int DET);
//void printYields(const TString &name, const TMatrixD &cs, const TMatrixD &csErr, const TMatrixD &csErrSyst);

void  applyUnfolding(const TMatrixD &vin, const TMatrixD &vinStatErr, const TMatrixD &vinSystErr,
             TMatrixD &vout, TMatrixD &voutStatErr, TMatrixD &voutSystErr);
void  applyUnfoldingToMc(int fsr); //TString fullUnfoldingConstFileName, TString fullMcRefYieldsFileName);
void  efficiencyCorrection(const TMatrixD &vin, const TMatrixD &vinStatErr, const TMatrixD &vinSystErr,
             TMatrixD &vout, TMatrixD &voutStatErr, TMatrixD &voutSystErr);
void  acceptanceCorrection(const TMatrixD &vin, const TMatrixD &vinStatErr, const TMatrixD &vinSystErr,
             TMatrixD &vout, TMatrixD &voutStatErr, TMatrixD &voutSystErr);
void  fsrCorrectionBase(const TMatrixD &vin, const TMatrixD &vinStatErr, const TMatrixD &vinSystErr,
			TMatrixD &vout, TMatrixD &voutStatErr, TMatrixD &voutSystErr, 
			const TString &fname, const TString &title);
void  fsrCorrection(const TMatrixD &vin, const TMatrixD &vinStatErr, const TMatrixD &vinSystErr,
             TMatrixD &vout, TMatrixD &voutStatErr, TMatrixD &voutSystErr);
void  fsrCorrectionSansAcceptance(const TMatrixD &vin, const TMatrixD &vinStatErr, const TMatrixD &vinSystErr,
				  TMatrixD &vout, TMatrixD &voutStatErr, TMatrixD &voutSystErr);

void  crossSections(const TMatrixD &vin, const TMatrixD &vinStatErr, const TMatrixD &vinSystErr,
		    TMatrixD &vout, TMatrixD &voutStatErr, TMatrixD &voutSystErr,
		    TMatrixD &voutNorm, TMatrixD &voutNormStatErr, TMatrixD &voutNormSystErr,
		    const TriggerSelection &triggers, const TString specTag="");

void  crossSectionsDET(const TMatrixD &vin, const TMatrixD &vinStatErr, const TMatrixD &vinSystErr,
		       TMatrixD &vout, TMatrixD &voutStatErr, TMatrixD &voutSystErr,
		       TMatrixD &voutNorm, TMatrixD &voutNormStatErr, 
		       TMatrixD &voutNormSystErr, 
		       const TriggerSelection &triggers, const TString specTag="");

void  postFsrCrossSections(const TMatrixD &vin, const TMatrixD &vinStatErr, const TMatrixD &vinSystErr,
			   TMatrixD &vout, TMatrixD &voutStatErr, TMatrixD &voutSystErr,
			   TMatrixD &voutNorm, TMatrixD &voutNormStatErr, TMatrixD &voutNormSystErr,
			   const TriggerSelection &triggers, const TString specTag="");

void  postFsrCrossSectionsDET(const TMatrixD &vin, const TMatrixD &vinStatErr, const TMatrixD &vinSystErr,
			      TMatrixD &vout, TMatrixD &voutStatErr, TMatrixD &voutSystErr,
			      TMatrixD &voutNorm, TMatrixD &voutNormStatErr, TMatrixD &voutNormSystErr,
			      const TriggerSelection &triggers, const TString specTag="");

void printTableForNotes(const TMatrixD &obs, const TMatrixD &obsErr, 
			const TMatrixD &unf, const TMatrixD &unfErr,
			const TMatrixD &ecor, const TMatrixD &ecorErr,
			const TMatrixD &acor, const TMatrixD &acorErr,
			const TMatrixD &fcor, const TMatrixD &fcorErr);

void printAllCorrections();
void printRelativeSystErrors();
void printRelativeSystErrorsPAS();
void printRelativeSystErrorsForCovarianceMatrix();

void plotAccEff(); // make a plot: A, MCeff, A*MCeff

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

TMatrixD unfoldedYields(DYTools::nMassBins,nMaxYBins);
TMatrixD unfoldedYieldsStatErr(DYTools::nMassBins,nMaxYBins);
TMatrixD unfoldedYieldsSystErr(DYTools::nMassBins,nMaxYBins);
  
// The constants below are obtained by running the full chain
// in 1-bin mode for 60-120 GeV range. They need to be updated
// whenever there are non-trivial chagnes to any part of the
// full chain. For more notes, see below the locations where
// these constants are used.
double normAbsStatErrPreFsrFullAcc  = 0.9;
double normAbsStatErrPostFsrFullAcc = 0.9;
double normAbsStatErrPreFsrDET      = 0.4;
double normAbsStatErrPostFsrDET     = 0.4;
//
double normAbsSystErrPreFsrFullAcc  = 7.3;
double normAbsSystErrPostFsrFullAcc = 7.0;
double normAbsSystErrPreFsrDET      = 3.5;
double normAbsSystErrPostFsrDET     = 3.4;

// ---------------------------------------------------------------
// Main function
// ---------------------------------------------------------------

void calcCrossSectionFsr(const TString conf) { //="../config_files/xsecCalc.conf")

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
      gTagDirYields = TString(line);
      state++;
    }else if(state==2){
      gTagDirConstants = TString(line);
      state++;
    } 
    else if (state==3) {
      gTagDirScaleFactorConstants = TString(line);
      if (PosOk(line,"hltEff")) {
	std::cout << "input file probably does not contain a line for tagDirESFConstants\n";
	assert(0);
      }
      state++;
    }
    else if (state==4) {
      triggerSetString = TString(line);
      break;
    }
  }
  ifs.close();
  if (!DYTools::checkTotalLumi(lumi)) return;

  gTagDirXSect = gTagDirScaleFactorConstants + fsrCorrectionName(fsrCorrection_BinByBin);
  CPlot::sOutDir = TString("plots_") + DYTools::analysisTag + 
    TString("_") +  gTagDirScaleFactorConstants;

  // Construct the trigger object
  TriggerSelection triggers(triggerSetString, true, 0);
  if (!triggers.isDefined()) {
    std::cout << "failed to construct trigger object from <" << triggerSetString << ">\n";
    assert ( triggers.isDefined() );
  }

  // Prepare global file names
  initGlobalFileNames(triggers,gTagDirYields,gTagDirConstants,gTagDirScaleFactorConstants,gTagDirXSect, fsrCorrection_BinByBin);

  if (DYTools::study2D==0) plotAccEff();

  // Do a closure test on MC
  if (doClosureTest) {
  applyUnfoldingToMc(0);

  if (0) if (fsrCorrection_BinByBin!=_fsrCorr_binByBin) {
    applyUnfoldingToMc(1);
    applyUnfoldingToMc(2);
  }

  if (1 || (fsrCorrection_BinByBin!=_fsrCorr_binByBin)) {
    applyUnfoldingToMc(3);
    applyUnfoldingToMc(4);
  }
  }


  TMatrixD signalYields(DYTools::nMassBins,nMaxYBins);
  TMatrixD signalYieldsStatErr(DYTools::nMassBins,nMaxYBins);
  TMatrixD signalYieldsSystErr(DYTools::nMassBins,nMaxYBins);
  
  //TMatrixD unfoldedYields(DYTools::nMassBins,nMaxYBins);
  //TMatrixD unfoldedYieldsStatErr(DYTools::nMassBins,nMaxYBins);
  //TMatrixD unfoldedYieldsSystErr(DYTools::nMassBins,nMaxYBins);
  
  TMatrixD effCorrectedYields(DYTools::nMassBins,nMaxYBins);
  TMatrixD effCorrectedYieldsStatErr(DYTools::nMassBins,nMaxYBins);
  TMatrixD effCorrectedYieldsSystErr(DYTools::nMassBins,nMaxYBins);
  
  TMatrixD accCorrectedYields(DYTools::nMassBins,nMaxYBins);
  TMatrixD accCorrectedYieldsStatErr(DYTools::nMassBins,nMaxYBins);
  TMatrixD accCorrectedYieldsSystErr(DYTools::nMassBins,nMaxYBins);
  
  TMatrixD preFsrYields(DYTools::nMassBins,nMaxYBins);
  TMatrixD preFsrYieldsStatErr(DYTools::nMassBins,nMaxYBins);
  TMatrixD preFsrYieldsSystErr(DYTools::nMassBins,nMaxYBins);
  
  TMatrixD preFsrSansAccYields(DYTools::nMassBins,nMaxYBins);
  TMatrixD preFsrSansAccYieldsStatErr(DYTools::nMassBins,nMaxYBins);
  TMatrixD preFsrSansAccYieldsSystErr(DYTools::nMassBins,nMaxYBins);
  
  TMatrixD absCrossSection(DYTools::nMassBins,nMaxYBins);
  TMatrixD absCrossSectionStatErr(DYTools::nMassBins,nMaxYBins);
  TMatrixD absCrossSectionSystErr(DYTools::nMassBins,nMaxYBins);
  
  TMatrixD absCrossSectionDET(DYTools::nMassBins,nMaxYBins);
  TMatrixD absCrossSectionStatErrDET(DYTools::nMassBins,nMaxYBins);
  TMatrixD absCrossSectionSystErrDET(DYTools::nMassBins,nMaxYBins);
  
  TMatrixD relCrossSection(DYTools::nMassBins,nMaxYBins);
  TMatrixD relCrossSectionStatErr(DYTools::nMassBins,nMaxYBins);
  TMatrixD relCrossSectionSystErr(DYTools::nMassBins,nMaxYBins);

  TMatrixD relCrossSectionDET(DYTools::nMassBins,nMaxYBins);
  TMatrixD relCrossSectionStatErrDET(DYTools::nMassBins,nMaxYBins);
  TMatrixD relCrossSectionSystErrDET(DYTools::nMassBins,nMaxYBins);

  TMatrixD absPostFsrCrossSection(DYTools::nMassBins,nMaxYBins);
  TMatrixD absPostFsrCrossSectionStatErr(DYTools::nMassBins,nMaxYBins);
  TMatrixD absPostFsrCrossSectionSystErr(DYTools::nMassBins,nMaxYBins);
  
  TMatrixD absPostFsrCrossSectionDET(DYTools::nMassBins,nMaxYBins);
  TMatrixD absPostFsrCrossSectionStatErrDET(DYTools::nMassBins,nMaxYBins);
  TMatrixD absPostFsrCrossSectionSystErrDET(DYTools::nMassBins,nMaxYBins);
  
  TMatrixD relPostFsrCrossSection(DYTools::nMassBins,nMaxYBins);
  TMatrixD relPostFsrCrossSectionStatErr(DYTools::nMassBins,nMaxYBins);
  TMatrixD relPostFsrCrossSectionSystErr(DYTools::nMassBins,nMaxYBins);

  TMatrixD relPostFsrCrossSectionDET(DYTools::nMassBins,nMaxYBins);
  TMatrixD relPostFsrCrossSectionStatErrDET(DYTools::nMassBins,nMaxYBins);
  TMatrixD relPostFsrCrossSectionSystErrDET(DYTools::nMassBins,nMaxYBins);

  systBackgrBeforeUnfolding = 0;
  systEfficiency = 0;
  systOthers = 0;

  // Read data yields from file (background subtraction is already done)
  std::cout << "read bg-subtracted data yields\n";
  readData(signalYields, signalYieldsStatErr, signalYieldsSystErr);

  saveYields(signalYields,signalYieldsStatErr,signalYieldsSystErr, triggers,"debug_raw",0);

  // Apply unfolding
  std::cout << "apply unfolding\n";
  applyUnfolding(signalYields, signalYieldsStatErr, signalYieldsSystErr,
		 unfoldedYields, unfoldedYieldsStatErr, unfoldedYieldsSystErr);
  saveYields(unfoldedYields,unfoldedYieldsStatErr,unfoldedYieldsSystErr, triggers,"debug_unf",0);

  // Apply efficiency correction
  std::cout << "efficiency correction\n";
  efficiencyCorrection(unfoldedYields, unfoldedYieldsStatErr, unfoldedYieldsSystErr,
		       effCorrectedYields, effCorrectedYieldsStatErr, effCorrectedYieldsSystErr);
  saveYields(effCorrectedYields,effCorrectedYieldsStatErr,effCorrectedYieldsSystErr, triggers, "debug_effCorr",0);

  // Acceptance corrections
  std::cout << "acceptance correction\n";
  acceptanceCorrection(effCorrectedYields, effCorrectedYieldsStatErr, effCorrectedYieldsSystErr,
		       accCorrectedYields, accCorrectedYieldsStatErr, accCorrectedYieldsSystErr);
  saveYields(accCorrectedYields,accCorrectedYieldsStatErr,accCorrectedYieldsSystErr, triggers, "debug_effAccCorr",0);
  
  // FSR corrections
  std::cout << "fsr correction\n";
  fsrCorrection(accCorrectedYields, accCorrectedYieldsStatErr, accCorrectedYieldsSystErr,
		preFsrYields, preFsrYieldsStatErr, preFsrYieldsSystErr);
  saveYields(preFsrYields,preFsrYieldsStatErr,preFsrYieldsSystErr,
	     triggers,"debug_preFSR",0);

  // Alternative (for DET shapes): all corrections except for acceptance correction 
  fsrCorrectionSansAcceptance(effCorrectedYields, effCorrectedYieldsStatErr, effCorrectedYieldsSystErr,
			      preFsrSansAccYields, preFsrSansAccYieldsStatErr, preFsrSansAccYieldsSystErr);
  saveYields(preFsrSansAccYields,preFsrSansAccYieldsStatErr,preFsrSansAccYieldsSystErr,
	     triggers,"debug_preFSRDET",0);
 
  // Calculate absolute and relative cross-sections
  std::cout << "absolute and relative cross sections" << std::endl;
  std::cout << "1. crossSections" << std::endl;
  if (DYTools::study2D==0) {
    crossSections(preFsrYields, preFsrYieldsStatErr, preFsrYieldsSystErr,
		  absCrossSection, absCrossSectionStatErr, absCrossSectionSystErr,
		  relCrossSection, relCrossSectionStatErr, relCrossSectionSystErr,
		  triggers,"");
  }
  else {
    std::cout << "...skipped\n";
  }

  // Calculate absolute and relative cross-sections DET (shapes, no Acc, but with FSR)
  std::cout << "2. crossSectionsDET" << std::endl;
  crossSectionsDET(preFsrSansAccYields, preFsrSansAccYieldsStatErr, preFsrSansAccYieldsSystErr,
		   absCrossSectionDET, absCrossSectionStatErrDET, absCrossSectionSystErrDET,
		   relCrossSectionDET, relCrossSectionStatErrDET, relCrossSectionSystErrDET,
		   triggers,"");

  // Also, calculate absolute and relative cross-sections for post-FSR stage
  std::cout << "3. postFsrCrossSection" << std::endl;
  postFsrCrossSections(accCorrectedYields, accCorrectedYieldsStatErr, accCorrectedYieldsSystErr,
		absPostFsrCrossSection, absPostFsrCrossSectionStatErr, absPostFsrCrossSectionSystErr,
	       relPostFsrCrossSection, relPostFsrCrossSectionStatErr, relPostFsrCrossSectionSystErr,
		       triggers, "");

  // calculate absolute and relative cross-sections DET (shapes, no Acc, no FSR corrections)
  std::cout << "4. postFsrCrossSectionsDET" << std::endl;
  postFsrCrossSectionsDET(effCorrectedYields, effCorrectedYieldsStatErr, effCorrectedYieldsSystErr,
			  absPostFsrCrossSectionDET, absPostFsrCrossSectionStatErrDET, absPostFsrCrossSectionSystErrDET,
			  relPostFsrCrossSectionDET, relPostFsrCrossSectionStatErrDET, relPostFsrCrossSectionSystErrDET,
			  triggers, "");

  // Output
  std::cout << "\nprint output\n" << std::endl;
  if (printAllCorrectionTable) printAllCorrections();
  if (printRelativeSystErrTable){
    printRelativeSystErrors();
    printRelativeSystErrorsPAS();
  }
  if (printRelativeSystErrTableForCovarianceMatrix) printRelativeSystErrorsForCovarianceMatrix();

  if (callPrintTableForNotes) {
  std::cout << "\nprintTablesForNotes\n" << std::endl;
  printTableForNotes(signalYields      , signalYieldsStatErr, 
		     unfoldedYields    , unfoldedYieldsStatErr,
		     effCorrectedYields, effCorrectedYieldsStatErr,
		     accCorrectedYields, accCorrectedYieldsStatErr,
		     preFsrYields      , preFsrYieldsStatErr);
 
  }

  latexPrintoutCrossSection(signalYields,       signalYieldsStatErr, 
		            unfoldedYields,     unfoldedYieldsStatErr,
		            effCorrectedYields, effCorrectedYieldsStatErr,
		            accCorrectedYields, accCorrectedYieldsStatErr,
		            preFsrYields      , preFsrYieldsStatErr ,
                            relCrossSection,           relCrossSectionStatErr, 
                                                       relCrossSectionSystErr,
                            relCrossSectionDET,        relCrossSectionStatErrDET, 
                                                       relCrossSectionSystErrDET,
                            relPostFsrCrossSection,    relPostFsrCrossSectionStatErr, 
                                                       relPostFsrCrossSectionSystErr,
                            relPostFsrCrossSectionDET, relPostFsrCrossSectionStatErrDET, 
                                                       relPostFsrCrossSectionSystErrDET,
                            absCrossSection,           absCrossSectionStatErr, 
                                                       absCrossSectionSystErr,
                            absCrossSectionDET,        absCrossSectionStatErrDET, 
                                                       absCrossSectionSystErrDET,
                            absPostFsrCrossSection,    absPostFsrCrossSectionStatErr, 
                                                       absPostFsrCrossSectionSystErr,
                            absPostFsrCrossSectionDET, absPostFsrCrossSectionStatErrDET, 
                                                       absPostFsrCrossSectionSystErrDET,
                            "CrossSection/calcCrossSectionFsr.C");

  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////
  //draw and save plots
 
  PlotMatrixVariousBinning(relCrossSection, "relative_CS", "LEGO2", 0, "Pre FSR All Phase Space", 1);
  PlotMatrixVariousBinning(relCrossSectionDET, "relative_CS_DET", "LEGO2", 0, "Pre FSR Detector Phase space", 1);
  PlotMatrixVariousBinning(relPostFsrCrossSection, "relative_postFSR_CS", "LEGO2", 0, "Post FSR All Phase Space", 1);
  PlotMatrixVariousBinning(relPostFsrCrossSectionDET, "relative_postFSR_CS_DET", "LEGO2", 0, "Post FSR Detector Phase space", 1);


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

  std::cout << "Load data yields from <" << fnameDataYields << ">" << std::endl;
  TFile fileYields   (fnameDataYields);
  TMatrixD *YieldsSignalPtr       = (TMatrixD *)fileYields.FindObjectAny("YieldsSignal");
  TMatrixD *YieldsSignalErrPtr    = (TMatrixD *)fileYields.FindObjectAny("YieldsSignalErr");
  TMatrixD *YieldsSignalSystErrPtr= (TMatrixD *)fileYields.FindObjectAny("YieldsSignalSystErr");
  TVectorD *MassBinLimitsForYieldsPtr = (TVectorD *)fileYields.FindObjectAny("massBinLimits");
  if (!MassBinLimitsForYieldsPtr) MassBinLimitsForYieldsPtr = (TVectorD *)fileYields.FindObjectAny("massBinning");
  TVectorD *YBinCountsForYieldsPtr = (TVectorD *)fileYields.FindObjectAny("rapidityCounts");
  fileYields.Close();

  if (!YieldsSignalPtr || !YieldsSignalErrPtr || !YieldsSignalSystErrPtr ||
      !MassBinLimitsForYieldsPtr || !YBinCountsForYieldsPtr) {
    std::cout << "failed to get at least one of the required objects from file <" << fnameDataYields << ">\n";
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
  if (checkResult) checkResult=unfolding::checkBinningRanges(MassBinLimitsForYields,YBinCountsForYields,fnameDataYields);
  if (checkResult) checkResult=unfolding::checkRangesMY(YieldsSignal,"YieldsSignal");
  if( !checkResult ){
    printf("ERROR: inconsistent binning in the inputs\n");
    assert(0);
  }else
    printf("readData: Binning in the inputs is consistent\n");

  // Prepare output yields and errors
  for(int i=0; i<DYTools::nMassBins; i++){
    for (int yi=0; yi<DYTools::nYBins[i]; yi++) {
      v[i][yi] = YieldsSignal[i][yi];
      vErr1[i][yi] = YieldsSignalErr[i][yi];
      // Background systematics should be already in, add 
      // energy scale systematics
      vErr2[i][yi] = YieldsSignalSystErr[i][yi];
      systBackgrBeforeUnfolding[i][yi] = YieldsSignalSystErr[i][yi];
    }
  }

  //printYields(fnameDataYields, v,vErr1,vErr2);

  return;
}

//-----------------------------------------------------------------

void saveYields(const TMatrixD &vin, const TMatrixD &vinStatErr, const TMatrixD &vinSystErr,
		const TriggerSelection &triggers, const TString specTag, int DET)
{

  TMatrixD zero=vin;
  zero=0;
  TMatrixD vout=zero;
  TMatrixD voutStatErr=zero;
  TMatrixD voutSystErr=zero;
  TMatrixD voutNorm=zero;
  TMatrixD voutNormStatErr=zero;
  TMatrixD voutNormSystErr=zero;

  //printYields(specTag,vin,vinStatErr,vinSystErr);

  if (DET==0) {
    crossSections(vin,vinStatErr,vinSystErr,
		  vout,voutStatErr,voutSystErr,
		  voutNorm,voutNormStatErr,voutNormSystErr,
		  triggers,specTag);
  }
  else {
    crossSectionsDET(vin,vinStatErr,vinSystErr,
		     vout,voutStatErr,voutSystErr,
		     voutNorm,voutNormStatErr,voutNormSystErr,
		     triggers,specTag);
  }
  return;
}

//-----------------------------------------------------------------
// Unfold
//-----------------------------------------------------------------
void  applyUnfolding(const TMatrixD &vinM, const TMatrixD &vinStatErrM, const TMatrixD &vinSystErrM,
             TMatrixD &voutM, TMatrixD &voutStatErrM, TMatrixD &voutSystErrM)
{

  int nUnfoldingBins = DYTools::getTotalNumberOfBins();
  TVectorD vin(nUnfoldingBins),vinStatErr(nUnfoldingBins),vinSystErr(nUnfoldingBins);
  TVectorD vout(nUnfoldingBins),voutStatErr(nUnfoldingBins),voutSystErr(nUnfoldingBins);

  // First, propagate through unfolding the signal yields with stat and syst errors
  assert(unfolding::unfold(vinM, voutM, fnameUnfoldingConstants, vin, vout)==1);
  assert(unfolding::propagateErrorThroughUnfolding(vinStatErrM,voutStatErrM, fnameUnfoldingConstants, vinStatErr, voutStatErr)==1);
  assert(unfolding::propagateErrorThroughUnfolding(vinSystErrM,voutSystErrM, fnameUnfoldingConstants, vinSystErr, voutSystErr)==1);

  // Second, propagate separately systematic error components that need it.
  // These are already included in the total systematic error above in vinSystErr,
  // however we do it separately so that we can quote the breakdown in the
  // table of systematic errors
  TMatrixD systBackgrM(DYTools::nMassBins,nMaxYBins);
  TVectorD systBackgrBeforeUnfoldingV(nUnfoldingBins), systBackgrV(nUnfoldingBins);
  unfolding::propagateErrorThroughUnfolding(systBackgrBeforeUnfolding, systBackgrM, fnameUnfoldingConstants, systBackgrBeforeUnfoldingV,systBackgrV);

  // The electron energy scale systematics that is loaded here
  // is estimated on the unfolded yields. So we read it in at this time
  // and will add to the outgoing total systematic below in this function.
  TVectorD systEscaleV(nUnfoldingBins);
  systEscaleV=0;
  if (includeEScaleSystematics) {
    TFile fileEscaleSystematics(fnameEscaleSystErrors);
    if( ! fileEscaleSystematics.IsOpen()){
      printf("ERROR: required file with escale errors %s is not found!\n",
	     fnameEscaleSystErrors.Data());
      assert(0);
    }
    assert(unfolding::checkBinningArrays(fileEscaleSystematics));
    TVectorD *escaleSystematicsPercentPtr
      = (TVectorD *)fileEscaleSystematics.FindObjectAny("escaleSystPercentFI");
    assert(escaleSystematicsPercentPtr);
    assert(unfolding::checkRangesFI(*escaleSystematicsPercentPtr,"escaleSystPercent"));
    TVectorD escaleSystematicsPercent = *escaleSystematicsPercentPtr;
  
    systEscaleV=0;
    for(int idx=0; idx<nUnfoldingBins; ++idx) {
      systEscaleV[idx] = (escaleSystematicsPercent[idx]/100.0) * vout[idx];
    }
  }

  // Pool together the unfolding systematics and add it to the total systematics
  TVectorD systUnfoldingV(nUnfoldingBins);
  systUnfoldingV=0;
  if (includeUnfoldingSystematics) {
    unfolding::calculateTotalUnfoldingSystErrorFlat(vin, systUnfoldingV, 
						    fnameUnfoldingConstants,
						    fnameUnfoldingSystErrors);
  }

  // Add unfolding and escale systematics to the total systematic error
  for(int i=0; i<nUnfoldingBins; i++){
    voutSystErr[i] = sqrt( voutSystErr[i]*voutSystErr[i] 
			   + systUnfoldingV[i]*systUnfoldingV[i]
			   + systEscaleV[i]*systEscaleV[i]);
  }
  unfolding::deflattenMatrix(voutSystErr, voutSystErrM);

  // After propagating through unfolding all errors that we had on yields before 
  // unfolding we can compute the relative errors of each kind. While unfolding
  // changes relative errors, all subsequent manipulations do not, so we can 
  // save the errors here.
  for(int i=0; i<DYTools::nMassBins; i++){
    for (int yi=0; yi<DYTools::nYBins[i]; yi++) {
      int idx=DYTools::findIndexFlat(i,yi);
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
      for (int yi=0; yi<DYTools::nYBins[i]; ++yi) {
      int idx=DYTools::findIndexFlat(i,yi);
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

// -------------------------------------------------------------------------

inline  double asym(double x, double y) { return (x-y)/(x+y); }

// -------------------------------------------------------------------------

void applyUnfoldingToMc(int fsr) { //TString fullUnfoldingConstFileName, TString fullMcRefYieldsFileName)
  std::string line(80,'-');
  std::cout << line << "\n entered applyUnfoldingToMc(fsr=" << fsr << ")\n";


  if (fsr && (fsrCorrection_BinByBin==_fsrCorr_binByBin)) {
    std::cout << "applyUnfoldingToMc: fsr flag is set on, but fsrCorrection is bin by bin\n";
    return;
    //assert(0);
  }
  TString unfoldingConstFileName;
  TString mcRefYieldsFileName;
  TString iniArrName, finArrName;

  switch(fsr) {
  case 0: 
    unfoldingConstFileName=fnameUnfoldingConstants;
    mcRefYieldsFileName=fnameMcReferenceYields;
    iniArrName="yieldsMcPostFsrGenFIArray";
    finArrName="yieldsMcPostFsrRecFIArray";
    break;
  case 1:
  case 3:
    unfoldingConstFileName = fnameFsrCorrectionConstantsUnf;
    mcRefYieldsFileName = fnameMcReferenceYieldsFsr;
    iniArrName="yieldsMcPreFsrGenFIArray";
    finArrName="yieldsMcPostFsrGenFIArray";
    break;
  case 2:
  case 4:
    unfoldingConstFileName = fnameFsrCorrectionDETConstantsUnf;
    mcRefYieldsFileName = fnameMcReferenceYieldsFsrDET;
    iniArrName="yieldsMcPreFsrGenDETFIArray";
    finArrName="yieldsMcPostFsrGenDETFIArray";
    break;
  default:
    std::cout << "applyUnfoldingToMc(" << fsr << "): code is not ready for this value of fsr\n";
    assert(0);
  }

  HERE("Load MC reference yields");

  std::cout << "mcRefYieldsFileName=<" << mcRefYieldsFileName << ">\n";
  TFile fileMcRef(mcRefYieldsFileName);
  if (!fileMcRef.IsOpen()) {
    std::cout << "failed to open a file <" << mcRefYieldsFileName << ">\n";
    assert(fileMcRef.IsOpen());
  }
  TVectorD *yieldsMcGenPtr = (TVectorD *)fileMcRef.FindObjectAny(iniArrName);
  TVectorD *yieldsMcRecPtr = (TVectorD *)fileMcRef.FindObjectAny(finArrName);
  if (!yieldsMcRecPtr || !yieldsMcGenPtr) {
    std::cout << "null pointers from <" << mcRefYieldsFileName << ">\n";
    std::cout << "  was looking for <" << iniArrName << ">\n";
    std::cout << "  and <" << finArrName << ">\n";
    assert(0);
  }
  TVectorD yieldsMcRecV = *yieldsMcRecPtr;
  TVectorD yieldsMcGenV = *yieldsMcGenPtr;

  //std::cout << "yieldsMcPostFsrRec:\n"; yieldsMcRecV.Print();
  //std::cout << "yieldsMcPostFsrGen:\n"; yieldsMcGenV.Print();

  int nUnfoldingBins = DYTools::getTotalNumberOfBins();
  TVectorD dNdMmcCheckVgen(nUnfoldingBins), dNdMmcCheckVreco(nUnfoldingBins);
  dNdMmcCheckVgen = 0; 
  dNdMmcCheckVreco = 0;

  if ((fsrCorrection_BinByBin!=int(_fsrCorr_unf)) ||
      ((fsrCorrection_BinByBin==int(_fsrCorr_unf)) && (fsr!=2) && (fsr!=4))) {
    unfolding::unfold(yieldsMcRecV, dNdMmcCheckVgen, unfoldingConstFileName);
    unfolding::unfoldTrueToReco(yieldsMcGenV, dNdMmcCheckVreco, unfoldingConstFileName);
  }
  else {
    unfolding::unfoldFSR(yieldsMcRecV, dNdMmcCheckVgen, unfoldingConstFileName, fnameFsrDETcorrFactors);
    unfolding::unfoldFSRTrueToReco(yieldsMcGenV, dNdMmcCheckVreco, unfoldingConstFileName, fnameFsrDETcorrFactors);
  }

  int saveTestRes=0;
  if (fsrCorrection_BinByBin!=_fsrCorr_binByBin) {
    switch(fsr) {
    case 0:
      {
	saveTestRes=1;
	TriggerSelection triggers("Full2011_hltEffOld",true,0);
	TMatrixD mIni(DYTools::nMassBins,nMaxYBins);
	mIni=0;
	TMatrixD mIniErr=mIni, mIniSystErr=mIni;
	TMatrixD mFin=mIni, mFinErr=mIni, mFinSystErr=mIni;
	assert(unfolding::deflattenMatrix(yieldsMcRecV,mFin));
	
	int det=1;
	TString specTag=TString("_test_") + finArrName;
	specTag.ReplaceAll("FIArray","_orig");
	saveYields(mFin,mFinErr,mFinSystErr, triggers,specTag, det);

	assert(unfolding::deflattenMatrix(yieldsMcGenV,mIni));
	mIniErr=0; mIniSystErr=0;
	specTag=TString("_test_") + iniArrName;
	specTag.ReplaceAll("FIArray","_orig");
	saveYields(mIni,mIniErr,mIniSystErr, triggers,specTag, det);

	specTag=TString("_test_") + iniArrName;
	specTag.ReplaceAll("FIArray","_from_unf_test");
	assert(unfolding::deflattenMatrix(dNdMmcCheckVgen,mIni));
	saveYields(mIni,mIniErr,mIniSystErr, triggers,specTag, det);
      }
      break;
    case 3:
    case 4: {
      saveTestRes=1;
      TriggerSelection triggers("Full2011_hltEffOld",true,0);
      TMatrixD mIni(DYTools::nMassBins,nMaxYBins);
      mIni=0;
      TMatrixD mIniErr=mIni, mIniSystErr=mIni;
      TMatrixD mFin=mIni, mFinErr=mIni, mFinSystErr=mIni;
      assert(unfolding::deflattenMatrix(yieldsMcRecV,mFin));
      if (fsr==3) fsrCorrection(mFin,mFinErr,mFinSystErr, mIni,mIniErr,mIniSystErr);
      else fsrCorrectionSansAcceptance(mFin,mFinErr,mFinSystErr, mIni,mIniErr,mIniSystErr);
      int det=(fsr==3) ? 0 : 1;
      TString specTag=TString("_fsrTest_") + finArrName;
      specTag.ReplaceAll("FIArray","_orig");
      saveYields(mFin,mFinErr,mFinSystErr, triggers,specTag, det);

      specTag=TString("_fsrTest_") + iniArrName;
      specTag.ReplaceAll("FIArray","_from_fnc");
      saveYields(mIni,mIniErr,mIniSystErr, triggers,specTag, det);

      assert(unfolding::deflattenMatrix(yieldsMcGenV,mIni));
      mIniErr=0; mIniSystErr=0;
      specTag=TString("_fsrTest_") + iniArrName;
      specTag.ReplaceAll("FIArray","_orig");
      saveYields(mIni,mIniErr,mIniSystErr, triggers,specTag, det);

      specTag=TString("_fsrTest_") + iniArrName;
      specTag.ReplaceAll("FIArray","_from_unf_test");
      assert(unfolding::deflattenMatrix(dNdMmcCheckVgen,mIni));
      saveYields(mIni,mIniErr,mIniSystErr, triggers,specTag, det);
    }
      break;
    default:
      saveTestRes=0;
    }
  }


  // Report results
  if (1) {
    printf("\nUNFOLD: Check on the MC, yields reco->gen:\n");
    printf(" mass      rapidity    yieldsReco   yieldsGen     unfolded  asym(%%)\n");
    for(int i=0; i<DYTools::nMassBins; i++){
      double *rapidityBinLimits=DYTools::getYBinLimits(i);
      for (int yi=0; yi<DYTools::nYBins[i]; ++yi) {
	int idx=DYTools::findIndexFlat(i,yi);
	printf("%4.0f-%4.0f  %4.2lf-%4.2lf  %10.0f  %10.0f    %10.0f   %6.2f\n",
	       DYTools::massBinLimits[i],DYTools::massBinLimits[i+1],
	       rapidityBinLimits[yi],rapidityBinLimits[yi+1],
	       yieldsMcRecV[idx],yieldsMcGenV[idx],
	       dNdMmcCheckVgen[idx],
	       100*asym(yieldsMcGenV[idx],dNdMmcCheckVgen[idx])
	       );
      }
      delete rapidityBinLimits;
    }

    printf("\nUNFOLD: Check on the MC, yields gen->reco:\n");
    printf(" mass      rapidity    yieldsGen   yieldsReco    unfolded  asym(%%)\n");
    for(int i=0; i<DYTools::nMassBins; i++){
      double *rapidityBinLimits=DYTools::getYBinLimits(i);
      for (int yi=0; yi<DYTools::nYBins[i]; ++yi) {
        int idx=DYTools::findIndexFlat(i,yi);
        printf("%4.0f-%4.0f  %4.2lf-%4.2lf  %10.0f   %10.0f    %10.0f  %6.2lf\n",
               DYTools::massBinLimits[i],DYTools::massBinLimits[i+1],
               rapidityBinLimits[yi],rapidityBinLimits[yi+1],
	       yieldsMcGenV[idx],
               yieldsMcRecV[idx],
               dNdMmcCheckVreco[idx],
	       100*asym(yieldsMcRecV[idx],dNdMmcCheckVreco[idx])
	       );
      }
      delete rapidityBinLimits;
    }
  }

  fileMcRef.Close();
  if (yieldsMcGenPtr) delete yieldsMcGenPtr;
  if (yieldsMcRecPtr) delete yieldsMcRecPtr;

  std::cout << " leaving applyUnfoldingToMc(fsr=" << fsr << ")\n" << line << "\n";

}

// -------------------------------------------------------------------------

void  efficiencyCorrection(const TMatrixD &vin, const TMatrixD &vinStatErr, const TMatrixD &vinSystErr,
             TMatrixD &vout, TMatrixD &voutStatErr, TMatrixD &voutSystErr)
{
  const int nUnfoldingBins= DYTools::getTotalNumberOfBins();

  // Read efficiency constants
  HERE("Efficiency: Load constants");
    
  TFile fileConstants(fnameEfficiencyConstants);
  TMatrixD* efficiencyArrayPtr    = (TMatrixD *)fileConstants.FindObjectAny("efficiencyArray");
  TMatrixD* efficiencyErrArrayPtr = (TMatrixD *)fileConstants.FindObjectAny("efficiencyErrArray");

  if (!efficiencyArrayPtr || !efficiencyErrArrayPtr) {
    std::cout << "at least one needed object is not present in <" << fnameEfficiencyConstants << ">\n";
    assert(0);
  }

  TFile fileScaleConstants(fnameScaleFactorConstants);
  TVectorD* rhoDataMcPtr    = (TVectorD *)fileScaleConstants.FindObjectAny("scaleFactorFlatIdxArray");
  TVectorD* rhoDataMcErrPtr = (TVectorD *)fileScaleConstants.FindObjectAny("scaleFactorErrFlatIdxArray");

  if (!rhoDataMcPtr || !rhoDataMcErrPtr) {
    std::cout << "a least one needed object is not present in <" << fnameScaleFactorConstants << ">\n";
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

  if (0) { // localDebug
    TMatrixD effScaleF(DYTools::nMassBins,nMaxYBins);
    TMatrixD effScaleFErr(DYTools::nMassBins,nMaxYBins);
    for (int mi=0, idx=0; mi<DYTools::nMassBins; ++mi) {
      for (int yi=0; yi<DYTools::nYBins[mi]; ++yi, ++idx) {
	effScaleF(mi,yi) = efficiencyArray[mi][yi] * rhoDataMc[idx];
	effScaleFErr(mi,yi) = sqrt(efficiencyErrArray[mi][yi]*efficiencyErrArray[mi][yi] + rhoDataMcErr[idx]*rhoDataMcErr[idx]);
      }
    }
    TString chkFName="effScale_" + DYTools::analysisTag + ".root";
    TFile chk(chkFName,"RECREATE");
    effScaleF.Write("effScaleF");
    effScaleFErr.Write("effScaleFErr");
    chk.Close();
  }

  // Apply the correction
  TMatrixD systErrorPropagated(DYTools::nMassBins,nMaxYBins);
  TMatrixD systErrorAdditional(DYTools::nMassBins,nMaxYBins);
  systErrorAdditional = 0;

  for(int mi=0, idx=0; mi<DYTools::nMassBins; mi++){
    for (int yi=0; yi<DYTools::nYBins[mi]; ++yi,++idx) {
      double effFactor = efficiencyArray[mi][yi] * rhoDataMc[idx];
      double effErr = effFactor
	* sqrt( efficiencyErrArray[mi][yi]*efficiencyErrArray[mi][yi]/efficiencyArray[mi][yi]/efficiencyArray[mi][yi]
		+ rhoDataMcErr[idx]*rhoDataMcErr[idx]/rhoDataMc[idx]/rhoDataMc[idx]);
      vout[mi][yi]        = vin[mi][yi] / effFactor;
      // Statistical error propagated
      voutStatErr[mi][yi] = vinStatErr[mi][yi] / effFactor;
      // Old systematic error, propagated
      systErrorPropagated[mi][yi] = vinSystErr[mi][yi] / effFactor;
      // Extra systematic error due to the errors on the efficiency and scale factor
      systErrorAdditional[mi][yi] = (vin[mi][yi]/effFactor)*(effErr/effFactor);
      systEfficiency[mi][yi] = systErrorAdditional[mi][yi]/vout[mi][yi];
      voutSystErr[mi][yi] = 
	sqrt(systErrorPropagated[mi][yi]*systErrorPropagated[mi][yi]
	     + systErrorAdditional[mi][yi]*systErrorAdditional[mi][yi]);
    }
  }

  if (printEfficiencyTable) {
  printf("\nEfficiency: Results for the data, yields:\n");
  printf("                after unfolding            eff. factors,%%       rho(data/mc)         efficiency-corrected        syst-err-eff, %%\n");
  for(int mi=0, idx=0; mi<DYTools::nMassBins; mi++){
    double *rapidityBinLimits=DYTools::getYBinLimits(mi);
    for (int yi=0; yi<DYTools::nYBins[mi]; ++yi,++idx) {
      printf("%4.0f-%4.0f  %4.2lf-%4.2lf  %8.1f +- %7.1f +- %7.1f   %5.2f +- %5.2f      %4.3f +- %4.3f      %8.1f +- %7.1f +- %7.1f        %4.2f\n",
	     DYTools::massBinLimits[mi],DYTools::massBinLimits[mi+1],
	     rapidityBinLimits[yi],rapidityBinLimits[yi+1],
	     vin[mi][yi], vinStatErr[mi][yi], vinSystErr[mi][yi],
	     efficiencyArray[mi][yi]*100, efficiencyErrArray[mi][yi]*100, 
	     rhoDataMc[idx], rhoDataMcErr[idx],
	     vout[mi][yi], voutStatErr[mi][yi], voutSystErr[mi][yi],
	     systErrorAdditional[mi][yi]*100.0/vout[mi][yi]);
    }
    delete rapidityBinLimits;
  }
  }

  fileConstants.Close();
  fileScaleConstants.Close();

  return;

}

// -------------------------------------------------------------------------

void  acceptanceCorrection(const TMatrixD &vin, const TMatrixD &vinStatErr, const TMatrixD &vinSystErr,
			   TMatrixD &vout, TMatrixD &voutStatErr, TMatrixD &voutSystErr){
  int nUnfoldingBins = DYTools::getTotalNumberOfBins();

  // Read efficiency constants
  printf("Acceptance: Load constants\n"); fflush(stdout);
    
  TFile fileConstants(fnameAcceptanceConstants);
  if (!fileConstants.IsOpen()) {
    std::cout << "failed to open acceptance systematics file <" << fnameAcceptanceConstants << ">\n";
    assert(0);
  }
  TMatrixD *acceptanceMatrixPtr    = (TMatrixD *)fileConstants.FindObjectAny("acceptanceMatrix");
  TMatrixD *acceptanceErrMatrixPtr = (TMatrixD *)fileConstants.FindObjectAny("acceptanceErrMatrix");
  if (!acceptanceMatrixPtr || !acceptanceErrMatrixPtr) {
    std::cout << "at least one object from file <" << fnameAcceptanceConstants << "> is null\n";
    assert(0);
  }

  TMatrixD acceptanceMatrix = *acceptanceMatrixPtr;
  TMatrixD acceptanceErrMatrix = *acceptanceErrMatrixPtr;

  TVectorD acceptanceTheoryErrArray(nUnfoldingBins);
  acceptanceTheoryErrArray=0;

  if (DYTools::study2D==0) {
    TFile fileSystematics(fnameAcceptanceSystematics);
			
    TVectorD *acceptanceTheoryErrArrayPtr = (TVectorD *)fileSystematics.FindObjectAny("acceptanceTheoryErrArray");
    if (!acceptanceTheoryErrArrayPtr) {
      std::cout << "failed to get object from <" << fnameAcceptanceSystematics << ">\n";
      if (0) {
	assert(0);
      }
      else {
	std::cout << "\n\tfailing acceptanceTheoryErrArray from <" << fnameAcceptanceSystematics << ">\n";
	vout=0; voutStatErr=0; voutSystErr=0;
	if (DYTools::study2D==0) assert(0);
      }
    }
    acceptanceTheoryErrArray = *acceptanceTheoryErrArrayPtr;
    fileSystematics.Close();
  }

  TFile fileAccFSRSyst(fnameAcceptanceFSRSystematics);
  TMatrixD *acceptanceFSRErrMatrixPtr = (TMatrixD *)fileAccFSRSyst.FindObjectAny("accSystPercent");
  if (!acceptanceFSRErrMatrixPtr) {
    std::cout << "failed to get object from <" << fnameAcceptanceFSRSystematics << ">\n";
    if (0) {
      assert(0);
    }
    else {
      std::cout << "\n\tfailing acceptanceCorrection (FSR systematics)\n";
      vout=0; voutStatErr=0; voutSystErr=0;
      //if (DYTools::study2D==0) assert(0);
    }
  }
  TMatrixD acceptanceFSRErrMatrix=acceptanceMatrix;
  if (acceptanceFSRErrMatrixPtr) {
    acceptanceFSRErrMatrix=*acceptanceFSRErrMatrixPtr;
  }
  else acceptanceFSRErrMatrix=0;

  // Check that the binning is consistent
  bool checkResult = true;
  if (checkResult) checkResult=unfolding::checkRangesMY(acceptanceMatrix,"acceptanceMatrix");
  //if( acceptanceMatrix.GetNoElements() != DYTools::nMassBins) checkResult = false;
  if(( acceptanceTheoryErrArray.GetNoElements() != nUnfoldingBins) &&
     ( acceptanceTheoryErrArray.GetNoElements() != DYTools::nUnfoldingBinsMax) ) {
    checkResult = false;
  }

  if( !checkResult ){
    printf("Acceptance: ERROR: inconsistent binning in the inputs\n");
    std::cout << "acceptanceTheoryErrArray.size=" << acceptanceTheoryErrArray.GetNoElements() << " (allowed=" << nUnfoldingBins << " or " << DYTools::nUnfoldingBinsMax << ")\n";
    assert(0);
  }else
    printf("Acceptance: Binning in the inputs is consistent\n");

  // Apply the correction
  TMatrixD systErrorPropagated(DYTools::nMassBins,nMaxYBins);
  TMatrixD systErrorAdditional(DYTools::nMassBins,nMaxYBins);
  TMatrixD systErrorTheory(DYTools::nMassBins,nMaxYBins);
  systErrorAdditional = 0;

  for(int i=0, idx=0; i<DYTools::nMassBins; i++){
    for (int yi=0; yi<DYTools::nYBins[i]; ++yi, ++idx) {
      double accFactor = acceptanceMatrix[i][yi];
      double accErr    = acceptanceErrMatrix[i][yi];
      double accThErr  = acceptanceTheoryErrArray[idx];
      double accFSRErr = acceptanceFSRErrMatrix[i][yi]/100;
      systAccTheoryRelative[i][yi]=acceptanceTheoryErrArray[idx];
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
    for(int i=0, idx=0; i<DYTools::nMassBins; i++){
      double *rapidityBinLimits=DYTools::getYBinLimits(i);
      for (int yi=0; yi<DYTools::nYBins[i]; ++yi, ++idx) {
	printf("%4.0f-%4.0f  %4.2lf-%4.2lf  %8.1f +- %7.1f +- %7.1f   %5.2f +- %5.2f      %9.1f +- %8.1f +- %8.1f        %4.2f     %4.2f\n",
	       DYTools::massBinLimits[i],DYTools::massBinLimits[i+1],
	       rapidityBinLimits[yi],rapidityBinLimits[yi+1],
	       vin[i][yi], vinStatErr[i][yi], vinSystErr[i][yi],
	       acceptanceMatrix[i][yi]*100, acceptanceErrMatrix[i][yi]*100, 
	       vout[i][yi], voutStatErr[i][yi], voutSystErr[i][yi],
	       systErrorAdditional[i][yi]*100/vout[i][yi],100*acceptanceTheoryErrArray[idx]);
      }
    }
  }

  fileConstants.Close();

  return;
}

// -------------------------------------------------------------------------

void  fsrCorrectionViaUnfolding(const TMatrixD &vinM, const TMatrixD &vinStatErrM, const TMatrixD &vinSystErrM,
				TMatrixD &voutM, TMatrixD &voutStatErrM, TMatrixD &voutSystErrM,
				const TString &unfFileName, const TString &name)
{
  int nUnfoldingBins = DYTools::getTotalNumberOfBins();
  TVectorD vin(nUnfoldingBins),vinStatErr(nUnfoldingBins),vinSystErr(nUnfoldingBins);
  TVectorD vout(nUnfoldingBins),voutStatErr(nUnfoldingBins),voutSystErr(nUnfoldingBins);

  // First, propagate through unfolding the signal yields with stat and syst errors
  int res=0;

  switch(fsrCorrection_BinByBin) {
  case _fsrCorr_binByBin: res=0; break;
  case _fsrCorr_unfPure:
  case _fsrCorr_unfMdf:
  case _fsrCorr_unfGood:
    res=( (unfolding::unfold(vinM, voutM, unfFileName, vin, vout) == 1 ) &&
	  (unfolding::propagateErrorThroughUnfolding(vinStatErrM,voutStatErrM, unfFileName, vinStatErr, voutStatErr) == 1 ) &&
	  (unfolding::propagateErrorThroughUnfolding(vinSystErrM,voutSystErrM, unfFileName, vinSystErr, voutSystErr) == 1) ) ? 1:0;
    break;
  case _fsrCorr_unf:
    res=( (unfolding::unfoldFSR(vinM, voutM, unfFileName, fnameFsrDETcorrFactors, vin, vout) == 1 ) &&
	  (unfolding::propagateErrorThroughFsrUnfolding(vinStatErrM,voutStatErrM, unfFileName, fnameFsrDETcorrFactors, vinStatErr, voutStatErr) == 1) &&
	  (unfolding::propagateErrorThroughFsrUnfolding(vinSystErrM,voutSystErrM, unfFileName, fnameFsrDETcorrFactors, vinSystErr, voutSystErr) ==1 ) ) ? 1:0;
    break;
  default:
    std::cout << "fsrCorrectionViaUnfolding is not ready for this _fsrCorr_* case\n";
    res=0;
  }
  if (!res) {
    std::cout << "fsrCorrectionViaUnfolding: unfolding procedure failed for name=<" << name << ">\n";
    assert(0);
  }
}

// -------------------------------------------------------------------------

void  fsrCorrectionBase(const TMatrixD &vin, const TMatrixD &vinStatErr, const TMatrixD &vinSystErr,
			TMatrixD &vout, TMatrixD &voutStatErr, TMatrixD &voutSystErr,
			const TString &fname, const TString &title)
{
  // Read efficiency constants
  std::cout << "FsrCorrectionBase(" << title << "): Load constants" << std::endl;
    
  TFile fileConstants(fname);
  assert(fileConstants.IsOpen());
  TMatrixD *fsrCorrectionMatrixPtr  = (TMatrixD *)fileConstants.FindObjectAny("fsrCorrectionMatrix");
  TMatrixD *fsrCorrectionErrMatrixPtr = (TMatrixD *)fileConstants.FindObjectAny("fsrCorrectionErrMatrix");
  assert(fsrCorrectionMatrixPtr && fsrCorrectionErrMatrixPtr);
  TMatrixD fsrCorrectionMatrix= *fsrCorrectionMatrixPtr;
  TMatrixD fsrCorrectionErrMatrix= *fsrCorrectionErrMatrixPtr;

  printSanityCheck(fsrCorrectionMatrix,fsrCorrectionErrMatrix,title);

  // Check that the binning is consistent
  bool checkResult = true;
  if (checkResult) checkResult=unfolding::checkRangesMY(fsrCorrectionMatrix,"fsrCorrectionMatrix");
  //if( fsrCorrectionArray.GetNoElements() != DYTools::nMassBins) checkResult = false;
  if( !checkResult ){
    printf("FsrCorrectionBase: ERROR: inconsistent binning in the inputs\n");
    assert(0);
  }else
    printf("FsrCorrectionBase: Binning in the inputs is consistent\n");

  // Apply the correction

  TMatrixD systErrorPropagated(DYTools::nMassBins,nMaxYBins);
  for(int i=0; i<DYTools::nMassBins; i++){
    for (int yi=0; yi<DYTools::nYBins[i]; ++yi) {
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
  for(int i=0; i<DYTools::nMassBins; i++){
    double *rapidityBinLimits=DYTools::getYBinLimits(i);
    for (int yi=0; yi<DYTools::nYBins[i]; ++yi) {
      printf("%4.0f-%4.0f  %4.2f-%4.2f  %9.1f +- %8.1f +- %8.1f   %4.3f +- %4.3f      %9.1f +- %8.1f +- %8.1f        %5.2f\n",
	     DYTools::massBinLimits[i],DYTools::massBinLimits[i+1],
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

// -------------------------------------------------------------------------

void  fsrCorrection(const TMatrixD &vin, const TMatrixD &vinStatErr, const TMatrixD &vinSystErr,
		    TMatrixD &vout, TMatrixD &voutStatErr, TMatrixD &voutSystErr)
{
  if (fsrCorrection_BinByBin==_fsrCorr_binByBin) {
    fsrCorrectionBase(vin,vinStatErr,vinSystErr, vout,voutStatErr,voutSystErr,
		      fnameFsrCorrectionConstantsBbB,"FsrCorrectionBbB");
  }
  else {
    fsrCorrectionViaUnfolding(vin,vinStatErr,vinSystErr, vout,voutStatErr,voutSystErr,
			      fnameFsrCorrectionConstantsUnf, "FsrCorrectionUnf");
  }

  return;
}

// -------------------------------------------------------------------------

void  fsrCorrectionSansAcceptance(const TMatrixD &vin, const TMatrixD &vinStatErr, const TMatrixD &vinSystErr,
				  TMatrixD &vout, TMatrixD &voutStatErr, TMatrixD &voutSystErr)
{
  if (fsrCorrection_BinByBin==_fsrCorr_binByBin) {
    fsrCorrectionBase(vin,vinStatErr,vinSystErr, vout,voutStatErr,voutSystErr,
		      fnameFsrCorrectionSansAccConstantsBbB, "FsrCorrectionSansAcc");
  }
  else {
    fsrCorrectionViaUnfolding(vin,vinStatErr,vinSystErr, vout,voutStatErr,voutSystErr,
			      fnameFsrCorrectionDETConstantsUnf, "FsrCorrectionDET");
  }

  return;
}

// -------------------------------------------------------------------------

void  crossSections(const TMatrixD &vin, const TMatrixD &vinStatErr, const TMatrixD &vinSystErr,
		    TMatrixD &vout, TMatrixD &voutStatErr, TMatrixD &voutSystErr,
		    TMatrixD &voutNorm, TMatrixD &voutNormStatErr, TMatrixD &voutNormSystErr,
		    const TriggerSelection &triggers, const TString specTag)
{
  assert(triggers.isDefined()); // eliminate compiler warning

  // Find absolute cross-section
  for(int i=0; i<DYTools::nMassBins; i++) {
    for (int yi=0; yi<DYTools::nYBins[i]; ++yi) {
      vout[i][yi] = vin[i][yi] / lumi;
      voutStatErr[i][yi] = vinStatErr[i][yi] / lumi;
      voutSystErr[i][yi] = vinSystErr[i][yi] / lumi;
      //std::cout << "i=" << i << ", yi=" << yi << ": vinSystErr = " << vinSystErr[i][yi] << "\n";
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
  int lowYBins=DYTools::nYBins[low];

  for( int i=low; i<=high; i++){
    if (lowYBins!=DYTools::nYBins[i]) {
      std::cout << "y binning error in crossSections. Correct the code: wrong assumption\n";
      assert(0);
    }
    for (int yi=0; yi<DYTools::nYBins[i]; ++yi) {
      xsecReference += vout[i][yi];
      xsecReferenceStatErr += voutStatErr[i][yi] * voutStatErr[i][yi];
      xsecReferenceSystErr += voutSystErr[i][yi] * voutSystErr[i][yi];
      //std::cout << "i=" << i << ", yi=" << yi << ": { xsecReferenceSystErr += " << voutSystErr[i][yi] << " } = " << xsecReferenceSystErr << "\n";
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

  bool use1binErrors = false;
  if( use1binErrorsForNorm == 1 && DYTools::study2D==0 && specTag == "" ) {
    use1binErrors = true;
    xsecReferenceStatErr = normAbsStatErrPreFsrFullAcc;
    xsecReferenceSystErr = normAbsSystErrPreFsrFullAcc;
  }

  // Find normalized cross-section
  for(int i=0; i<DYTools::nMassBins; i++) {
    for (int yi=0; yi<DYTools::nYBins[i]; ++yi) {
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
    //printf("                    absolute                       normalized +- stat +- sys (total)           (1/sigma)(1/dM)norm +-stat +-syst (total) \n");
    printf("#                   (1/sigma)(1/dM)norm +-stat +-syst  +-total \n");
  }
  for(int i=0; i<DYTools::nMassBins; i++){
    double *rapidityBinLimits=DYTools::getYBinLimits(i);
    for (int yi=0; yi<DYTools::nYBins[i]; ++yi) {
      double binw = DYTools::massBinLimits[i+1] - DYTools::massBinLimits[i];
      normXSec[i][yi]=voutNorm[i][yi];
      normXSecErrStat[i][yi]=voutNormStatErr[i][yi];
      normXSecErrSyst[i][yi]=voutNormSystErr[i][yi];
      normXSecErr[i][yi]=sqrt( SQR(voutNormStatErr[i][yi]) + SQR(voutNormSystErr[i][yi]) );
      
      normXSecByBin[i][yi]=voutNorm[i][yi]/binw;
      normXSecErrByBinStat[i][yi]=voutNormStatErr[i][yi]/binw;
      normXSecErrByBinSyst[i][yi]=voutNormSystErr[i][yi]/binw;
      normXSecErrByBin[i][yi]=sqrt( SQR(voutNormStatErr[i][yi]) + SQR(voutNormSystErr[i][yi]) )/binw;

      if (printPreFSRCrossSectionTable && (DYTools::study2D==0)) {
	if (specTag.Length() && (i==0) && (yi==0)) std::cout << "special tag=<" << specTag << ">\n";
	printf("%4.0f %4.0f      %1.8e %1.8e %1.8e   %1.8e     \n",
	       DYTools::massBinLimits[i],DYTools::massBinLimits[i+1],
	       voutNorm[i][yi]/binw, voutNormStatErr[i][yi]/binw, voutNormSystErr[i][yi]/binw,
	       normXSecErrByBin[i][yi]
	       );
      }

      if (0)
      if (printPreFSRCrossSectionTable) {
	if (specTag.Length() && (i==0) && (yi==0)) std::cout << "special tag=<" << specTag << ">\n";
	printf("%4.0f-%4.0f  %4.2f-%4.2f    %6.1f +- %4.1f +- %4.1f      %1.6f +- %1.6f +- %1.6f  ( %1.6f )     %1.8f +- %1.8f +- %1.8f  ( %1.8f )    \n",
	       DYTools::massBinLimits[i],DYTools::massBinLimits[i+1],
	       rapidityBinLimits[yi],rapidityBinLimits[yi+1],
	       vout[i][yi], voutStatErr[i][yi], voutSystErr[i][yi],
	       voutNorm[i][yi], voutNormStatErr[i][yi], voutNormSystErr[i][yi],
	       normXSecErr[i][yi],
	       voutNorm[i][yi]/binw, voutNormStatErr[i][yi]/binw, voutNormSystErr[i][yi]/binw,
	       normXSecErrByBin[i][yi]
	       );
      }
    }
    delete rapidityBinLimits;
  }

  gSystem->mkdir(pathXSect,kTRUE);
  TString extraTag=DYTools::analysisTag;
  if (specTag.Length()) extraTag.Append(TString("_") + specTag);
  TString xSecResultFileName(pathXSect+TString("xSec_results_") + 
			     extraTag + TString(".root"));
  std::cout << "xSecResultFileName= " << xSecResultFileName << "\n";
  TFile fa(xSecResultFileName,"recreate");
  vin.Write("XSec");
  vinStatErr.Write("XSecErr");
  vinSystErr.Write("XSecSystErr");
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


  printf("\nPre FSR cross-section in the Z peak from %3.0f to %3.0f:\n",
	 DYTools::massBinLimits[low], DYTools::massBinLimits[high+1]);
  printf("           %9.1f +- %8.1f +- %6.1f \n",
	 xsecReference, xsecReferenceStatErr, xsecReferenceSystErr);
  if( use1binErrors ){
    printf("   INFO: the stat and syst errors come from 1-bin 60-120 GeV measurement\n");
    printf("         hardwired in the code. Make sure the numbers are up to date.\n");
  }else{
    printf("   INFO: the stat and syst errors come from the sum of 13 bins in the range 60-120 GeV\n");
  }
//     printf("check %f %f\n", xsecReferenceStatErr, xsecReferenceSystErr);

  return;
}

// -------------------------------------------------------------------------

void  crossSectionsDET(const TMatrixD &vin, const TMatrixD &vinStatErr, const TMatrixD &vinSystErr,
		       TMatrixD &vout, TMatrixD &voutStatErr, TMatrixD &voutSystErr,
		       TMatrixD &voutNorm, TMatrixD &voutNormStatErr, 
		       TMatrixD &voutNormSystErr,
		       const TriggerSelection &triggers, const TString specTag)
{
  assert(triggers.isDefined()); // eliminate compiler warning

  // Find absolute cross-section
  for(int i=0; i<DYTools::nMassBins; i++) {
    for (int yi=0; yi<DYTools::nYBins[i]; ++yi) {
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
    for (int yi=0; yi<DYTools::nYBins[i]; ++yi) {
      xsecReference += vout[i][yi];
      xsecReferenceStatErr += voutStatErr[i][yi] * voutStatErr[i][yi];
      xsecReferenceSystErr += voutSystErr[i][yi] * voutSystErr[i][yi];
    }
  }

  xsecReferenceStatErr = sqrt(xsecReferenceStatErr);
  xsecReferenceSystErr = sqrt(xsecReferenceSystErr);

  bool use1binErrors = false;
  if( use1binErrorsForNorm == 1 && DYTools::study2D==0 && specTag == "" ) {
    use1binErrors = true;
    xsecReferenceStatErr = normAbsStatErrPreFsrDET;
    xsecReferenceSystErr = normAbsSystErrPreFsrDET;
  }

  // Find normalized cross-section
  for(int i=0; i<DYTools::nMassBins; i++) {
    for (int yi=0; yi<DYTools::nYBins[i]; ++yi) {
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
    if (specTag.Length()) std::cout << "special tag=<" << specTag << ">\n";
    printf("                    absolute                      normalized +-stat +-sys (total)\n");
    for(int i=0; i<DYTools::nMassBins; i++){
      double *rapidityBinLimits=DYTools::getYBinLimits(i);
      for (int yi=0; yi<DYTools::nYBins[i]; ++yi) {
	printf("%4.0f-%4.0f  %4.2f-%4.2f    %6.1f +- %4.1f +- %4.1f      %1.6f +- %1.6f +- %1.6f  ( %1.6f )\n",
	       //    printf("%4.0f-%4.0f   %9.1f +- %8.1f +- %8.1f   %1.6f +- %1.6f +- %1.6f  ( %1.6f )\n",
	       DYTools::massBinLimits[i],DYTools::massBinLimits[i+1],
	       rapidityBinLimits[yi], rapidityBinLimits[yi+1],
	       vout[i][yi], voutStatErr[i][yi], voutSystErr[i][yi],
	       voutNorm[i][yi], voutNormStatErr[i][yi], voutNormSystErr[i][yi],
	       sqrt(SQR(voutNormStatErr[i][yi]) + SQR(voutNormSystErr[i][yi])) );
      }
      delete rapidityBinLimits;
    }
  }
  printf("\nPreFsr cross-section in the Z peak from %3.0f to %3.0f:\n",
	 DYTools::massBinLimits[low], DYTools::massBinLimits[high+1]);
  if (specTag.Length()) std::cout << "special tag=<" << specTag << ">\n";
  printf("           %9.1f +- %8.1f +- %6.1f \n",
	 xsecReference, xsecReferenceStatErr, xsecReferenceSystErr);
  if( use1binErrors ){
    printf("   INFO: the stat and syst errors come from 1-bin 60-120 GeV measurement\n");
    printf("         hardwired in the code. Make sure the numbers are up to date.\n");
  }else{
    printf("   INFO: the stat and syst errors come from the sum of 13 bins in the range 60-120 GeV\n");
  }


  gSystem->mkdir(pathXSect,kTRUE);
  TString extraTag=DYTools::analysisTag;
  if (specTag.Length()) extraTag.Append(TString("_") + specTag);
  TString xSecResultFileName(pathXSect+TString("xSecDET_results_") + 
			     extraTag + TString(".root"));
  std::cout << "xSecDETResultFileName= " << xSecResultFileName << "\n";

  TMatrixD normXSec=voutNorm;
  TMatrixD normXSecErr=voutNormStatErr; // has to be calculated
  normXSecErr=0;
  for (int iM=0; iM<DYTools::nMassBins; ++iM) {
    for (int iY=0; iY<DYTools::nYBins[iM]; ++iY) {
      normXSecErr[iM][iY]=sqrt( SQR(voutNormStatErr[iM][iY]) + SQR(voutNormSystErr[iM][iY]) );
    }
  }
  TMatrixD normXSecErrStat=voutNormStatErr;
  TMatrixD normXSecErrSyst=voutNormSystErr;

  TFile fa(xSecResultFileName,"recreate");
  vin.Write("XSec");
  vinStatErr.Write("XSecErr");
  vinSystErr.Write("XSecSystErr");
  normXSec.Write("normXSec");
  normXSecErr.Write("normXSecErr");
  normXSecErrStat.Write("normXSecErrStat");
  normXSecErrSyst.Write("normXSecErrSyst");
  //normXSecByBin.Write("normXSecByBin");
  //normXSecErrByBin.Write("normXSecErrByBin");
  //normXSecErrByBinStat.Write("normXSecErrByBinStat");
  //normXSecErrByBinSyst.Write("normXSecErrByBinSyst");
  unfolding::writeBinningArrays(fa);
  fa.Close();

  return;
}

// -------------------------------------------------------------------------

void  postFsrCrossSections(const TMatrixD &vin, const TMatrixD &vinStatErr, const TMatrixD &vinSystErr,
		    TMatrixD &vout, TMatrixD &voutStatErr, TMatrixD &voutSystErr,
			   TMatrixD &voutNorm, TMatrixD &voutNormStatErr, TMatrixD &voutNormSystErr,
			   const TriggerSelection &triggers, const TString specTag)
{

  // Find absolute cross-section
  for(int i=0; i<DYTools::nMassBins; i++) {
    for (int yi=0; yi<DYTools::nYBins[i]; ++yi) {
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
    for (int yi=0; yi<DYTools::nYBins[i]; ++yi) {
      xsecReference += vout[i][yi];
      xsecReferenceStatErr += voutStatErr[i][yi] * voutStatErr[i][yi];
      xsecReferenceSystErr += voutSystErr[i][yi] * voutSystErr[i][yi];
    }
  }

  xsecReferenceStatErr = sqrt(xsecReferenceStatErr);
  xsecReferenceSystErr = sqrt(xsecReferenceSystErr);

  bool use1binErrors = false;
  if( use1binErrorsForNorm == 1 && DYTools::study2D==0 && specTag == "" ) {
    use1binErrors = true;
    xsecReferenceStatErr = normAbsStatErrPostFsrFullAcc;
    xsecReferenceSystErr = normAbsSystErrPostFsrFullAcc;
  }

  // Find normalized cross-section
  for(int i=0; i<DYTools::nMassBins; i++) {
    for (int yi=0; yi<DYTools::nYBins[i]; ++yi) {
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
    for(int i=0; i<DYTools::nMassBins; i++){
      double *rapidityBinLimits=DYTools::getYBinLimits(i);
      for (int yi=0; yi<DYTools::nYBins[i]; ++yi) {
	printf("%4.0f-%4.0f  %4.2f-%4.2f    %6.1f +- %4.1f +- %4.1f      %1.6f +- %1.6f +- %1.6f  ( %1.6f )\n",
	       //    printf("%4.0f-%4.0f   %9.1f +- %8.1f +- %8.1f   %1.6f +- %1.6f +- %1.6f  ( %1.6f )\n",
	       DYTools::massBinLimits[i],DYTools::massBinLimits[i+1],
	       rapidityBinLimits[yi], rapidityBinLimits[yi+1],
	       vout[i][yi], voutStatErr[i][yi], voutSystErr[i][yi],
	       voutNorm[i][yi], voutNormStatErr[i][yi], voutNormSystErr[i][yi],
	       sqrt(SQR(voutNormStatErr[i][yi]) + SQR(voutNormSystErr[i][yi])) );
      }
      delete rapidityBinLimits;
    }
  }
  printf("\nPostFsr cross-section in the Z peak from %3.0f to %3.0f:\n",
	 DYTools::massBinLimits[low], DYTools::massBinLimits[high+1]);
  printf("           %9.1f +- %8.1f +- %6.1f \n",
	 xsecReference, xsecReferenceStatErr, xsecReferenceSystErr);
//     printf("check %f %f\n", xsecReferenceStatErr, xsecReferenceSystErr);
  if( use1binErrors ){
    printf("   INFO: the stat and syst errors come from 1-bin 60-120 GeV measurement\n");
    printf("         hardwired in the code. Make sure the numbers are up to date.\n");
  }else{
    printf("   INFO: the stat and syst errors come from the sum of 13 bins in the range 60-120 GeV\n");
  }


  gSystem->mkdir(pathXSect,kTRUE);
  TString extraTag=DYTools::analysisTag;
  if (specTag.Length()) extraTag.Append(TString("_") + specTag);
  TString xSecResultFileName(pathXSect+TString("xSecPostFsr_results_") + 
			     extraTag +
			     triggers.triggerConditionsName() + TString(".root"));
  std::cout << "xSecDETResultFileName= " << xSecResultFileName << "\n";

  TMatrixD normXSec=voutNorm;
  TMatrixD normXSecErr=voutNormStatErr; // has to be calculated
  normXSecErr=-1;
  TMatrixD normXSecErrStat=voutNormStatErr;
  TMatrixD normXSecErrSyst=voutNormSystErr;

  TFile fa(xSecResultFileName,"recreate");
  normXSec.Write("normXSec");
  normXSecErr.Write("normXSecErr");
  normXSecErrStat.Write("normXSecErrStat");
  normXSecErrSyst.Write("normXSecErrSyst");
  //normXSecByBin.Write("normXSecByBin");
  //normXSecErrByBin.Write("normXSecErrByBin");
  //normXSecErrByBinStat.Write("normXSecErrByBinStat");
  //normXSecErrByBinSyst.Write("normXSecErrByBinSyst");
  unfolding::writeBinningArrays(fa);
  fa.Close();


  return;
}

// -------------------------------------------------------------------------

void  postFsrCrossSectionsDET(const TMatrixD &vin, const TMatrixD &vinStatErr, const TMatrixD &vinSystErr,
			      TMatrixD &vout, TMatrixD &voutStatErr, TMatrixD &voutSystErr,
			      TMatrixD &voutNorm, TMatrixD &voutNormStatErr, TMatrixD &voutNormSystErr,
			      const TriggerSelection &triggers, const TString specTag)
{
  
  // Find absolute cross-section
  for(int i=0; i<DYTools::nMassBins; i++) {
    for (int yi=0; yi<DYTools::nYBins[i]; ++yi) {
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
    for ( int yi=0; yi<DYTools::nYBins[i]; ++yi) {
      xsecReference += vout[i][yi];
      xsecReferenceStatErr += voutStatErr[i][yi] * voutStatErr[i][yi];
      xsecReferenceSystErr += voutSystErr[i][yi] * voutSystErr[i][yi];
    }
  }
  xsecReferenceStatErr = sqrt(xsecReferenceStatErr);
  xsecReferenceSystErr = sqrt(xsecReferenceSystErr);

  bool use1binErrors = false;
  if( use1binErrorsForNorm == 1 && DYTools::study2D==0 && specTag == "" ) {
    use1binErrors = true;
    xsecReferenceStatErr = normAbsStatErrPostFsrDET;
    xsecReferenceSystErr = normAbsSystErrPostFsrDET;
  }

  // Find normalized cross-section
  for(int i=0; i<DYTools::nMassBins; i++) {
    for (int yi=0; yi<DYTools::nYBins[i]; ++yi) {
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
    for(int i=0; i<DYTools::nMassBins; i++){
      double *rapidityBinLimits=DYTools::getYBinLimits(i);
      for (int yi=0; yi<DYTools::nYBins[i]; ++yi) {
	printf("%4.0f-%4.0f  %4.2lf-%4.2lf    %6.1f +- %4.1f +- %4.1f      %1.6f +- %1.6f +- %1.6f  ( %1.6f )\n",
	       DYTools::massBinLimits[i],DYTools::massBinLimits[i+1],
	       rapidityBinLimits[yi], rapidityBinLimits[yi+1],
	       vout[i][yi], voutStatErr[i][yi], voutSystErr[i][yi],
	       voutNorm[i][yi], voutNormStatErr[i][yi], voutNormSystErr[i][yi],
	       sqrt(SQR(voutNormStatErr[i][yi]) + SQR(voutNormSystErr[i][yi])) );
      }
      delete rapidityBinLimits;
    }
  }

  printf("\nPostFsrDET cross-section in the Z peak from %3.0f to %3.0f:\n",
	 DYTools::massBinLimits[low], DYTools::massBinLimits[high+1]);
  printf("           %9.1f +- %8.1f +- %6.1f \n",
	 xsecReference, xsecReferenceStatErr, xsecReferenceSystErr);
  if( use1binErrors ){
    printf("   INFO: the stat and syst errors come from 1-bin 60-120 GeV measurement\n");
    printf("         hardwired in the code. Make sure the numbers are up to date.\n");
  }else{
    printf("   INFO: the stat and syst errors come from the sum of 13 bins in the range 60-120 GeV\n");
  }

  gSystem->mkdir(pathXSect,kTRUE);
  TString extraTag=DYTools::analysisTag;
  if (specTag.Length()) extraTag.Append(TString("_") + specTag);
  TString xSecResultFileName(pathXSect+TString("xSecPostFsrDET_results_") + 
			     extraTag +
			     triggers.triggerConditionsName() + TString(".root"));
  std::cout << "xSecDETResultFileName= " << xSecResultFileName << "\n";

  TMatrixD normXSec=voutNorm;
  TMatrixD normXSecErr=voutNormStatErr; // has to be calculated
  normXSecErr=-1;
  TMatrixD normXSecErrStat=voutNormStatErr;
  TMatrixD normXSecErrSyst=voutNormSystErr;

  TFile fa(xSecResultFileName,"recreate");
  normXSec.Write("normXSec");
  normXSecErr.Write("normXSecErr");
  normXSecErrStat.Write("normXSecErrStat");
  normXSecErrSyst.Write("normXSecErrSyst");
  //normXSecByBin.Write("normXSecByBin");
  //normXSecErrByBin.Write("normXSecErrByBin");
  //normXSecErrByBinStat.Write("normXSecErrByBinStat");
  //normXSecErrByBinSyst.Write("normXSecErrByBinSyst");
  unfolding::writeBinningArrays(fa);
  fa.Close();

  return;
}

// -------------------------------------------------------------------------

void printTableForNotes(const TMatrixD& obs, const TMatrixD& obsErr, 
			const TMatrixD& unf, const TMatrixD& unfErr,
			const TMatrixD& ecor, const TMatrixD& ecorErr,
			const TMatrixD& acor, const TMatrixD& acorErr,
			const TMatrixD& fcor, const TMatrixD& fcorErr)
{

  printf("\n\nLatex table for notes\n");
  printf("               obs-bg                  unfolded                 eff-corrected                acc-corrected              fsr-corrected\n");
  for(int i=0; i<DYTools::nMassBins; i++){
    for (int yi=0; yi<DYTools::nYBins[i]; ++yi) {
      printf("$%4.0f-%4.0f$ &", DYTools::massBinLimits[i],DYTools::massBinLimits[i+1]);
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

// -------------------------------------------------------------------------

void printAllCorrections(){

  TFile fileConstantsEff(fnameEfficiencyConstants);
  TMatrixD *efficiencyArrayPtr    = (TMatrixD *)fileConstantsEff.FindObjectAny("efficiencyArray");
  TMatrixD *efficiencyErrArrayPtr = (TMatrixD *)fileConstantsEff.FindObjectAny("efficiencyErrArray");
  assert(efficiencyArrayPtr); 
  assert(efficiencyErrArrayPtr);
  fileConstantsEff.Close();

  TFile fileScaleConstants(fnameScaleFactorConstants);
  TVectorD *rhoDataMcPtr    = (TVectorD *)fileScaleConstants.FindObjectAny("scaleFactorFlatIdxArray");
  TVectorD *rhoDataMcErrPtr = (TVectorD *)fileScaleConstants.FindObjectAny("scaleFactorErrFlatIdxArray");
  assert(rhoDataMcPtr);
  assert(rhoDataMcErrPtr);
  fileScaleConstants.Close();

  TFile fileConstantsAcc(fnameAcceptanceConstants);
  TMatrixD *acceptanceMatrixPtr    = NULL;
  TMatrixD *acceptanceErrMatrixPtr = NULL;
  if (fileConstantsAcc.IsOpen()) {
    acceptanceMatrixPtr=(TMatrixD *)fileConstantsAcc.FindObjectAny("acceptanceMatrix");
    acceptanceErrMatrixPtr=(TMatrixD *)fileConstantsAcc.FindObjectAny("acceptanceErrMatrix");
    assert(acceptanceMatrixPtr);
    assert(acceptanceErrMatrixPtr);
    fileConstantsAcc.Close();
  }
  
  TFile fileConstantsFsr(fnameFsrCorrectionConstantsBbB);
  TMatrixD *fsrCorrectionArrayPtr    = (TMatrixD *)fileConstantsFsr.FindObjectAny("fsrCorrectionMatrix");
  TMatrixD *fsrCorrectionErrArrayPtr = (TMatrixD *)fileConstantsFsr.FindObjectAny("fsrCorrectionErrMatrix");
  assert(fsrCorrectionArrayPtr);
  assert(fsrCorrectionErrArrayPtr);

  TFile fileConstantsFsrSansAcc(fnameFsrCorrectionSansAccConstantsBbB);
  TMatrixD *fsrCorrectionSansAccArrayPtr    = (TMatrixD *)fileConstantsFsrSansAcc.FindObjectAny("fsrCorrectionMatrix");
  TMatrixD *fsrCorrectionSansAccErrArrayPtr = (TMatrixD *)fileConstantsFsrSansAcc.FindObjectAny("fsrCorrectionErrMatrix");
  assert(fsrCorrectionSansAccArrayPtr);
  assert(fsrCorrectionSansAccErrArrayPtr);

  TMatrixD zero= *efficiencyArrayPtr;
  zero=0;
  TMatrixD efficiencyArray= *efficiencyArrayPtr;
  TMatrixD efficiencyErrArray= *efficiencyErrArrayPtr;
  TVectorD rhoDataMc= *rhoDataMcPtr;
  TVectorD rhoDataMcErr= *rhoDataMcErrPtr;
  TMatrixD acceptanceMatrix= (acceptanceMatrixPtr) ? (*acceptanceMatrixPtr) : zero;
  TMatrixD acceptanceErrMatrix= (acceptanceErrMatrixPtr) ? (*acceptanceErrMatrixPtr) : zero;
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
  for(int i=0, idx=0; i<DYTools::nMassBins; i++){
    double *rapidityBinLimits=DYTools::getYBinLimits(i);
    for (int yi=0; yi<DYTools::nYBins[i]; ++yi, idx++) {
      double effFactor = efficiencyArray[i][yi] * rhoDataMc[idx];
//     double effErr = effFactor
//       * sqrt( efficiencyErrArray[idx]*efficiencyErrArray[idx]/efficiencyArray[idx]/efficiencyArray[idx]
// 	      + rhoDataMcErr[idx]*rhoDataMcErr[idx]/rhoDataMc[idx]/rhoDataMc[idx]);

      double accFactor = acceptanceMatrix[i][yi];
      double accErr    = acceptanceErrMatrix[i][yi];

      double accEff = accFactor * effFactor;
//     double accEffErr = accEff*sqrt( (effErr/effFactor)*(effErr/effFactor) 
// 				    + (accErr/accFactor)*(accErr/accFactor));

      double fsrFactor = fsrCorrectionArray[i][yi];
      double fsrErr    = fsrCorrectionErrArray[i][yi];

      double fsrFactorSansAcc = fsrCorrectionSansAccArray[i][yi];
      double fsrErrSansAcc    = fsrCorrectionSansAccErrArray[i][yi];

      sprintf(buf,"$ %4.0f-%4.0f $&", DYTools::massBinLimits[i],DYTools::massBinLimits[i+1]);
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

// -------------------------------------------------------------

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
  for(int i=0; i<DYTools::nMassBins; i++){
    double *rapidityBinLimits=DYTools::getYBinLimits(i);
    for (int yi=0; yi<DYTools::nYBins[i]; ++yi) {
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
      sprintf(buf,"$%4.0f-%4.0f $&$ %4.2f-%4.2f $&", DYTools::massBinLimits[i],DYTools::massBinLimits[i+1],rapidityBinLimits[yi],rapidityBinLimits[yi+1]);
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

void printRelativeSystErrorsPAS(){

  // The modeling error is just as for muons (as an approximation)
  double systModelingRelative[40] = {
    9.700000e-02 ,  3.100000e-02 ,  1.900000e-02 ,  7.000000e-03 ,  5.000000e-03 ,
    3.000000e-03 ,  1.000000e-03 ,  1.000000e-03 ,  2.000000e-03 ,  2.000000e-03 ,
    3.000000e-03 ,  2.000000e-03 ,  2.000000e-03 ,  2.000000e-03 ,  1.000000e-03 ,
    2.000000e-03 ,  2.000000e-03 ,  3.000000e-03 ,  3.000000e-03 ,  5.000000e-03 ,
    6.000000e-03 ,  6.000000e-03 ,  5.000000e-03 ,  6.000000e-03 ,  7.000000e-03 ,
    1.000000e-02 ,  1.100000e-02 ,  1.100000e-02 ,  1.100000e-02 ,  1.100000e-02 ,
    1.000000e-02 ,  1.500000e-02 ,  1.400000e-02 ,  1.300000e-02 ,  1.500000e-02 ,
    1.200000e-02 ,  1.300000e-02 ,  2.100000e-02 ,  2.400000e-02 ,  3.100000e-02 
  };
  
  std::string fileName="tbl_relativeSystErrorsPAS.tex";
  std::ofstream fout;
  fout.open(fileName.c_str());
  fout << "\n\nLatex table of relative systematic errors  in percent for PAS/paper\n";
  fout << "\\begin{table}[tbh]\n";
  fout << "\\begin{center}\\begin{tabular}{ccccccccc}\n";
  fout << "\\hline\n";
  fout << " Mass range & rapidity range   &   Escale  &   Eff.   &   Bkg    &    Unfol &    sum   & Acc(PDF+stat)%%  & modeling\n";   
  fout << "\\\\ \\hline\n";

  char buf[200];
  for(int i=0; i<DYTools::nMassBins; i++){
    double *rapidityBinLimits=DYTools::getYBinLimits(i);
    for (int yi=0; yi<DYTools::nYBins[i]; ++yi) {
      // Factor out theory error from the total acceptance error
      double systAcceptanceExpRelative = systAcceptanceRelative[i][yi];
      // The "sum" contains only the experimental systematic errors
      double sum = sqrt(systEscaleRelative[i][yi]*systEscaleRelative[i][yi]
			+ systEfficiency[i][yi]*systEfficiency[i][yi]
			+ systBackgrRelative[i][yi]*systBackgrRelative[i][yi]
			+ systUnfoldRelative[i][yi]*systUnfoldRelative[i][yi]);
      sprintf(buf,"$%4.0f-%4.0f $&$ %4.2f-%4.2f $&", DYTools::massBinLimits[i],DYTools::massBinLimits[i+1],rapidityBinLimits[yi],rapidityBinLimits[yi+1]);
      fout << buf;
      sprintf(buf,"$  %5.1f  $&$ %5.1f $&$ %5.1f $&$ %5.1f $&$ %5.1f $&$ %5.1f $&$ %5.1f $ \\\\" ,
	      100*systEscaleRelative[i][yi], 
	      100*systEfficiency[i][yi], 
	      100*systBackgrRelative[i][yi], 
	      100*systUnfoldRelative[i][yi],
	      100*sum,
	      100*systAcceptanceExpRelative,
	      100*systModelingRelative[i]
	      );
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

// -------------------------------------------------------------

void printRelativeSystErrorsForCovarianceMatrix(){
  if (DYTools::study2D==1) {
    std::cout << "\n\n\tprintRelativeSystErrorsForCovarianceMatrix is not for 2D\n";
    return;
  }

  TFile fileConstants(fnameFsrCorrectionConstantsBbB);
  assert(fileConstants.IsOpen());
  TMatrixD *fsrCorrectionMatrixPtr  = (TMatrixD *)fileConstants.FindObjectAny("fsrCorrectionMatrix");
  TMatrixD *fsrCorrectionErrMatrixPtr = (TMatrixD *)fileConstants.FindObjectAny("fsrCorrectionErrMatrix");
  assert(fsrCorrectionMatrixPtr && fsrCorrectionErrMatrixPtr);
  TMatrixD fsrCorrectionMatrix= *fsrCorrectionMatrixPtr;
  TMatrixD fsrCorrectionErrMatrix= *fsrCorrectionErrMatrixPtr;
  fileConstants.Close();


  if (1) {
  std::string fileName=Form("tbl_relativeSystErrorsForCovariance%s.dat",DYTools::analysisTag.Data());
  std::ofstream fout;
  fout.open(fileName.c_str());
  fout << "#preFSR-full acceptance\n";
  fout << "#\n#Table of relative statistical and systematic errors  in percent for covariance matrix\n#\n";
  fout << "# Statistical errors is on yields. No systematics included\n";
  fout << "# Sum = sqrt(Stat^2 + Escale^2 + Eff^2 + Bkgr^2 + Unfol^2 + fsr^2 )\n";
  fout << "# Mass idx    stat.err     Escale           Eff.            Unfol             Bkg            FSR             sum\n"; 
  fout << "#\n";

  char buf[200];
  for(int i=0; i<DYTools::nMassBins; i++){
    //double *rapidityBinLimits=DYTools::getYBinLimits(i);
    for (int yi=0; yi<DYTools::nYBins[i]; ++yi) {
      double fsrErrRelative=fsrCorrectionErrMatrix[i][yi]/fsrCorrectionMatrix[i][yi];
      // The "sum" contains only the experimental systematic errors
      double sum = sqrt(systEscaleRelative[i][yi]*systEscaleRelative[i][yi]
			+ systEfficiency[i][yi]*systEfficiency[i][yi]
			+ systBackgrRelative[i][yi]*systBackgrRelative[i][yi]
			+ systUnfoldRelative[i][yi]*systUnfoldRelative[i][yi]
			+ fsrErrRelative*fsrErrRelative
			);
      sprintf(buf,"%d   %14.11f   %14.11f   %14.11f   %14.11f   %14.11f   %14.11f   %14.11f",
	      i,
	      100*unfoldedYieldsStatErr[i][yi]/unfoldedYields[i][yi],
	      100*systEscaleRelative[i][yi],
	      100*systEfficiency[i][yi],
	      100*systUnfoldRelative[i][yi],
	      100*systBackgrRelative[i][yi],
	      100*fsrErrRelative,
	      100*sum
	      );
      fout << buf << "\n";
    }
    //delete rapidityBinLimits;
  }
  fout.close();
  std::cout << "file <" << fileName << "> created\n";
  }

  if (1) {
  std::string fileName=Form("tbl_relativeSystErrorsForCovariance%s_accOnly.dat",DYTools::analysisTag.Data());
  std::ofstream fout;
  fout.open(fileName.c_str());
  fout << "#preFSR-full acceptance\n";
  fout << "#\n#Table of systematic errors on acceptance.\n";
  fout << "#\n";
  fout << " Mass idx    accExperiment   accTheory   sum\n";

  char buf[200];
  for(int i=0; i<DYTools::nMassBins; i++){
    //double *rapidityBinLimits=DYTools::getYBinLimits(i);
    for (int yi=0; yi<DYTools::nYBins[i]; ++yi) {
      // Factor out theory error from the total acceptance error
      double systAcceptanceExpRelative 
	= sqrt( systAcceptanceRelative[i][yi]*systAcceptanceRelative[i][yi]
		- systAccTheoryRelative[i][yi]*systAccTheoryRelative[i][yi]);
      sprintf(buf,"%d   %14.11f   %14.11f   %14.11f",
	      i,
	      100*systAcceptanceExpRelative,
	      100*systAccTheoryRelative[i][yi],
	      100*systAcceptanceRelative[i][yi]
	      );
      fout << buf << "\n";
    }
    //delete rapidityBinLimits;
  }
  fout.close();
  std::cout << "file <" << fileName << "> created\n";
  }

  std::cout << "HERE" << std::endl;
  return;
}

// --------------------------------------------------------------
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

// --------------------------------------------------------------

void plotAccEff() {
  TFile fileConstantsAcc(fnameAcceptanceConstants);
  if (!fileConstantsAcc.IsOpen()) {
    std::cout << "failed to open acceptance systematics file <" << fnameAcceptanceConstants << ">\n";
    assert(0);
  }
  TMatrixD *acceptanceMatrixPtr    = (TMatrixD *)fileConstantsAcc.FindObjectAny("acceptanceMatrix");
  TMatrixD *acceptanceErrMatrixPtr = (TMatrixD *)fileConstantsAcc.FindObjectAny("acceptanceErrMatrix");
  if (!acceptanceMatrixPtr || !acceptanceErrMatrixPtr) {
    std::cout << "at least one object from file <" << fnameAcceptanceConstants << "> is null\n";
    assert(0);
  }
  fileConstantsAcc.Close();

  TMatrixD acceptanceMatrix = *acceptanceMatrixPtr;
  TMatrixD acceptanceErrMatrix = *acceptanceErrMatrixPtr;

  TFile fileConstants(fnameEfficiencyConstants);
  TMatrixD* efficiencyMatrixPtr    = (TMatrixD *)fileConstants.FindObjectAny("efficiencyArray");
  TMatrixD* efficiencyErrMatrixPtr = (TMatrixD *)fileConstants.FindObjectAny("efficiencyErrArray");
  fileConstants.Close();

  if (!efficiencyMatrixPtr || !efficiencyErrMatrixPtr) {
    std::cout << "at least one needed object is not present in <" << fnameEfficiencyConstants << ">\n";
    assert(0);
  }
  TMatrixD efficiencyMatrix= *efficiencyMatrixPtr;
  TMatrixD efficiencyErrMatrix= *efficiencyErrMatrixPtr;

  const int iYBin=0;
  const int perMassBinWidth=0;
  const int perRapidityBinWidth=0;
 
  TCanvas *canv=MakeCanvas("canvAccEff","canvAccEff",800,800);

  TH1F* hAcc=extractMassDependence("hAcc","",
				   acceptanceMatrix, acceptanceErrMatrix,
				   iYBin,
				   perMassBinWidth,perRapidityBinWidth);
     
  TH1F* hEff=extractMassDependence("hEff","",
				   efficiencyMatrix,efficiencyErrMatrix,
				   iYBin,
				   perMassBinWidth,perRapidityBinWidth);
     
  InitHist(hAcc,"M_{ee} [GeV]","",kBlack);
  TH1F *hAccEff=(TH1F*) hAcc->Clone("hAccEff");
  hAccEff->Multiply(hEff);
  CPlot cp("cpAccEff","","M_{ee} [GeV]",""); 
  cp.SetLogx();
  cp.SetYRange(0.,1.2);
  cp.AddTextBox("#gamma*/Z#rightarrow ee",0.25,0.75,0.5,0.91, 0,kBlack,kWhite);
  cp.AddTextBox("at #sqrt{s}=7 TeV, CMS Simulation", 0.55,0.75,0.92,0.91, 0,kBlack,kWhite);
  hAcc->GetXaxis()->SetMoreLogLabels();
  hAcc->GetXaxis()->SetNoExponent();
  hAcc->SetMarkerStyle(20);
  hEff->SetMarkerStyle(26);
  hAccEff->SetMarkerStyle(25);
  //printHisto(hAcc);
  hAcc->SetDirectory(0);
  hEff->SetDirectory(0);
  hAccEff->SetDirectory(0);
  cp.AddHist1D(hAcc,"A","LP",kBlue,1,1,1);
  cp.AddHist1D(hEff,"#varepsilon","LP",kRed,1,0,1);
  cp.AddHist1D(hAccEff,"A#times#varepsilon","LP",kBlack,1,0,1);
  cp.Draw(canv,false,"png",0);
  cp.TransLegend(0.1,-0.5);
  canv->Update();
  SaveCanvas(canv,canv->GetName());
}


// --------------------------------------------------------------


void initGlobalFileNames(const TriggerSelection &triggers, const TString &tagDirYields, const TString &tagDirConstants, const TString &tagDirScaleFactorConstants, const TString &tagDirXSect, TFsrCorrectionType_t fsrBinByBin) {
  const TString fnameEnd(DYTools::analysisTag + TString(".root"));
  const TString fnameEndWithPU(DYTools::analysisTag + TString("_PU") + TString(".root"));
  pathYields=TString("../root_files/yields/") + tagDirYields +
    TString("/");
  pathSystematics=TString("../root_files/systematics/") + tagDirConstants +
    TString("/");
  pathConstants= TString("../root_files/constants/") + tagDirConstants +
    TString("/");
  pathScaleFactors= TString("../root_files/constants/") + tagDirScaleFactorConstants + TString("/");
  pathXSect= TString("../root_files/") + tagDirXSect + TString("/");


  // Data yields
  fnameDataYields= pathYields + TString("yields_bg-subtracted") + fnameEnd;

  // Vectors for MC closure test
  if (fsrCorrection_BinByBin==_fsrCorr_binByBin) {
    fnameMcReferenceYields= 
      pathConstants + TString("yields_MC_unfolding_reference_");
  }
  else {
    fnameMcReferenceYields= pathConstants + TString("yields_detResponse");
    if (useExactVectorsForMcClosureTest) fnameMcReferenceYields.Append("Exact");
  }
  fnameMcReferenceYields.Append( fnameEndWithPU );
  
  fnameMcReferenceYieldsFsr=pathConstants + TString("yields_fsr");
  if (useExactVectorsForMcClosureTest) fnameMcReferenceYieldsFsr.Append("Exact");
  else fnameMcReferenceYieldsFsr.Append("Good");
  fnameMcReferenceYieldsFsr.Append(fnameEndWithPU);

  fnameMcReferenceYieldsFsrDET=pathConstants + TString("yields_fsrDET");
  if (useExactVectorsForMcClosureTest) fnameMcReferenceYieldsFsrDET.Append("exact");
  else fnameMcReferenceYieldsFsrDET.Append("_good");
  fnameMcReferenceYieldsFsrDET.Append(fnameEndWithPU);


// This file contains unfolding matrix 
  TString unfFileStart=pathConstants;
  if (fsrBinByBin!=_fsrCorr_binByBin) {
    if (useExactVectorsForMcClosureTest) unfFileStart.Append("detResponseExact_");
    else unfFileStart.Append("detResponse_");
  }
  fnameUnfoldingConstants=unfFileStart + TString("unfolding_constants") + fnameEndWithPU;

// Contains relative unfolding systematic errors
  fnameUnfoldingSystErrors=pathSystematics +
    TString("unfolding_systematics") + fnameEnd;

// Contains relative escale systematic errors
  fnameEscaleSystErrors= pathSystematics + 
    TString("escale_systematics") + fnameEnd;

  // efficiency factors
  fnameEfficiencyConstants = pathConstants +
    TString("event_efficiency_constants") + fnameEnd;

  fnameScaleFactorConstants = pathScaleFactors +
    TString("scale_factors_") + DYTools::analysisTag +
    TString("_") + triggers.triggerConditionsName() +
    TString("_PU") +
    TString(".root");

  // acceptance factors
  fnameAcceptanceConstants = pathConstants +
    TString("acceptance_constants") + fnameEnd;
  fnameAcceptanceSystematics = pathSystematics +
    TString("theoretical_uncertainties.root");
  fnameAcceptanceFSRSystematics = pathSystematics + 
    TString("acceptance_FSR_systematics") + fnameEnd;

  // FSR correction
  fnameFsrCorrectionConstantsBbB = pathConstants + 
    TString("fsr_constants_") + fnameEnd;
  fnameFsrCorrectionSansAccConstantsBbB = pathConstants +
    TString("fsr_constants_") + DYTools::analysisTag + TString("_sans_acc.root");
  switch(fsrCorrection_BinByBin) {
  case _fsrCorr_binByBin: break;
  case _fsrCorr_unfPure:
  case _fsrCorr_unf:
    fnameFsrCorrectionConstantsUnf = pathConstants +
      TString("fsrExact_unfolding_constants") + fnameEnd;
    fnameFsrCorrectionDETConstantsUnf = pathConstants +
      TString("fsrDET_unfolding_constants") + fnameEnd;
    break;
  case _fsrCorr_unfMdf:
    fnameFsrCorrectionConstantsUnf = pathConstants +
      TString("fsrGood_unfolding_constants") + fnameEnd;
    fnameFsrCorrectionDETConstantsUnf = pathConstants +
      TString("fsrDET_Mdf_unfolding_constants") + fnameEnd;
    break;
  case _fsrCorr_unfGood:
    fnameFsrCorrectionConstantsUnf = pathConstants +
      TString("fsrGood_unfolding_constants") + fnameEndWithPU;
    fnameFsrCorrectionDETConstantsUnf = pathConstants +
      TString("fsrDET_good_unfolding_constants") + fnameEndWithPU;
    break;
  default:
    std::cout << "code not ready\n";
    assert(0);
  }
  fnameFsrDETcorrFactors = pathConstants +
    TString("yields_fsrCorrFactors") + fnameEnd;
}



// --------------------------------------------------------------
