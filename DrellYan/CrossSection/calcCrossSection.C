#include "TVectorD.h"
#include "TFile.h"
#include "TString.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include "../Include/DYTools.hh"

// This global constants will be filled from 
// the configuration file. This is not great C++...
TString tagDirYields = "";
TString tagDirConstants = "";
Double_t lumi = 0;

const TString fileDataYields        ("yields_bg-subtracted.root");

const TString fileMcReferenceYields ("yields_MC_unfolding_reference.root");

const TString fileUnfoldingConstants("unfolding_constants.root");

const TString fileEfficiencyConstants("event_efficiency_constants.root");

const TString fileScaleFactorConstants("scale_factors.root");

const TString fileAcceptanceConstants("acceptance_constants.root");

const TString fileFsrCorrectionConstants("fsr_constants.root");

const TString fileFsrCorrectionSansAccConstants("fsr_constants_sans_acc.root");

// Forward declarations
void readData(TVectorD &v, TVectorD &vErr1, TVectorD &vErr2);

void  unfold(TVectorD &vin, TVectorD &vinStatErr, TVectorD &vinSystErr,
             TVectorD &vout, TVectorD &voutStatErr, TVectorD &voutSystErr);
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

// The arrays below contain the estimate of relative systematic
// error in percent obtained elsewhere. The calculation of these
// errors is done outside of the main scripts that calculate
// the unfolding constants and energy scale corrections. We want
// to change this in the future from hardwired text constants here
// to reading ROOT files.
const double unfoldingSystematicsPercent[nMassBins] = 
  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

const double escaleSystematicsPercent[nMassBins] = 
  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

// Some global arrays for convenience.
// These will contain errors on R
// (after unfolding, etc).
// This is messy, until the next re-write
TVectorD systBackgrBeforeUnfolding(DYTools::nMassBins);
TVectorD systEscaleBeforeUnfolding(DYTools::nMassBins);
TVectorD systBackgr(DYTools::nMassBins);
TVectorD systEscale(DYTools::nMassBins);
TVectorD systEfficiency(DYTools::nMassBins);
TVectorD systUnfolding(DYTools::nMassBins);
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
  systBackgr = 0;
  systEscale = 0;
  systEfficiency = 0;
  systUnfolding = 0;
  systOthers = 0;

  // Read data yields from file (background subtraction is already done)
  readData(signalYields, signalYieldsStatErr, signalYieldsSystErr);

  // Apply unfolding
  unfold(signalYields, signalYieldsStatErr, signalYieldsSystErr,
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
  return;
}

//-----------------------------------------------------------------
// Unfold
//-----------------------------------------------------------------
void  unfold(TVectorD &vin, TVectorD &vinStatErr, TVectorD &vinSystErr,
             TVectorD &vout, TVectorD &voutStatErr, TVectorD &voutSystErr)
{

  // Read unfolding constants
  printf("unfold: Load constants\n"); fflush(stdout);
    
  TFile fileConstants(TString("../root_files/constants/")+tagDirConstants+TString("/")+fileUnfoldingConstants);
  TMatrixD DetResponse             = *(TMatrixD *)fileConstants.FindObjectAny("DetResponse");
  TMatrixD DetInvertedResponse     = *(TMatrixD *)fileConstants.FindObjectAny("DetInvertedResponse");
  TMatrixD DetInvertedResponseErr  = *(TMatrixD *)fileConstants.FindObjectAny("DetInvertedResponseErr");
  TVectorD BinLimitsArray          = *(TVectorD *)fileConstants.FindObjectAny("BinLimitsArray");

  printf("Load MC reference yields\n"); fflush(stdout);
  TFile fileMcRef(TString("../root_files/constants/")+tagDirConstants+TString("/")+fileMcReferenceYields);
  TVectorD yieldsMcFsrOfRec        = *(TVectorD *)fileMcRef.FindObjectAny("yieldsMcFsrOfRec");
  TVectorD yieldsMcRec             = *(TVectorD *)fileMcRef.FindObjectAny("yieldsMcRec");

 // Check that the binning is consistent
  bool checkResult = true;
  int nBins = BinLimitsArray.GetNoElements()-1;
  if( DetResponse.GetNrows()    != nMassBins ) checkResult = false;
  if( BinLimitsArray.GetNoElements() != nMassBins+1) checkResult = false;
  for(int i=0; i<nBins+1; i++){
    if( BinLimitsArray[i] != massBinLimits[i] )
      checkResult = false;
  }
  if( !checkResult ){
    printf("unfold: ERROR: inconsistent binning in the inputs\n");
    assert(0);
  }else
    printf("unfold: Binning in the inputs is consistent\n");

  // Apply unfolding matrix
  TVectorD systErrorPreviousPropagated(nBins);
  TVectorD systErrorAdditional(nBins);
  TVectorD dNdMmcCheck(nBins);
  // Initialize all elements. Stat error is propagated, 
  // but systematic error has the propagated component and the new additions
  vout = 0;
  voutStatErr = 0;
  voutSystErr = 0;
  systErrorPreviousPropagated = 0;
  systErrorAdditional = 0;
  dNdMmcCheck = 0;
  systBackgr = 0;
  systEscale = 0;
  for(int i=0; i<nBins; i++){
    for(int j=0; j<nBins; j++){
      vout[i] += DetInvertedResponse(j,i) * vin[j];
      voutStatErr                [i] += pow( DetInvertedResponse   (j,i) * vinStatErr[j], 2);

      systErrorPreviousPropagated[i] += pow( DetInvertedResponse   (j,i) * vinSystErr[j], 2);
      systErrorAdditional        [i] += pow( DetInvertedResponseErr(j,i) * vin[j], 2);

      systBackgr[i] += pow( DetInvertedResponse   (j,i) * systBackgrBeforeUnfolding[j], 2);
      systEscale[i] += pow( DetInvertedResponse   (j,i) * systEscaleBeforeUnfolding[j], 2);

      dNdMmcCheck[i] += DetInvertedResponse(j,i) * yieldsMcRec[j];
    }
    voutStatErr[i] = sqrt(voutStatErr[i]);
    systBackgr[i] = sqrt(systBackgr[i])/vout[i];
    systEscale[i] = sqrt(systEscale[i])/vout[i];

    systErrorPreviousPropagated[i] = sqrt(systErrorPreviousPropagated[i]);
    // the "additional" systematic error below comes from error propagation
    // of the errors on the unfolding matrix elements
    systErrorAdditional[i]         = sqrt(systErrorAdditional[i]);
    // Add externally computed systematics from other sources for this step
    double unfoldingSystematics = (unfoldingSystematicsPercent[i]/100.0)*vout[i];
    systErrorAdditional[i] = sqrt( systErrorAdditional[i]*systErrorAdditional[i] + 
				   unfoldingSystematics*unfoldingSystematics);
    systUnfolding[i] = systErrorAdditional[i]/vout[i];

    // The total systematic error ends up being from 3 sources:
    //  - original systematic error on the signal yields before unfolding, propagated through unfolding
    //  - the systematic error due to uncertainty in the elements of the unfolding matrix due
    //         to limites statistics
    //  - the systematic error due to uncertainty in extra smearing applied to MC
    //         when the unfolding matrix is extracted
    voutSystErr[i] = sqrt(systErrorPreviousPropagated[i]*systErrorPreviousPropagated[i] 
			  + systErrorAdditional[i]*systErrorAdditional[i]);
  }
  
  // Report results
  printf("\nUNFOLD: Check on the MC, yields:\n");
  for(int i=0; i<nBins; i++){
    printf("%4.0f-%4.0f   %10.0f    %10.0f\n",
	   BinLimitsArray[i],BinLimitsArray[i+1],
	   yieldsMcFsrOfRec[i],
	   dNdMmcCheck[i]);
  }
  printf("\nUNFOLD: Results for the data, yields:\n");
  printf("                   yields observed        after unfolding                     syst-err-unfolding,%%\n");
  for(int i=0; i<nBins; i++){
    printf("%4.0f-%4.0f   %8.1f +- %6.1f +- %5.1f       %8.1f +- %7.1f +- %6.1f       %7.1f\n",
	   BinLimitsArray[i],BinLimitsArray[i+1],
	   vin[i], vinStatErr[i], vinSystErr[i],
	   vout[i], voutStatErr[i], voutSystErr[i],
	   systErrorAdditional[i]*100.0/vout[i]);
  }

  fileMcRef.Close();
  fileConstants.Close();
  return;
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
  int low = 9;
  int high = 22;
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
    printf("%4.0f-%4.0f   %9.1f +- %8.1f +- %6.1f   %9.4f +- %8.4f +- %8.4f  ( %8.4f )\n",
	   massBinLimits[i],massBinLimits[i+1],
	   vout[i], voutStatErr[i], voutSystErr[i],
	   voutNorm[i], voutNormStatErr[i], voutNormSystErr[i],
	   sqrt(voutNormStatErr[i] * voutNormStatErr[i] + voutNormSystErr[i] * voutNormSystErr[i]) );
  }
  printf("\nPre FSR cross-section in the Z peak from %3f to %3f:\n",
	 massBinLimits[low], massBinLimits[high]);
  printf("           %9.1f +- %8.1f +- %6.1f \n",
	 xsecReference, xsecReferenceStatErr, xsecReferenceSystErr);
    printf("check %f %f\n", xsecReferenceStatErr, xsecReferenceSystErr);

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
  int low = 9;
  int high = 22;
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
    printf("%4.0f-%4.0f   %9.1f +- %8.1f +- %6.1f   %9.4f +- %8.4f +- %8.4f  ( %8.4f )\n",
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
  int low = 9;
  int high = 22;
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
    printf("%4.0f-%4.0f   %9.1f +- %8.1f +- %6.1f   %9.4f +- %8.4f +- %8.4f  ( %8.4f )\n",
	   massBinLimits[i],massBinLimits[i+1],
	   vout[i], voutStatErr[i], voutSystErr[i],
	   voutNorm[i], voutNormStatErr[i], voutNormSystErr[i],
	   sqrt(voutNormStatErr[i] * voutNormStatErr[i] + voutNormSystErr[i] * voutNormSystErr[i]) );
  }
  printf("\nPostFsr cross-section in the Z peak from %3f to %3f:\n",
	 massBinLimits[low], massBinLimits[high]);
  printf("           %9.1f +- %8.1f +- %6.1f \n",
	 xsecReference, xsecReferenceStatErr, xsecReferenceSystErr);
    printf("check %f %f\n", xsecReferenceStatErr, xsecReferenceSystErr);

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
  int low = 9;
  int high = 22;
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
    printf("%4.0f-%4.0f   %9.1f +- %8.1f +- %8.1f   %9.4f +- %8.4f +- %8.4f  ( %8.4f )\n",
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
    double effErr = effFactor
      * sqrt( efficiencyErrArray[i]*efficiencyErrArray[i]/efficiencyArray[i]/efficiencyArray[i]
	      + rhoDataMcErr[i]*rhoDataMcErr[i]/rhoDataMc[i]/rhoDataMc[i]);

    double accFactor = acceptanceArray[i];
    double accErr    = acceptanceErrArray[i];

    double accEff = accFactor * effFactor;
    double accEffErr = accEff*sqrt( (effErr/effFactor)*(effErr/effFactor) 
				    + (accErr/accFactor)*(accErr/accFactor));

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
    double sum = sqrt(systEscale[i]*systEscale[i]
		      + systEfficiency[i]*systEfficiency[i]
		      + systBackgr[i]*systBackgr[i]
		      + systUnfolding[i]*systUnfolding[i]);
    printf("%4.0f-%4.0f &", massBinLimits[i],massBinLimits[i+1]);
    printf("   $%5.1f$ &  $%5.1f$ &  $%5.1f$ &  $%5.1f$ &  $%5.1f$", 
	   100*systEscale[i], 
	   100*systEfficiency[i], 
	   100*systBackgr[i], 
	   100*systUnfolding[i],
	   100*sum);
    printf("\n");
  }
}

