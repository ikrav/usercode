#include "TVectorD.h"
#include "TFile.h"
#include "TString.h"
#include "TSystem.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "../Include/DYTools.hh"
#include "../Include/UnfoldingTools.hh"

using std::string;
using std::stringstream;

// The global variables will be used in several functions:
TString tagDirYields = "";
TString tagDirConstants = "";
Double_t lumi = 0;

// load in flat index format
int readData(const TString &fname, TVectorD &vFI, TVectorD &vErr1FI, 
	     TVectorD &vErr2FI, int debug=0);

void  applyUnfoldingLocal(TVectorD &vinFI, TVectorD &voutFI, bool ifSeed, int seed, int reweightInt);

const TString fileEnd( DYTools::analysisTag + TString(".root") );
const TString fileDataYields        (TString("yields_bg-subtracted") + fileEnd);

const TString fileUnfoldingConstantsBaseSeed(
     TString("unfolding_constants_seed_") + DYTools::analysisTag + TString("_") );
const TString fileUnfoldingConstantsBaseReweight(
     TString("unfolding_constants_reweight_") + DYTools::analysisTag + TString("_") );

const int seedFirst = 1001;
const int seedLast = 1020;


//
//  Main code
//

void calcUnfoldingSystematics(const TString conf){

  // check whether it is a calculation
  if (conf.Contains("_DebugRun_")) {
    std::cout << "calcUnfoldingSystematics: _DebugRun_ detected. Terminating the script\n";
    return;
  }

  // First, read the configuration file. The configuration file
  // is the same as the one used for the cross section calculation
  // script, we need to know the location of data yields and
  // unfolding matrices

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

  int nUnfoldingBins = DYTools::getTotalNumberOfBins();

  TVectorD signalYields(nUnfoldingBins);
  TVectorD signalYieldsStatErr(nUnfoldingBins);
  TVectorD signalYieldsSystErr(nUnfoldingBins);

  TVectorD unfoldedYields(nUnfoldingBins);

  TVectorD unfoldedYieldsMean(nUnfoldingBins);
  TVectorD unfoldedYieldsRMS(nUnfoldingBins);
  TVectorD unfoldedYieldsSquaredMean(nUnfoldingBins);

  unfoldedYieldsMean = 0;
  unfoldedYieldsRMS = 0;
  unfoldedYieldsSquaredMean = 0;

  // Read data yields from file
  TString dataYieldsFileName=TString("../root_files/yields/") +
    tagDirYields + TString("/") + fileDataYields;
  int debugRead=0;
  if (readData(dataYieldsFileName,signalYields, 
	       signalYieldsStatErr, signalYieldsSystErr, debugRead) != 1) {
    std::cout << "failed to load yields from a file <" 
	      << dataYieldsFileName << ">\n";
    throw 2;
  }


/////////////////////////////////
//calculate smearing systematics 
/////////////////////////////////

  int nseeds = 0;

  for(int i=seedFirst; i<=seedLast; i++){
    nseeds++;
    applyUnfoldingLocal(signalYields, unfoldedYields, 1, i, 100);
    for(int idx = 0; idx < nUnfoldingBins; idx++){
      unfoldedYieldsMean[idx] += unfoldedYields[idx];
      unfoldedYieldsSquaredMean[idx] += unfoldedYields[idx]*unfoldedYields[idx];
    }
  }


  // Final calculation of the mean and RMS for Smearing
  TVectorD unfoldingSystPercentSmear(nUnfoldingBins);
  for(int idx = 0; idx < nUnfoldingBins; idx++){
    unfoldedYieldsMean[idx] = unfoldedYieldsMean[idx]/double(nseeds);
    unfoldedYieldsSquaredMean[idx] = 
      unfoldedYieldsSquaredMean[idx]/double(nseeds);
    unfoldedYieldsRMS[idx] = 
      sqrt(unfoldedYieldsSquaredMean[idx] - 
	   unfoldedYieldsMean[idx]*unfoldedYieldsMean[idx]);
    unfoldingSystPercentSmear[idx] = 
      unfoldedYieldsRMS[idx]*100.0/unfoldedYieldsMean[idx];
  }
  
  
/////////////////////////////////
//calculate Fsr systematics 
/////////////////////////////////

  TVectorD unfoldedYieldsFsrMax(nUnfoldingBins);
  TVectorD unfoldedYieldsFsrMin(nUnfoldingBins);
  TVectorD unfoldedYieldsFsrErr(nUnfoldingBins);
  applyUnfoldingLocal(signalYields, unfoldedYieldsFsrMax, 0, 1000, 105);
  applyUnfoldingLocal(signalYields, unfoldedYieldsFsrMin, 0, 1000, 95);

  TVectorD unfoldingSystPercentFsr(nUnfoldingBins); 

  for(int idx = 0; idx < nUnfoldingBins; idx++) {
    unfoldedYieldsFsrErr[idx] = 
      fabs(unfoldedYieldsFsrMax[idx]-unfoldedYieldsFsrMin[idx]) /
      (unfoldedYieldsFsrMax[idx]+unfoldedYieldsFsrMin[idx]);
    unfoldingSystPercentFsr[idx] = unfoldedYieldsFsrErr[idx]*100.0;
  }
 

/////////////////////////////////
//calculate total systematics 
/////////////////////////////////
    TVectorD unfoldingSystPercent(nUnfoldingBins); 
    for(int idx = 0; idx < nUnfoldingBins; idx++) {
      unfoldingSystPercent[idx] = 
	sqrt(
	     unfoldingSystPercentSmear[idx]*unfoldingSystPercentSmear[idx] +
	     unfoldingSystPercentFsr[idx]*unfoldingSystPercentFsr[idx] 
	     );
    }

//printing out to the screen

   printf("mass     rapidity   mean-unfolded   RMS-unfolded   rel-error     rel-err-percent-Smear     rel-err-percent-Fsr      rel-err-percent-total \n");
   for(int i=0, idx=0; i<DYTools::nMassBins; i++){
    double *rapidityBinLimits=DYTools::getYBinLimits(i);
    for (int yi=0; yi<DYTools::nYBins[i]; ++yi, ++idx) {
      printf("%4.0f-%4.0f  %4.2lf-%4.2lf  %7.1f      %7.1f          %6.4f             %6.1f                 %6.3f                       %6.1f\n",
	     DYTools::massBinLimits[i],DYTools::massBinLimits[i+1],
	     rapidityBinLimits[yi],rapidityBinLimits[yi+1],
	     unfoldedYieldsMean[idx], unfoldedYieldsRMS[idx],
	     unfoldedYieldsRMS[idx]/unfoldedYieldsMean[idx],
	     unfoldingSystPercentSmear[idx], unfoldingSystPercentFsr[idx], 
	     unfoldingSystPercent[idx]);
      //unfoldedYieldsRMS[idx]*100.0/unfoldedYieldsMean[idx]);
    }
  }
  


  // Store constants in the file
  TString outputDir(TString("../root_files/systematics/")+tagDirConstants);
  gSystem->mkdir(outputDir,kTRUE);
  TString unfoldingSystFileName(outputDir+TString("/unfolding_systematics") 
				+ DYTools::analysisTag + TString(".root"));

  TFile fa(unfoldingSystFileName,"recreate");
  unfolding::writeBinningArrays(fa);
  unfoldedYieldsMean.Write("unfoldedYieldsMeanFI");
  unfoldedYieldsRMS.Write("unfoldedYieldsRMSFI");
  unfoldingSystPercentSmear.Write("unfoldingSystPercentSmearFI");
  unfoldingSystPercentFsr.Write("unfoldingSystPercentFsrFI");
  unfoldingSystPercent.Write("unfoldingSystPercentFI");
  fa.Close();

  return;
}

//-----------------------------------------------------------------
// Read data
//-----------------------------------------------------------------
int readData(const TString &fname, TVectorD &vFI, TVectorD &vErrFI1, TVectorD &vErrFI2, int debug){

  if (debug) std::cout << "Load data yields from <" << fname << ">" << std::endl;
  TFile fileYields   (fname);
  if (!fileYields.IsOpen()) {
    std::cout << "failed to open a file <" << fname << ">\n";
    return 0;
  }
  int res=unfolding::checkBinningArrays(fileYields);
  if (res!=1) return res;
  TMatrixD *YieldsSignalPtr = (TMatrixD *)fileYields.FindObjectAny("YieldsSignal");
  TMatrixD *YieldsSignalErrPtr    = (TMatrixD *)fileYields.FindObjectAny("YieldsSignalErr");
  TMatrixD *YieldsSignalSystErrPtr= (TMatrixD *)fileYields.FindObjectAny("YieldsSignalSystErr");
  fileYields.Close();

  if (!YieldsSignalPtr || !YieldsSignalErrPtr || !YieldsSignalSystErrPtr) {
    std::cout << "failed to get at least one of the required objects from file <" 
	      << fileDataYields << ">\n";
    assert(0);
  }

  TMatrixD YieldsSignal = *YieldsSignalPtr;
  TMatrixD YieldsSignalErr = *YieldsSignalErrPtr;

  if (1) {
    std::cout << "readData: YieldsSignal from <" << fname << ">\n";
    YieldsSignal.Print();
  }

  // Check that the binning is consistent
  res= ( unfolding::checkRangesMY(YieldsSignal,"YieldsSignal from file") == 1 ) &&
    ( unfolding::checkRangesMY(YieldsSignalErr,"YieldsSignalErr from file") == 1 );
  if (res!=1) return 0;

  res = ( unfolding::flattenMatrix(YieldsSignal, vFI) &&
	  unfolding::flattenMatrix(YieldsSignalErr, vErrFI1) );
  vErrFI2=0;

  return 1;
}


//-----------------------------------------------------------------
// Unfold
//-----------------------------------------------------------------
void  applyUnfoldingLocal(TVectorD &vin, TVectorD &vout, bool ifSeed, int seed, int reweightInt)
//if ifSeed==1, smearing systematics
//if ifSeed==0, Fsr 
//reweightInt = 95%, 105%
{

  // Read unfolding constants
  printf("unfold: Load constants\n"); fflush(stdout);


  // Construct file names
  TString fullUnfoldingConstFileName = TString("../root_files/systematics/")+tagDirConstants+TString("/");
  if (ifSeed){
    fullUnfoldingConstFileName += fileUnfoldingConstantsBaseSeed;
    fullUnfoldingConstFileName += seed;
    fullUnfoldingConstFileName += ".root";
    printf("Apply unfolding using unfolding matrix from %s\n", fullUnfoldingConstFileName.Data());
  }
 else{
    fullUnfoldingConstFileName += fileUnfoldingConstantsBaseReweight;    
    fullUnfoldingConstFileName += reweightInt;
    fullUnfoldingConstFileName += ".root";
    printf("Apply unfolding using unfolding matrix from %s\n", fullUnfoldingConstFileName.Data());
 }
    
  if ( unfolding::unfold(vin, vout, fullUnfoldingConstFileName) != 1 ) {
    std::cout << "failed to unfold using matrix from <" << fullUnfoldingConstFileName << ">\n";
    assert(0);
  }

  // Print the result. Mainly for debugging purposes
  if (1) {
    printf("\nUNFOLD: Results for the data, yields:\n");
    printf("                   yields observed        after unfolding            \n");
    for(int i=0, idx=0; i<DYTools::nMassBins; i++){
      double *rapidityBinLimits=DYTools::getYBinLimits(i);
      for (int yi=0; yi<DYTools::nYBins[i]; ++yi, ++idx) {
	printf("%4.0f-%4.0f %4.2lf-%4.2lf   %8.1f       %8.1f\n",
	       DYTools::massBinLimits[i], DYTools::massBinLimits[i+1],
	       rapidityBinLimits[yi],rapidityBinLimits[yi+1],
	       vin[idx], vout[idx]);
      }
      delete rapidityBinLimits;
    }
    printf("\n");
  }

  return;
}

