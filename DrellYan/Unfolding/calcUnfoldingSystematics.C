#include "TVectorD.h"
#include "TFile.h"
#include "TString.h"
#include "TSystem.h"
#include <iostream>
#include <fstream>
#include <sstream>

#include "../Include/DYTools.hh"
#include "../Include/UnfoldingTools.hh"

// The global variables will be used in several functions:
TString tagDirYields = "";
TString tagDirConstants = "";
Double_t lumi = 0;

void readData(TVectorD &v, TVectorD &vErr1, TVectorD &vErr2);

void  applyUnfoldingShort(TVectorD &vin, TVectorD &vout, int seed);

const TString fileDataYields        ("yields_bg-subtracted.root");

const TString fileUnfoldingConstantsBase("unfolding_constants_seed_");
const int seedFirst = 1001;
const int seedLast = 1020;

void calcUnfoldingSystematics(const TString conf){

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


  TVectorD signalYields(nMassBins);
  TVectorD signalYieldsStatErr(nMassBins);
  TVectorD signalYieldsSystErr(nMassBins);

  TVectorD unfoldedYields(nMassBins);

  TVectorD unfoldedYieldsMean(nMassBins);
  TVectorD unfoldedYieldsRMS(nMassBins);
  TVectorD unfoldedYieldsSquaredMean(nMassBins);

  unfoldedYieldsMean = 0;
  unfoldedYieldsRMS = 0;
  unfoldedYieldsSquaredMean = 0;

  // Read data yields from file
  readData(signalYields, signalYieldsStatErr, signalYieldsSystErr);
  
  int nseeds = 0;
  for(int i=seedFirst; i<=seedLast; i++){
    nseeds++;
    applyUnfoldingShort(signalYields, unfoldedYields, i);
    for(int iMassBin = 0; iMassBin < nMassBins; iMassBin++){
      unfoldedYieldsMean[iMassBin] += unfoldedYields[iMassBin];
      unfoldedYieldsSquaredMean[iMassBin] += unfoldedYields[iMassBin]*unfoldedYields[iMassBin];
    }
  }
  // Final calculation of the mean and RMS
  TVectorD unfoldingSystPercent(nMassBins);
  for(int iMassBin = 0; iMassBin < nMassBins; iMassBin++){
    unfoldedYieldsMean[iMassBin] = unfoldedYieldsMean[iMassBin]/(1.0*nseeds);
    unfoldedYieldsSquaredMean[iMassBin] = unfoldedYieldsSquaredMean[iMassBin]/(1.0*nseeds);
    unfoldedYieldsRMS[iMassBin] = sqrt(unfoldedYieldsSquaredMean[iMassBin] - 
				unfoldedYieldsMean[iMassBin]*unfoldedYieldsMean[iMassBin]);
    unfoldingSystPercent[iMassBin] = unfoldedYieldsRMS[iMassBin]*100.0/unfoldedYieldsMean[iMassBin];
  }
  
  printf("               mean-unfolded   RMS-unfolded   rel-error     rel-error-in-percent\n");
  for(int i=0; i<nMassBins; i++){
    printf("%4.0f-%4.0f   %7.1f      %7.1f          %6.4f        %6.1f\n",
	   massBinLimits[i],massBinLimits[i+1],
	   unfoldedYieldsMean[i], unfoldedYieldsRMS[i],
	   unfoldedYieldsRMS[i]/unfoldedYieldsMean[i],
	   unfoldingSystPercent[i]);
// 	   unfoldedYieldsRMS[i]*100.0/unfoldedYieldsMean[i]);
  }
  
  // Store constants in the file
  TString outputDir(TString("../root_files/systematics/")+tagDirConstants);
  gSystem->mkdir(outputDir,kTRUE);
  TString unfoldingSystFileName(outputDir+TString("/unfolding_systematics.root"));

  TFile fa(unfoldingSystFileName,"recreate");
  unfoldingSystPercent.Write("unfoldingSystPercent");
  fa.Close();

  return;
}

//-----------------------------------------------------------------
// Read data
//-----------------------------------------------------------------
void readData(TVectorD &v, TVectorD &vErr1, TVectorD &vErr2){

  printf("Load data yields\n"); fflush(stdout);
  TFile fileYields   (TString("../root_files/yields/")+tagDirYields+TString("/")+fileDataYields);
  TVectorD YieldsSignal       = *(TVectorD *)fileYields.FindObjectAny("YieldsSignal");
  TVectorD YieldsSignalErr    = *(TVectorD *)fileYields.FindObjectAny("YieldsSignalErr");
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
    // Presently zero, but in the future non-zero
    // (e.g. background systematics, or energy scale correction effect on yields)
    vErr2[i] = 0;
  } 

  fileYields.Close();
  return;
}

//-----------------------------------------------------------------
// Unfold
//-----------------------------------------------------------------
void  applyUnfoldingShort(TVectorD &vin, TVectorD &vout, int seed)
{

  // Read unfolding constants
  printf("unfold: Load constants\n"); fflush(stdout);


  // Construct file names
  TString fullUnfoldingConstFileName = TString("../root_files/systematics/")
    +tagDirConstants+TString("/")+fileUnfoldingConstantsBase;
  fullUnfoldingConstFileName += seed;
  fullUnfoldingConstFileName += ".root";
  printf("Apply unfolding using unfolding matrix from %s\n", fullUnfoldingConstFileName.Data());
    
  unfolding::unfold(vin, vout, fullUnfoldingConstFileName);

  // Print the result. Mainly for debugging purposes
  printf("\nUNFOLD: Results for the data, yields:\n");
  printf("                   yields observed        after unfolding            \n");
  for(int i=0; i<DYTools::nMassBins; i++){
    printf("%4.0f-%4.0f   %8.1f       %8.1f\n",
	   DYTools::massBinLimits[i], DYTools::massBinLimits[i+1],
	   vin[i], vout[i]);
  }
  printf("\n");

  return;
}

