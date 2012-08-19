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

//void readData(TVectorD &v, TVectorD &vErr1, TVectorD &vErr2);

void  extractAcceptance(TMatrixD &vout, int reweightInt);

const TString fileEnd(DYTools::analysisTag + TString(".root"));
const TString fileDataYields(TString("yields_bg-subtracted") + fileEnd);
const TString fileAcceptanceConstantsBaseReweight(TString("acceptance_constants_reweight_") + 
					DYTools::analysisTag + TString("_"));


void calcAcceptanceSystematics(const TString conf){

  // check whether it is a calculation
  if (conf.Contains("_DebugRun_")) {
    std::cout << "calcAcceptanceSystematics: _DebugRun_ detected. Terminating the script\n";
    return;
  }

  // normal calculation

  // First, read the configuration file. The configuration file
  // is the same as the one used for the cross section calculation
  // script, we need to know the location of data yields and
  // unfolding matrices

  ifstream ifs;
  ifs.open(conf.Data());
  assert(ifs.is_open());
  std::string line;
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

  /////////////////////////////////
  //calculate Fsr systematics 
  /////////////////////////////////

  TMatrixD accFsrMax(DYTools::nMassBins,DYTools::nYBinsMax);
  TMatrixD accFsrMin(DYTools::nMassBins,DYTools::nYBinsMax);
  TMatrixD accSystPercentFsr(DYTools::nMassBins,DYTools::nYBinsMax);

  extractAcceptance(accFsrMax, 105);
  extractAcceptance(accFsrMin, 95);

  for(int i = 0; i < DYTools::nMassBins; i++){
    for (int yi = 0; yi < DYTools::nYBins[i]; ++yi) {
      accSystPercentFsr[i][yi] = 100.0*fabs(accFsrMax[i][yi]-accFsrMin[i][yi])/
	(accFsrMax[i][yi]+accFsrMin[i][yi]);
    }
  }

  //printing out to the screen

  printf("mass     rapidity    <acc>  rel-err-percent-Fsr\n");
  for(int i=0; i<DYTools::nMassBins; i++){
    double *rapidityBinLimits=DYTools::getYBinLimits(i);
    for (int yi=0; yi<DYTools::nYBins[i]; ++yi) {
      printf("%4.0f-%4.0f  %4.2f-%4.2f  %1.3f   %1.3f \n", 
	     DYTools::massBinLimits[i], DYTools::massBinLimits[i+1],
	     rapidityBinLimits[yi], rapidityBinLimits[yi+1],
	     0.5*(accFsrMax[i][yi]+accFsrMin[i][yi]), accSystPercentFsr[i][yi]);
    }
    delete rapidityBinLimits;
  }


  // Store constants in the file
  TString outputDir(TString("../root_files/systematics/")+tagDirConstants);
  gSystem->mkdir(outputDir,kTRUE);
  TString acceptanceSystFileName(outputDir+TString("/acceptance_FSR_systematics") + fileEnd);

  TFile fa(acceptanceSystFileName,"recreate");
  accSystPercentFsr.Write("accSystPercent");
  accFsrMin.Write("accFSRSystMin");
  accFsrMax.Write("accFsrSystMax");
  unfolding::writeBinningArrays(fa);
  fa.Close();

  return;
}


//-----------------------------------------------------------------
// Acceptance
//-----------------------------------------------------------------
void  extractAcceptance(TMatrixD &vout, int reweightInt)
{
  if ((vout.GetNrows()!=DYTools::nMassBins) || (vout.GetNcols()!=DYTools::nYBinsMax)) {
    std::cout << "extractAcceptance: vout should be of size [DYTools::nMassBins]x[DYTools::nYBinsMax]=[nUnfoldingBins]\n";
    assert(0);
  }

  // Read unfolding constants
  printf("acc: Load constants\n"); fflush(stdout);


  // Construct file names
  TString fullAcceptanceConstFileName = TString("../root_files/systematics/")+tagDirConstants+TString("/");
  fullAcceptanceConstFileName += fileAcceptanceConstantsBaseReweight;    
  fullAcceptanceConstFileName += reweightInt;
  fullAcceptanceConstFileName += ".root";
  printf("Apply acceptace using constants from %s\n", fullAcceptanceConstFileName.Data());
 
  TFile fAcc(fullAcceptanceConstFileName);
  unfolding::checkBinningArrays(fAcc);
  TMatrixD* acc = (TMatrixD*)fAcc.Get("acceptanceMatrix");
  assert(acc);

  for (int i=0; i<DYTools::nMassBins; i++){
    for (int yi=0; yi<DYTools::nYBins[i]; ++yi) {
      vout[i][yi]=(*acc)[i][yi];
    }
  }

  delete acc;
  return;
}

