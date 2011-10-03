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

void  extractAcceptance(TVectorD &vout, int reweightInt);

const TString fileDataYields("yields_bg-subtracted.root");

const TString fileAcceptanceConstantsBaseReweight("acceptance_constants_reweight_");

void calcAcceptanceSystematics(const TString conf){

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

  TVectorD accFsrMax(nMassBins);
  TVectorD accFsrMin(nMassBins);
  TVectorD accSystPercentFsr(nMassBins);


  extractAcceptance(accFsrMax, 105);
  extractAcceptance(accFsrMin, 95);

  for(int i = 0; i < nMassBins; i++){
    accSystPercentFsr[i] = 100.0*fabs(accFsrMax[i]-accFsrMin[i])/(accFsrMax[i]+accFsrMin[i]);
  }

  //printing out to the screen

  printf("            <acc>  rel-err-percent-Fsr\n");
  for(int i=0; i<nMassBins; i++){
    printf("%4.0f-%4.0f   %1.3f   %1.3f \n", 
           massBinLimits[i], massBinLimits[i+1],  
           (accFsrMax[i]+accFsrMin[i])/2, accSystPercentFsr[i]);
  }


  // Store constants in the file
  TString outputDir(TString("../root_files/systematics/")+tagDirConstants);
  gSystem->mkdir(outputDir,kTRUE);
  TString unfoldingSystFileName(outputDir+TString("/acceptance_systematics.root"));

  TFile fa(unfoldingSystFileName,"recreate");
  accSystPercentFsr.Write("accSystPercent");
  fa.Close();

  return;
}


//-----------------------------------------------------------------
// Acceptance
//-----------------------------------------------------------------
void  extractAcceptance(TVectorD &vout, int reweightInt)
{

  // Read unfolding constants
  printf("acc: Load constants\n"); fflush(stdout);


  // Construct file names
  TString fullAcceptanceConstFileName = TString("../root_files/systematics/")+tagDirConstants+TString("/");
  fullAcceptanceConstFileName += fileAcceptanceConstantsBaseReweight;    
  fullAcceptanceConstFileName += reweightInt;
  fullAcceptanceConstFileName += ".root";
  printf("Apply acceptace using constants from %s\n", fullAcceptanceConstFileName.Data());
 
  TFile fAcc(fullAcceptanceConstFileName);
  TVectorD* acc = (TVectorD*)fAcc.Get("acceptanceArray");
  double* accArr=acc->GetMatrixArray();

  for (int i=0; i<DYTools::nMassBins; i++){
   vout[i]=accArr[i];
  }  

  return;
}

