#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TCanvas.h>                // class for drawing
#include <TH1F.h>                   // 1D histograms
#include <TGraphErrors.h>           // graphs
#include <TProfile.h>               // profile histograms
#include <TBenchmark.h>             // class to track macro running statistics
#include <TLorentzVector.h>         // class for Lorentz vector computations
#include <TVectorD.h>
#include <TGraphErrors.h>
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <vector>                   // STL vector class
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings

#include "../Include/CPlot.hh"          // helper class for plots
#include "../Include/MitStyleRemix.hh"  // style settings for drawing
#include "../Include/MyTools.hh"        // miscellaneous helper functions
#include "../Include/DYTools.hh"

// define classes and constants to read in ntuple
#include "../Include/EWKAnaDefs.hh"
#include "../Include/TGenInfo.hh"
#endif

//=== FUNCTION DECLARATIONS ======================================================================================

//=== MAIN MACRO =================================================================================================

void TheoryErrors(const TString input)
{
  gBenchmark->Start("TheoryErrors");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================
  
  //Bool_t doSave  = false;    // save plots?
  TString format = "png";   // output file format
  
  vector<TString> fnamev;   // file names   
  vector<TString> labelv;   // legend label
  vector<Int_t>   colorv;   // color in plots
  vector<Int_t>   linev;    // line style
  vector<Double_t> xsecv;
  vector<Double_t> lumiv;
  TString          dirTag;

  //Double_t massLow  = DYTools::massBinLimits[0];
  //Double_t massHigh = DYTools::massBinLimits[nMassBins];
  
  ifstream ifs;
  ifs.open(input.Data());
  //assert(ifs.is_open());
  string line;
  Int_t state=0;
  while(getline(ifs,line)) {
    if(line[0]=='#') continue;
    if(state == 0){
      dirTag = TString(line);
      state++;
      continue;
    }else{
      string fname;
      Int_t color, linesty;
      stringstream ss(line);
      Double_t xsec;
      ss >> fname >> xsec >> color >> linesty;
      string label = line.substr(line.find('@')+1);
      fnamev.push_back(fname);
      labelv.push_back(label);
      colorv.push_back(color);
      linev.push_back(linesty);
      xsecv.push_back(xsec);
      lumiv.push_back(0);
    }
  }
  ifs.close();
  
  //const Double_t kGAP_LOW  = 1.4442;
  //const Double_t kGAP_HIGH = 1.566;
  
  //--------------------------------------------------------------------------------------------------------------
  // Positive and negative errors on acceptance
  //==============================================================================================================

const double accErrpv[40]=
              {0.0202589,0.0166771,0.0162125,0.0158289,0.0152556,0.0150349,0.0152044,0.0150686,0.0148437,0.0147129
              ,0.0143715,0.0139860,0.0135677,0.0133275,0.0129153,0.0128559,0.0127118,0.0127453,0.0127401,0.0127825
              ,0.0126450,0.0127544,0.0127726,0.0127865,0.0126603,0.0124779,0.0123707,0.0121409,0.0116166,0.0114122
              ,0.0107999,0.0101864,0.0093395,0.0086239,0.0074198,0.0063859,0.0054075,0.0044680,0.0031472,0.0020292};
const double  accErrmv[40]=
              {0.0287934,0.0261841,0.0270555,0.0272950,0.0264058,0.0254006,0.0251779,0.0244064,0.0239047,0.0236810
              ,0.0231674,0.0229479,0.0224777,0.0223520,0.0218297,0.0216087,0.0212085,0.0206289,0.0200064,0.0197212
              ,0.0190722,0.0185641,0.0181809,0.0177539,0.0171240,0.0164836,0.0159473,0.0155092,0.0145961,0.0140976
              ,0.0131865,0.0120779,0.0107776,0.0095906,0.0076742,0.0063116,0.0051286,0.0042477,0.0032955,0.0018471};



  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================
 
  TVectorD accErrv  (DYTools::nMassBins);
  TVectorD accErrPlusv   (DYTools::nMassBins);
  TVectorD accErrMinusv  (DYTools::nMassBins);

  accErrv   = 0;
  for(int i=0; i<DYTools::nMassBins; i++){
      accErrPlusv[i]=accErrpv[i];
      accErrMinusv[i]=accErrmv[i];
      //accErrv[i] = sqrt(accErrpv[i]*accErrpv[i]+accErrmv[i]*accErrmv[i])/2;
      accErrv[i] = (accErrpv[i]+accErrmv[i])/2;
  }

   // Store constants in the file
  TString outputDir(TString("../root_files/systematics/")+dirTag);
  gSystem->mkdir(outputDir,kTRUE);
  TString accConstFileName(outputDir+TString("/theoretical_uncertainties.root"));

  TFile fa(accConstFileName,"recreate");
  accErrPlusv.Write("acceptanceTheoryErrPlusArray");
  accErrMinusv.Write("acceptanceTheoryErrMinusArray");
  accErrv.Write("acceptanceTheoryErrArray");
  fa.Close();

  //--------------------------------------------------------------------------------------------------------------
  // Summary print out
  //==============================================================================================================
  cout << endl;
  cout << "*--------------------------------------------------" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  cout << endl; 

  printf(" mass bin  Bin     Errp      Errm      ErrAverage\n");
  for(int i=0; i<DYTools::nMassBins; i++){
    printf(" %4.0f-%4.0f  %i  +%1.7f  -%1.7f  %1.7f \n",
	   DYTools::massBinLimits[i], DYTools::massBinLimits[i+1],
	   i, accErrPlusv[i], accErrMinusv[i], accErrv[i]);
  }
  cout << endl;
}
