#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TChain.h>
#include <TCanvas.h>                // class for drawing
#include <TH1F.h>                   // 1D histograms
#include <TGraphErrors.h>           // graphs
#include <TProfile.h>               // profile histograms
#include <TBenchmark.h>             // class to track macro running statistics
#include <TLorentzVector.h>         // class for Lorentz vector computations
#include <TGraphAsymmErrors.h>
#include <TClonesArray.h>
#include <TMatrixD.h>

#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <vector>                   // STL vector class
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings

#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooBreitWigner.h"
#include "RooExponential.h"
#include "RooCBShape.h"
#include "RooGaussian.h"
#include "RooGenericPdf.h"
#include "RooHistPdf.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooFFTConvPdf.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooDataSet.h"
#include "RooFitResult.h"

#include "CPlot.hh"          // helper class for plots
#include "DYTools.hh"

#endif


void printCorrelations(ostream& os, RooFitResult *res);

void fitMass(TTree *passTree, TTree *failTree, TString cut, int mode, 
     double &efficiency, double &efficiencyErrHigh, double &efficiencyErrLow, 
	     TPad *passPad, TPad *failPad, TFile *plotsRootFile,
	     ofstream &fitLog, int NsetBins, 
	     bool isRECO, const char* setBinsType, TString dirTag);

void fitMassWithTemplates(TTree *passTree, TTree *failTree, TString cut, 
			  int mode, 
			  double &efficiency, 
			  double &efficiencyErrHigh, double &efficiencyErrLow, 
			  TPad *passPad, TPad *failPad, TFile *plotsRootFile,
			  ofstream &fitLog, 
			  TH1F *templatePass, TH1F *templateFail, 
			  bool isRECO, const char* setBinsType, 
			  TString dirTag, const TString &picFileExtraTag );
