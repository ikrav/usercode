#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TChain.h>
#include <TString.h>
#include <TCanvas.h>                // class for drawing

#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <vector>                   // STL vector class
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings

#include "fitFunctionsCore.hh"

#endif

//void measurePassAndFail(double &signal, double &signalErr, 
// double &efficiency, double &efficiencyErr,
// TTree *passTree, TTree *failTree,TCanvas *passCanvas, TCanvas *failCanvas,
//const char* setBinsType);

void measureEfficiency(TTree *passTree, TTree *failTree, 
	       int method, int etBinning, int etaBinning, TCanvas *canvas, 
		       ofstream &effOutput, ofstream &fitLog, 
		       bool useTemplates, TFile *templatesFile, 
		       TFile *resultsRootFile, TFile *plotsRootFile,
		       int NsetBins, DYTools::TEfficiencyKind_t effType, 
		       const char* setBinsType, 
		       TString dirTag, const TString &picFileExtraTag, 
		       int puBin=-1); // puBin is important for the fit

void measureEfficiencyPU(TTree *passTreeFull, TTree *failTreeFull, 
		 int method, int etBinning, int etaBinning, TCanvas *canvas, 
			 ofstream &effOutput, ofstream &fitLog, 
			 bool useTemplates, TFile *templatesFile, 
			 const TString &resultRootFileBase,
			 int NsetBins, DYTools::TEfficiencyKind_t effType,
			 const char* setBinsType, 
			 TString dirTag, const TString &picFileExtraTag, 
			 int puDependence=0);

void measureEfficiencyCountAndCount(TTree *passTree, TTree *failTree, 
			    int etBinning, int etaBinning, 
			    TCanvas *canvas, ofstream &effOutput, 
			    bool saveResultsToRootFile, TFile *resultsRootFile,
				    TFile *plotsRootFile, 
				    DYTools::TEfficiencyKind_t effType);

void measureEfficiencyWithFit(TTree *passTree, TTree *failTree, 
			      int method, int etBinning, int etaBinning, TCanvas *canvas, 
			      ofstream &effOutput, ofstream &fitLog, 
			      bool useTemplates, TFile *templatesFile, 
			      TFile *resultsRootFile, TFile *plotsRootFile,
			      int NsetBins, DYTools::TEfficiencyKind_t effType,
			      const char* setBinsType,
			      TString dirTag, const TString &picFileExtraTag, int puBin=-1);

int getTemplateBin(int etBin, int etaBin, int etaBinning);


TString getTemplateName(int etBin, int etaBin, const char *pass_fail_str, 
			int puBin=-1);

TH1F * getPassTemplate(int etBin, int etaBin, int etaBinning, TFile *file, 
		       int puBin=-1);

TH1F * getFailTemplate(int etBin, int etaBin, int etaBinning, TFile *file, 
		       int puBin=-1);
