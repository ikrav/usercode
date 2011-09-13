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

void measurePassAndFail(double &signal, double &signalErr, double &efficiency, double &efficiencyErr,TTree *passTree, TTree *failTree,TCanvas *passCanvas, TCanvas *failCanvas,char* setBinsType);

void measureEfficiency(TTree *passTree, TTree *failTree, int method, int etBinning, int etaBinning, TCanvas *canvas, ofstream &effOutput, ofstream &fitLog, bool useTemplates, TFile *templatesFile, TFile *resultsRootFile,int NsetBins, bool isRECO, char* setBinsType, TString dirTag);

void measureEfficiencyCountAndCount(TTree *passTree, TTree *failTree, int etBinning, int etaBinning, TCanvas *canvas, ofstream &effOutput, bool saveResultsToRootFile, TFile *resultsRootFile);

void measureEfficiencyWithFit(TTree *passTree, TTree *failTree, int method, int etBinning, int etaBinning, TCanvas *canvas, ofstream &effOutput, ofstream &fitLog, bool useTemplates, TFile *templatesFile, TFile *resultsRootFile, int NsetBins, bool isRECO, char* setBinsType, TString dirTag);

int getTemplateBin(int etBin, int etaBin, int etaBinning);

TH1F * getPassTemplate(int etBin, int etaBin, int etaBinning, TFile *file);

TH1F * getFailTemplate(int etBin, int etaBin, int etaBinning, TFile *file);
