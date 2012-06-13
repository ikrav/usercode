#ifndef plotFunctions_HH
#define plotFunctions_HH

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TTree.h>                  // class to access ntuples
#include <TCanvas.h>                // class for drawing
#include <TH1F.h>                   // 1D histograms
#include <TH2D.h>                   // 2D histograms
#include <TStyle.h>                   
#include <TMatrixD.h>
#include <TString.h>


#include <vector>                   // STL vector class

#include "../Include/CSample.hh"        // helper class for organizing input ntuple files
#include "../Include/DYTools.hh"        // helper class for organizing input ntuple files
#include "../Include/MyTools.hh"        // helper class for organizing input ntuple files
     
#endif

//for prepareYields
void DrawMassPeak(vector<TH1F*> hMassv, vector<CSample*> samplev, vector<TString> snamev, TH1F* hMassDibosons, bool hasData, 
		  bool mergeDibosons, TString labelDibosons, Int_t colorDibosons, Double_t lumi, char* lumitext, bool actualBinning,
		  TFile *histoFile);

//for prepareYields
void DrawFlattened(vector<TMatrixD*> yields, vector<TMatrixD*> yieldsSumw2, vector<CSample*> samplev, vector<TString> snamev, bool hasData, 
                   bool mergeDibosons, TString labelDibosons, Int_t colorDibosons, Double_t lumi, char* lumitext,
		   TFile *histoFile);

//for prepareYields
void Draw6Canvases(vector<TMatrixD*> yields, vector<TMatrixD*> yieldsSumw2,
		   vector<CSample*> samplev, vector<TString> snamev, 
		   bool hasData, double dataOverMc, double* dataOverMcEachBin, bool normEachBin, bool singleCanvas,
		   TFile *histoFile);
//		   bool hasData, double dataOverMc, double* dataOverMcEachBin, bool normEachBin=1, bool singleCanvas=0,
//		   TFile *histoFile=NULL);

//for prepareYields
void SetSomeHistAttributes (TH1F* hist, TString samplename);

//for subtractBackground
//void PlotBkgMatrix(TMatrixD bkgRatesPercent);

//for acceptance
//void PlotAccMatrix(TMatrixD accv);

//for everything
Int_t minMutualMultiple();
Int_t minMutualMultipleTwo(Int_t n1, Int_t n2);
TMatrixD AdjustMatrixBinning(TMatrixD matrUsualBinning);
void PlotMatrixVariousBinning(TMatrixD matr, TString name, TString drawOption, TFile *histoFile);

//for cross section
void RShapePlot (TMatrixD relCrossSection, TMatrixD relCrossSectionStatErr, 
TMatrixD relCrossSectionDET, TMatrixD relCrossSectionStatErrDET, 
TMatrixD relPostFsrCrossSection, TMatrixD relPostFsrCrossSectionStatErr, 
TMatrixD relPostFsrCrossSectionDET, TMatrixD relPostFsrCrossSectionStatErrDET);

void RShapeDrawAndSave(Int_t n, double* x,double* ex,double* y1,double* ey1,double* y2,double* ey2,double* y3,double* ey3,double* y4,double* ey4, TString name);


#endif
