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
#include <TLatex.h>


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
		   bool hasData, double dataOverMc, double* dataOverMcEachBin, 
		   bool normEachBin, int singleCanvas,
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
void Plot1D(const TMatrixD &matrValues, const TMatrixD &matrErrors, TString name, TString title);
Int_t minMutualMultiple();
Int_t minMutualMultipleTwo(Int_t n1, Int_t n2);
TMatrixD AdjustMatrixBinning(const TMatrixD &matrUsualBinning);
void PlotMatrixVariousBinning(const TMatrixD &matr, TString name, TString drawOption, TFile *histoFile, bool SetLogz=0);
void PlotMatrixVariousBinning(const TMatrixD &matr, TString name, TString drawOption, TFile *histoFile, TString title, bool SetLogz);


//for cross section
void RShapePlot (TMatrixD relCrossSection, TMatrixD relCrossSectionStatErr, 
TMatrixD relCrossSectionDET, TMatrixD relCrossSectionStatErrDET, 
TMatrixD relPostFsrCrossSection, TMatrixD relPostFsrCrossSectionStatErr, 
TMatrixD relPostFsrCrossSectionDET, TMatrixD relPostFsrCrossSectionStatErrDET);

void RShapeDrawAndSave(Int_t n, double* x,double* ex,double* y1,double* ey1,double* y2,double* ey2,double* y3,double* ey3,double* y4,double* ey4, TString name, double mass1, double mass2);


// -----------------------------------------------------------
//   Utility plots
// -----------------------------------------------------------

void PlotMatrixMYSlices(const std::vector<int> &indices, int functionOfRapidity,
			const std::vector<TMatrixD> &matrV,
			const std::vector<TMatrixD> &matrErrV,
			const std::vector<TString> &labelV,
			TString name, TString drawOption, TFile *histoFile, 
			TString title, TString yAxisLabel="Count",
			int ncolX=-1, int ncolY=-1,
			double yAxisMin=99999, double yAxisMax=-99999);

void PlotMatrixMYSlices(const std::vector<int> &indices, int functionOfRapidity,
			const std::vector<TVectorD> &matrFIV,
			const std::vector<TVectorD> &matrErrFIV,
			const std::vector<TString> &labelV,
			TString name, TString drawOption, TFile *histoFile, 
			TString title, TString yAxisLabel="Count",
			int ncolX=-1, int ncolY=-1,
			double yAxisMin=99999, double yAxisMax=-99999);

void PrintMatrixMYSlices(const std::vector<int> &indices, int functionOfRapidity,
			 const std::vector<TMatrixD> &matrV,
			 const std::vector<TMatrixD> &matrErrV,
			 const std::vector<TString> &labelV);

void PrintMatrixMYSlices(const std::vector<int> &indices, int functionOfRapidity,
		       const std::vector<TVectorD> &matrFIV,
		       const std::vector<TVectorD> &matrErrFIV,
		       const std::vector<TString> &labelV);

// ------------------------------------------------------------------

template<class Container_t>
inline
void PlotMatrixMYSlices(int idx, int functionOfRapidity,
			const Container_t &arr,
			const TString label,
			const TString yAxisLabel="Count") {
  std::vector<int> indices; indices.push_back(idx);
  std::vector<Container_t> arrV; arrV.push_back(arr);
  Container_t err=arr; err=0;
  std::vector<Container_t> errV; errV.push_back(err);
  std::vector<TString> labelV; labelV.push_back(label);
  PlotMatrixMYSlices(indices,functionOfRapidity,arrV,errV,labelV,
		     label, "hist", NULL, "hist", yAxisLabel);
}

// ------------------------------------------------------------------

template<class Container_t>
inline
void PrintMatrixMYSlices(int idx, int functionOfRapidity,
			 const Container_t &arr,
			 const TString label) {
  std::vector<int> indices; indices.push_back(idx);
  std::vector<Container_t> arrV; arrV.push_back(arr);
  Container_t err=arr; err=0;
  std::vector<Container_t> errV; errV.push_back(err);
  std::vector<TString> labelV; labelV.push_back(label);
  PrintMatrixMYSlices(indices,functionOfRapidity,arrV,errV,labelV);
}

// ------------------------------------------------------------------

#endif
