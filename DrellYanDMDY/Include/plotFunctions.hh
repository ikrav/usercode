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

void PlotMatrixVariousBinning(TMatrixD matr, TString name, TString drawOption, TFile *histoFile, const TString &saveDir="plots")
{
//acceptance, bkgRatesPercent, LEGO2,COLZ 
   TMatrixD matrDraw = AdjustMatrixBinning(matr);
   gStyle->SetPalette(1);

   TCanvas* canv=new TCanvas(name,name);
   canv->SetFillColor(0);
   canv->SetLeftMargin     (0.1);
   canv->SetRightMargin    (0.1);

   int nM=matr.GetNrows();
   int nY=matr.GetNcols();

   TString histName = "Hist" + name;
   TH2D* Hist= new TH2D(histName,name, nM, DYTools::massBinLimits, nY, DYTools::yRangeMin, DYTools::yRangeMax);
   for (int i=0; i<nM; i++)
     for (int j=0; j<nY; j++)
       {
         Hist->SetBinContent(i+1,j+1,matrDraw(i,j));
       }
  
   Hist->SetStats(0);

   Hist->SetTitle(name);
   TAxis* xAx=Hist->GetXaxis();
   xAx->SetTitle("mass, GeV");
   //xAx->SetMoreLogLabels(kTRUE);
   TAxis* yAx=Hist->GetYaxis();
   yAx->SetTitle("|Y|");
   canv->SetLogx();

  Hist->Draw(drawOption);
  SaveCanvas(canv, name, saveDir);
  if (histoFile) canv->Write();
}

TMatrixD AdjustMatrixBinning(TMatrixD matrUsualBinning)
{
  TMatrixD matrOut(DYTools::nMassBins,minMutualMultiple());
  for (int i=0; i<DYTools::nMassBins; i++)
    {
      int nTheSameCells=minMutualMultiple()/DYTools::nYBins[i];
      for (int j=0; j<DYTools::nYBins[i]; j++)
        for (int l=0; l<nTheSameCells; l++)          
           matrOut(i,j*nTheSameCells+l)=matrUsualBinning(i,j);     
    }
  return matrOut;
}

Int_t minMutualMultiple()
{
  int mult=nYBins[0];
  for (int i=1; i<nMassBins; i++)
  {
    mult=minMutualMultipleTwo(mult,nYBins[i]);

  }
  return mult;
}

Int_t minMutualMultipleTwo(Int_t n1, Int_t n2)
{
  Int_t subMultiple=1;
  Int_t min; 
  if (n1<n2) min=n1;
  else if (n2<=n1) min=n2;
  for (int i=1; i<=min; i++)
    {
      if (((n1%i)==0) && ((n2%i)==0))
        subMultiple=i;
    }
  return subMultiple*(n1/subMultiple)*(n2/subMultiple);
}
