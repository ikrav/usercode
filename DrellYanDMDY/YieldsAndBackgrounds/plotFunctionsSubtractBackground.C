#include "../Include/plotFunctions.hh"
#include "../Include/MyTools.hh"
#include "../Include/DYTools.hh"

#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH2D.h>
#include <TAxis.h>

void PlotBkgMatrix(TMatrixD bkgRatesPercent)
{
   gStyle->SetPalette(1);
   TString name="bkgRatesPercent"; 

   TCanvas* canv=new TCanvas(name,name);
   canv->SetFillColor(0);
   canv->SetLeftMargin     (0.1);
   canv->SetRightMargin    (0.1);

   int nM=bkgRatesPercent.GetNrows();
   int nY=bkgRatesPercent.GetNcols();
 
   TH2D* bkgHist= new TH2D("bkgHist","bkgHist", nM, massBinLimits2D, nY, DYTools::yRangeMin, DYTools::yRangeMax);
   for (int i=0; i<nM; i++)
     for (int j=0; j<nY; j++)
       {
         bkgHist->SetBinContent(i+1,j+1,bkgRatesPercent(i,j));
       }
  
   bkgHist->SetStats(0);

   bkgHist->SetTitle(name);
   TAxis* xAx=bkgHist->GetXaxis();
   xAx->SetTitle("mass, GeV");
   xAx->SetMoreLogLabels(kTRUE);
   TAxis* yAx=bkgHist->GetYaxis();
   yAx->SetTitle("|Y|");
   canv->SetLogx();

  bkgHist->Draw("COLZ");
  SaveCanvas(canv, name);
  
}
