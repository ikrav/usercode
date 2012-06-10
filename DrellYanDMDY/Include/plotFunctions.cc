#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TTree.h>                  // class to access ntuples
#include <TCanvas.h>                // class for drawing
#include <TH1F.h>                   // 1D histograms
#include <TH2D.h>                   // 2D histograms
#include <TStyle.h>                   
#include <TMatrixD.h>
#include <TString.h>
#include <TMultiGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>


#include <vector>                   // STL vector class

#include "../Include/CSample.hh"        // helper class for organizing input ntuple files
#include "../Include/DYTools.hh"        // helper class for organizing input ntuple files
#include "../Include/MyTools.hh"        // helper class for organizing input ntuple files
#include "plotFunctions.hh"
     
#endif

void PlotMatrixVariousBinning(TMatrixD matr, TString name, TString drawOption, TFile *histoFile, const TString &saveDir)
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

  void RShapePlot (TMatrixD relCrossSection, TMatrixD relCrossSectionStatErr, 
TMatrixD relCrossSectionDET, TMatrixD relCrossSectionStatErrDET, 
TMatrixD relPostFsrCrossSection, TMatrixD relPostFsrCrossSectionStatErr, 
TMatrixD relPostFsrCrossSectionDET, TMatrixD relPostFsrCrossSectionStatErrDET)
{

      // create first graph
   //Pre-FSR -All phase space
   double x[nMassBins];
   double ex[nMassBins];
   /*
   double y1[nMassBins];
   double ey1[nMassBins];
   double y2[nMassBins];
   double ey2[nMassBins];
   double y3[nMassBins];
   double ey3[nMassBins];
   double y4[nMassBins];
   double ey4[nMassBins];
   */
   double* y1;
   double* ey1;
   double* y2;
   double* ey2;
   double* y3;
   double* ey3;
   double* y4;
   double* ey4;

   TString saveDir="plots"+ DYTools::analysisTag;
   Int_t n;
   if (DYTools::study2D==0)
     {
       for (int i=0; i<nMassBins; i++)
         {
           x[i]=(massBinLimits[i+1]+massBinLimits[i])/2;
           ex[i]=(massBinLimits[i+1]-massBinLimits[i])/2;
         }
       y1=relCrossSection.GetMatrixArray();
       ey1=relCrossSectionStatErr.GetMatrixArray();
       y2=relCrossSectionDET.GetMatrixArray();
       ey2=relCrossSectionStatErrDET.GetMatrixArray();
       y3=relPostFsrCrossSection.GetMatrixArray();
       ey3=relPostFsrCrossSectionStatErr.GetMatrixArray();
       y4=relPostFsrCrossSectionDET.GetMatrixArray();
       ey4=relPostFsrCrossSectionStatErrDET.GetMatrixArray();

       n = DYTools::nMassBins;
       RShapeDrawAndSave(n,x,ex,y1,ey1,y2,ey2,y3,ey3,y4,ey4,"RShape1D",saveDir);
     }
   else if (DYTools::study2D==1)
     {
       y1  = new double[DYTools::nYBinsMax];
       ey1 = new double[DYTools::nYBinsMax];
       y2  = new double[DYTools::nYBinsMax];
       ey2 = new double[DYTools::nYBinsMax];
       y3  = new double[DYTools::nYBinsMax];
       ey3 = new double[DYTools::nYBinsMax];
       y4  = new double[DYTools::nYBinsMax];
       ey4 = new double[DYTools::nYBinsMax];
/*
        for (int i=0; i<DYTools::nMassBins; i++)
         {
           n = DYTools::nYBins[i];
           for (int j=0; j<DYTools::nYBins[i]; j++)
             {
               x[j]=yRangeMin+(j+0.5)*(yRangeMax-yRangeMin)/n;
               ex[j]=0.5*(yRangeMax-yRangeMin)/n;
               y1[j]=relCrossSection(i,j); 
               ey1[j]=relCrossSectionStatErr(i,j);
               y2[j]=relCrossSectionDET(i,j);
               ey2[j]=relCrossSectionStatErrDET(i,j);
               y3[j]=relPostFsrCrossSection(i,j);
               ey3[j]=relPostFsrCrossSectionStatErr(i,j);
               y4[j]=relPostFsrCrossSectionDET(i,j);
               ey4[j]=relPostFsrCrossSectionStatErrDET(i,j);
             }
           
           TString name2D="RShape2D_mass_";
           name2D+=massBinLimits[i];
           name2D+="-";
           name2D+=massBinLimits[i+1];
           RShapeDrawAndSave(n,x,ex,y1,ey1,y2,ey2,y3,ey3,y4,ey4,name2D,saveDir);
         }
*/
       delete y1;
       delete ey1;
       delete y2;
       delete ey2;
       delete y3;
       delete ey3;
       delete y4;
       delete ey4;

     }
}

void RShapeDrawAndSave(Int_t n, double* x,double* ex,double* y1,double* ey1,double* y2,double* ey2,double* y3,double* ey3,double* y4,double* ey4, TString name, const TString &saveDir)
{

   TCanvas* canv=new TCanvas(name,name);
   canv->SetGrid();
   canv->SetLogx(1);
   canv->SetLogy(1);
   canv->SetFillColor(0);
   // draw a frame to define the range
   TMultiGraph *mg = new TMultiGraph();
   TGraphErrors *gr1 = new TGraphErrors(n,x,y1,ex,ey1);
   gr1->SetName("PreFSR (d#sigma /dM)/ (d#sigma /dM)_{z}");
   gr1->SetTitle("Pre-FSR (d#sigma /dM)/ (d#sigma /dM)_{z}");
   gr1->SetFillColor(0);
   gr1->SetMarkerColor(kBlack);
   gr1->SetMarkerStyle(20);
   gr1->SetMarkerSize(1.0);

   TGraphErrors *gr2 = new TGraphErrors(n,x,y2,ex,ey2);
   gr2->SetName("PreFSR DET (d#sigma /dM)/ (d#sigma /dM)_{z}");
   gr2->SetTitle("Pre-FSR DET (d#sigma /dM)/ (d#sigma /dM)_{z}");
   gr2->SetFillColor(0);
   gr2->SetMarkerColor(kBlue);
   gr2->SetMarkerStyle(20);
   gr2->SetMarkerSize(1.0);
   gr2->SetLineColor(kBlue);

   TGraphErrors *gr3 = new TGraphErrors(n,x,y3,ex,ey3);
   gr3->SetName("PostFSR (d#sigma /dM)/ (d#sigma /dM)_{z}");
   gr3->SetTitle("Post-FSR (d#sigma /dM)/ (d#sigma /dM)_{z}");
   gr3->SetFillColor(0);
   gr3->SetMarkerColor(kRed);
   gr3->SetMarkerStyle(20);
   gr3->SetMarkerSize(1.0);
   gr3->SetLineColor(kRed);

   TGraphErrors *gr4 = new TGraphErrors(n,x,y4,ex,ey4);
   gr4->SetName("PostFSR DET(d#sigma /dM)/ (d#sigma /dM)_{z}");
   gr4->SetTitle("Post-FSR DEt(d#sigma /dM)/ (d#sigma /dM)_{z}");
   gr4->SetFillColor(0);
   gr4->SetMarkerColor(kGreen);
   gr4->SetMarkerStyle(20);
   gr4->SetMarkerSize(1.0);
   gr4->SetLineColor(kGreen);

   mg->Add(gr4);
   mg->Add(gr3);
   mg->Add(gr2);
   mg->Add(gr1);
   mg->Draw("ap");
   TAxis* yax=mg->GetYaxis();
   yax->SetRangeUser(5e-6,2);
   mg->GetXaxis()->SetTitle("M_{ee}");
   mg->GetYaxis()->SetTitle("1/#sigma_{z} d#sigma /dM");
   mg->GetYaxis()->SetTitleOffset(1.20);
   //mg->SetName("1/#sigma_z d#sigma /dM");

   TLegend *leg = new TLegend(.60,.55,.95,.90);
   leg->AddEntry(gr1,"Pre FSR All Phase Space");
   leg->AddEntry(gr2,"Pre FSR Detector Phase space");
   leg->AddEntry(gr3,"Post FSR All Phase Space");
   leg->AddEntry(gr4,"Post FSR Detector Phase space");
   leg->Draw();
   
   
   SaveCanvas(canv, name, saveDir);

}
