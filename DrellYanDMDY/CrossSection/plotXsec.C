//// Run by following commands
////$ root -l 
////root [0] .L plotXsec.C+
////root [1] plotXsec();
////////////////////////
#include "TCanvas.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TAxis.h"
#include "TAxis.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH1D.h"
#include "TFile.h"
#include "TPad.h"
#include "TObject.h"
#include "TMath.h"
#include "TLatex.h"
#include "THStack.h"
#include "TLegend.h"
#include "TLine.h"
#include "TFrame.h"
#include "TString.h"

#include <sstream>
#include <iostream>
#include <fstream>

#include "../Include/MitStyleRemix.hh"  // style settings for drawing
#include "../Include/DYTools.hh"
#include "../Include/UnfoldingTools.hh"
#include "../Include/MyTools.hh"        // miscellaneous helper functions
#include "../Include/TriggerSelection.hh"

void readData(TVectorD &v, TVectorD &vErr1, TVectorD &vErr2, const TriggerSelection &triggers);
void readTh(TVectorD &vTh, TVectorD &vThErr, const TriggerSelection &triggers);
const int nMassBinTh=518;

// Forward declarations
//void setHistAttributes(bool doIt, TH1F *hist, int fillColor, int lineColor,	Width_t lineWidth);

// Main function

void plotXsec(){
  TriggerSelection triggers("Full2011_hltEffOld",true,0);
  TVectorD xSec(DYTools::nMassBins);
  TVectorD xSecErr(DYTools::nMassBins);
  TVectorD xSecErrSyst(DYTools::nMassBins);

  TVectorD BinLimitsForXsec(DYTools::nMassBins+1);
  for(int ibin=0; ibin<=DYTools::nMassBins; ibin++){
    BinLimitsForXsec[ibin] = DYTools::massBinLimits[ibin];
  }

  TVectorD xSecTh(nMassBinTh);
  TVectorD xSecThErr(nMassBinTh);

//-----------------------------------------------------------------
// Read data
//-----------------------------------------------------------------
  readData(xSec, xSecErr, xSecErrSyst, triggers);

  // Create a canvas with pads
  TCanvas *c1 = MakeCanvas("cXsec","cXsec",600,600);
  //c1->SetGrid();
  c1->SetLogx(1);
  c1->SetFillColor(0);
  c1->SetLogy(1);
  c1->SetTickx(1);
  c1->SetTicky(1);

  // draw a frame to define the range
  TMultiGraph *mg = new TMultiGraph();
        // create first graph
  const Int_t n1 = DYTools::nMassBins, n2=nMassBinTh;
  Double_t x1[n1], x3[n2], ex1[n1], ex3[n2] ;
  Double_t y1[n1], y2[n2],  y3[n2], ey1[n1], ey2[n2], ey3[n2];
////  printf("Bin  MeanMass  Width        Xsec         ErrXsec\n");
  //Int_t i1=-1, i2=0;
  for (int iL=0; iL < n1; iL++ )
  {
      x1[iL]=(BinLimitsForXsec[iL+1]+BinLimitsForXsec[iL])/2;
      ex1[iL]=(BinLimitsForXsec[iL+1]-BinLimitsForXsec[iL])/2;
      y1[iL]=fabs(xSec[iL]);
      ey1[iL]=xSecErr[iL];
////      printf(" %4.2f     %4.2f       %1.7f      %1.7f\n",
////	   x1[iL], ex1[iL], y1[iL],ey1[iL]);
  }
  Double_t peak_bin_width = 1.;
  Double_t peak_val_theory = 1009.0;
  Double_t mass_xlow4[n2+1];
  Double_t mxl = 14.0;
  for( int i = 0; i < n2+1; i++ ) {

    if     ( i >=   0 && i <  11 ) {mxl += 1.0;}
    else if( i >=  11 && i <  18 ) {mxl += 5.0;}
    else if( i >=  18 && i < 118 ) {mxl += 1.0;}
    else if( i >= 118 && i < 340 ) {mxl += 2.0;}
    else if( i >= 340 && i < n2)   {mxl += 5.0;}
    else if( i == n2)              {mxl = 1500; }
     mass_xlow4[i] = mxl;
////     if (i > 10  && i < 330) continue;
////     cout <<"i= "<<i <<" mxl= "<<mass_xlow4[i] << endl;
  }
  readTh(xSecTh, xSecThErr, triggers);
  for (int iL2=0; iL2 < n2; iL2++ )
  {   
      x3[iL2]=(mass_xlow4[iL2+1]+mass_xlow4[iL2])/2;
      ex3[iL2]=mass_xlow4[iL2+1]-mass_xlow4[iL2];
      y2[iL2]=xSecTh[iL2];
      ey2[iL2]=xSecThErr[iL2];
      y3[iL2]=y2[iL2]*peak_bin_width/(peak_val_theory*ex3[iL2]);//divide to the Z peak
      ey3[iL2]=sqrt(pow(ey2[iL2]/(peak_val_theory*ex3[iL2]),2)+pow(y2[iL2]/(peak_val_theory*2*ex3[iL2]),2));
  }

   TGraphErrors *gr1 = new TGraphErrors(n1,x1,y1,ex1,ey1);

   gr1->SetFillColor(0);
   gr1->SetMarkerColor(kBlack);
   gr1->SetMarkerStyle(20);
   gr1->SetMarkerSize(1.0);
   gr1->SetLineColor(kBlack);

   TGraphErrors *gr2 = new TGraphErrors(n2,x3,y3);

   gr2->SetFillColor(0);
   gr2->SetMarkerColor(kBlue);
   //gr2->SetMarkerStyle(20);
   gr2->SetLineWidth(2);
   gr2->SetLineColor(kBlue);
////   gr2->SetFillColor(kBlue);
////   gr2->SetFillStyle(3001);

   mg->Add(gr2,"L");
   mg->Add(gr1,"Psame");
   //mg->Add(gr2);
   //mg->Add(gr1);
   //gr1->Draw("ap");
   mg->Draw("A");
   mg->GetYaxis()->SetRangeUser(5e-11,1.0);
   mg->GetXaxis()->SetRangeUser(15,1500);
   mg->GetXaxis()->SetTitle("M_{ee} [GeV]");
   mg->GetYaxis()->SetTitle("1/#sigma_{z} d#sigma/dM [GeV^{-1}]");

   TLatex *cmsText = new TLatex();
   cmsText->SetTextFont(42);
   cmsText->SetTextSize(0.055);
   cmsText->SetTextAlign(31);
   cmsText->SetNDC();
   cmsText->SetText(0.93, 0.94, "CMS Preliminary");

   TLatex *chText = new TLatex();
   chText->SetTextFont(42);
   chText->SetTextSize(0.055);
   chText->SetTextAlign(33);
   chText->SetNDC();
   chText->SetText(0.90, 0.80, "#gamma*/Z #rightarrow ee");

   TLatex *lumiText = new TLatex();
   lumiText->SetTextFont(42);
   lumiText->SetTextSize(0.04);
   lumiText->SetTextAlign(33);
   lumiText->SetNDC();
   lumiText->SetText(0.91, 0.90, "4.7 fb^{-1} at #sqrt{s} = 7 TeV");

   TLegend *leg = new TLegend(.20,.20,.60,.30);
   leg->AddEntry(gr1,"Data ee 4.7 fb^{-1} 2011 ");
   leg->AddEntry(gr2,"NNLO, FEWZ+MSTW08");
   leg->SetTextSize(0.03);
   leg->SetFillColor(0);
   leg->SetLineColor(0);
   leg->SetShadowColor(0);
   leg->Draw();

   cmsText->Draw();
   chText->Draw();
   lumiText->Draw();

   SaveCanvas(c1,"cXsec");
  return;
}


void readData(TVectorD &v, TVectorD &vErr1, TVectorD &vErr2, const TriggerSelection &triggers){

  //printf("Load data yields\n"); fflush(stdout);
  std::cout << "Load data yields" << std::endl;
  TString xSecResultFileName(TString("../root_files/xSec_results_") + 
			     DYTools::analysisTag + TString("_") +
		   triggers.triggerConditionsName() + TString(".root"));
  std::cout << "xSecResultFileName= " << xSecResultFileName << "\n";

  //TFile fileXsecResult   (TString("../root_files/xSec_results.root"));
  TFile fileXsecResult   (xSecResultFileName);
  TMatrixD xSecM          = *(TMatrixD *)fileXsecResult.FindObjectAny("normXSecByBin");
  TMatrixD xSecErr1M      = *(TMatrixD *)fileXsecResult.FindObjectAny("normXSecErrByBin");
  TMatrixD xSecErr2M      = *(TMatrixD *)fileXsecResult.FindObjectAny("normXSecErrByBinSyst");

  const int nUnfoldingBins=DYTools::getTotalNumberOfBins();
  TVectorD xSec(nUnfoldingBins);
  TVectorD xSecErr1(nUnfoldingBins);
  TVectorD xSecErr2(nUnfoldingBins);

  std::cout << "nUnfoldingBins=" << nUnfoldingBins << ", nMassBins=" << DYTools::nMassBins << "\n";

  unfolding::flattenMatrix(xSecM, xSec);
  unfolding::flattenMatrix(xSecErr1M, xSecErr1);
  unfolding::flattenMatrix(xSecErr1M, xSecErr1);

  // Check that the binning is consistent
  bool checkResult = true;
  if( v.GetNoElements() != DYTools::nMassBins ) checkResult = false;
  if( !checkResult ){
    printf("ERROR: Data inconsistent binning in the inputs\n");
    assert(0);
  }else
    printf("readData: Binning in the inputs is consistent\n");

  // Prepare output yields and errors
  for(int i=0; i<DYTools::nMassBins; i++){
    v[i] = xSec[i];
    vErr1[i] = xSecErr1[i];
    vErr2[i] = xSecErr2[i];
  }

  fileXsecResult.Close();
  return;
}
void readTh(TVectorD &v, TVectorD &vErr, const TriggerSelection &triggers){

  //printf("Load data yields\n"); fflush(stdout);
  //TFile fileXsecTh   (TString("../root_files/xSecTh_results.root"));

  TString xSecThResultFileName(TString("../root_files/xSecTh_results_") +
			       DYTools::analysisTag + TString("_") +
		      triggers.triggerConditionsName() + TString(".root"));
  std::cout << "Load theory predictions\n";
  std::cout << "xSecThResultFileName=" << xSecThResultFileName << std::endl;

  TFile fileXsecTh   (xSecThResultFileName);
  TVectorD xSecTh          = *(TVectorD *)fileXsecTh.FindObjectAny("XSecTh");
  TVectorD xSecThErr      = *(TVectorD *)fileXsecTh.FindObjectAny("XSecThErr");

  // Check that the binning is consistent
  bool checkResult = true;
  if( v.GetNoElements() != nMassBinTh ) checkResult = false;
  if( !checkResult ){
    printf("ERROR: Th inconsistent binning in the inputs\n");
    assert(0);
  }else
    printf("readTh: Binning in the inputs is consistent\n");

  // Prepare output yields and errors
  for(int i=0; i<nMassBinTh; i++){
    v[i] = xSecTh[i];
    vErr[i] = xSecThErr[i];
  }

  fileXsecTh.Close();
  return;
}
