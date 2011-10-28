#include <TH1F.h>
#include <TFile.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TObject.h>
#include <TMath.h>
#include <TLatex.h>
#include <THStack.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TFrame.h>

#include <sstream>
#include <iostream>
#include <fstream>

#include "../Include/MitStyleRemix.hh"  // style settings for drawing

// Forward declarations
void setHistAttributes(bool doIt, TH1F *hist, int fillColor, int lineColor,
		       Width_t lineWidth);

// Main function

void plotYields(const TString conf = "data.conf", bool plotRatio = true){

  // Find the file with histograms to plot.

  // First find the directory
  //
  // parse .conf file
  //
  TString inputDir;
  Double_t lumi;
  Bool_t doWeight;
  ifstream ifs;
  ifs.open(conf.Data());
  assert(ifs.is_open());
  string line;
  while(getline(ifs,line)) {
    if(line[0]=='#') continue;
    stringstream ss1(line); ss1 >> lumi;
    getline(ifs,line);
    stringstream ss2(line); ss2 >> doWeight;
    getline(ifs,line);
    inputDir = TString(line);
    break;
  }
  ifs.close();

  inputDir.ReplaceAll("selected_events","yields");

  TFile *file = new TFile(inputDir+TString("/massHist.root"));

  TH1F *data  = (TH1F*) file->Get("data"); 
  TH1F *zee   = (TH1F*) file->Get("zee");   bool zeeMc = true;
  TH1F *ztt   = (TH1F*) file->Get("ztt");   bool zttMc = true;
  TH1F *qcd   = (TH1F*) file->Get("qcd");   bool qcdMc = true;
  TH1F *ttbar = (TH1F*) file->Get("ttbar"); bool ttbarMc = true;
  TH1F *wjets = (TH1F*) file->Get("wjets"); bool wjetsMc = true;
  TH1F *ww   = (TH1F*) file->Get("ww");     bool wwMc = true;
  TH1F *wz   = (TH1F*) file->Get("wz");     bool wzMc = true;
  TH1F *zz   = (TH1F*) file->Get("zz");     bool zzMc = true;

  // Make sure that all MC predictions that are expected according
  // to boolean flags above are present. Data has to be present always.
  if( !data) return;
  if( !zee   && zeeMc   ) return;
  if( !ztt   && zttMc   ) return;
  if( !qcd   && qcdMc   ) return;
  if( !ttbar && ttbarMc ) return;
  if( !wjets && wjetsMc ) return;
  if( !ww    && wwMc    ) return;
  if( !wz    && wzMc    ) return;
  if( !zz    && zzMc    ) return;

  // Merge EWK backgrounds  
  TH1F *ewk = (TH1F*) data->Clone("ewk");
  ewk->Reset();
  ewk->SetFillStyle(1001);
  bool ewkMc = (wjetsMc || wwMc || wzMc || zzMc);
  if(wjetsMc)
    ewk->Add(wjets);
  if(wwMc)
    ewk->Add(ww);
  if(wzMc)
    ewk->Add(wz);
  if(zzMc)
    ewk->Add(zz);

  // Create a canvas with pads
  TCanvas *c1 = MakeCanvas("c1","c1",600,600);

  TPad *pad1 = 0;
  TPad *pad2 = 0;
  if(!plotRatio){
    pad1 = (TPad*)c1->GetPad(0);
    pad1->SetRightMargin(0.07); // default is 0.05
  }else{
    c1->Divide(1,2,0,0);
    pad1 = (TPad*)c1->GetPad(1);
    pad1->SetPad(0,0.3,1.0,1.0);
    pad1->SetLeftMargin(0.18);
    pad1->SetTopMargin(0.08);
    pad1->SetRightMargin(0.07);
    pad1->SetBottomMargin(0.01); // All X axis labels and titles are thus cut off

    pad2 = (TPad*)c1->GetPad(2);
    pad2->SetPad(0,0,1.0,0.28);
    pad2->SetLeftMargin(0.18);
    pad2->SetTopMargin(0.01);
    pad2->SetRightMargin(0.07);
    pad2->SetBottomMargin(0.45);
    
  }

  // ===================================
  // Work with the yields distributions
  // ===================================
  pad1->cd();
  pad1->SetLogx(1);
  pad1->SetLogy(1);
  pad1->SetTicky();
  pad1->SetFrameLineWidth(2);

  // Gautier's salamandar pallet
  setHistAttributes(zeeMc  , zee  , kOrange-2 , kOrange+3 , 1);
  setHistAttributes(zttMc  , ztt  , kOrange+7 , kOrange+3, 1);
  setHistAttributes(ewkMc  , ewk  , kOrange+10, kOrange+3, 1);
  setHistAttributes(qcdMc  , qcd  , kViolet-5 , kViolet+3, 1);
  setHistAttributes(ttbarMc, ttbar, kRed+2    , kRed+4   , 1);
  THStack *mc = new THStack("mc","");
  if(qcdMc)
    mc->Add(qcd);
  if(ttbarMc)
    mc->Add(ttbar);
  if(ewkMc)
    mc->Add(ewk);
  if(zttMc)
    mc->Add(ztt);
  if(zeeMc)
    mc->Add(zee);

  // Get the total MC histogram
  TList *list = mc->GetHists();
  TH1F *totalMc = (TH1F*)data->Clone("totalMc");
  totalMc->Reset();
  TIter iter(list);
  TH1F *iterHist;
  while (( iterHist = (TH1F*)iter()))
    totalMc->Add(iterHist);

  // Have to draw stack to have axes defined
  mc->Draw("hist");
  mc->GetXaxis()->SetMoreLogLabels();
  mc->GetXaxis()->SetNoExponent();
  mc->GetXaxis()->SetTitle("M(ee) [GeV]");
  mc->GetYaxis()->SetTitle("entries per bin");
  mc->GetYaxis()->SetTitleOffset(1.2);

  // Custom labels for the X axis
  // Get the original labels out of the way
  mc->GetXaxis()->SetLabelOffset(99);
  const int nMassValues = 7;
  int massValues[nMassValues] = 
    {15, 30, 60, 120, 240, 600, 1500};
  double xLabelNdc[nMassValues];
  double norm = pad1->GetX2() - pad1->GetX1();
  for(int i=0; i<nMassValues; i++){
    // NOTE: the log is needed only if X axis
    // is displayed on the log scale
    xLabelNdc[i] = (TMath::Log10(massValues[i]) - pad1->GetX1())/norm;
  }
  double yLabelNdc = pad1->GetBottomMargin() - 0.02;
  TText xLabelText;
  xLabelText.SetTextFont(42);
  xLabelText.SetTextAlign(23);  

  // This is the *approximate* low edge of the Y axis for
  // the stack histogram. The exact low limit, when in log mode,
  // appears to be rather difficult to set.
  mc->SetMinimum(1.5);

  // Other labels on the plot
  TLatex *cmsText = new TLatex();
  cmsText->SetTextFont(42);
  cmsText->SetTextSize(0.055);
  cmsText->SetTextAlign(31);
  cmsText->SetNDC();
  cmsText->SetText(0.93, 0.94, "CMS Preliminary");

  TLatex *lumiText = new TLatex();
  lumiText->SetTextFont(42);
  lumiText->SetTextSize(0.04);
  lumiText->SetTextAlign(33);
  lumiText->SetNDC();
  lumiText->SetText(0.91, 0.90, "1.1 fb^{-1} at #sqrt{s} = 7 TeV");

  // Put together the legend
  
  TLegend *legend = new TLegend(0.70,0.50, 0.90,0.80 );
  legend->SetTextFont(42);
  legend->SetTextSize(0.03);
  legend->SetLineColor(0);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);
  legend->AddEntry(data,"  data","PL");
  legend->AddEntry(data,"      ","0"); // skip a line
  legend->AddEntry(zee,"#gamma*/Z #rightarrow ee","F");
  legend->AddEntry(ztt,"#gamma*/Z #rightarrow #tau#tau","F");
  legend->AddEntry(ewk,"EWK","F");
  legend->AddEntry(ttbar,"t#bar{t}","F");
  legend->AddEntry(qcd,"QCD","F");

  // Draw all objects
  mc->Draw("hist");
  data->Draw("pe,same");
  totalMc->Draw("same,hist,bar");
  
  for(int i=0; i<nMassValues; i++){
    TString thisLabel = "";
    thisLabel += massValues[i];
    xLabelText.DrawTextNDC(xLabelNdc[i],yLabelNdc,thisLabel);
  }

  cmsText->Draw();
  lumiText->Draw();
  legend->Draw();

  if(plotRatio){
    // ===================================
    // Work with the ratio distribution
    // ===================================
    
    pad2->cd();
    pad2->SetLogx();
    pad2->SetFrameLineWidth(2);
    
    // Find the ratio histogram
    TH1F *ratio = (TH1F*)data->Clone("ratio");
    ratio->Reset();
    ratio->Divide(data,totalMc,1.0,1.0);
    
    // Set attributes
    ratio->SetLineColor(kBlack);
    ratio->SetLineWidth(1);
    ratio->SetMarkerSize(1.0);
    
    double ymin = 0.3;
    double ymax = 1.7;
    ratio->GetYaxis()->SetRangeUser(ymin,ymax);
    ratio->GetYaxis()->SetNdivisions(3);
    ratio->GetYaxis()->SetLabelSize(0.15);
    ratio->GetXaxis()->SetLabelOffset(99);
    ratio->GetXaxis()->SetTickLength(0.10);
    ratio->GetXaxis()->SetTitle(mc->GetXaxis()->GetTitle());
    ratio->GetXaxis()->SetTitleSize(3 * mc->GetXaxis()->GetTitleSize());
    ratio->SetMarkerColor(kOrange+7);
    ratio->SetMarkerStyle(kFullCircle);
    TH1F *ratioClone = (TH1F*)ratio->Clone("ratioClone");
    ratioClone->SetMarkerColor(kBlack);
    ratioClone->SetMarkerStyle(kOpenCircle);
    
    TLine *lineAtOne = new TLine(15.0, 1.0, 1500, 1.0);
    lineAtOne->SetLineStyle(kDashed);
    lineAtOne->SetLineWidth(1);
    lineAtOne->SetLineColor(kBlue);
    
    TLine *verticalLines[nMassValues];
    for(int i=0; i<nMassValues; i++){
      verticalLines[i] = new TLine(massValues[i],ymin,massValues[i],ymax);
      verticalLines[i]->SetLineWidth(1);
      verticalLines[i]->SetLineColor(15);
    }

    ratio->Draw("PE");
    ratioClone->Draw("PE,same");
    
    lineAtOne->Draw();
    
    // Adjust for the smaller pad of the ratio plot
    yLabelNdc = pad2->GetBottomMargin() - 0.06;
    xLabelText.SetTextSize(3 * xLabelText.GetTextSize());
    for(int i=0; i<nMassValues; i++){
      TString thisLabel = "";
      thisLabel += massValues[i];
      xLabelText.DrawTextNDC(xLabelNdc[i],yLabelNdc,thisLabel);
    }

    for(int i=0; i<nMassValues; i++)
      verticalLines[i]->Draw();

    pad2->RedrawAxis();
    pad2->GetFrame()->Draw();


  } // end if plot ratio

  return;
}

void setHistAttributes(bool doIt, TH1F *hist, int fillColor, int lineColor,
		       Width_t lineWidth){

  if(doIt){
    hist->SetFillColor(fillColor);
    hist->SetLineColor(lineColor);
    hist->SetLineWidth(lineWidth);
  }
  return;
}
