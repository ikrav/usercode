// This script creates a shape_weights.root file, defining the ratio
// between the MC signal distribution and Data-MCbackgrounds
// The file is used in plotDYUnfoldingMatrix.C with the study flag
// DYTOOLS::ESCALE_RESIDUAL 

#include <TROOT.h>
#include <TH1F.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TSystem.h>
#include "../Include/MitStyleRemix.hh"
#include "../Include/DYTools.hh"

void create_mass_shape_weights(TString dataSelectionTag="DY_m10+pr+a05+o03+pr_4680pb", int plotDistributions=0)
{

  TString sourceFile=TString("../root_files/yields/") + dataSelectionTag + TString("/massHist.root");

  TFile f(sourceFile);
  if (!f.IsOpen()) {
    std::cout << "failed to open a file <" << sourceFile << ">\n";
    std::cout << "Possible causes: \n";
    std::cout << "\t1. Check the dataSelectionTag <" << dataSelectionTag << ">\n";
    std::cout << "\t2. Have you run preparation steps?\n";
    throw 2;
  }
  TH1F * data = (TH1F*)f.Get("data");
  TH1F * ttbar= (TH1F*)f.Get("ttbar");
  TH1F * wjets= (TH1F*)f.Get("wjets");
  TH1F * ww   = (TH1F*)f.Get("ww");
  TH1F * wz   = (TH1F*)f.Get("wz");
  TH1F * zz   = (TH1F*)f.Get("zz");
  TH1F * ztt  = (TH1F*)f.Get("ztt");
  TH1F * qcd  = (TH1F*)f.Get("qcd");
  TH1F * zee  = (TH1F*)f.Get("zee");

  TH1F *signal = (TH1F*)data->Clone("signal");

  signal->Add(ttbar,-1.0);
  signal->Add(wjets,-1.0);
  signal->Add(ww,-1.0);
  signal->Add(wz,-1.0);
  signal->Add(zz,-1.0);
  signal->Add(ztt,-1.0);
  signal->Add(qcd,-1.0);

  TCanvas *c0=(!plotDistributions) ? NULL : MakeCanvas("c0","c0",800,800);
  if (c0) { 
    zee->DrawCopy("a hist"); signal->DrawCopy("same pe");
    c0->SetLogx(1); c0->SetLogy(1); c0->Update(); 
    c0->SaveAs("figInitialMass.png");
  }
  else {
    zee->Draw("hist");
    signal->Draw("same,pe");
  }

  int ifirst = zee->GetXaxis()->FindBin(60.001);
  int ilast = zee->GetXaxis()->FindBin(119.999);
  
  double areaSignal = signal->Integral(ifirst,ilast);
  double areaZee    = zee->Integral(ifirst,ilast);
  zee->Scale(areaSignal/areaZee);

  TCanvas *c1=(!plotDistributions) ? NULL : MakeCanvas("c1","c1",800,800);
  if (c1) { 
    zee->DrawCopy("a hist"); signal->DrawCopy("same pe");
    c1->SetLogx(1); c1->SetLogy(1);
    c1->Update(); 
    c1->SaveAs("figMassAreaMatched.png");
  }
  else {
    zee->Draw("hist");
    signal->Draw("same,pe");
  }

  /*
  const int nWiderBins = 30;
  double widerBins[nWiderBins+1] = 
    {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 
     81, 86, 91, 96, 101, 106, 110, 115, 120, 126, 141, 
     160, 185, 220, 320, 510, 1500}; 
  TH1F *signal_wider = new TH1F("signal_wider","",nWiderBins, widerBins);
  TH1F *zee_wider = new TH1F("zee_wider","",nWiderBins, widerBins);
  
  for(int i=0; i<signal->GetNbinsX(); i++){
    double x = signal->GetBinCenter(i);
    signal_wider->Fill(x, signal->GetBinContent(i));
    zee_wider->Fill(x, zee->GetBinContent(i));
  }
  */
  TH1F* signal_wider=(TH1F*)signal->Clone("signal_wider");
  TH1F* zee_wider=(TH1F*)zee->Clone("zee_wider");

  TH1F *weights = (TH1F*)signal_wider->Clone("weights");
  weights->Divide(zee_wider);
  weights->GetXaxis()->SetTitle("m(e^{+}e^{-}) [GeV/c^{2}]");
  weights->GetYaxis()->SetTitle("data/MC ratio");


  // Reset by hand the first 3 points where normally
  //the fewz correction is needed.
  weights->SetBinContent(1,1.0);
  weights->SetBinContent(2,1.0);
  weights->SetBinContent(3,1.0);

  TCanvas *cWeights=(!plotDistributions) ? NULL : MakeCanvas("cWeights","cWeights",800,800);
  weights->DrawCopy("AP");
  weights->SetDirectory(0);
  if (cWeights) { cWeights->SetLogx(1); cWeights->Update(); cWeights->SaveAs("figWeights.png"); }

  f.Close();

  for (int i=0; i<2; ++i) {
    TString outPath=TString("../root_files/constants/") + dataSelectionTag;
    if (i==1) outPath=TString("../root_files/yields/") + dataSelectionTag;
    gSystem->mkdir(outPath,true);
    TString fname=outPath + TString("/shape_weights.root");
    TFile f2(fname,"recreate");
    weights->Write();
    f2.Close();
    std::cout << "File <" << fname << "> created (copy " << (i+1) << ")\n";
  }
  return;
}
