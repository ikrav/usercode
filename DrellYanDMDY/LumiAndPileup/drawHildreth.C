//
// This script draws an overlay of the pile-up distributions
// used in the primary part of the analysis: for signal data
// and MC.
//

void drawHildreth(){

  TH1F *h1;
  
  TFile f1("../root_files/pileup/mcPileupHildreth_full2011_20121110_repacked.root");
  h1 = (TH1F*)f1.Get("pileup_simulevel_mc");
  TH1F *hmc = (TH1F*)h1->Clone("hmc");
  hmc->SetDirectory(0);
  f1.Close();

  TFile f2("../root_files/pileup/dataPileupHildreth_full2011_20121110_repacked_default.root");
  h1 = (TH1F*)f2.Get("pileup_lumibased_data");
  TH1F *hdataReco = (TH1F*)h1->Clone("hdataReco");
  hdataReco->SetDirectory(0);
  f2.Close();

  hdataReco->Scale(1.0/hdataReco->GetSumOfWeights());
  hmc->Scale(1.0/hmc->GetSumOfWeights());

  hdataReco->GetXaxis()->SetRangeUser(0.0,40.0);


  hmc->SetLineColor(kOrange+7);
  hmc->SetFillColor(kOrange+7);
  hmc->SetFillStyle(1001);
  
  hdataReco->GetXaxis()->SetTitle("pile-up interactions");
  hdataReco->GetXaxis()->SetTitleOffset(1.0);

  TCanvas *c = MakeCanvas("c","c");

  hdataReco->Draw("hist");

  hmc->Draw("hist,same");

  hdataReco->Draw("same,hist");
  
  TLegend *leg = new TLegend(0.6,0.7,0.85, 0.85);
  leg->AddEntry(hdataReco, "2011 data","l");
  leg->AddEntry(hmc, "Fall 2011 MC", "f");
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw();

}

