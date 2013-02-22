//
// This script draws the plot of the pile-up distributions
// from data: the primary one, and the ones obtained with
// +- 5 % variation of the inelastic pp cross section at given energy.
//

void drawHildrethSyst(){

  TH1F *h1;
  
  TFile f1("../root_files/pileup/dataPileupHildreth_full2011_20121110_repacked_default.root");
  h1 = (TH1F*)f1.Get("pileup_lumibased_data");
  TH1F *hdata = (TH1F*)h1->Clone("hdata");
  hdata->SetDirectory(0);
  f1.Close();

  TFile f2("../root_files/pileup/dataPileupHildreth_full2011_20130214_repacked_plus5percent.root");
  h1 = (TH1F*)f2.Get("pileup_lumibased_data");
  TH1F *hdataUp = (TH1F*)h1->Clone("hdataUp");
  hdataUp->SetDirectory(0);
  f2.Close();

  TFile f3("../root_files/pileup/dataPileupHildreth_full2011_20130214_repacked_minus5percent.root");
  h1 = (TH1F*)f3.Get("pileup_lumibased_data");
  TH1F *hdataDown = (TH1F*)h1->Clone("hdataDown");
  hdataDown->SetDirectory(0);
  f3.Close();


  hdataUp->Scale(1.0/hdataUp->GetSumOfWeights());
  hdataDown  ->Scale(1.0/hdataDown  ->GetSumOfWeights());
  hdata ->Scale(1.0/hdata ->GetSumOfWeights());

  hdata->SetLineWidth(2);
  hdataUp->SetLineWidth(2);
  hdataDown->SetLineWidth(2);

  hdataDown->GetXaxis()->SetRangeUser(0.0,30.0);

  hdataUp->SetLineColor(kBlue);
  hdataDown->SetLineColor(kRed);
  
  hdataDown->GetXaxis()->SetTitle("pile-up interactions");
  hdataDown->GetXaxis()->SetTitleOffset(1.0);

  TCanvas *c = MakeCanvas("c","c");

  hdataDown->Draw("hist");
  hdataUp->Draw("hist,same");
  hdata->Draw("hist,same");

  
  TLegend *leg = new TLegend(0.6,0.7,0.85, 0.85);
  leg->AddEntry(hdata     , "standard #sigma_{pp}","l");
  leg->AddEntry(hdataUp   , "#sigma_{pp} + 5%","l");
  leg->AddEntry(hdataDown , "#sigma_{pp} - 5%","l");
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw();

}

