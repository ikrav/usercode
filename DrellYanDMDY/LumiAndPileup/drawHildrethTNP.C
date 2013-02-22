//
// This script draws an overlay of the pile-up distributions
// used in the analysis for TnP samples: the distribution for MC, 
// the three different distributions for the data (since the
// required triggers and prescales are different for different 
// efficiency types in tag and probe

void drawHildrethTNP(){

  TH1F *h1;
  
  TFile f1("../root_files/pileup/mcPileupHildreth_full2011_20121110_repacked.root");
  h1 = (TH1F*)f1.Get("pileup_simulevel_mc");
  TH1F *hmc = (TH1F*)h1->Clone("hmc");
  hmc->SetDirectory(0);
  f1.Close();

  TFile f2("../root_files/pileup/dataPileupHildreth_full2011_TnP_RECO_20121118_repacked_default.root");
  h1 = (TH1F*)f2.Get("pileup_lumibased_data");
  TH1F *hdataReco = (TH1F*)h1->Clone("hdataReco");
  hdataReco->SetDirectory(0);
  f2.Close();

  TFile f3("../root_files/pileup/dataPileupHildreth_full2011_TnP_ID_20121118_repacked_default.root");
  h1 = (TH1F*)f3.Get("pileup_lumibased_data");
  TH1F *hdataId = (TH1F*)h1->Clone("hdataId");
  hdataId->SetDirectory(0);
  f3.Close();

  TFile f4("../root_files/pileup/dataPileupHildreth_full2011_TnP_HLT_20121118_repacked_default.root");
  h1 = (TH1F*)f4.Get("pileup_lumibased_data");
  TH1F *hdataHlt = (TH1F*)h1->Clone("hdataHlt");
  hdataHlt->SetDirectory(0);
  f4.Close();

  hdataReco->Scale(1.0/hdataReco->GetSumOfWeights());
  hdataId  ->Scale(1.0/hdataId  ->GetSumOfWeights());
  hdataHlt ->Scale(1.0/hdataHlt ->GetSumOfWeights());
  hmc->Scale(1.0/hmc->GetSumOfWeights());

  hdataHlt->GetXaxis()->SetRangeUser(0.0,40.0);

  hdataReco->SetLineColor(kBlue);
  hdataId->SetLineColor(kGreen);

  hmc->SetLineColor(kOrange+7);
  hmc->SetFillColor(kOrange+7);
  hmc->SetFillStyle(1001);
  
  hdataHlt->GetXaxis()->SetTitle("pile-up interactions");
  hdataHlt->GetXaxis()->SetTitleOffset(1.0);

  TCanvas *c = MakeCanvas("c","c");

  hdataHlt->Draw("hist");

  hmc->Draw("hist,same");

  hdataHlt->Draw("same,hist");
  hdataReco->Draw("hist,same");
  hdataId->Draw("hist,same");
  
  TLegend *leg = new TLegend(0.5,0.6,0.9, 0.9);
  leg->AddEntry(hdataReco, "TnP 2011 data for #epsilon_{RECO}","l");
  leg->AddEntry(hdataId  , "TnP 2011 data for #epsilon_{ID}","l");
  leg->AddEntry(hdataHlt, "TnP 2011 data for #epsilon_{HLT}","l");
  leg->AddEntry(hmc, "Fall 201111 TnP MC", "f");
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw();

}

