#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"

void prepareMCPileupReference(){

  TFile f("/home/hep/ikrav/work/ntuples/DrellYan_8TeV_53X_local/s12-zeem20-v7a_tight-loose_skim.root");
  TTree *tree = (TTree*)f.Get("Events");

  TH1F *pileup = new TH1F("pileup","",100,0,100);

  tree->Draw("Info.nPU>>pileup");

  pileup->SetDirectory(0);
  pileup->Print();

  TFile fout("mcPileupHildreth_full2012.root", "recreate");
  fout.cd();
  pileup->Write();
  fout.Close();
}
