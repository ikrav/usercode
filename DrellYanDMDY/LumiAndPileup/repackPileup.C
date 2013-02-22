//
// This script reads in the pile-up histograms created by the
// recommended CMS procedure (which comes as 50 bins from 0 to 50)
// and repackages them into histograms with binning of our choosing
// and with the desired names for the DrellYanDMDY package use.
//
// The file names need to be edited by hand in the if statements below

#include <TFile.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TString.h>
#include <assert.h>

enum SAMPLE {
  GENERAL_DATA = 0,
  GENERAL_MC   = 1,
  TNP_RECO     = 2,
  TNP_ID       = 3,
  TNP_HLT      = 4
};

enum SYSTEMATICS {
  DEFAULT = 0,
  PLUS_5_PERCENT = 1,
  MINUS_5_PERCENT = 2
};

void repackPileup(int sample, int systematics_flag = DEFAULT){

  int nbins;
  double lowEdge, upEdge;
  int nbinsnew = 46;
  double lowEdgeNew = -0.5;
  double upEdgeNew = 45.5;

  TString fname;
  TString fnameNew;
  TString outputHistName;
  if( sample == GENERAL_DATA ){
    fname    = "../root_files/pileup/tmp/dataPileupHildreth_full2011_20130214";
    fnameNew = "../root_files/pileup/dataPileupHildreth_full2011_20130214_repacked";
    outputHistName = "pileup_lumibased_data";
  }else if( sample == GENERAL_MC ){
    fname    = "../root_files/pileup/tmp/mcPileupHildreth_full2011_20130214";
    fnameNew = "../root_files/pileup/mcPileupHildreth_full2011_20130214_repacked";
    outputHistName = "pileup_simulevel_mc";
  }else if( sample == TNP_RECO ){
    fname    = "../root_files/pileup/tmp/dataPileupHildreth_full2011_TnP_RECO_20130214";
    fnameNew = "../root_files/pileup/dataPileupHildreth_full2011_TnP_RECO_20130214_repacked";
    outputHistName = "pileup_lumibased_data";
  }else if( sample == TNP_ID ){
    fname    = "../root_files/pileup/tmp/dataPileupHildreth_full2011_TnP_ID_20130214";
    fnameNew = "../root_files/pileup/dataPileupHildreth_full2011_TnP_ID_20130214_repacked";
    outputHistName = "pileup_lumibased_data";
  }else if( sample == TNP_HLT ){
    fname    = "../root_files/pileup/tmp/dataPileupHildreth_full2011_TnP_HLT_20130214";
    fnameNew = "../root_files/pileup/dataPileupHildreth_full2011_TnP_HLT_20130214_repacked";
    outputHistName = "pileup_lumibased_data";
  }
  if( systematics_flag == PLUS_5_PERCENT ){
    fname += "_plus5percent";
    fnameNew += "_plus5percent";
  }else if( systematics_flag == MINUS_5_PERCENT ){
    fname += "_minus5percent";
    fnameNew += "_minus5percent";
  }else{
    // no modifer for systematics_flag == DEFAULT
  }
  fname += ".root";
  fnameNew += ".root";

 // Read input histogram file and find the histogram
  TFile f1(fname);
  if( !f1.IsOpen())
    assert(0);

  TH1D *hist1 = (TH1D*)f1.Get("pileup");
  if(hist1 == 0)
    assert(0);
  
  // Make sure the binning is expected:
  nbins = hist1->GetNbinsX();
  lowEdge = hist1->GetXaxis()->GetBinLowEdge(1);
  upEdge  = hist1->GetXaxis()->GetBinUpEdge(nbins);
  if( ! (nbins == 50 && lowEdge == 0.0 && upEdge == 50.0) )
    assert(0);

  // Create and fill a new histogram
  TH1F *hist1new = new TH1F(outputHistName, "", nbinsnew, lowEdgeNew, upEdgeNew);
  hist1new->SetDirectory(0);
  for(int i=1; i<=nbinsnew; i++){
    hist1new->SetBinContent(i, hist1->GetBinContent(i));
  }
  f1.Close();
  
  // Write the result out
  TFile f1new(fnameNew, "recreate");
  hist1new->Write();
  f1new.Close();

  return;
}

