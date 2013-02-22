#include "TString.h"
#include "TFile.h"
#include "TH1D.h"

enum MODE {
  DO_ID  = 0,
  DO_HLT = 1
};

enum SYSTEMATICS {
  DEFAULT = 0,
  PLUS_5_PERCENT = 1,
  MINUS_5_PERCENT = 2
};

  
void combineTwoPileups(int mode, int syst = DEFAULT){
  
  TString f1, f2, f3;
  double lumi1, lumi2;
  
  if(mode == DO_ID){
    f1 = "../root_files/pileup/tmp/MyDataPileupHistogram_sc";
    f2 = "../root_files/pileup/tmp/MyDataPileupHistogram_ele";
    f3 = "../root_files/pileup/tmp/dataPileupHildreth_full2011_TnP_ID_20130214";
    lumi1 = 420.306;
    lumi2 = 3780.07;
  }else if( mode == DO_HLT ){
    f1 = "../root_files/pileup/tmp/MyDataPileupHistogram_sc";
    f2 = "../root_files/pileup/tmp/MyDataPileupHistogram_ele";
    f3 = "../root_files/pileup/tmp/dataPileupHildreth_full2011_TnP_HLT_20130214";
    lumi1 = 420.306;
    lumi2 = 955.21;
  }else 
    assert (0);

  if(syst == PLUS_5_PERCENT ){
    f1 += "_plus5percent";
    f2 += "_plus5percent";
    f3 += "_plus5percent";
  }else if(syst == MINUS_5_PERCENT ){
    f1 += "_minus5percent";
    f2 += "_minus5percent";
    f3 += "_minus5percent";
  }else{
    // for syst == DEFAULT do nothing
  }
  f1 += ".root";
  f2 += ".root";
  f3 += ".root";
  
    
  TFile file1(f1);
  if( !file1.IsOpen() )
    assert(0);
  TH1D *pileup1 = (TH1D*)file1.Get("pileup");
  pileup1->Sumw2();

  TFile file2(f2);
  if( !file2.IsOpen() )
    assert(0);
  TH1D *pileup2 = (TH1D*)file2.Get("pileup");
  pileup2->Sumw2();

  double w1 = lumi1/(lumi1 + lumi2);
  double w2 = lumi2/(lumi1 + lumi2);

  pileup1->Scale(1.0/pileup1->GetSumOfWeights());
  pileup2->Scale(1.0/pileup2->GetSumOfWeights());

  TH1D *pileup3 = (TH1D*)pileup1->Clone("pileup");
  pileup3->Reset();
  pileup3->Add(pileup1, pileup2, w1, w2);
  
  pileup3->Print();

  TFile file3(f3, "recreate");
  file3.cd();
  pileup3->Write();
  file3.Close();

}
