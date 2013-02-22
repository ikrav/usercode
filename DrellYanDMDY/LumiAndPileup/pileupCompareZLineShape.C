// This script draws an overlay of the MC mass spectrum
// prepared with full selection and different pile-up scenarios

const int nbins = 40;
const double limits[nbins+1] = 
  {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 
   81, 86, 91, 96, 101, 106, 110, 115, 120, 126, 133, 141, 
   150, 160, 171, 185, 200, 220, 243, 273, 320, 380, 440, 
   510, 600, 1000, 1500}; // 40 bins

void pileupCompareZLineShape(){

  TString dirDef = "/home/hep/ikrav/releases/another_UserCode_v2/UserCode/ikravchenko/DrellYanDMDY/root_files/yields/DY_m10+pr+a05+o03+pr_4839pb/";
  TString dirPlus = "/home/hep/ikrav/releases/another_UserCode_v2_pileup_syst/UserCode/ikravchenko/DrellYanDMDY/root_files/yields/DY_m10+pr+a05+o03+pr_4839pb/";
  TString dirMinus = "/home/hep/ikrav/releases/another_UserCode_v2_pileup_syst_minus/UserCode/ikravchenko/DrellYanDMDY/root_files/yields/DY_m10+pr+a05+o03+pr_4839pb/";
  TString fileName = "yields1D.root";

  TString fullFileNameDef = dirDef + fileName;
  TString fullFileNamePlus = dirPlus + fileName;
  TString fullFileNameMinus = dirMinus + fileName;

  TFile f1(fullFileNameDef);
  TMatrixD *mdef = (TMatrixD*)f1.Get("yields_zee");

  TFile f2(fullFileNamePlus);
  TMatrixD *mplus = (TMatrixD*)f2.Get("yields_zee");

  TFile f3(fullFileNameMinus);
  TMatrixD *mminus = (TMatrixD*)f3.Get("yields_zee");
  
  TH1F *hdef = new TH1F("hdef","",nbins,limits);
  TH1F *hplus = new TH1F("hplus","",nbins,limits);
  TH1F *hminus = new TH1F("hminus","",nbins,limits);

  hdef->SetDirectory(0);
  hplus->SetDirectory(0);
  hminus->SetDirectory(0);

  hplus->SetLineColor(kBlue);
  hminus->SetLineColor(kRed);

  for(int i=0; i<nbins; i++){
    hdef->SetBinContent(i+1, (*mdef)(i,0));
    hplus->SetBinContent(i+1, (*mplus)(i,0));
    hminus->SetBinContent(i+1, (*mminus)(i,0));    
  }

  TCanvas *c1 = MakeCanvas("c1","c1");
  c1->SetLogx();
  hdef->Draw("hist");
  hplus->Draw("same");
  hminus->Draw("same");
};

