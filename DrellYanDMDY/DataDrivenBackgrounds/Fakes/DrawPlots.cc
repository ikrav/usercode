#include <TROOT.h>
#include <TH1.h>
#include <TString.h>
#include <TFile.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TVectorT.h>
#include <TMatrixT.h>
#include <boost/shared_ptr.hpp>
#include <iostream>
#include <iomanip> //setprecision
#include <math.h>
#include <vector>

using std::vector;
using std::cout;
using std::setprecision;
using boost::shared_ptr;

const int eTBins(15);
const int nMassBins = 13;
const double eTBinLimits[eTBins+1] =  {0,5,10,15,20,25,30,35,40,50,60,70,90,150,200,300};

int main(int argc, char** argv){

  gROOT->SetStyle("Plain");
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.18);
  //gStyle->SetHistFillColor(4);

  TApplication graphicsPlease("graphicsPlease", &argc, argv); 

  TString inName = TString("plots.root");
  TFile* infile = new TFile(inName);

  TH1F* fRate_reg = (TH1F*)infile->Get("fRate_reg");
  TH1F* fRate_reg2 = (TH1F*)infile->Get("fRate_reg2");
  TH1F* fRate_reg3 = (TH1F*)infile->Get("fRate_reg3");

  TH1F* rapMass20to30 = (TH1F*)infile->Get("rapMass20to30");
  TH1F* rapMass30to45 = (TH1F*)infile->Get("rapMass30to45");
  TH1F* rapMass45to60 = (TH1F*)infile->Get("rapMass45to60");
  TH1F* rapMass60to120 = (TH1F*)infile->Get("rapMass60to120");
  TH1F* rapMass120to200 = (TH1F*)infile->Get("rapMass120to200");
  TH1F* rapMass200to1500 = (TH1F*)infile->Get("rapMass200to1500");

  fRate_reg->SetMarkerColor(kGreen);
  fRate_reg->SetLineColor(kGreen);
  fRate_reg2->SetMarkerColor(kRed);
  fRate_reg2->SetLineColor(kRed);
  fRate_reg3->SetMarkerColor(kBlue);
  fRate_reg3->SetLineColor(kBlue);

  fRate_reg->SetXTitle("electron E_{T} GeV/c");
  fRate_reg->SetYTitle("Fake Rate");

  TCanvas* c1 = new TCanvas("c1","Fake rates v #{eta}",0,0,700,500);
  c1->cd();
    
  fRate_reg->Draw("");
  fRate_reg2->Draw("same");
  fRate_reg3->Draw("same");

  TLegend *leg1 = new TLegend(0.20 ,0.68,0.45,0.85);
  leg1->AddEntry(fRate_reg, "#eta < 1.4442", "lp");
  leg1->AddEntry(fRate_reg2, "1.566< #eta < 2.0", "lp");
  leg1->AddEntry(fRate_reg3, "2.0 < #eta < 2.5", "lp");
  
  leg1->Draw();

  c1->SaveAs("FakeRateOrig.pdf");
  c1->SaveAs("FakeRate.C");

  //Plot obtained fake rate from data

   TH1F* testEGsfF = (TH1F*)infile->Get("testEGsfF"); 
   TCanvas* c2 = new TCanvas("c2","Faked Electrons qcdhist",0,0,700,500);
   c2->cd();
   testEGsfF->SetMarkerColor(kBlue);
   testEGsfF->SetLineColor(kBlue);
   testEGsfF->SetXTitle("e^{+}e^{-} mass GeV/c^{2}");
   testEGsfF->SetYTitle("Number of Events");
   testEGsfF->SetFillColor(4);
   testEGsfF->Draw("e");
   testEGsfF->SetFillColor(4);

   //c2->SaveAs("FakeData.png");
   //c2->SaveAs("FakeData.C");

   //Plot MC estimation

   // inName = TString("bgplotsNOLeadingCut.root");
   inName = TString("bgplots.root");
   TFile* inFile2 = new TFile(inName);
   
   TH1F* BkgEGsfF = (TH1F*) inFile2->Get("BkgEGsfF");
   TH1F* BkgRapMass20to30 = (TH1F*)inFile2->Get("BkgRapMass20to30");
   TH1F* BkgRapMass30to45 = (TH1F*)inFile2->Get("BkgRapMass30to45");
   TH1F* BkgRapMass45to60 = (TH1F*)inFile2->Get("BkgRapMass45to60");
   TH1F* BkgRapMass60to120 = (TH1F*)inFile2->Get("BkgRapMass60to120");
   TH1F* BkgRapMass120to200 = (TH1F*)inFile2->Get("BkgRapMass120to200");
   TH1F* BkgRapMass200to1500 = (TH1F*)inFile2->Get("BkgRapMass200to1500");



   TCanvas* c3 = new TCanvas("c3","Faked Electrons",0,0,700,500);
   c3->cd();
   BkgEGsfF->SetMarkerColor(kRed);
   BkgEGsfF->SetLineColor(kRed);
   BkgEGsfF->SetXTitle("e^{+}e^{-} mass GeV/c^{2}");
   BkgEGsfF->SetYTitle("Number of Events");
   BkgEGsfF->Draw("e");

   cout << "Total number of MC events is: " << BkgEGsfF->Integral(1,40) << "\n";
   c3->SaveAs("BgFakeData.png");


   //Clone signal histos
   TH1F* fakeH = (TH1F*)testEGsfF->Clone("fakeH");

   TH1F* fake20to30 = (TH1F*)rapMass20to30->Clone("fake20to30");
   TH1F* fake30to45 = (TH1F*)rapMass30to45->Clone("fake30to45");
   TH1F* fake45to60 = (TH1F*)rapMass45to60->Clone("fake45to60");
   TH1F* fake60to120 = (TH1F*)rapMass60to120->Clone("fake60to120");
   TH1F* fake120to200 = (TH1F*)rapMass120to200->Clone("fake120to200");
   TH1F* fake200to1500 = (TH1F*)rapMass200to1500->Clone("fake200to1500");


   //lets put these in a vector to help with performing repeat tasks


   vector<TH1F*> fakeHistos;
   fakeHistos.push_back(rapMass20to30);
   fakeHistos.push_back(rapMass30to45);
   fakeHistos.push_back(rapMass45to60);
   fakeHistos.push_back(rapMass60to120);
   fakeHistos.push_back(rapMass120to200);
   fakeHistos.push_back(rapMass200to1500);

   //Plot the difference of data - MC background


   //calculate errors for subtracted histo

   unsigned int numXbins = fakeH->GetNbinsX() + 1;

   for (unsigned int i =1; i < numXbins; ++i){
     double x = BkgEGsfF->GetBinContent(i);
     double y =  testEGsfF->GetBinContent(i);
     double err = sqrt(x+y);
     fakeH->SetBinError(i, err);
     //cout << "x: " << x << " y: " << y << "err: " << err << "\n";
   }

   //Calculate errors for subtracting rapidity plots

   numXbins = fake20to30->GetNbinsX() + 1;
   const unsigned int numXbinsLower = fake200to1500->GetNbinsX() + 1;
   for (unsigned int i =1; i < numXbins; ++i){
     double x = BkgRapMass20to30->GetBinContent(i);
     double y =  rapMass20to30->GetBinContent(i);
     double err = sqrt(x+y);
     fake20to30->SetBinError(i, err);
     x = BkgRapMass30to45->GetBinContent(i);
     y =  rapMass30to45->GetBinContent(i);
     err = sqrt(x+y);
     fake30to45->SetBinError(i, err);
     x = BkgRapMass45to60->GetBinContent(i);
     y =  rapMass45to60->GetBinContent(i);
     err = sqrt(x+y);
     fake45to60->SetBinError(i, err);
     x = BkgRapMass60to120->GetBinContent(i);
     y =  rapMass60to120->GetBinContent(i);
     err = sqrt(x+y);
     fake60to120->SetBinError(i, err);
     x = BkgRapMass120to200->GetBinContent(i);
     y =  rapMass120to200->GetBinContent(i);
     err = sqrt(x+y);
     fake120to200->SetBinError(i, err);
     
     if (numXbins < numXbinsLower){
       x = BkgRapMass200to1500->GetBinContent(i);
       y =  rapMass200to1500->GetBinContent(i);
       err = sqrt(x+y);
       fake200to1500->SetBinError(i, err);
     }
     //cout << "x: " << x << " y: " << y << "err: " << err << "\n";
   }
   

   //apply corrections (subtract backgrounds)  
   fakeH->Add(BkgEGsfF,-1.0);
   fake20to30->Add(BkgRapMass20to30,-1.0);
   fake30to45->Add(BkgRapMass30to45,-1.0);
   fake45to60->Add(BkgRapMass45to60,-1.0);
   fake60to120->Add(BkgRapMass60to120,-1.0);
   fake120to200->Add(BkgRapMass120to200,-1.0);
   fake200to1500->Add(BkgRapMass200to1500,-1.0);

   //pretify plots
   fake20to30->SetMarkerColor(kBlue);
   fake30to45->SetMarkerColor(kBlue);
   fake45to60->SetMarkerColor(kBlue);
   fake60to120->SetMarkerColor(kBlue);
   fake120to200->SetMarkerColor(kBlue);
   fake200to1500->SetMarkerColor(kBlue);

   fake20to30->SetXTitle("e^{+}e^{-} rapidity");
   fake30to45->SetXTitle("e^{+}e^{-} rapidity");
   fake45to60->SetXTitle("e^{+}e^{-} rapidity");
   fake60to120->SetXTitle("e^{+}e^{-} rapidity");
   fake120to200->SetXTitle("e^{+}e^{-} rapidity");
   fake200to1500->SetXTitle("e^{+}e^{-} rapidity");

   /*
    fake20to30->
   fake30to45->
   fake45to60->
   fake60to120->
   fake120to200->
   fake200to1500->
   */

   TCanvas* c4 = new TCanvas("c4","Faked Electrons (background subtracted)",0,0,700,500);
   c4->cd();
   TPad *systpad = new TPad("systpad","systpad",0.1,0.1,0.95,0.995);
   systpad->Draw();
   systpad->SetLogy();
   systpad->SetLogx();
   systpad->cd();
   //fakeH->SetMarkerColor(kGreen);
   fakeH->SetLineColor(kGreen);
   fakeH->SetMarkerStyle(22);
   fakeH->SetXTitle("e^{+}e^{-} mass GeV/c^{2}");
   fakeH->SetYTitle("Number of Events");
   fakeH->Draw("e");
   
   cout << "Total number of final events is: " << fakeH->Integral(1,40) << "\n";

   c4->SaveAs("FakeDist.C");
   c4->SaveAs("FakeDist.png");

   TCanvas* c5_2d = new TCanvas("c5_2d","Faked Electrons rapidity (background subtracted)",0,0,700,500);
   c5_2d->cd();
   TPad *rap2dpad = new TPad("rap2dpad","rap2dpad",0.1,0.1,0.95,0.995);
   rap2dpad->Draw();
   //systpad->SetLogy();
   //systpad->SetLogx();
   rap2dpad->Divide(2,3);

   rap2dpad->cd(1);
   fake20to30->Draw("e");
   rap2dpad->cd(2);
   fake30to45->Draw("e");
   rap2dpad->cd(3);
   fake45to60->Draw("e");
   rap2dpad->cd(4);
   fake60to120->Draw("e");
   rap2dpad->cd(5);
   fake120to200->Draw("e");
   rap2dpad->cd(6);
   fake200to1500->Draw("e");
   c5_2d->SaveAs("rap.png");

   //print out numbers and associated errors

   cout << "Total number of events is: " << fakeH->Integral() << "\n";

   //Plot leading trigger cut fake rate

   TH1F* fRate = (TH1F*) infile->Get("fRate_reg");
   TH1F* fRate2 = (TH1F*) infile->Get("fRate_reg2");
   TH1F* fRate3 = (TH1F*) infile->Get("fRate_reg3");
   TCanvas* c5 = new TCanvas("c5","Faked Rate",0,0,700,500);
   c5->cd();
 
   fRate->SetMarkerColor(kGreen);
   fRate->SetLineColor(kGreen);
   fRate2->SetMarkerColor(kRed);
   fRate2->SetLineColor(kRed);
   fRate3->SetMarkerColor(kBlue);
   fRate3->SetLineColor(kBlue);

   fRate2->SetXTitle("electron E_{T} GeV/c");
   fRate2->SetYTitle("Fake Rate");
   fRate2->Draw("");
   fRate->Draw("same");
   fRate3->Draw("same");
   //fRate->Draw("same");

   TLegend *leg2 = new TLegend(0.55 ,0.68,0.80,0.85);
   //leg2->AddEntry(fRate, "all #eta", "lp");
   leg2->AddEntry(fRate_reg, "#eta < 1.4442", "lp");
   leg2->AddEntry(fRate_reg2, "1.566< #eta < 2.0", "lp");
   leg2->AddEntry(fRate_reg3, "2.0 < #eta < 2.5", "lp");
   leg2->Draw();
   c5->SaveAs("FakeRate.png");

   //==================
   // systematics plots
   //===================

   //c9->SaveAs("eEt.png");

   //Plot fake rate percentage errors
   
   cout << "fake 1: " << BkgEGsfF->Integral(1,40) << "\n";
   //cout << "fake 2: " << qcdHist4_syst->Integral(1,40) << "\n";

   //Fill matrix with 2D info

   TMatrixT<double> fakeBackgroundFromData(6,25);//hardwired, need to change this
   TMatrixT<double> fakeBackgroundFromDataError(6,25);
   TMatrixT<double> fakeBackgroundFromDataErrorSyst(6,25);

   unsigned int numXbins_rap25 = 25;//hardwired need to change this
   for (unsigned int i = 0; i < 5; ++i){
     for (unsigned int j = 0; j < numXbins_rap25; ++j){
       fakeBackgroundFromData[i][j] = fakeHistos.at(i)->GetBinContent(j+1);
       fakeBackgroundFromDataError[i][j] = fakeHistos.at(i)->GetBinError(j+1);
       fakeBackgroundFromDataErrorSyst[i][j] = (0.5 * fakeHistos.at(i)->GetBinContent(j+1));
     }
   }

   unsigned int numXbins_rap10 = 10;//hardwired need to change this
   for (unsigned int j = 0; j < numXbins_rap10; ++j){
     fakeBackgroundFromData[5][j] = fake200to1500->GetBinContent(j+1);
     fakeBackgroundFromDataError[5][j] = fake200to1500->GetBinError(j+1);
     fakeBackgroundFromDataErrorSyst[5][j] = (0.5 * fake200to1500->GetBinContent(j+1));
   }

   /*

   //fill vectors with fake rate info

   TVectorT<double> fakeBackgroundFromData;
   TVectorT<double> fakeBackgroundFromDataError;
   TVectorT<double> fakeBackgroundFromDataErrorSyst;

   int nBinsV(40);//number of bins
 
   //store central value, stat error and syst error of ee estimation
   fakeBackgroundFromData.ResizeTo(nBinsV);
   fakeBackgroundFromDataError.ResizeTo(nBinsV);
   fakeBackgroundFromDataErrorSyst.ResizeTo(nBinsV);

   double error(0);
   double binVal(0);
   for (int i = 0; i < nBinsV; ++i){
     binVal = fakeH->GetBinContent(i+1);
     error = fakeH->GetBinError(i+1);
     if (binVal < 0){//No negative values
       binVal = 0;
       error = 1.15; //use limits on 0 events, ie. 1.15 for 1 sigma
     }
     fakeBackgroundFromData[i] = binVal;     
     if (error == 0.0) error = 1.15;
     fakeBackgroundFromDataError[i] = error;
     fakeBackgroundFromDataErrorSyst[i] = 0.5 * binVal;
   }

   */

   TString outFileName = "fakeBkgDataPoints.root";
   TFile *outFile = new TFile(outFileName,"RECREATE");
   outFile->cd();
   fakeBackgroundFromData.Write("fakeBackgroundFromData");
   fakeBackgroundFromDataError.Write("fakeBackgroundFromDataError");
   fakeBackgroundFromDataErrorSyst.Write("fakeBackgroundFromDataErrorSyst");
   outFile->Write();

   graphicsPlease.Run();

}
