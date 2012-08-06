#include "../Include/CPlot.hh"
#include "../Include/MitStyleRemix.hh"
#include "../Include/ElectronEnergyScale.hh"
//#include <TCanvas.h>
//#include <TLegend.h>

void plotTwoESFSetsDemo() {
  CPlot::sOutDir = "plots";


  ElectronEnergyScale esf1("Date20120101_default");
  ElectronEnergyScale esf2("Date20120802_default");

  if (!esf1.isInitialized() || 
      !esf2.isInitialized()
      ) {
    esf1.print();
    esf2.print();
    std::cout << "failed to prepare esf's\n";
    return;
  }

  TCanvas *esfCanvas= MakeCanvas("esfCanvas","esfCanvas",1200,600);
  esfCanvas->Divide(2);

  const char *label1="esf20120101";
  const char *label2="esf20120802";
  TH1F *h1Scale=esf1.createScaleHisto("esf20120101");
  TH1F *h1Smear=esf1.createSmearHisto("esf20120101",0);
  TH1F *h2Scale=esf2.createScaleHisto("esf20120802");
  TH1F *h2Smear=esf2.createSmearHisto("esf20120802",0);
  assert(h1Scale); assert(h1Smear);
  assert(h2Scale); assert(h2Smear);

  int ci1=kGreen+2; 
  h1Scale->SetLineColor(ci1); h1Scale->SetMarkerColor(ci1);
  h1Smear->SetLineColor(ci1); h1Smear->SetMarkerColor(ci1);
  h1Scale->SetMarkerStyle(24); h1Smear->SetMarkerStyle(24);
  h1Scale->GetXaxis()->SetTitleOffset(1.100);
  h1Smear->GetXaxis()->SetTitleOffset(1.100);
  h1Scale->GetYaxis()->SetTitleOffset(1.40);
  h1Smear->GetYaxis()->SetTitleOffset(1.40);

  int ci2=46; 
  h2Scale->SetLineColor(ci2); h2Scale->SetMarkerColor(ci2);
  h2Smear->SetLineColor(ci2); h2Smear->SetMarkerColor(ci2);
  h2Scale->SetMarkerStyle(20); h2Smear->SetMarkerStyle(20);


  CPlot cpScale("cpScale"," ", "#eta","data scaling parameter");
  CPlot cpSmear("cpSmear"," ", "#eta","MC smearing pamareter");

  cpScale.AddHist1D(h1Scale, label1, "LPE", ci1, 1,0,1);
  cpScale.AddHist1D(h2Scale, label2, "LPE same", ci2, 1,0,1);
  cpSmear.AddHist1D(h1Smear, label1, "LPE", ci1, 1,0,1);
  cpSmear.AddHist1D(h2Smear, label2, "LPE same", ci2, 1,0,1);

  cpScale.SetYRange(0.95,1.1);  cpScale.TransLegend(-0.25,0);
  cpSmear.SetYRange(0. ,2.5);   cpSmear.TransLegend(-0.25,0);

  cpScale.Draw(esfCanvas,0,"png",1);
  cpSmear.Draw(esfCanvas,0,"png",2);
  esfCanvas->Update();
}
