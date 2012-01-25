
#include "../../Include/ElectronEnergyScale.hh"
//#include "../Interface/ElectronEnergyScaleAdv.hh"
#include <TCanvas.h>
#include <TLegend.h>

void eesSmearEventDemo(double eta1=0.0, double eta2=2.1, double mass=83.) {
  std::cout << "\nElectronEnergyScale event smear demo\n";
  ElectronEnergyScale esc("Date20120101_default");
  //esc.init(ElectronEnergyScale::UNCORRECTED);

  int seed = 0;
  esc.randomizeSmearingWidth(seed);

  int eta1Bin=esc.getEtaBinIdx(eta1);
  int eta2Bin=esc.getEtaBinIdx(eta2);

  TCanvas *c = new TCanvas("c","c",600,600);
  TString hSpikesName="hSpikes  (calibrationSet=";
  hSpikesName.Append(esc.calibrationSetName());
  hSpikesName.Append(")");
  TH1F *hSpikes=new TH1F("hSpikes",hSpikesName.Data(),100,60.,120.);
  TH1F *hSpikesShifted=new TH1F("hSpikesShifted","hSpikesShifted",100,60.,120.);
  TH1F *hSmeared=new TH1F("hSmeared","hSmeared",100,60,120);
  TH1F *hSmearedRnd=new TH1F("hSmearedRnd","hSmearedRnd",100,60,120);

  int ci=0;
  ci=kBlue+2;
  hSpikes->SetLineColor(ci); hSpikes->SetMarkerColor(ci);
  hSpikes->GetXaxis()->SetTitle("mass (GeV)");
  hSpikes->GetYaxis()->SetTitle("event count");
  ci=38;
  hSpikesShifted->SetLineColor(31); hSpikesShifted->SetMarkerColor(31);
  ci=43;
  hSmeared->SetLineColor(ci); hSmeared->SetMarkerColor(ci);
  ci=46;
  hSmearedRnd->SetLineColor(ci); hSmearedRnd->SetMarkerColor(ci);

  const double weight=1.;
  hSpikes->Fill(mass);
  hSpikesShifted->Fill(mass + esc.generateMCSmear(eta1,eta2));
  esc.addSmearedWeight(hSmeared,eta1Bin,eta2Bin,mass ,weight);
  esc.addSmearedWeightRandomized(hSmearedRnd,eta1Bin,eta2Bin,mass, weight);

  hSpikes->Draw();
  hSpikesShifted->Draw("same");
  hSmeared->Draw("LP same");
  hSmearedRnd->Draw("L same");

  std::cout << "hSpikes->Integral()=" << hSpikes->Integral() << "\n";
  std::cout << "hSpikesShifted->Integral()=" << hSpikesShifted->Integral() << "\n";
  std::cout << "hSmeared->Integral()=" << hSmeared->Integral() << "\n";
  std::cout << "hSmearedRnd->Integral()=" << hSmearedRnd->Integral() << "\n";

  TLegend *leg= new TLegend(0.62,0.75, 0.90,0.85);
  leg->AddEntry(hSpikes, "original event","l");
  leg->AddEntry(hSpikesShifted,"randomly shifted event","l");
  leg->AddEntry(hSmeared,"smeared event","lp");
  leg->AddEntry(hSmearedRnd,"randomized smeared event","l");
  leg->Draw();

  c->Update();
  std::cout << "\n";
}
