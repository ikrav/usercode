
#include "../../Include/ElectronEnergyScale.hh"
//#include "../Interface/ElectronEnergyScaleAdv.hh"
#include <TCanvas.h>
#include <TLegend.h>

void eesSmearEventDemo(double eta1=0.0, double eta2=2.1, double mass=83., UInt_t nEvents=1) {
  std::cout << "\nElectronEnergyScale event smear demo\n";
  //ElectronEnergyScale esc("Date20120101_default");
  ElectronEnergyScale esc("UNCORRECTED");
  //ElectronEnergyScale esc;
  //esc.init("FileGauss /home/andriusj/testESF_3EB3EENegs_Gauss.inp",1);
  esc.init("FileBreitWigner /home/andriusj/testESF_6binNegs_BreitWigner.inp",1);
  //esc.init("FileVoigt /home/andriusj/testESF_6binNegs_Voigtian.inp",1);
  esc.print();
  //esc.init(ElectronEnergyScale::UNCORRECTED);

  int seed = 0;
  int randomized=0;
  if (esc.calibrationSetFunctionName().CompareTo("Gauss")==0) {
    randomized=1;
    std::cout << "randomizing smearing width for a Gaussian\n";
    esc.randomizeSmearingWidth(seed);
  }

  int eta1Bin=esc.getEtaBinIdx(eta1);
  int eta2Bin=esc.getEtaBinIdx(eta2);

  TCanvas *c = new TCanvas("c","c",600,600);
  double xmin=60;
  double xmax=120;
  TString hSpikesName="hSpikes  (calibrationSet=";
  hSpikesName.Append(esc.calibrationSetName());
  hSpikesName.Append(")");
  TH1F *hSpikes=new TH1F("hSpikes",hSpikesName.Data(),100,xmin,xmax);
  TH1F *hSpikesShifted=new TH1F("hSpikesShifted","hSpikesShifted",100,xmin,xmax);
  TH1F *hSmeared=new TH1F("hSmeared","hSmeared",100,xmin,xmax);
  TH1F *hSmearedRnd=new TH1F("hSmearedRnd","hSmearedRnd",100,xmin,xmax);

  int ci=0;
  ci=kBlue+2;
  hSpikes->SetLineColor(ci); hSpikes->SetMarkerColor(ci);
  hSpikes->GetXaxis()->SetTitle("mass (GeV)");
  hSpikes->GetYaxis()->SetTitle("event count");
  ci=38;
  hSpikesShifted->SetLineColor(ci); hSpikesShifted->SetMarkerColor(ci);
  ci=43;
  hSmeared->SetLineColor(ci); hSmeared->SetMarkerColor(ci);
  ci=46;
  hSmearedRnd->SetLineColor(ci); hSmearedRnd->SetMarkerColor(ci);

  const double weight=1.;
  for (UInt_t ev=0; ev<nEvents; ++ev) {
    hSpikes->Fill(mass);
    hSpikesShifted->Fill(mass + esc.generateMCSmear(eta1,eta2));
    if (randomized) esc.addSmearedWeightRandomized(hSmearedRnd,eta1Bin,eta2Bin,mass, weight);
  }
  //std::cout << "eta1Bin=" << eta1Bin << ", eta2Bin=" << eta2Bin << "\n";
  esc.addSmearedWeight(hSmeared,eta1Bin,eta2Bin,mass ,nEvents*weight);
  double hsmF=nEvents/hSmeared->Integral(); hSmeared->Scale(hsmF); std::cout << "hSmeared is rescaled by " << hsmF << "\n";

  if (0) {
    hSpikes->Draw();
    hSpikesShifted->Draw("same");
  }
  else {
    hSpikesShifted->Draw();
    hSpikes->Draw("same");
  }
  hSmeared->Draw("LP same");
  if (randomized) hSmearedRnd->Draw("L same");

  std::cout << "hSpikes->Integral()=" << hSpikes->Integral() << ", mean=" << hSpikes->GetMean() << ", rms=" << hSpikes->GetRMS() << "\n";
  std::cout << "hSpikesShifted->Integral()=" << hSpikesShifted->Integral() << ", mean=" << hSpikesShifted->GetMean() << ", rms=" << hSpikesShifted->GetRMS() << "\n";
  std::cout << "hSmeared->Integral()=" << hSmeared->Integral() << ", mean=" << hSmeared->GetMean() << ", rms=" << hSmeared->GetRMS() << "\n";
  std::cout << "hSmearedRnd->Integral()=" << hSmearedRnd->Integral() << ", mean=" << hSmearedRnd->GetMean() << ", rms=" << hSmearedRnd->GetRMS() << "\n";

  TLegend *leg= new TLegend(0.62,0.75, 0.90,0.85);
  leg->AddEntry(hSpikes, "original event","l");
  leg->AddEntry(hSpikesShifted,"randomly shifted event","l");
  leg->AddEntry(hSmeared,"smeared event","lp");
  if (randomized) leg->AddEntry(hSmearedRnd,"randomized smeared event","l");
  leg->Draw();

  c->Update();
  std::cout << "\n";
}
