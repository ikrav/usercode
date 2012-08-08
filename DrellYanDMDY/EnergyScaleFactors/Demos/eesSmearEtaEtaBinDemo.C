
#include "../Include/ElectronEnergyScale.hh"
//#include "../ElectronEnergyScaleAdv.hh"
#include <TCanvas.h>
#include <TLegend.h>
#include <TText.h>
#include "fileLocation.hh"


void DrawShadedText(TText *label, Double_t x, Double_t y, const char *text, int color=kBlack) {
  label->SetTextColor(kGray); label->DrawTextNDC(x+0.005,y-0.005,text);
  label->SetTextColor(color); label->DrawTextNDC(x,y,text);
}


// --------------------------------------------

void eesSmearEtaEtaBinDemo(double eta1=0.0, double eta2=2.1) {
  std::cout << "\nElectronEnergyScale smearing demo for (eta1,eta2) bin, defined by (" << eta1 << "," << eta2 << ") values\n";
  ElectronEnergyScale esc("Date20120802_default");
  //ElectronEnergyScale esc("UNCORRECTED");
  //esc.init("FileVoigt ../root_files/constants/EScale/testESF_6binNegs_Voigtian_20120802.inp",1);
  esc.init("FileGauss ../../root_files/constants/EScale/testESF_4EB3EENegs_Gauss_20120802.inp",1);
  //esc.setCalibrationSet(ElectronEnergyScale::UNCORRECTED);
  esc.print();
  if (!esc.isInitialized()) return;

  TString eemFileExp, eemFileMC;
  setEEMFileLocation(eemFileExp,eemFileMC);

  double x_extend=0;
  double xmin=60-x_extend, xmax=120+x_extend;
  std::vector<std::vector<double>*> mcDataAll,expDataAll;
  if (!esc.loadEEMFile( eemFileMC, mcDataAll,xmin,xmax) ||
      !esc.loadEEMFile(eemFileExp,expDataAll,xmin,xmax)) {
    std::cout << "failed to load data\n";
    return;
  }

  int etaEtaIdx=esc.getEtaEtaIdx(eta1,eta2);
  const int eta1Bin=esc.getEtaBinIdx(eta1);
  const int eta2Bin=esc.getEtaBinIdx(eta2);
  const double expEnScaleFactor=sqrt(esc.getEnergyScaleCorrection(eta1)*esc.getEnergyScaleCorrection(eta2));
  //std::cout << "eta1Bin=" << eta1Bin << ", eta2Bin=" << eta2Bin << ", expEnScaleFactor=" << expEnScaleFactor << "\n";
  printf("eta1=%4.2lf (bin #%d), eta2=%4.2lf (bin #%d), expEnScaleFactor=%7.5lf (sqrt(%7.5lf*%7.5lf))\n",eta1,eta1Bin,eta2,eta2Bin,expEnScaleFactor,esc.getEnergyScaleCorrection(eta1),esc.getEnergyScaleCorrection(eta2));
  const std::vector<double> *mcData = mcDataAll[etaEtaIdx];
  const std::vector<double> *expData= expDataAll[etaEtaIdx];

  TCanvas *c = new TCanvas("c","c",800,800);
  c->Divide(2,2);


  TString hName="Energy scale corrected (calibrationSet=";
  hName.Append(esc.calibrationSetName());
  hName.Append(")");
  const double massMin=xmin;
  const double massMax=xmax;
  const int numMassBins=100; //int(massMax-massMin+1e-3);
  const int subdivisionForSmearing=10;
  int ci=kBlack;

  TH1::SetDefaultSumw2 ();
  TH1F *hMCSmeared= new TH1F("hMCSmeared","hMCSmeared",numMassBins,massMin,massMax);
  ci=kBlue; hMCSmeared->SetLineColor(ci); hMCSmeared->SetMarkerColor(ci);
  hMCSmeared->SetMarkerSize(1);
  TH1F *hTemp= new TH1F("hTemp","hTemp", subdivisionForSmearing*numMassBins, massMin, massMax);
  ci=kPink; hTemp->SetLineColor(ci); hTemp->SetMarkerColor(ci);
  TH1F *hMCRaw= new TH1F("hMCRaw","hMCRaw",numMassBins,massMin,massMax);
  ci=38; hMCRaw->SetLineColor(ci); hMCRaw->SetMarkerColor(ci);
  hMCRaw->GetXaxis()->SetTitle("mass [GeV]");
  TH1F *hMCShifted = new TH1F("hMCShifted","hMCShifted",numMassBins,massMin,massMax);
  ci=47; hMCShifted->SetLineColor(ci); hMCShifted->SetMarkerColor(ci);
  TH1F *hExpScaled= new TH1F("hExpScaled","hExpScaled",numMassBins,massMin,massMax);
  ci=46; hExpScaled->SetLineColor(ci); hExpScaled->SetMarkerColor(ci);
  hExpScaled->SetMarkerSize(1);
  hExpScaled->GetXaxis()->SetTitle("mass [GeV]");
  TH1F *hExpRaw= new TH1F("hExpRaw","hExpRaw",numMassBins,massMin,massMax);
  ci=kBlack; hExpRaw->SetLineColor(ci); hExpRaw->SetMarkerColor(ci);
  hExpRaw->SetMarkerSize(1);
  hExpRaw->GetXaxis()->SetTitle("mass [GeV]");

  for (unsigned int i=0; i<mcData->size(); ++i) {
    hMCRaw->Fill((*mcData)[i]);
  }
  for (unsigned int i=0; i<expData->size(); ++i) {
    hExpRaw->Fill((*expData)[i]);
  }
  for (unsigned int i=0; i<mcData->size(); ++i) {
    hTemp->Fill((*mcData)[i]);
  }
  for (unsigned int i=0; i<mcData->size(); ++i) {
    hMCShifted->Fill((*mcData)[i]+esc.generateMCSmear(eta1,eta2));
  }
  double scale=hExpRaw->Integral()/hMCRaw->Integral();
  std::cout << "exp/MC scale=" << scale << "\n";
  hMCRaw->Scale(scale); //hExpRaw->Integral()/hMCRaw->Integral());
  hTemp->Scale(scale); //hExpRaw->Integral()/hTemp->Integral());
  hMCShifted->Scale(scale);

  for (unsigned int i=0; i<expData->size(); ++i) {
    hExpScaled->Fill((*expData)[i] * expEnScaleFactor);
  }
  esc.smearDistribution(hMCSmeared,eta1Bin,eta2Bin,hTemp);
  hMCSmeared->Scale(hExpScaled->GetMaximum()/hMCSmeared->GetMaximum());

  TH1F* hMCSmearedRescaled=(TH1F*)hMCSmeared->Clone("hMCSmearedRescaled");
  //hMCSmearedRescaled->Scale(hExpScaled->GetMaximum()/hMCSmeared->GetMaximum());
  ci=30; hMCSmearedRescaled->SetLineColor(ci); hMCSmearedRescaled->SetMarkerColor(ci);
  hMCSmearedRescaled->SetMarkerSize(1);

  std::cout << "hMCRaw->Integral    =" << hMCRaw->Integral() << "\n";
  std::cout << "hTemp->Integral     =" << hTemp->Integral() << "\n";
  std::cout << "hMCSmeared->Integral=" << hMCSmeared->Integral() << "\n";
  std::cout << "hMCShifted->Integral=" << hMCShifted->Integral() << "\n";
  std::cout << "hExpRaw->Integral   =" << hExpRaw->Integral() << "\n";
  std::cout << "hExpScaled->Integral=" << hExpScaled->Integral() << "\n";
  std::cout << "rescaling hMCSmeared\n";
  hMCSmearedRescaled->Scale(hExpScaled->Integral()/hMCSmeared->Integral());
  std::cout << "hMCSmearedRescaled->Integral=" << hMCSmeared->Integral() << "\n";  

  c->cd(1);
  TLegend *leg1= new TLegend(0.7,0.75, 0.93,0.85);
  hExpRaw->Draw("LPE");
  hMCRaw->Draw("hist same");
  leg1->AddEntry(hExpRaw,"original data","pe");
  leg1->AddEntry(hMCRaw,"original MC","l");
  leg1->Draw();
  TText *label1=new TText();
  DrawShadedText(label1,0.5,0.93,"original distributions");

  c->cd(2);
  TLegend *leg2= new TLegend(0.65,0.75, 0.96,0.85);
  hExpScaled->Draw("LPE");
  leg2->AddEntry(hExpScaled,"en.scaled data","pe");
  //hMCSmeared->Draw("hist same"); leg2->AddEntry(hMCSmeared,"smeared MC","l");
  hMCSmearedRescaled->Draw("hist same"); leg2->AddEntry(hMCSmearedRescaled,"smeared MC, scaled","l");
  leg2->Draw();
  DrawShadedText(label1,0.5,0.93,"modified distributions",46);

  c->cd(3);
  TLegend *leg3= new TLegend(0.7,0.75, 0.93,0.9);
  hMCRaw->Draw("hist");
  hMCSmeared->Draw("lp same");
  //hMCShifted->Draw("hist same"); leg3->AddEntry(hMCShifted,"shifted MC (var)","l");
  leg3->AddEntry(hMCRaw,"original MC","l");
  leg3->AddEntry(hMCSmeared,"smeared MC","lp");
  //  hTemp->Draw("lp same");  leg3->AddEntry(hTemp,"temp MC","pe");
  leg3->Draw();
  DrawShadedText(label1,0.5,0.93,"MC distributions");

  c->cd(4);
  TLegend *leg4= new TLegend(0.7,0.75, 0.93,0.85);
  hExpScaled->Draw("LPE");
  hExpRaw->Draw("hist same");
  leg4->AddEntry(hExpScaled,"scaled data","pe");
  leg4->AddEntry(hExpRaw,"original data","l");
  leg4->Draw();
  DrawShadedText(label1,0.5,0.93,"data distributions");

  c->Update();
  std::cout << "\n\n";
  std::cout << " ! shifted MC is denoted as \"(var)\", since the distribution\n";
  std::cout << " ! will change each time you call the subroutine\n";
  std::cout << "\n";
}
