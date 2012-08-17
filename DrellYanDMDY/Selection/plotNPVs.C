#include <TROOT.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>

#include "../Include/CPlot.hh"
#include "../Include/MitStyleRemix.hh"
#include "../Include/InputFileMgr.hh"
#include "../Include/DYTools.hh"
#include "../Include/MyTools.hh"
#include "../Include/PUReweight.hh"


void plotNPVs(const TString inpFile) {
  InitialInputMgr_t inpMgr;
  if (!inpMgr.Load(inpFile)) return;

  PUReweight_t puMgr;
  if (!puMgr.setDefaultFile(inpMgr.dirTag(),DYTools::analysisTag_USER,0)) return;
  

  TCanvas *c=new TCanvas("nPVs","nPVs",800,600);// MakeCanvas("nPVs","nPVs",800,600);

  if (!puMgr.setReference("hNGoodPV_data") ||
      !puMgr.setActiveSample("hNGoodPV_zee")) return;
  TH1F *hData=(TH1F*)puMgr.getHRef()->Clone("hpv_data");
  hData->SetDirectory(0);
  hData->Scale(1/hData->Integral());
  std::cout << "hData->Integral=" << hData->Integral() << "\n";
  TH1F *hZee=(TH1F*)puMgr.getHActive()->Clone("hpv_zee");
  hZee->SetDirectory(0);
  hZee->Scale(1/hZee->Integral());
  hZee->SetMarkerStyle(24);

  CPlot *cp=new CPlot("cp","","nGoodPVs","count (a.u.)");
  cp->AddHist1D(hData,"data","lpe",kBlue);
  cp->AddHist1D(hZee,"MC signal","lpe same",kBlack);
  cp->SetYRange(0,0.2);

  if (1) {
    hData->SetTitle("");
    hData->GetXaxis()->SetTitle("nGoodPVs");
    hData->GetYaxis()->SetTitle("probability");
    hData->GetYaxis()->SetTitleOffset(1.2);
    hData->GetXaxis()->SetTitleOffset(1.);
    hData->SetMaximum(0.2);
    hData->Draw("lpe"); // set axes
    hZee->Draw("lpe same");

    TLegend *fLeg = new TLegend(0.6,0.64,0.93,0.9);
    fLeg->AddEntry(hData,"data","PL");
    fLeg->AddEntry(hZee,"MC signal","PL");
    fLeg->Draw();

  }
  else {
    cp->Draw(c,false,"png");
  }

  c->Update();
  SaveCanvas(c,"NPVs");

  return;
}
