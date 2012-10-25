#include <TCanvas.h>
#include <TAxis.h>
#include <TGaxis.h>
#include <TH1F.h>
#include <TPad.h>
#include <TFile.h>
#include <TLatex.h>
#include <TString.h>

#include <sstream>
#include <iostream>
#include <fstream>

#include "../Include/MitStyleRemix.hh"  // style settings for drawing
#include "../Include/DYTools.hh"
#include "../Include/UnfoldingTools.hh"
#include "../Include/MyTools.hh"        // miscellaneous helper functions
#include "../Include/TriggerSelection.hh"
#include "../Include/InputFileMgr.hh"
#include "../Include/ComparisonPlot.hh"
#include "../Include/classXSect.h"

//
// User defined constants
//

const int landscape=1;  // page orientation for 2D plots
const int plot_theoryCT10=1; // 2D theory prediction
const int plot_theoryMSTW2008=1; // 2D theory prediction

const int showPoints=0; // whether to show 1D theory prediction in 40 bins

//
// Main function
//

void plotXsec(const TString xsecConfFile, const TString xSecKindString, 
	      const TString crossSectionSet="_fsrUnfGood"){

  // --------------------------------------------------------------
  //           Process input
  // --------------------------------------------------------------

  XSecInputFileMgr_t inpMgr;
  if (!inpMgr.Load(xsecConfFile)) {
    std::cout << "failed to load file <" << xsecConfFile << ">\n";
    return;
  }

  DYTools::TCrossSectionKind_t csKind=DYTools::_cs_None;
  if (xSecKindString.Contains("auto") ||
      xSecKindString.Contains("default")) {
    if (DYTools::study2D==1) csKind=DYTools::_cs_preFsrDetNorm;
    else csKind=DYTools::_cs_preFsrNorm;
  }
  else {
    csKind=DetermineCrossSectionKind(xSecKindString);  // auto-stop on failure
  }

  // --------------------------------------------------------------
  //           Prepare variables
  // --------------------------------------------------------------

  TGaxis::SetMaxDigits(3);
  CPlot::sOutDir = TString("plots_") + DYTools::analysisTag +
    TString("_") + inpMgr.evtEffScaleTag();

  TriggerSelection triggers(inpMgr.triggerSetName(),true,0);
  
  ComparisonPlot_t::TRatioType_t ratio=
      ComparisonPlot_t::_ratioPlain;  //   #1/#2
    //ComparisonPlot_t::_ratioRel;  //   (#1-#2)/#1
  TString ratioYLabel=(ratio==ComparisonPlot_t::_ratioPlain) ?
    "data/theory" : "(th-dt)/th";
  if (DYTools::study2D) {
    ratioYLabel=(ratio==ComparisonPlot_t::_ratioPlain) ?
      "theory/data" : "(dt-th)/dt";
  }

  ComparisonPlot_t compPlot(ratio, "check","",
    (DYTools::study2D) ? "rapidity |y|" : "m_{ee} [GeV]",
			    CrossSectionKindLongName(csKind),
			    ratioYLabel);

  //compPlot.SetRatioNdivisions(805);
  //compPlot.SetPrintValues(1);
  //compPlot.SetPrintRatios(1);
  compPlot.SetPrintRatioNames(1);

  int cWidth=600, cHeight=700;
  if (DYTools::study2D) {
    //double rdy=0.15; compPlot.SetRatioYRange(1-rdy,1+rdy);
    cWidth=400*(2+landscape); cHeight=400*(3-landscape);
  }
  else {
    compPlot.SetLogx(1);
    compPlot.SetLogy(1);
    //double rdy=0.15; compPlot.SetRatioYRange(1-rdy,1+rdy);
    cWidth=600; cHeight=700;
  }

  TString dataFName=TString("../root_files/") + 
    inpMgr.evtEffScaleTag() +
    crossSectionSet +
    TString("/xSec");
  if (!csPreFsr(csKind)) dataFName.Append("PostFsr");
  if (csInDET(csKind)) dataFName.Append("DET");
  dataFName.Append("_results_"); 
  dataFName.Append(DYTools::analysisTag);
  if (!csPreFsr(csKind)) dataFName.Append(triggers.triggerConditionsName());
  dataFName.Append(".root");
  

  // --------------------------------------------------------------
  //      Prepare plots
  // --------------------------------------------------------------

  TString canvasName=TString("cXsec_") + CrossSectionKindName(csKind) + 
    TString("_") + DYTools::analysisTag;
  TCanvas *canvas=MakeCanvas(canvasName,canvasName,cWidth,cHeight);


  std::vector<ComparisonPlot_t*> compPlotsV;
  //std::vector<TH1F*> dataHistos;
  std::vector<TH1F*> histos;

  if (DYTools::study2D) { // 2D

    compPlot.Prepare6Pads(canvas,landscape);

    vector<TLatex*> massLabels;
    prepareMassRanges(massLabels,0.55,0.17, kBlue+2);

    // 

    // data
    for (int i=0; i<6; ++i) {
      const int iM=i+1; // skip 1st mass bin
      ComparisonPlot_t *cp=new ComparisonPlot_t(compPlot,Form("compPlot_%d",i+1),"");
      compPlotsV.push_back(cp);
      std::vector<TH1F*> localHV;

      if ((csKind==DYTools::_cs_preFsrDet) || (csKind==DYTools::_cs_preFsrDetNorm)) {
	if (plot_theoryCT10) {
	  TH1F* hTh=readTh2D_CT10(csKind, iM);
	  hTh->SetMarkerStyle(27);
	  hTh->SetMarkerSize(0.8);
	  cp->AddHist1D(hTh,"CTEQ10W","L",kGreen+1,1,0,1);
	  histos.push_back(hTh);
	  localHV.push_back(hTh);
	}
	if (plot_theoryMSTW2008) {
	  TH1F* hTh=readTh2D_MSTW2008(csKind, iM);
	  hTh->SetMarkerStyle(24);
	  hTh->SetMarkerSize(0.8);
	  cp->AddHist1D(hTh,"MSTW2008","L",kBlue,1,0,1);
	  histos.push_back(hTh);
	  localHV.push_back(hTh);
	}
      }

      // data
      TH1F* hData=readData(triggers,csKind,iM,dataFName);
      histos.push_back(hData);
      localHV.push_back(hData);
      cp->AddHist1D(hData,"data","P",kBlack,1,0,1);
      cp->SetRefIdx(hData);

      for (unsigned int ii=0; ii<localHV.size(); ++ii) {
	localHV[ii]->GetYaxis()->SetTitleOffset(1.2);
      }

      // Plotting and beautifying
      int subpad=-1;
      cp->Draw6(canvas,landscape,i+1,false,"png",&subpad);
      canvas->cd(subpad);
      massLabels[iM]->Draw();
      if (i==2) {
	cp->AddTextCMSPreliminary();
	cp->AddTextLumi(0.55,0.27);
      }
    }

  }
  else { // 1D

    canvas->Divide(1,2);
    compPlot.PreparePads(canvas);

    // theory
    if ((csKind==DYTools::_cs_preFsr) || (csKind==DYTools::_cs_preFsrNorm)) {
      TH1F* hTh=readTh1D_MSTW2008(csKind);
      removeError(hTh);
      compPlot.AddHist1D(hTh,"NNLO, FEWZ+MSTW08","L",kBlue,1,0,0);
      compPlot.SkipInRatioPlots(hTh); // do not consider for ratios
      histos.push_back(hTh);
    }
    // theory: special set to calculate ratios
    if (csKind==DYTools::_cs_preFsrNorm) {
      int showLabel=-1;
      TString draw_opt=(showPoints) ? "P" : "P skip";
      TH1F* hTh=readTh1D_MSTW2008(csKind,"_NNLO","default",_th2011_nnlo);
      hTh->SetMarkerSize(0.8);
      compPlot.AddHist1D(hTh,"NNLO, FEWZ+MSTW08 (40bins)",draw_opt,kBlue,1,0,showLabel);
      compPlot.SetRefIdx(hTh); // use to calculate ratios
      histos.push_back(hTh);
    }

    // data 1D
    TH1F* hData=readData(triggers,csKind,-1, dataFName);
    histos.push_back(hData);
    compPlot.AddHist1D(hData,"data","P",kBlack,1,0,1);

  }


  for (unsigned int i=0; i<histos.size(); i++) {
    histos[i]->SetDirectory(0);
  }

  // --------------------------------------------------------------
  //      Plot
  // --------------------------------------------------------------

  if (DYTools::study2D) {
  }
  else { // 1D
    compPlot.Draw(canvas,false,"png");
    canvas->cd(1);
    compPlot.AddTextCMSPreliminary();
    compPlot.AddTextLumi(0.55,0.27);
  }

  canvas->Update();
  canvas->cd();
  SaveCanvas(canvas,canvas->GetName());
  
  return;
}


