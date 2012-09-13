#include "TCanvas.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH1D.h"
#include "TFile.h"
#include "TPad.h"
#include "TObject.h"
#include "TMath.h"
#include "TLatex.h"
#include "THStack.h"
#include "TLegend.h"
#include "TLine.h"
#include "TFrame.h"
#include "TString.h"
#include <TRandom.h>

#include <sstream>
#include <iostream>
#include <fstream>

#include "../Include/MitStyleRemix.hh"  // style settings for drawing
#include "../Include/DYTools.hh"
#include "../Include/DYToolsUI.hh"
#include "../Include/UnfoldingTools.hh"
#include "../Include/MyTools.hh"        // miscellaneous helper functions
#include "../Include/TriggerSelection.hh"
#include "../Include/plotFunctions.hh"
#include "../Include/ComparisonPlot.hh"

#include "../Include/classXSect.h"

// -----------------------------------------------------

// constants

DYTools::TCrossSectionKind_t kind=DYTools::_cs_preFsrNorm;


int c_plotTheoryWithFEWZ=1*kBlack;
int c_plotTheoryWithFEWZfine=1*kBlack;
int c_plotTheoryNoFEWZfine=1*(48+1);
int c_plotTheoryNoFEWZ=1*47;
int c_plotTheory2011=1*(kBlue+1);

int testTheoryPredictions=0;

// Forward declarations

void plotXsec_data_vs_theory1D_work(TString dataFName="", TString dataLabel="", std::vector<TString> *extra_dataFNamesV=NULL, std::vector<TString> *extra_labelsV=NULL);

void PlotTheoryPredictions(ComparisonPlot_t &compPlot, const TriggerSelection &triggers, DYTools::TCrossSectionKind_t theKind, int iMassBin, std::vector<TH1F*> &theories, int testFlag);

// -----------------------------------------------------

// Main function

const int testFsrUnfolding=0;
//int useFEWZ=0;
TCanvas *canvas=NULL;

void plotXsec_study_data_vs_theory1D(){
  //TGaxis ::fgMaxDigits=3;
  TGaxis::SetMaxDigits(3);

  TString dataFName="";
  TString dataLabel="data";
  dataFName="../root_files/DY_m10+pr+a05+o03+pr_4680pb/xSecDET_results_1D__fsrUnf.root";
  
  dataFName="../root_files/date_20120825-et6-eta5/xSec_results_1D__fsrUnf.root";
  //dataFName="../root_files/date_20120825-et6-eta5/xSecDET_results_1D__fsrUnf.root";
  dataFName="../root_files/date_20120825-et7altB-eta8altNegs/xSec_results_1D__fsrUnf.root";
  //dataFName="../root_files/date_20120825/xSecDET_results_2D_Full2011_hltEffOld.root";
  if (testTheoryPredictions!=0) dataFName="";
  std::vector<TString> dataFNamesV, labelsV;

  if (0) {
    dataFNamesV.push_back("../root_files/date_20120825-et6-eta5/xSec_results_1D__fsrUnf.root");
    labelsV.push_back("et6-eta5");
    dataFNamesV.push_back("../root_files/date_20120825-et6altB-eta4test/xSec_results_1D__fsrUnf.root");
    labelsV.push_back("et6altB-eta4test");
    dataFNamesV.push_back("../root_files/date_20120825-et7altB-eta8altNegs/xSec_results_1D__fsrUnf.root");
    labelsV.push_back("et7altB-eta8altNegs");
  }

  // Test FSR unfolding
  if (0 && testFsrUnfolding) {
    const char *preFsrDet="xSecDET_results_2D_yieldsMcPreFsrGenDET_orig_fsrUnf.root";
    const char *postFsrDet="xSecDET_results_2D_yieldsMcPostFsrGenDET_orig_fsrUnf.root";
    const char *postFsrDet_unfolded="xSecDET_results_2D_yieldsMcPreFsrGenDET_from_unf_test_fsrUnf.root";
    const char *postFsr_unfolded="xSec_results_2D_yieldsMcPreFsrGen_from_unf_test_fsrUnf.root";
    const char *postFsrDet_unfoldedTest="xSecDET_results_2D_yieldsMcPreFsrGenDET_from_unf_test_fsrUnf.root";
    const char *preFsrDet_fromFnc="xSecDET_results_2D_yieldsMcPreFsrGenDET_from_fnc_fsrUnf.root";
    if (1) {
      const char* path="../root_files/date_20120825-et7altB-eta8alt/";
      const char *label="et7altB-eta8alt";
      //dataFNamesV.push_back(Form("%s%s",path,postFsrDet));
      //labelsV.push_back(Form("%s (postFsrDet)",label));
      //dataFNamesV.push_back(Form("%s%s",path,preFsrDet));
      //labelsV.push_back(Form("%s (preFsrDet)",label));
      //dataFNamesV.push_back(Form("%s%s",path,preFsrDet_fromFnc));
      //labelsV.push_back(Form("%s (fnc.unf.postFsrDet)",label));
      //dataFNamesV.push_back(Form("%s%s",path,postFsrDet_unfolded));
      //labelsV.push_back(Form("%s (unf.postFsrDet)",label));
      //dataFNamesV.push_back(Form("%s%s",path,postFsrDet_unfoldedTest));
      //labelsV.push_back(Form("%s (unf.postFsrDet.tst)",label));

      //dataFName=Form("%s%s",path,postFsrDet_unfolded);
      //dataLabel="unf.postFsrDet";
      dataFName=Form("%s%s",path,postFsr_unfolded);
      dataLabel="unf.postFsr";
    }
  }

  plotXsec_data_vs_theory1D_work(dataFName,dataLabel,
				    &dataFNamesV,&labelsV);

  if (canvas && 0) {
    TString canvFName=Form("fig-xSec-2D-data_%s",CrossSectionKindName(kind).Data());
    std::cout << "canvFName=<" << canvFName << "\n";
    SaveCanvas(canvas,canvFName);
  }

}

// -----------------------------------------------------

void plotXsec_data_vs_theory1D_work(TString dataFName, TString dataLabel, std::vector<TString> *dataFNamesV, std::vector<TString> *labelsV){

  const int testMassBin=0;
  int iMassBin=-1;

  HERE("entered *_work");

  TriggerSelection triggers("Full2011_hltEffOld",true,0);

//-----------------------------------------------------------------
// Read data
//-----------------------------------------------------------------

  // Create a canvas with pads

  TString compPlotYAxisLabel="y axis";
  switch (kind) {
  case DYTools::_cs_preFsrDetNorm: 
    compPlotYAxisLabel="normalized cross section (DET)";
    break;
  case DYTools::_cs_preFsrDet:
    compPlotYAxisLabel="counts";
    break;
  case DYTools::_cs_preFsr:
    compPlotYAxisLabel="counts in full space";
    break;
  case DYTools::_cs_preFsrNorm: 
    compPlotYAxisLabel="normalized cross section";
    break;
  default:
    compPlotYAxisLabel="y axis";
  }

  ComparisonPlot_t compPlot(ComparisonPlot_t::_ratioRel,
			    "check","",
			    (DYTools::study2D) ? "rapidity |y|" : "m_{ee} [GeV]",
			    compPlotYAxisLabel,"(dt-th)/th");
  compPlot.ErrorsOnRatios(0);
  compPlot.SetLogx(1);
  compPlot.SetLogy(1);
  //compPlot.SetYRange(5e-11,1.);
  //compPlot.SetRatioYRange(-0.06,0.06);
  compPlot.SetRatioNdivisions(805);
  compPlot.SetPrintRatioNames(1);

  const int cWidth=600;
  const int cHeight=700;
  if (!canvas) {
    TString canvasName=TString("cXsec_") + DYTools::analysisTag;
    canvas=MakeCanvas(canvasName, canvasName, cWidth,cHeight);
    canvas->Divide(1,2);
    compPlot.PreparePads(canvas);
  }

  //kind=DYTools::_cs_preFsr;
  TH1F* hData=NULL;
  if (dataFName.Length()) {
    HERE(dataFName);
    hData= readData(triggers, kind,iMassBin, dataFName);
    //removeError(hData);

    hData->GetXaxis()->SetLabelOffset(1.5);
  }

  std::vector<TH1F*> hDataV;
  if (dataFNamesV) {
    HERE("load dataFNamesV");
    hDataV.reserve(dataFNamesV->size());
    for (unsigned int i=0; i<dataFNamesV->size(); ++i) {
      TH1F* histo=readData(triggers,kind,iMassBin,(*dataFNamesV)[i]);
      TString  hName= histo->GetName();
      histo->SetName( hName + (*labelsV)[i] );
      hDataV.push_back(histo);
      if (testMassBin) printHisto(histo);
    }
  }

  std::vector<TH1F*> theories;
 
  PlotTheoryPredictions(compPlot,triggers,kind,iMassBin,theories, testTheoryPredictions);

  if (dataFName.Length()) {
    printf("Data cross section\n");
    printHisto(std::cout,hData);
  }

  if (dataFName.Length()) {
    hData->SetMarkerSize(0.8);
    compPlot.AddHist1D(hData,dataLabel,"L",kBlack,1,0,1);
  }

  for (unsigned int i=0; i<hDataV.size(); ++i) {
    hDataV[i]->SetMarkerStyle(21+i);
    hDataV[i]->SetMarkerSize(0.8);
    compPlot.AddHist1D(hDataV[i],(*labelsV)[i],"L",getStandardColor(i),1,0,1);
  }

  int subpad1=1;
  compPlot.Draw(canvas,false,"png");
  canvas->cd(subpad1);
  compPlot.AddTextCMSPreliminary();
  //compPlot.AddTextLumi(0.55,0.27);

  canvas->Update();
  canvas->cd();
  
  return;
}

// ------------------------------------------------------------

void PlotTheoryPredictions(ComparisonPlot_t &compPlot, const TriggerSelection &triggers, DYTools::TCrossSectionKind_t theKind, int iMassBin, std::vector<TH1F*> &theories, int testFlag) {

  int plotTheoryWithFEWZ=c_plotTheoryWithFEWZ;
  int plotTheoryWithFEWZfine=c_plotTheoryWithFEWZfine;
  int plotTheoryNoFEWZfine=c_plotTheoryNoFEWZfine;
  int plotTheoryNoFEWZ=c_plotTheoryNoFEWZ;
  int plotTheory2011=c_plotTheory2011;
  int plotTheoryWithFEWZ_2MCfiles=kOrange+1;
  int plotTheoryWithFEWZ_2MCfilesFine=kOrange+1;
  int plotTheoryWithFEWZ_2MCfiles_debug=0;

  int extraFlag=0;
  TH1F *theory2011=NULL;
  TH1F *theoryFEWZ=NULL;
  TH1F *theoryWithFEWZfine=NULL;
  TH1F *theoryNoFEWZ=NULL;
  TH1F *theoryNoFEWZfine=NULL;
  TH1F *theoryWithFEWZ_2MCfiles=NULL;
  TH1F *theoryWithFEWZ_2MCfilesFine=NULL;
  TH1F *theoryWithFEWZ_2MCfiles_debug=NULL;
  //TH1F *theoryWithFEWZ_2MCfilesFine_debug=NULL;

  switch(testFlag) {
  case 0 :
    plotTheoryWithFEWZfine *=0;
    plotTheoryWithFEWZ *=1;
    plotTheoryNoFEWZfine *=0;
    plotTheoryNoFEWZ *=1;
    plotTheory2011 *=1;
    plotTheoryWithFEWZ_2MCfiles *=0;
    plotTheoryWithFEWZ_2MCfilesFine *=0;
    plotTheoryWithFEWZ_2MCfiles_debug *=0;
    break;
  case 1: 
    plotTheoryWithFEWZ_2MCfiles=0;
    plotTheoryWithFEWZ_2MCfilesFine=0;
    break;
  case 2:
    // no modifications;
    break;
  case 3: 
    plotTheoryNoFEWZ=0;
    plotTheoryNoFEWZfine=0;
    plotTheory2011=0;
    break;
  case 4: 
    plotTheoryNoFEWZ=0;
    plotTheoryNoFEWZfine=0;
    plotTheory2011=0;
    plotTheoryWithFEWZ_2MCfiles *= 1;
    plotTheoryWithFEWZ_2MCfilesFine=0;
    plotTheoryWithFEWZ_2MCfiles_debug=kOrange+5;
    break;
  }


  // ------------- 2MC files
  if (plotTheoryWithFEWZ_2MCfilesFine) {
    theoryWithFEWZ_2MCfilesFine=readTh(triggers,theKind,iMassBin,1,extraFlag_2MCfilesFineGrid);
    AppendToHistoName(theoryWithFEWZ_2MCfilesFine,"_withFEWZ_2MCfilesFineGrid");
    //printHisto(theoryWithFEWZ_2MCfilesFine,1);
    int showLabel=(plotTheoryWithFEWZ_2MCfilesFine==0) ? 1:-1;
    compPlot.AddHist1D(theoryWithFEWZ_2MCfilesFine,"powheg (2MC files)","L",plotTheoryWithFEWZ_2MCfilesFine,1,0,showLabel);
    compPlot.SkipInRatioPlots(theoryWithFEWZ_2MCfilesFine);
  }

  if (plotTheoryWithFEWZ_2MCfiles) {
    theoryWithFEWZ_2MCfiles=readTh(triggers,theKind,iMassBin,1,extraFlag_2MCfiles);
    AppendToHistoName(theoryWithFEWZ_2MCfiles,"_withFEWZ_2MCfiles");
    //printHisto(theoryFEWZ);
    compPlot.AddHist1D(theoryWithFEWZ_2MCfiles,"powheg + FEWZ (2MC files)","P",plotTheoryWithFEWZ_2MCfiles,1,0,1);
  }

  if (plotTheoryWithFEWZ_2MCfiles_debug) {
    theoryWithFEWZ_2MCfiles_debug=readTh(triggers,theKind,iMassBin,1,extraFlag_2MCfiles_debug);
    AppendToHistoName(theoryWithFEWZ_2MCfiles_debug,"_withFEWZ_2MCfiles_debug");
    //printHisto(theoryFEWZ);
    compPlot.AddHist1D(theoryWithFEWZ_2MCfiles_debug,"powheg + FEWZ (2MC files, small count)","P",plotTheoryWithFEWZ_2MCfiles_debug,1,0,1);
  }

  // Theory with FEWZ
  if (plotTheoryWithFEWZfine) {
    theoryWithFEWZfine=readTh(triggers,theKind,iMassBin,1,extraFlag_fineGrid);
    AppendToHistoName(theoryWithFEWZfine,"_withFEWZ_fineGrid");
    //printHisto(theoryWithFEWZfine,1);
    int showLabel=(plotTheoryWithFEWZ==0) ? 1:-1;
    compPlot.AddHist1D(theoryWithFEWZfine,"powheg","L",plotTheoryWithFEWZfine,1,0,showLabel);
    compPlot.SkipInRatioPlots(theoryWithFEWZfine);
  }

  if (plotTheoryWithFEWZ) {
    theoryFEWZ=readTh(triggers,theKind,iMassBin,1,extraFlag);
    AppendToHistoName(theoryFEWZ,"_withFEWZ");
    //printHisto(theoryFEWZ);
    compPlot.AddHist1D(theoryFEWZ,"powheg + FEWZ","P",plotTheoryWithFEWZ,1,0,1);
  }

  if (plotTheoryNoFEWZfine) {
    theoryNoFEWZfine=readTh(triggers,theKind,iMassBin,0,extraFlag_fineGrid);
    AppendToHistoName(theoryNoFEWZfine,"_noFEWZ_fineGrid");
    //printHisto(theoryNoFEWZfine,1);
    removeError(theoryNoFEWZfine);
    int showLabel=(plotTheoryNoFEWZ==0) ? 1:-1;
    compPlot.AddHist1D(theoryNoFEWZfine,"powheg","L",plotTheoryNoFEWZfine,1,0,showLabel);
    compPlot.SkipInRatioPlots(theoryNoFEWZfine);
  }

  if (plotTheoryNoFEWZ) {
    theoryNoFEWZ=readTh(triggers,theKind,iMassBin,0,extraFlag);
    AppendToHistoName(theoryNoFEWZ,"_noFEWZ");
    //printHisto(theoryNoFEWZ);
    compPlot.AddHist1D(theoryNoFEWZ,"powheg","P",plotTheoryNoFEWZ,1,0,1);
    if (!plotTheory2011) compPlot.SetRefIdx(theoryNoFEWZ);
  }
  if (plotTheory2011) {
    if ((theKind!=DYTools::_cs_preFsr) && (theKind!=DYTools::_cs_preFsrNorm)) {
      std::cout << "\n\t theory2011 is only preFSR absolute or normalized\n";
    }
    else {
      theory2011=readTh2011(theKind);
      AppendToHistoName(theory2011,"_powheg2011");
      //printHisto(theory2011);
      removeError(theory2011);
      compPlot.AddHist1D(theory2011,"theory 2011","L",plotTheory2011,1,0,1);
      compPlot.SkipInRatioPlots(theory2011);
    }
  }

  if (theoryNoFEWZ) compPlot.SetRefIdx(theoryNoFEWZ);
  else if (theoryFEWZ) compPlot.SetRefIdx(theoryFEWZ);

  if (theory2011) theories.push_back(theory2011);
  if (theoryFEWZ) theories.push_back(theoryFEWZ);
  if (theoryWithFEWZfine) theories.push_back(theoryWithFEWZfine);
  if (theoryNoFEWZ) theories.push_back(theoryNoFEWZ);
  if (theoryNoFEWZfine) theories.push_back(theoryNoFEWZfine);
  if (theoryWithFEWZ_2MCfiles) theories.push_back(theoryWithFEWZ_2MCfiles);
  if (theoryWithFEWZ_2MCfilesFine) theories.push_back(theoryWithFEWZ_2MCfilesFine);
  if (theoryWithFEWZ_2MCfiles_debug) theories.push_back(theoryWithFEWZ_2MCfiles_debug);
}


// ------------------------------------------------------------


// ------------------------------------------------------------
