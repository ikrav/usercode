
#include <TFile.h>
#include <TChain.h>
#include <TString.h>
#include <TBenchmark.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TRandom.h>
#include <TVectorD.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TTimeStamp.h>
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <vector>                   // STL vector class
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings

#include <vector>

#include "../Include/CPlot.hh"
#include "../Include/MitStyleRemix.hh"

#include "../Include/DYTools.hh"
#include "../Include/EleIDCuts.hh"

#include "../Include/EWKAnaDefs.hh"
#include "../Include/TGenInfo.hh"
#include "../Include/TEventInfo.hh"
#include "../Include/TDielectron.hh"
#include "../Include/TElectron.hh"
#include "../Include/TVertex.hh"
#include "../Include/TriggerSelection.hh"
#include "../Include/cutFunctions.hh"
#include "../EventScaleFactors/tnpSelectEvents.hh"

#include "../Include/InputFileMgr.hh"
#include "../Include/EventSelector.hh"

using namespace mithep;

const int NEffTypes=3;
// Declaration of arrays into which efficiencies will be loaded
typedef double EffArray_t[NEffTypes][DYTools::nEtBins5][DYTools::nEtaBins5]; // largest storage

template T SQR(const T& x) { return x*x; }


//=== FUNCTION DECLARATIONS ======================================================================================

int fillEfficiencyConstants(const TnPInputFileMgr_t &mcMgr, 
	 const TnPInputFileMgr_t &dataMgr, const TriggerSelection &triggers );
int fillOneEfficiency(const TnPInputFileMgr_t &mgr, const TString filename, 
  UInt_t kindIdx, vector<TMatrixD*> &data, vector<TMatrixD*> &dataErrLo, 
  vector<TMatrixD*> &dataErrHi, vector<TMatrixD*> &dataAvgErr);


Bool_t matchedToGeneratorLevel(const TGenInfo *gen, 
  const TDielectron *dielectron);
int createSelectionFile(const MCInputFileMgr_t &mcMgr, 
    const TString &outSkimFName, TriggerSelection &triggers, int debugMode);

double findEventScaleFactor(const esfSelectEvent_t &data);
double findEventScaleFactor(int kind, const esfSelectEvent_t &data);
double findScaleFactor(int kind, double scEt, double scEta);
double findScaleFactor(int kind, int etBin, int etaBin);

double findEventScaleFactorSmeared(const esfSelectEvent_t &data, int iexp);
double findEventScaleFactorSmeared(int kind, const esfSelectEvent_t &data, 
    int iexp);
double findScaleFactorSmeared(int kind, double scEt, double scEta, int iexp);
double findScaleFactorSmeared(int kind, int etBin, int etaBin, 
    const EffArray_t &dataRndWeight, const EffArray_t &mcRndWeight);

void drawEfficiencies();
void drawEfficiencyGraphs(TGraphErrors *grData, TGraphErrors *grMc,
			  TString yAxisTitle, TString text, TString plotName);
//template<class Graph_t>
//void drawEfficiencyGraphsAsymmErrors(Graph_t *grData, Graph_t *grMc,
//			  TString yAxisTitle, TString text, TString plotName);
void drawScaleFactors();
void drawScaleFactorGraphs(TGraphErrors *gr, TString yAxisTitle, TString text, 
			   TString plotName);

void drawEventScaleFactors(TVectorD scaleRecoV, TVectorD scaleRecoErrV,
			   TVectorD scaleIdV , TVectorD scaleIdErrV ,
			   TVectorD scaleHltV, TVectorD scaleHltErrV,
			   TVectorD scaleV   , TVectorD scaleErrV    );
void drawEventScaleFactorGraphs(TGraphErrors *gr, TString yAxisTitle, 
				TString plotName);

double errOnRatio(double a, double da, double b, double db);

void PrintEffInfoLines(const char *msg, int effKind, int effMethod, 
     int binCount, const double *eff, const double *effErr);

//=== Constants ==========================

const bool savePlots = true;

// File names for efficiency measurements from tag and probe
TString          dirTag;

// Define the method used to obtain the efficiencies
//const int dataRecoEffMethod = FITnFIT;
//const int mcRecoEffMethod   = COUNTnCOUNT;
//const int dataIdEffMethod  = FITnFIT;
//const int mcIdEffMethod    = COUNTnCOUNT;
//const int dataHltEffMethod = COUNTnCOUNT;
//const int mcHltEffMethod   = COUNTnCOUNT;

//int dataEffMethods[NEffTypes], mcEffMethods[NEffTypes];


vector<TMatrixD*> dataEff,mcEff;
vector<TMatrixD*> dataEffErrLo,mcEffErrLo;
vector<TMatrixD*> dataEffErrHi,mcEffErrHi;
vector<TMatrixD*> dataEffAvgErr,mcEffAvgErr;

// Global variables
const int nexp = 100;
EffArray_t ro_Data[nexp], ro_MC[nexp];


int etBinning=-1;
int etBinCount=0;
double *etBinLimits=NULL;

int etaBinning=-1;
int etaBinCount=0;
double *etaBinLimits=NULL;


//=== MAIN MACRO =================================================================================================

void calcEventEff(const TString mcInputFile, const TString tnpDataInputFile, 
    const TString tnpMCInputFile, TString triggerSetString, int selectEvents, 
    int debugMode=0)
{

//  ---------------------------------
//       Preliminary checks
//  ---------------------------------

  // verify whether it was a compilation check
  if (mcInputFile.Contains("_DebugRun_") || 
      triggerSetString.Contains("_DebugRun_")) {
    std::cout << "calcEventEff: _DebugRun_ detected. Terminating the script\n";
    return;
  }

  // fast check
  // Construct the trigger object
  TriggerSelection triggers(triggerSetString, false, 0); 
  assert ( triggers.isDefined() );

  if (debugMode) std::cout << "\n\n\tDEBUG MODE is ON\n\n";

//  ---------------------------------
//         Normal execution
//  ---------------------------------

  gBenchmark->Start("calcEventEff");
  

  CPlot::sOutDir = "plots";

  MCInputFileMgr_t mcMgr;
  TnPInputFileMgr_t tnpDataMgr,tnpMCMgr;

  if (!mcMgr.Load(mcInputFile) ||
      !tnpMCMgr.Load(tnpMCInputFile) ||
      !tnpDataMgr.Load(tnpDataInputFile)) {
    return;
  }
  if (!tnpMCMgr.hasSameBinCounts(tnpDataMgr)) {
    cout << "Files tnpMCInputFile=<" << tnpMCInputFile 
	 << ">, tnpDataInputFile=<" << tnpDataInputFile 
	 << "> have different bin counts:\n";
    cout << "MC   input: " << tnpMCMgr;
    cout << "Data input: " << tnpDataMgr;
    return;
  }
  dirTag=tnpMCMgr.dirTag();

  etBinning=tnpMCMgr.etBinsKind();
  etBinCount=DYTools::getNEtBins(etBinning);
  etBinLimits=DYTools::getEtBinLimits(etBinning);
  
  etaBinning=tnpMCMgr.etaBinsKind();
  etaBinCount=DYTools::getNEtaBins(etaBinning);
  etaBinLimits=DYTools::getEtaBinLimits(etaBinning);
  

  TString selectEventsFName=TString("../root_files/tag_and_probe/") +
    dirTag + TString("/eventSFSelectEvents") + analysisTag + TString(".root");
  
  if (selectEvents && 
      !createSelectionFile(mcMgr, selectEventsFName, triggers, debugMode)) {
    std::cout << "failed to create selection file <" 
	      << selectEventsFName << ">\n";
    return;
  }

  // Read efficiency constants from ROOT files
  // This has to be done AFTER configuration file is parsed
  if (!fillEfficiencyConstants( tnpMCMgr, tnpDataMgr, triggers )) {
    return;
  }

  TH1F *hScale = new TH1F("hScale", "", 150, 0.0, 1.5);
  TH1F *hScaleReco = new TH1F("hScaleReco", "", 150, 0.0, 1.5);
  TH1F *hScaleId  = new TH1F("hScaleId" , "", 150, 0.0, 1.5);
  TH1F *hScaleHlt = new TH1F("hScaleHlt", "", 150, 0.0, 1.5);
  std::vector<TH1F*> hScaleEffV;
  hScaleEffV.reserve(3);
  hScaleEffV.push_back(hScaleReco);
  hScaleEffV.push_back(hScaleId);
  hScaleEffV.push_back(hScaleHlt);

  TH1F *hZpeakEt = new TH1F("hZpeakEt", "", etBinCount, etBinLimits);
  vector<TH1F*> hLeadingEtV;
  vector<TH1F*> hTrailingEtV;
  vector<TH1F*> hElectronEtV;

  vector<TH1F*> hScaleV;
  vector<TH1F*> hScaleRecoV;
  vector<TH1F*> hScaleIdV;
  vector<TH1F*> hScaleHltV;
  vector<vector<TH1F*>*> hScaleEffVV;
  hScaleEffVV.reserve(3);
  hScaleEffVV.push_back(&hScaleRecoV);
  hScaleEffVV.push_back(&hScaleIdV);
  hScaleEffVV.push_back(&hScaleHltV);

  for(int i=0; i<nMassBins; i++){
    TString base = "hScaleV_";
    base += i;
    hScaleV.push_back(new TH1F(base,base,150,0.0,1.5));
    hScaleRecoV.push_back(new TH1F(base+TString("_reco"),
				   base+TString("_reco"),150,0.0,1.5));
    hScaleIdV .push_back(new TH1F(base+TString("_id" ),
				  base+TString("_id" ),150,0.0,1.5));
    hScaleHltV.push_back(new TH1F(base+TString("_hlt"),
				  base+TString("_hlt"),150,0.0,1.5));
    base = "hLeadingEt_";
    base += i;
    hLeadingEtV.push_back(new TH1F(base,base,etBinCount, etBinLimits));
    base = "hTrailingEt_";
    base += i;
    hTrailingEtV.push_back(new TH1F(base,base,etBinCount, etBinLimits));
    base = "hElectronEt_";
    base += i;
    hElectronEtV.push_back(new TH1F(base,base,etBinCount, etBinLimits));
  }

  // Create Gaussian-distributed random offsets for each pseudo-experiment
  /*
  for(int i=0; i<nexp; i++){
    ro_D_B_reco[i] = gRandom->Gaus(0.0,1.0);
    ro_D_E_reco[i] = gRandom->Gaus(0.0,1.0);
    ro_M_B_reco[i] = gRandom->Gaus(0.0,1.0);
    ro_M_E_reco[i] = gRandom->Gaus(0.0,1.0);
    
    ro_D_B_id[i] = gRandom->Gaus(0.0,1.0);
    ro_D_E_id[i] = gRandom->Gaus(0.0,1.0);
    ro_M_B_id[i] = gRandom->Gaus(0.0,1.0);
    ro_M_E_id[i] = gRandom->Gaus(0.0,1.0);
    
    ro_D_B_hlt[i] = gRandom->Gaus(0.0,1.0);
    ro_D_E_hlt[i] = gRandom->Gaus(0.0,1.0);
    ro_M_B_hlt[i] = gRandom->Gaus(0.0,1.0);
    ro_M_E_hlt[i] = gRandom->Gaus(0.0,1.0);
  }
  */
  int debug_pseudo_exps=0;
  for (int i=0; i<nexp; i++) {
    EffArray_t *arr= & ro_Data[i];
    for (int kind=0; kind<NEffTypes; ++kind) {
      for (int iEt=0; iEt<DYTools::nEtBins5; ++iEt) {
	for (int iEta=0; iEta<DYTools::nEtaBins5; ++iEta) {
	  (*arr)[kind][iEt][iEta]=
	    (debug_pseudo_exps) ? ((kind+1)*100 + (iEt+1)*10 + iEta+1) :
	    gRandom->Gaus(0.0,1.0);
	}
      }
    }
    arr= & ro_MC[i];
    for (int kind=0; kind<NEffTypes; ++kind) {
      for (int iEt=0; iEt<DYTools::nEtBins5; ++iEt) {
	for (int iEta=0; iEta<DYTools::nEtaBins5; ++iEta) {
	  (*arr)[kind][iEt][iEta]=
	    (debug_pseudo_exps) ? -((kind+1)*100 + (iEt+1)*10 + iEta+1) :
	    gRandom->Gaus(0.0,1.0);
	}
      }
    }
  }

  // Create container for data for error estimates based on pseudo-experiments
  TH1F *systScale[nMassBins][nexp];
  TH1F *systScaleReco[nMassBins][nexp];
  TH1F *systScaleId[nMassBins][nexp];
  TH1F *systScaleHlt[nMassBins][nexp];
  for(int i=0; i<nMassBins; i++)
    for(int j=0; j<nexp; j++){
      TString base = "hScaleM_mass";
      base += i;
      base += "_exp";
      base += j;
      systScale[i][j] = new TH1F(base,base,150,0.0,1.5);
      systScaleReco[i][j] = new TH1F(base+TString("_reco"),
				     base+TString("_reco"),150,0.0,1.5);
      systScaleId [i][j] = new TH1F(base+TString("_id" ),
				    base+TString("_id" ),150,0.0,1.5);
      systScaleHlt[i][j] = new TH1F(base+TString("_hlt"),
				    base+TString("_hlt"),150,0.0,1.5);
    }
  
  TFile *skimFile=new TFile(selectEventsFName);
  if (!skimFile || !skimFile->IsOpen()) {
    std::cout << "failed to open file <" << selectEventsFName << ">\n";
    assert(0);
  }
  TTree *skimTree = (TTree*)skimFile->Get("Events");
  assert(skimTree);
  esfSelectEvent_t selData;
  selData.setBranchAddress(skimTree);

  std::cout << "there are " << skimTree->GetEntries() 
	    << " entries in the <" << selectEventsFName << "> file\n";
  for (UInt_t ientry=0; ientry<skimTree->GetEntries(); ++ientry) {
    if (debugMode && (ientry>10000)) break;
    if ( ientry%10000 == 0 ) std::cout << "ientry=" << ientry << "\n";

    skimTree->GetEntry(ientry);

    double scaleFactor = findEventScaleFactor(selData);
    double scaleFactorReco = sqrt(findEventScaleFactor(0,selData));
    double scaleFactorId  = sqrt(findEventScaleFactor(1,selData));
    double scaleFactorHlt = sqrt(findEventScaleFactor(2,selData));
    double weight=selData.weight;
    hScale->Fill(scaleFactor, weight);
    hScaleReco->Fill( scaleFactorReco, weight);
    hScaleId ->Fill( scaleFactorId, weight);
    hScaleHlt->Fill( scaleFactorHlt, weight);
    // Use generator-level post-FSR mass 
    int ibin = findMassBin(selData.genMass);
    hScaleRecoV[ibin]->Fill( scaleFactorReco, weight);
    hScaleIdV [ibin]->Fill( scaleFactorId, weight);
    hScaleHltV[ibin]->Fill( scaleFactorHlt, weight);
    hScaleV   [ibin]->Fill( scaleFactor, weight);
    
    hLeadingEtV [ibin]->Fill( selData.et_1, weight);
    hTrailingEtV[ibin]->Fill( selData.et_2, weight);
    hElectronEtV[ibin]->Fill( selData.et_1, weight);
    hElectronEtV[ibin]->Fill( selData.et_2, weight);
    if( selData.insideMassWindow(60,120) ) {
      hZpeakEt->Fill(selData.et_1, weight);
      hZpeakEt->Fill(selData.et_2, weight);
    }

    // Acumulate pseudo-experiments for error estimate
    for(int iexp = 0; iexp<nexp; iexp++){
      scaleFactor = findEventScaleFactorSmeared(selData, iexp);
      scaleFactorReco = sqrt(findEventScaleFactorSmeared(0,selData,iexp));
      scaleFactorId  = sqrt(findEventScaleFactorSmeared(1,selData,iexp));
      scaleFactorHlt = sqrt(findEventScaleFactorSmeared(2,selData,iexp));
      systScale   [ibin][iexp]->Fill(scaleFactor, weight);
      systScaleReco[ibin][iexp]->Fill(scaleFactorReco, weight);
      systScaleId[ibin][iexp]->Fill(scaleFactorId, weight);
      systScaleHlt[ibin][iexp]->Fill(scaleFactorHlt, weight);
    }
    
    // 	if(scaleFactor>1.3)
    // 	  printf("  leading:   %f    %f      trailing:   %f   %f     mass: %f\n",
// 		 leading->scEt, leading->scEta, trailing->scEt, trailing->scEta, dielectron->mass);
    

  } // end loop over selected events
    
  delete skimTree;
  delete skimFile;
  std::cout << "loop over selected events done\n";
  
  // Calculate errors on the scale factors
  // The "Mean" are the mean among all pseudo-experiments, very close to the primary scale factor values
  TVectorD scaleMeanV(nMassBins);
  TVectorD scaleMeanErrV(nMassBins);
  TVectorD scaleMeanRecoV(nMassBins);
  TVectorD scaleMeanRecoErrV(nMassBins);
  TVectorD scaleMeanIdV(nMassBins);
  TVectorD scaleMeanIdErrV(nMassBins);
  TVectorD scaleMeanHltV(nMassBins);
  TVectorD scaleMeanHltErrV(nMassBins);
  // Put into these vectors the content of the mean of the primary scale factor distributions
  TVectorD scaleV(nMassBins);
  TVectorD scaleRecoV(nMassBins);
  TVectorD scaleIdV(nMassBins);
  TVectorD scaleHltV(nMassBins);
  for(int ibin = 0; ibin < nMassBins; ibin++){
    scaleMeanV[ibin] = 0;
    scaleMeanErrV[ibin] = 0;
    scaleMeanRecoV[ibin] = 0;
    scaleMeanRecoErrV[ibin] = 0;
    scaleMeanIdV[ibin] = 0;
    scaleMeanIdErrV[ibin] = 0;
    scaleMeanHltV[ibin] = 0;
    scaleMeanHltErrV[ibin] = 0;
    for(int iexp = 0; iexp < nexp; iexp++){
      scaleMeanV[ibin] += systScale[ibin][iexp]->GetMean();
      scaleMeanErrV[ibin] += SQR(systScale[ibin][iexp]->GetMean());

      scaleMeanRecoV[ibin] += systScaleReco[ibin][iexp]->GetMean();
      scaleMeanRecoErrV[ibin] += SQR(systScaleReco[ibin][iexp]->GetMean());

      scaleMeanIdV[ibin] += systScaleId[ibin][iexp]->GetMean();
      scaleMeanIdErrV[ibin] += SQR(systScaleId[ibin][iexp]->GetMean());

      scaleMeanHltV[ibin] += systScaleHlt[ibin][iexp]->GetMean();
      scaleMeanHltErrV[ibin] += SQR(systScaleHlt[ibin][iexp]->GetMean());
    }
    scaleRecoV[ibin] = hScaleRecoV[ibin]->GetMean();
    scaleIdV [ibin] = hScaleIdV [ibin] ->GetMean();
    scaleHltV[ibin] = hScaleHltV[ibin]->GetMean();
    scaleV   [ibin] = hScaleV   [ibin]->GetMean();

    scaleMeanV[ibin] = scaleMeanV[ibin]/double(nexp);
    scaleMeanErrV[ibin] = sqrt( scaleMeanErrV[ibin] / double(nexp) 
				- scaleMeanV[ibin]*scaleMeanV[ibin] ); 
				
    scaleMeanRecoV[ibin] = scaleMeanRecoV[ibin]/double(nexp);
    scaleMeanRecoErrV[ibin] = sqrt( scaleMeanRecoErrV[ibin] / double(nexp)
				- scaleMeanRecoV[ibin]*scaleMeanRecoV[ibin] ); 
				
    scaleMeanIdV[ibin] = scaleMeanIdV[ibin]/double(nexp);
    scaleMeanIdErrV[ibin] = sqrt( scaleMeanIdErrV[ibin] / double(nexp) 
				- scaleMeanIdV[ibin]*scaleMeanIdV[ibin] ); 
				
    scaleMeanHltV[ibin] = scaleMeanHltV[ibin]/double(nexp);
    scaleMeanHltErrV[ibin] = sqrt( scaleMeanHltErrV[ibin] / double(nexp) 
				- scaleMeanHltV[ibin]*scaleMeanHltV[ibin] ); 
				
  }

  // Store constants in the file
  TString outputDir(TString("../root_files/constants/")+dirTag);
  gSystem->mkdir(outputDir,kTRUE);
  TString sfConstFileName(outputDir+TString("/scale_factors_") + analysisTag + 
			  TString("_") +
			  triggers.triggerConditionsName() + 
			  TString(".root"));

  TFile fa(sfConstFileName, "recreate");
  scaleV.Write("scaleFactorArray");
  scaleMeanErrV.Write("scaleFactorErrArray");
  fa.Close();

  //
  // Some plots
  //

  TCanvas *c1 = new TCanvas("c1","c1",10,10,500,500);
  c1->Divide(2,2);
 
  c1->cd(1);
  hScale->Draw();

  c1->cd(2);
  hScaleReco->Draw();

  c1->cd(3);
  hScaleId->Draw();

  c1->cd(4);
  hScaleHlt->Draw();

  std::cout << "\nScale factors as a function of mass bin\n";
  std::cout << "    mass          rho_reco          rho_id"
	    << "       rho_hlt        rho_total\n";
  std::string format1=
    std::string("   %3.0f - %3.0f     %5.3f +- %5.3f    %5.3f +- %5.3f")+
    std::string("     %5.3f +- %5.3f     %5.3f +- %5.3f\n");
  for(int i=0; i<nMassBins; i++){
    printf(format1.c_str(),
	   massBinLimits[i], massBinLimits[i+1],
	   hScaleRecoV[i]->GetMean(), scaleMeanRecoErrV[i],
	   hScaleIdV[i]->GetMean(),  scaleMeanIdErrV[i],
	   hScaleHltV[i]->GetMean(), scaleMeanHltErrV[i],
	   hScaleV[i]->GetMean()   , scaleMeanErrV[i]);
  }

  drawEfficiencies();
  drawScaleFactors();
  drawEventScaleFactors(scaleRecoV, scaleMeanRecoErrV,
			scaleIdV , scaleMeanIdErrV ,
			scaleHltV, scaleMeanHltErrV,
			scaleV   , scaleMeanErrV    );

  //
  // Make plots of Et spectra
  //
  // Normalize first
  for(int i=0; i<nMassBins; i++){
    printf("Total events in mass bin %3d     %10.0f\n", 
	   i, hLeadingEtV[i]->GetSumOfWeights());
    hLeadingEtV  [i]->Sumw2();
    hTrailingEtV [i]->Sumw2();
    hElectronEtV [i]->Sumw2();
    hLeadingEtV [i]->Scale(1.0/hLeadingEtV[i]->GetSumOfWeights());
    hTrailingEtV[i]->Scale(1.0/hTrailingEtV[i]->GetSumOfWeights());
    hElectronEtV[i]->Scale(1.0/hElectronEtV[i]->GetSumOfWeights());
  }
  printf("Total events around Z peak    %10.0f\n", 
	 hZpeakEt->GetSumOfWeights()/2.0);
  hZpeakEt->Sumw2();
  hZpeakEt->Scale(1.0/hZpeakEt->GetSumOfWeights());

  TCanvas *c3 = MakeCanvas("ET canvas 1", "ET canvas 1");
  CPlot etplot1("etplot_bin0", "","E_{T} [GeV]", "N_{ele}, normalized");
  int cbin = 0;
  TString label = "mass bin 0";
  etplot1.SetLogx();
  etplot1.AddHist1D(hElectronEtV[cbin] , label , "PE", kRed);
  etplot1.AddHist1D(hZpeakEt       , "60-120 mass range", "hist,f", kBlack);
  etplot1.SetYRange(0.0,0.6);
  hElectronEtV[cbin]->GetXaxis()->SetMoreLogLabels();
  hElectronEtV[cbin]->GetXaxis()->SetNoExponent();
  hZpeakEt->SetFillStyle(3001);
  hElectronEtV[cbin]->SetMarkerColor(kRed);
  etplot1.Draw(c3,savePlots,"png");

  TCanvas *c4 = MakeCanvas("ET canvas 2", "ET canvas 2");
  CPlot etplot2("etplot_bin3", "","E_{T} [GeV]", "N_{ele}, normalized");
  cbin = 3;
  label = "mass bin 3";
  etplot2.SetLogx();
  etplot2.AddHist1D(hElectronEtV[cbin] , label , "PE", kRed);
  etplot2.AddHist1D(hZpeakEt       , "60-120 mass range", "hist,f", kBlack);
  etplot2.SetYRange(0.0,0.6);
  hElectronEtV[cbin]->GetXaxis()->SetMoreLogLabels();
  hElectronEtV[cbin]->GetXaxis()->SetNoExponent();
  hZpeakEt->SetFillStyle(3001);
  hZpeakEt->SetMarkerColor(kBlack);
  hElectronEtV[cbin]->SetMarkerColor(kRed);
  etplot2.Draw(c4,savePlots,"png");

  TCanvas *c5 = MakeCanvas("ET canvas 3", "ET canvas 3");
  CPlot etplot3("etplot_bin5", "","E_{T} [GeV]", "N_{ele}, normalized");
  cbin = 5;
  label = "mass bin 5";
  etplot3.SetLogx();
  etplot3.AddHist1D(hElectronEtV[cbin] , label , "PE", kRed);
  etplot3.AddHist1D(hZpeakEt       , "60-120 mass range", "hist,f", kBlack);
  etplot3.SetYRange(0.0,0.6);
  hElectronEtV[cbin]->GetXaxis()->SetMoreLogLabels();
  hElectronEtV[cbin]->GetXaxis()->SetNoExponent();
  hZpeakEt->SetFillStyle(3001);
  hZpeakEt->SetMarkerColor(kBlack);
  hElectronEtV[cbin]->SetMarkerColor(kRed);
  etplot3.Draw(c5,savePlots,"png");

  gBenchmark->Show("calcEventEff");
  return;
}

// ------------------------------------------------------------

Bool_t matchedToGeneratorLevel(const TGenInfo *gen, 
			       const TDielectron *dielectron){

  Bool_t result = kTRUE;
  // In the generator branch of this ntuple, first particle is always
  // negative, and second always positive. In the Dielectron block
  // of the ntuple, the first particle is always the one with larger Pt.
  double dR1=999, dR2=999;
  TLorentzVector v1reco, v2reco, v1gen, v2gen;
  v1reco.SetPtEtaPhiM(dielectron->pt_1, dielectron->eta_1, 
		      dielectron->phi_1, 0.000511);
  v2reco.SetPtEtaPhiM(dielectron->pt_2, dielectron->eta_2, 
		      dielectron->phi_2, 0.000511);
  v1gen .SetPtEtaPhiM(gen->pt_1, gen->eta_1, gen->phi_1, 0.000511);
  v2gen .SetPtEtaPhiM(gen->pt_2, gen->eta_2, gen->phi_2, 0.000511);
  if( dielectron->q_1 < 0 ){
    dR1 = v1reco.DeltaR(v1gen);
    dR2 = v2reco.DeltaR(v2gen);
  }else{
    dR1 = v1reco.DeltaR(v2gen);
    dR2 = v2reco.DeltaR(v1gen);
  }
  // Require that both are within loose dR of 0.4, otherwise bail out
  if( fabs(dR1) > 0.4 || fabs(dR2) > 0.4 ) result = kFALSE; 
  
  return result;
}

// ------------------------------------------------------------

int createSelectionFile(const MCInputFileMgr_t &mcMgr, 
     const TString &outSkimFName, TriggerSelection &triggers, int debugMode) {
  int eventsInNtuple = 0;
  double weightedEventsInNtuple = 0;
  int eventsAfterTrigger = 0;
  int totalCand = 0;
  int totalCandInMassWindow = 0;
  int totalCandInEtaAcceptance = 0;
  int totalCandEtAbove10GeV = 0;
  int totalCandMatchedToGen = 0;
  int totalCandFullSelection = 0;

#ifdef esfSelectEventsIsObject
  esfSelectEvent_t::Class()->IgnoreTObjectStreamer();
#endif
  esfSelectEvent_t selData;
  TFile *skimFile=new TFile(outSkimFName,"recreate");
  cout << "createSelectionFileName=<" << outSkimFName << ">\n";
  if (!skimFile || !skimFile->IsOpen()) {
    cout << "createSelectionFile: failed to create a file <" 
	 << outSkimFName << ">\n";
    return 0;
  }
  TTree *skimTree= new TTree("Events","Events");
  assert(skimTree);
  selData.createBranches(skimTree);

  double lumi0=1.;

  // Loop over files
  for(UInt_t ifile=0; ifile<mcMgr.size(); ifile++){

    //
    // Access samples and fill histograms
    //  
    TFile *infile = 0;
    TTree *eventTree = 0;
        
    // Data structures to store info from TTrees
    mithep::TEventInfo *info = new mithep::TEventInfo();
    mithep::TGenInfo   *gen  = new mithep::TGenInfo();
    TClonesArray *dielectronArr = new TClonesArray("mithep::TDielectron");
    TClonesArray *pvArr   = new TClonesArray("mithep::TVertex");
    
    // Read input file
    TString fname=mcMgr.fileName(ifile);
    cout << "Processing " << fname << "..." << endl;
    infile = new TFile(fname); 
    assert(infile);
    
    // Get the TTrees
    eventTree = (TTree*)infile->Get("Events"); assert(eventTree);

    // Find weight for events for this file
    // The first file in the list comes with weight 1,
    // all subsequent ones are normalized to xsection and luminosity
    double xsec=mcMgr.xsec(ifile);
    AdjustXSectionForSkim(infile,xsec,eventTree->GetEntries(),1);
    double lumi = eventTree->GetEntries()/xsec;
    if (ifile==0) lumi0=lumi;
    double weight = lumi0/lumi;
    cout << "       -> sample weight is " << weight << endl;

    // Set branch address to structures that will store the info  
    eventTree->SetBranchAddress("Info",&info);
    TBranch *infoBr       = eventTree->GetBranch("Info");
    eventTree->SetBranchAddress("Dielectron",&dielectronArr);
    TBranch *dielectronBr = eventTree->GetBranch("Dielectron");
    eventTree->SetBranchAddress("Gen",&gen);
    TBranch *genBr = eventTree->GetBranch("Gen");
    eventTree->SetBranchAddress("PV", &pvArr); 
    TBranch *pvBr         = eventTree->GetBranch("PV");

    // loop over events    
    eventsInNtuple         += eventTree->GetEntries();
    weightedEventsInNtuple += weight * eventTree->GetEntries();
    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
//       for(UInt_t ientry=0; ientry<100000; ientry++) { // This is for faster turn-around in testing
      if (debugMode && (ientry>100000)) break;
      genBr->GetEntry(ientry);
      infoBr->GetEntry(ientry);
      
      /* old trigger defs
      // For EPS2011 for both data and MC (starting from Summer11 production)
      // we use an OR of the twi triggers below. Both are unpresecaled.
      ULong_t eventTriggerBit = kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL 
	| kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL;
      ULong_t leadingTriggerObjectBit = kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele1Obj
	| kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele1Obj;
      ULong_t trailingTriggerObjectBit = kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele2Obj
	| kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele2Obj;
      */
      ULong_t eventTriggerBit = triggers.getEventTriggerBit(info->runNum);
      ULong_t leadingTriggerObjectBit  = 
	triggers.getLeadingTriggerObjectBit(info->runNum);
      ULong_t trailingTriggerObjectBit = 
	triggers.getTrailingTriggerObjectBit(info->runNum);
      
      if(!(info->triggerBits & eventTriggerBit)) continue;  // no trigger accept? Skip to next event...
      eventsAfterTrigger++;
      
      // check the number of goodPVs
      pvBr->GetEntry(ientry);
      int nGoodPV=countGoodVertices(pvArr);

      // loop through dielectrons
      dielectronArr->Clear();
      dielectronBr->GetEntry(ientry);    

      for(Int_t i=0; i<dielectronArr->GetEntriesFast(); i++) {
	
	totalCand++;
	const mithep::TDielectron *dielectron = 
	  (mithep::TDielectron*)((*dielectronArr)[i]);

	// Consider only events in the mass range of interest
	// Use generator level post-FSR mass.
 	if( gen->mass < massBinLimits[0] || 
	    gen->mass > massBinLimits[nMassBins]) continue;
	totalCandInMassWindow++;

	// Exclude ECAL gap region (should already be done for ntuple, but just to make sure...)
	if((fabs(dielectron->scEta_1)>kECAL_GAP_LOW) &&
	   (fabs(dielectron->scEta_1)<kECAL_GAP_HIGH)) continue;
	if((fabs(dielectron->scEta_2)>kECAL_GAP_LOW) &&
	   (fabs(dielectron->scEta_2)<kECAL_GAP_HIGH)) continue;
	// ECAL acceptance cut on supercluster Et
	if((fabs(dielectron->scEta_1) > 2.5)       || 
	   (fabs(dielectron->scEta_2) > 2.5)) continue;  // outside eta range? Skip to next event...
	totalCandInEtaAcceptance++;
	// None of the electrons should be below 10 GeV
	if((dielectron->scEt_1 < 10)            || 
	   (dielectron->scEt_2 < 10))	      continue;  // below supercluster ET cut? Skip to next event...
	totalCandEtAbove10GeV++;

	// For MC-only, do generator level matching
	if( ! matchedToGeneratorLevel(gen, dielectron) ) continue;
	totalCandMatchedToGen++;

	// ECAL-driven reconstruction is not required.

	// No cut on opposite charges to avoid systematics related to charge mis-ID	
	//  	if( (dielectron->q_1 == dielectron->q_2 )) continue;

	TElectron *ele1 = extractElectron(dielectron, 1);
	TElectron *ele2 = extractElectron(dielectron, 2);

	// ID cuts
	if( !( passSmurf(ele1) && passSmurf(ele2) ) ) continue;

	// ET and trigger cut on the leading electron
	const TElectron *leading = ele1;
	const TElectron *trailing = ele2;
	if( ele1->scEt < ele2->scEt ){
	  leading = ele2;
	  trailing = ele1;
	}

	// The individual electron trigger match should be done
	// exactly the same way as it is done in the signal selection.
	// (At the moment of this writing it is a bit different, needs to
	// be fixed).
	if( !( leading->scEt  > 20 && 
	       (leading ->hltMatchBits & leadingTriggerObjectBit) ) ) continue;
	if( !( trailing->scEt > 10 && 
	      (trailing->hltMatchBits & trailingTriggerObjectBit) ) ) continue;
	totalCandFullSelection++;

	selData.assign(gen->mass, dielectron->mass, 
		       leading->scEt, leading->scEta,
		       trailing->scEt, trailing->scEta,
		       weight,
		       nGoodPV
		       );
	skimTree->Fill();
 
	// 	  printf("  leading:   %f    %f      trailing:   %f   %f     mass: %f\n",
	// 		 leading->scEt, leading->scEta, trailing->scEt, trailing->scEta, dielectron->mass);

	
      } // end loop over dielectrons
    } // end loop over events
    
    delete infile;
    infile = 0;
    eventTree = 0;
    delete gen;
    delete info;
    delete dielectronArr;
    delete pvArr;
  } // end loop over files
  
  skimFile->cd();
  skimTree->Write();
  skimFile->Close();
  delete skimFile;

  printf("Total events in ntuple(s)                          %15d\n",
	 eventsInNtuple);
  printf("    number of weighted events in ntuple            %17.1lf\n",
	 weightedEventsInNtuple);
  printf("    events after event level trigger cut           %15d\n",
	 eventsAfterTrigger);
  printf("\nTotal candidates (no cuts)                       %15d\n",
	 totalCand);
  printf("        candidates in 15-600 mass window           %15d\n",
	 totalCandInMassWindow);
  printf("        candidates with eta 0-1.4442, 1.566-2.5    %15d\n",
	 totalCandInEtaAcceptance);
  printf("        candidates, both electrons above 10 GeV    %15d\n",
	 totalCandEtAbove10GeV);
  printf("        candidates matched to GEN level (if MC)    %15d\n",
	 totalCandMatchedToGen);
  printf("        candidates, full selection                 %15d\n",
	 totalCandFullSelection);
  return 1;
}


// ------------------------------------------------------------

double findEventScaleFactor(const esfSelectEvent_t &data) {

  double esf1=1.0;
  int etBin1 = findEtBin(data.et_1, etBinning);
  int etaBin1 = findEtaBin(data.eta_1, etaBinning);
  if ((etBin1!=-1) && (etaBin1!=-1)) {
    esf1=
      findScaleFactor(0, etBin1, etaBin1) *
      findScaleFactor(1, etBin1, etaBin1) *
      findScaleFactor(2, etBin1, etaBin1);
  }

  double esf2=1.0;
  int etBin2 = findEtBin(data.et_2, etBinning);
  int etaBin2 = findEtaBin(data.eta_2, etaBinning);
  if ((etBin2!=-1) && (etaBin2!=-1)) {
    esf2=
      findScaleFactor(0, etBin2, etaBin2) *
      findScaleFactor(1, etBin2, etaBin2) *
      findScaleFactor(2, etBin2, etaBin2);
  }

  return esf1*esf2;
}

// ------------------------------------------------------------

double findEventScaleFactor(int kind, const esfSelectEvent_t &data) {

  double esf1=1.0;
  int etBin1 = findEtBin(data.et_1, etBinning);
  int etaBin1 = findEtaBin(data.eta_1, etaBinning);
  if ((etBin1!=-1) && (etaBin1!=-1)) {
    esf1= findScaleFactor(kind, etBin1, etaBin1);
  }

  double esf2=1.0;
  int etBin2 = findEtBin(data.et_2, etBinning);
  int etaBin2 = findEtaBin(data.eta_2, etaBinning);
  if ((etBin2!=-1) && (etaBin2!=-1)) {
    esf2= findScaleFactor(kind, etBin2, etaBin2);
  }

  return esf1*esf2;
}

// --------------------------------------

double findScaleFactor(int kind, double scEt, double scEta) {

  int etBin = findEtBin(scEt, etBinning);
  int etaBin = findEtaBin(scEta, etaBinning);
  return findScaleFactor(kind,etBin,etaBin);

}

// --------------------------------------

double findScaleFactor(int kind, int etBin, int etaBin) {

  double result = 0;
  if( (etBin == -1) || (etaBin == -1)) {
    // Found bin outside of calibrated range, return 1.0
    result = 1.0;
    return result;
  }

  result = (*dataEff[kind])[etBin][etaBin] / (*mcEff[kind])[etBin][etaBin];
  return result;
}

// ---------------------- all scale factors smeared ------------------------------

// ------------------------------------------------------------

double findEventScaleFactorSmeared(const esfSelectEvent_t &data, int iexp) {

  double esf1=1.0;
  int etBin1 = findEtBin(data.et_1, etBinning);
  int etaBin1 = findEtaBin(data.eta_1, etaBinning);
  if ((etBin1!=-1) && (etaBin1!=-1)) {
    esf1=
      findScaleFactorSmeared(0, etBin1, etaBin1, iexp) *
      findScaleFactorSmeared(1, etBin1, etaBin1, iexp) *
      findScaleFactorSmeared(2, etBin1, etaBin1, iexp);
  }

  double esf2=1.0;
  int etBin2 = findEtBin(data.et_2, etBinning);
  int etaBin2 = findEtaBin(data.eta_2, etaBinning);
  if ((etBin2!=-1) && (etaBin2!=-1)) {
    esf2=
      findScaleFactorSmeared(0, etBin2, etaBin2, iexp) *
      findScaleFactorSmeared(1, etBin2, etaBin2, iexp) *
      findScaleFactorSmeared(2, etBin2, etaBin2, iexp);
  }

  return esf1*esf2;
}

// ------------------------------------------------------------

double findEventScaleFactorSmeared(int kind, const esfSelectEvent_t &data,
				   int iexp) {

  double esf1=1.0;
  int etBin1 = findEtBin(data.et_1, etBinning);
  int etaBin1 = findEtaBin(data.eta_1, etaBinning);
  if ((etBin1!=-1) && (etaBin1!=-1)) {
    esf1= findScaleFactorSmeared(kind, etBin1, etaBin1, iexp);
  }

  double esf2=1.0;
  int etBin2 = findEtBin(data.et_2, etBinning);
  int etaBin2 = findEtaBin(data.eta_2, etaBinning);
  if ((etBin2!=-1) && (etaBin2!=-1)) {
    esf2= findScaleFactorSmeared(kind, etBin2, etaBin2, iexp);
  }

  return esf1*esf2;
}

// --------------------------------------

double findScaleFactorSmeared(int kind, double scEt, double scEta, int iexp) {

  int etBin = findEtBin(scEt, etBinning);
  int etaBin = findEtaBin(scEta, etaBinning);
  return findScaleFactorSmeared(kind,etBin,etaBin,ro_Data[iexp],ro_MC[iexp]);

}

// --------------------------------------

double findScaleFactorSmeared(int kind, int etBin, int etaBin,
	   const EffArray_t &dataRndWeight, const EffArray_t &mcRndWeight) {

  double result = 0;
  if( (etBin == -1) || (etaBin == -1)) {
    // Found bin outside of calibrated range, return 1.0
    result = 1.0;
    return result;
  }

  double effData=
    (*dataEff[kind])[etBin][etaBin] + 
    dataRndWeight[kind][etBin][etaBin] * (*dataEffAvgErr[kind])[etBin][etaBin];
  if (effData>100.) effData=100.;
  double effMC=
    (*mcEff[kind])[etBin][etaBin] + 
    mcRndWeight[kind][etBin][etaBin] * (*mcEffAvgErr[kind])[etBin][etaBin];
  if (effMC>100.) effMC=100.;
  return (effData/effMC);
}

// --------------------------------------
// --------------------------------------


 void drawEfficiencyGraphs(TGraphErrors *grData, TGraphErrors *grMc,
			   TString yAxisTitle, TString text, TString plotName){
   
   // Generate "random" canvas name
//    TTimeStamp time;
//    TString c = "c";
//    c += time.GetNanoSec();
   TString c = plotName;
//    printf("Canvas name %s\n", c.Data());

   TCanvas *c2 = MakeCanvas(c,c);
   CPlot plot1(c,"","E_{T} [GeV]", yAxisTitle);
   plot1.SetLogx(); 
   plot1.AddGraph(grData,"data","PE2", kBlue);
   plot1.AddGraph(grMc  ,"MC"  ,"PE", kBlack);
   plot1.Draw(c2);
   plot1.SetYRange(0.2,1.1);
   grData->GetXaxis()->SetTitle("E_{T} [GeV]");
   grData->GetXaxis()->SetMoreLogLabels();
   grData->GetXaxis()->SetNoExponent();
   grData->SetFillStyle(3001);
   grData->SetFillColor(kBlue);
   grMc->SetMarkerStyle(24);
   plot1.TransLegend(0.0, -0.5);
   plot1.AddTextBox(text, 0.6,0.4,0.87,0.55, 0);

   TLine *line = new TLine(10,1.0,500,1.0);
   line->SetLineStyle(kDashed);
   line->Draw("same");

   plot1.Draw(c2, savePlots, "png");

   return;
 }

// -------------------------------------------------------------------------
/*
template<class Graph_t>
void drawEfficiencyGraphsAsymmErrors(Graph_t *grData, Graph_t *grMc,
                     TString yAxisTitle, TString text, TString plotName){
  
   // Generate "random" canvas name
//    TTimeStamp time;
//    TString c = "c";
//    c += time.GetNanoSec();
   TString c = plotName;
//    printf("Canvas name %s\n", c.Data());

   TCanvas *c2 = MakeCanvas(c,c);
   CPlot plot1(c,"","E_{T} [GeV]", yAxisTitle);
   plot1.SetLogx(); 
   plot1.AddGraph(grData,"data","PE2", kBlue);
   plot1.AddGraph(grMc  ,"MC"  ,"PE", kBlack);
   plot1.Draw(c2);
   plot1.SetYRange(0.2,1.1);
   grData->GetXaxis()->SetTitle("E_{T} [GeV]");
   grData->GetXaxis()->SetMoreLogLabels();
   grData->GetXaxis()->SetNoExponent();
   grData->SetFillStyle(3001);
   grData->SetFillColor(kBlue);
   grMc->SetMarkerStyle(24);
   plot1.TransLegend(0.0, -0.5);
   plot1.AddTextBox(text, 0.6,0.4,0.8,0.5, 0);
   plot1.Draw(c2, savePlots, "png");

   TLine *line = new TLine(10,1.0,500,1.0);
   line->SetLineStyle(kDashed);
   line->Draw("same");

   return;
 }

// -------------------------------------------------------------------------

template<class Graph_t>
void drawEfficiencyGraphsAsymmErrorsPU(Graph_t *grData, Graph_t *grMc,
				     const TString &yAxisTitle, 
				     const TString &text, 
				     const TString &plotName){
  
   // Generate "random" canvas name
//    TTimeStamp time;
//    TString c = "c";
//    c += time.GetNanoSec();
   TString c = plotName;
//    printf("Canvas name %s\n", c.Data());

   TString xAxisTitle = "nGoodPV";
   TCanvas *c2 = MakeCanvas(c,c);
   CPlot plot1(c,"",xAxisTitle, yAxisTitle);
   //plot1.SetLogx(); 
   plot1.AddGraph(grData,"data","PE2", kBlue);
   plot1.AddGraph(grMc  ,"MC"  ,"PE", kBlack);
   plot1.Draw(c2);
   plot1.SetYRange(0.2,1.1);
   grData->GetXaxis()->SetTitle(xAxisTitle);
   //grData->GetXaxis()->SetMoreLogLabels();
   //grData->GetXaxis()->SetNoExponent();
   grData->SetFillStyle(3001);
   grData->SetFillColor(kBlue);
   grMc->SetMarkerStyle(24);
   plot1.TransLegend(0.0, -0.5);
   plot1.AddTextBox(text, 0.6,0.4,0.8,0.5, 0);
   plot1.Draw(c2, savePlots, "png");

   TLine *line = new TLine(0,1.0,DYTools::nPVLimits[DYTools::nPVBinCount],1.0);
   line->SetLineStyle(kDashed);
   line->Draw("same");

   return;
 }
*/

// -------------------------------------------------------------------------

void drawScaleFactorGraphs(TGraphErrors *gr, TString yAxisTitle, TString text,
			   TString plotName){
   
   // Generate "random" canvas name
//    TTimeStamp time;
//    TString c = "c";
//    c += time.GetNanoSec();
  TString c = plotName;
//   printf("Canvas name %s\n", c.Data());

   TCanvas *c2 = MakeCanvas(c,c);
   CPlot plot1(c,"","E_{T} [GeV]", yAxisTitle);
   plot1.SetLogx(); 
   plot1.AddGraph(gr,"present E_{T}/#eta dep.","PE", kBlack);
   plot1.Draw(c2);
   plot1.SetYRange(0.5,1.5);
   gr->GetXaxis()->SetTitle("E_{T} [GeV]");
   gr->GetXaxis()->SetMoreLogLabels();
   gr->GetXaxis()->SetNoExponent();
   plot1.AddTextBox(text, 0.6,0.35,0.87,0.5, 0);
   
   TLine *line = new TLine(10,1.0,500,1.0);
   line->SetLineStyle(kDashed);
   line->Draw("same");

   plot1.Draw(c2,savePlots, "png");

   return;
 }

// -------------------------------------------------------------------------

void drawEventScaleFactorGraphs(TGraphErrors *gr, TString yAxisTitle, 
				TString plotName){
  
  // Generate "random" canvas name
//   TTimeStamp time;
//   TString c = "c";
//   c += time.GetNanoSec();
  TString c = plotName;
//   printf("Canvas name %s\n", c.Data());
  
  TCanvas *c2 = MakeCanvas(c,c);
  CPlot plot1(c,"","m(e^{+}e^{-}) [GeV]", yAxisTitle);
  plot1.SetLogx(); 
  plot1.AddGraph(gr,"present E_{T}/#eta dep.","PE", kBlack);
  plot1.Draw(c2);
  plot1.SetYRange(0.0,1.5);
  gr->GetXaxis()->SetTitle("m(e^{+}e^{-}) [GeV]");
  gr->GetXaxis()->SetMoreLogLabels();
  gr->GetXaxis()->SetNoExponent();
//   cout << "From main progam CPlot::sOutDir " << CPlot::sOutDir << endl;
  CPlot::sOutDir = "plots";
  
  TLine *line = new TLine(15,1.0,500,1.0);
  line->SetLineStyle(kDashed);
  line->Draw("same");

  plot1.Draw(c2, savePlots, "png");

  return;
 }

// -------------------------------------------------------------------------

void drawEfficiencies(){
  // Make graphs
  double x[etBinCount];
  double dx[etBinCount];
  for(int i=0; i<etBinCount; i++){
    x[i]  = 0.5*(etBinLimits[i] + etBinLimits[i+1]);
    dx[i] = 0.5*(etBinLimits[i+1] - etBinLimits[i]);
  }


  double effData[etBinCount],effDataErr[etBinCount];
  double effMC[etBinCount],effMCErr[etBinCount];
  char bufEta[30];
  char plotLabel[30];

  for (int kind=0; kind<3; ++kind) {
    for (int iEta=0; iEta<etaBinCount; ++iEta) {
      sprintf(bufEta,"_abs_eta_%5.3lf__%5.3lf",
	      etaBinLimits[iEta],etaBinLimits[iEta+1]);
      TString etaStr=bufEta;
      etaStr.ReplaceAll(".","_");
      sprintf(plotLabel,"%5.3lf < |#eta| < %5.3lf",
	      etaBinLimits[iEta],etaBinLimits[iEta+1]);

      for (int iEt=0; iEt<etBinCount; ++iEt) {
	effData[iEt]= (*dataEff[kind])[iEt][iEta];
	effDataErr[iEt]= (*dataEffAvgErr[kind])[iEt][iEta];
      }
      for (int iEt=0; iEt<etBinCount; ++iEt) {
	effMC[iEt]= (*mcEff[kind])[iEt][iEta];
	effMCErr[iEt]= (*mcEffAvgErr[kind])[iEt][iEta];
      }
      TGraphErrors *grDataEff 
	= new TGraphErrors(etBinCount, x,  effData, dx, effDataErr);
  
      TGraphErrors *grMcEff 
	= new TGraphErrors(etBinCount, x,  effMC, dx, effMCErr);
  
      // Draw all graphs
      TString effName=EfficiencyKindName(DYTools::TEfficiencyKind_t(kind));
      TString ylabel=TString("efficiency_{") + effName + TString("}");
      TString plotName = TString("plot_eff_") + analysisTag + TString("_") +
	effName + etaStr;
      
      drawEfficiencyGraphs(grDataEff, grMcEff,
 		       ylabel, plotLabel, plotName);

      //delete grMcEff;
      //delete grDataEff;
    }
  }

  return;
}

// -------------------------------------------------------------------------

void drawScaleFactors(){

  double x[etBinCount];
  double dx[etBinCount];

  for(int i=0; i<etBinCount; i++){
    x[i]  = (etBinLimits[i] + etBinLimits[i+1])/2.0;
    dx[i] = (etBinLimits[i+1] - etBinLimits[i])/2.0;
  }


  double scale[etBinCount];
  double scaleErr[etBinCount];
  char bufEta[30];
  char plotLabel[30];

  for (int kind=0; kind<3; ++kind) {
    for (int iEta=0; iEta<etaBinCount; ++iEta) {
      sprintf(bufEta,"_abs_eta_%5.3lf__%5.3lf",
	      etaBinLimits[iEta],etaBinLimits[iEta+1]);
      TString etaStr=bufEta;
      etaStr.ReplaceAll(".","_");
      sprintf(plotLabel,"%5.3lf < |#eta| < %5.3lf",
	      etaBinLimits[iEta],etaBinLimits[iEta+1]);

      for (int iEt=0; iEt<etBinCount; ++iEt) {
	scale[iEt]= (*dataEff[kind])[iEt][iEta] / (*mcEff[kind])[iEt][iEta];
	scaleErr[iEt]= errOnRatio( (*dataEff[kind])[iEt][iEta], 
				   (*dataEffAvgErr[kind])[iEt][iEta],
				   (*mcEff[kind])[iEt][iEta], 
				   (*mcEffAvgErr[kind])[iEt][iEta] );
      }

      TGraphErrors *grScaleFactor
	= new TGraphErrors(etBinCount, x, scale, dx, scaleErr);
      TString effName=EfficiencyKindName(DYTools::TEfficiencyKind_t(kind));
      TString ylabel=TString("scale factor ") + effName;
      TString plotName = TString("plot_scale_") + analysisTag + TString("_") +
	effName + etaStr;

      drawScaleFactorGraphs(grScaleFactor, ylabel, plotLabel, plotName);
      //delete grScaleFactor;
    }
  }
}

// -------------------------------------------------------------------------

double errOnRatio(double a, double da, double b, double db){

  double result = 0;
  if(a == 0 || b == 0)
    return result;
 
  result = (a/b)*sqrt( (da/a)*(da/a) + (db/b)*(db/b) );
  return result;
}

// -------------------------------------------------------------------------

void drawEventScaleFactors(TVectorD scaleRecoV, TVectorD scaleRecoErrV,
			   TVectorD scaleIdV , TVectorD scaleIdErrV ,
			   TVectorD scaleHltV, TVectorD scaleHltErrV,
			   TVectorD scaleV   , TVectorD scaleErrV    )
{

  // repackage into arrays
  double x[nMassBins];
  double dx[nMassBins];
  double scaleRecoA   [nMassBins];
  double scaleIdA    [nMassBins];
  double scaleHltA   [nMassBins];
  double scaleA      [nMassBins];
  double scaleRecoErrA[nMassBins];
  double scaleIdErrA [nMassBins];
  double scaleHltErrA[nMassBins];
  double scaleErrA   [nMassBins];
  for(int i=0; i<nMassBins; i++){
    x[i] = (massBinLimits[i] + massBinLimits[i+1])/2.0;
    dx[i]= (massBinLimits[i+1] - massBinLimits[i])/2.0;
    scaleRecoA       [i] = scaleRecoV   [i];
    scaleIdA        [i] = scaleIdV    [i];
    scaleHltA       [i] = scaleHltV   [i];
    scaleA          [i] = scaleV      [i];
    scaleRecoErrA    [i] = scaleRecoErrV[i];
    scaleIdErrA     [i] = scaleIdErrV [i];
    scaleHltErrA    [i] = scaleHltErrV[i];
    scaleErrA       [i] = scaleErrV   [i];
  }

  TGraphErrors *grScale = 
    new TGraphErrors(nMassBins, x, scaleA, dx, scaleErrA);
  TGraphErrors *grScaleReco = 
    new TGraphErrors(nMassBins, x, scaleRecoA, dx, scaleRecoErrA);
  TGraphErrors *grScaleId  = 
    new TGraphErrors(nMassBins, x, scaleIdA , dx, scaleIdErrA );
  TGraphErrors *grScaleHlt = 
    new TGraphErrors(nMassBins, x, scaleHltA, dx, scaleHltErrA);

  TString plotName;
  TString plotNameBase = 
    TString("plot_event_scale_") + analysisTag + TString("_");
  plotName = plotNameBase + TString("reco");
  drawEventScaleFactorGraphs(grScaleReco, "RECO scale factor" , plotName);
  plotName = "plot_event_scale_id";
  plotName = plotNameBase + TString("id");
  drawEventScaleFactorGraphs(grScaleId , "ID scale factor"   , plotName);
  plotName = plotNameBase + TString("hlt");
  drawEventScaleFactorGraphs(grScaleHlt, "HLT scale factor"  , plotName);
  plotName = plotNameBase + TString("full");
  drawEventScaleFactorGraphs(grScale   , "event scale factor", plotName);

}

// -------------------------------------------------------------------------

// This method reads all ROOT files that have efficiencies from
// tag and probe in TMatrixD form and converts the matrices into 
// more simple arrays.
int fillEfficiencyConstants(  const TnPInputFileMgr_t &mcMgr, 
      const TnPInputFileMgr_t &dataMgr, const TriggerSelection &triggers ) {

  TString fnStart="efficiency_TnP_"; //+ analysisTag;
  TString fnEnd=".root";

  if (dataEff.size()) dataEff.clear();
  if (dataEffErrLo.size()) dataEffErrLo.clear();
  if (dataEffErrHi.size()) dataEffErrHi.clear();
  if (mcEff.size()) mcEff.clear();
  if (mcEffErrLo.size()) mcEffErrLo.clear();
  if (mcEffErrHi.size()) mcEffErrHi.clear();
  dataEff.reserve(3); dataEffErrLo.reserve(3); dataEffErrHi.reserve(3);
  mcEff.reserve(3); mcEffErrLo.reserve(3); mcEffErrHi.reserve(3);

  int res=1;
  for (int kind=0; res && (kind<3); ++kind) {
    TString dataFName=fnStart + 
      getLabel(DATA,TEfficiencyKind_t(kind),dataMgr.effCalcMethod(kind),
	       etBinning, etaBinning, triggers)
      + fnEnd;
    res=fillOneEfficiency(dataMgr, dataFName, kind, 
			  dataEff, dataEffErrLo, dataEffErrHi, dataEffAvgErr);
  }
  for (int kind=0; res && (kind<3); ++kind) {
    TString mcFName=fnStart + 
      getLabel(MC,TEfficiencyKind_t(kind),mcMgr.effCalcMethod(kind),
	       etBinning, etaBinning, triggers)
      + fnEnd;
    res=fillOneEfficiency(mcMgr, mcFName, kind, 
			  mcEff, mcEffErrLo, mcEffErrHi, mcEffAvgErr);
  }
  if (res!=1) std::cout << "Error in fillEfficiencyConstants\n"; 
  else std::cout << "fillEfficiencyConstants ok\n";
  return res;
}

// -------------------------------------------------------------------------

int fillOneEfficiency(const TnPInputFileMgr_t &mgr, const TString filename, 
   UInt_t kindIdx, vector<TMatrixD*> &effV, vector<TMatrixD*> &errLoV, 
   vector<TMatrixD*> &errHiV, vector<TMatrixD*> &avgErrV) {

  TFile f(TString("../root_files/tag_and_probe/")+mgr.dirTag()+TString("/")+
	  filename);
  if(!f.IsOpen()) assert(0);
  std::cout << "reading <" << filename << ">\n";

  TMatrixD *effMatrix        = (TMatrixD*)f.Get("effArray2D");
  TMatrixD *effMatrixErrLow  = (TMatrixD*)f.Get("effArrayErrLow2D");
  TMatrixD *effMatrixErrHigh = (TMatrixD*)f.Get("effArrayErrHigh2D");
  f.Close();

  // Make sure that the objects are present
  if( !(effMatrix && effMatrixErrLow && effMatrixErrHigh) ) assert(0);

  // Make sure that there are only two eta bins and appropriate number of ET bins
  if( effMatrix->GetNcols() != etaBinCount ) {
    std::cout << "The number of eta bins stored in constants file ("
	      << effMatrix->GetNcols() << ") is not " <<
	      << etaBinCount << "\n";
    return 0;
  }
  if( effMatrix->GetNrows() != etBinCount ) {
    std::cout << "The number of ET bins stored in constants file ("
	      << effMatrix->GetNrows() << ") is different form expected ("
	      << etBinCount << ")\n";
    return 0;
  }

  if ((effV.size()!=kindIdx) ||
      (errLoV.size()!=kindIdx) ||
      (errHiV.size()!=kindIdx) ||
      (avgErrV.size()!=kindIdx)) {
    cout << "Error: effV.size=" << effV.size() 
	 << ", errLoV.size=" << errLoV.size() 
	 << ", errHiV.size=" << errHiV.size() 
	 << ", avgErrV.size=" << avgErrV.size() 
	 << ", kindIdx=" <<  kindIdx << "\n";
    return 0;
  }

  effV.push_back(effMatrix);
  errLoV.push_back(effMatrixErrLow);
  errHiV.push_back(effMatrixErrHigh);
  TMatrixD* avgErr=(TMatrixD*)effMatrixErrLow->Clone("avgErr");
  for (int col=0; col<effMatrixErrLow->GetNcols(); ++col) {
    for (int row=0; row<effMatrixErrLow->GetNrows(); ++row) {
      (*avgErr)[row][col]=
	0.5*((*effMatrixErrLow)[row][col] + (*effMatrixErrHigh)[row][col]);
    }
  }
  avgErrV.push_back(avgErr);

  return 1;
}


// -------------------------------------------------------------------------

void PrintEffInfoLines(const char *msg, int effKind, int effMethod, 
		       int binCount, const double *eff, const double *effErr) {
  std::cout << "PrintEffInfoLines(" << ((msg) ? msg : "<null>") 
	    << ", effKind=" << effKind << ", effMethod=" << effMethod 
	    << ", binCount=" << binCount << "\n";
  for (int i=0; i<binCount; i++) {
    printf("  i=%d  eff=%6.4lf effErr=%8.6e\n",i,eff[i],effErr[i]);
  }
  std::cout << std::endl;
  return;
}


// -------------------------------------------------------------------------

