#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TChain.h>
#include <TCanvas.h>                // class for drawing
#include <TH1F.h>                   // 1D histograms
#include <TGraphErrors.h>           // graphs
#include <TProfile.h>               // profile histograms
#include <TBenchmark.h>             // class to track macro running statistics
#include <TLorentzVector.h>         // class for Lorentz vector computations
#include <TGraphAsymmErrors.h>
#include <TClonesArray.h>
#include <TMatrixD.h>
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <vector>                   // STL vector class
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings

#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooBreitWigner.h"
#include "RooExponential.h"
#include "RooCBShape.h"
#include "RooGenericPdf.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooFFTConvPdf.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooGaussian.h"

#include "../Include/CPlot.hh"          // helper class for plots
#include "../Include/MitStyleRemix.hh"  // style settings for drawing

// define classes and constants to read in ntuple
#include "../Include/EWKAnaDefs.hh"
#include "../Include/TGenInfo.hh"
#include "../Include/TEventInfo.hh"
#include "../Include/TDielectron.hh"
#include "../Include/TElectron.hh"

#include "../Include/DYTools.hh"

#include "../Include/EleIDCuts.hh"

#endif

using namespace mithep;


//=== COMMON CONSTANTS ===========================================================================================

const Double_t kGAP_LOW  = 1.4442;
const Double_t kGAP_HIGH = 1.566;

//=== FUNCTION DECLARATIONS ======================================================================================

Bool_t matchedToGeneratorLevel(const TGenInfo *gen, const TDielectron *dielectron);

void measurePassAndFail(double &signal, double &signalErr, double &efficiency, double &efficiencyErr,
			TTree *passTree, TTree *failTree,
			TCanvas *passCanvas, TCanvas *failCanvas);
void measureEfficiency(TTree *passTree, TTree *failTree, 
		       int method, int etBinning, int etaBinning, 
		       TCanvas *canvas, ofstream &effOutput, ofstream &fitLog,
		       bool useTemplates, TFile *templatesFile, TFile *resultsRootFile);
void measureEfficiencyCountAndCount(TTree *passTree, TTree *failTree, 
				    int etBinning, int etaBinning, 
				    TCanvas *canvas, ofstream &effOutput, bool saveToRootFile, TFile *resultsRootFile);
void measureEfficiencyWithFit(TTree *passTree, TTree *failTree, 
			      int method, int etBinning, int etaBinning, 
			      TCanvas *canvas, ofstream &effOutput, ofstream &fitLog,
			      bool useTemplates, TFile *templatesFile, TFile *resultsRootFile);
void fitMass(TTree *passTree, TTree *failTree, 
	     TString cut, int mode,
	     double &efficiency, double &efficiencyErrHigh, double &efficiencyErrLow,
	     TPad *passPad, TPad *failPad, ofstream &fitLog);

void fitMassWithTemplates(TTree *passTree, TTree *failTree, 
			  TString cut, int mode,
			  double &efficiency, double &efficiencyErrHigh, double &efficiencyErrLow, 
			  TPad *passPad, TPad *failPad, ofstream &fitLog,
			  TH1F *templatePass, TH1F *templateFail);

bool isTag(const TElectron *electron, UInt_t trigger);
bool passID(const TElectron *electron);

void printCorrelations(ostream& os, RooFitResult *res);
TString getLabel(int sample, int effType, int method, 
		 int etBinning, int etaBinning);

int getTemplateBin(int etBin, int etaBin, int etaBinning);
TH1F * getPassTemplate(int etBin, int etaBin, int etaBinning, TFile *file);
TH1F * getFailTemplate(int etBin, int etaBin, int etaBinning, TFile *file);

//=== MAIN MACRO =================================================================================================

void eff_IdHlt(const TString configFile) 
{

  gBenchmark->Start("eff_IdHlt");
  
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================
  
  Double_t massLow  = 60;
  Double_t massHigh = 120;

  // Read in the configuratoin file
  TString sampleTypeString = "";
  TString effTypeString    = "";
  TString calcMethodString = "";
  TString etBinningString  = "";
  TString etaBinningString = "";
  TString dirTag;
  vector<TString> ntupleFileNames;
  ifstream ifs;
  ifs.open(configFile.Data());
  assert(ifs.is_open());
  string line;
  Int_t state=0;
  while(getline(ifs,line)) {
    if(line[0]=='#') continue;
    if(line[0]=='%') break;
    if(state==0){
      // Read 1st line of content: data or MC?
      sampleTypeString = TString(line);
      state++;
    }else if(state==1) {
      // Read 2d content line: efficiency type string
      effTypeString = TString(line);
      state++;
    }else if(state==2) {
      // Read 3d content line: fitting mode
      calcMethodString = TString(line);
      state++;
    }else if(state==3) {
      // Read 4th content line: SC ET binning
      etBinningString = TString(line);
      state++;
    }else if(state==4) {
      // Read 5th content line: SC eta binning
      etaBinningString = TString(line);
      state++;
    }else if(state==5) {
      // Read 5th content line: SC eta binning
      dirTag = TString(line);
      state++;
    }else if(state==6) {
      ntupleFileNames.push_back(TString(line));
    }
  }
  
  int calcMethod = 0;
  if(calcMethodString == "COUNTnCOUNT")
    calcMethod = COUNTnCOUNT;
  else if(calcMethodString == "COUNTnFIT")
    calcMethod = COUNTnFIT;
  else if(calcMethodString == "FITnFIT")
    calcMethod = FITnFIT;
  else
    assert(0);
  printf("Efficiency calculation method: %s\n", calcMethodString.Data());

  int effType = 0;
  if(effTypeString == "ID")
    effType = ID;
  else if(effTypeString == "HLT")
    effType = HLT;
  else
    assert(0);
  printf("Efficiency type to measure: %s\n", effTypeString.Data());

  int etBinning = 0;
  if(etBinningString == "ETBINS1")
    etBinning = ETBINS1;
  else if(etBinningString == "ETBINS5")
    etBinning = ETBINS5;
  else
    assert(0);
  printf("SC ET binning: %s\n", etBinningString.Data());

  int etaBinning = 0;
  if(etaBinningString == "ETABINS1")
    etaBinning = ETABINS1;
  else if(etaBinningString == "ETABINS2")
    etaBinning = ETABINS2;
  else
    assert(0);
  printf("SC eta binning: %s\n", etaBinningString.Data());

  int sample;
  if(sampleTypeString == "DATA")
    sample = DATA;
  else if(sampleTypeString == "MC")
    sample = MC;
  else
    assert(0);
  printf("Sample: %s\n", sampleTypeString.Data());

  // The label is a string that contains the fields that are passed to
  // the function below, to be used to name files with the output later.
  TString label = getLabel(sample, effType, calcMethod, etBinning, etaBinning);

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================
  
  //  
  // Set up histograms
  //
  TH1F* hMass           = new TH1F("hMass"     ,"",30,massLow,massHigh);
  TH1F* hMassTotal      = new TH1F("hMassTotal","",30,massLow,massHigh);
  TH1F* hMassPass       = new TH1F("hMassPass" ,"",30,massLow,massHigh);
  TH1F* hMassFail       = new TH1F("hMassFail" ,"",30,massLow,massHigh);

  // Save MC templates if sample is MC
  TString tagAndProbeDir(TString("../root_files/tag_and_probe/")+dirTag);
  gSystem->mkdir(tagAndProbeDir,kTRUE);

  TFile *templatesFile = 0;
  vector<TH1F*> hPassTemplateV;
  vector<TH1F*> hFailTemplateV;
  if( sample != DATA) {
    // For simulation, we will be saving templates
    TString labelMC = "";
    labelMC += label.Data();
    // For the file name of templates
    labelMC.Replace(labelMC.Index("count-count_"),12,"");
    TString templatesLabel = tagAndProbeDir + TString("/mass_templates_")+labelMC+TString(".root");
    templatesFile = new TFile(templatesLabel,"recreate");
    for(int i=0; i<getNEtBins(etBinning); i++){
      for(int j=0; j<getNEtaBins(etaBinning); j++){
	TString hname = "hMassTemplate_Et";
	hname += i;
	hname += "_eta";
	hname += j;
	hPassTemplateV.push_back(new TH1F(hname+TString("_pass"),"",60,massLow,massHigh));
	hFailTemplateV.push_back(new TH1F(hname+TString("_fail"),"",60,massLow,massHigh));
      }
    }
  } else {
    // For data, we will be using templates,
    // however, if the request is COUNTnCOUNT, do nothing
    if( calcMethod != COUNTnCOUNT ){
      TString labelMC = "";
      labelMC += label.Data();
      // Find the corresponding label of MC templates. For this,
      // strip calculation method and replace data with MC in the label
      labelMC.Replace(labelMC.Index("data"),4,"mc");
      if( calcMethod == COUNTnFIT)
	labelMC.Replace(labelMC.Index("count-fit_"),10,"");
      if( calcMethod == FITnFIT)
	labelMC.Replace(labelMC.Index("fit-fit_"),8,"");
      TString templatesLabel = tagAndProbeDir+TString("/mass_templates_")+labelMC+TString(".root");
      templatesFile = new TFile(templatesLabel);
      if( ! templatesFile->IsOpen() )
	assert(0);
    }
  }

  // This file can be utilized in the future, but for now
  // opening it just removes complaints about memory resident
  // trees. No events are actually written.
  TFile *selectedEventsFile = new TFile("selectedEventsFile.root","recreate");

  TTree *passTree = new TTree("passTree","passTree");
  Double_t storeMass, storeEt, storeEta;
  passTree->Branch("mass",&storeMass,"mass/D");
  passTree->Branch("et",&storeEt  ,"et/D");
  passTree->Branch("eta",&storeEta ,"eta/D");

  TTree *failTree = new TTree("failTree","failTree");
  failTree->Branch("mass",&storeMass,"mass/D");
  failTree->Branch("et",&storeEt  ,"et/D");
  failTree->Branch("eta",&storeEta ,"eta/D");

  int nDivisions = getNEtBins(etBinning)*getNEtaBins(etaBinning);
  double ymax = 800;
  if(nDivisions <4 )
    ymax = nDivisions * 200;
  TCanvas *c1 = MakeCanvas("c1","c1", 600, ymax);
  c1->Divide(2,nDivisions);

  int eventsInNtuple = 0;
  int eventsAfterTrigger = 0;
  int totalCand = 0;
  int totalCandInMassWindow = 0;
  int totalCandInEtaAcceptance = 0;
  int totalCandEtAbove10GeV = 0;
  int totalCandMatchedToGen = 0;

  // Loop over files
  for(UInt_t ifile=0; ifile<ntupleFileNames.size(); ifile++){

    //
    // Access samples and fill histograms
    //  
    TFile *infile = 0;
    TTree *eventTree = 0;
        
    // Data structures to store info from TTrees
    mithep::TEventInfo *info = new mithep::TEventInfo();
    mithep::TGenInfo   *gen  = new mithep::TGenInfo();
    TClonesArray *dielectronArr = new TClonesArray("mithep::TDielectron");
    
    // Read input file
    cout << "Processing " << ntupleFileNames[ifile] << "..." << endl;
    infile = new TFile(ntupleFileNames[ifile]); 
    assert(infile);
    
    // Get the TTrees
    eventTree = (TTree*)infile->Get("Events"); assert(eventTree);

    // Set branch address to structures that will store the info  
    eventTree->SetBranchAddress("Info",&info);                TBranch *infoBr       = eventTree->GetBranch("Info");
    eventTree->SetBranchAddress("Dielectron",&dielectronArr); TBranch *dielectronBr = eventTree->GetBranch("Dielectron");
    TBranch *genBr = 0;
    if(sample != DATA){
      eventTree->SetBranchAddress("Gen",&gen);
      genBr = eventTree->GetBranch("Gen");
    }

    // loop over events    
    eventsInNtuple += eventTree->GetEntries();
    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
//       for(UInt_t ientry=0; ientry<100000; ientry++) { // This is for faster turn-around in testing
       
      if(sample != DATA)
	genBr->GetEntry(ientry);
      
      // Check that the whole event has fired the appropriate trigger
      infoBr->GetEntry(ientry);
      
      // For EPS2011 for both data and MC (starting from Summer11 production)
      // we use a special trigger for tag and probe that has second leg
      // unbiased with cuts at HLT
      UInt_t eventTriggerBit = kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30
	| kHLT_Ele32_CaloIdL_CaloIsoVL_SC17;
      // The tag trigger bit matches the "electron" of the trigger we
      // use for this tag and probe study: electron+sc
      UInt_t tagTriggerObjectBit = kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_EleObj
	| kHLT_Ele32_CaloIdL_CaloIsoVL_SC17_EleObj;
      // The probe trigger, however, is any of possibilities used in
      // the trigger that is used in the main analysis
      UInt_t probeTriggerObjectBit = kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele1Obj
	| kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele1Obj
	| kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele2Obj
	| kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele2Obj;

      if(!(info->triggerBits & eventTriggerBit)) continue;  // no trigger accept? Skip to next event... 
      eventsAfterTrigger++;
      
      // loop through dielectrons
      dielectronArr->Clear();
      dielectronBr->GetEntry(ientry);    
      for(Int_t i=0; i<dielectronArr->GetEntriesFast(); i++) {
	
	totalCand++;
	const mithep::TDielectron *dielectron = (mithep::TDielectron*)((*dielectronArr)[i]);
	
	// Tag and probe is done around the Z peak
	if((dielectron->mass < massLow) || (dielectron->mass > massHigh)) continue;
	totalCandInMassWindow++;
	//
	// Exclude ECAL gap region (should already be done for ntuple, but just to make sure...)
	if((fabs(dielectron->scEta_1)>kGAP_LOW) && (fabs(dielectron->scEta_1)<kGAP_HIGH)) continue;
	if((fabs(dielectron->scEta_2)>kGAP_LOW) && (fabs(dielectron->scEta_2)<kGAP_HIGH)) continue;
	// ECAL acceptance cut on supercluster Et
	if((fabs(dielectron->scEta_1) > 2.5)       || (fabs(dielectron->scEta_2) > 2.5)) continue;  // outside eta range? Skip to next event...
	totalCandInEtaAcceptance++;
	// None of the electrons should be below 10 GeV
	if((dielectron->scEt_1 < 10)               || (dielectron->scEt_2 < 10))	      continue;  // below supercluster ET cut? Skip to next event...
	totalCandEtAbove10GeV++;
	
	// Next, we will do a loose kinematic matching to generator level
	// info. 
	// For the data, this is not needed and not done. We take all
	// candidates, and take care of background by fitting.
	// For MC, however, we do not fit, but count pass/fail events.
	// So we need to make sure there is no background. However, even
	// in the signal Z->ee MC sample there jets and therefore fake
	// electrons. So we drop all candidates that do not have both leptons
	// matched.
	// 
	if( sample != DATA )
	  if( ! matchedToGeneratorLevel(gen, dielectron) ) continue;
	totalCandMatchedToGen++;

	// ECAL driven: this condition is NOT applied	

	// Preliminary selection is complete. Now work on tags and probes.
	
	TElectron *ele1 = extractElectron(dielectron, 1);
	TElectron *ele2 = extractElectron(dielectron, 2);
	bool isTag1 = isTag(ele1, tagTriggerObjectBit);
	bool isTag2 = isTag(ele2, tagTriggerObjectBit);
	
	// Any electron that made it here is eligible to be a probe
	// for ID cuts.
	bool isIDProbe1     = true;
	bool isIDProbe2     = true;
	bool isIDProbePass1 = passID(ele1);
	bool isIDProbePass2 = passID(ele2);
	
	// Probes for HLT cuts:

	bool isHLTProbe1     = passID(ele1);
	bool isHLTProbe2     = passID(ele2);
	bool isHLTProbePass1 = ( isHLTProbe1 && (ele1 ->hltMatchBits & probeTriggerObjectBit) );
	bool isHLTProbePass2 = ( isHLTProbe2 && (ele2 ->hltMatchBits & probeTriggerObjectBit) );

	// 
	//  Apply tag and probe, and accumulate counters or histograms
	//       
	
	bool isProbe1     = false;
	bool isProbe2     = false;
	bool isProbePass1 = false;
	bool isProbePass2 = false;
	if( effType == ID ){
	  isProbe1     = isIDProbe1;
	  isProbe2     = isIDProbe2;
	  isProbePass1 = isIDProbePass1;
	  isProbePass2 = isIDProbePass2;
	}else if( effType == HLT ){
	  isProbe1     = isHLTProbe1;
	  isProbe2     = isHLTProbe2;
	  isProbePass1 = isHLTProbePass1;
	  isProbePass2 = isHLTProbePass2;
	}else {
	  printf("ERROR: unknown efficiency type requested\n");
	}

	storeMass = dielectron->mass;
	// First electron is the tag, second is the probe
	if( isTag1 && isProbe2){
	  // total probes
	  hMassTotal->Fill(dielectron->mass);
	  storeEt   = dielectron->scEt_2;
	  storeEta  = dielectron->scEta_2;
	  int templateBin = getTemplateBin( findEtBin(storeEt,etBinning),
					    findEtaBin(storeEta,etaBinning),
					    etaBinning);
	  if( isProbePass2 ){
	    // passed
	    hMassPass->Fill(dielectron->mass);
	    passTree->Fill();
	    if(sample != DATA && templateBin != -1)
	      hPassTemplateV[templateBin]->Fill(dielectron->mass);
	  }else{
	    // fail
	    hMassFail->Fill(dielectron->mass);
	    failTree->Fill();
	    if(sample != DATA && templateBin != -1)
	      hFailTemplateV[templateBin]->Fill(dielectron->mass);
	  }
	}
	// Second electron is the tag, first is the probe
	if( isTag2 && isProbe1 ){
	  // total probes
	  hMassTotal->Fill(dielectron->mass);
	  storeEt   = dielectron->scEt_1;
	  storeEta  = dielectron->scEta_1;
	  int templateBin = getTemplateBin( findEtBin(storeEt,etBinning),
					    findEtaBin(storeEta,etaBinning),
					    etaBinning);
	  if( isProbePass1 ){
	    // passed
	    hMassPass->Fill(dielectron->mass);
	    passTree->Fill();
	    if(sample != DATA && templateBin != -1)
	      hPassTemplateV[templateBin]->Fill(dielectron->mass);
	  }else{
	    // fail
	    hMassFail->Fill(dielectron->mass);
	    failTree->Fill();
	    if(sample != DATA && templateBin != -1)
	      hFailTemplateV[templateBin]->Fill(dielectron->mass);
	  }
	}
	
	// In case the full selection is applied:
	//       if( !(isTag1 && ele2_passID) && !(isTag2 && ele1_passID) ) continue;
	if( !(isTag1 && isIDProbePass2) && !(isTag2 && isIDProbePass1) ) continue;
	//       if( !(isTag1) && !(isTag2) ) continue;
	// Fill histogram
	hMass->Fill(dielectron->mass);
	
      } // end loop over dielectron candidates
    } // end loop over events
  
    delete infile;
    infile=0;
    eventTree=0;
    
    delete gen;
    delete info;
    delete dielectronArr;
  } // end loop over files

  //
  // Efficiency analysis
  //
  
//   printf("Number of regular candidates:      %15.0f\n", hMass->GetSumOfWeights());
  printf("Total events in ntuple                                       %15d\n",eventsInNtuple);
  printf("    events after event level trigger cut                     %15d\n",eventsAfterTrigger);
  printf("\nTotal candidates (no cuts)                                   %15d\n",totalCand);
  printf("        candidates in 60-120 mass window                     %15d\n",totalCandInMassWindow);
  printf("        candidates witheta 0-1.4442, 1.566-2.5               %15d\n",totalCandInEtaAcceptance);
  printf("        candidates, both electrons above 10 GeV              %15d\n",totalCandEtAbove10GeV);
  printf("        candidates matched to GEN level (if MC)              %15d\n",totalCandMatchedToGen);
  printf("\nNumber of probes, total                                      %15.0f\n", hMassTotal->GetSumOfWeights());
  printf("Number of probes, passed                                     %15.0f\n", hMassPass->GetSumOfWeights());
  printf("Number of probes, failed                                     %15.0f\n", hMassFail->GetSumOfWeights());


  // Human-readbale text file to store measured efficiencies
  TString reslog = tagAndProbeDir+TString("/efficiency_TnP_")+label+TString(".txt");
  ofstream effOutput;
  effOutput.open(reslog);
  // Print into the results file the header.
  effOutput << "Efficiency calculation method: " << calcMethodString.Data() << endl;
  effOutput << "Efficiency type to measure: " << effTypeString.Data() << endl;
  effOutput << "SC ET binning: " << etBinningString.Data() << endl;
  effOutput << "SC eta binning: " << etaBinningString.Data() << endl;
  effOutput << "Sample: " << sampleTypeString.Data() << endl;
  effOutput << "Files processed: " << endl;
  for(UInt_t i=0; i<ntupleFileNames.size(); i++)
    effOutput << "   " << ntupleFileNames[i].Data() << endl;
  
  // ROOT file to store measured efficiencies in ROOT format
  TString resroot = tagAndProbeDir+TString("/efficiency_TnP_")+label+TString(".root");
  TFile *resultsRootFile = new TFile(resroot,"recreate");

  // Fit log 
  TString fitlogname = tagAndProbeDir+TString("/efficiency_TnP_")+label+TString("_fitlog.dat");
  ofstream fitLog;
  fitLog.open(fitlogname);

  //
  //  Find efficiency
  //
  bool useTemplates = false;
  if(sample == DATA && effType == ID && (calcMethod == COUNTnFIT || FITnFIT) )
    useTemplates = true;

  measureEfficiency(passTree, failTree,
		    calcMethod, etBinning, etaBinning, c1, effOutput, fitLog,
		    useTemplates, templatesFile, resultsRootFile);

  effOutput.close();
  fitLog.close();
  TString command = "cat ";
  command += reslog;
  system(command.Data());

  TString fitpicname = tagAndProbeDir+TString("/efficiency_TnP_")+label+TString("_fit.png");
  c1->SaveAs(fitpicname);

  // Save MC templates
  if(sample != DATA){
    templatesFile->cd();
    for(int i=0; i<getNEtBins(etBinning); i++){
      for(int j=0; j<getNEtaBins(etaBinning); j++){
	int templateBin = getTemplateBin( i, j, etaBinning);
	hPassTemplateV[templateBin]->Write();
	hFailTemplateV[templateBin]->Write();
      }
    }
    templatesFile->Close();
  }

  gBenchmark->Show("eff_IdHlt");
  
  
}


//=== FUNCTION DEFINITIONS ======================================================================================


Bool_t matchedToGeneratorLevel(const TGenInfo *gen, const TDielectron *dielectron){

  Bool_t result = kTRUE;
  // In the generator branch of this ntuple, first particle is always
  // negative, and second always positive. In the Dielectron block
  // of the ntuple, the first particle is always the one with larger Pt.
  double dR1=999, dR2=999;
  TLorentzVector v1reco, v2reco, v1gen, v2gen;
  v1reco.SetPtEtaPhiM(dielectron->pt_1, dielectron->eta_1, dielectron->phi_1, 0.000511);
  v2reco.SetPtEtaPhiM(dielectron->pt_2, dielectron->eta_2, dielectron->phi_2, 0.000511);
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

void measurePassAndFail(double &signal, double &signalErr, double &efficiency, double &efficiencyErr,
			TTree *passTree, TTree *failTree,
			TCanvas *passCanvas, TCanvas *failCanvas){

  // Define data sets
  RooRealVar mass("mass","mass",60, 120);
  mass.setBins(120);
  RooRealVar pt ("et" ,"et" ,10.0, 1000);
  RooRealVar eta("eta","eta",-10, 10);
  RooDataSet *dataPass = new RooDataSet("dataPass","dataPass",passTree,RooArgSet(mass,pt,eta));

  RooDataSet *dataFail = new RooDataSet("dataFail","dataFail",RooArgList(mass,pt,eta),Import(*failTree));
  RooCategory probeType("probeType","probeType");
  probeType.defineType("pass");
  probeType.defineType("fail");
  RooAbsData *data;
  // If needed do binned fit
  bool unbinnedFit = false;
  if(unbinnedFit){
    data = new RooDataSet("data","data",mass,Index(probeType),
			  Import("pass",*dataPass), Import("fail",*dataFail));
  }else{
    RooDataHist *dataPassBinned = dataPass->binnedClone("dataPassBinned","dataPassBinned");
    RooDataHist *dataFailBinned = dataFail->binnedClone("dataFailBinned","dataFailBinned");
    data = new RooDataHist("data","data",mass,Index(probeType),
			   Import("pass",*dataPassBinned), Import("fail",*dataFailBinned));    
  }
  // Define the PDFs
  //
  // Common pieces
  //
  // True signal model
  RooRealVar zMass ("zMass" ,"zMass" ,91.188);
  RooRealVar zWidth("zWidth","zWidth",2.495);
  RooBreitWigner bwPdf("bwPdf","bwPdf",mass,zMass,zWidth);
  // Total signal events and efficiency
  RooRealVar nsignal("nsignal","nsignal",1000,0.0,1.0e7);
  RooRealVar eff    ("eff"    ,"eff"    ,0.7,0.0,1.0);

  //  PDF for the PASS sample
  // 
  // Background
  RooRealVar lambdaBgPass("lambdaBgPass","lambdaBgPass",-0.1, -0.5, 0.0);
  RooExponential bgPassPdf("bgPassPdf","bgPassPdf",mass,lambdaBgPass);
  // Signal
  //     - resolution function
  RooRealVar cbMeanPass("cbMeanPass","cbMeanPass"   ,0.0, -10.0,10.0);
  RooRealVar cbWidthPass("cbWidthPass","cbWidthPass",1.0,  0.1, 5.0);
  RooRealVar cbAlphaPass("cbAlphaPass","cbAlphaPass",5.0,  0.0,20.0);
  RooRealVar cbNPass("cbNPass","cbNPass"            ,1.0,  0.0,10.0);
  RooCBShape cbPassPdf("cbPassPdf","cbPassPdf",mass,cbMeanPass,cbWidthPass,cbAlphaPass,cbNPass);
  //     - realistic model
  mass.setBins(10000,"fft");
  RooFFTConvPdf signalPassPdf("signalPassPdf","signalPassPdf",mass, bwPdf, cbPassPdf);
  // Combine signal and background
  RooFormulaVar nsigPass("nsigPass","nsigPass","@0*@1",RooArgList(nsignal,eff));
  RooRealVar nbgPass ("nbgPass" ,"nbgPass" ,1,0.0,1.0e5);
  RooAddPdf passPdf("passPdf","passPdf",RooArgList(signalPassPdf,bgPassPdf), RooArgList(nsigPass,nbgPass));
  
  //
  //  PDF for the FAIL sample
  // 
  // Background
  RooRealVar lambdaBgFail("lambdaBgFail","lambdaBgFail",-0.1, -0.5, 0.0);
  RooExponential bgFailPdf("bgFailPdf","bgFailPdf",mass,lambdaBgFail);
  // Signal
  //     - resolution function
  RooRealVar cbMeanFail("cbMeanFail","cbMeanFail"   ,0.0, -10.0,10.0);
  RooRealVar cbWidthFail("cbWidthFail","cbWidthFail",1.0,   0.1, 5.0);
  RooRealVar cbAlphaFail("cbAlphaFail","cbAlphaFail",5.0,   0.0,20.0);
  RooRealVar cbNFail("cbNFail","cbNFail"            ,1.0,   0.0,10.0);
  RooCBShape cbFailPdf("cbFailPdf","cbFailPdf",mass,cbMeanFail,cbWidthFail,cbAlphaFail,cbNFail);
  //     - realistic model
  RooFFTConvPdf signalFailPdf("signalFailPdf","signalFailPdf",mass, bwPdf, cbFailPdf);
  // Combine signal and background
  RooFormulaVar nsigFail("nsigFail","nsigFail","@0*(1.0-@1)",RooArgList(nsignal,eff));
  RooRealVar nbgFail ("nbgFail" ,"nbgFail" ,1,0.0,1.0e5);
  RooAddPdf failPdf("failPdf","failPdf",RooArgList(signalFailPdf,bgFailPdf), RooArgList(nsigFail,nbgFail));
  
  // Combine pass and fail
  RooSimultaneous fullPdf("fullPdf","fullPdf",probeType);
  fullPdf.addPdf(passPdf,"pass");
  fullPdf.addPdf(failPdf,"fail");

  // Do the fit

  // Start with a reasonable point and do rough approximation first
  double total = dataPass->numEntries()+dataFail->numEntries();
  nsignal.setVal( 0.99*total);
  eff.setVal(0.90);
  nbgPass.setVal(0.01*total);
  nbgFail.setVal(0.01*total);
  cbAlphaPass.setVal(1.0);
  cbAlphaFail.setVal(0.5);
  cbNPass    .setVal(5.0);
  cbNFail    .setVal(5.0);
  cbAlphaPass.setConstant(kTRUE);
  cbAlphaFail.setConstant(kTRUE);
  cbNPass    .setConstant(kTRUE);
  cbNFail    .setConstant(kTRUE);
  RooFitResult *result = fullPdf.fitTo(*data,Extended(kTRUE),Save());

  // Release shape parameters and refine the fit
  cbAlphaPass.setConstant(kFALSE);
  cbAlphaFail.setConstant(kFALSE);
  cbNPass    .setConstant(kFALSE);
  cbNFail    .setConstant(kFALSE);
  result = fullPdf.fitTo(*data,Extended(kTRUE),Save());

  cout << "Fit status 1st iteration " << result->status() << endl;
//   if(!result->status()){
//     result = fullPdf.fitTo(*data,Extended(kTRUE),Save());
//     cout << "Fit status 2d iteration " << result->status() << endl;
//   }

  // Plot
  passCanvas->cd();
  passCanvas->SetWindowPosition(0,0);
  passCanvas->Draw();
  RooPlot *framePass = mass.frame();
  dataPass->plotOn(framePass);
  passPdf.plotOn(framePass);
  passPdf.plotOn(framePass,Components("bgPassPdf"),LineStyle(kDashed));
  framePass->Draw();

  failCanvas->cd();
  failCanvas->SetWindowPosition(0+ failCanvas->GetWindowWidth(),0);
  failCanvas->Draw();
  RooPlot *frameFail = mass.frame();
  dataFail->plotOn(frameFail);
  failPdf.plotOn(frameFail);
  failPdf.plotOn(frameFail,Components("bgFailPdf"),LineStyle(kDashed));
  frameFail->Draw();

  signal        = nsignal.getVal();
  signalErr     = nsignal.getError();
  efficiency    = eff.getVal();
  efficiencyErr = eff.getError();

  return;
}

bool isTag(const TElectron *electron, UInt_t trigger){

  bool elePassID  = passID(electron);
  bool elePassHLT =  (electron ->hltMatchBits & trigger);

  bool result = ( elePassID && elePassHLT && (electron->scEt > 20) );

  return result;
}

bool passID(const TElectron *electron){

  bool result = passSmurf(electron);
  return result;
}

void measureEfficiency(TTree *passTree, TTree *failTree, 
		       int method, int etBinning, int etaBinning, 
		       TCanvas *canvas, ofstream &effOutput, ofstream &fitLog,
		       bool useTemplates, TFile *templatesFile, TFile *resultsRootFile){

  // For COUNTnCOUNT method we should write to root file results
  // from measureEfficiencyCountAndCount routine, otherwise
  // from measureEfficiencyWithFit routine.
  bool saveCountingToRootFile = true;
  if( method == COUNTnFIT || method == FITnFIT )
    saveCountingToRootFile = false;

  // Always report counting method results
  measureEfficiencyCountAndCount(passTree, failTree, etBinning, etaBinning, 
				 canvas, effOutput, saveCountingToRootFile, resultsRootFile);

  if( method == COUNTnFIT || method == FITnFIT )
    measureEfficiencyWithFit(passTree, failTree, 
			     method, etBinning, etaBinning, 
			     canvas, effOutput, fitLog,
			     useTemplates, templatesFile, resultsRootFile);
  
  return;
}

void measureEfficiencyCountAndCount(TTree *passTree, TTree *failTree, 
				    int etBinning, int etaBinning, 
				    TCanvas *canvas, ofstream &effOutput, 
				    bool saveResultsToRootFile, TFile *resultsRootFile){

  int nEt                = getNEtBins(etBinning);
  const double *limitsEt = getEtBinLimits(etBinning);

  int nEta                = getNEtaBins(etaBinning);
  const double *limitsEta = getEtaBinLimits(etaBinning);
  printf("eta bins %d\n", nEta);

  TMatrixD effArray2D(nEt, nEta);
  TMatrixD effArrayErrLow2D(nEt, nEta);
  TMatrixD effArrayErrHigh2D(nEt, nEta);
 
  effOutput << endl;
  effOutput << "Efficiency, counting method:\n";  
  effOutput << "     SC ET         SC eta           efficiency             pass         fail\n";
  for(int j=0; j<nEta; j++){
    for(int i=0; i<nEt; i++){
      double effCount, effErrLowCount, effErrHighCount;
      TString etCut = TString::Format(" ( et >=%6.1f && et <%6.1f ) ",
				      limitsEt[i], limitsEt[i+1]);
      TString etaCut = TString::Format(" ( abs(eta) >= %5.3f && abs(eta) < %5.3f ) ",
				       limitsEta[j], limitsEta[j+1]);
      TString cut = etCut + TString(" && ") + etaCut;
      //cout << cut.Data() << endl;
      double probesPass = passTree->GetEntries(cut);
      double probesFail = failTree->GetEntries(cut);
      DYTools::calcEfficiency( probesPass, probesPass+probesFail, DYTools::EFF_CLOPPER_PEARSON,
			       effCount, effErrLowCount, effErrHighCount);
      char strOut[200];
      sprintf(strOut, "   %3.0f - %3.0f   %5.3f - %5.3f   %5.1f +%5.1f -%5.1f    %10.0f  %10.0f\n",
	     limitsEt[i], limitsEt[i+1],
	     limitsEta[j], limitsEta[j+1],
	     effCount*100, effErrLowCount*100, effErrHighCount*100,
	     probesPass, probesFail);
      effOutput << strOut;
      canvas->cd(1 + 2*(i + j*nEt) + 0);
      passTree->Draw("mass",cut);
      canvas->cd(1 + 2*(i + j*nEt) + 1);
      failTree->Draw("mass",cut);
      canvas->Update();
      effArray2D(i,j) = effCount;
      effArrayErrLow2D(i,j) = effErrLowCount;
      effArrayErrHigh2D(i,j) = effErrHighCount;
    }
  }
  effOutput << endl;

  if(saveResultsToRootFile){
    if(resultsRootFile && resultsRootFile->IsOpen()){
      resultsRootFile->cd();
      effArray2D.Write("effArray2D");
      effArrayErrLow2D.Write("effArrayErrLow2D");
      effArrayErrHigh2D.Write("effArrayErrHigh2D");    
      resultsRootFile->Close();
    }else assert(0);
  }

  return;
}

void measureEfficiencyWithFit(TTree *passTree, TTree *failTree, 
			      int method, int etBinning, int etaBinning, 
			      TCanvas *canvas, ofstream &effOutput, ofstream &fitLog,
			      bool useTemplates, TFile *templatesFile, TFile *resultsRootFile){
  
  int nEt                = getNEtBins(etBinning);
  const double *limitsEt = getEtBinLimits(etBinning);

  int nEta                = getNEtaBins(etaBinning);
  const double *limitsEta = getEtaBinLimits(etaBinning);
  printf("eta bins %d\n", nEta);
  
  TMatrixD effArray2D(nEt, nEta);
  TMatrixD effArrayErrLow2D(nEt, nEta);
  TMatrixD effArrayErrHigh2D(nEt, nEta);
 
  effOutput << endl;
  effOutput << "Efficiency, Count+Fit method:\n";  
  effOutput << "     SC ET         SC eta                efficiency             pass         fail\n";
  for(int j=0; j<nEta; j++){
    for(int i=0; i<nEt; i++){
      TString etCut = TString::Format(" ( et >=%6.1f && et <%6.1f ) ",
				      limitsEt[i], limitsEt[i+1]);
      TString etaCut = TString::Format(" ( abs(eta) >= %5.3f && abs(eta) < %5.3f ) ",
				       limitsEta[j], limitsEta[j+1]);
      TString cut = etCut + TString(" && ") + etaCut;
      //cout << cut.Data() << endl;
      double probesPass = passTree->GetEntries(cut);
      double probesFail = failTree->GetEntries(cut);
      TPad *passPad = (TPad*)canvas->GetPad(1 + 2*(i + j*nEt) + 0);
      TPad *failPad = (TPad*)canvas->GetPad(1 + 2*(i + j*nEt) + 1);
      double efficiency, efficiencyErrHi, efficiencyErrLo;
      printf("\n ==\n");
      char strOut[200];
      sprintf(strOut," ==   Start fitting Et: %3.0f - %3.0f  and eta:  %5.3f - %5.3f \n",
	     limitsEt[i], limitsEt[i+1],
	     limitsEta[j], limitsEta[j+1]);
      printf("%s",strOut);
      printf(" ==\n\n");
      fitLog << endl << strOut << endl;
      // In case templates are used, find the right templates
      TH1F *templatePass = getPassTemplate(i,j,etaBinning, templatesFile);
      TH1F *templateFail = getFailTemplate(i,j,etaBinning, templatesFile);
      if(!useTemplates){
	fitMass(passTree, failTree, cut, method, efficiency, efficiencyErrHi, efficiencyErrLo, 
		passPad, failPad, fitLog);
      }else{
	printf("\nMASS TEMPLATES ARE USED IN THE FIT\n\n");
	fitMassWithTemplates(passTree, failTree, cut, method, 
			     efficiency, efficiencyErrHi, efficiencyErrLo,
		passPad, failPad, fitLog, templatePass, templateFail);
      }
      sprintf(strOut, "   %3.0f - %3.0f   %5.3f - %5.3f   %5.1f +%5.1f -%5.1f        %10.0f  %10.0f\n",
	      limitsEt[i], limitsEt[i+1],
	      limitsEta[j], limitsEta[j+1],
	      efficiency*100, efficiencyErrHi*100, efficiencyErrLo*100,
	      probesPass, probesFail);
      effOutput << strOut;
      effArray2D(i,j) = efficiency;
      effArrayErrLow2D(i,j) = efficiencyErrLo;
      effArrayErrHigh2D(i,j) = efficiencyErrHi;
    }
  }
  effOutput << endl;

  if(resultsRootFile && resultsRootFile->IsOpen()){
    resultsRootFile->cd();
    effArray2D.Write("effArray2D");
    effArrayErrLow2D.Write("effArrayErrLow2D");
    effArrayErrHigh2D.Write("effArrayErrHigh2D");    
    resultsRootFile->Close();
  }else assert(0);

  return;
}

void fitMass(TTree *passTree, TTree *failTree, 
	     TString cut, int mode,
	     double &efficiency, double &efficiencyErrHi, double &efficiencyErrLo,
	     TPad *passPad, TPad *failPad, ofstream &fitLog){
  
  // Define data sets
  RooRealVar mass("mass","mass",60, 120);
  mass.setBins(120);
  RooRealVar et ("et" ,"et" ,10.0, 1000);
  RooRealVar eta("eta","eta",-10, 10);
  RooFormulaVar rooCut("rooCut","rooCut",cut,RooArgSet(et,eta));
  RooDataSet  *dataUnbinnedPass = new RooDataSet("dataUnbinnedPass","dataUnbinnedPass",
						 passTree,RooArgSet(mass,et,eta), rooCut);
  RooDataSet  *dataUnbinnedFail = new RooDataSet("dataUnbinnedFail","dataUnbinnedFail",
						 failTree,RooArgSet(mass,et,eta), rooCut);
  RooDataHist *dataBinnedPass   = dataUnbinnedPass->binnedClone("dataBinnedPass","dataBinnedPass");
  RooDataHist *dataBinnedFail   = dataUnbinnedFail->binnedClone("dataBinnedFail","dataBinnedFail");
  RooCategory probeType("probeType","probeType");
  probeType.defineType("pass");
  probeType.defineType("fail");
  RooAbsData *data;

  // If needed do binned fit
  bool unbinnedFit = true;
  if(unbinnedFit){
    data = new RooDataSet("data","data",mass,Index(probeType),
			  Import("pass",*dataUnbinnedPass), 
			  Import("fail",*dataUnbinnedFail));
    cout << endl << "Setting up UNBINNED fit" << endl << endl;
  }else{
    data = new RooDataHist("data","data",mass,Index(probeType),
			   Import("pass",*dataBinnedPass), 
			   Import("fail",*dataBinnedFail));    
    cout << endl << "Setting up BINNED fit" << endl << endl;
  }
  // Define the PDFs
  //
  // Common pieces
  //
  // True signal model
  RooRealVar zMass ("zMass" ,"zMass" ,91.188);
  RooRealVar zWidth("zWidth","zWidth",2.495);
  RooBreitWigner bwPdf("bwPdf","bwPdf",mass,zMass,zWidth);
  // Total signal events and efficiency
  RooRealVar nsignal("nsignal","nsignal",1000,0.0,1.0e7);
  RooRealVar eff    ("eff"    ,"eff"    ,0.7,0.0,1.0);

  //  PDF for the PASS sample
  // 
  // Background
  RooRealVar lambdaBgPass("lambdaBgPass","lambdaBgPass",-0.1, -0.5, 0.0);
  RooExponential bgPassPdf("bgPassPdf","bgPassPdf",mass,lambdaBgPass);
  // Signal
  //     - resolution function
  RooRealVar cbMeanPass("cbMeanPass","cbMeanPass"   ,0.0, -10.0,10.0);
  RooRealVar cbWidthPass("cbWidthPass","cbWidthPass",1.0,  0.1, 5.0);
  RooRealVar cbAlphaPass("cbAlphaPass","cbAlphaPass",5.0,  0.0,20.0);
  RooRealVar cbNPass("cbNPass","cbNPass"            ,1.0,  0.0,10.0);
  RooCBShape cbPassPdf("cbPassPdf","cbPassPdf",mass,cbMeanPass,cbWidthPass,cbAlphaPass,cbNPass);
  //     - realistic model
  mass.setBins(10000,"fft");
  RooFFTConvPdf signalPassPdf("signalPassPdf","signalPassPdf",mass, bwPdf, cbPassPdf);
  // Combine signal and background
  RooFormulaVar nsigPass("nsigPass","nsigPass","@0*@1",RooArgList(nsignal,eff));
  RooRealVar nbgPass ("nbgPass" ,"nbgPass" ,1,0.0,1.0e5);
  RooAbsPdf *passPdf;
  if( mode == COUNTnFIT ){
    RooGenericPdf * simpleSignal = new RooGenericPdf("simpleSignal","simpleSignal",
						     "1.0",RooArgList());
    RooExtendPdf * simpleSignalExtended = new RooExtendPdf("passPdf", "passPdf",
							   *simpleSignal, nsigPass);
    passPdf = simpleSignalExtended;
  }else if( mode == FITnFIT ){
    passPdf = new RooAddPdf("passPdf","passPdf",RooArgList(signalPassPdf,bgPassPdf), RooArgList(nsigPass,nbgPass));
  }else{
    printf("ERROR: inappropriate mode requested\n");
    return;
  }
  //
  //  PDF for the FAIL sample
  // 
  // Background
  RooRealVar lambdaBgFail("lambdaBgFail","lambdaBgFail",-0.1, -0.5, 0.0);
  RooExponential bgFailPdf("bgFailPdf","bgFailPdf",mass,lambdaBgFail);
  // Signal
  //     - resolution function
  RooRealVar cbMeanFail("cbMeanFail","cbMeanFail"   ,0.0, -10.0,10.0);
  RooRealVar cbWidthFail("cbWidthFail","cbWidthFail",1.0,   0.1, 5.0);
  RooRealVar cbAlphaFail("cbAlphaFail","cbAlphaFail",5.0,   0.0,20.0);
  RooRealVar cbNFail("cbNFail","cbNFail"            ,1.0,   0.0,10.0);
  RooCBShape cbFailPdf("cbFailPdf","cbFailPdf",mass,cbMeanFail,cbWidthFail,cbAlphaFail,cbNFail);
  //     - realistic model
  RooFFTConvPdf signalFailPdf("signalFailPdf","signalFailPdf",mass, bwPdf, cbFailPdf);
  // Combine signal and background
  RooFormulaVar nsigFail("nsigFail","nsigFail","@0*(1.0-@1)",RooArgList(nsignal,eff));
  RooRealVar nbgFail ("nbgFail" ,"nbgFail" ,1,0.0,1.0e5);
  RooAddPdf failPdf("failPdf","failPdf",RooArgList(signalFailPdf,bgFailPdf), RooArgList(nsigFail,nbgFail));
  
  // Combine pass and fail
  RooSimultaneous fullPdf("fullPdf","fullPdf",probeType);
  fullPdf.addPdf(*passPdf,"pass");
  fullPdf.addPdf(failPdf,"fail");

  
  // Do the fit
  // Start with a reasonable point and do rough approximation first
  double total = dataUnbinnedPass->numEntries() + dataUnbinnedFail->numEntries();
  nsignal.setVal( 0.99*total);
  eff.setVal(0.90);
  if( mode == FITnFIT ){
    nbgPass.setVal(0.01*total);
    cbAlphaPass.setVal(1.0);
    cbNPass    .setVal(5.0);
    cbAlphaPass.setConstant(kTRUE);
    cbNPass    .setConstant(kTRUE);
  }
  nbgFail.setVal(0.01*total);
  cbAlphaFail.setVal(0.5);
  cbNFail    .setVal(5.0);
  cbAlphaFail.setConstant(kTRUE);
  cbNFail    .setConstant(kTRUE);
  RooFitResult *result = fullPdf.fitTo(*data,
				       Extended(kTRUE),
				       Save());

  // Release shape parameters and refine the fit
  if( mode == FITnFIT ){
    cbAlphaPass.setConstant(kFALSE);
    cbNPass    .setConstant(kFALSE);
  }
  cbAlphaFail.setConstant(kFALSE);
  cbNFail    .setConstant(kFALSE);
  result = fullPdf.fitTo(*data,
			 Extended(kTRUE),
			 Minos(RooArgSet(eff)),
			 Save());
  // If minos fails, refit without minos
  if((fabs(eff.getErrorLo())<5e-5) || (eff.getErrorHi()<5e-5)){
    cout << "MINOS FAILS" << endl;
    result = fullPdf.fitTo(*data,
			   Extended(kTRUE),
			   Save());
  }

  efficiency       = eff.getVal();
  efficiencyErrHi  = eff.getErrorHi();
  efficiencyErrLo  = fabs(eff.getErrorLo());

  // Draw fit results
  passPad->cd();
  passPad->Clear();
  RooPlot *framePass = mass.frame();
  dataUnbinnedPass->plotOn(framePass);
  if(mode == FITnFIT){
    passPdf->plotOn(framePass);
    passPdf->plotOn(framePass,Components("bgPassPdf"),LineStyle(kDashed));
  }
  framePass->Draw();
  passPad->Update();

  failPad->cd();
  failPad->Clear();
  RooPlot *frameFail = mass.frame();
  dataUnbinnedFail->plotOn(frameFail);
  failPdf.plotOn(frameFail);
  failPdf.plotOn(frameFail,Components("bgFailPdf"),LineStyle(kDashed));
  frameFail->Draw();
  failPad->Update();

  // Print fit outcome into fit log
  result->printStream(fitLog,RooPrintable::kValue,RooPrintable::kVerbose);
  fitLog << endl;
  printCorrelations(fitLog, result);
  
  return;
}

void printCorrelations(ostream& os, RooFitResult *res)
{
  ios_base::fmtflags flags = os.flags();
  const RooArgList *parlist = res->correlation("eff");
  
  os << "  Correlation Matrix" << endl;
  os << " --------------------" << endl;
  for(Int_t i=0; i<parlist->getSize(); i++) {
    for(Int_t j=0; j<parlist->getSize(); j++) 
      os << "  " << setw(7) << setprecision(4) << fixed << res->correlationMatrix()(i,j);    
    os << endl;
  }
  os.flags(flags);
}

TString getLabel(int sample, int effType, int method, 
		 int etBinning, int etaBinning){

  TString label = "";

  if(sample == DATA)
    label += "data";
  else if(sample == MC)
    label += "mc";
  else
    assert(0);

  if( effType == ID )
    label += "_id";
  else if(effType == HLT )
    label += "_hlt";
  else
    assert(0);

  if(method == COUNTnCOUNT)
    label += "_count-count";
  else if( method == COUNTnFIT ) 
    label += "_count-fit";
  else if( method == FITnFIT ) 
    label += "_fit-fit";
  else
    assert(0);

  label += "_bins-et";
  label += getNEtBins(etBinning);
  label += "-eta";
  label += getNEtaBins(etaBinning);

  return label;
}

int getTemplateBin(int etBin, int etaBin, int etaBinning){

  int templateBin = -1;

  if( etBin != -1 && etaBin != -1)
    templateBin = etBin * getNEtaBins(etaBinning) + etaBin;

  return templateBin;

}

TH1F * getPassTemplate(int etBin, int etaBin, int etaBinning, TFile *file){

  TH1F *hist = 0;
  if(file == 0)
    return hist;

  int templateBin = getTemplateBin(etBin, etaBin, etaBinning);
  if( templateBin == -1 )
    return hist;

  TString name = "hMassTemplate_Et";
  name += etBin;
  name += "_eta";
  name += etaBin;
  name += "_pass";

  hist = (TH1F*)file->Get(name);
  return hist;
}

TH1F * getFailTemplate(int etBin, int etaBin, int etaBinning, TFile *file){

  TH1F *hist = 0;
  if(file == 0)
    return hist;

  int templateBin = getTemplateBin(etBin, etaBin, etaBinning);
  if( templateBin == -1 )
    return hist;

  TString name = "hMassTemplate_Et";
  name += etBin;
  name += "_eta";
  name += etaBin;
  name += "_fail";

  hist = (TH1F*)file->Get(name);
  return hist;
}


// Alternative fit model
void fitMassWithTemplates(TTree *passTree, TTree *failTree, 
			  TString cut, int mode,
			  double &efficiency, double &efficiencyErrHi, double &efficiencyErrLo, 
			  TPad *passPad, TPad *failPad, ofstream &fitLog,
			  TH1F *templatePass, TH1F *templateFail){
  
  // Define data sets
  RooRealVar mass("mass","mass",60, 120);
  mass.setBins(30);
  RooRealVar et ("et" ,"et" ,10.0, 1000);
  RooRealVar eta("eta","eta",-10, 10);
  RooFormulaVar rooCut("rooCut","rooCut",cut,RooArgSet(et,eta));
  RooDataSet  *dataUnbinnedPass = new RooDataSet("dataUnbinnedPass","dataUnbinnedPass",
						 passTree,RooArgSet(mass,et,eta), rooCut);
  RooDataSet  *dataUnbinnedFail = new RooDataSet("dataUnbinnedFail","dataUnbinnedFail",
						 failTree,RooArgSet(mass,et,eta), rooCut);
  RooDataHist *dataBinnedPass   = dataUnbinnedPass->binnedClone("dataBinnedPass","dataBinnedPass");
  RooDataHist *dataBinnedFail   = dataUnbinnedFail->binnedClone("dataBinnedFail","dataBinnedFail");
  RooCategory probeType("probeType","probeType");
  probeType.defineType("pass");
  probeType.defineType("fail");
  RooAbsData *data;

  // If needed do binned fit
  bool unbinnedFit = true;
  if(unbinnedFit){
    data = new RooDataSet("data","data",mass,Index(probeType),
			  Import("pass",*dataUnbinnedPass), 
			  Import("fail",*dataUnbinnedFail));
    cout << endl << "Setting up UNBINNED fit" << endl << endl;
  }else{
    data = new RooDataHist("data","data",mass,Index(probeType),
			   Import("pass",*dataBinnedPass), 
			   Import("fail",*dataBinnedFail));    
    cout << endl << "Setting up BINNED fit" << endl << endl;
  }
  // Define the PDFs
  //
  // Common pieces
  //
  // Total signal events and efficiency
  RooRealVar nsignal("nsignal","nsignal",1000,0.0,1.0e7);
  RooRealVar eff    ("eff"    ,"eff"    ,0.7,0.0,1.0);

  //  PDF for the PASS sample
  // 
  // Background
  RooRealVar lambdaBgPass("lambdaBgPass","lambdaBgPass",-0.1, -0.5, 0.0);
  RooExponential bgPassPdf("bgPassPdf","bgPassPdf",mass,lambdaBgPass);
  // Signal
  //     - resolution function
  RooRealVar resMeanPass("resMeanPass","cbMeanPass"   ,0.0, -5.0, 5.0);
  RooRealVar resSigma   ("resSigma"   ,"resSigma  "   ,1.0, 0.1, 6.0);
  RooGaussian resPassPdf("resPassPdf","resPassPdf", mass, resMeanPass, resSigma);
//   RooRealVar cbMeanPass("cbMeanPass","cbMeanPass"   ,0.0, -5.0, 5.0);
//   RooRealVar cbWidthPass("cbWidthPass","cbWidthPass",1.0,  0.1, 6.0);
//   RooRealVar cbAlphaPass("cbAlphaPass","cbAlphaPass",5.0,  0.0,20.0);
//   RooRealVar cbNPass("cbNPass","cbNPass"            ,1.0,  0.0,10.0);
//   RooCBShape cbPassPdf("cbPassPdf","cbPassPdf",mass,cbMeanPass,cbWidthPass,cbAlphaPass,cbNPass);
  //      - mc template
  RooDataHist rooTemplatePass("rooTemplatePass","rooTemplatePass",RooArgList(mass),templatePass);
  RooHistPdf templatePassPdf("templatePassPdf","templatePassPdf",RooArgSet(mass),rooTemplatePass);
  //     - realistic model
  mass.setBins(10000,"fft");
//   RooFFTConvPdf signalPassPdf("signalPassPdf","signalPassPdf",mass, templatePassPdf, cbPassPdf);
  RooFFTConvPdf signalPassPdf("signalPassPdf","signalPassPdf",mass, templatePassPdf, resPassPdf);
  // Combine signal and background
  RooFormulaVar nsigPass("nsigPass","nsigPass","@0*@1",RooArgList(nsignal,eff));
  RooRealVar nbgPass ("nbgPass" ,"nbgPass" ,1,0.0,1.0e5);
  RooAbsPdf *passPdf;
  if( mode == COUNTnFIT ){
    RooGenericPdf * simpleSignal = new RooGenericPdf("simpleSignal","simpleSignal",
						     "1.0",RooArgList());
    RooExtendPdf * simpleSignalExtended = new RooExtendPdf("passPdf", "passPdf",
							   *simpleSignal, nsigPass);
    passPdf = simpleSignalExtended;
  }else if( mode == FITnFIT ){
    passPdf = new RooAddPdf("passPdf","passPdf",RooArgList(signalPassPdf,bgPassPdf), RooArgList(nsigPass,nbgPass));
  }else{
    printf("ERROR: inappropriate mode requested\n");
    return;
  }
  //
  //  PDF for the FAIL sample
  // 
  // Background
  RooRealVar lambdaBgFail("lambdaBgFail","lambdaBgFail",-0.1, -0.5, 0.0);
  RooExponential bgFailPdf("bgFailPdf","bgFailPdf",mass,lambdaBgFail);
  // Signal
  //     - resolution function
  RooRealVar resMeanFail("resMeanFail","cbMeanFail"   ,0.0, -5.0, 5.0);
  RooGaussian resFailPdf("resFailPdf","resFailPdf", mass, resMeanFail, resSigma);
//   RooRealVar cbMeanFail("cbMeanFail","cbMeanFail"   ,0.0,  -5.0, 5.0);
//   RooRealVar cbWidthFail("cbWidthFail","cbWidthFail",1.0,   0.1, 6.0);
//   RooRealVar cbAlphaFail("cbAlphaFail","cbAlphaFail",5.0,   0.0,20.0);
//   RooRealVar cbNFail("cbNFail","cbNFail"            ,1.0,   0.0,10.0);
//   RooCBShape cbFailPdf("cbFailPdf","cbFailPdf",mass,cbMeanFail,cbWidthFail,cbAlphaFail,cbNFail);
  //      - mc template
  RooDataHist rooTemplateFail("rooTemplateFail","rooTemplateFail",RooArgList(mass),templateFail);
  RooHistPdf templateFailPdf("templateFailPdf","templateFailPdf",RooArgSet(mass),rooTemplateFail);
  //     - realistic model
//   RooFFTConvPdf signalFailPdf("signalFailPdf","signalFailPdf",mass, templateFailPdf, cbFailPdf);
  RooFFTConvPdf signalFailPdf("signalFailPdf","signalFailPdf",mass, templateFailPdf, resFailPdf);
  // Combine signal and background
  RooFormulaVar nsigFail("nsigFail","nsigFail","@0*(1.0-@1)",RooArgList(nsignal,eff));
  RooRealVar nbgFail ("nbgFail" ,"nbgFail" ,1,0.0,1.0e5);
  RooAddPdf failPdf("failPdf","failPdf",RooArgList(signalFailPdf,bgFailPdf), RooArgList(nsigFail,nbgFail));
  
  // Combine pass and fail
  RooSimultaneous fullPdf("fullPdf","fullPdf",probeType);
  fullPdf.addPdf(*passPdf,"pass");
  fullPdf.addPdf(failPdf,"fail");

  
  // Do the fit
  double total = dataUnbinnedPass->numEntries() + dataUnbinnedFail->numEntries();
  nsignal.setVal( 0.99*total);
  eff.setVal(0.90);
  if( mode == FITnFIT ){
    nbgPass.setVal(0.01*total);
  }
  nbgFail.setVal(0.01*total);
  RooFitResult *result = fullPdf.fitTo(*data,
				       Extended(kTRUE),
				       Minos(RooArgSet(eff)),
				       Save());
  // If minos fails, refit without minos
  if((fabs(eff.getErrorLo())<5e-5) || (eff.getErrorHi()<5e-5)){
    cout << "MINOS FAILS" << endl;
    result = fullPdf.fitTo(*data,
			   Extended(kTRUE),
			   Save());
  }

  efficiency     = eff.getVal();
  efficiencyErrHi  = eff.getErrorHi();
  efficiencyErrLo  = fabs(eff.getErrorLo());

  // Draw fit results
  passPad->cd();
  passPad->Clear();
  RooPlot *framePass = mass.frame();
  dataUnbinnedPass->plotOn(framePass);
  if(mode == FITnFIT){
    passPdf->plotOn(framePass);
    passPdf->plotOn(framePass,Components("bgPassPdf"),LineStyle(kDashed));
  }
  framePass->Draw();
  passPad->Update();

  failPad->cd();
  failPad->Clear();
  RooPlot *frameFail = mass.frame();
  dataUnbinnedFail->plotOn(frameFail);
  failPdf.plotOn(frameFail);
  failPdf.plotOn(frameFail,Components("bgFailPdf"),LineStyle(kDashed));
  frameFail->Draw();
  failPad->Update();

  // Print fit outcome into fit log
  result->printStream(fitLog,RooPrintable::kValue,RooPrintable::kVerbose);
  fitLog << endl;
  printCorrelations(fitLog, result);
  
  return;
}

