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
#include <TRandom.h>
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
#include "../Include/TriggerSelection.hh"

#include "../Include/cutFunctions.hh"
#include "../Include/fitFunctions.hh"
#include "../Include/fitFunctionsCore.hh"

#endif

using namespace mithep;


//=== COMMON CONSTANTS ===========================================================================================

const Double_t kGAP_LOW  = 1.4442;
const Double_t kGAP_HIGH = 1.566;

//=== FUNCTION DECLARATIONS ======================================================================================

//=== MAIN MACRO =================================================================================================

void eff_IdHlt(const TString configFile, TString triggerSetString) 
{

  //  ---------------------------------
  //       Preliminary checks
  //  ---------------------------------

  // verify whether it was a compilation check
  if (configFile.Contains("_DebugRun_") || triggerSetString.Contains("_DebugRun_")) {
    std::cout << "eff_IdHlt: _DebugRun_ detected. Terminating the script\n";
    return;
  }

  // fast check
  TriggerConstantSet triggerSet=DetermineTriggerSet(triggerSetString);  
  assert ( triggerSet != TrigSet_UNDEFINED );

  //  ---------------------------------
  //         Normal execution
  //  ---------------------------------

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
  if (!ifs.is_open()) {
    std::cout << "tried to open the configuration file <" << configFile << ">\n";
    assert(ifs.is_open());
  }
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

  // Construct the trigger object
  TriggerSelection triggers(triggerSetString, (sample==DATA)?true:false, 0);
  if (effType==HLT) std::cout << "\tHLT efficiency calculation method " << triggers.hltEffCalcName() << ", triggerSet=" << triggers.triggerSetName() << "\n";
  else triggers.hltEffCalcMethod(HLTEffCalc_2011Old);

  TRandom *rnd= new TRandom();
  rnd->SetSeed(0); 

 // The label is a string that contains the fields that are passed to
  // the function below, to be used to name files with the output later.
  TString label = getLabel(sample, effType, calcMethod, etBinning, etaBinning, triggers);

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
    TString labelMC = getLabel(-1111, effType, calcMethod, etBinning, etaBinning, triggers);
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
      TString labelMC = getLabel(-1111, effType, calcMethod, etBinning, etaBinning, triggers);
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
  if (!selectedEventsFile) {
    assert(0);
  }

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
  TCanvas *c1 = MakeCanvas("c1","c1", 600, int(ymax));
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

    // check whether the file is suitable for the requested run range
    UInt_t runNumMin = UInt_t(eventTree->GetMinimum("runNum"));
    UInt_t runNumMax = UInt_t(eventTree->GetMaximum("runNum"));
    std::cout << "runNumMin=" << runNumMin << ", runNumMax=" << runNumMax << "\n";
    if (!triggers.validRunRange(runNumMin,runNumMax)) {
      std::cout << "... file contains uninteresting run range\n";
      continue;
    }

    // Define other branches
    eventTree->SetBranchAddress("Dielectron",&dielectronArr); TBranch *dielectronBr = eventTree->GetBranch("Dielectron");
    TBranch *genBr = 0;
    if(sample != DATA){
      eventTree->SetBranchAddress("Gen",&gen);
      genBr = eventTree->GetBranch("Gen");
    }

    // loop over events    
    eventsInNtuple += eventTree->GetEntries();
     for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
       //  for(UInt_t ientry=0; ientry<200000; ientry++) { // This is for faster turn-around in testing
       
      if(sample != DATA)
	genBr->GetEntry(ientry);
      
      // Check that the whole event has fired the appropriate trigger
      infoBr->GetEntry(ientry);
      
      /* Old code
      // For EPS2011 for both data and MC (starting from Summer11 production)
      // we use a special trigger for tag and probe that has second leg
      // unbiased with cuts at HLT
      ULong_t eventTriggerBit = kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30
	| kHLT_Ele32_CaloIdL_CaloIsoVL_SC17;
      // The tag trigger bit matches the "electron" of the trigger we
      // use for this tag and probe study: electron+sc
      ULong_t tagTriggerObjectBit = kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_EleObj
	| kHLT_Ele32_CaloIdL_CaloIsoVL_SC17_EleObj;
      // The probe trigger, however, is any of possibilities used in
      // the trigger that is used in the main analysis
      ULong_t probeTriggerObjectBit = kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele1Obj
	| kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele1Obj
	| kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele2Obj
	| kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele2Obj;
      */

      ULong_t eventTriggerBit= triggers.getEventTriggerBit_TagProbe(info->runNum);

      if(!(info->triggerBits & eventTriggerBit)) continue;  // no trigger accept? Skip to next event... 
      eventsAfterTrigger++;

      ULong_t tagTriggerObjectBit= triggers.getTagTriggerObjBit(info->runNum);
      ULong_t probeTriggerObjectBit_Tight= triggers.getProbeTriggerObjBit_Tight(info->runNum);
      ULong_t probeTriggerObjectBit_Loose= triggers.getProbeTriggerObjBit_Loose(info->runNum);
      ULong_t probeTriggerObjectBit= probeTriggerObjectBit_Tight | probeTriggerObjectBit_Loose;

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
	  if( ! dielectronMatchedToGeneratorLevel(gen, dielectron) ) continue;
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
	  if (triggers.useRandomTagTnPMethod(info->runNum)) {
	    if (rnd->Uniform() <= 0.5) {
	      // tag is 1st electron
	      if (!isTag1) continue;
	      isTag2=0; // ignore whether ele2 can be a tag
	    }
	    else {
	      if (!isTag2) continue;
	      isTag1=0; // ignore whether ele1 can be a tag
	    }
	  }
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

  int NsetBins=120;
  bool isRECO=0;
  const char* setBinsType="cache";

  measureEfficiency(passTree, failTree,
		    calcMethod, etBinning, etaBinning, c1, effOutput, fitLog,
		    useTemplates, templatesFile, resultsRootFile,
		    NsetBins, isRECO, setBinsType, 
		    dirTag, triggers.triggerSetName());

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


