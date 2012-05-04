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
#include "RooGaussian.h"
#include "RooGenericPdf.h"
#include "RooHistPdf.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooFFTConvPdf.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooDataSet.h"
#include "RooFitResult.h"

#include "../Include/CPlot.hh"          // helper class for plots
#include "../Include/MitStyleRemix.hh"  // style settings for drawing

// define classes and constants to read in ntuple
#include "../Include/EWKAnaDefs.hh"
#include "../Include/TGenInfo.hh"
#include "../Include/TEventInfo.hh"
#include "../Include/TDielectron.hh"
#include "../Include/TElectron.hh"
#include "../Include/TPhoton.hh"
#include "../Include/TVertex.hh"
#include "../Include/DYTools.hh"
#include "../Include/EleIDCuts.hh"
#include "../Include/TriggerSelection.hh"

#include "../Include/cutFunctions.hh"
#include "../Include/fitFunctions.hh"
#include "../Include/fitFunctionsCore.hh"
#include "../EventScaleFactors/tnpSelectEvents.hh"

#include "../Include/EventSelector.hh"

#endif

using namespace mithep;


//=== COMMON CONSTANTS ===========================================================================================


//=== FUNCTION DECLARATIONS ======================================================================================

//=== MAIN MACRO =================================================================================================

void eff_Reco(const TString configFile, const TString effTypeString, const TString triggerSetString, int debugMode=0) 
{

  //  ---------------------------------
  //       Preliminary checks
  //  ---------------------------------

  // verify whether it was a compilation check
  if (configFile.Contains("_DebugRun_") || triggerSetString.Contains("_DebugRun_")) {
    std::cout << "eff_Reco: _DebugRun_ detected. Terminating the script\n";
    return;
  }

  if (!effTypeString.Contains("RECO")) {
    std::cout << "eff_Reco: effTypeString should be \"RECO\"\n";
    return;
  }

  // fast check
  // Construct the trigger object
  TriggerSelection triggers(triggerSetString, true, 0); // later calls actOnData
  assert ( triggers.isDefined() );

  if (debugMode) std::cout << "\n\n\tDEBUG MODE is ON\n\n";

  //  ---------------------------------
  //         Normal execution
  //  ---------------------------------

  using namespace mithep; 
 
  gBenchmark->Start("eff_Reco");
  

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================
  
  Double_t massLow  = 60;
  Double_t massHigh = 120;

  // Read in the configuratoin file
  TString sampleTypeString = "";
  TString calcMethodString = "";
  TString etBinningString  = "";
  TString etaBinningString = "";
  TString dirTag;
  vector<TString> ntupleFileNames;
  ifstream ifs;
  ifs.open(configFile.Data());
  if (!ifs.is_open()) {
    std::cout << "configFile=" << configFile << "\n";
    assert(ifs.is_open());
  }
  string line;
  Int_t state=0;
  Int_t subState=0;
  while(getline(ifs,line)) {
    if(line[0]=='#') continue;
    if(line[0]=='%') break;
    if(state==0){
      // Read 1st line of content: data or MC?
      sampleTypeString = TString(line);
      state++;
    }else if(state==1) {
      // Read 2d content line: efficiency fitting mode
      size_t pos=line.find(':');
      if (pos==string::npos) {
	std::cout << "expected format is EFFICIENCY:fitting_mode\n";
	return;
      }
      subState++;
      if (line.find(effTypeString)!=string::npos) {
	calcMethodString = TString(line.c_str()+pos+1);
      }
      if (subState==3) state++;
    }else if(state==2) {
      // Read 3rd content line: SC ET binning
      etBinningString = TString(line);
      state++;
    }else if(state==3) {
      // Read 4th content line: SC eta binning
      etaBinningString = TString(line);
      state++;
    }else if(state==4) {
      // Read 5th content line: directory tag
      dirTag = TString(line);
      state++;
    }else if(state==5) {
      ntupleFileNames.push_back(TString(line));
    }
  }
  ifs.close();

  if (state!=5) {
    std::cout << "failed to read input file\n";
    return;
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
  if(effTypeString == "RECO")
    effType = RECO;
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
  else if(etaBinningString == "ETABINS5")
    etaBinning = ETABINS5;
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

  // Correct the trigger object
  triggers.actOnData((sample==DATA)?true:false);

  // The label is a string that contains the fields that are passed to
  // the function below, to be used to name files with the output later.
  TString label = getLabel(sample, effType, calcMethod, etBinning, etaBinning, triggers);

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================
  
#ifdef tnpSelectEventsIsObject
  tnpSelectEvent_t::Class()->IgnoreTObjectStreamer();
#endif

  //  
  // Set up histograms
  //
  //   TH1F* hMass           = new TH1F("hMass"     ,"",30,massLow,massHigh);
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
    // For data, we will be using templates
    // however, if the request is COUNTnCOUNT, do nothing
    if( calcMethod != COUNTnCOUNT ){
      TString labelMC = getLabel(-1111, effType, calcMethod, etBinning, etaBinning, triggers);
      TString templatesLabel = tagAndProbeDir+TString("/mass_templates_")+labelMC+TString(".root");
      templatesFile = new TFile(templatesLabel);
      if( ! templatesFile->IsOpen() ) {
	std::cout << "templatesFile name " << templatesLabel << "\n";
	assert(0);
      }
    }
  }

  // This file is utilized by fit_EffReco
  TString uScore="_";
  TString selectEventsFName=tagAndProbeDir + TString("/selectEvents_") 
    + analysisTag + uScore
    + sampleTypeString + uScore +
    + effTypeString + uScore +  triggers.triggerSetName() + TString(".root");
  std::cout << "selectEventsFName=<" << selectEventsFName << ">\n"; 
  TFile *selectedEventsFile = new TFile(selectEventsFName,"recreate");
  if(!selectedEventsFile) 
    assert(0);

  tnpSelectEvent_t storeData;
  const int new_store_data_code=1;
  TTree *passTree = new TTree("passTree","passTree");
  assert(passTree);
  Double_t storeMass, storeEt, storeEta;
  UInt_t storeNGoodPV;
  if (new_store_data_code) {
    storeData.createBranches(passTree);
  }
  else {
    passTree->Branch("mass",&storeMass,"mass/D");
    passTree->Branch("et",&storeEt  ,"et/D");
    passTree->Branch("eta",&storeEta ,"eta/D");
    passTree->Branch("nGoodPV",&storeNGoodPV,"nGoodPV/i");
  }
  

  TTree *failTree = new TTree("failTree","failTree");
  assert(failTree);
  if (new_store_data_code) {
    storeData.createBranches(failTree);
  }
  else {
    failTree->Branch("mass",&storeMass,"mass/D");
    failTree->Branch("et",&storeEt  ,"et/D");
    failTree->Branch("eta",&storeEta ,"eta/D");
    failTree->Branch("nGoodPV",&storeNGoodPV,"nGoodPV/i");
  }


  int eventsInNtuple = 0;
  int eventsAfterTrigger = 0;
  int tagCand = 0;
  int tagCandPassEt = 0;
  int tagCandPassEta = 0;
  int tagCandGenMatched = 0;
  int tagCandEcalDriven = 0;
  int tagCandFinalCount = 0;
  int numTagProbePairs = 0;
  int numTagProbePairsPassEt = 0;
  int numTagProbePairsPassEta = 0;
  int numTagProbePairsGenMatched = 0;
  int numTagProbePairsInMassWindow = 0;
  
  // Loop over files
  for(UInt_t ifile=0; ifile<ntupleFileNames.size(); ifile++){

    //
    // Access samples and fill histograms
    //  
    TFile *infile = 0;
    TTree *eventTree = 0;
        
    // Data structures to store info from TTrees
    mithep::TEventInfo *info = new mithep::TEventInfo();
    mithep::TGenInfo      *gen  = new mithep::TGenInfo();
    TClonesArray *scArr   = new TClonesArray("mithep::TPhoton");
    TClonesArray *eleArr  = new TClonesArray("mithep::TElectron");
    TClonesArray *pvArr   = new TClonesArray("mithep::TVertex");
    
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
    eventTree->SetBranchAddress("Photon"  ,&scArr); 
    eventTree->SetBranchAddress("Electron",&eleArr); 
    eventTree->SetBranchAddress("PV", &pvArr); 
    TBranch *electronBr   = eventTree->GetBranch("Electron");
    TBranch *photonBr     = eventTree->GetBranch("Photon");
    TBranch *pvBr         = eventTree->GetBranch("PV");
    assert(electronBr); assert(photonBr);
    assert(pvBr);
    TBranch *genBr = 0;
    if(sample != DATA){
      eventTree->SetBranchAddress("Gen",&gen);
      genBr = eventTree->GetBranch("Gen");
    }

    // loop over events    
    eventsInNtuple += eventTree->GetEntries();
    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
      //for(UInt_t ientry=0; ientry<1000; ientry++) { 
      if (debugMode && (ientry>100000)) break;  // This is for faster turn-around in testing
      
      if(sample != DATA)
	genBr->GetEntry(ientry);
      eleArr->Clear();
      electronBr->GetEntry(ientry);
      scArr->Clear();
      photonBr->GetEntry(ientry);
      
      // Check that the whole event has fired the appropriate trigger
      infoBr->GetEntry(ientry);

      /*  Old code
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
      */
      ULong_t eventTriggerBit= triggers.getEventTriggerBit_SCtoGSF(info->runNum);
      ULong_t tagTriggerObjectBit= triggers.getLeadingTriggerObjectBit_SCtoGSF(info->runNum);

      if(!(info->triggerBits & eventTriggerBit)) continue;  // no trigger accept? Skip to next event... 
      eventsAfterTrigger++;
      
      // Loop over the tag electrons
      for(int iele = 0; iele < eleArr->GetEntriesFast(); iele++){
	
	const mithep::TElectron *electron = (mithep::TElectron*)((*eleArr)[iele]);
	tagCand++;

	// All cuts for the tag electron should be applied here
	if(electron->scEt<20) continue;
	tagCandPassEt++;

	bool isBele = isBarrel(electron->scEta);
	bool isEele = isEndcap(electron->scEta);
	if( ! isBele && ! isEele) continue;
	tagCandPassEta++;
	
	if( sample != DATA)
	  if( ! electronMatchedToGeneratorLevel(gen, electron) ) continue;
	
	tagCandGenMatched++;

	// ECAL driven: this condition is NOT applied	
	
	if( !isTag( electron, tagTriggerObjectBit) ) continue;

	tagCandFinalCount++;
	
	
	// Loop over superclusters in this event: the probes
	// Note: each supercluster has a number assigned: scID,
	// and each electron object from TElectron collection has
	// the field scID that tells from which supercluster this electron
	// comes from. That allows to make the match between the 
	// object in TPhoton collection and TElectron collection to
	// find superclusters reconstructed as electrons.
	
	for(int isc = 0; isc < scArr->GetEntriesFast(); isc++){
	  
	  const TPhoton *sc = (TPhoton*)((*scArr)[isc]);
	  // Avoid probe that is same as tag
	  if( sc->scID == electron->scID ) continue;

	  numTagProbePairs++;
	  // Apply probe cuts
	  if(sc->scEt < 10) continue;
	  numTagProbePairsPassEt++;

	  bool isBsc = isBarrel(sc->scEta);
	  bool isEsc = isEndcap(sc->scEta);
	  if( ! isBsc && ! isEsc) continue;
	  numTagProbePairsPassEta++;

	  if( sample != DATA)
	    if( ! scMatchedToGeneratorLevel(gen, sc) ) continue;
	  numTagProbePairsGenMatched++;
	  
	  // Find mass of the electron-supercluster pair
	  TLorentzVector ele4V, sc4V, dycand4V;
	  ele4V.SetPtEtaPhiM(electron->pt, electron->eta, electron->phi, 0.000511);
	  sc4V .SetPtEtaPhiM(sc->pt, sc->eta, sc->phi, 0.000511);
	  dycand4V = ele4V + sc4V;
	  double mass = dycand4V.M();	  
	  // Tag and probe is done around the Z peak
	  if((mass < massLow) || (mass > massHigh)) continue;
	  numTagProbePairsInMassWindow++;

	  // The probes are fully selected at this point.

	  // Loop over electron collection again to find a match to this supercluster
	  // Match only to ECAL-driven GSF electrons
	  const TElectron *electronMatch = 0;
	  for(int iele2 = 0; iele2 < eleArr->GetEntriesFast(); iele2++){
	    const TElectron *electron2 = (TElectron*)((*eleArr)[iele2]);
	    if( sc->scID == electron2->scID )
	      if( electron2->isEcalDriven )
		electronMatch = electron2;
	  } // end loop over electrons searching for SC match
	  
	  // get the number of goodPVs
	  pvBr->GetEntry(ientry);
	  storeNGoodPV=0;
	  for(Int_t ipv=0; ipv<pvArr->GetEntriesFast(); ipv++) {
	    const mithep::TVertex *pv = (mithep::TVertex*)((*pvArr)[ipv]);
	    if(pv->nTracksFit                        < 1)  continue;
	    if(pv->ndof                              < 4)  continue;
	    if(fabs(pv->z)                           > 24) continue;
	    if(sqrt((pv->x)*(pv->x)+(pv->y)*(pv->y)) > 2)  continue;
	    storeNGoodPV++;
	  }

	  // total probes
	  hMassTotal->Fill(mass);
	  storeMass = mass;
	  storeEt   = sc->scEt;
	  storeEta  = sc->scEta;
	  if (new_store_data_code) {
	    storeData.assign(mass,sc->scEt,sc->scEta,storeNGoodPV);
	  }
	  int templateBin = getTemplateBin( findEtBin(sc->scEt,etBinning),
					    findEtaBin(sc->scEta,etaBinning),
					    etaBinning);
	  if( electronMatch != 0 ){
	    // supercluster has match in reconstructed electrons: "pass"
	    hMassPass->Fill(mass);
	    passTree->Fill();
	    if(sample != DATA && templateBin != -1)
	      hPassTemplateV[templateBin]->Fill(mass);
	  }else{
	    // supercluster is not reconstructed as an electron
	    hMassFail->Fill(mass);
	    failTree->Fill();
	    if(sample != DATA && templateBin != -1)
	      hFailTemplateV[templateBin]->Fill(mass);
	  }
	  
	  } // end loop over superclusters - probes	  
      } // end loop over electrons - tags      
    } // end loop over events
    
    delete infile;
    infile=0;
    eventTree=0;
    
    delete gen;
    delete info;
    delete eleArr;
    delete scArr;
  } // end loop over files

  // save the selected trees
  selectedEventsFile->cd();
  passTree->Write();
  failTree->Write();
  selectedEventsFile->Write();

  //
  // Efficiency analysis
  //
  
  //   printf("Number of regular candidates:      %15.0f\n", hMass->GetSumOfWeights());
  printf("Total events in ntuple                                       %15d\n",eventsInNtuple);
  printf("    events after event level trigger cut                     %15d\n",eventsAfterTrigger);
  printf("\nTotal electron tag candidates (no cuts)                      %15d\n",tagCand);
  printf("                 tag candidates Et>20                        %15d\n",tagCandPassEt);
  printf("                 tag candidates, eta in acceptance           %15d\n",tagCandPassEta);
  printf("                 tag candidates, matched to GEN (if MC)      %15d\n",tagCandGenMatched);
  printf("                 tag candidates, ECAL driven                 %15d\n",tagCandEcalDriven);
  printf("                 tag candidates, full selection(ID,HLT)      %15d\n",tagCandFinalCount);

  printf("\nTotal tag(electron)-probe(supercluster) pairs                %15d\n",numTagProbePairs);
  printf("               probe Et>10                                   %15d\n",numTagProbePairsPassEt);
  printf("               probe eta in acceptance                       %15d\n",numTagProbePairsPassEta);
  printf("               probe matched to GEN (if MC)                  %15d\n",numTagProbePairsGenMatched);
  printf("               tag-probe mass in 60-120 GeV window           %15d\n",numTagProbePairsInMassWindow);

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
  TString fitlogname = TString("results_unsorted/efficiency_TnP_")+label+TString("_fitlog.dat");
  ofstream fitLog;
  fitLog.open(fitlogname);

  //
  //  Find efficiency
  //
  bool useTemplates = false;
  if(sample == DATA)
    useTemplates = true;

  int NsetBins=30;
  bool isRECO=1;
  const char* setBinsType="cache";
  
  int nDivisions = getNEtBins(etBinning)*getNEtaBins(etaBinning);
  double ymax = 800;
  if(nDivisions <4 )
    ymax = nDivisions * 200;
  TCanvas *c1 = MakeCanvas("c1","c1", 600, (int)ymax);
  c1->Divide(2,nDivisions);

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
  //c1->Update();
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

    
  selectedEventsFile->Close();
  std::cout << "selectedEventsFile <" << selectEventsFName << "> saved\n";
 
  gBenchmark->Show("eff_Reco");
  
  
}


