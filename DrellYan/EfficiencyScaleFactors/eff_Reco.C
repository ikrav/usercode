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

#include "../Include/DYTools.hh"

#include "../Include/EleIDCuts.hh"

#endif

using namespace mithep;


//=== COMMON CONSTANTS ===========================================================================================

//=== FUNCTION DECLARATIONS ======================================================================================

Bool_t dielectronMatchedToGeneratorLevel(const TGenInfo *gen, const TDielectron *dielectron);

Bool_t electronMatchedToGeneratorLevel(const TGenInfo *gen, const TElectron *electron);

Bool_t scMatchedToGeneratorLevel(const TGenInfo *gen, const TPhoton *sc);

void measurePassAndFail(double &signal, double &signalErr, double &efficiency, double &efficiencyErr,
			TTree *passTree, TTree *failTree,
			TCanvas *passCanvas, TCanvas *failCanvas);
void measureEfficiency(TTree *passTree, TTree *failTree, 
		       int method, int etBinning, int etaBinning, 
		       TCanvas *canvas, ofstream &effOutput, ofstream &fitLog,
		       bool useTemplates, TFile *templatesFile, TFile *resultsRootFile);
void measureEfficiencyCountAndCount(TTree *passTree, TTree *failTree, 
				    int etBinning, int etaBinning, 
				    TCanvas *canvas, ofstream &effOutput, bool saveResultsToRootFile, TFile *resultsRootFile);
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

void eff_Reco(const TString configFile) 
{

  gBenchmark->Start("eff_Reco");
  
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
  if(effTypeString == "GSF")
    effType = GSF;
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
    // For data, we will be using templates
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
  if(!selectedEventsFile) 
    assert(0);

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
    
    // Read input file
    cout << "Processing " << ntupleFileNames[ifile] << "..." << endl;
    infile = new TFile(ntupleFileNames[ifile]); 
    assert(infile);
    
    // Get the TTrees
    eventTree = (TTree*)infile->Get("Events"); assert(eventTree);

    // Set branch address to structures that will store the info  
    eventTree->SetBranchAddress("Info",&info);                TBranch *infoBr       = eventTree->GetBranch("Info");
    eventTree->SetBranchAddress("Photon"  ,&scArr); 
    eventTree->SetBranchAddress("Electron",&eleArr); 
    TBranch *electronBr   = eventTree->GetBranch("Electron");
    TBranch *photonBr     = eventTree->GetBranch("Photon");
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
      eleArr->Clear();
      electronBr->GetEntry(ientry);
      scArr->Clear();
      photonBr->GetEntry(ientry);
      
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

      if(!(info->triggerBits & eventTriggerBit)) continue;  // no trigger accept? Skip to next event... 
      eventsAfterTrigger++;
      
      // Loop over the tag electrons
      for(int iele = 0; iele < eleArr->GetEntriesFast(); iele++){
	const TElectron *electron = (TElectron*)((*eleArr)[iele]);
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
	  
	  // total probes
	  hMassTotal->Fill(mass);
	  storeMass = mass;
	  storeEt   = sc->scEt;
	  storeEta  = sc->scEta;
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

  gBenchmark->Show("eff_Reco");
  
  
}


//=== FUNCTION DEFINITIONS ======================================================================================

Bool_t dielectronMatchedToGeneratorLevel(const TGenInfo *gen, const TDielectron *dielectron){

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

Bool_t electronMatchedToGeneratorLevel(const TGenInfo *gen, const TElectron *electron){

  Bool_t result = kTRUE;
  // In the generator branch of this ntuple, first particle is always
  // negative, and second always positive. 
  double dR=999;
  TLorentzVector vreco, v1gen, v2gen;
  vreco.SetPtEtaPhiM(electron->pt, electron->eta, electron->phi, 0.000511);
  v1gen .SetPtEtaPhiM(gen->pt_1, gen->eta_1, gen->phi_1, 0.000511);
  v2gen .SetPtEtaPhiM(gen->pt_2, gen->eta_2, gen->phi_2, 0.000511);
  if( electron->q < 0 ){
    dR = vreco.DeltaR(v1gen);
  }else{
    dR = vreco.DeltaR(v2gen);
  }
  // Require that both are within loose dR of 0.4, otherwise bail out
  if( fabs(dR) > 0.4 ) result = kFALSE; 
  
  return result;
}

Bool_t scMatchedToGeneratorLevel(const TGenInfo *gen, const TPhoton *sc){

  Bool_t result = kTRUE;
  // We do not know which of the gen electrons possibly
  // produced this supercluster, so we check both.
  double dR1=999, dR2=999;
  TLorentzVector vreco, v1gen, v2gen;
  vreco.SetPtEtaPhiM(sc->pt, sc->eta, sc->phi, 0.000511);
  v1gen .SetPtEtaPhiM(gen->pt_1, gen->eta_1, gen->phi_1, 0.000511);
  v2gen .SetPtEtaPhiM(gen->pt_2, gen->eta_2, gen->phi_2, 0.000511);
  dR1 = vreco.DeltaR(v1gen);
  dR2 = vreco.DeltaR(v2gen);
  // Require that at least one is  within loose dR of 0.4, otherwise bail out
  if( fabs(dR1) > 0.4 && fabs(dR2) > 0.4 ) result = kFALSE; 
  
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
  mass.setBins(10000,"cache");
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
  effOutput << "     SC ET         SC eta           efficiency             pass         fail\n";
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
      double efficiency, efficiencyErrHigh, efficiencyErrLow;
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
	fitMass(passTree, failTree, cut, method, 
		efficiency, efficiencyErrHigh, efficiencyErrLow,
		passPad, failPad, fitLog);
      }else{
	printf("\nMASS TEMPLATES ARE USED IN THE FIT\n\n");
	fitMassWithTemplates(passTree, failTree, cut, method, 
			     efficiency, efficiencyErrHigh, efficiencyErrLow,
		passPad, failPad, fitLog, templatePass, templateFail);
      }
      sprintf(strOut, "   %3.0f - %3.0f   %5.3f - %5.3f   %5.1f +%5.1f -%5.1f        %10.0f  %10.0f\n",
	      limitsEt[i], limitsEt[i+1],
	      limitsEta[j], limitsEta[j+1],
	      efficiency*100, efficiencyErrHigh*100, efficiencyErrLow*100,
	      probesPass, probesFail);
//       sprintf(strOut, "   %3.0f - %3.0f   %5.3f - %5.3f   %5.1f +-%5.1f           %10.0f  %10.0f\n",
// 	      limitsEt[i], limitsEt[i+1],
// 	      limitsEta[j], limitsEta[j+1],
// 	      efficiency*100, efficiencyErr*100,
// 	      probesPass, probesFail);
      effOutput << strOut;
      effArray2D(i,j) = efficiency;
      effArrayErrLow2D(i,j) = efficiencyErrLow;
      effArrayErrHigh2D(i,j) = efficiencyErrHigh;
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
	     double &efficiency, double &efficiencyErrHigh, double &efficiencyErrLow,  
	     TPad *passPad, TPad *failPad, ofstream &fitLog){
  
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
  RooRealVar cbMeanPass("cbMeanPass","cbMeanPass"   ,0.0, -5.0, 5.0);
  RooRealVar cbWidthPass("cbWidthPass","cbWidthPass",1.0,  0.1, 6.0);
  RooRealVar cbAlphaPass("cbAlphaPass","cbAlphaPass",5.0,  0.0,20.0);
  RooRealVar cbNPass("cbNPass","cbNPass"            ,1.0,  0.0,10.0);
  RooCBShape cbPassPdf("cbPassPdf","cbPassPdf",mass,cbMeanPass,cbWidthPass,cbAlphaPass,cbNPass);
  //     - realistic model
  mass.setBins(10000,"cache");
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
  RooRealVar cbMeanFail("cbMeanFail","cbMeanFail"   ,0.0,  -5.0, 5.0);
  RooRealVar cbWidthFail("cbWidthFail","cbWidthFail",1.0,   0.1, 6.0);
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

  efficiency     = eff.getVal();
  efficiencyErrHigh  = eff.getErrorHi();
  efficiencyErrLow   = fabs(eff.getErrorLo());

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

  if( effType == GSF )
    label += "_gsf";
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
			  double &efficiency, double &efficiencyErrHigh, double &efficiencyErrLow, 
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
  mass.setBins(10000,"cache");
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
  // The limits for the "fail" come from stating at the fit results without
  // splitting into Et bins. In some cases, the peak is not there at all.
  RooRealVar resMeanFail("resMeanFail","cbMeanFail"   ,0.0, -3.0, 3.0);
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
  efficiencyErrHigh  = eff.getErrorHi();
  efficiencyErrLow   = fabs(eff.getErrorLo());
//   efficiencyErr  = eff.getError();

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

