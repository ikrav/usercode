#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TClonesArray.h>
#include <TString.h>
#include <TBenchmark.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

// define structures to read in ntuple
#include "../Include/EWKAnaDefs.hh"
#include "../Include/TEventInfo.hh"
#include "../Include/TGenInfo.hh"
#include "../Include/TMuon.hh"
#include "../Include/TElectron.hh"
#include "../Include/TDielectron.hh"
#include "../Include/TPhoton.hh"
#include "../Include/TJet.hh"
#include "../Include/TVertex.hh"

#include "../Include/DYTools.hh"
#include "../Include/EleIDCuts.hh"
#endif

// Main macro function
//--------------------------------------------------------------------------------------------------
void SkimNtuples(const TString input = "skim.input") 
{
  gBenchmark->Start("SkimNtuples");
  
  TString outfilename;          // output of skimming 
  vector<TString> infilenames;  // list input ntuple files to be skimmed
  TString sample;
  
  // 
  // parse input file
  //  
  ifstream ifs;
  ifs.open(input.Data()); 
  assert(ifs.is_open());
  string line;
  // First line should be DATA or SIGNALMC or BGMC
  getline(ifs,line); 
  sample = line;
  // Second line is the OUTPUT skim file name
  getline(ifs,line); 
  outfilename = line;
  // All subsequent lines are names of INPUT root files
  while(getline(ifs,line)) { infilenames.push_back(line); }
  ifs.close();
  
  bool isGenPresent = true;
  if( sample == "DATA" || sample == "BGMC")
    isGenPresent = false;
  else if( sample == "SIGNALMC" )
    isGenPresent = true;
  else{
    printf("Unknown sample type: use DATA or SIGNALMC or BGMC only in the input configuration file.\n");
    return;
  }

  TTree::SetMaxTreeSize(kMaxLong64);
  
  // Don't write TObject part of the objects
  mithep::TEventInfo::Class()->IgnoreTObjectStreamer();
  mithep::TGenInfo::Class()->IgnoreTObjectStreamer();
  mithep::TElectron::Class()->IgnoreTObjectStreamer();
  mithep::TDielectron::Class()->IgnoreTObjectStreamer();
  mithep::TMuon::Class()->IgnoreTObjectStreamer();
  mithep::TJet::Class()->IgnoreTObjectStreamer();
  mithep::TPhoton::Class()->IgnoreTObjectStreamer();
  mithep::TVertex::Class()->IgnoreTObjectStreamer();
  
  // Data structures to store info from TTrees
  mithep::TEventInfo *info    = new mithep::TEventInfo();
  mithep::TGenInfo *gen    = new mithep::TGenInfo();
  TClonesArray *electronArr   = new TClonesArray("mithep::TElectron");
  TClonesArray *dielectronArr = new TClonesArray("mithep::TDielectron");
  TClonesArray *muonArr       = new TClonesArray("mithep::TMuon");
  TClonesArray *pfJetArr      = new TClonesArray("mithep::TJet");
  TClonesArray *photonArr     = new TClonesArray("mithep::TPhoton");
  TClonesArray *pvArr         = new TClonesArray("mithep::TVertex");
  
  UInt_t nInputEvts = 0;
  UInt_t nPassEvts  = 0;
  
  TFile* outfile = new TFile(outfilename, "RECREATE");
  
  //
  // Initialize data trees and structs
  // 
  TTree *outEventTree = new TTree("Events","Events"); 
  outEventTree->Branch("Info",       &info);
  if( isGenPresent )
    outEventTree->Branch("Gen",       &gen);
  outEventTree->Branch("Electron",   &electronArr);
  outEventTree->Branch("Dielectron", &dielectronArr);
  outEventTree->Branch("Muon",       &muonArr);
  outEventTree->Branch("PFJet",      &pfJetArr);
  outEventTree->Branch("Photon",     &photonArr);
  outEventTree->Branch("PV",         &pvArr);

  for(UInt_t ifile=0; ifile<infilenames.size(); ifile++) {
    cout << "Skimming " << infilenames[ifile] << "..." << endl;
    TFile *infile = new TFile(infilenames[ifile]);
    assert(infile);
    
    TTree *eventTree = (TTree*)infile->Get("Events");
    assert(eventTree);
    
    // Set branch address to structures that will store the info  
    eventTree->SetBranchAddress("Info",       &info);          TBranch *infoBr       = eventTree->GetBranch("Info");
    TBranch *genBr = 0;
    if(isGenPresent){
      eventTree->SetBranchAddress("Gen" ,       &gen);
      genBr        = eventTree->GetBranch("Gen");
      if( !genBr ){
	printf("MC info is not found in signal MC file\n");
	assert(0);
      }
    }
    eventTree->SetBranchAddress("Electron",   &electronArr);   TBranch *electronBr   = eventTree->GetBranch("Electron");
    eventTree->SetBranchAddress("Dielectron", &dielectronArr); TBranch *dielectronBr = eventTree->GetBranch("Dielectron");
    eventTree->SetBranchAddress("Muon",       &muonArr);       TBranch *muonBr       = eventTree->GetBranch("Muon");
    eventTree->SetBranchAddress("PFJet",      &pfJetArr);      TBranch *pfJetBr      = eventTree->GetBranch("PFJet");
    eventTree->SetBranchAddress("Photon",     &photonArr);     TBranch *photonBr     = eventTree->GetBranch("Photon");
    eventTree->SetBranchAddress("PV",         &pvArr);         TBranch *pvBr         = eventTree->GetBranch("PV");
    
    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) { 
//     for(UInt_t ientry=0; ientry< 100000; ientry++) { // For testing
      infoBr->GetEntry(ientry);
      if( isGenPresent)
	genBr->GetEntry(ientry);

      electronArr->Clear();   
      dielectronArr->Clear(); 
      muonArr->Clear();       
      pfJetArr->Clear();      
      photonArr->Clear();     
      pvArr->Clear();         
      
      nInputEvts++;
            
      Bool_t keep = kFALSE;

      dielectronBr->GetEntry(ientry);

      // Require at least one dielectron
      if(dielectronArr->GetEntriesFast() > 0) {
	// Apply cuts on the dielectron
	int nDielectronsPass = 0;
	for(Int_t i=0; i<dielectronArr->GetEntriesFast(); i++) {
	  mithep::TDielectron *dielectron = (mithep::TDielectron*)((*dielectronArr)[i]);
	  // Require at least one dielectron above 19 GeV and the other above 9 GeV
	  bool etCut = false;
	  if( (dielectron->scEt_1 > 19 && dielectron->scEt_2 > 9) ||
	      (dielectron->scEt_1 > 9 && dielectron->scEt_2 > 19) )
	    etCut = true;
	  // Require at least one dielectron to pass full ID
	  bool idCut = false;
	  if( passSmurf(extractElectron(dielectron,1)) ||
	      passSmurf(extractElectron(dielectron,2)) )
	    idCut = true;
	  if( etCut && idCut )
	    nDielectronsPass++;
	}
	if(nDielectronsPass > 0)
	  keep = kTRUE;
      }
      
      if(keep) {
	// Some of the objects are dropped for skimmed events
// 	electronBr->GetEntry(ientry);
// 	muonBr->GetEntry(ientry);
// 	pfJetBr->GetEntry(ientry);
// 	photonBr->GetEntry(ientry);
// Fill only those branches that we want to write out
	pvBr->GetEntry(ientry);
        nPassEvts++;
      }else{
// Clear branches which we filled to make the pass/fail check,
// but do not want to write out for failed events
	dielectronArr->Clear();
      }

      outEventTree->Fill();

    }
  }
  
  outfile->Write();
  outfile->Close();
  
  delete info;
  delete electronArr;
  delete dielectronArr;
  delete muonArr;
  delete pfJetArr;
  delete photonArr;
  delete pvArr;
    
  std::cout << outfilename << " created!" << std::endl;
  std::cout << " >>> Events processed: " << nInputEvts << std::endl;
  std::cout << " >>>   Events passing: " << nPassEvts << std::endl;
  
  gBenchmark->Show("SkimNtuples");
}  
