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
#include <sstream>

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
void TrimNtuples(const TString input = "trim.input") 
{
  gBenchmark->Start("TrimNtuples");
  
  TString outfilename;          // output of skimming 
  vector<TString> infilenames;  // list input ntuple files to be skimmed
  double massMin,massMax;
  
  // 
  // parse input file
  //  
  ifstream ifs;
  ifs.open(input.Data()); 
  assert(ifs.is_open());
  string line;
  // First line should contain two numbers: min and max mass of virtual Z
  getline(ifs,line); 
  stringstream ss(line);
  ss >> massMin >> massMax;
  // Second line is the OUTPUT skim file name
  getline(ifs,line); 
  outfilename = line;
  // All subsequent lines are names of INPUT root files
  while(getline(ifs,line)) { infilenames.push_back(line); }
  ifs.close();
  
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
    eventTree->SetBranchAddress("Gen" ,       &gen);           TBranch *genBr        = eventTree->GetBranch("Gen");
    eventTree->SetBranchAddress("Electron",   &electronArr);   TBranch *electronBr   = eventTree->GetBranch("Electron");
    eventTree->SetBranchAddress("Dielectron", &dielectronArr); TBranch *dielectronBr = eventTree->GetBranch("Dielectron");
    eventTree->SetBranchAddress("Muon",       &muonArr);       TBranch *muonBr       = eventTree->GetBranch("Muon");
    eventTree->SetBranchAddress("PFJet",      &pfJetArr);      TBranch *pfJetBr      = eventTree->GetBranch("PFJet");
    eventTree->SetBranchAddress("Photon",     &photonArr);     TBranch *photonBr     = eventTree->GetBranch("Photon");
    eventTree->SetBranchAddress("PV",         &pvArr);         TBranch *pvBr         = eventTree->GetBranch("PV");
    
    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) { 
//     for(UInt_t ientry=0; ientry< 100000; ientry++) { // For testing
      infoBr->GetEntry(ientry);
      genBr->GetEntry(ientry);

      electronArr->Clear();   
      dielectronArr->Clear(); 
      muonArr->Clear();       
      pfJetArr->Clear();      
      photonArr->Clear();     
      pvArr->Clear();         
      
      nInputEvts++;
            
      Bool_t keep = kFALSE;

      if(gen->vmass >= massMin && gen->vmass < massMax )
	keep = kTRUE;

      if(keep) {
	dielectronBr->GetEntry(ientry);
 	electronBr->GetEntry(ientry);
 	muonBr->GetEntry(ientry);
 	pfJetBr->GetEntry(ientry);
 	photonBr->GetEntry(ientry);
	pvBr->GetEntry(ientry);
	outEventTree->Fill();
        nPassEvts++;
      }

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
  
  gBenchmark->Show("TrimNtuples");
}  
