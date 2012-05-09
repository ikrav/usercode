//================================================================================================
//
// applyJSONfilter.C
//
//  * processes the provided configuration data files, producing 
//    JSON-filtered files
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TBenchmark.h>             // class to track macro running statistics
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TBranch.h>                // class to access branches in TTree
#include <TClonesArray.h>           // ROOT array class

/*
#include <TCanvas.h>                // class for drawing
#include <TH1F.h>                   // 1D histograms
#include <TH2D.h>
#include <TLorentzVector.h>         // 4-vector class
#include <TVector3.h>               // 3D vector class
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings
*/


// define structures to read in ntuple
#include "../Include/TEventInfo.hh"
#include "../Include/TGenInfo.hh"
#include "../Include/TDielectron.hh"
#include "../Include/TJet.hh"
#include "../Include/TVertex.hh"
#include "../Include/TPhoton.hh"
#include "../Include/TMuon.hh"

// lumi section selection with JSON files
#include "../Include/TriggerSelection.hh"
#include "../Include/JsonParser.hh"

// input file processor 
#include "../Include/InputFileMgr.hh"

#endif


//=== MAIN MACRO =================================================================================================

void applyJSONandTriggerFilter(const TString &conf, std::string triggerSelectionString="Full2011", int applyJSON=0)
{  
  gBenchmark->Start("applyJSONandTriggerFilter");

  
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================

  // Defined triggers
  TriggerSelection requiredTriggers(triggerSelectionString, true, 0);
  ULong_t allEventTriggerBits=requiredTriggers.getCombinedEventTriggerBit();
  std::string evtTrigStr = PrintBitsToStdString(allEventTriggerBits,"_evtTrig_",1);
  std::cout << " ** all trigger bits=<" << evtTrigStr << ">\n";
  std::cout << " ** applyJSON=" << applyJSON << "\n";
 
  // load input file
  InitialInputMgr_t mgr;
  if (!mgr.Load(conf)) {
    std::cout << "failed to load input file <" << conf << ">\n";
    return;
  }

  Bool_t hasData = (mgr.sampleInfo(0)->fnamev.size()>0);
  if (!hasData) {
    std::cout << "the provided configuration file has no data\n";
    return;
  }

  //--------------------------------------------------------------------------------------------------------------
  // Main filtering code 
  //==============================================================================================================  
  
  const TString ntupDir = mgr.outputDir() + TString("/jsonTrigNtuples"); gSystem->mkdir(ntupDir,kTRUE);
 

  TTree::SetMaxTreeSize(kMaxLong64);
 
  //
  // Don't write TObject part of the objects
  //
  TDescriptiveInfo_t::Class()->IgnoreTObjectStreamer();
  mithep::TEventInfo::Class()->IgnoreTObjectStreamer();
  mithep::TGenInfo::Class()->IgnoreTObjectStreamer();
  mithep::TElectron::Class()->IgnoreTObjectStreamer();
  mithep::TDielectron::Class()->IgnoreTObjectStreamer();
  mithep::TMuon::Class()->IgnoreTObjectStreamer();
  mithep::TJet::Class()->IgnoreTObjectStreamer();
  mithep::TPhoton::Class()->IgnoreTObjectStreamer();
  mithep::TVertex::Class()->IgnoreTObjectStreamer();

 
  //
  // Access samples and fill histograms
  //  
  TFile *infile=0;
  TTree *eventTree=0;  
  
  // Data structures to store info from TTrees
  mithep::TEventInfo *info    = new mithep::TEventInfo();
  mithep::TGenInfo *gen       = new mithep::TGenInfo();
  // TClonesArray *caloJetArr    = new TClonesArray("mithep::TJet");
  // TClonesArray *trackJetArr   = new TClonesArray("mithep::TJet");
  TClonesArray *pfJetArr      = new TClonesArray("mithep::TJet");
  TClonesArray *dielectronArr = new TClonesArray("mithep::TDielectron");
  TClonesArray *electronArr = new TClonesArray("mithep::TElectron");
  TClonesArray *photonArr= new TClonesArray("mithep::TPhoton");
  TClonesArray *pvArr         = new TClonesArray("mithep::TVertex");
  UInt_t origNumEntries = 0;
  
  TDescriptiveInfo_t *description= new TDescriptiveInfo_t();

  //
  // loop over samples
  //
  for(UInt_t isam=0; isam<mgr.sampleCount(); isam++) {
    //if (isam==0) continue;
    //if (isam!=mgr.sampleCount()-1) continue;
    if (isam>0) continue;

    const CSample* samp = mgr.sampleInfo(isam);
    const UInt_t nfiles = samp->fnamev.size();    
    for(UInt_t ifile=0; ifile<nfiles; ifile++) {
      cout << "Processing " << samp->fnamev[ifile] << "... "; cout.flush();

      //if (ifile>0) break; // debug

      infile = new TFile(samp->fnamev[ifile]);
      assert(infile);

      Bool_t hasJSON = kFALSE;
      // mithep::RunLumiRangeMap rlrm;
      JsonParser jsonParser;
      if((samp->jsonv.size()>0) && (samp->jsonv[ifile].CompareTo("NONE")!=0)) { 
        hasJSON = kTRUE;
	// rlrm.AddJSONFile(samp->jsonv[ifile].Data()); 
	std::cout << "JSON file " << samp->jsonv[ifile] << "\n";
        jsonParser.Initialize(samp->jsonv[ifile].Data()); 
      }
      //if (!hasJSON) {
      //std::cout << "this script cannot apply JSON filter, if no JSON file is provided\n";
      //return;
      //}   
      

      // Get the TTree
      eventTree = (TTree*)infile->Get("Events"); assert(eventTree);

      // Set branch address to structures that will store the info  
      eventTree->SetBranchAddress("Info",       &info);   
      TBranch *infoBr       = eventTree->GetBranch("Info");
      eventTree->SetBranchAddress("Dielectron", &dielectronArr); 
      //TBranch *dielectronBr = eventTree->GetBranch("Dielectron");
      eventTree->SetBranchAddress("Electron", &electronArr); 
      eventTree->SetBranchAddress("Photon", &photonArr);
      // eventTree->SetBranchAddress("CaloJet",    &caloJetArr);    TBranch *caloJetBr    = eventTree->GetBranch("CaloJet");
      // eventTree->SetBranchAddress("TrackJet",   &trackJetArr);   TBranch *trackJetBr   = eventTree->GetBranch("TrackJet");
      eventTree->SetBranchAddress("PFJet",      &pfJetArr);
      //      TBranch *pfJetBr      = eventTree->GetBranch("PFJet");
      eventTree->SetBranchAddress("PV",         &pvArr);
      //         TBranch *pvBr         = eventTree->GetBranch("PV");

      // Generator information is present only for signal MC, which should be the last entry
      TBranch *genBr = 0;
      if( isam == mgr.sampleCount()-1 ){  
	eventTree->SetBranchAddress("Gen",&gen);
	genBr = eventTree->GetBranch("Gen");
      }
      

      // Prepare a filtered file 
      TString okFileName=samp->fnamev[ifile];
      std::string extension=evtTrigStr;
      if (applyJSON && (isam==0)) extension += "_json";
      extension+=".root";
      okFileName.Replace(okFileName.Index(".root"),5,extension.c_str());
      std::cout << "output file name=<" << okFileName << ">\n";
      TFile *okFile = new TFile(okFileName,"RECREATE");
      TTree *descriptionTree = new TTree("Description","description");
      descriptionTree->Branch("description","TDescriptiveInfo_t",&description);
      descriptionTree->Branch("origNumEntries",&origNumEntries,"origNumEntries/i",sizeof(origNumEntries));
      TTree *okTree = new TTree("Events","Events");
      okTree->Branch("Info","mithep::TEventInfo",&info);
      okTree->Branch("Dielectron","TClonesArray(mithep::TDielectron)",&dielectronArr,320000,0);
      okTree->Branch("Electron","TClonesArray(mithep::TElectron)",&electronArr,320000,0);
      okTree->Branch("Photon","TClonesArray(mithep::TPhoton)",&photonArr,320000,0);
      okTree->Branch("PFJet","TClonesArray(mithep::TJet)",&pfJetArr,320000,0);
      okTree->Branch("PV", "TClonesArray(mithep::TVertex)",&pvArr,320000,0);

      if( isam == mgr.sampleCount()-1 ) okTree->Branch("Gen","mithep::TGenInfo",&gen);

      //
      //   Fill the description tree
      //
      {
	// get number of entries
	origNumEntries = eventTree->GetEntries();
	// describe
	std::vector<std::string> *lines= & description->_info;
	lines->clear();
	std::string s;
	s="# Initial file: "; s+=samp->fnamev[ifile].Data();
	lines->push_back(s);
	s="# Json file applied: "; 
	if (applyJSON && samp->jsonv[ifile].Length()) s+=samp->jsonv[ifile].Data(); else s+="(none)";
	lines->push_back(s);
	char buf[100];
	sprintf(buf,"# Event trigger bits applied (%lu): ",allEventTriggerBits); 
	s=buf; s+=evtTrigStr;
	lines->push_back(s);

	// fill
	descriptionTree->Fill();
      }
      


      // loop through events
      const UInt_t nEvents=eventTree->GetEntries();
      UInt_t nsel=0;
      std::cout << "nEvents=" << nEvents << std::endl;
      for(UInt_t ientry=0; ientry<nEvents; ientry++) {       
	
	infoBr->GetEntry(ientry);

	//
	//   JSON filtering
	//
        if(applyJSON && hasJSON && !jsonParser.HasRunLumi(info->runNum, info->lumiSec)) {
	  //std::cout << " JSON rejected the run\n";
	  continue;  // not certified run? Skip to next event...
	}
        
	//
	//   Event trigger filtering
	//
	if ((info->triggerBits & allEventTriggerBits)==0) {
	  continue;
	}
	
	//
	//  filters passed. Get all events and save to another file
	//
	eventTree->GetEntry(ientry);
	nsel++;
	okTree->Fill();
      }

      delete infile;
      infile=0; eventTree=0;
      okFile->Write();
      okFile=0;
      std::cout << "selected " << nsel << "/" << nEvents << " (" << 0.1*trunc(nsel*1000./double(nEvents)) << "\%) events\n";
      std::cout << "file " << okFileName << " created\n";
    }
  }

  gBenchmark->Show("applyJSONandTriggerFilter");
  return;
}
