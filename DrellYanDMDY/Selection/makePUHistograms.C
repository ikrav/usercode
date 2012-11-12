//================================================================================================
//
// Z->e e selection macro
//
//  * plots distributions associated with selected events
//  * prints list of selected events from data
//  * outputs ROOT files of events passing selection for each sample, 
//    which can be processed by plotSelect.C
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TBranch.h>                // class to access branches in TTree
#include <TClonesArray.h>           // ROOT array class
#include <TCanvas.h>                // class for drawing
#include <TH1F.h>                   // 1D histograms
#include <TH2D.h>
#include <TBenchmark.h>             // class to track macro running statistics
#include <TLorentzVector.h>         // 4-vector class
#include <TVector3.h>               // 3D vector class
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings

using namespace std;

#include "../Include/CPlot.hh"          // helper class for plots
#include "../Include/MitStyleRemix.hh"  // style settings for drawing
#include "../Include/MyTools.hh"        // miscellaneous helper functions
#include "../Include/CSample.hh"        // helper class for organizing input ntuple files
#include "../Include/DYTools.hh"
#include "../Include/DYToolsUI.hh"
#include "../Include/TriggerSelection.hh"

// define structures to read in ntuple
#include "../Include/EWKAnaDefs.hh"
#include "../Include/TEventInfo.hh"
#include "../Include/TVertex.hh"

// lumi section selection with JSON files
#include "../Include/JsonParser.hh"

#include "../Include/EventSelector.hh"

#endif

// define structure for output ntuple
#include "../Include/ZeeData.hh"


//=== MAIN MACRO =================================================================================================

void makePUHistograms(const TString conf, 
		      const TString triggerSetString="Full2011DatasetTriggers", 
		      int debugMode=0) 
{  
  gBenchmark->Start("makePUHistograms");

  // fast check
  TriggerConstantSet triggerSet=DetermineTriggerSet(triggerSetString);  
  assert ( triggerSet != TrigSet_UNDEFINED );

  // Construct the trigger object
  TriggerSelection requiredTriggers(triggerSetString, true, 0);
  assert(requiredTriggers.isDefined());

  if (debugMode) std::cout << "\n\n\tDEBUG MODE is ON\n\n";

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================
  
  TString  outputDir;         // output directory
  Double_t lumi;              // luminosity (pb^-1)
  Bool_t   doWeight;          // weight events?
  TString  escaleTag;         // Energy scale calibrations tag
  TString  format;            // plot format

  vector<TString>  snamev;    // sample name (for output file)  
  vector<CSample*> samplev;   // data/MC samples
  Bool_t hasData=false;
    
  //
  // parse .conf file
  //
  ifstream ifs;
  ifs.open(conf.Data());
  assert(ifs.is_open());
  string line;
  Int_t state=0;
  string generateEEMFile;
  int generateEEMFEWZFile=0;
  while(getline(ifs,line)) {
    if ((line[0]=='#') && (line[1]=='$') && (line[2]=='$')) {
      if (line.find("generate_EEM_files=") != string::npos) {
	generateEEMFile=line.substr(line.find('=')+1);
	generateEEMFEWZFile=(line.find("FEWZ") != string::npos) ? 1:0;
	std::cout << "\n\tEEM files will be generated, tag=<" 
		  << generateEEMFile 
		  << ">, with_FEWZ_weights=" << generateEEMFEWZFile << "\n\n";
	continue;
      }
    }
    if(line[0]=='#') continue;
    if(line[0]=='%') { 
      state++; 
      continue; 
    }
    if(line[0]=='$') {
      samplev.push_back(new CSample());
      stringstream ss(line);
      string chr;
      string sname;
      Int_t color;
      ss >> chr >> sname >> color;
      string label = line.substr(line.find('@')+1);
      snamev.push_back(sname);
      samplev.back()->label = label;
      samplev.back()->color = color;
      continue;
    }

    if(state==0) {  // general settings
      stringstream ss1(line); ss1 >> lumi;
      getline(ifs,line);
      stringstream ss2(line); ss2 >> doWeight;
      getline(ifs,line);
      outputDir = TString(line);
      getline(ifs,line);
      stringstream ss3(line); ss3 >> escaleTag;
      getline(ifs,line);
      format = TString(line);
      
    } else if(state==1) {  // define data sample
      string fname;
      Double_t xsec;
      string json;
      stringstream ss(line);
      ss >> fname >> xsec >> json;
      samplev.back()->fnamev.push_back(fname);
      samplev.back()->xsecv.push_back(xsec);
      samplev.back()->jsonv.push_back(json);
      hasData=true;
    
    } else if(state==2) {  // define MC samples
      string fname;
      Double_t xsec;
      stringstream ss(line);
      ss >> fname >> xsec;
      samplev.back()->fnamev.push_back(fname);
      samplev.back()->xsecv.push_back(xsec);
    }
  }
  ifs.close();


  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
  
  //
  // Set up histograms
  //
  vector<TH1F*> hNGoodPV_signal_trigger_v;
  vector<TH1F*> hNGoodPV_tnp_ele17_sc8_v;
  vector<TH1F*> hNGoodPV_tnp_ele17_ele8_v;
  vector<TH1F*> hNGoodPV_tnp_all_triggers_v;

  char hname[100];
  for(UInt_t isam=0; isam<samplev.size(); isam++) {

    sprintf(hname, "hNGoodPV_signal_trigger_v_%s", snamev[isam].Data());
    hNGoodPV_signal_trigger_v.push_back(new TH1F(hname,"",47,-0.5,46.5));        
    hNGoodPV_signal_trigger_v[isam]->Sumw2();

    sprintf(hname, "hNGoodPV_tnp_ele17_sc8_v_%s", snamev[isam].Data());
    hNGoodPV_tnp_ele17_sc8_v.push_back(new TH1F(hname,"",47,-0.5,46.5));        
    hNGoodPV_tnp_ele17_sc8_v[isam]->Sumw2();

    sprintf(hname, "hNGoodPV_tnp_ele17_ele8_v_%s", snamev[isam].Data());
    hNGoodPV_tnp_ele17_ele8_v.push_back(new TH1F(hname,"",47,-0.5,46.5));        
    hNGoodPV_tnp_ele17_ele8_v[isam]->Sumw2();

    sprintf(hname, "hNGoodPV_tnp_all_triggers_v_%s", snamev[isam].Data());
    hNGoodPV_tnp_all_triggers_v.push_back(new TH1F(hname,"",47,-0.5,46.5));        
    hNGoodPV_tnp_all_triggers_v[isam]->Sumw2();

  }

  //
  // Access samples and fill histograms
  //  
  TFile *infile=0;
  TTree *eventTree=0;  
  
  // Data structures to store info from TTrees
  mithep::TEventInfo *info    = new mithep::TEventInfo();
  TClonesArray *pvArr         = new TClonesArray("mithep::TVertex");
  
  //
  // loop over samples
  //
  for(UInt_t isam=0; isam<samplev.size(); isam++) {        
    
    //
    // loop through files
    //
    CSample* samp = samplev[isam];
    const UInt_t nfiles = samplev[isam]->fnamev.size();    
    for(UInt_t ifile=0; ifile<nfiles; ifile++) {
      cout << "Processing " << samp->fnamev[ifile] << "... "; cout.flush();
      infile = new TFile(samp->fnamev[ifile]);
      assert(infile);
    
      Bool_t hasJSON = kFALSE;
      // mithep::RunLumiRangeMap rlrm;
      JsonParser jsonParser;
      if((samp->jsonv.size()>0) && 
	 samp->jsonv[ifile].Length() &&
	 (samp->jsonv[ifile].CompareTo("NONE")!=0)) { 
        hasJSON = kTRUE;
	std::cout << "JSON file <" << samp->jsonv[ifile] << ">\n";
        jsonParser.Initialize(samp->jsonv[ifile].Data()); 
      }
      
      // Get the TTree
      eventTree = (TTree*)infile->Get("Events"); assert(eventTree);
      
      // Set branch address to structures that will store the info  
      eventTree->SetBranchAddress("Info",       &info);          TBranch *infoBr       = eventTree->GetBranch("Info");
      eventTree->SetBranchAddress("PV",         &pvArr);         TBranch *pvBr         = eventTree->GetBranch("PV");

      // Determine maximum number of events to consider
      // *** CASES ***
      // <> lumi < 0                             => use all events in the sample
      // <> xsec = 0                             => for data (use all events)
      // <> lumi > 0, xsec > 0, doWeight = true  => use all events and scale to lumi
      // <> lumi > 0, xsec > 0, doWeight = false => compute expected number of events
      UInt_t maxEvents = eventTree->GetEntries();
      Double_t weight = 1;
      if(lumi>0) {
        Double_t xsec = samp->xsecv[ifile];
	if(xsec>0) { 

	  if(doWeight) { weight = lumi*xsec/(Double_t)eventTree->GetEntries(); } 
	  else         { maxEvents = (UInt_t)(lumi*xsec); } 
	}       
      }  
      if(maxEvents > eventTree->GetEntries()) {
        cout << "Not enough events for " << lumi << " pb^-1 in file: " << samp->fnamev[ifile];
        return;
      }
      samp->weightv.push_back(weight);
     
      // loop through events
      std::cout << "numEntries = " << eventTree->GetEntries() << std::endl;
      for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
	if (debugMode && (ientry>10000)) break; // debug option
	if(ientry >= maxEvents) break;
	
	infoBr->GetEntry(ientry);

        if(hasJSON && !jsonParser.HasRunLumi(info->runNum, info->lumiSec)) continue;  // not certified run? Skip to next event...
	
	// Configure the object for trigger matching	
	bool isData = ((isam == 0) && hasData);
	requiredTriggers.actOnData(isData);
	ULong_t analysisTriggerBit = requiredTriggers.getEventTriggerBit(info->runNum);
	// We retrieve the trigger bit for the "id" (a bit more triggers than for "hlt")
	// indicated by the flag true in the call below.
	ULong_t tnpTriggerBit       = requiredTriggers.getEventTriggerBit_TagProbe(info->runNum, true);

	// Count number of good PVs
	pvArr->Clear();
	pvBr->GetEntry(ientry);
	UInt_t nGoodPV=0;
	nGoodPV=countGoodVertices(pvArr);
	
	// The primary signal trigger
	if( (info->triggerBits & analysisTriggerBit) ){
	  hNGoodPV_signal_trigger_v[isam]->Fill(nGoodPV, weight);
	}
	
	// The primary tag and probe trigger
	if( (info->triggerBits & tnpTriggerBit) ){
	  hNGoodPV_tnp_all_triggers_v[isam]->Fill(nGoodPV, weight);
	}
	
	// One of the two tag and probe triggers, the one for RECO, by itself
	if( (info->triggerBits & Triggers2011::kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30) ){
	  hNGoodPV_tnp_ele17_sc8_v[isam]->Fill(nGoodPV, weight);
	}
	
	// The other of the two tag and probe triggers, used for HLT and ID, by itself
	if( (info->triggerBits & Triggers2011::kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30) ){
	  hNGoodPV_tnp_ele17_ele8_v[isam]->Fill(nGoodPV, weight);
	}
	
      }
      delete infile;
      infile=0, eventTree=0;
    }
    std::cout << "next sample" << std::endl;
    
  }
  delete info;
  delete pvArr;

  TFile fout(TString("./referencePVHistogramByTrigger.root"),"recreate");
  for(UInt_t isam=0; isam<samplev.size(); isam++) {        
    hNGoodPV_signal_trigger_v[isam]->Write();
    hNGoodPV_tnp_ele17_sc8_v [isam]->Write();
    hNGoodPV_tnp_ele17_ele8_v[isam]->Write();
    hNGoodPV_tnp_all_triggers_v[isam]->Write();
  }
  fout.Close();

  gBenchmark->Show("makePUHistograms");       
} 
