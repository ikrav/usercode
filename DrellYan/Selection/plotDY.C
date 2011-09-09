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
#include <TBenchmark.h>             // class to track macro running statistics
#include <TLorentzVector.h>         // 4-vector class
#include <TVector3.h>               // 3D vector class
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings

#include "../Include/CPlot.hh"          // helper class for plots
#include "../Include/MitStyleRemix.hh"  // style settings for drawing
#include "../Include/MyTools.hh"        // miscellaneous helper functions
#include "../Include/CSample.hh"        // helper class for organizing input ntuple files

// define structures to read in ntuple
#include "../Include/EWKAnaDefs.hh"
#include "../Include/TEventInfo.hh"
#include "../Include/TDielectron.hh"
#include "../Include/TJet.hh"
#include "../Include/TVertex.hh"

// lumi section selection with JSON files
#include "../Include/JsonParser.hh"

// Helper functions for Electron ID selection
#include "../Include/EleIDCuts.hh" 

// define structure for output ntuple
#include "../Include/ZeeData.hh"

#include "../Include/ElectronEnergyScale.hh" //energy scale correction

#endif

//=== FUNCTION DECLARATIONS ======================================================================================

// generate web page
void makeHTML(const TString outDir);

// fill ntuple of selected events
void fillData(ZeeData *data, const mithep::TEventInfo *info, const mithep::TDielectron *dielectron,
              const UInt_t npv, const UInt_t njets, const Double_t weight);

// print event dump
void eventDump(ofstream &ofs, const mithep::TDielectron *dielectron, 
               const UInt_t runNum, const UInt_t lumiSec, const UInt_t evtNum, 
	       const UInt_t triggerObj1, const UInt_t triggerObj2);


//=== MAIN MACRO =================================================================================================

void plotDY(const TString conf) 
{  
  gBenchmark->Start("plotDY");

  
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================
  
  TString  outputDir;         // output directory
  Double_t lumi;              // luminosity (pb^-1)
  Bool_t   doWeight;          // weight events?
  TString  format;            // plot format

  vector<TString>  snamev;    // sample name (for output file)  
  vector<CSample*> samplev;   // data/MC samples
    
  //
  // parse .conf file
  //
  ifstream ifs;
  ifs.open(conf.Data());
  assert(ifs.is_open());
  string line;
  Int_t state=0;
  while(getline(ifs,line)) {
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

  // sOutDir is a static data member in the CPlot class.
  // There is a strange crash of the whole ROOT session well after
  // this script is executed when one attempts to exit ROOT, with 
  // a dump of memory map. This happens only on UNL Tier3, but
  // there is no sign of a problem on any other computer.
  //   The consequence of this variable is not set is that the plots
  // will be created in the local directory rather than the
  // one configured through sOutDir.
//   CPlot::sOutDir        = outputDir + TString("/plots");   gSystem->mkdir(CPlot::sOutDir,kTRUE);

  const TString ntupDir = outputDir + TString("/ntuples"); gSystem->mkdir(ntupDir,kTRUE);
  
  Bool_t hasData = (samplev[0]->fnamev.size()>0);

  const Double_t kGAP_LOW  = 1.4442;
  const Double_t kGAP_HIGH = 1.566;
    
  //
  // Canvas dimensions
  //
  Int_t canw=800, canh=600;

  
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
  
  //
  // Set up histograms
  //
  vector<TH1F*> hMassv, hMass2v, hMass3v;
  vector<TH1F*> hMassBBv, hMassBEv, hMassEEv;
  vector<TH1F*> hPtv, hPt2v, hyv, hPhiv;   
  
  vector<TH1F*> hElePtv, hEleEtav, hElePhiv;
  vector<TH1F*> hTrkIso03v, hEmIso03v, hHadIso03v;

  vector<TH1F*> hCaloMexv, hCaloMeyv, hCaloMetv, hCaloSumEtv, hCaloMetPhiv;
  vector<TH1F*> hTCMexv, hTCMeyv, hTCMetv, hTCSumEtv, hTCMetPhiv;
  vector<TH1F*> hPFMexv, hPFMeyv, hPFMetv, hPFSumEtv, hPFMetPhiv;

//   vector<TH1F*> hNCaloJetsv, hCaloJetEtv, hCaloJetEtav, hCaloJetPhiv; 
//   vector<TH1F*> hNTrackJetsv, hTrackJetEtv, hTrackJetEtav, hTrackJetPhiv; 
  vector<TH1F*> hNPFJetsv, hPFJetEtv, hPFJetEtav, hPFJetPhiv;
  vector<TH1F*> hNTracks0v;
  vector<TH1F*> hNCaloTowers0v;
  vector<TH1F*> hNPVv, hNGoodPVv, hGoodPVNTracksv, hGoodPVChi2v, hGoodPVNdofv, hGoodPVSumPtv;
  
  vector<TH1F*> hDeltaEtaInBv, hDeltaEtaInEv;
  vector<TH1F*> hDeltaPhiInBv, hDeltaPhiInEv;
  
  vector<Double_t> nSelv, nSelVarv;  
  vector<Double_t> nPosSSv;
  vector<Double_t> nNegSSv;
  
  TH1F* hGoodPVDz = new TH1F("hGoodPVDz","",200,-40,40);
  
  UInt_t nProcessedEvents=0;
  TH1F* hTrigger = new TH1F("hTrigger","",32,-0.5,31.5);
  
  char hname[100];
  for(UInt_t isam=0; isam<samplev.size(); isam++) {

    sprintf(hname,"hMass_%i",isam);   hMassv.push_back(new TH1F(hname,"",30,60,120));   hMassv[isam]->Sumw2();
    sprintf(hname,"hMass2_%i",isam);  hMass2v.push_back(new TH1F(hname,"",35,20,160));  hMass2v[isam]->Sumw2();
    sprintf(hname,"hMass3_%i",isam);  hMass3v.push_back(new TH1F(hname,"",50,0,500));   hMass3v[isam]->Sumw2();
    sprintf(hname,"hMassBB_%i",isam); hMassBBv.push_back(new TH1F(hname,"",30,60,120)); hMassBBv[isam]->Sumw2();
    sprintf(hname,"hMassBE_%i",isam); hMassBEv.push_back(new TH1F(hname,"",30,60,120)); hMassBEv[isam]->Sumw2();
    sprintf(hname,"hMassEE_%i",isam); hMassEEv.push_back(new TH1F(hname,"",30,60,120)); hMassEEv[isam]->Sumw2();
    sprintf(hname,"hPt_%i",isam);     hPtv.push_back(new TH1F(hname,"",50,0,50));       hPtv[isam]->Sumw2();
    sprintf(hname,"hPt2_%i",isam);    hPt2v.push_back(new TH1F(hname,"",50,0,200));     hPt2v[isam]->Sumw2();
    sprintf(hname,"hy_%i",isam);      hyv.push_back(new TH1F(hname,"",20,-3,3));        hyv[isam]->Sumw2();
    sprintf(hname,"hPhi_%i",isam);    hPhiv.push_back(new TH1F(hname,"",20,-3.2,3.2));  hPhiv[isam]->Sumw2();
    
    sprintf(hname,"hElePt_%i",isam);      hElePtv.push_back(new TH1F(hname,"",50,0,150));     hElePtv[isam]->Sumw2();
    sprintf(hname,"hEleEta_%i",isam);     hEleEtav.push_back(new TH1F(hname,"",20,-3,3));     hEleEtav[isam]->Sumw2();
    sprintf(hname,"hElePhi_%i",isam);     hElePhiv.push_back(new TH1F(hname,"",20,-3.2,3.2)); hElePhiv[isam]->Sumw2(); 
    sprintf(hname,"hTrkIso03_%i",isam);   hTrkIso03v.push_back(new TH1F(hname,"",60,0,3));    hTrkIso03v[isam]->Sumw2();
    sprintf(hname,"hEmIso03_%i",isam);    hEmIso03v.push_back(new TH1F(hname,"",60,0,3));     hEmIso03v[isam]->Sumw2();
    sprintf(hname,"hHadIso03_%i",isam);   hHadIso03v.push_back(new TH1F(hname,"",60,0,3));    hHadIso03v[isam]->Sumw2();

    sprintf(hname,"hCaloMex_%i",isam);    hCaloMexv.push_back(new TH1F(hname,"",40,-50,50));      hCaloMexv[isam]->Sumw2();
    sprintf(hname,"hCaloMey_%i",isam);    hCaloMeyv.push_back(new TH1F(hname,"",40,-50,50));      hCaloMeyv[isam]->Sumw2();
    sprintf(hname,"hCaloMet_%i",isam);    hCaloMetv.push_back(new TH1F(hname,"",20,0,60));        hCaloMetv[isam]->Sumw2();
    sprintf(hname,"hCaloSumEt_%i",isam);  hCaloSumEtv.push_back(new TH1F(hname,"",30,0,600));     hCaloSumEtv[isam]->Sumw2();
    sprintf(hname,"hCaloMetPhi_%i",isam); hCaloMetPhiv.push_back(new TH1F(hname,"",20,-3.2,3.2)); hCaloMetPhiv[isam]->Sumw2();
    sprintf(hname,"hTCMex_%i",isam);      hTCMexv.push_back(new TH1F(hname,"",40,-50,50));        hTCMexv[isam]->Sumw2();
    sprintf(hname,"hTCMey_%i",isam);      hTCMeyv.push_back(new TH1F(hname,"",40,-50,50));        hTCMeyv[isam]->Sumw2();
    sprintf(hname,"hTCMet_%i",isam);      hTCMetv.push_back(new TH1F(hname,"",20,0,60));          hTCMetv[isam]->Sumw2();
    sprintf(hname,"hTCSumEt_%i",isam);    hTCSumEtv.push_back(new TH1F(hname,"",30,0,600));       hTCSumEtv[isam]->Sumw2();
    sprintf(hname,"hTCMetPhi_%i",isam);   hTCMetPhiv.push_back(new TH1F(hname,"",20,-3.2,3.2));   hTCMetPhiv[isam]->Sumw2();
    sprintf(hname,"hPFMex_%i",isam);      hPFMexv.push_back(new TH1F(hname,"",40,-50,50));        hPFMexv[isam]->Sumw2();
    sprintf(hname,"hPFMey_%i",isam);      hPFMeyv.push_back(new TH1F(hname,"",40,-50,50));        hPFMeyv[isam]->Sumw2();
    sprintf(hname,"hPFMet_%i",isam);      hPFMetv.push_back(new TH1F(hname,"",20,0,60));          hPFMetv[isam]->Sumw2();
    sprintf(hname,"hPFSumEt_%i",isam);    hPFSumEtv.push_back(new TH1F(hname,"",30,0,600));       hPFSumEtv[isam]->Sumw2();
    sprintf(hname,"hPFMetPhi_%i",isam);   hPFMetPhiv.push_back(new TH1F(hname,"",20,-3.2,3.2));   hPFMetPhiv[isam]->Sumw2(); 

//     sprintf(hname,"hNCaloJets_%i",isam);   hNCaloJetsv.push_back(new TH1F(hname,"",10,-0.5,9.5));   hNCaloJetsv[isam]->Sumw2();
//     sprintf(hname,"hCaloJetEt_%i",isam);   hCaloJetEtv.push_back(new TH1F(hname,"",50,0,300));      hCaloJetEtv[isam]->Sumw2();
//     sprintf(hname,"hCaloJetEta_%i",isam);  hCaloJetEtav.push_back(new TH1F(hname,"",40,-5,5));      hCaloJetEtav[isam]->Sumw2();
//     sprintf(hname,"hCaloJetPhi_%i",isam);  hCaloJetPhiv.push_back(new TH1F(hname,"",20,-3.2,3.2));  hCaloJetPhiv[isam]->Sumw2();        
//     sprintf(hname,"hNTrackJets_%i",isam);  hNTrackJetsv.push_back(new TH1F(hname,"",10,-0.5,9.5));  hNTrackJetsv[isam]->Sumw2();
//     sprintf(hname,"hTrackJetEt_%i",isam);  hTrackJetEtv.push_back(new TH1F(hname,"",50,0,300));     hTrackJetEtv[isam]->Sumw2();
//     sprintf(hname,"hTrackJetEta_%i",isam); hTrackJetEtav.push_back(new TH1F(hname,"",40,-5,5));     hTrackJetEtav[isam]->Sumw2();
//     sprintf(hname,"hTrackJetPhi_%i",isam); hTrackJetPhiv.push_back(new TH1F(hname,"",20,-3.2,3.2)); hTrackJetPhiv[isam]->Sumw2();    
    sprintf(hname,"hNPFJets_%i",isam);     hNPFJetsv.push_back(new TH1F(hname,"",10,-0.5,9.5));     hNPFJetsv[isam]->Sumw2();
    sprintf(hname,"hPFJetEt_%i",isam);     hPFJetEtv.push_back(new TH1F(hname,"",50,0,300));        hPFJetEtv[isam]->Sumw2();
    sprintf(hname,"hPFJetEta_%i",isam);    hPFJetEtav.push_back(new TH1F(hname,"",40,-5,5));        hPFJetEtav[isam]->Sumw2();
    sprintf(hname,"hPFJetPhi_%i",isam);    hPFJetPhiv.push_back(new TH1F(hname,"",20,-3.2,3.2));    hPFJetPhiv[isam]->Sumw2();
    
    sprintf(hname,"hNPV_%i",isam);           hNPVv.push_back(new TH1F(hname,"",15,-0.5,14.5));            hNPVv[isam]->Sumw2();
//     sprintf(hname,"hNGoodPV_%i",isam);       hNGoodPVv.push_back(new TH1F(hname,"",15,-0.5,14.5));        hNGoodPVv[isam]->Sumw2();
    sprintf(hname,"hNGoodPV_%s",snamev[isam].Data());       hNGoodPVv.push_back(new TH1F(hname,"",15,-0.5,14.5));        hNGoodPVv[isam]->Sumw2();
    sprintf(hname,"hGoodPVNTracks_%i",isam); hGoodPVNTracksv.push_back(new TH1F(hname,"",50,-0.5,149.5)); hGoodPVNTracksv[isam]->Sumw2();
    sprintf(hname,"hGoodPVChi2_%i",isam);    hGoodPVChi2v.push_back(new TH1F(hname,"",50,0,300));         hGoodPVChi2v[isam]->Sumw2();
    sprintf(hname,"hGoodPVNdof_%i",isam);    hGoodPVNdofv.push_back(new TH1F(hname,"",50,-0.5,299.5));    hGoodPVNdofv[isam]->Sumw2();
    sprintf(hname,"hGoodPVSumPt_%i",isam);   hGoodPVSumPtv.push_back(new TH1F(hname,"",50,0,300));        hGoodPVSumPtv[isam]->Sumw2();
    
    sprintf(hname,"hNTracks0_%i",isam);     hNTracks0v.push_back(new TH1F(hname,"",30,-0.5,299.5));  hNTracks0v[isam]->Sumw2();
    sprintf(hname,"hNCaloTowers0_%i",isam); hNCaloTowers0v.push_back(new TH1F(hname,"",40,100,900)); hNCaloTowers0v[isam]->Sumw2();
      
    sprintf(hname,"hDeltaEtaInB_%i",isam); hDeltaEtaInBv.push_back(new TH1F(hname,"",50,-0.02,0.02)); hDeltaEtaInBv[isam]->Sumw2();
    sprintf(hname,"hDeltaPhiInB_%i",isam); hDeltaPhiInBv.push_back(new TH1F(hname,"",50,-0.1,0.1));   hDeltaPhiInBv[isam]->Sumw2();
    sprintf(hname,"hDeltaEtaInE_%i",isam); hDeltaEtaInEv.push_back(new TH1F(hname,"",50,-0.02,0.02)); hDeltaEtaInEv[isam]->Sumw2();
    sprintf(hname,"hDeltaPhiInE_%i",isam); hDeltaPhiInEv.push_back(new TH1F(hname,"",50,-0.2,0.2));   hDeltaPhiInEv[isam]->Sumw2();
    
    nSelv.push_back(0);
    nSelVarv.push_back(0);
    nPosSSv.push_back(0);
    nNegSSv.push_back(0);    
  }
  
  //
  // Access samples and fill histograms
  //  
  TFile *infile=0;
  TTree *eventTree=0;  
  
  // Data structures to store info from TTrees
  mithep::TEventInfo *info    = new mithep::TEventInfo();
  // TClonesArray *caloJetArr    = new TClonesArray("mithep::TJet");
  // TClonesArray *trackJetArr   = new TClonesArray("mithep::TJet");
  TClonesArray *pfJetArr      = new TClonesArray("mithep::TJet");
  TClonesArray *dielectronArr = new TClonesArray("mithep::TDielectron");
  TClonesArray *pvArr         = new TClonesArray("mithep::TVertex");
  
  //
  // Set up event dump to file
  //
  ofstream evtfile;
  char evtfname[100];    
  sprintf(evtfname,"%s/events.txt",outputDir.Data());
  evtfile.open(evtfname);
  assert(evtfile.is_open());
  
  //
  // loop over samples
  //
  for(UInt_t isam=0; isam<samplev.size(); isam++) {        
    if(isam==0 && !hasData) continue;
    
    //
    // Set up output ntuple file for the sample
    //
    TString outName = ntupDir + TString("/") + snamev[isam] + TString("_select.root");
    TFile *outFile = new TFile(outName,"RECREATE");
    TTree *outTree = new TTree("Events","Events");
    ZeeData data;
    outTree->Branch("Events", &data.runNum, 
"runNum/i:evtNum:lumiSec:nTracks0:nCaloTowers0:nPV:nJets:caloMEx/F:caloMEy:caloSumET:tcMEx:tcMEy:tcSumET:pfMEx:pfMEy:pfSumET:mass:pt:y:phi:pt_1:eta_1:phi_1:scEt_1:scEta_1:scPhi_1:hltMatchBits_1/i:q_1/I:pt_2/F:eta_2:phi_2:scEt_2:scEta_2:scPhi_2:hltMatchBits_2/i:q_2/I:weight/F");
    
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
      if((samp->jsonv.size()>0) && (samp->jsonv[ifile].CompareTo("NONE")!=0)) { 
        hasJSON = kTRUE;
	// rlrm.AddJSONFile(samp->jsonv[ifile].Data()); 
        jsonParser.Initialize(samp->jsonv[ifile].Data()); 
      }
      
      // Get the TTree
      eventTree = (TTree*)infile->Get("Events"); assert(eventTree);

      // Set branch address to structures that will store the info  
      eventTree->SetBranchAddress("Info",       &info);          TBranch *infoBr       = eventTree->GetBranch("Info");
      eventTree->SetBranchAddress("Dielectron", &dielectronArr); TBranch *dielectronBr = eventTree->GetBranch("Dielectron");
      // eventTree->SetBranchAddress("CaloJet",    &caloJetArr);    TBranch *caloJetBr    = eventTree->GetBranch("CaloJet");
      // eventTree->SetBranchAddress("TrackJet",   &trackJetArr);   TBranch *trackJetBr   = eventTree->GetBranch("TrackJet");
      eventTree->SetBranchAddress("PFJet",      &pfJetArr);      TBranch *pfJetBr      = eventTree->GetBranch("PFJet");
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
      Double_t nsel=0, nselvar=0;
      for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {       
	if(ientry >= maxEvents) break;
	
	infoBr->GetEntry(ientry);
		
	// mithep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);      
	// if(hasJSON && !rlrm.HasRunLumi(rl)) continue;  // not certified run? Skip to next event...
        if(hasJSON && !jsonParser.HasRunLumi(info->runNum, info->lumiSec)) continue;  // not certified run? Skip to next event...
        
	// For EPS2011 for both data and MC (starting from Summer11 production)
	// we use an OR of the twi triggers below. Both are unpresecaled.
	UInt_t eventTriggerBit = kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL 
	  | kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL;
	UInt_t leadingTriggerObjectBit = kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele1Obj
	  | kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele1Obj;
	UInt_t trailingTriggerObjectBit = kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele2Obj
	  | kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele2Obj;
	if(isam==0) {
	  nProcessedEvents++;
	  for(UInt_t ibit=0; ibit<32; ibit++) {
            UInt_t itrig = 1<<ibit;
	    if(info->triggerBits & itrig)
	      hTrigger->Fill(ibit);
	  }
 	}
	
        if(!(info->triggerBits & eventTriggerBit)) continue;  // no trigger accept? Skip to next event...                                   
	dielectronArr->Clear(); 
	dielectronBr->GetEntry(ientry);	
	// loop through dielectrons
        for(Int_t i=0; i<dielectronArr->GetEntriesFast(); i++) {
	  mithep::TDielectron *dielectron = (mithep::TDielectron*)((*dielectronArr)[i]);
          // Exclude ECAL gap region and cut out of acceptance electrons
          if((fabs(dielectron->scEta_1)>kGAP_LOW) && (fabs(dielectron->scEta_1)<kGAP_HIGH)) continue;
          if((fabs(dielectron->scEta_2)>kGAP_LOW) && (fabs(dielectron->scEta_2)<kGAP_HIGH)) continue;
          if((fabs(dielectron->scEta_1) > 2.5)       || (fabs(dielectron->scEta_2) > 2.5))       continue;  // outside eta range? Skip to next event...
	  //
	  // Energy scale corrections for data
	  // NOTE: the electrons and dielectron 4-vectors are updated, the supercluster quantities are not
	  //
          Double_t scEt1 = dielectron->scEt_1;
	  Double_t scEt2 = dielectron->scEt_2;
	  // Electron energy scale correction
          if(isam==0) {            	    
    	    double corr1 = escale::findEnergyScaleCorrection(dielectron->scEta_1);
    	    double corr2 = escale::findEnergyScaleCorrection(dielectron->scEta_2);
	    scEt1 = dielectron->scEt_1 * corr1;
	    scEt2 = dielectron->scEt_2 * corr2;

            TLorentzVector ele1; 
	    ele1.SetPtEtaPhiM(dielectron->pt_1,dielectron->eta_1,dielectron->phi_1,0.000511);
	    ele1 *= corr1;
	    dielectron->pt_1  = ele1.Pt();
	    dielectron->eta_1 = ele1.Eta();
	    dielectron->phi_1 = ele1.Phi();
            
	    TLorentzVector ele2; 
	    ele2.SetPtEtaPhiM(dielectron->pt_2,dielectron->eta_2,dielectron->phi_2,0.000511);
	    ele2 *= corr2;
	    dielectron->pt_2  = ele2.Pt();
	    dielectron->eta_2 = ele2.Eta();
	    dielectron->phi_2 = ele2.Phi();
	    
	    TLorentzVector vDiEle = ele1+ele2;            
	    dielectron->mass = vDiEle.M();
	    dielectron->pt   = vDiEle.Pt();
	    dielectron->y    = vDiEle.Rapidity();
	    dielectron->phi  = vDiEle.Phi(); 
          }
       	  
	  // requirements on BOTH electrons
	  // For DY ET cuts are asymmetric:
       	  if( ! ( (scEt1>20 && scEt2>10) || (scEt1>10 && scEt2>20) ) ) continue;
	  // Both electrons must match trigger objects. At least one ordering
	  // must match
	  if( ! ( 
		 (dielectron->hltMatchBits_1 & leadingTriggerObjectBit && 
		  dielectron->hltMatchBits_2 & trailingTriggerObjectBit )
		  ||
		 (dielectron->hltMatchBits_1 & trailingTriggerObjectBit && 
		  dielectron->hltMatchBits_2 & leadingTriggerObjectBit ) ) ) continue;
	  // Other cuts to both electrons

	  // The Smurf electron ID package is the same as used in HWW analysis
	  // and contains cuts like VBTF WP80 for pt>20, VBTF WP70 for pt<10
	  // with some customization, plus impact parameter cuts dz and dxy
	  if(!passSmurf(dielectron)) continue;  
	  
	  hMass2v[isam]->Fill(dielectron->mass,weight);
	  hMass3v[isam]->Fill(dielectron->mass,weight);
	  
          // mass window 
          if((dielectron->mass < 10) || (dielectron->mass > 1000)) continue;
	  
	  
          /******** We have a Z candidate! HURRAY! ********/
	  
          if(dielectron->q_1 == dielectron->q_2) {
	    if(dielectron->q_1 > 0) nPosSSv[isam] += weight;
	    else                    nNegSSv[isam] += weight;
	  }
	  
	  // event printout
          if(isam==0)
            eventDump(evtfile, dielectron, info->runNum, info->lumiSec, info->evtNum, 
		      leadingTriggerObjectBit, trailingTriggerObjectBit);
	  
	  //
	  // Fill histograms
	  // 
	  hMassv[isam]->Fill(dielectron->mass,weight);
          hPtv[isam]  ->Fill(dielectron->pt,  weight);
	  hPt2v[isam] ->Fill(dielectron->pt,  weight);
          hyv[isam]   ->Fill(dielectron->y,   weight);
          hPhiv[isam] ->Fill(dielectron->phi, weight);
      
          hElePtv[isam] ->Fill(dielectron->pt_1, weight); hElePtv[isam] ->Fill(dielectron->pt_2, weight);
          hEleEtav[isam]->Fill(dielectron->eta_1,weight); hEleEtav[isam]->Fill(dielectron->eta_2,weight);
          hElePhiv[isam]->Fill(dielectron->phi_1,weight); hElePhiv[isam]->Fill(dielectron->phi_2,weight);
	
          const Bool_t isB1 = fabs(dielectron->scEta_1)<kGAP_LOW;
	  const Bool_t isB2 = fabs(dielectron->scEta_2)<kGAP_LOW;
	  
	  if(isB1) {
	    hDeltaEtaInBv[isam]->Fill(dielectron->deltaEtaIn_1,weight);
	    hDeltaPhiInBv[isam]->Fill(dielectron->deltaPhiIn_1,weight);
	  } else {
	    hDeltaEtaInEv[isam]->Fill(dielectron->deltaEtaIn_1,weight);
	    hDeltaPhiInEv[isam]->Fill(dielectron->deltaPhiIn_1,weight);
	  }
	  if(isB2) {
	    hDeltaEtaInBv[isam]->Fill(dielectron->deltaEtaIn_2,weight);
	    hDeltaPhiInBv[isam]->Fill(dielectron->deltaPhiIn_2,weight);
	  } else {
	    hDeltaEtaInEv[isam]->Fill(dielectron->deltaEtaIn_2,weight);
	    hDeltaPhiInEv[isam]->Fill(dielectron->deltaPhiIn_2,weight);
	  }
	  
	  if(isB1 && isB2)        { hMassBBv[isam]->Fill(dielectron->mass,weight); } 
	  else if(!isB1 && !isB2) { hMassBEv[isam]->Fill(dielectron->mass,weight); } 
	  else                    { hMassEEv[isam]->Fill(dielectron->mass,weight);}
	  
 	  const Double_t jetEtMin = 15;
	  
// 	  UInt_t ncalojets=0;
// 	  caloJetArr->Clear();
// 	  caloJetBr->GetEntry(ientry);	  
// 	  for(Int_t ijet=0; ijet<caloJetArr->GetEntriesFast(); ijet++) {
// 	    const mithep::TJet *jet = (mithep::TJet*)((*caloJetArr)[ijet]);
// 	    if(toolbox::deltaR(dielectron->eta_1,dielectron->phi_1,jet->eta,jet->phi)<0.5) continue;
// 	    if(toolbox::deltaR(dielectron->eta_2,dielectron->phi_2,jet->eta,jet->phi)<0.5) continue;
// 	    hCaloJetEtv[isam]->Fill(jet->et,weight);
// 	    hCaloJetEtav[isam]->Fill(jet->eta,weight);
// 	    hCaloJetPhiv[isam]->Fill(jet->phi,weight);
// 	    if(jet->et > jetEtMin) ncalojets++;
// 	  }
// 	  hNCaloJetsv[isam]->Fill(ncalojets,weight);
	  
// 	  UInt_t ntkjets=0;
//           trackJetArr->Clear();
// 	  trackJetBr->GetEntry(ientry);
// 	  for(Int_t ijet=0; ijet<trackJetArr->GetEntriesFast(); ijet++) {
// 	    const mithep::TJet *jet = (mithep::TJet*)((*trackJetArr)[ijet]);
// 	    if(toolbox::deltaR(dielectron->eta_1,dielectron->phi_1,jet->eta,jet->phi)<0.5) continue;
// 	    if(toolbox::deltaR(dielectron->eta_2,dielectron->phi_2,jet->eta,jet->phi)<0.5) continue;
// 	    hTrackJetEtv[isam]->Fill(jet->et,weight);
// 	    hTrackJetEtav[isam]->Fill(jet->eta,weight);
// 	    hTrackJetPhiv[isam]->Fill(jet->phi,weight);
// 	    if(jet->et > jetEtMin) ntkjets++;
// 	  }
// 	  hNTrackJetsv[isam]->Fill(ntkjets,weight);
          	  
 	  UInt_t npfjets=0;	  
 	  pfJetArr->Clear();
 	  pfJetBr->GetEntry(ientry);	  
 	  for(Int_t ijet=0; ijet<pfJetArr->GetEntriesFast(); ijet++) {
 	    const mithep::TJet *jet = (mithep::TJet*)((*pfJetArr)[ijet]);
 	    if(toolbox::deltaR(dielectron->eta_1,dielectron->phi_1,jet->eta,jet->phi)<0.5) continue;
 	    if(toolbox::deltaR(dielectron->eta_2,dielectron->phi_2,jet->eta,jet->phi)<0.5) continue;
 	    hPFJetEtv[isam]->Fill(jet->pt,weight);
 	    hPFJetEtav[isam]->Fill(jet->eta,weight);
 	    hPFJetPhiv[isam]->Fill(jet->phi,weight);
 	    if(jet->pt > jetEtMin) npfjets++;
 	  }
 	  hNPFJetsv[isam]->Fill(npfjets,weight);
	  
	  pvArr->Clear();
          pvBr->GetEntry(ientry);
          hNPVv[isam]->Fill(pvArr->GetEntriesFast(),weight);
          UInt_t nGoodPV=0;
          for(Int_t ipv=0; ipv<pvArr->GetEntriesFast(); ipv++) {
            const mithep::TVertex *pv = (mithep::TVertex*)((*pvArr)[ipv]);                  
            if(pv->nTracksFit                        < 1)  continue;
            if(pv->ndof                              < 4)  continue;
            if(fabs(pv->z)                           > 24) continue;
            if(sqrt((pv->x)*(pv->x)+(pv->y)*(pv->y)) > 2)  continue;
            nGoodPV++;
            hGoodPVNTracksv[isam]->Fill(pv->nTracksFit,weight);
            hGoodPVChi2v[isam]   ->Fill(pv->chi2,weight);
            hGoodPVNdofv[isam]   ->Fill(pv->ndof,weight);
            hGoodPVSumPtv[isam]  ->Fill(pv->sumPt,weight);
            if(isam!=0) continue;
            for(Int_t jpv=ipv+1; jpv<pvArr->GetEntriesFast(); jpv++) {
              const mithep::TVertex *pv2 = (mithep::TVertex*)((*pvArr)[jpv]);
              if(pv2->nTracksFit                           < 1)  continue;
              if(pv2->ndof                                 < 4)  continue;
              if(fabs(pv2->z)                              > 24) continue;
              if(sqrt((pv2->x)*(pv2->x)+(pv2->y)*(pv2->y)) > 2)  continue;
              hGoodPVDz->Fill(pv->z - pv2->z,weight);
            }
          }
          hNGoodPVv[isam]->Fill(nGoodPV,weight);
	  
	  // fill ntuple data
	  fillData(&data, info, dielectron, pvArr->GetEntriesFast(), npfjets, weight);
	  outTree->Fill();
	  
	  nsel    += weight;
	  nselvar += weight*weight;

      
//           TVector3 met;
      
//           if(info->caloMEx!=0 || info->caloMEy!=0) {       
//             hCaloMexv[isam]->Fill(info->caloMEx,weight);
//             hCaloMeyv[isam]->Fill(info->caloMEy,weight);
//             met.SetXYZ(info->caloMEx, info->caloMEy, 0);
//             hCaloMetv[isam]   ->Fill(met.Perp(),weight);
//             hCaloMetPhiv[isam]->Fill(met.Phi(),weight);
//             hCaloSumEtv[isam] ->Fill(info->caloSumET,weight);
//           }
      
//           if(info->tcMEx!=0 || info->tcMEy!=0) {       
//             hTCMexv[isam]->Fill(info->tcMEx,weight);
//             hTCMeyv[isam]->Fill(info->tcMEy,weight);
//             met.SetXYZ(info->tcMEx, info->tcMEy, 0);
//             hTCMetv[isam]   ->Fill(met.Perp(),weight);
//             hTCMetPhiv[isam]->Fill(met.Phi(),weight);
//             hTCSumEtv[isam] ->Fill(info->tcSumET,weight);
//           }
      
//           if(info->pfMEx!=0 || info->pfMEy!=0) {       
//             hPFMexv[isam]->Fill(info->pfMEx,weight);
//             hPFMeyv[isam]->Fill(info->pfMEy,weight);
//             met.SetXYZ(info->pfMEx, info->pfMEy, 0);
//             hPFMetv[isam]   ->Fill(met.Perp(),weight);
//             hPFMetPhiv[isam]->Fill(met.Phi(),weight);
//             hPFSumEtv[isam] ->Fill(info->pfSumET,weight);
//           }
	  
// 	  hNTracks0v[isam]    ->Fill(info->nTracks0,weight);
//           hNCaloTowers0v[isam]->Fill(info->nCaloTowers0,weight);
        }	 
      }           
      cout << nsel << " +/- " << sqrt(nselvar) << " events" << endl;
      nSelv[isam]    += nsel;
      nSelVarv[isam] += nselvar;
      delete infile;
      infile=0, eventTree=0;
    }
    outFile->Write();
    delete outTree;
    outFile->Close();        
    delete outFile;
  }
  delete info;
  delete dielectronArr;
  delete pvArr;
  evtfile.close();


  // Write useful histograms
  TString outNamePV = outputDir + TString("/npv.root");
  TFile *outFilePV = new TFile(outNamePV,"RECREATE");
  outFilePV->cd();
  for(UInt_t isam=0; isam<samplev.size(); isam++) {
    hNGoodPVv[isam]->Write();
  }
  outFilePV->Close();

  //--------------------------------------------------------------------------------------------------------------
  // Make plots
  //==============================================================================================================
  TCanvas *c = MakeCanvas("c","c",canw,canh);

  printf("Make plots\n");fflush(stdout);
  // string buffers
  char ylabel[100];   // y-axis label
  char lumitext[50];
  if(lumi>0) {
    if(lumi<1) { sprintf(lumitext,"#int#font[12]{L}dt = %.0f nb^{-1}",1000.*lumi); }
    else       { sprintf(lumitext,"#int#font[12]{L}dt = %.3g pb^{-1}",lumi); }
  }
      
  if(hasData) {
    CPlot plotTrigger("trigger","","trigger bit","Fraction of Processed Events");
    hTrigger->Scale(1.0/(Double_t)nProcessedEvents);
    plotTrigger.AddHist1D(hTrigger);
    plotTrigger.Draw(c,kFALSE,format);
  }
  
  // scale factor for yield in MC to equal yield in data
  Double_t mcscale=1;
  if(hasData) {
    Double_t numer = nSelv[0];
    Double_t denom = 0;
    for(UInt_t isam=1; isam<samplev.size(); isam++)
      denom += nSelv[isam];
    mcscale = (denom>0) ? numer/denom : 1.0;
  }
  
  printf("Plot dielectron mass\n");fflush(stdout);
  // dielectron mass
  sprintf(ylabel,"Events / %.1f GeV/c^{2}",hMassv[0]->GetBinWidth(1));
  CPlot plotMass("mass","","m(e^{+}e^{-}) [GeV/c^{2}]",ylabel);
  if(hasData) { plotMass.AddHist1D(hMassv[0],samplev[0]->label,"E"); }
  printf("  debug1\n");fflush(stdout);
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    hMassv[isam]->Scale(mcscale);
    plotMass.AddToStack(hMassv[isam],samplev[isam]->label,samplev[isam]->color);
  }
  printf("  debug2\n");fflush(stdout);
  if(samplev.size()>5)
    plotMass.SetLegend(0.75,0.55,0.98,0.9);
  else
    plotMass.TransLegend(0.1,0);
  if(lumi>0) plotMass.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
  plotMass.Draw(c,kFALSE,format);
  
  printf("  debug3\n");fflush(stdout);
  plotMass.SetName("masslog");
  printf("  debug4\n");fflush(stdout);
  plotMass.SetLogy();
  printf("  debug5\n");fflush(stdout);
  if( plotMass.GetStack() != NULL)
    plotMass.SetYRange((1e-4)*(plotMass.GetStack()->GetMaximum()),10.*(plotMass.GetStack()->GetMaximum()));  
  printf("  debug6\n");fflush(stdout);
  plotMass.Draw(c,kFALSE,format);

  printf("Make many other plots\n");fflush(stdout);
  bool zillionPlots = false;
  if(zillionPlots){

    sprintf(ylabel,"Events / %.1f GeV/c^{2}",hMass2v[0]->GetBinWidth(1));
    CPlot plotMass2("mass2","","m(e^{+}e^{-}) [GeV/c^{2}]",ylabel);
    if(hasData) { plotMass2.AddHist1D(hMass2v[0],samplev[0]->label,"E"); }
    for(UInt_t isam=1; isam<samplev.size(); isam++) {
      hMass2v[isam]->Scale(mcscale);
      plotMass2.AddToStack(hMass2v[isam],samplev[isam]->label,samplev[isam]->color);
    }
    if(samplev.size()>5)
      plotMass2.SetLegend(0.75,0.55,0.98,0.9);
    else
      plotMass2.TransLegend(0.1,0);
    if(lumi>0) plotMass2.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
    plotMass2.SetYRange((1e-4)*(plotMass2.GetStack()->GetMaximum()),10.*(plotMass2.GetStack()->GetMaximum()));
    plotMass2.SetLogy();
    plotMass2.Draw(c,kFALSE,format);
    
  sprintf(ylabel,"Events / %.1f GeV/c^{2}",hMass3v[0]->GetBinWidth(1));
  CPlot plotMass3("mass3","","m(e^{+}e^{-}) [GeV/c^{2}]",ylabel);
  if(hasData) { plotMass3.AddHist1D(hMass3v[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    hMass3v[isam]->Scale(mcscale);
    plotMass3.AddToStack(hMass3v[isam],samplev[isam]->label,samplev[isam]->color);
  }
  if(samplev.size()>5)
    plotMass3.SetLegend(0.75,0.55,0.98,0.9);
  else
    plotMass3.TransLegend(0.1,0);
  if(lumi>0) plotMass3.AddTextBox(lumitext,0.45,0.85,0.65,0.8,0);
  plotMass3.SetYRange((1e-4)*(plotMass3.GetStack()->GetMaximum()),10.*(plotMass3.GetStack()->GetMaximum()));
  plotMass3.SetLogy();
  plotMass3.Draw(c,kFALSE,format);

  sprintf(ylabel,"Events / %.1f GeV/c^{2}",hMassBBv[0]->GetBinWidth(1));
  CPlot plotMassBB("massBB","barrel-barrel","m(e^{+}e^{-}) [GeV/c^{2}]",ylabel);
  if(hasData) { plotMassBB.AddHist1D(hMassBBv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    hMassBBv[isam]->Scale(mcscale);
    plotMassBB.AddToStack(hMassBBv[isam],samplev[isam]->label,samplev[isam]->color);
  }
  if(samplev.size()>5)
    plotMassBB.SetLegend(0.75,0.55,0.98,0.9);
  else
    plotMassBB.TransLegend(0.1,0);
  if(lumi>0) plotMassBB.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
  plotMassBB.Draw(c,kFALSE,format);  
  
  sprintf(ylabel,"Events / %.1f GeV/c^{2}",hMassBEv[0]->GetBinWidth(1));
  CPlot plotMassBE("massBE","barrel-endcap","m(e^{+}e^{-}) [GeV/c^{2}]",ylabel);
  if(hasData) { plotMassBE.AddHist1D(hMassBEv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    hMassBEv[isam]->Scale(mcscale);
    plotMassBE.AddToStack(hMassBEv[isam],samplev[isam]->label,samplev[isam]->color);
  }
  if(samplev.size()>5)
    plotMassBE.SetLegend(0.75,0.55,0.98,0.9);
  else
    plotMassBE.TransLegend(0.1,0);
  if(lumi>0) plotMassBE.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
  plotMassBE.Draw(c,kFALSE,format);
  
  sprintf(ylabel,"Events / %.1f GeV/c^{2}",hMassEEv[0]->GetBinWidth(1));
  CPlot plotMassEE("massEE","endcap-endcap","m(e^{+}e^{-}) [GeV/c^{2}]",ylabel);
  if(hasData) { plotMassEE.AddHist1D(hMassEEv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    hMassEEv[isam]->Scale(mcscale);
    plotMassEE.AddToStack(hMassEEv[isam],samplev[isam]->label,samplev[isam]->color);
  }
  if(samplev.size()>5)
    plotMassEE.SetLegend(0.75,0.55,0.98,0.9);
  else
    plotMassEE.TransLegend(0.1,0);
  if(lumi>0) plotMassEE.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
  plotMassEE.Draw(c,kFALSE,format);
    
  // dielectron pT
  sprintf(ylabel,"Events / %.1f GeV/c",hPtv[0]->GetBinWidth(1));
  CPlot plotPt("pt","","p_{T}(e^{+}e^{-}) [GeV/c]",ylabel);
  if(hasData) { plotPt.AddHist1D(hPtv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    hPtv[isam]->Scale(mcscale);
    plotPt.AddToStack(hPtv[isam],samplev[isam]->label,samplev[isam]->color);
  }
  if(lumi>0) plotPt.AddTextBox(lumitext,0.41,0.85,0.61,0.8,0);
  if(samplev.size()>5)
    plotPt.SetLegend(0.75,0.55,0.98,0.9);
  else
    plotPt.TransLegend(0.1,0);
  plotPt.Draw(c,kFALSE,format);
  
  sprintf(ylabel,"Events / %.1f GeV/c",hPt2v[0]->GetBinWidth(1));
  CPlot plotPt2("ptlog","","p_{T}(#mu^{+}#mu^{-}) [GeV/c]",ylabel);
  if(hasData) { plotPt2.AddHist1D(hPt2v[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    hPt2v[isam]->Scale(mcscale);
    plotPt2.AddToStack(hPt2v[isam],samplev[isam]->label,samplev[isam]->color);
  }
  if(lumi>0) plotPt2.AddTextBox(lumitext,0.41,0.85,0.61,0.8,0);
  if(samplev.size()>5)
    plotPt2.SetLegend(0.75,0.55,0.98,0.9);
  else
    plotPt2.TransLegend(0.1,0);
  plotPt2.SetLogy();
  plotPt2.SetYRange((1e-4)*(plotPt2.GetStack()->GetMaximum()),10.*(plotPt2.GetStack()->GetMaximum())); 
  plotPt2.Draw(c,kFALSE,format);
    
  // dielectron eta
  sprintf(ylabel,"Events / %.2f",hyv[0]->GetBinWidth(1));
  CPlot ploty("y","","y(e^{+}e^{-})",ylabel);
  if(hasData) { ploty.AddHist1D(hyv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    hyv[isam]->Scale(mcscale);
    ploty.AddToStack(hyv[isam],samplev[isam]->label,samplev[isam]->color);
  }
  if(lumi>0) ploty.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
  if(samplev.size()>5)
    ploty.SetLegend(0.75,0.55,0.98,0.9);
  else
    ploty.TransLegend(0.1,0);
  ploty.SetYRange(0,2.0*(ploty.GetStack()->GetMaximum()));
  ploty.Draw(c,kFALSE,format);
    
  // dielectron phi
  sprintf(ylabel,"Events / %.2f",hPhiv[0]->GetBinWidth(1));
  CPlot plotPhi("phi","","#phi(e^{+}e^{-})",ylabel);
  if(hasData) { plotPhi.AddHist1D(hPhiv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    hPhiv[isam]->Scale(mcscale);
    plotPhi.AddToStack(hPhiv[isam],samplev[isam]->label,samplev[isam]->color);
  }
  if(lumi>0) plotPhi.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
  plotPhi.SetYRange(0,1.8*(plotPhi.GetStack()->GetMaximum()));
  if(samplev.size()>5)
    plotPhi.SetLegend(0.75,0.55,0.98,0.9);
  else
    plotPhi.TransLegend(0.1,0);
  plotPhi.Draw(c,kFALSE,format);
    
  // electron pT
  sprintf(ylabel,"Events / %.1f GeV/c",hElePtv[0]->GetBinWidth(1));
  CPlot plotElePt("elept","","p_{T}(e) [GeV/c]",ylabel);
  if(hasData) { plotElePt.AddHist1D(hElePtv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    hElePtv[isam]->Scale(mcscale);
    plotElePt.AddToStack(hElePtv[isam],samplev[isam]->label,samplev[isam]->color);
  }
  if(lumi>0) plotElePt.AddTextBox(lumitext,0.45,0.85,0.65,0.8,0);
  if(samplev.size()>5)
    plotElePt.SetLegend(0.75,0.55,0.98,0.9);
  else
    plotElePt.TransLegend(0.1,0);
  plotElePt.Draw(c,kFALSE,format);
  
  plotElePt.SetName("eleptlog");  
  plotElePt.SetLogy();
  plotElePt.SetYRange((1e-4)*(plotElePt.GetStack()->GetMaximum()),10.*(plotElePt.GetStack()->GetMaximum())); 
  plotElePt.Draw(c,kFALSE,format);
    
  // electron eta
  sprintf(ylabel,"Events / %.2f",hEleEtav[0]->GetBinWidth(1));
  CPlot plotEleEta("eleeta","","#eta(e) [GeV/c]",ylabel);
  if(hasData) { plotEleEta.AddHist1D(hEleEtav[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    hEleEtav[isam]->Scale(mcscale);
    plotEleEta.AddToStack(hEleEtav[isam],samplev[isam]->label,samplev[isam]->color);
  }
  if(lumi>0) plotEleEta.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
  if(samplev.size()>5)
    plotEleEta.SetLegend(0.75,0.55,0.98,0.9);
  else
    plotEleEta.TransLegend(0.1,0);
  plotEleEta.SetYRange(0,1.5*(plotEleEta.GetStack()->GetMaximum()));
  plotEleEta.Draw(c,kFALSE,format);
    
  // electron phi
  sprintf(ylabel,"Events / %.2f",hElePhiv[0]->GetBinWidth(1));
  CPlot plotElePhi("elephi","","#phi(e) [GeV/c]",ylabel);
  if(hasData) { plotElePhi.AddHist1D(hElePhiv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    hElePhiv[isam]->Scale(mcscale);
    plotElePhi.AddToStack(hElePhiv[isam],samplev[isam]->label,samplev[isam]->color);
  }
  if(lumi>0) plotElePhi.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
  if(samplev.size()>5)
    plotElePhi.SetLegend(0.75,0.55,0.98,0.9);
  else
    plotElePhi.TransLegend(0.1,0);
  plotElePhi.SetYRange(0,2.0*(plotElePhi.GetStack()->GetMaximum()));
  plotElePhi.Draw(c,kFALSE,format);

  //
  // MET and SumET
  // 
  sprintf(ylabel,"Events / %.1f GeV",hCaloMexv[0]->GetBinWidth(1));
  CPlot plotCaloMex("calomex","","Calo #slash{E}_{x} [GeV]",ylabel);
  if(hasData) { plotCaloMex.AddHist1D(hCaloMexv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    hCaloMexv[isam]->Scale(mcscale);
    plotCaloMex.AddToStack(hCaloMexv[isam],samplev[isam]->label,samplev[isam]->color);
  }
  if(lumi>0) plotCaloMex.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
  plotCaloMex.TransLegend(0.1,0);
  plotCaloMex.Draw(c,kFALSE,format);
  
  sprintf(ylabel,"Events / %.1f GeV",hCaloMeyv[0]->GetBinWidth(1));
  CPlot plotCaloMey("calomey","","Calo #slash{E}_{y} [GeV]",ylabel);
  if(hasData) { plotCaloMey.AddHist1D(hCaloMeyv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    hCaloMeyv[isam]->Scale(mcscale);
    plotCaloMey.AddToStack(hCaloMeyv[isam],samplev[isam]->label,samplev[isam]->color);
  }
  if(lumi>0) plotCaloMey.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
  plotCaloMey.TransLegend(0.1,0);
  plotCaloMey.Draw(c,kFALSE,format);
    
  sprintf(ylabel,"Events / %.1f GeV",hCaloMetv[0]->GetBinWidth(1));
  CPlot plotCaloMet("calomet","","Calo #slash{E}_{T} [GeV]",ylabel);
  if(hasData) { plotCaloMet.AddHist1D(hCaloMetv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    hCaloMetv[isam]->Scale(mcscale);
    plotCaloMet.AddToStack(hCaloMetv[isam],samplev[isam]->label,samplev[isam]->color);
  }
  if(lumi>0) plotCaloMet.AddTextBox(lumitext,0.35,0.85,0.55,0.8,0);
  plotCaloMet.TransLegend(0.05,0);
  plotCaloMet.Draw(c,kFALSE,format);

  sprintf(ylabel,"Events / %.1f GeV",hCaloMetPhiv[0]->GetBinWidth(1));
  CPlot plotCaloMetPhi("calometphi","","Calo #phi(#slash{E}_{T})",ylabel);
  if(hasData) { plotCaloMetPhi.AddHist1D(hCaloMetPhiv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    hCaloMetPhiv[isam]->Scale(mcscale);
    plotCaloMetPhi.AddToStack(hCaloMetPhiv[isam],samplev[isam]->label,samplev[isam]->color);
  }
  if(lumi>0) plotCaloMetPhi.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
  plotCaloMetPhi.SetYRange(0,1.8*(plotCaloMetPhi.GetStack()->GetMaximum()));
  plotCaloMetPhi.TransLegend(0.1,0);
  plotCaloMetPhi.Draw(c,kFALSE,format);
  
  sprintf(ylabel,"Events / %.1f GeV/c",hCaloSumEtv[0]->GetBinWidth(1));
  CPlot plotCaloSumEt("calosumet","","Calo #Sigma_{}E_{T} [GeV]",ylabel);
  if(hasData) { plotCaloSumEt.AddHist1D(hCaloSumEtv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    hCaloSumEtv[isam]->Scale(mcscale);
    plotCaloSumEt.AddToStack(hCaloSumEtv[isam],samplev[isam]->label,samplev[isam]->color);
  }
  if(lumi>0) plotCaloSumEt.AddTextBox(lumitext,0.7,0.4,0.9,0.45,0);
  if(samplev.size()>5)
    plotCaloSumEt.SetLegend(0.7,0.55,0.9,0.9);
  else
    plotCaloSumEt.TransLegend(0.1,0);
  plotCaloSumEt.Draw(c,kFALSE,format);

  sprintf(ylabel,"Events / %.1f GeV",hTCMexv[0]->GetBinWidth(1));
  CPlot plotTCMex("tcmex","","TC #slash{E}_{x} [GeV]",ylabel);
  if(hasData) { plotTCMex.AddHist1D(hTCMexv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    hTCMexv[isam]->Scale(mcscale);
    plotTCMex.AddToStack(hTCMexv[isam],samplev[isam]->label,samplev[isam]->color);
  }
  if(lumi>0) plotTCMex.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
  plotTCMex.TransLegend(0.1,0);
  plotTCMex.Draw(c,kFALSE,format);
  
  sprintf(ylabel,"Events / %.1f GeV",hTCMeyv[0]->GetBinWidth(1));
  CPlot plotTCMey("tcmey","","TC #slash{E}_{y} [GeV]",ylabel);
  if(hasData) { plotTCMey.AddHist1D(hTCMeyv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    hTCMeyv[isam]->Scale(mcscale);
    plotTCMey.AddToStack(hTCMeyv[isam],samplev[isam]->label,samplev[isam]->color);
  }
  if(lumi>0) plotTCMey.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
  plotTCMey.TransLegend(0.1,0);
  plotTCMey.Draw(c,kFALSE,format);

  sprintf(ylabel,"Events / %.1f GeV",hTCMetv[0]->GetBinWidth(1));
  CPlot plotTCMet("tcmet","","TC #slash{E}_{T} [GeV]",ylabel);
  if(hasData) { plotTCMet.AddHist1D(hTCMetv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    hTCMetv[isam]->Scale(mcscale);
    plotTCMet.AddToStack(hTCMetv[isam],samplev[isam]->label,samplev[isam]->color);
  }
  if(lumi>0) plotTCMet.AddTextBox(lumitext,0.35,0.85,0.55,0.8,0);
  plotTCMet.TransLegend(0.05,0);
  plotTCMet.Draw(c,kFALSE,format);

  sprintf(ylabel,"Events / %.1f GeV",hTCMetPhiv[0]->GetBinWidth(1));
  CPlot plotTCMetPhi("tcmetphi","","TC #phi(#slash{E}_{T})",ylabel);
  if(hasData) { plotTCMetPhi.AddHist1D(hTCMetPhiv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    hTCMetPhiv[isam]->Scale(mcscale);
    plotTCMetPhi.AddToStack(hTCMetPhiv[isam],samplev[isam]->label,samplev[isam]->color);
  }
  if(lumi>0) plotTCMetPhi.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
  plotTCMetPhi.SetYRange(0,1.8*(plotTCMetPhi.GetStack()->GetMaximum()));
  plotTCMetPhi.TransLegend(0.1,0);
  plotTCMetPhi.Draw(c,kFALSE,format);
  
  sprintf(ylabel,"Events / %.1f GeV/c",hTCSumEtv[0]->GetBinWidth(1));
  CPlot plotTCSumEt("tcsumet","","TC #Sigma_{}E_{T} [GeV]",ylabel);
  if(hasData) { plotTCSumEt.AddHist1D(hTCSumEtv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    hTCSumEtv[isam]->Scale(mcscale);
    plotTCSumEt.AddToStack(hTCSumEtv[isam],samplev[isam]->label,samplev[isam]->color);
  }
  if(lumi>0) plotTCSumEt.AddTextBox(lumitext,0.7,0.4,0.9,0.45,0);
  if(samplev.size()>5)
    plotTCSumEt.SetLegend(0.7,0.55,0.9,0.9);
  else
    plotTCSumEt.TransLegend(0.1,0);
  plotTCSumEt.Draw(c,kFALSE,format);
  
  sprintf(ylabel,"Events / %.1f GeV",hPFMexv[0]->GetBinWidth(1));
  CPlot plotPFMex("pfmex","","PF #slash{E}_{x} [GeV]",ylabel);
  if(hasData) { plotPFMex.AddHist1D(hPFMexv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    hPFMexv[isam]->Scale(mcscale);
    plotPFMex.AddToStack(hPFMexv[isam],samplev[isam]->label,samplev[isam]->color);
  }
  if(lumi>0) plotPFMex.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
  plotPFMex.TransLegend(0.1,0);
  plotPFMex.Draw(c,kFALSE,format);
  
  sprintf(ylabel,"Events / %.1f GeV",hPFMeyv[0]->GetBinWidth(1));
  CPlot plotPFMey("pfmey","","PF #slash{E}_{y} [GeV]",ylabel);
  if(hasData) { plotPFMey.AddHist1D(hPFMeyv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    hPFMeyv[isam]->Scale(mcscale);
    plotPFMey.AddToStack(hPFMeyv[isam],samplev[isam]->label,samplev[isam]->color);
  }
  if(lumi>0) plotPFMey.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
  plotPFMey.TransLegend(0.1,0);
  plotPFMey.Draw(c,kFALSE,format);
    
  sprintf(ylabel,"Events / %.1f GeV",hPFMetv[0]->GetBinWidth(1));
  CPlot plotPFMet("pfmet","","PF #slash{E}_{T} [GeV]",ylabel);
  if(hasData) { plotPFMet.AddHist1D(hPFMetv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    hPFMetv[isam]->Scale(mcscale);
    plotPFMet.AddToStack(hPFMetv[isam],samplev[isam]->label,samplev[isam]->color);
  }
  if(lumi>0) plotPFMet.AddTextBox(lumitext,0.35,0.85,0.55,0.8,0);
  plotPFMet.TransLegend(0.05,0);
  plotPFMet.Draw(c,kFALSE,format);
  
  sprintf(ylabel,"Events / %.1f GeV",hPFMetPhiv[0]->GetBinWidth(1));
  CPlot plotPFMetPhi("pfmetphi","","PF #phi(#slash{E}_{T})",ylabel);
  if(hasData) { plotPFMetPhi.AddHist1D(hPFMetPhiv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    hPFMetPhiv[isam]->Scale(mcscale);
    plotPFMetPhi.AddToStack(hPFMetPhiv[isam],samplev[isam]->label,samplev[isam]->color);
  }
  if(lumi>0) plotPFMetPhi.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
  plotPFMetPhi.SetYRange(0,1.8*(plotPFMetPhi.GetStack()->GetMaximum()));
  plotPFMetPhi.TransLegend(0.1,0);
  plotPFMetPhi.Draw(c,kFALSE,format);
  
  sprintf(ylabel,"Events / %.1f GeV/c",hPFSumEtv[0]->GetBinWidth(1));
  CPlot plotPFSumEt("pfsumet","","PF #Sigma_{}E_{T} [GeV]",ylabel);
  if(hasData) { plotPFSumEt.AddHist1D(hPFSumEtv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    hPFSumEtv[isam]->Scale(mcscale);
    plotPFSumEt.AddToStack(hPFSumEtv[isam],samplev[isam]->label,samplev[isam]->color);
  }
  if(lumi>0) plotPFSumEt.AddTextBox(lumitext,0.7,0.4,0.9,0.45,0);
  if(samplev.size()>5) 
    plotPFSumEt.SetLegend(0.7,0.55,0.9,0.9);
  else 
    plotPFSumEt.TransLegend(0.1,0);
  plotPFSumEt.Draw(c,kFALSE,format);


  //
  // Jet Multiplicity
  //
//   CPlot plotNCaloJets("ncalojets","Jet E_{T} > 15 GeV","N_{calo jets}","Events");
//   if(hasData) { plotNCaloJets.AddHist1D(hNCaloJetsv[0],samplev[0]->label,"E"); }
//   for(UInt_t isam=1; isam<samplev.size(); isam++) {
//     hNCaloJetsv[isam]->Scale(mcscale);
//     plotNCaloJets.AddToStack(hNCaloJetsv[isam],samplev[isam]->label,samplev[isam]->color);
//   }
//   if(lumi>0) plotNCaloJets.AddTextBox(lumitext,0.43,0.85,0.63,0.8,0);
//   plotNCaloJets.TransLegend(0.08,0);
//   plotNCaloJets.SetLogy();
//   plotNCaloJets.Draw(c,kFALSE,format);

//   CPlot plotNTrackJets("ntrackjets","Jet E_{T} > 15 GeV","N_{track jets}","Events");
//   if(hasData) { plotNTrackJets.AddHist1D(hNTrackJetsv[0],samplev[0]->label,"E"); }
//   for(UInt_t isam=1; isam<samplev.size(); isam++) {
//     hNTrackJetsv[isam]->Scale(mcscale);
//     plotNTrackJets.AddToStack(hNTrackJetsv[isam],samplev[isam]->label,samplev[isam]->color);
//   }
//   if(lumi>0) plotNTrackJets.AddTextBox(lumitext,0.43,0.85,0.63,0.8,0);
//   plotNTrackJets.TransLegend(0.08,0);
//   plotNTrackJets.SetLogy();
//   plotNTrackJets.Draw(c,kFALSE,format);
  
  CPlot plotNPFJets("npfjets","Jet E_{T} > 15 GeV","N_{PF jets}","Events");
  if(hasData) { plotNPFJets.AddHist1D(hNPFJetsv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    hNPFJetsv[isam]->Scale(mcscale);
    plotNPFJets.AddToStack(hNPFJetsv[isam],samplev[isam]->label,samplev[isam]->color);
  }
  if(lumi>0) plotNPFJets.AddTextBox(lumitext,0.43,0.85,0.63,0.8,0);
  plotNPFJets.TransLegend(0.08,0);
  plotNPFJets.SetLogy();
  plotNPFJets.Draw(c,kFALSE,format);

  //
  // Calo Jet kinematics
  //
//   sprintf(ylabel,"Events / %.1f GeV",hCaloJetEtv[0]->GetBinWidth(1));
//   CPlot plotCaloJetEt("calojetet","Jet E_{T} > 15 GeV","Calo Jet E_{T} [GeV]",ylabel);
//   if(hasData) { plotCaloJetEt.AddHist1D(hCaloJetEtv[0],samplev[0]->label,"E"); }
//   for(UInt_t isam=1; isam<samplev.size(); isam++) {
//     hCaloJetEtv[isam]->Scale(mcscale);
//     plotCaloJetEt.AddToStack(hCaloJetEtv[isam],samplev[isam]->label,samplev[isam]->color);
//   }
//   if(lumi>0) plotCaloJetEt.AddTextBox(lumitext,0.41,0.85,0.61,0.8,0);
//   if(samplev.size()>5) 
//     plotCaloJetEt.SetLegend(0.7,0.55,0.9,0.9);
//   else 
//     plotCaloJetEt.TransLegend(0.1,0);  
//   plotCaloJetEt.Draw(c,kFALSE,format);
  
//   plotCaloJetEt.SetName("calojetetlog");  
//   plotCaloJetEt.SetLogy();  
//   plotCaloJetEt.SetYRange((1e-4)*(plotCaloJetEt.GetStack()->GetMaximum()),10.*(plotCaloJetEt.GetStack()->GetMaximum()));
//   plotCaloJetEt.Draw(c,kFALSE,format);
    
//   sprintf(ylabel,"Events / %.2f",hCaloJetEtav[0]->GetBinWidth(1));
//   CPlot plotCaloJetEta("calojeteta","Jet E_{T} > 15 GeV","Calo Jet #eta",ylabel);
//   if(hasData) { plotCaloJetEta.AddHist1D(hCaloJetEtav[0],samplev[0]->label,"E"); }
//   for(UInt_t isam=1; isam<samplev.size(); isam++) {
//     hCaloJetEtav[isam]->Scale(mcscale);
//     plotCaloJetEta.AddToStack(hCaloJetEtav[isam],samplev[isam]->label,samplev[isam]->color);
//   }
//   if(lumi>0) plotCaloJetEta.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
//   if(samplev.size()>5) 
//     plotCaloJetEta.SetLegend(0.7,0.55,0.9,0.9);
//   else 
//     plotCaloJetEta.TransLegend(0.1,0);  
//   plotCaloJetEta.SetYRange(0,1.5*(plotCaloJetEta.GetStack()->GetMaximum()));
//   plotCaloJetEta.Draw(c,kFALSE,format);
    
//   sprintf(ylabel,"Events / %.2f",hCaloJetPhiv[0]->GetBinWidth(1));
//   CPlot plotCaloJetPhi("calojetphi","Jet E_{T} > 15 GeV","Calo Jet #phi",ylabel);
//   if(hasData) { plotCaloJetPhi.AddHist1D(hCaloJetPhiv[0],samplev[0]->label,"E"); }
//   for(UInt_t isam=1; isam<samplev.size(); isam++) {
//     hCaloJetPhiv[isam]->Scale(mcscale);
//     plotCaloJetPhi.AddToStack(hCaloJetPhiv[isam],samplev[isam]->label,samplev[isam]->color);
//   }
//   if(lumi>0) plotCaloJetPhi.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
//   plotCaloJetPhi.SetYRange(0,1.8*(plotCaloJetPhi.GetStack()->GetMaximum()));
//   if(samplev.size()>5) 
//     plotCaloJetPhi.SetLegend(0.7,0.55,0.9,0.9);
//   else 
//     plotCaloJetPhi.TransLegend(0.1,0);
//   plotCaloJetPhi.Draw(c,kFALSE,format);

  //
  // Track Jet kinematics
  //
//   sprintf(ylabel,"Events / %.1f GeV",hTrackJetEtv[0]->GetBinWidth(1));
//   CPlot plotTrackJetEt("trackjetet","Jet E_{T} > 15 GeV","Track Jet E_{T} [GeV]",ylabel);
//   if(hasData) { plotTrackJetEt.AddHist1D(hTrackJetEtv[0],samplev[0]->label,"E"); }
//   for(UInt_t isam=1; isam<samplev.size(); isam++) {
//     hTrackJetEtv[isam]->Scale(mcscale);
//     plotTrackJetEt.AddToStack(hTrackJetEtv[isam],samplev[isam]->label,samplev[isam]->color);
//   }
//   if(lumi>0) plotTrackJetEt.AddTextBox(lumitext,0.41,0.85,0.61,0.8,0);
//   if(samplev.size()>5) 
//     plotTrackJetEt.SetLegend(0.7,0.55,0.9,0.9);
//   else 
//     plotTrackJetEt.TransLegend(0.1,0);  
//   plotTrackJetEt.Draw(c,kFALSE,format);
  
//   plotTrackJetEt.SetName("trackjetetlog");  
//   plotTrackJetEt.SetLogy();
//   plotTrackJetEt.SetYRange((1e-4)*(plotTrackJetEt.GetStack()->GetMaximum()),10.*(plotTrackJetEt.GetStack()->GetMaximum()));  
//   plotTrackJetEt.Draw(c,kFALSE,format);
    
//   sprintf(ylabel,"Events / %.2f",hTrackJetEtav[0]->GetBinWidth(1));
//   CPlot plotTrackJetEta("trackjeteta","Jet E_{T} > 15 GeV","Track Jet #eta",ylabel);
//   if(hasData) { plotTrackJetEta.AddHist1D(hTrackJetEtav[0],samplev[0]->label,"E"); }
//   for(UInt_t isam=1; isam<samplev.size(); isam++) {
//     hTrackJetEtav[isam]->Scale(mcscale);
//     plotTrackJetEta.AddToStack(hTrackJetEtav[isam],samplev[isam]->label,samplev[isam]->color);
//   }
//   if(lumi>0) plotTrackJetEta.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
//   if(samplev.size()>5) 
//     plotTrackJetEta.SetLegend(0.7,0.55,0.9,0.9);
//   else 
//     plotTrackJetEta.TransLegend(0.1,0);
//   plotTrackJetEta.SetYRange(0,1.5*(plotTrackJetEta.GetStack()->GetMaximum()));
//   plotTrackJetEta.Draw(c,kFALSE,format);
    
//   sprintf(ylabel,"Events / %.2f",hTrackJetPhiv[0]->GetBinWidth(1));
//   CPlot plotTrackJetPhi("trackjetphi","Jet E_{T} > 15 GeV","Track Jet #phi",ylabel);
//   if(hasData) { plotTrackJetPhi.AddHist1D(hTrackJetPhiv[0],samplev[0]->label,"E"); }
//   for(UInt_t isam=1; isam<samplev.size(); isam++) {
//     hTrackJetPhiv[isam]->Scale(mcscale);
//     plotTrackJetPhi.AddToStack(hTrackJetPhiv[isam],samplev[isam]->label,samplev[isam]->color);
//   }
//   if(lumi>0) plotTrackJetPhi.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
//   plotTrackJetPhi.SetYRange(0,1.8*(plotTrackJetPhi.GetStack()->GetMaximum()));
//   if(samplev.size()>5) 
//     plotTrackJetPhi.SetLegend(0.7,0.55,0.9,0.9);
//   else 
//     plotTrackJetPhi.TransLegend(0.1,0);
//   plotTrackJetPhi.Draw(c,kFALSE,format);

  //
  // PF Jet kinematics
  //
  sprintf(ylabel,"Events / %.1f GeV",hPFJetEtv[0]->GetBinWidth(1));
  CPlot plotPFJetEt("pfjetet","Jet E_{T} > 15 GeV","PF Jet E_{T} [GeV]",ylabel);
  if(hasData) { plotPFJetEt.AddHist1D(hPFJetEtv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    hPFJetEtv[isam]->Scale(mcscale);
    plotPFJetEt.AddToStack(hPFJetEtv[isam],samplev[isam]->label,samplev[isam]->color);
  }
  if(lumi>0) plotPFJetEt.AddTextBox(lumitext,0.41,0.85,0.61,0.8,0);
  if(samplev.size()>5) 
    plotPFJetEt.SetLegend(0.7,0.55,0.9,0.9);
  else 
    plotPFJetEt.TransLegend(0.1,0);  
  plotPFJetEt.Draw(c,kFALSE,format);
  
  plotPFJetEt.SetName("pfjetetlog");  
  plotPFJetEt.SetLogy();
  plotPFJetEt.SetYRange((1e-4)*(plotPFJetEt.GetStack()->GetMaximum()),10.*(plotPFJetEt.GetStack()->GetMaximum()));  
  plotPFJetEt.Draw(c,kFALSE,format);
    
  sprintf(ylabel,"Events / %.2f",hPFJetEtav[0]->GetBinWidth(1));
  CPlot plotPFJetEta("pfjeteta","Jet E_{T} > 15 GeV","PF Jet #eta",ylabel);
  if(hasData) { plotPFJetEta.AddHist1D(hPFJetEtav[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    hPFJetEtav[isam]->Scale(mcscale);
    plotPFJetEta.AddToStack(hPFJetEtav[isam],samplev[isam]->label,samplev[isam]->color);
  }
  if(lumi>0) plotPFJetEta.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
  if(samplev.size()>5) 
    plotPFJetEta.SetLegend(0.7,0.55,0.9,0.9);
  else 
    plotPFJetEta.TransLegend(0.1,0);
  plotPFJetEta.SetYRange(0,1.5*(plotPFJetEta.GetStack()->GetMaximum()));
  plotPFJetEta.Draw(c,kFALSE,format);
    
  sprintf(ylabel,"Events / %.2f",hPFJetPhiv[0]->GetBinWidth(1));
  CPlot plotPFJetPhi("pfjetphi","Jet E_{T} > 15 GeV","PF Jet #phi",ylabel);
  if(hasData) { plotPFJetPhi.AddHist1D(hPFJetPhiv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    hPFJetPhiv[isam]->Scale(mcscale);
    plotPFJetPhi.AddToStack(hPFJetPhiv[isam],samplev[isam]->label,samplev[isam]->color);
  }
  if(lumi>0) plotPFJetPhi.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
  plotPFJetPhi.SetYRange(0,1.8*(plotPFJetPhi.GetStack()->GetMaximum()));
  if(samplev.size()>5) 
    plotPFJetPhi.SetLegend(0.7,0.55,0.9,0.9);
  else 
    plotPFJetPhi.TransLegend(0.1,0);
  plotPFJetPhi.Draw(c,kFALSE,format);
  

  //
  // Primary Vertex
  //
  CPlot plotNPV("npv","","N_{PV}","Events");
  if(hasData) { plotNPV.AddHist1D(hNPVv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    hNPVv[isam]->Scale(mcscale);
    plotNPV.AddToStack(hNPVv[isam],samplev[isam]->label,samplev[isam]->color);
  }
  if(lumi>0) plotNPV.AddTextBox(lumitext,0.45,0.85,0.65,0.8,0);
  plotNPV.TransLegend(0.1,0);
  plotNPV.SetLogy();
  plotNPV.Draw(c,kFALSE,format);

  CPlot plotNGoodPV("ngoodpv","Good PVs","N_{PV}","Events");
  if(hasData) { plotNGoodPV.AddHist1D(hNGoodPVv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    hNGoodPVv[isam]->Scale(mcscale);
    plotNGoodPV.AddToStack(hNGoodPVv[isam],samplev[isam]->label,samplev[isam]->color);
  }
  if(lumi>0) plotNGoodPV.AddTextBox(lumitext,0.45,0.85,0.65,0.8,0);
  plotNGoodPV.TransLegend(0.1,0);
  plotNGoodPV.SetLogy();
  plotNGoodPV.Draw(c,kFALSE,format);

  sprintf(ylabel,"Events / %i",(Int_t)hGoodPVNTracksv[0]->GetBinWidth(1));    
  CPlot plotGoodPVNTracks("goodpvntracks","Good PVs","PV N_{tracks}",ylabel);
  if(hasData) { plotGoodPVNTracks.AddHist1D(hGoodPVNTracksv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    hGoodPVNTracksv[isam]->Scale(mcscale);
    plotGoodPVNTracks.AddToStack(hGoodPVNTracksv[isam],samplev[isam]->label,samplev[isam]->color);
  }
  if(lumi>0) plotGoodPVNTracks.AddTextBox(lumitext,0.41,0.85,0.61,0.8,0);
  plotGoodPVNTracks.TransLegend(0.08,0);
  plotGoodPVNTracks.Draw(c,kFALSE,format);
  
  sprintf(ylabel,"Events / %.1f",hGoodPVChi2v[0]->GetBinWidth(1));
  CPlot plotGoodPVChi2("goodpvchi2","Good PVs","PV #chi^{2}",ylabel);
  if(hasData) { plotGoodPVChi2.AddHist1D(hGoodPVChi2v[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    hGoodPVChi2v[isam]->Scale(mcscale);
    plotGoodPVChi2.AddToStack(hGoodPVChi2v[isam],samplev[isam]->label,samplev[isam]->color);
  }
  if(lumi>0) plotGoodPVChi2.AddTextBox(lumitext,0.41,0.85,0.61,0.8,0);
  plotGoodPVChi2.TransLegend(0.08,0);
  plotGoodPVChi2.Draw(c,kFALSE,format);

  sprintf(ylabel,"Events / %i",(Int_t)hGoodPVNdofv[0]->GetBinWidth(1));
  CPlot plotGoodPVNdof("goodpvndof","Good PVs","PV N_{DOF}",ylabel);
  if(hasData) { plotGoodPVNdof.AddHist1D(hGoodPVNdofv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    hGoodPVNdofv[isam]->Scale(mcscale);
    plotGoodPVNdof.AddToStack(hGoodPVNdofv[isam],samplev[isam]->label,samplev[isam]->color);
  }
  if(lumi>0) plotGoodPVNdof.AddTextBox(lumitext,0.41,0.85,0.61,0.8,0);
  plotGoodPVNdof.TransLegend(0.08,0);
  plotGoodPVNdof.Draw(c,kFALSE,format);
  
  sprintf(ylabel,"Events / %.1f GeV/c",hGoodPVSumPtv[0]->GetBinWidth(1));
  CPlot plotGoodPVSumPt("goodpvsumpt","Good PVs","PV #Sigma^{}p_{T} [GeV/c]",ylabel);
  if(hasData) { plotGoodPVSumPt.AddHist1D(hGoodPVSumPtv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    hGoodPVSumPtv[isam]->Scale(mcscale);
    plotGoodPVSumPt.AddToStack(hGoodPVSumPtv[isam],samplev[isam]->label,samplev[isam]->color);
  }
  if(lumi>0) plotGoodPVSumPt.AddTextBox(lumitext,0.41,0.85,0.61,0.8,0);
  plotGoodPVSumPt.TransLegend(0.08,0);
  plotGoodPVSumPt.Draw(c,kFALSE,format);

  if(hasData) {
    sprintf(ylabel,"Events / %.1f",hGoodPVDz->GetBinWidth(1));
    CPlot plotGoodPVDz("goodpvdz","Good PV pairs","#Deltaz [cm]",ylabel);
    plotGoodPVDz.AddHist1D(hGoodPVDz);
    if(lumi>0) plotGoodPVDz.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
    plotGoodPVDz.ShowStats();
    plotGoodPVDz.Draw(c,kFALSE,format);
  }  
      
  //
  // Track Multiplicity
  //
  sprintf(ylabel,"Events / %.0f GeV",hNTracks0v[0]->GetBinWidth(1));
  CPlot plotNTracks0("ntracks0","","N_{tracks}",ylabel);
  if(hasData) { plotNTracks0.AddHist1D(hNTracks0v[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    hNTracks0v[isam]->Scale(mcscale);
    plotNTracks0.AddToStack(hNTracks0v[isam],samplev[isam]->label,samplev[isam]->color);
  }
  if(lumi>0) plotNTracks0.AddTextBox(lumitext,0.45,0.85,0.65,0.8,0);
  plotNTracks0.TransLegend(0.1,0);
  plotNTracks0.Draw(c,kFALSE,format);

  //
  // Calorimeter Tower Multiplicity
  //
  sprintf(ylabel,"Events / %.0f GeV",hNCaloTowers0v[0]->GetBinWidth(1));
  CPlot plotNCaloTowers0("ncalotowers0","","N_{calo towers}",ylabel);
  if(hasData) { plotNCaloTowers0.AddHist1D(hNCaloTowers0v[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    hNCaloTowers0v[isam]->Scale(mcscale);
    plotNCaloTowers0.AddToStack(hNCaloTowers0v[isam],samplev[isam]->label,samplev[isam]->color);
  }
  if(lumi>0) plotNCaloTowers0.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
  plotNCaloTowers0.TransLegend(0.08,0);
  plotNCaloTowers0.Draw(c,kFALSE,format);


  // DeltaEta_In (barrel)
  sprintf(ylabel,"Events / %.3f",hDeltaEtaInBv[0]->GetBinWidth(1));
  CPlot plotDeltaEtaInB("deltaetainB","barrel","#Delta^{}#eta_{in}",ylabel);
  if(hasData) { plotDeltaEtaInB.AddHist1D(hDeltaEtaInBv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    hDeltaEtaInBv[isam]->Scale(mcscale);
    plotDeltaEtaInB.AddToStack(hDeltaEtaInBv[isam],samplev[isam]->label,samplev[isam]->color);
  }
  if(lumi>0) plotDeltaEtaInB.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
  plotDeltaEtaInB.TransLegend(0.1,0);
  plotDeltaEtaInB.Draw(c,kFALSE,format);
    
  // DeltaPhi_In (barrel)
  sprintf(ylabel,"Events / %.3f",hDeltaPhiInBv[0]->GetBinWidth(1));
  CPlot plotDeltaPhiInB("deltaphiinB","barrel","#Delta^{}#phi_{in}",ylabel);  
  if(hasData) { plotDeltaPhiInB.AddHist1D(hDeltaPhiInBv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    hDeltaPhiInBv[isam]->Scale(mcscale);
    plotDeltaPhiInB.AddToStack(hDeltaPhiInBv[isam],samplev[isam]->label,samplev[isam]->color);
  }
  if(lumi>0) plotDeltaPhiInB.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
  plotDeltaPhiInB.TransLegend(0.05,0);
  plotDeltaPhiInB.Draw(c,kFALSE,format);

  // DeltaEta_In (endcap)
  sprintf(ylabel,"Events / %.3f",hDeltaEtaInEv[0]->GetBinWidth(1));
  CPlot plotDeltaEtaInE("deltaetainE","endcap","#Delta^{}#eta_{in}",ylabel);
  if(hasData) { plotDeltaEtaInE.AddHist1D(hDeltaEtaInEv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    hDeltaEtaInEv[isam]->Scale(mcscale);
    plotDeltaEtaInE.AddToStack(hDeltaEtaInEv[isam],samplev[isam]->label,samplev[isam]->color);
  }
  if(lumi>0) plotDeltaEtaInE.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
  plotDeltaEtaInE.TransLegend(0.1,0);
  plotDeltaEtaInE.Draw(c,kFALSE,format);
    
  // DeltaPhi_In (endcap)
  sprintf(ylabel,"Events / %.3f",hDeltaPhiInEv[0]->GetBinWidth(1));
  CPlot plotDeltaPhiInE("deltaphiinE","endcap","#Delta^{}#phi_{in}",ylabel);  
  if(hasData) { plotDeltaPhiInE.AddHist1D(hDeltaPhiInEv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    hDeltaPhiInEv[isam]->Scale(mcscale);
    plotDeltaPhiInE.AddToStack(hDeltaPhiInEv[isam],samplev[isam]->label,samplev[isam]->color);
  }
  if(lumi>0) plotDeltaPhiInE.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
  plotDeltaPhiInE.TransLegend(0.05,0);
  plotDeltaPhiInE.Draw(c,kFALSE,format);

  } // end if zillionPlots

  //--------------------------------------------------------------------------------------------------------------
  // Summary print out
  //============================================================================================================== 
  ofstream txtfile;
  char txtfname[100];    
  sprintf(txtfname,"%s/summary.txt",outputDir.Data());
  txtfile.open(txtfname);
  assert(txtfile.is_open());
  txtfile << "*" << endl;
  txtfile << "* SUMMARY" << endl;
  txtfile << "*--------------------------------------------------" << endl;
  txtfile << endl;

  txtfile << "  L_int = " << lumi << "/pb" << endl;
  txtfile << endl;
  
  if(hasData) {
    txtfile << "   Data: " << setprecision(1) << fixed << nProcessedEvents << " events processed!" << endl;
    txtfile << "         " << setprecision(1) << fixed << nSelv[0] << " Z events!" << endl;
    txtfile << "         " << setprecision(1) << fixed << nPosSSv[0] << " SS (+) events!" << endl;
    txtfile << "         " << setprecision(1) << fixed << nNegSSv[0] << " SS (-) events!" << endl;
    for(UInt_t ifile=0; ifile<samplev[0]->fnamev.size(); ifile++)
      txtfile << "     " << samplev[0]->fnamev[ifile] << endl;
      txtfile << endl;
  } 
  
  if(samplev.size()>1) {
    txtfile << "   MC:" << endl;
    for(UInt_t isam=1; isam<samplev.size(); isam++) {      
      for(UInt_t ifile=0; ifile<samplev[isam]->fnamev.size(); ifile++) {
        if(ifile==0) {
          txtfile << setw(10) << snamev[isam];
          txtfile << setw(10) << setprecision(2) << fixed << nSelv[isam] << " +/- ";
          txtfile << setw(5) << setprecision(2) << fixed << sqrt(nSelVarv[isam]);
          txtfile << "   " << "SS (+) = " << setw(5) << setprecision(3) << nPosSSv[isam];
	  txtfile << "   " << "SS (-) = " << setw(5) << setprecision(3) << nNegSSv[isam];
          txtfile << "   " << samplev[isam]->fnamev[ifile] << endl;
        } else {
          txtfile << setw(48) << "" << "   " << samplev[isam]->fnamev[ifile] << endl;
        }
      }
      txtfile << endl;
    }
  }
  txtfile.close();

  // make webpage
  makeHTML(outputDir);
         
  cout << endl;
  cout << " <> Output saved in " << outputDir << "/" << endl;
  cout << endl;
        
  gBenchmark->Show("plotDY");       
} 


//=== FUNCTION DEFINITIONS ======================================================================================

//--------------------------------------------------------------------------------------------------
void makeHTML(const TString outDir)
{
  ofstream htmlfile;
  char htmlfname[100];
  sprintf(htmlfname,"%s/plots.html",outDir.Data());
  htmlfile.open(htmlfname);
  htmlfile << "<!DOCTYPE html" << endl;
  htmlfile << "    PUBLIC \"-//W3C//DTD HTML 3.2//EN\">" << endl;
  htmlfile << "<html>" << endl;
  htmlfile << "<head><title>Zee</title></head>" << endl;
  htmlfile << "<body bgcolor=\"EEEEEE\">" << endl;
  htmlfile << "<h3 style=\"text-align:left; color:DD6600;\">Zee</h3>" << endl;

  htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl;  
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/mass.png\"><img src=\"plots/mass.png\" alt=\"plots/mass.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/masslog.png\"><img src=\"plots/masslog.png\" alt=\"plots/masslog.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/mass2.png\"><img src=\"plots/mass2.png\" alt=\"plots/mass2.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/mass3.png\"><img src=\"plots/mass3.png\" alt=\"plots/mass3.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/massBB.png\"><img src=\"plots/massBB.png\" alt=\"plots/massBB.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/massBE.png\"><img src=\"plots/massBE.png\" alt=\"plots/massBE.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/massEE.png\"><img src=\"plots/massEE.png\" alt=\"plots/massEE.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/pt.png\"><img src=\"plots/pt.png\" alt=\"plots/pt.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/ptlog.png\"><img src=\"plots/ptlog.png\" alt=\"plots/ptlog.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/y.png\"><img src=\"plots/y.png\" alt=\"plots/y.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/phi.png\"><img src=\"plots/phi.png\" alt=\"plots/phi.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "</table>" << endl;
  htmlfile << "<hr />" << endl;

  htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl;  
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/elept.png\"><img src=\"plots/elept.png\" alt=\"plots/elept.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/eleptlog.png\"><img src=\"plots/eleptlog.png\" alt=\"plots/eleptlog.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/eleeta.png\"><img src=\"plots/eleeta.png\" alt=\"plots/eleeta.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/elephi.png\"><img src=\"plots/elephi.png\" alt=\"plots/elephi.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "</tr>" << endl; 
  htmlfile << "</table>" << endl;
  htmlfile << "<hr />" << endl;

  htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl;  
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/calomex.png\"><img src=\"plots/calomex.png\" alt=\"plots/calomex.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/tcmex.png\"><img src=\"plots/tcmex.png\" alt=\"plots/tcmex.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/pfmex.png\"><img src=\"plots/pfmex.png\" alt=\"plots/pfmex.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/calomey.png\"><img src=\"plots/calomey.png\" alt=\"plots/calomey.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/tcmey.png\"><img src=\"plots/tcmey.png\" alt=\"plots/tcmey.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/pfmey.png\"><img src=\"plots/pfmey.png\" alt=\"plots/pfmey.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/calomet.png\"><img src=\"plots/calomet.png\" alt=\"plots/calomet.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/tcmet.png\"><img src=\"plots/tcmet.png\" alt=\"plots/tcmet.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/pfmet.png\"><img src=\"plots/pfmet.png\" alt=\"plots/pfmet.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/calometphi.png\"><img src=\"plots/calometphi.png\" alt=\"plots/calometphi.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/tcmetphi.png\"><img src=\"plots/tcmetphi.png\" alt=\"plots/tcmetphi.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/pfmetphi.png\"><img src=\"plots/pfmetphi.png\" alt=\"plots/pfmetphi.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl; 
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/calosumet.png\"><img src=\"plots/calosumet.png\" alt=\"plots/calosumet.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/tcsumet.png\"><img src=\"plots/tcsumet.png\" alt=\"plots/tcsumet.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/pfsumet.png\"><img src=\"plots/pfsumet.png\" alt=\"plots/pfsumet.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl; 
  htmlfile << "</table>" << endl;
  htmlfile << "<hr />" << endl;

  htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl;  
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/ntracks0.png\"><img src=\"plots/ntracks0.png\" alt=\"plots/ntracks0.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/ncalotowers0.png\"><img src=\"plots/ncalotowers0.png\" alt=\"plots/ncalotowers0.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/trigger.png\"><img src=\"plots/trigger.png\" alt=\"plots/trigger.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/npv.png\"><img src=\"plots/npv.png\" alt=\"plots/npv.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/ngoodpv.png\"><img src=\"plots/ngoodpv.png\" alt=\"plots/ngoodpv.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/goodpvdz.png\"><img src=\"plots/goodpvdz.png\" alt=\"plots/goodpvdz.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/goodpvntracks.png\"><img src=\"plots/goodpvntracks.png\" alt=\"plots/goodpvntracks.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/goodpvchi2.png\"><img src=\"plots/goodpvchi2.png\" alt=\"plots/goodpvchi2.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/goodpvndof.png\"><img src=\"plots/goodpvndof.png\" alt=\"plots/goodpvndof.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/goodpvsumpt.png\"><img src=\"plots/goodpvsumpt.png\" alt=\"plots/goodpvsumpt.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "</tr>" << endl; 
  htmlfile << "<tr>" << endl;
//   htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/ncalojets.png\"><img src=\"plots/ncalojets.png\" alt=\"plots/ncalojets.png\" width=\"100%\"></a></td>" << endl;
//   htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/ntrackjets.png\"><img src=\"plots/ntrackjets.png\" alt=\"plots/ntrackjets.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/npfjets.png\"><img src=\"plots/npfjets.png\" alt=\"plots/npfjets.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl; 
  htmlfile << "</table>" << endl;
  htmlfile << "<hr />" << endl;  

  htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl;  
  htmlfile << "<tr>" << endl;
//   htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/calojetet.png\"><img src=\"plots/calojetet.png\" alt=\"plots/calojetet.png\" width=\"100%\"></a></td>" << endl;
//   htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/calojetetlog.png\"><img src=\"plots/calojetetlog.png\" alt=\"plots/calojetetlog.png\" width=\"100%\"></a></td>" << endl;
//   htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/calojeteta.png\"><img src=\"plots/calojeteta.png\" alt=\"plots/calojeteta.png\" width=\"100%\"></a></td>" << endl;
//   htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/calojetphi.png\"><img src=\"plots/calojetphi.png\" alt=\"plots/calojetphi.png\" width=\"100%\"></a></td>" << endl;
//   htmlfile << "</tr>" << endl;
//   htmlfile << "<tr>" << endl;
//   htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/trackjetet.png\"><img src=\"plots/trackjetet.png\" alt=\"plots/trackjetet.png\" width=\"100%\"></a></td>" << endl;
//   htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/trackjetetlog.png\"><img src=\"plots/trackjetetlog.png\" alt=\"plots/trackjetetlog.png\" width=\"100%\"></a></td>" << endl;
//   htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/trackjeteta.png\"><img src=\"plots/trackjeteta.png\" alt=\"plots/trackjeteta.png\" width=\"100%\"></a></td>" << endl;
//   htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/trackjetphi.png\"><img src=\"plots/trackjetphi.png\" alt=\"plots/trackjetphi.png\" width=\"100%\"></a></td>" << endl;
//   htmlfile << "</tr>" << endl; 
//   htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/pfjetet.png\"><img src=\"plots/pfjetet.png\" alt=\"plots/pfjetet.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/pfjetetlog.png\"><img src=\"plots/pfjetetlog.png\" alt=\"plots/pfjetetlog.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/pfjeteta.png\"><img src=\"plots/pfjeteta.png\" alt=\"plots/pfjeteta.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/pfjetphi.png\"><img src=\"plots/pfjetphi.png\" alt=\"plots/pfjetphi.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "</table>" << endl; 
  htmlfile << "<hr />" << endl;
  
  htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl;  
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/deltaetainB.png\"><img src=\"plots/deltaetainB.png\" alt=\"plots/deltaetainB.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/deltaetainE.png\"><img src=\"plots/deltaetainE.png\" alt=\"plots/deltaetainE.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/deltaphiinB.png\"><img src=\"plots/deltaphiinB.png\" alt=\"plots/deltaphiinB.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/deltaphiinE.png\"><img src=\"plots/deltaphiinE.png\" alt=\"plots/deltaphiinE.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "</tr>" << endl; 
  htmlfile << "</table>" << endl;
  htmlfile << "<hr />" << endl;
      
  htmlfile << "</body>" << endl;
  htmlfile << "</html>" << endl;
  htmlfile.close();
}

//--------------------------------------------------------------------------------------------------
void fillData(ZeeData *data, const mithep::TEventInfo *info, const mithep::TDielectron *dielectron, 
              const UInt_t npv, const UInt_t njets, const Double_t weight)
{
  data->runNum         = info->runNum;
  data->evtNum         = info->evtNum;
  data->lumiSec        = info->lumiSec;
//   data->nTracks0       = info->nTracks0;
//   data->nCaloTowers0   = info->nCaloTowers0;
  data->nPV            = npv;
  data->nJets          = njets;                                        
//   data->caloMEx        = info->caloMEx;
//   data->caloMEy        = info->caloMEy;
//   data->caloSumET      = info->caloSumET;
//   data->tcMEx          = info->tcMEx;
//   data->tcMEy          = info->tcMEy;
//   data->tcSumET        = info->tcSumET;
//   data->pfMEx          = info->pfMEx;
//   data->pfMEy          = info->pfMEy;
  data->pfSumET        = info->pfSumET;
  data->mass           = dielectron->mass;
  data->pt             = dielectron->pt;
  data->y              = dielectron->y;
  data->phi            = dielectron->phi; 
  data->pt_1           = dielectron->pt_1;
  data->eta_1          = dielectron->eta_1;
  data->phi_1          = dielectron->phi_1;
  data->scEt_1         = dielectron->scEt_1;
  data->scEta_1        = dielectron->scEta_1;
  data->scPhi_1        = dielectron->scPhi_1;
  data->hltMatchBits_1 = dielectron->hltMatchBits_1;
  data->q_1            = dielectron->q_1;
  data->pt_2           = dielectron->pt_2;
  data->eta_2          = dielectron->eta_2;
  data->phi_2          = dielectron->phi_2;
  data->scEt_2         = dielectron->scEt_2;
  data->scEta_2        = dielectron->scEta_2;
  data->scPhi_2        = dielectron->scPhi_2;
  data->hltMatchBits_2 = dielectron->hltMatchBits_2;
  data->q_2            = dielectron->q_2;
  data->weight         = weight;
}

//--------------------------------------------------------------------------------------------------
void eventDump(ofstream &ofs, const mithep::TDielectron *dielectron, 
               const UInt_t runNum, const UInt_t lumiSec, const UInt_t evtNum, 
	       const UInt_t triggerObj1, const UInt_t triggerObj2)
{
  ofs << endl;
  ofs << "Run:" << runNum;
  ofs << "  Lumi:" << lumiSec;
  ofs << "  Event:" << evtNum;
  ofs << "  mass: " << dielectron->mass;
  ofs << "  pt: " << dielectron->pt << endl;
  
  ofs << "----------+-----------+-----------+-----------+-----------+-----------+-------------+------------+------------+-----------+------" << endl;
  ofs << "  SC ET   |  SC eta   |   SC phi  | trkiso/pt | emiso/pt  | hadiso/pt | sigiEtaiEta |    deta    |    dphi    |    H/E    | HLT" << endl;
  ofs << "----------+-----------+-----------+-----------+-----------+-----------+-------------+------------+------------+-----------+------" << endl;
      
  ofs << setw(9) << dielectron->scEt_1 << " |";
  ofs << setw(10) << dielectron->scEta_1 << " |";
  ofs << setw(10) << dielectron->scPhi_1 << " |";
  ofs << setw(10) << dielectron->trkIso03_1/dielectron->pt_1 << " |";
  ofs << setw(10) << dielectron->emIso03_1/dielectron->pt_1 << " |";
  ofs << setw(10) << dielectron->hadIso03_1/dielectron->pt_1 << " |";
  ofs << setw(12) << dielectron->sigiEtaiEta_1 << " |";
  ofs << setw(12) << dielectron->deltaEtaIn_1 << "|";
  ofs << setw(12) << dielectron->deltaPhiIn_1 << "|";
  ofs << setw(10) << dielectron->HoverE_1 << " |";
  if(dielectron->hltMatchBits_1 & triggerObj1)
    ofs << " LEAD" << endl; 
  else if(dielectron->hltMatchBits_1 & triggerObj2)
    ofs << " TRAIL" << endl;
  else
    ofs << " NOMAT" << endl;
    
  ofs << setw(9) << dielectron->scEt_2 << " |";
  ofs << setw(10) << dielectron->scEta_2 << " |";
  ofs << setw(10) << dielectron->scPhi_2 << " |";
  ofs << setw(10) << dielectron->trkIso03_2/dielectron->pt_2 << " |";
  ofs << setw(10) << dielectron->emIso03_2/dielectron->pt_2 << " |";
  ofs << setw(10) << dielectron->hadIso03_2/dielectron->pt_2 << " |";
  ofs << setw(12) << dielectron->sigiEtaiEta_2 << " |";
  ofs << setw(12) << dielectron->deltaEtaIn_2 << "|";
  ofs << setw(12) << dielectron->deltaPhiIn_2 << "|";
  ofs << setw(10) << dielectron->HoverE_2 << " |";
  if(dielectron->hltMatchBits_2 & triggerObj1)
    ofs << " LEAD" << endl; 
  else if(dielectron->hltMatchBits_2 & triggerObj2)
    ofs << " TRAIL" << endl;
  else
    ofs << " NOMAT" << endl;
}    

