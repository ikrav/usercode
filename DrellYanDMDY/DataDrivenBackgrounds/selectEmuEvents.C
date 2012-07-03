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

#include "../Include/CPlot.hh"          // helper class for plots
#include "../Include/MitStyleRemix.hh"  // style settings for drawing
#include "../Include/MyTools.hh"        // miscellaneous helper functions
#include "../Include/CSample.hh"        // helper class for organizing input ntuple files
#include "../Include/DYTools.hh"
#include "../Include/TriggerSelection.hh"

// define structures to read in ntuple
#include "../Include/EWKAnaDefs.hh"
#include "../Include/TEventInfo.hh"
#include "../Include/TGenInfo.hh"
#include "../Include/TElectron.hh"
#include "../Include/TMuon.hh"
#include "../Include/TVertex.hh"

// lumi section selection with JSON files
#include "../Include/JsonParser.hh"

// Helper functions for Electron ID selection
#include "../Include/EleIDCuts.hh" 

// define structure for output ntuple
#include "../Include/XemuData.hh"

#include "../Include/ElectronEnergyScale.hh" //energy scale correction
#include "../Include/EtaEtaMass.hh" // EtaEtaMassData_t definition

#endif

//=== FUNCTION DECLARATIONS ======================================================================================

// fill ntuple of selected events
void fillData(XemuData *data, const mithep::TEventInfo *info, const mithep::TElectron *electron, 
	      const mithep::TMuon *muon, const TLorentzVector &emu4v,
              const UInt_t npv, const UInt_t njets, const Double_t weight);

// print event dump
void eventDump(ofstream &ofs, const mithep::TElectron *electron, 
               const UInt_t runNum, const UInt_t lumiSec, const UInt_t evtNum, 
	       const Float_t &mass,
	       const UInt_t triggerObj1, const UInt_t triggerObj2);


template < typename T > inline T highbit(T& t);
template < typename T > std::ostream& bin(T& value, std::ostream &o); 

//=== MAIN MACRO =================================================================================================

void selectEmuEvents(const TString conf) 
{  
  gBenchmark->Start("selectEvents");

  
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
	std::cout << "\n\tEEM files will be generated, tag=<" << generateEEMFile << ">, with_FEWZ_weights=" << generateEEMFEWZFile << "\n\n";
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


  // 
  // Set up energy scale corrections
  //
  ElectronEnergyScale escale(escaleTag);
  assert(escale.isInitialized());
  escale.print();


  // sOutDir is a static data member in the CPlot class.
  // There is a strange crash of the whole ROOT session well after
  // this script is executed when one attempts to exit ROOT, with 
  // a dump of memory map. This happens only on UNL Tier3, but
  // there is no sign of a problem on any other computer.
  //   The consequence of this variable is not set is that the plots
  // will be created in the local directory rather than the
  // one configured through sOutDir.
//   CPlot::sOutDir        = outputDir + TString("/plots");   gSystem->mkdir(CPlot::sOutDir,kTRUE);

  const TString ntupDir = outputDir + TString("/ntuples/EMU"); gSystem->mkdir(ntupDir,kTRUE);
  
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
  vector<TH1F*> hMassv, hMass2v, hMass3v, hMass4v;
  vector<TH1F*> hyv;   

  vector<TH1F*> hNGoodPVv;
  
  vector<Double_t> nSelv, nSelVarv;  
  vector<Double_t> nPosSSv;
  vector<Double_t> nNegSSv;
    
  UInt_t nProcessedEvents=0;
//   TH1F* hTrigger = new TH1F("hTrigger","",32,-0.5,31.5);
  
  char hname[100];
  for(UInt_t isam=0; isam<samplev.size(); isam++) {

    sprintf(hname,"hMass_%i",isam);   hMassv.push_back(new TH1F(hname,"",30,60,120));   hMassv[isam]->Sumw2();
    sprintf(hname,"hMass2_%i",isam);  hMass2v.push_back(new TH1F(hname,"",35,20,160));  hMass2v[isam]->Sumw2();
    sprintf(hname,"hMass3_%i",isam);  hMass3v.push_back(new TH1F(hname,"",50,0,500));   hMass3v[isam]->Sumw2();
    sprintf(hname,"hy_%i",isam);      hyv.push_back(new TH1F(hname,"",20,-3,3));        hyv[isam]->Sumw2();
    sprintf(hname,"hNGoodPV_%s",snamev[isam].Data());       hNGoodPVv.push_back(new TH1F(hname,"",15,-0.5,14.5));        hNGoodPVv[isam]->Sumw2();
    sprintf(hname,"hMass4_%i",isam);   hMass4v.push_back(new TH1F(hname,"",275,0,1100));   hMass4v[isam]->Sumw2();
    
    nSelv.push_back(0);
    nSelVarv.push_back(0);
    nPosSSv.push_back(0);
    nNegSSv.push_back(0);    
  }
  
  // 
  // Read weights from a file
  //
  //const bool useFewzWeights = true;
  const bool cutZPT100 = true;
  if(cutZPT100)
    cout << "NOTE: in MC, for Z/gamma* PT>100 the FEWZ weights for 80<PT<100 GeV are used!" << endl;
  TH2D *weights[DYTools::nMassBins];
  TH2D *weightErrors[DYTools::nMassBins];
  TFile fweights("../root_files/fewz/weights_stepwise_prec10-5_fine12.root");
  if( !fweights.IsOpen() ) assert(0);
  for(int i=0; i<DYTools::nMassBins; i++){
    TString hnames = TString::Format("weight_%02d",i+1);
    weights[i] = (TH2D*)fweights.Get(hnames);
    hnames = TString::Format("h_weighterror_%02d",i+1);
    weightErrors[i] = (TH2D*)fweights.Get(hnames);
  }


  //
  // Access samples and fill histograms
  //  
  TFile *infile=0;
  TTree *eventTree=0;  
  
  // Data structures to store info from TTrees
  mithep::TEventInfo *info    = new mithep::TEventInfo();
  //mithep::TGenInfo *gen       = new mithep::TGenInfo();
  TClonesArray *muonArr = new TClonesArray("mithep::TMuon");
  TClonesArray *electronArr = new TClonesArray("mithep::TElectron");
  TClonesArray *pvArr         = new TClonesArray("mithep::TVertex");
  //EtaEtaMassData_t *eem = new EtaEtaMassData_t();
  
  //
  // Set up event dump to file
  //
  ofstream evtfile;
  char evtfname[100];    
  sprintf(evtfname,"%s/events.txt",outputDir.Data());
  printf("Checking that file can be written: %s\n",evtfname);
  evtfile.open(evtfname);
  assert(evtfile.is_open());
  
  //
  // loop over samples
  //
  for(UInt_t isam=0; isam<samplev.size(); isam++) {        
    if(isam==0 && !hasData) continue;
    
    /*  DONT USE THIS. REMOVE WHEN SURE IT IS NOT REQUIRED

    //
    // Set up output (eta,eta,mass) EEM file, if needed
    //
    TString outEEMName;
    TFile *eemFile=NULL;
    TTree *eemTree=NULL;
    if (generateEEMFile.size()) {
      outEEMName = ntupDir + TString("/") + snamev[isam] + TString("_") + TString(generateEEMFile.c_str()) + TString("_EtaEtaM.root");
      eemFile = new TFile(outEEMName,"RECREATE");
      eemTree = new TTree("Data","Data");
      assert(eemTree);
      eemTree->Branch("Data","EtaEtaMassData_t",&eem);
    }
    */
    //
    // Set up output ntuple file for the sample
    //
    TString outName = ntupDir + TString("/") + snamev[isam] + TString("_select.root");
    TFile *outFile = new TFile(outName,"RECREATE");
    TTree *outTree = new TTree("Events","Events");
    XemuData data;
    /*outTree->Branch("Events", &data.runNum, 
      "runNum/i:evtNum:lumiSec:nTracks0:nCaloTowers0:nPV:nJets:caloMEx/F:caloMEy:caloSumET:tcMEx:tcMEy:tcSumET:pfMEx:pfMEy:pfSumET:mass:pt:y:phi:pt_1:eta_1:phi_1:scEt_1:scEta_1:scPhi_1:hltMatchBits_1/i:q_1/I:pt_2/F:eta_2:phi_2:scEt_2:scEta_2:scPhi_2:hltMatchBits_2/i:q_2/I:weight/F");*/
    outTree->Branch("mass", &(data.mass), "data.mass/F");
    outTree->Branch("weight", &data.weight, "data.weight/F");
    outTree->Branch("pt_e", &(data.pt_e), "data.pt_e/F");   
    outTree->Branch("rapidity", &(data.rapidity), "data.rapidity/F");

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
	std::cout << "JSON file " << samp->jsonv[ifile] << "\n";
        jsonParser.Initialize(samp->jsonv[ifile].Data()); 
      }
      
      // Get the TTree
      eventTree = (TTree*)infile->Get("Events"); assert(eventTree);

      // Set branch address to structures that will store the info  
      eventTree->SetBranchAddress("Info",       &info);          TBranch *infoBr       = eventTree->GetBranch("Info");
      eventTree->SetBranchAddress("Electron" ,  &electronArr  ); TBranch *electronBr   = eventTree->GetBranch("Electron");
      eventTree->SetBranchAddress("Muon"     ,  &muonArr      ); TBranch *muonBr       = eventTree->GetBranch("Muon");
      eventTree->SetBranchAddress("PV",         &pvArr);         TBranch *pvBr         = eventTree->GetBranch("PV");
      // Generator information is present only for MC. Moreover, we
      // need to look it up only for signal MC in this script

      //REMOVE THIS *genBr STUFF BELOW
      /*
      TBranch *genBr = 0;
      if( snamev[isam] == "zee" ){
	eventTree->SetBranchAddress("Gen",&gen);
	genBr = eventTree->GetBranch("Gen");
      }
      */

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
        //if( snamev[isam] == "zee" )
	//genBr->GetEntry(ientry);

        if(hasJSON && !jsonParser.HasRunLumi(info->runNum, info->lumiSec)) continue;  // not certified run? Skip to next event...

	// Configure the object for trigger matching	
	//bool isData = (isam == 0 && hasData);

	//ADD MY CODE TO TRIGGERCONSTANT CLASS AND CALL AS BELOW
	/*	TriggerConstantSet constantsSet = Full2011DatasetTriggers; // Enum from TriggerSelection.hh
	TriggerSelection requiredTriggers(constantsSet, isData, info->runNum);
	ULong_t eventTriggerBit = requiredTriggers.getEventTriggerBit();
	ULong_t leadingTriggerObjectBit = requiredTriggers.getLeadingTriggerObjectBit();
	ULong_t trailingTriggerObjectBit = requiredTriggers.getTrailingTriggerObjectBit();*/

        /*ULong_t eventTriggerBit = kHLT_Mu17_Ele8_CaloIdL_MuObj | kHLT_Mu17_Ele8_CaloIdL_EGObj 
	  | kHLT_Mu8_Ele17_CaloIdL_MuObj | kHLT_Mu8_Ele17_CaloIdL_EGObj;*/

      ULong_t eventTriggerBit = (kHLT_Mu17_Ele8_CaloIdL | kHLT_Mu8_Ele17_CaloIdL | kHLT_Mu15_Photon20_CaloIdL | kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL);
      /*      ULong_t leadingTriggerObjectBit = kHLT_Mu17_Ele8_CaloIdL_MuObj | kHLT_Mu17_Ele8_CaloIdL_EGObj
	| kHLT_Mu8_Ele17_CaloIdL_MuObj | kHLT_Mu8_Ele17_CaloIdL_EGObj;
      ULong_t trailingTriggerObjectBit = kHLT_Mu17_Ele8_CaloIdL_MuObj | kHLT_Mu17_Ele8_CaloIdL_EGObj
      | kHLT_Mu8_Ele17_CaloIdL_MuObj | kHLT_Mu8_Ele17_CaloIdL_EGObj;*/
  
      ULong_t electronTriggerObjectBit = (kHLT_Mu17_Ele8_CaloIdL_EGObj | kHLT_Mu8_Ele17_CaloIdL_EGObj | kHLT_Mu15_Photon20_CaloIdL_EGObj | kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_EGObj);
      ULong_t muonTriggerObjectBit =  (kHLT_Mu17_Ele8_CaloIdL_MuObj | kHLT_Mu8_Ele17_CaloIdL_MuObj | kHLT_Mu15_Photon20_CaloIdL_MuObj | kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_MuObj);
      // Apply trigger cut at the event level        	
      if(!(info->triggerBits & eventTriggerBit)) continue;  // no trigger accept? Skip to next event...                                   

	/*cout << "trigger bit info is: ";
	bin(info->triggerBits,std::cout);
	std::cout << "\n";
        cout << "emu bits are: ";
	bin(eventTriggerBit,std::cout);
	std::cout << "\n";*/

        electronArr->Clear(); 
	electronBr->GetEntry(ientry);	
	muonArr->Clear(); 
	muonBr->GetEntry(ientry);
        // loop through electrons
        for(Int_t i=0; i<electronArr->GetEntriesFast(); i++) {
	  mithep::TElectron *electron = (mithep::TElectron*)((*electronArr)[i]);
	  // energy scale correction is done only for the data
	  if(isam==0)
	    //need to look at energy scale correction for single electron
	    //correctEnergyScale(electron, info);
	  
	  //
	  // Apply electron cuts
	  //
	  if(    electron->scEt < 20      )   continue;
	  if((fabs(electron->scEta)>kGAP_LOW) && (fabs(electron->scEta)<kGAP_HIGH)) continue;
	  if( fabs(electron->scEta) > 2.5 )   continue;  // outside eta range? Skip to next event...
	  if( !(electron->isEcalDriven)   )   continue;  // not ECAL seeded electrons? Skip to next event...
	  if( ! passSmurf(electron)        )   continue;  
	  if( !(electron->hltMatchBits & electronTriggerObjectBit)) continue;  // electron matched to HLT object? 
          //cout << "electron hlt match: " <<  electron->hltMatchBits << "\n";
	  if( fabs(electron->d0)>0.02 ) continue;
	  //if( fabs(electron->dz)>1.0 ) continue;//smurf cut is actually tighter <0.1
	  
	  // loop through muons
	  for(Int_t j=0; j<muonArr->GetEntriesFast(); j++) {
	    mithep::TMuon *muon = (mithep::TMuon*)((*muonArr)[j]);
            //cout << "muon hlt match: " <<  muon->hltMatchBits << "\n";

	    //
	    // Apply muon cuts
	    //
	    // Pt and eta
	    //if( muon->pt < 10 )              continue;
            if( muon->pt < 16 )              continue;
	    if( fabs(muon->eta > 2.4 ) )     continue;
	    // ID, quality, etc
	    if(!(muon->typeBits & kGlobal) ) continue;
	    if((muon->nTkHits < 11) )        continue;
	    if((fabs(muon->d0) > 0.2) )      continue;
	    if((muon->nPixHits < 1) )        continue;
	    if((muon->nSeg < 2) )            continue;
	    if((muon->nValidHits < 1) )      continue;
	    // not isolated? Skip to next event...
	    Double_t iso = (muon->trkIso03 + muon->emIso03 + muon->hadIso03)/(muon->pt);
	    if((iso > 0.15) )                   continue;
	    // no trigger matches? Skip to next event...
	    if(!(muon->hltMatchBits & muonTriggerObjectBit) ) continue; // Muon matched to triggerObject bit
	    //cout << "muon hlt match: " <<  muon->hltMatchBits << "\n";
	    if( fabs(muon->d0)>0.02 ) continue;
	    //if( fabs(muon->dz)>1.0 ) continue;

	    //
	    // Apply muon cuts
	    //
          
	    bool isOppositeSign = true;
	    if(electron->q == muon->q) {
	      isOppositeSign = false;
	      if(electron->q > 0) nPosSSv[isam] += weight;
	      else                nNegSSv[isam] += weight;
	    }

	    // Keep only opposite sign e-mu pairs
	    if( !isOppositeSign) continue;

	    /******** We have an e-mu candidate! HURRAY! ********/
            nsel    += weight;
	    nselvar += weight*weight;

	    // Calculate e-mu invariant mass

            TLorentzVector electron4v, muon4v, emu4v; 
            electron4v.SetPtEtaPhiM(electron->pt,electron->eta,electron->phi,0.000511);
            muon4v      .SetPtEtaPhiM(muon->pt    ,muon->eta    ,muon->phi    ,0.105658);
	    emu4v = electron4v + muon4v;
	    double mass = emu4v.M();
	  
	  
	    //DO NOT USE EVENT DUMP, UNLESS i MODIFY CODE TO TAKE EMU PAIRS. EASY TO DO
	    // event printout
	    /*if(isam==0)
	      eventDump(evtfile, dielectron, info->runNum, info->lumiSec, info->evtNum, 
	      leadingTriggerObjectBit, trailingTriggerObjectBit);*/
	  
	    //
	    // Fill histograms
	    // 
	    hMassv[isam]->Fill(mass,weight);
            hMass4v[isam]->Fill(mass,weight);
	    
	    pvArr->Clear();
	    pvBr->GetEntry(ientry);
	    UInt_t nGoodPV=0;
	    for(Int_t ipv=0; ipv<pvArr->GetEntriesFast(); ipv++) {
	      const mithep::TVertex *pv = (mithep::TVertex*)((*pvArr)[ipv]);                  
	      if(pv->nTracksFit                        < 1)  continue;
	      if(pv->ndof                              < 4)  continue;
	      if(fabs(pv->z)                           > 24) continue;
	      if(sqrt((pv->x)*(pv->x)+(pv->y)*(pv->y)) > 2)  continue;
	      nGoodPV++;
	    }
	    hNGoodPVv[isam]->Fill(nGoodPV,weight);
	  
	    // fill ntuple data
	    double weightSave = weight;
	    //eem->weight(weight);
	  
	    // Note: we do not need jet count at the moment. It can be found
	    // by looping over PFJets list if needed. See early 2011 analysis.
	    int njets = -1;
	    fillData(&data, info, electron, muon, emu4v, pvArr->GetEntriesFast(), njets, weightSave);
	    outTree->Fill();

	  
	    nsel    += weight;
	    nselvar += weight*weight;
	  }
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
  
  plotMass.SetName("masslog");
  plotMass.SetLogy();
  if( plotMass.GetStack() != NULL)
    plotMass.SetYRange((1e-7)*(plotMass.GetStack()->GetMaximum()),10.*(plotMass.GetStack()->GetMaximum()));  
  plotMass.Draw(c,kFALSE,format);

  //=======

  TCanvas *c2 = MakeCanvas("c2","c",canw,canh);

  // dielectron mass
  sprintf(ylabel,"Events / %.1f GeV/c^{2}",hMass4v[0]->GetBinWidth(1));
  CPlot plotMass2("mass","","m(e^{+}e^{-}) [GeV/c^{2}]",ylabel);
  if(hasData) { plotMass2.AddHist1D(hMass4v[0],samplev[0]->label,"E"); }
  printf("  debug1\n");fflush(stdout);
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    hMass4v[isam]->Scale(mcscale);
    plotMass2.AddToStack(hMass4v[isam],samplev[isam]->label,samplev[isam]->color);
  }
  printf("  debug2\n");fflush(stdout);
  if(samplev.size()>5)
    plotMass2.SetLegend(0.75,0.55,0.98,0.9);
  else
    plotMass2.TransLegend(0.1,0);
  if(lumi>0) plotMass2.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
  plotMass2.Draw(c2,kFALSE,format);
  
  plotMass2.SetName("masslog");
  plotMass2.SetLogy();
  if( plotMass2.GetStack() != NULL)
    plotMass2.SetYRange((1e-7)*(plotMass2.GetStack()->GetMaximum()),10.*(plotMass2.GetStack()->GetMaximum()));  
  plotMass2.Draw(c2,kFALSE,format);
    
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

  cout << endl;
  cout << " <> Output saved in " << outputDir << "/" << endl;
  cout << endl;
        
  gBenchmark->Show("selectEvents");       
} 


//=== FUNCTION DEFINITIONS ======================================================================================

//--------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------
void fillData(XemuData *data, const mithep::TEventInfo *info, const mithep::TElectron *electron,
	      const mithep::TMuon *muon, const TLorentzVector &emu4v,
              const UInt_t npv, const UInt_t njets, const Double_t weight)
{
  data->runNum         = info->runNum;
  data->evtNum         = info->evtNum;
  data->lumiSec        = info->lumiSec;
  data->nPV            = npv;
  data->nJets          = njets;                                        
  data->pfSumET        = info->pfSumET;
  data->mass           = emu4v.M();
  data->pt_e           = electron->pt;
  data->eta_e          = electron->eta;
  data->phi_e          = electron->phi;
  data->scEt_e         = electron->scEt;
  data->scEta_e        = electron->scEta;
  data->scPhi_e        = electron->scPhi;
  data->hltMatchBits_e = electron->hltMatchBits;
  data->q_e            = electron->q;
  data->pt_mu           = muon->pt;
  data->eta_mu          = muon->eta;
  data->phi_mu          = muon->phi;
  data->hltMatchBits_mu = muon->hltMatchBits;
  data->q_mu            = muon->q;
  data->weight         = weight;
  data->rapidity       = emu4v.Rapidity();
}

//--------------------------------------------------------------------------------------------------
void eventDump(ofstream &ofs, const mithep::TElectron *electron, 
               const UInt_t runNum, const UInt_t lumiSec, const UInt_t evtNum,
	       const Float_t &mass, 
	       const UInt_t triggerObj1, const UInt_t triggerObj2)
{
  ofs << endl;
  ofs << "Run:" << runNum;
  ofs << "  Lumi:" << lumiSec;
  ofs << "  Event:" << evtNum;
  ofs << "  mass: " << mass;
  ofs << "  pt: " << electron->pt << endl;
  
  ofs << "----------+-----------+-----------+-----------+-----------+-----------+-------------+------------+------------+-----------+------" << endl;
  ofs << "  SC ET   |  SC eta   |   SC phi  | trkiso/pt | emiso/pt  | hadiso/pt | sigiEtaiEta |    deta    |    dphi    |    H/E    | HLT" << endl;
  ofs << "----------+-----------+-----------+-----------+-----------+-----------+-------------+------------+------------+-----------+------" << endl;
      
  ofs << setw(9) << electron->scEt << " |";
  ofs << setw(10) << electron->scEta << " |";
  ofs << setw(10) << electron->scPhi << " |";
  ofs << setw(10) << electron->trkIso03/electron->pt << " |";
  ofs << setw(10) << electron->emIso03/electron->pt << " |";
  ofs << setw(10) << electron->hadIso03/electron->pt << " |";
  ofs << setw(12) << electron->sigiEtaiEta << " |";
  ofs << setw(12) << electron->deltaEtaIn << "|";
  ofs << setw(12) << electron->deltaPhiIn << "|";
  ofs << setw(10) << electron->HoverE << " |";
  if(electron->hltMatchBits & triggerObj1)
    ofs << " LEAD" << endl; 
  else if(electron->hltMatchBits & triggerObj2)
    ofs << " TRAIL" << endl;
  else
    ofs << " NOMAT" << endl;
    
  ofs << setw(9) << electron->scEt << " |";
  ofs << setw(10) << electron->scEta << " |";
  ofs << setw(10) << electron->scPhi << " |";
  ofs << setw(10) << electron->trkIso03/electron->pt << " |";
  ofs << setw(10) << electron->emIso03/electron->pt << " |";
  ofs << setw(10) << electron->hadIso03/electron->pt << " |";
  ofs << setw(12) << electron->sigiEtaiEta << " |";
  ofs << setw(12) << electron->deltaEtaIn << "|";
  ofs << setw(12) << electron->deltaPhiIn << "|";
  ofs << setw(10) << electron->HoverE << " |";
  if(electron->hltMatchBits & triggerObj1)
    ofs << " LEAD" << endl; 
  else if(electron->hltMatchBits & triggerObj2)
    ofs << " TRAIL" << endl;
  else
    ofs << " NOMAT" << endl;
}    

//===================================
//Code to print out numbers as binary
//===================================
template < typename T >

inline T highbit(T& t)
{
  return t = (((T)(-1)) >> 1) + 1;
}

//============================ 

template < typename T >

std::ostream& bin(T& value, std::ostream &o)
{
  for ( T bit = highbit(bit); bit; bit >>= 1 )
    {
      o << ( ( value & bit ) ? '1' : '0' );
    }
  return o;
}
