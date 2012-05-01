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
#include "../Include/TGenInfo.hh"
#include "../Include/TDielectron.hh"
#include "../Include/TVertex.hh"

// lumi section selection with JSON files
#include "../Include/JsonParser.hh"

// Helper functions for Electron ID selection
#include "../Include/EleIDCuts.hh" 

#include "../Include/ElectronEnergyScale.hh" //energy scale correction
#include "../Include/EtaEtaMass.hh" // EtaEtaMassData_t definition

#include "../Include/EventSelector.hh"
#include "../Include/FEWZ.hh"

#endif

// define structure for output ntuple
#include "../Include/ZeeData.hh"


//=== FUNCTION DECLARATIONS ======================================================================================

// fill ntuple of selected events
#ifdef ZeeData_is_TObject
// Function is defined in fillData
//void fillData(ZeeData_t *data, const mithep::TEventInfo *info, const mithep::TDielectron *dielectron,
//              const UInt_t npv, const UInt_t nGoodPV, const UInt_t njets, const Double_t weight);
#else
void fillData(ZeeData *data, const mithep::TEventInfo *info, const mithep::TDielectron *dielectron,
              const UInt_t npv, const UInt_t nGoodPV, const UInt_t njets, const Double_t weight);
#endif

// print event dump
void eventDump(ofstream &ofs, const mithep::TDielectron *dielectron, 
               const UInt_t runNum, const UInt_t lumiSec, const UInt_t evtNum, 
	       const UInt_t triggerObj1, const UInt_t triggerObj2);


//=== MAIN MACRO =================================================================================================

void selectEvents(const TString conf, const TString triggerSetString="Full2011DatasetTriggers", DYTools::TSystematicsStudy_t runMode=DYTools::NORMAL, int debugMode=0) 
{  
  gBenchmark->Start("selectEvents");

  // fast check
  TriggerConstantSet triggerSet=DetermineTriggerSet(triggerSetString);  
  assert ( triggerSet != TrigSet_UNDEFINED );

  // Construct the trigger object
  TriggerSelection requiredTriggers(triggerSetString, true, 0);
  assert(requiredTriggers.isDefined());

  if (debugMode) std::cout << "\n\n\tDEBUG MODE is ON\n\n";

  std::cout << "\n\nRun mode: " << SystematicsStudyName(runMode) << "\n";
  switch(runMode) {
  case DYTools::NORMAL:
  case DYTools::ESCALE_STUDY:
    break;
  default:
    std::cout << "selectEvents is not ready for runMode=" << SystematicsStudyName(runMode) << "\n";
    throw 2;
  }
  
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
  if (runMode==DYTools::ESCALE_STUDY) { 
    std::cout << "Since runMode=EScale_Study, changing escaleTag=<" << escaleTag << "> to 'UNCORRECTED'\n";
    escaleTag="UNCORRECTED";
  }
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

  const TString ntupDir = outputDir + TString("/ntuples"); gSystem->mkdir(ntupDir,kTRUE);
  
  Bool_t hasData = (samplev[0]->fnamev.size()>0);

  //const Double_t kGAP_LOW  = 1.4442;
  //const Double_t kGAP_HIGH = 1.566;
    
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
    sprintf(hname,"hNGoodPV_%s",snamev[isam].Data());       hNGoodPVv.push_back(new TH1F(hname,"",35,-0.5,34.5));        hNGoodPVv[isam]->Sumw2();
    
    nSelv.push_back(0);
    nSelVarv.push_back(0);
    nPosSSv.push_back(0);
    nNegSSv.push_back(0);    
  }
  
  // 
  // Read weights from a file
  //
  const bool useFewzWeights = true;
  const bool cutZPT100 = true;
  FEWZ_t fewz(useFewzWeights,cutZPT100);
  if (useFewzWeights && !fewz.isInitialized()) {
    std::cout << "failed to prepare FEWZ correction\n";
    throw 2;
  }
  const int new_fewz_code=1;
  TH2D *weights[DYTools::nMassBins];   // temporary
  TH2D *weightErrors[DYTools::nMassBins]; // temporary
  if (!new_fewz_code) {
  if(cutZPT100)
    cout << "NOTE: in MC, for Z/gamma* PT>100 the FEWZ weights for 80<PT<100 GeV are used!" << endl;
  TFile fweights("../root_files/fewz/weights_stepwise_prec10-5_fine12.root");
  if( !fweights.IsOpen() ) assert(0);
  for(int i=0; i<DYTools::nMassBins; i++){
    TString hnames = TString::Format("weight_%02d",i+1);
    weights[i] = (TH2D*)fweights.Get(hnames);
    hnames = TString::Format("h_weighterror_%02d",i+1);
    weightErrors[i] = (TH2D*)fweights.Get(hnames);
  }
  }


  //
  // Access samples and fill histograms
  //  
  TFile *infile=0;
  TTree *eventTree=0;  
  
  // Data structures to store info from TTrees
  mithep::TEventInfo *info    = new mithep::TEventInfo();
  mithep::TGenInfo *gen       = new mithep::TGenInfo();
  TClonesArray *dielectronArr = new TClonesArray("mithep::TDielectron");
  TClonesArray *pvArr         = new TClonesArray("mithep::TVertex");
  EtaEtaMassData_t *eem = new EtaEtaMassData_t();
  
#ifdef ZeeData_is_TObject
  ZeeData_t::Class()->IgnoreTObjectStreamer();
#endif

  //
  // Set up event dump to file
  //
  const int dump_events=0;
  ofstream evtfile;
  char evtfname[100];
  sprintf(evtfname,"%s/events.txt",outputDir.Data());
  if (dump_events) {
    evtfile.open(evtfname);
    assert(evtfile.is_open());
  }
  
  //
  // loop over samples
  //
  for(UInt_t isam=0; isam<samplev.size(); isam++) {        
    if(isam==0 && !hasData) continue;
    
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

    // Define dielectron selector
    DielectronSelector_t eeSelector(DielectronSelector_t::_selectDefault,
				    &escale);

    //
    // Set up output ntuple file for the sample
    //
    TString outName = ntupDir + TString("/") + snamev[isam] + analysisTag_USER + TString("_select.root");
    if (runMode==DYTools::ESCALE_STUDY) {
      if (isam==0) {
	outName.Replace(outName.Index("_select."),sizeof("_select."),"_select_uncorrected.");
      }
      else {
	std::cout << "... runMode=EScale_Study, skipping the non-data files\n";
	break;
      }
    }
    TFile *outFile = new TFile(outName,"RECREATE");
    TTree *outTree = new TTree("Events","Events");
#ifdef ZeeData_is_TObject
    ZeeData_t *data=new ZeeData_t();
    outTree->Branch("Events","ZeeData_t",&data);
#else
    ZeeData *data=new ZeeData();
    outTree->Branch("Events", &data.runNum, 
    "runNum/i:evtNum:lumiSec:nTracks0:nCaloTowers0:nPV:nGoodPV:nJets:caloMEx/F:caloMEy:caloSumET:tcMEx:tcMEy:tcSumET:pfMEx:pfMEy:pfSumET:mass:pt:y:phi:pt_1:eta_1:phi_1:scEt_1:scEta_1:scPhi_1:hltMatchBits_1/i:q_1/I:pt_2/F:eta_2:phi_2:scEt_2:scEta_2:scPhi_2:hltMatchBits_2/i:q_2/I:weight/F"
    //"runNum/i:evtNum:lumiSec:nTracks0:nCaloTowers0:nPV:nJets:caloMEx/F:caloMEy:caloSumET:tcMEx:tcMEy:tcSumET:pfMEx:pfMEy:pfSumET:mass:pt:y:phi:pt_1:eta_1:phi_1:scEt_1:scEta_1:scPhi_1:hltMatchBits_1/i:q_1/I:pt_2/F:eta_2:phi_2:scEt_2:scEta_2:scPhi_2:hltMatchBits_2/i:q_2/I:weight/F"
    );
#endif

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
      eventTree->SetBranchAddress("Dielectron", &dielectronArr); TBranch *dielectronBr = eventTree->GetBranch("Dielectron");
      eventTree->SetBranchAddress("PV",         &pvArr);         TBranch *pvBr         = eventTree->GetBranch("PV");
      // Generator information is present only for MC. Moreover, we
      // need to look it up only for signal MC in this script
      TBranch *genBr = 0;
      if( snamev[isam] == "zee" ){
	eventTree->SetBranchAddress("Gen",&gen);
	genBr = eventTree->GetBranch("Gen");
      }
      
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
	  // if this is a spec.skim file, rescale xsec
	  const int new_adjustment_code=1;
	  if (new_adjustment_code) {
	  AdjustXSectionForSkim(infile,xsec,eventTree->GetEntries(),1);
	  }
	  else {
	  TTree *descrTree=(TTree*)infile->Get("Description");
	  if (descrTree) {
	    UInt_t origNumEntries=0;
	    descrTree->SetBranchAddress("origNumEntries",&origNumEntries);
	    descrTree->GetEntry(0);
	    if (origNumEntries>0) {
	      Double_t factor=eventTree->GetEntries()/double(origNumEntries);
	      std::cout << " -> rescaling xsec by " << factor << " due to skimming\n";
	      xsec*=factor;
	    }
	    delete descrTree;
	  }
	  else std::cout << "descrTree not found\n";
	  // proceed
	  }
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
      std::cout << "numEntries = " << eventTree->GetEntries() << std::endl;
      for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
	if (debugMode && (ientry>10000)) break; // debug option
	if(ientry >= maxEvents) break;
	
	infoBr->GetEntry(ientry);
      if( snamev[isam] == "zee" )
	  genBr->GetEntry(ientry);

	// Load FEWZ weights for signal MC
	double fewz_weight = 1.0;
	if(( snamev[isam] == "zee" ) && useFewzWeights) {
	  if (new_fewz_code) {
	    fewz_weight=fewz.getWeight(gen->vmass,gen->vpt,gen->vy);
	  }
	  else {
	  int ibinPreFsr = DYTools::findMassBin(gen->vmass);
	  // If mass is larger than the highest bin boundary
	  // (last bin), use the last bin.
	  if(ibinPreFsr == -1 && gen->vmass >= massBinLimits[nMassBins] )
	    ibinPreFsr = nMassBins-1;
	  // Find FEWZ-powheg reweighting factor 
	  // that depends on pre-FSR Z/gamma* rapidity, pt, and mass
	  if(useFewzWeights){
	    if(ibinPreFsr != -1 && ibinPreFsr < DYTools::nMassBins){
	      int ptBin = weights[ibinPreFsr]->GetXaxis()->FindBin( gen->vpt );
	      int yBin = weights[ibinPreFsr]->GetYaxis()->FindBin( gen->vy );
	      // In case if pt or y are outside of the weight maps,
	      // set them to the closest bin.
	      if(ptBin == weights[ibinPreFsr]->GetNbinsX() + 1)
		ptBin = weights[ibinPreFsr]->GetNbinsX();
	      if(ptBin == 0)
		ptBin = 1;
	      if(yBin == weights[ibinPreFsr]->GetNbinsY() + 1)
		yBin = weights[ibinPreFsr]->GetNbinsY();
	      if(yBin == 0)
		yBin = 1;
	      // Apply PT cut if needed
	      if( cutZPT100 ) 
		if( ptBin == weights[ibinPreFsr]->GetNbinsX() )
		  ptBin = weights[ibinPreFsr]->GetNbinsX() - 1;
	      fewz_weight = weights[ibinPreFsr]->GetBinContent( ptBin, yBin);
	    }else{
	      // Error printout is commented out: the maps now go down
	      // to 15 GeV. Events with generator level mass below 15 GeV 
	      // may contribute to 
	      // reconstructed events when reconstructed mass is above 15 GeV
// 	      cout << "Error: vmass outside of FEWZ weight maps, vmass=" 
// 		   << gen->vmass << endl;
	    }
	  }
	  }
	}


	// mithep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);      
	// if(hasJSON && !rlrm.HasRunLumi(rl)) continue;  // not certified run? Skip to next event...
        if(hasJSON && !jsonParser.HasRunLumi(info->runNum, info->lumiSec)) continue;  // not certified run? Skip to next event...

	// Configure the object for trigger matching	
	bool isData = (isam == 0 && hasData);
	//TriggerConstantSet constantsSet = Full2011DatasetTriggers; // Enum from TriggerSelection.hh
	//TriggerSelection requiredTriggers(constantsSet, isData, info->runNum);
	requiredTriggers.actOnData(isData);
	ULong_t eventTriggerBit = requiredTriggers.getEventTriggerBit(info->runNum);
	ULong_t leadingTriggerObjectBit = requiredTriggers.getLeadingTriggerObjectBit(info->runNum);
	ULong_t trailingTriggerObjectBit = requiredTriggers.getTrailingTriggerObjectBit(info->runNum);
	// Apply trigger cut at the event level        	
        if(!(info->triggerBits & eventTriggerBit)) continue;  // no trigger accept? Skip to next event...                                   

	dielectronArr->Clear(); 
	dielectronBr->GetEntry(ientry);	
	// loop through dielectrons
        for(Int_t i=0; i<dielectronArr->GetEntriesFast(); i++) {
	  mithep::TDielectron *dielectron = (mithep::TDielectron*)((*dielectronArr)[i]);

	  const int ee_selection_new_code=1;
	  if (ee_selection_new_code) {
	    // Keep the EEM values before any changes
	    eem->Assign(dielectron->scEta_1,dielectron->scEta_2,dielectron->mass);
	    
	    // Determine correction type
	    DielectronSelector_t::TEScaleCorrection_t escaleCorrType=
	      DielectronSelector_t::_escaleNone;
	    if (isam==0) {
	      switch (runMode) {
	      case DYTools::NORMAL: 
		escaleCorrType=DielectronSelector_t::_escaleData;
		break;
	      case DYTools::ESCALE_STUDY:
		escaleCorrType=DielectronSelector_t::_escaleUncorrected;
		break;
	      default:
		std::cout << "does not know what escale to apply for runMode=" << SystematicsStudyName(runMode) << "\n";
		throw 2;
	      }
	    }
	    if (!eeSelector(dielectron,
			    escaleCorrType,
			    leadingTriggerObjectBit,
			    trailingTriggerObjectBit)) continue;
	  
	    hMass2v[isam]->Fill(dielectron->mass,weight);
	    hMass3v[isam]->Fill(dielectron->mass,weight);
	  }
	  else {

          // Exclude ECAL gap region and cut out of acceptance electrons
          if((fabs(dielectron->scEta_1)>kECAL_GAP_LOW) && (fabs(dielectron->scEta_1)<kECAL_GAP_HIGH)) continue;
          if((fabs(dielectron->scEta_2)>kECAL_GAP_LOW) && (fabs(dielectron->scEta_2)<kECAL_GAP_HIGH)) continue;
          if((fabs(dielectron->scEta_1) > 2.5)       || (fabs(dielectron->scEta_2) > 2.5))       continue;  // outside eta range? Skip to next event...

	  // Keep the EEM values before any changes
	  eem->Assign(dielectron->scEta_1,dielectron->scEta_2,dielectron->mass);

	  //
	  // Energy scale corrections for data
	  // NOTE: the electrons and dielectron 4-vectors are updated, the supercluster quantities are not
	  //
          Double_t scEt1 = dielectron->scEt_1;
	  Double_t scEt2 = dielectron->scEt_2;
	  // Electron energy scale correction
          if(isam==0) {            	    
	    double corr1 = escale.getEnergyScaleCorrection(dielectron->scEta_1);
	    double corr2 = escale.getEnergyScaleCorrection(dielectron->scEta_2);
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
	  
          // loose mass window 
          if( dielectron->mass < 10 ) continue;
	  } // ee_selection_new_code
	  
	  
          /******** We have a Z candidate! HURRAY! ********/
	  
          if(dielectron->q_1 == dielectron->q_2) {
	    if(dielectron->q_1 > 0) nPosSSv[isam] += weight;
	    else                    nNegSSv[isam] += weight;
	  }
	  
	  // event printout
          if((isam==0) && evtfile.is_open())
            eventDump(evtfile, dielectron, info->runNum, info->lumiSec, info->evtNum, 
		      leadingTriggerObjectBit, trailingTriggerObjectBit);
	  
	  //
	  // Fill histograms
	  // 
	  hMassv[isam]->Fill(dielectron->mass,weight);
	  
	  pvArr->Clear();
          pvBr->GetEntry(ientry);
          UInt_t nGoodPV=0;
	  const int new_good_pv_code=1;
	  if (new_good_pv_code) {
	    nGoodPV=countGoodVertices(pvArr);
	  }
	  else {
          for(Int_t ipv=0; ipv<pvArr->GetEntriesFast(); ipv++) {
            const mithep::TVertex *pv = (mithep::TVertex*)((*pvArr)[ipv]);                  
            if(pv->nTracksFit                        < 1)  continue;
            if(pv->ndof                              < 4)  continue;
            if(fabs(pv->z)                           > 24) continue;
            if(sqrt((pv->x)*(pv->x)+(pv->y)*(pv->y)) > 2)  continue;
            nGoodPV++;
          }
	  }
          hNGoodPVv[isam]->Fill(nGoodPV,weight);
	  
	  // fill ntuple data
	  double weightSave = weight;
	  eem->weight(weight);
	  if( snamev[isam] == "zee" ) {
	    weightSave *= fewz_weight;
	    if (generateEEMFEWZFile) eem->weight(weightSave);
	  }
	  // Note: we do not need jet count at the moment. It can be found
	  // by looping over PFJets list if needed. See early 2011 analysis.
	  int njets = -1;
	  fillData(data, info, dielectron, pvArr->GetEntriesFast(), nGoodPV, njets, weightSave);
	  outTree->Fill();
	  if (eemTree) {
	    eemTree->Fill();
	    //std::cout << "store eem=" << (*eem) << "\n"; 
	  }
	  
	  nsel    += weight;
	  nselvar += weight*weight;
      
        }	 
      }
      cout << nsel << " +/- " << sqrt(nselvar) << " events" << endl;
      nSelv[isam]    += nsel;
      nSelVarv[isam] += nselvar;
      delete infile;
      infile=0, eventTree=0;
    }
    std::cout << "next file" << std::endl;
    outFile->Write();
    delete outTree;
    outFile->Close();        
    delete outFile;
    if (eemFile) {
      eemFile->Write();
      delete eemTree;
      eemFile->Close();
      delete eemFile;
    }

    eeSelector.printCounts(std::cout);
  }
  delete info;
  delete dielectronArr;
  delete pvArr;
  if (evtfile.is_open()) evtfile.close();


  // Write useful histograms
  // npv.root
  TString outNamePV = outputDir + TString("/npv") + analysisTag_USER + TString(".root");
  TFile *outFilePV = new TFile(outNamePV,"RECREATE");
  outFilePV->cd();
  for(UInt_t isam=0; isam<samplev.size(); isam++) {
    hNGoodPVv[isam]->Write();
  }
  outFilePV->Close();

  //--------------------------------------------------------------------------------------------------------------
  // Make plots
  //==============================================================================================================

  TString outFNameHistos = outputDir + TString("/selectEvents") + analysisTag_USER + TString("-plots.root");
  TFile *outFileHistos = new TFile(outFNameHistos,"RECREATE");
  outFileHistos->cd();
  for(UInt_t isam=0; isam<samplev.size(); isam++) hNGoodPVv[isam]->Write();

  TString canvasName="selectEvents" + study2Dstr;
  TCanvas *c = MakeCanvas(canvasName,canvasName,canw,canh);

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
    plotMass.SetYRange((1e-4)*(plotMass.GetStack()->GetMaximum()),10.*(plotMass.GetStack()->GetMaximum()));  
  plotMass.Draw(c,kFALSE,format);

  c->Write();
  outFileHistos->Close();

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
  cout << "file <" << txtfname << "> created\n";

  cout << endl;
  cout << " <> Output saved in " << outputDir << "/" << endl;
  cout << endl;
        
  gBenchmark->Show("selectEvents");       
} 


//=== FUNCTION DEFINITIONS ======================================================================================

//--------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------
#ifdef ZeeData_is_TObject
  // the function is defined in ZeeData.hh
#else
void fillData(ZeeData *data, const mithep::TEventInfo *info, const mithep::TDielectron *dielectron, 
              const UInt_t npv, const UInt_t nGoodPV, const UInt_t njets, const Double_t weight)
{
  data->runNum         = info->runNum;
  data->evtNum         = info->evtNum;
  data->lumiSec        = info->lumiSec;
  data->nPV            = npv;
  data->nGoodPV        = nGoodPV;
  data->nJets          = njets;                                        
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
#endif

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

