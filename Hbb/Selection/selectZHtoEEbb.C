//================================================================================================
//
// ZH->e e b b selection macro
//
//  * plots distributions associated with selected events
//  * prints list of selected events from data
//  * outputs ROOT files of events passing selection for each sample, 
//    which can be processed later by another script
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
#include "../Include/CSample.hh"        // helper class for organizing input ntuple files

// define structures to read in ntuple
#include "../Include/HiggsAnaDefs.hh"
#include "../Include/TEventInfo.hh"
#include "../Include/TElectron.hh"
#include "../Include/TJet.hh"
#include "../Include/TVertex.hh"
#include "../Include/MyTools.hh"

// lumi section selection with JSON files, stand alone from MIT framework
#include "../Include/JsonParser.hh"

// Helper functions for Electron ID selection
#include "../Include/EleIDCuts.hh" 

// define structure for output ntuple
#include "../Include/ZHtoEEbbData.hh"
#endif

//=== FUNCTION DECLARATIONS ======================================================================================

// fill ntuple of selected events
void fillData(ZHtoEEbbData *data, const mithep::TEventInfo *info, TLorentzVector ZMomentum,
	      const mithep::TElectron *electron1, const mithep::TElectron *electron2,
	      const mithep::TJet *jet1, const mithep::TJet *jet2, double dijetMass, double dijetPt,
              const UInt_t npv, const UInt_t njets, const Double_t weight);

// Electron cut values
const double cutEleETMin  = 20.0;
const double cutEleEtaMax = 2.5;
const double cutEleIsoMax = 0.1;

// Z cut values
const double cutZMassMin  = 75.0;
const double cutZMassMax  = 105.0;
const double cutZPtMin    = 50.0;

// Jet cut values
const double cutJetPTMin  = 20.0;
const double cutJetEtaMax = 2.5;
const double cutJetTrackCountMin = 2;
const double cutEMFractionMin = 0.01;
const double cutHadFractionMin = 0.01;
const double csv1cut = 0.85;
const double csv2cut = 0.55;

// Higgs (di-b-jet) cut values)
// Mass cuts for 115 GeV Higgs
const double cutHMassMin =  95.0;
const double cutHMassMax = 125.0;
const double cutHPtMin    = 50.0;

// Other
const double cutDRJetLeptonMin = 0.3;
const double cutDPhiZHMin = 2.95;
const int cutAdditionalCentralJetsMax = 1;

//=== MAIN MACRO =================================================================================================

void selectZHtoEEbb(const TString conf) 
{  
  gBenchmark->Start("selectZHtoEEbb");

  
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

  CPlot::sOutDir        = outputDir + TString("/plots");   gSystem->mkdir(CPlot::sOutDir,kTRUE);
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
  vector<TH1F*> hZMassv;
  vector<TH1F*> hHMassv;
  vector<TH1F*> hZPtv;   
  vector<TH1F*> hHPtv;   
  
  vector<TH1F*> hNGoodPVv;
  
  vector<Double_t> nSelv, nSelVarv;  
  
  UInt_t nProcessedEvents=0;

  char hname[100];
  for(UInt_t isam=0; isam<samplev.size(); isam++) {
    sprintf(hname,"hZMass_%i",isam);   hZMassv.push_back(new TH1F(hname,"",30,60,120));   hZMassv[isam]->Sumw2();
    sprintf(hname,"hHMass_%i",isam);   hHMassv.push_back(new TH1F(hname,"",50,50,150));   hHMassv[isam]->Sumw2();
    sprintf(hname,"hZPt_%i",isam);     hZPtv.push_back(new TH1F(hname,"",50,0,200));       hZPtv[isam]->Sumw2();
    sprintf(hname,"hHPt_%i",isam);     hHPtv.push_back(new TH1F(hname,"",50,0,300));       hHPtv[isam]->Sumw2();
    sprintf(hname,"hNGoodPV_%s",snamev[isam].Data());       hNGoodPVv.push_back(new TH1F(hname,"",15,-0.5,14.5));        hNGoodPVv[isam]->Sumw2();
    
    nSelv.push_back(0);
    nSelVarv.push_back(0);
  }
  
  //
  // Access samples and fill histograms
  //  
  TFile *infile=0;
  TTree *eventTree=0;  
  
  // Data structures to store info from TTrees
  mithep::TEventInfo *info    = new mithep::TEventInfo();
  TClonesArray *pfJetArr      = new TClonesArray("mithep::TJet");
  TClonesArray *electronArr   = new TClonesArray("mithep::TElectron");
  TClonesArray *pvArr         = new TClonesArray("mithep::TVertex");
  
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
    ZHtoEEbbData data;
    outTree->Branch("Events",&data.runNum,
		    "runNum/i:evtNum:lumiSec:nPV:nJets:caloMEx/F:caloMEy:caloSumET:tcMEx:tcMEy:tcSumET:pfMEx:pfMEy:pfSumET:diEMass:diEPt:diEY:diEPhi:pt_1:eta_1:phi_1:scEt_1:scEta_1:scPhi_1:hltMatchBits_1/i:q_1/I:pt_2/F:eta_2:phi_2:scEt_2:scEta_2:scPhi_2:hltMatchBits_2/i:q_2/I:diJetMass/F:diJetPt:jPt_1:jEta_1:jPhi_1:jMass_1:tche_1:tchp_1:jCSV_1:mcFlavor_1/I:jHltMatchBits_1:jPt_2/F:jEta_2:jPhi_2:jMass_2:tche_2:tchp_2:jCSV_2:mcFlavor_2/I:jHltMatchBits_2:weight/F");    

    //
    // loop through files
    //
    CSample* samp = samplev[isam];
    const UInt_t nfiles = samplev[isam]->fnamev.size();    
    for(UInt_t ifile=0; ifile<nfiles; ifile++) {
      cout << "Processing " << samp->fnamev[ifile] << "... "; cout.flush();
      infile = new TFile(samp->fnamev[ifile]);
      assert(infile);
    
      // Note: commented out way is the one to use with MitXxx packages
      // in CMSSW, the JsonParser is for decoupled from CMSSW use
      Bool_t hasJSON = kFALSE;
//       mithep::RunLumiRangeMap rlrm;
      JsonParser jsonParser;
      if((samp->jsonv.size()>0) && (samp->jsonv[ifile].CompareTo("NONE")!=0)) { 
        hasJSON = kTRUE;
//         rlrm.AddJSONFile(samp->jsonv[ifile].Data()); 
        jsonParser.Initialize(samp->jsonv[ifile].Data()); 
      }
      
      // Get the TTree
      eventTree = (TTree*)infile->Get("Events"); assert(eventTree);

      // Set branch address to structures that will store the info  
      eventTree->SetBranchAddress("Info",       &info);          TBranch *infoBr       = eventTree->GetBranch("Info");
      eventTree->SetBranchAddress("Electron",   &electronArr);   TBranch *electronBr   = eventTree->GetBranch("Electron");
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
//       for(UInt_t ientry=0; ientry<100000; ientry++) {  // for debugging
	if(ientry >= maxEvents) break;
	
	infoBr->GetEntry(ientry);

	// As above: commented out is the method to use within CMSSW,
	// and JsonParser is a stand-alone tool
//         mithep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);      
//         if(hasJSON && !rlrm.HasRunLumi(rl)) continue;  // not certified run? Skip to next event...
        if(hasJSON && !jsonParser.HasRunLumi(info->runNum, info->lumiSec)) continue;  // not certified run? Skip to next event...
        
	// The double electron trigger bit definitions
	UInt_t eventTriggerBit          = kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL;
	UInt_t leadingTriggerObjectBit  = kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele1Obj;
	UInt_t trailingTriggerObjectBit = kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele2Obj;
	if(isam==0) {
	  nProcessedEvents++;
 	}
	
        if(!(info->triggerBits & eventTriggerBit)) continue;  // no trigger accept? Skip to next event...                                   

	// Count number of good primary vertices
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

	electronArr->Clear(); 
	electronBr->GetEntry(ientry);	

	// 
	// Loop over electrons
	//
        for(Int_t i=0; i<electronArr->GetEntriesFast(); i++) {
	  mithep::TElectron *electron1 = (mithep::TElectron*)((*electronArr)[i]);

	  for(Int_t j=0; j<electronArr->GetEntriesFast(); j++) {
	    if( i==j ) continue; 
	    mithep::TElectron *electron2 = (mithep::TElectron*)((*electronArr)[j]);

	    // Exclude ECAL gap region and cut out of acceptance electrons
	    if((fabs(electron1->scEta)>kGAP_LOW) && (fabs(electron1->scEta)<kGAP_HIGH)) continue;
	    if((fabs(electron2->scEta)>kGAP_LOW) && (fabs(electron2->scEta)<kGAP_HIGH)) continue;
	    if((fabs(electron1->scEta) > cutEleEtaMax) || (fabs(electron2->scEta) > cutEleEtaMax))       continue;  // outside eta range? Skip to next event...
	    //
	    // ET thresholds for electrons
	    if( ! (electron1->scEt > cutEleETMin && electron2->scEt > cutEleETMin) ) continue;
	    // Both electrons must match trigger objects. At least one ordering
	    // must match
	    if( ! ( 
		   (electron1->hltMatchBits & leadingTriggerObjectBit && 
		    electron2->hltMatchBits & trailingTriggerObjectBit )
		   ||
		   (electron1->hltMatchBits & trailingTriggerObjectBit && 
		    electron2->hltMatchBits & leadingTriggerObjectBit ) ) ) continue;
	    
	    // Other cuts to both electrons
	    
	    // Apply electron ID to both electrons
	    bool useSmurf = false;
	    if(useSmurf){
	      // The Smurf electron ID package as used in HWW analysis.
	      // It contains cuts like VBTF WP80 for pt>20, VBTF WP70 for pt<10
	      // with some customization, plus impact parameter cuts dz and dxy
	      if(! ( passSmurf(electron1) && passSmurf(electron2)) ) continue;  
	    }else{
	      // Use ID cuts as in AN 11-240.
	      // The note says "WP95", which is interpreted here is 2011 WP95
	      // ID cuts listed on https://twiki.cern.ch/twiki/bin/view/CMS/VbtfEleID2011
	      // NOTE: the fBrem cut is not applied, IP cuts are not applied,
	      // conversion rejection is not applied. Not clear if these
	      // are applied in 11-240.
	      if( ! ( passWP95ID2011(electron1) && passWP95ID2011(electron2)) ) continue;
	      // Cut on relative combined isolation of each electron
	      // as in 11-240
	      double iso1 = (electron1->trkIso03 + electron1->emIso03 + electron1->hadIso03)/electron1->pt;
	      double iso2 = (electron2->trkIso03 + electron2->emIso03 + electron2->hadIso03)/electron2->pt;
	      if( !( iso1 < cutEleIsoMax && iso2 < cutEleIsoMax ) ) continue;
	    }
	    
	    TLorentzVector lep1Momentum, lep2Momentum;
	    lep1Momentum.SetPtEtaPhiM(electron1->pt, electron1->eta, electron1->phi, 0.000511);
	    lep2Momentum.SetPtEtaPhiM(electron2->pt, electron2->eta, electron2->phi, 0.000511);
	    TLorentzVector ZMomentum = lep1Momentum + lep2Momentum;

	    // mass window for Z candidate
	    if((ZMomentum.M() < cutZMassMin) || (ZMomentum.M() > cutZMassMax)) continue;
	    
	    
	    // Selection of Z candidate is complete. 
	    //
	    // 	  printf("DEBUG: Z candidate found\n");
	    
	    // 
	    // Find a pair of candidate b-jets with highest combined
	    // probability.
	    //
	    double bestSignificance = -1.0;
	    const mithep::TJet *bjet1 = 0;
	    const mithep::TJet *bjet2 = 0;
	    int totalCentralJets = 0;
	    pfJetArr->Clear();
	    pfJetBr->GetEntry(ientry);
	    for(Int_t ijet=0; ijet<pfJetArr->GetEntriesFast(); ijet++) {
	      const mithep::TJet *jet1 = (mithep::TJet*)((*pfJetArr)[ijet]);
	      for(Int_t jjet=0; jjet<pfJetArr->GetEntriesFast(); jjet++) {
		if( ijet == jjet ) continue;
		const mithep::TJet *jet2 = (mithep::TJet*)((*pfJetArr)[jjet]);
		
		// Cuts on jet kinematics
		if( !( fabs(jet1->eta) < cutJetEtaMax && fabs(jet2->eta) < cutJetEtaMax) ) continue;
		// 	      printf("DEBUG: pair of central jets found\n");
		// Jet Pt cuts (JEC corrections are already applied before)
		// NOT CLEAR: do leptons from Z need to be excluded from jets
		// in addition to min separation cut below? Description in 11-240 not totally clear.
		if( !( jet1->pt > cutJetPTMin && jet2->pt > cutJetPTMin ) ) continue;
		// 	      printf("DEBUG: Pt is high enough\n");
		
		// Cuts on charged particle multiplicity and EM, Had energy fractions
		if( !(jet1->nCharged >= cutJetTrackCountMin && jet2->nCharged >= cutJetTrackCountMin ) ) continue;
		// Not clear in 11-240 if charged and neutral are cut on
		// separately or added
		if( !( (jet1->chgEMfrac + jet1->neuEMfrac) > cutEMFractionMin 
		       && (jet2->chgEMfrac + jet2->neuEMfrac) > cutEMFractionMin ) ) continue;
		if( !( (jet1->chgHadrfrac + jet1->neuHadrfrac) > cutHadFractionMin 
		       && (jet2->chgHadrfrac + jet2->neuHadrfrac) > cutHadFractionMin ) ) continue;
		
		// Make sure leptons from Z are not too close to jets
		if( toolbox::deltaR(jet1->eta, jet1->phi, electron1->scEta, electron1->scPhi) < cutDRJetLeptonMin ) continue;
		if( toolbox::deltaR(jet1->eta, jet1->phi, electron2->scEta, electron2->scPhi) < cutDRJetLeptonMin ) continue;
		if( toolbox::deltaR(jet2->eta, jet2->phi, electron1->scEta, electron1->scPhi) < cutDRJetLeptonMin ) continue;
		if( toolbox::deltaR(jet2->eta, jet2->phi, electron2->scEta, electron2->scPhi) < cutDRJetLeptonMin ) continue;
		// 	      printf("DEBUG: Not matched to leptons\n");

		double thisSignificance = jet1->csv + jet2->csv;
		if(thisSignificance > bestSignificance){
		  bjet1 = jet1;
		  bjet2 = jet2;
		  bestSignificance = thisSignificance;
		}
	      } // end inner loop over jets
	      
	      // Count jets toward jet veto
	      if( jet1->pt > cutJetPTMin && fabs(jet1->eta) < cutJetEtaMax
		  && toolbox::deltaR(jet1->eta, jet1->phi, electron1->scEta, electron1->scPhi) > cutDRJetLeptonMin
		  && toolbox::deltaR(jet1->eta, jet1->phi, electron2->scEta, electron2->scPhi) > cutDRJetLeptonMin 
		  ){
		totalCentralJets++;
	      } 
	      
	    } // end outer loop over jets
	    // 	  printf("DEBUG: N central jets %d\n", totalCentralJets);
	    
	    // Is there a good enough pair of jets?
	    if( !( bjet1 && bjet2 ) ) continue;
	    
	    // Cuts on the Higgs -> bjet1 bjet2
	    TLorentzVector b1Momentum, b2Momentum;
	    b1Momentum.SetPtEtaPhiM(bjet1->pt, bjet1->eta, bjet1->phi, bjet1->mass);
	    b2Momentum.SetPtEtaPhiM(bjet2->pt, bjet2->eta, bjet2->phi, bjet2->mass);
	    TLorentzVector HMomentum = b1Momentum + b2Momentum;

	    // Pt of the Higgs candidate
	    if( ! (HMomentum.Pt() > cutHPtMin ) ) continue;
	    // Pt of the Z candidate
	    if( !( ZMomentum.Pt() > cutZPtMin ) ) continue;
	    
	    // Apply b-tagging.
	    // First, order CSV of jets in magnitude
	    double csv1 = bjet1->csv;
	    double csv2 = bjet2->csv;
	    if(csv1 < csv2){
	      csv1 = bjet2->csv;
	      csv2 = bjet1->csv;
	    }
	    // Second, apply the cut on CSV
	    if( !( csv1 > csv1cut && csv2 > csv2cut ) ) continue;

	    // Z and H are back to back
	    if( ! ( fabs( ZMomentum.DeltaPhi(HMomentum) ) > cutDPhiZHMin ) ) continue;
	    
	    // Additional jet veto (no more than 1 additional jet)
	    if( ! (totalCentralJets - 2 <= cutAdditionalCentralJetsMax) ) continue;
	    
	    // The mass cut on dijet should be a sliding window if we follow AN 11-240
	    // The values below are for search at m_H = 115 GeV.
	    // NOTE: AN 11-240 does not make it clear whether this mass window
	    // should be applied to the pair of jets with highest b-jet probability,
	    // or whether the pair should be selected from those that passed this cut.
	    if( ! (HMomentum.M() > cutHMassMin && HMomentum.M() < cutHMassMax ) ) continue;
	    
	    //
	    // The ZH candidate for cut and count is identified
	    //
// 	    printf("DEBUG: ZH candidate found\n");
	    
	    //
	    // Fill histograms
	    // 
	    hZMassv[isam]->Fill(ZMomentum.M(),weight);
	    hHMassv[isam]->Fill(HMomentum.M(),weight);
	    hZPtv[isam]  ->Fill(ZMomentum.Pt(),  weight);
	    hHPtv[isam]  ->Fill(HMomentum.Pt(),  weight);
	    
	    //
	    // Fill ntuple data
	    //
	    fillData(&data, info, ZMomentum, electron1, electron2, bjet1, bjet2, HMomentum.M(), HMomentum.Pt(),
		     pvArr->GetEntriesFast(), totalCentralJets, weight);
	    outTree->Fill();
	    
	    nsel    += weight;
	    nselvar += weight*weight;
	    
	  } // end second lepton loop	 
	} // end first lepton loop
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
  delete electronArr;
  delete pvArr;

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
  
  // dielectron mass
  TCanvas *c1 = MakeCanvas("c1","c1",canw,canh);
  sprintf(ylabel,"Events / %.1f GeV/c^{2}",hZMassv[0]->GetBinWidth(1));
  CPlot plotZMass("mass","","m(e^{+}e^{-}) [GeV/c^{2}]",ylabel);
  if(hasData) { plotZMass.AddHist1D(hZMassv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    hZMassv[isam]->Scale(mcscale);
    plotZMass.AddToStack(hZMassv[isam],samplev[isam]->label,samplev[isam]->color);
  }
  if(samplev.size()>5)
    plotZMass.SetLegend(0.75,0.55,0.98,0.9);
  else
    plotZMass.TransLegend(0.1,0);
  if(lumi>0) plotZMass.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
  plotZMass.Draw(c1,kFALSE,format);
  
  // dijet mass
  TCanvas *c2 = MakeCanvas("c2","c2",canw,canh);
  sprintf(ylabel,"Events / %.1f GeV/c^{2}",hHMassv[0]->GetBinWidth(1));
  CPlot plotHMass("mass","","m(jj) [GeV/c^{2}]",ylabel);
  if(hasData) { plotHMass.AddHist1D(hHMassv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    hHMassv[isam]->Scale(mcscale);
    plotHMass.AddToStack(hHMassv[isam],samplev[isam]->label,samplev[isam]->color);
  }
  if(samplev.size()>5)
    plotHMass.SetLegend(0.75,0.55,0.98,0.9);
  else
    plotHMass.TransLegend(0.1,0);
  if(lumi>0) plotHMass.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
  plotHMass.Draw(c2,kFALSE,format);
  
  // dielectron pT
  TCanvas *c3 = MakeCanvas("c3","c3",canw,canh);
  sprintf(ylabel,"Events / %.1f GeV/c",hZPtv[0]->GetBinWidth(1));
  CPlot ploZPt("pt","","p_{T}(e^{+}e^{-}) [GeV/c]",ylabel);
  if(hasData) { ploZPt.AddHist1D(hZPtv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    hZPtv[isam]->Scale(mcscale);
    ploZPt.AddToStack(hZPtv[isam],samplev[isam]->label,samplev[isam]->color);
  }
  if(lumi>0) ploZPt.AddTextBox(lumitext,0.41,0.85,0.61,0.8,0);
  if(samplev.size()>5)
    ploZPt.SetLegend(0.75,0.55,0.98,0.9);
  else
    ploZPt.TransLegend(0.1,0);
  ploZPt.Draw(c3,kFALSE,format);

  // dijet pT
  TCanvas *c4 = MakeCanvas("c4","c4",canw,canh);
  sprintf(ylabel,"Events / %.1f GeV/c",hHPtv[0]->GetBinWidth(1));
  CPlot plotHPt("pt","","p_{T}(jj) [GeV/c]",ylabel);
  if(hasData) { plotHPt.AddHist1D(hHPtv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    hHPtv[isam]->Scale(mcscale);
    plotHPt.AddToStack(hHPtv[isam],samplev[isam]->label,samplev[isam]->color);
  }
  if(lumi>0) plotHPt.AddTextBox(lumitext,0.41,0.85,0.61,0.8,0);
  if(samplev.size()>5)
    plotHPt.SetLegend(0.75,0.55,0.98,0.9);
  else
    plotHPt.TransLegend(0.1,0);
  plotHPt.Draw(c4,kFALSE,format);

  TCanvas *c5 = MakeCanvas("c5","c5",canw,canh);
  CPlot plotNGoodPV("ngoodpv","Good PVs","N_{PV}","Events");
  if(hasData) { plotNGoodPV.AddHist1D(hNGoodPVv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    hNGoodPVv[isam]->Scale(mcscale);
    plotNGoodPV.AddToStack(hNGoodPVv[isam],samplev[isam]->label,samplev[isam]->color);
  }
  if(lumi>0) plotNGoodPV.AddTextBox(lumitext,0.45,0.85,0.65,0.8,0);
  plotNGoodPV.TransLegend(0.1,0);
//   plotNGoodPV.SetLogy();
  plotNGoodPV.Draw(c5,kFALSE,format);

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
    txtfile << "         " << setprecision(1) << fixed << nSelv[0] << " ZH events!" << endl;
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
        
  gBenchmark->Show("selectZHtoEEbb");       
} 


//=== FUNCTION DEFINITIONS ======================================================================================

//--------------------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------------------
void fillData(ZHtoEEbbData *data, const mithep::TEventInfo *info, TLorentzVector ZMomentum,
	      const mithep::TElectron *electron1, const mithep::TElectron *electron2, 
	      const mithep::TJet *jet1, const mithep::TJet *jet2, double dijetMass, double dijetPt,
              const UInt_t npv, const UInt_t njets, const Double_t weight)
{
  data->runNum         = info->runNum;
  data->evtNum         = info->evtNum;
  data->lumiSec        = info->lumiSec;
  data->nPV            = npv;
  data->nJets          = njets;                                        

  data->pfMET          = info->pfMET;
  data->pfMETphi       = info->pfMETphi;
  data->pfSumET        = info->pfSumET;

  // Z part
  data->diEMass           = ZMomentum.M();
  data->diEPt             = ZMomentum.Pt();
  data->diEY              = ZMomentum.Rapidity();
  data->diEPhi            = ZMomentum.Phi(); 

  data->pt_1           = electron1->pt;
  data->eta_1          = electron1->eta;
  data->phi_1          = electron1->phi;
  data->scEt_1         = electron1->scEt;
  data->scEta_1        = electron1->scEta;
  data->scPhi_1        = electron1->scPhi;
  data->hltMatchBits_1 = electron1->hltMatchBits;
  data->q_1            = electron1->q;

  data->pt_2           = electron2->pt;
  data->eta_2          = electron2->eta;
  data->phi_2          = electron2->phi;
  data->scEt_2         = electron2->scEt;
  data->scEta_2        = electron2->scEta;
  data->scPhi_2        = electron2->scPhi;
  data->hltMatchBits_2 = electron2->hltMatchBits;
  data->q_2            = electron2->q;

  // H part
  data->diJetMass      = dijetMass;
  data->diJetPt        = dijetPt;

  data->jPt_1           = jet1->pt          ;
  data->jEta_1          = jet1->eta         ;
  data->jPhi_1          = jet1->phi         ;
  data->jMass_1         = jet1->mass        ;
  data->tche_1          = jet1->tche        ;
  data->tchp_1          = jet1->tchp        ;
  data->jCSV_1          = -1.0              ; // Not available at the moment!
  data->mcFlavor_1      = jet1->mcFlavor    ;
  data->jHltMatchBits_1 = jet1->hltMatchBits;

  data->jPt_1          = jet2->pt          ;
  data->jEta_1         = jet2->eta         ;
  data->jPhi_1         = jet2->phi         ;
  data->jMass_1        = jet2->mass        ;
  data->tche_1         = jet2->tche        ;
  data->tchp_1         = jet2->tchp        ;
  data->jCSV_1         = -1.0              ; // Not available at the moment!
  data->mcFlavor_1     = jet2->mcFlavor    ;
  data->jHltMatchBits_1 =jet2->hltMatchBits;

  data->weight         = weight;
}

//--------------------------------------------------------------------------------------------------
