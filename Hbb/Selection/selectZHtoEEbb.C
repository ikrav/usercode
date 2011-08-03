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
#include "../Include/EWKAnaDefs.hh"
#include "../Include/TEventInfo.hh"
#include "../Include/TDielectron.hh"
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
void fillData(ZHtoEEbbData *data, const mithep::TEventInfo *info, const mithep::TDielectron *dielectron,
	      const mithep::TJet *jet1, const mithep::TJet *jet2, double dijetMass, double dijetPt,
              const UInt_t npv, const UInt_t njets, const Double_t weight);

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
    sprintf(hname,"hHMass_%i",isam);   hHMassv.push_back(new TH1F(hname,"",50,50,250));   hHMassv[isam]->Sumw2();
    sprintf(hname,"hZPt_%i",isam);     hZPtv.push_back(new TH1F(hname,"",50,0,100));       hZPtv[isam]->Sumw2();
    sprintf(hname,"hHPt_%i",isam);     hHPtv.push_back(new TH1F(hname,"",50,0,100));       hHPtv[isam]->Sumw2();
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
  TClonesArray *dielectronArr = new TClonesArray("mithep::TDielectron");
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
      eventTree->SetBranchAddress("Dielectron", &dielectronArr); TBranch *dielectronBr = eventTree->GetBranch("Dielectron");
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
	UInt_t eventTriggerBit = kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL 
	  | kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL;
	UInt_t leadingTriggerObjectBit = kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele1Obj
	  | kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele1Obj;
	UInt_t trailingTriggerObjectBit = kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele2Obj
	  | kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele2Obj;
	if(isam==0) {
	  nProcessedEvents++;
 	}
	
        if(!(info->triggerBits & eventTriggerBit)) continue;  // no trigger accept? Skip to next event...                                   

	dielectronArr->Clear(); 
	dielectronBr->GetEntry(ientry);	

	// 
	// Loop through dielectrons
	//
        for(Int_t i=0; i<dielectronArr->GetEntriesFast(); i++) {
	  mithep::TDielectron *dielectron = (mithep::TDielectron*)((*dielectronArr)[i]);

          // Exclude ECAL gap region and cut out of acceptance electrons
          if((fabs(dielectron->scEta_1)>kGAP_LOW) && (fabs(dielectron->scEta_1)<kGAP_HIGH)) continue;
          if((fabs(dielectron->scEta_2)>kGAP_LOW) && (fabs(dielectron->scEta_2)<kGAP_HIGH)) continue;
          if((fabs(dielectron->scEta_1) > 2.5)       || (fabs(dielectron->scEta_2) > 2.5))       continue;  // outside eta range? Skip to next event...
	  //
	  // ET thresholds for electrons
       	  if( ! (dielectron->scEt_1 > 20 && dielectron->scEt_2 > 20) ) continue;
	  // Both electrons must match trigger objects. At least one ordering
	  // must match
	  if( ! ( 
		 (dielectron->hltMatchBits_1 & leadingTriggerObjectBit && 
		  dielectron->hltMatchBits_2 & trailingTriggerObjectBit )
		  ||
		 (dielectron->hltMatchBits_1 & trailingTriggerObjectBit && 
		  dielectron->hltMatchBits_2 & leadingTriggerObjectBit ) ) ) continue;

	  // Other cuts to both electrons

	  // Apply electron ID to both electrons
	  bool useSmurf = false;
	  if(useSmurf){
	    // The Smurf electron ID package as used in HWW analysis.
	    // It contains cuts like VBTF WP80 for pt>20, VBTF WP70 for pt<10
	    // with some customization, plus impact parameter cuts dz and dxy
	    if(!passSmurf(dielectron)) continue;  
	  }else{
	    // Use ID cuts as in AN 11-240
	    // The second argument is for using relative combined isolation
	    // Note: relative combined isolation cut in this call is a bit
	    // looser than in 11-240. The values for other cuts haev nto been
	    // checked.
	    if( !passWP95(dielectron, kTRUE)) continue;
	  }

          // mass window for Z candidate
          if((dielectron->mass < 75) || (dielectron->mass > 105)) continue;
	  
	  TLorentzVector lep1Momentum, lep2Momentum;
	  lep1Momentum.SetPtEtaPhiM(dielectron->pt_1, dielectron->eta_1, dielectron->phi_1, 0.000511);
	  lep2Momentum.SetPtEtaPhiM(dielectron->pt_2, dielectron->eta_2, dielectron->phi_2, 0.000511);
	  TLorentzVector ZMomentum = lep1Momentum + lep2Momentum;

	  // Pt of the Z candidate
	  if( !( ZMomentum.Pt() > 100 ) ) continue;

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
	  double dRJetLeptonMin = 0.3;
	  int totalCentralJets = 0;
 	  pfJetArr->Clear();
 	  pfJetBr->GetEntry(ientry);
 	  for(Int_t ijet=0; ijet<pfJetArr->GetEntriesFast(); ijet++) {
 	    const mithep::TJet *jet1 = (mithep::TJet*)((*pfJetArr)[ijet]);
	    for(Int_t jjet=0; jjet<pfJetArr->GetEntriesFast(); jjet++) {
	      if( ijet == jjet ) continue;
	      const mithep::TJet *jet2 = (mithep::TJet*)((*pfJetArr)[jjet]);
	  
	      // Cuts on jet kinematics
	      if( !( fabs(jet1->eta) < 2.5 && fabs(jet2->eta) < 2.5) ) continue;
// 	      printf("DEBUG: pair of central jets found\n");
	      // Jet Pt cuts (JEC corrections are already applied before)
	      // NOT CLEAR: do leptons from Z need to be excluded from jets
	      // in addition to min separation cut below? Description in 11-240 not totally clear.
	      if( !( jet1->pt > 20 && jet2->pt > 20 ) ) continue;
// 	      printf("DEBUG: Pt is high enough\n");

	      // Cuts presently missing
	      // TRACK MULTIPLICITY >= 2
	      // FRACTION OF EM ENERGY >= 1%
	      // FRACTION OF HAD ENERGY >= 1%      

	      // Make sure leptons from Z are not too close to jets
	      if( toolbox::deltaR(jet1->eta, jet1->phi, dielectron->scEta_1, dielectron->scPhi_1) < dRJetLeptonMin ) continue;
	      if( toolbox::deltaR(jet1->eta, jet1->phi, dielectron->scEta_2, dielectron->scPhi_2) < dRJetLeptonMin ) continue;
	      if( toolbox::deltaR(jet2->eta, jet2->phi, dielectron->scEta_1, dielectron->scPhi_1) < dRJetLeptonMin ) continue;
	      if( toolbox::deltaR(jet2->eta, jet2->phi, dielectron->scEta_2, dielectron->scPhi_2) < dRJetLeptonMin ) continue;
// 	      printf("DEBUG: Not matched to leptons\n");

	      // b-tagging
	      bool useCSV = false;
	      bool isBTag1 = false;
	      bool isBTag2 = false;
	      if(useCSV){
		// The primary method is to be implemented once variables are available
		// in the ntuples. This is what 11-240 uses.
	      } else {
		// Use the impact parameter significance of the second most significant
		// track ("high efficiency")
		if(fabs(jet1->tche)>2)
		  isBTag1 = true;
		if(fabs(jet2->tche)>2)
		  isBTag2 = true;
	      }
	      if(! (isBTag1 && isBTag2 ) ) continue;
// 	      printf("DEBUG: B-tagged\n");
	      
	      // Save if this is highest significance pair so far
	      double thisSignificance = -1;
	      if(useCSV){
		// Here we should add CSV values of the two jets
	      } else {
		thisSignificance = fabs(jet1->tche) + fabs(jet1->tche);
	      }
	      if(thisSignificance > bestSignificance){
		bjet1 = jet1;
		bjet2 = jet2;
		bestSignificance = thisSignificance;
	      }
	    } // end inner loop over jets

	    // Count jets toward jet veto
	    if( jet1->pt > 20 && fabs(jet1->eta) < 2.5
		&& toolbox::deltaR(jet1->eta, jet1->phi, dielectron->scEta_1, dielectron->scPhi_1) > dRJetLeptonMin
		&& toolbox::deltaR(jet1->eta, jet1->phi, dielectron->scEta_2, dielectron->scPhi_2) > dRJetLeptonMin 
		){
	      totalCentralJets++;
	    } 

	  } // end outer loop over jets
// 	  printf("DEBUG: N central jets %d\n", totalCentralJets);

	  // Is there a good enough pair of b-jets?
	  if( !( bjet1 && bjet2 ) ) continue;
// 	  printf("DEBUG: two b-jets found\n");

	  // Cuts on the Higgs -> bjet1 bjet2
	  TLorentzVector b1Momentum, b2Momentum;
	  b1Momentum.SetPtEtaPhiM(bjet1->pt, bjet1->eta, bjet1->phi, bjet1->mass);
	  b2Momentum.SetPtEtaPhiM(bjet2->pt, bjet2->eta, bjet2->phi, bjet2->mass);
	  TLorentzVector HMomentum = b1Momentum + b2Momentum;
	  if( ! (HMomentum.Pt() > 100 ) ) continue;
	  // The mass cut on dijet should be a sliding window if we follow AN 11-240
	  // The values below are for search at m_H = 120 GeV.
	  // NOTE: AN 11-240 does not make it clear whether this mass window
	  // should be applied to the pair of jets with highest b-jet probability,
	  // or whether the pair should be selected from those that passed this cut.
	  if( ! (HMomentum.M() > 100 && HMomentum.M() < 130 ) ) continue;
	  
	  // Z and H are back to back
	  if( ! ( fabs( ZMomentum.DeltaPhi(HMomentum) ) > 2.90 ) ) continue;

	  // Additional jet veto (no more than 1 additional jet)
	  if( ! (totalCentralJets - 2 < 2) ) continue;
	  
	  //
	  // The ZH candidate for cut and count is identified
	  //
	  printf("DEBUG: ZH candidate found\n");

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
	  
	  //
	  // Fill histograms
	  // 
	  hZMassv[isam]->Fill(dielectron->mass,weight);
	  hHMassv[isam]->Fill(HMomentum.M(),weight);
          hZPtv[isam]  ->Fill(dielectron->pt,  weight);
          hHPtv[isam]  ->Fill(HMomentum.Pt(),  weight);
          hNGoodPVv[isam]->Fill(nGoodPV,weight);
	  
	  //
	  // Fill ntuple data
	  //
	  fillData(&data, info, dielectron, bjet1, bjet2, HMomentum.M(), HMomentum.Pt(),
		   pvArr->GetEntriesFast(), totalCentralJets, weight);
	  outTree->Fill();
	  
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
    outFile->Write();
    delete outTree;
    outFile->Close();        
    delete outFile;
  }
  delete info;
  delete dielectronArr;
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
  
  // dielectron mass
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
  plotZMass.Draw(c,kFALSE,format);
  
  // dijet mass
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
  plotHMass.Draw(c,kFALSE,format);
  
  // dielectron pT
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
  ploZPt.Draw(c,kFALSE,format);

  // dijet pT
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
  plotHPt.Draw(c,kFALSE,format);

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
void fillData(ZHtoEEbbData *data, const mithep::TEventInfo *info, const mithep::TDielectron *dielectron, 
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
  data->diEMass           = dielectron->mass;
  data->diEPt             = dielectron->pt;
  data->diEY              = dielectron->y;
  data->diEPhi            = dielectron->phi; 

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
