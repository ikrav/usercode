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
#include "../Include/TGenInfo.hh"
#include "../Include/TElectron.hh"
#include "../Include/TJet.hh"
#include "../Include/TVertex.hh"
#include "../Include/MyTools.hh"

// Helper functions for Electron ID selection
#include "../Include/EleIDCuts.hh" 

#endif

//=== FUNCTION DECLARATIONS ======================================================================================

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

void efficiencyZHtoEEbb(const TString input) 
{  
  gBenchmark->Start("efficiencyZHtoEEbb");

  
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================
  
  TString  format = "png";            // plot format

  vector<TString> fnamev;   // file names   
  vector<TString> labelv;   // legend label
  vector<Int_t>   colorv;   // color in plots
  vector<Int_t>   linev;    // line style
  vector<Double_t> xsecv;
  vector<Double_t> lumiv;

  ifstream ifs;
  ifs.open(input.Data());
  assert(ifs.is_open());
  string line;
  while(getline(ifs,line)) {
    if(line[0]=='#') continue;
    string fname;
    Int_t color, linesty;
    stringstream ss(line);
    Double_t xsec;
    ss >> fname >> xsec >> color >> linesty;
    string label = line.substr(line.find('@')+1);
    fnamev.push_back(fname);
    labelv.push_back(label);
    colorv.push_back(color);
    linev.push_back(linesty);
    xsecv.push_back(xsec);
    lumiv.push_back(0);
  }
  ifs.close();

  const Double_t kGAP_LOW  = 1.4442;
  const Double_t kGAP_HIGH = 1.566;
    
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
  
  //
  // Set up histograms and counters
  //

  double totalGenEvents     = 0;
  double passEventTrigger   = 0;
  double passLooseEE        = 0;
  double passInAcceptanceEE = 0;
  double passEETriggerMatch = 0;
  double passEEID           = 0;
  double passZMassWindow    = 0;
  double passDijetExists    = 0;
  // final selection
  double passHiggsPt        = 0;
  double passZPt            = 0;
  double passBTag1          = 0;
  double passBTag2          = 0;
  double passZHdphi         = 0;
  double passExtraJets      = 0;
  double passHMass          = 0;
  
  //
  // Access samples and fill histograms
  //  
  TFile *infile=0;
  TTree *eventTree=0;  
  
  // Data structures to store info from TTrees
  mithep::TEventInfo *info    = new mithep::TEventInfo();
  mithep::TGenInfo *gen       = new mithep::TGenInfo();
  TClonesArray *pfJetArr      = new TClonesArray("mithep::TJet");
  TClonesArray *electronArr   = new TClonesArray("mithep::TElectron");
  TClonesArray *pvArr         = new TClonesArray("mithep::TVertex");
  
  //
  // loop over files
  //
  const UInt_t nfiles = fnamev.size();    
  for(UInt_t ifile=0; ifile<nfiles; ifile++) {
    cout << "Processing " << fnamev[ifile] << "... "; cout.flush();
    infile = new TFile(fnamev[ifile]);
    assert(infile);
    
    // Get the TTree
    eventTree = (TTree*)infile->Get("Events"); assert(eventTree);

    // Set branch address to structures that will store the info  
    eventTree->SetBranchAddress("Info",       &info);          TBranch *infoBr       = eventTree->GetBranch("Info");
    eventTree->SetBranchAddress("Gen",&gen);                   TBranch *genBr = eventTree->GetBranch("Gen");
    eventTree->SetBranchAddress("Electron",   &electronArr);   TBranch *electronBr   = eventTree->GetBranch("Electron");
    eventTree->SetBranchAddress("PFJet",      &pfJetArr);      TBranch *pfJetBr      = eventTree->GetBranch("PFJet");
    eventTree->SetBranchAddress("PV",         &pvArr);         TBranch *pvBr         = eventTree->GetBranch("PV");

   // Find weight for events for this file
    // The first file in the list comes with weight 1,
    // all subsequent ones are normalized to xsection and luminosity
    lumiv[ifile] = eventTree->GetEntries()/xsecv[ifile];
    double scale = lumiv[0]/lumiv[ifile];
    cout << "       -> sample weight is " << scale << endl;
      
    // loop through events
    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {       
      //       for(UInt_t ientry=0; ientry<100000; ientry++) {  // for debugging
      
      infoBr->GetEntry(ientry);
      genBr->GetEntry(ientry);

      totalGenEvents += gen->weight * scale;

      // The double electron trigger bit definitions
      UInt_t eventTriggerBit          = kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL;
      UInt_t leadingTriggerObjectBit  = kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele1Obj;
      UInt_t trailingTriggerObjectBit = kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele2Obj;
	
      if(!(info->triggerBits & eventTriggerBit)) continue;  // no trigger accept? Skip to next event...                                   
      passEventTrigger  += gen->weight * scale;

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
	  
	  passLooseEE += gen->weight * scale;
	  // Exclude ECAL gap region and cut out of acceptance electrons
	  if((fabs(electron1->scEta)>kGAP_LOW) && (fabs(electron1->scEta)<kGAP_HIGH)) continue;
	  if((fabs(electron2->scEta)>kGAP_LOW) && (fabs(electron2->scEta)<kGAP_HIGH)) continue;
	  if((fabs(electron1->scEta) > cutEleEtaMax) || (fabs(electron2->scEta) > cutEleEtaMax))       continue;  // outside eta range? Skip to next event...
	  //
	  // ET thresholds for electrons
	  if( ! (electron1->scEt > cutEleETMin && electron2->scEt > cutEleETMin) ) continue;

	  passInAcceptanceEE  += gen->weight * scale;
	  // Both electrons must match trigger objects. At least one ordering
	  // must match
	  if( ! ( 
		 (electron1->hltMatchBits & leadingTriggerObjectBit && 
		  electron2->hltMatchBits & trailingTriggerObjectBit )
		 ||
		 (electron1->hltMatchBits & trailingTriggerObjectBit && 
		  electron2->hltMatchBits & leadingTriggerObjectBit ) ) ) continue;
	  
	  passEETriggerMatch += gen->weight * scale;
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
	  
	  passEEID += gen->weight * scale;
	  // mass window for Z candidate
	  if((ZMomentum.M() < cutZMassMin) || (ZMomentum.M() > cutZMassMax)) continue;
	  
	  passZMassWindow  += gen->weight * scale;

	  // Selection of Z candidate is complete. 
	  
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
	      // Jet Pt cuts (JEC corrections are already applied before)
	      // NOT CLEAR: do leptons from Z need to be excluded from jets
	      // in addition to min separation cut below? Description in 11-240 not totally clear.
	      if( !( jet1->pt > cutJetPTMin && jet2->pt > cutJetPTMin ) ) continue;
		
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
	  
	  // Is there a good enough pair of jets?
	  if( !( bjet1 && bjet2 ) ) continue;
	  passDijetExists += gen->weight * scale;

	  // Preselection is complete. Start final selection.	  

	  TLorentzVector b1Momentum, b2Momentum;
	  b1Momentum.SetPtEtaPhiM(bjet1->pt, bjet1->eta, bjet1->phi, bjet1->mass);
	  b2Momentum.SetPtEtaPhiM(bjet2->pt, bjet2->eta, bjet2->phi, bjet2->mass);
	  TLorentzVector HMomentum = b1Momentum + b2Momentum;
	  
	  // Pt of the Higgs candidate
	  if( ! (HMomentum.Pt() > cutHPtMin ) ) continue;
	  passHiggsPt += gen->weight * scale;
	  // Pt of the Z candidate
	  if( !( ZMomentum.Pt() > cutZPtMin ) ) continue;
	  passZPt += gen->weight * scale;

	  // Apply b-tagging.
	  // First, order CSV of jets in magnitude
	  double csv1 = bjet1->csv;
	  double csv2 = bjet2->csv;
	  if(csv1 < csv2){
	    csv1 = bjet2->csv;
	    csv2 = bjet1->csv;
	  }
	  // Second, apply the cut on CSV
	  if( !( csv1 > csv1cut) ) continue;
	  passBTag1 += gen->weight * scale;
	  if( !( csv2 > csv2cut ) ) continue;
	  passBTag2 += gen->weight * scale;

	  // Z and H are back to back
	  if( ! ( fabs( ZMomentum.DeltaPhi(HMomentum) ) > cutDPhiZHMin ) ) continue;
	  passZHdphi += gen->weight * scale;
	  
	  // Additional jet veto (no more than 1 additional jet)
	  if( ! (totalCentralJets - 2 <= cutAdditionalCentralJetsMax) ) continue;
	  passExtraJets += gen->weight * scale;

	  // The mass cut on dijet should be a sliding window if we follow AN 11-240
	  // The values below are for search at m_H = 115 GeV.
	  // NOTE: AN 11-240 does not make it clear whether this mass window
	  // should be applied to the pair of jets with highest b-jet probability,
	  // or whether the pair should be selected from those that passed this cut.
	  if( ! (HMomentum.M() > cutHMassMin && HMomentum.M() < cutHMassMax ) ) continue;
	  passHMass += gen->weight * scale;

	  //
	  // The ZH candidate for cut and count is identified
	  //
	    
	  //
	  // Fill histograms if needed
	  // 
	  
	} // end second lepton loop	 
      } // end first lepton loop
    }  
    delete infile;
    infile=0, eventTree=0;
  }
  delete info;
  delete electronArr;
  delete pvArr;

  //
  // Summary printout
  // 
  printf("\n");
  printf("Total generated events              %10.0f\n",totalGenEvents     );
  printf("Preselection:\n");
  printf("   events that pass trigger         %10.0f\n",passEventTrigger   );
  printf("   number of loose dielectrons      %10.0f\n",passLooseEE        );
  printf("   dielectrons passing ET/eta cuts  %10.0f\n",passInAcceptanceEE );
  printf("   dielectrons matched to trigger   %10.0f\n",passEETriggerMatch );
  printf("   dielectrons passing ID           %10.0f\n",passEEID           );
  printf("   dielectrons in Z mass window     %10.0f\n",passZMassWindow    );
  printf("   at least one dijet is found      %10.0f, %5.1f %%\n",passDijetExists, 100*passDijetExists/totalGenEvents    );
  printf("Full selection:\n");
  printf("   passed dijet Pt cut              %10.0f, %5.1f %%\n",passHiggsPt, 100*passHiggsPt/passDijetExists);
  printf("   passed Z Pt cut                  %10.0f, %5.1f %%\n",passZPt    , 100*passZPt/passHiggsPt        );
  printf("   pass 1st jet (stronger) btag     %10.0f, %5.1f %%\n",passBTag1  , 100*passBTag1/passZPt       );
  printf("   pass 2d jet (weaker) btag        %10.0f, %5.1f %%\n",passBTag2  , 100*passBTag2/passBTag1       );
  printf("   pass Z-H deltaPhi                %10.0f, %5.1f %%\n",passZHdphi , 100*passZHdphi/passBTag2       );
  printf("   pass extra jet veto              %10.0f, %5.1f %%\n",passExtraJets,100*passExtraJets/passZHdphi      );
  printf("   pass dijet mass window           %10.0f, %5.1f %%\n",passHMass  , 100*passHMass/passExtraJets        );
  printf("\n");
  printf("Total event efficiency: %5.2f %%\n", 100*passHMass/totalGenEvents);
  printf("\n");
  
  gBenchmark->Show("efficiencyZHtoEEbb");       
} 


//=== FUNCTION DEFINITIONS ======================================================================================

//--------------------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------------------
