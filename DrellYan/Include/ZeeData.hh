#ifndef ZEE_DATA_HH
#define ZEE_DATA_HH

struct ZeeData
{
  UInt_t  runNum;                          // run number in data
  UInt_t  evtNum;                          // event number in data
  UInt_t  lumiSec;                         // lumi section      
  UInt_t  nTracks0;                        // number of reconstructed tracks in event
  UInt_t  nCaloTowers0;                    // number of reconstructed calorimeter towers in event
  UInt_t  nPV;                             // number of valid reconstructed primary vertices in event                                          
  UInt_t  nGoodPV;                         // number of good reconstructed primary vertices in event
  UInt_t  nJets;                           // number of jets (with some requirements)
  Float_t caloMEx, caloMEy, caloSumET;     // calorimeter MET 
  Float_t tcMEx, tcMEy, tcSumET;           // track-corrected MET
  Float_t pfMEx, pfMEy, pfSumET;           // particle flow MET

  Float_t mass, pt, y, phi;                // dielectron kinematics    
  
  Float_t pt_1, eta_1, phi_1;              // leading electron
  Float_t scEt_1, scEta_1, scPhi_1;  
  UInt_t  hltMatchBits_1;
  Int_t   q_1;
  
  Float_t pt_2, eta_2, phi_2;              // lagging electron
  Float_t scEt_2, scEta_2, scPhi_2;
  UInt_t  hltMatchBits_2;
  Int_t   q_2;

  Float_t weight;                          // event weight  

};

// updated (added nGoodPV):
//"runNum/i:evtNum:lumiSec:nTracks0:nCaloTowers0:nPV:nGoodPV:nJets:caloMEx/F:caloMEy:caloSumET:tcMEx:tcMEy:tcSumET:pfMEx:pfMEy:pfSumET:mass:pt:y:phi:pt_1:eta_1:phi_1:scEt_1:scEta_1:scPhi_1:hltMatchBits_1/i:q_1/I:pt_2/F:eta_2:phi_2:scEt_2:scEta_2:scPhi_2:hltMatchBits_2/i:q_2/I:weight/F"

// outdated: 2012 Mar 05
//"runNum/i:evtNum:lumiSec:nTracks0:nCaloTowers0:nPV:nJets:caloMEx/F:caloMEy:caloSumET:tcMEx:tcMEy:tcSumET:pfMEx:pfMEy:pfSumET:mass:pt:y:phi:pt_1:eta_1:phi_1:scEt_1:scEta_1:scPhi_1:hltMatchBits_1/i:q_1/I:pt_2/F:eta_2:phi_2:scEt_2:scEta_2:scPhi_2:hltMatchBits_2/i:q_2/I:weight/F"

#endif
