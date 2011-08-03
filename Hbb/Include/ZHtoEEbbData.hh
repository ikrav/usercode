#ifndef ZHTOEEBB_DATA_HH
#define ZHTOEEBB_DATA_HH

struct ZHtoEEbbData
{
  UInt_t  runNum;                          // run number in data
  UInt_t  evtNum;                          // event number in data
  UInt_t  lumiSec;                         // lumi section      
  UInt_t  nPV;                             // number of valid reconstructed primary vertices in event                                          
  UInt_t  nJets;                           // number of jets (with some requirements)
  Float_t pfMET, pfMETphi, pfSumET;        // particle flow MET

  // Z part
  Float_t diEMass, diEPt, diEY, diEPhi;                // dielectron kinematics    
  
  Float_t pt_1, eta_1, phi_1;              // leading electron
  Float_t scEt_1, scEta_1, scPhi_1;  
  UInt_t  hltMatchBits_1;
  Int_t   q_1;
  
  Float_t pt_2, eta_2, phi_2;              // lagging electron
  Float_t scEt_2, scEta_2, scPhi_2;
  UInt_t  hltMatchBits_2;
  Int_t   q_2;

  // Higgs-> bb part
  Float_t diJetMass, diJetPt;         // kinematics of two b jets

  Float_t jPt_1, jEta_1, jPhi_1, jMass_1;   // leading jet
  Float_t tche_1, tchp_1, jCSV_1;
  Int_t mcFlavor_1;
  ULong_t jHltMatchBits_1;

  Float_t jPt_2, jEta_2, jPhi_2, jMass_2;   // trailing jet
  Float_t tche_2, tchp_2, jCSV_2;
  Int_t mcFlavor_2;
  ULong_t jHltMatchBits_2;

  Float_t weight;                          // event weight  
};

//"runNum/i:evtNum:lumiSec:nPV:nJets:caloMEx/F:caloMEy:caloSumET:tcMEx:tcMEy:tcSumET:pfMEx:pfMEy:pfSumET:diEMass:diEPt:diEY:diEPhi:pt_1:eta_1:phi_1:scEt_1:scEta_1:scPhi_1:hltMatchBits_1/i:q_1/I:pt_2/F:eta_2:phi_2:scEt_2:scEta_2:scPhi_2:hltMatchBits_2/i:q_2/I:diJetMass/F:diJetPt:jPt_1:jEta_1:jPhi_1:jMass_1:tche_1:tchp_1:jCSV_1:mcFlavor_1/I:jHltMatchBits_1:jPt_2/F:jEta_2:jPhi_2:jMass_2:tche_2:tchp_2:jCSV_2:mcFlavor_2/I:jHltMatchBits_2:weight/F"

#endif
