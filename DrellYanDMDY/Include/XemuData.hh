#ifndef XEMU_DATA_HH
#define XEMU_DATA_HH

struct XemuData
{
  UInt_t  runNum;                          // run number in data
  UInt_t  evtNum;                          // event number in data
  UInt_t  lumiSec;                         // lumi section      
  UInt_t  nPV;                             // number of valid reconstructed primary vertices in event                                          
  UInt_t  nGoodPV;                         // number of good reconstructed primary vertices in the event                                          
  UInt_t  nJets;                           // number of jets (with some requirements)

  Float_t pfSumET;

  Float_t mass, pt;                       // dilepton kinematics    
  Float_t rapidity;

  Float_t pt_mu, eta_mu, phi_mu;                   // muon
  UInt_t  hltMatchBits_mu;
  Int_t   q_mu;
  
  Float_t pt_e, eta_e, phi_e;                     // electron
  
  Float_t scEt_e, scEta_e, scPhi_e;
  UInt_t  hltMatchBits_e;
  Int_t   q_e;

  Float_t weight;                          // event weight  
};

//"runNum/i:evtNum:lumiSec:nPV:mass:pt:pt_mu:eta_mu:hltMatchBits_mu/i:q_mu/I:pt_e/F:eta_e:scEt_e:scEta_e:hltMatchBits_e/i:q_e/I:weight/F"

#endif
