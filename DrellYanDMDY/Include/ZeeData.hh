#ifndef ZEE_DATA_HH
#define ZEE_DATA_HH

#define ZeeData_is_TObject

#ifdef ZeeData_is_TObject
#include <TObject.h>
#include "TEventInfo.hh"
#include "TDielectron.hh"

class ZeeData_t : public TObject
#else
struct ZeeData
#endif
{
public:
  UInt_t  runNum;                          // run number in data
  UInt_t  evtNum;                          // event number in data
  UInt_t  lumiSec;                         // lumi section      
  UInt_t  nTracks0;                        // number of reconstructed tracks in event
  UInt_t  nCaloTowers0;                    // number of reconstructed calorimeter towers in event
  UInt_t  nPV;                             // number of valid reconstructed primary vertices in event                                          
  UInt_t  nGoodPV;                         // number of good PVs
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

#ifdef ZeeData_is_TObject
  ClassDef(ZeeData_t,1)
#endif
};

//"runNum/i:evtNum:lumiSec:nTracks0:nCaloTowers0:nPV:nGoodPV:nJets:caloMEx/F:caloMEy:caloSumET:tcMEx:tcMEy:tcSumET:pfMEx:pfMEy:pfSumET:mass:pt:y:phi:pt_1:eta_1:phi_1:scEt_1:scEta_1:scPhi_1:hltMatchBits_1/i:q_1/I:pt_2/F:eta_2:phi_2:scEt_2:scEta_2:scPhi_2:hltMatchBits_2/i:q_2/I:weight/F"
//"runNum/i:evtNum:lumiSec:nTracks0:nCaloTowers0:nPV:nJets:caloMEx/F:caloMEy:caloSumET:tcMEx:tcMEy:tcSumET:pfMEx:pfMEy:pfSumET:mass:pt:y:phi:pt_1:eta_1:phi_1:scEt_1:scEta_1:scPhi_1:hltMatchBits_1/i:q_1/I:pt_2/F:eta_2:phi_2:scEt_2:scEta_2:scPhi_2:hltMatchBits_2/i:q_2/I:weight/F"

#ifdef ZeeData_is_TObject
inline
void fillData(ZeeData_t *data, const mithep::TEventInfo *info, const mithep::TDielectron *dielectron, 
              const UInt_t npv, const UInt_t nGoodPV, const UInt_t njets, const Double_t weight) {
  data->runNum         = info->runNum;
  data->evtNum         = info->evtNum;
  data->lumiSec        = info->lumiSec;
  data->nTracks0       = 0;
  data->nCaloTowers0   = 0;
  data->nPV            = npv;
  data->nGoodPV        = nGoodPV;
  data->nJets          = njets;
  data->caloMEx        = 0;
  data->caloMEy        = 0;
  data->caloSumET      = 0;
  data->tcMEx          = 0; //info->trkMET * cos(info->trkMETphi);
  data->tcMEy          = 0; //info->trkMET * sin(info->trkMETphi);
  data->tcSumET        = info->trkSumET;
  data->pfMEx          = 0; //info->pfMET * cos(info->pfMETphi);
  data->pfMEy          = 0; //info->pfMET * sin(info->pfMETphi);
  data->pfSumET        = info->pfSumET;
  data->mass           = dielectron->mass;
  data->pt             = dielectron->pt;
  data->y              = dielectron->y;
  data->phi            = dielectron->phi; 
  // no endorsement that the 1st electron is the leading one
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

#endif
