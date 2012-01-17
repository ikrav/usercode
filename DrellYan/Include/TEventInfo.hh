#ifndef EWKANA_NTUPLER_TEVENTINFO_HH
#define EWKANA_NTUPLER_TEVENTINFO_HH

#include <TObject.h>

namespace mithep 
{
  class TEventInfo : public TObject
  {
    public:
      TEventInfo(){}
      ~TEventInfo(){}

      UInt_t  runNum; 			             // run number in data
      UInt_t  evtNum; 			             // event number in data
      UInt_t  lumiSec;			             // lumi section
      UInt_t  nPU, nPUminus, nPUplus;                // number of reconstructed pile up vertices in event (MC only)       
      Float_t pvx, pvy, pvz;		             // best primary vertex
      Float_t bsx, bsy, bsz;		             // beamspot
      Float_t pfMET, pfMETphi, pfSumET;	             // particle flow MET
      Float_t trkMET, trkMETphi, trkSumET;           // track MET
      Float_t rhoLowEta, rhoHighEta;                 // average energy density for isolation correction
      ULong_t triggerBits;		             // HLT trigger bits
      Bool_t  hasGoodPV;                             // event has a good PV?

    ClassDef(TEventInfo,2)
  };
}
#endif
