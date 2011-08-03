#ifndef MITHTT_NTUPLER_TEVENTINFO_HH
#define MITHTT_NTUPLER_TEVENTINFO_HH

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
      UInt_t  nPU;                                   // number of reconstructed pile up vertices in event (MC only)
      UInt_t  triggerBits;		             // HLT trigger bits 
      Float_t pvx, pvy, pvz;		             // best primary vertex
      Float_t bsx, bsy, bsz;		             // beamspot
      Float_t pfMET, pfMETphi, pfSumET;	             // particle flow MET
      Float_t trkMET, trkMETphi, trkSumET;           // track MET
      Float_t rho;                                   // average energy density for isolation correction
      Bool_t  hasGoodPV;                             // event has a good PV?

    ClassDef(TEventInfo,1)
  };
}
#endif
