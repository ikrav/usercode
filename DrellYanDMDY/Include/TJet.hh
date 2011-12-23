#ifndef EWKANA_NTUPLER_TJET_HH
#define EWKANA_NTUPLER_TJET_HH

#include <TObject.h>

namespace mithep 
{
  class TJet : public TObject
  {
    public:
      TJet(){}
      ~TJet(){}

      Float_t pt, eta, phi, mass;  // kinematics
      Float_t area;                // jet area
      Float_t tche;                // TrackCountingHighEfficiency b-tag discriminator
      Float_t tchp;                // TrackCountingHighPurity b-tag discriminator
      Int_t   mcFlavor;            // PDG ID of matched parton flavor
      ULong_t hltMatchBits;        // bits from matching with HLT primitives

    ClassDef(TJet,1)
  };
}
#endif
