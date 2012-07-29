#ifndef EWKANA_NTUPLER_TJET_HH
#define EWKANA_NTUPLER_TJET_HH

#include <TObject.h>

namespace mithep 
{
  class TJet : public TObject
  {
    public:
      TJet():
      pt(0), eta(0), phi(0), mass(0), rawPt(0), area(0), tche(0), tchp(0), dz(0), mcFlavor(0), hltMatchBits(0)
      {}
      ~TJet(){}

      Float_t pt, eta, phi, mass;  // kinematics
      Float_t rawPt;               // raw jet pT (before any correction)
      Float_t area;                // jet area
      Float_t tche;                // TrackCountingHighEfficiency b-tag discriminator
      Float_t tchp;                // TrackCountingHighPurity b-tag discriminator
      Float_t dz;                  // effective dz of the jet based on charged components
      Int_t   mcFlavor;            // PDG ID of matched parton flavor
      ULong_t hltMatchBits;        // bits from matching with HLT primitives

    ClassDef(TJet,2)
  };
}
#endif
