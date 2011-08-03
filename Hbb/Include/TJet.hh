#ifndef MITHTT_NTUPLER_TJET_HH
#define MITHTT_NTUPLER_TJET_HH

#include <TObject.h>

namespace mithep 
{
  class TJet : public TObject
  {
    public:
      TJet(){}
      ~TJet(){}

      Float_t pt, eta, phi, mass;  // kinematics
      Float_t unc;                 // energy scale uncertainty
      Float_t area;                // jet area
      Float_t tche;                // TrackCountingHighEfficiency b-tag discriminator
      Float_t tchp;                // TrackCountingHighPurity b-tag discriminator
      Float_t csv;                 // Combined secondary vertex BJet Tags Disc
      Float_t csvMva;			// Combined Secondary Vertex MVA BJet Tags Disc
      UInt_t nCharged;			// ChargedMultiplicity
      Float_t chgEMfrac, neuEMfrac;	// charged and neutral EM energy fractions
      Float_t chgHadrfrac, neuHadrfrac; // charged and neutral Hadr energy fractions
      Int_t   mcFlavor;			// PDG ID of matched parton flavor
      UInt_t  hltMatchBits;		// bits from matching with HLT primitives

    ClassDef(TJet,4)
  };
}
#endif
