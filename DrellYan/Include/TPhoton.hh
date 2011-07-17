#ifndef EWKANA_NTUPLER_TPHOTON_HH
#define EWKANA_NTUPLER_TPHOTON_HH

#include <TObject.h>

namespace mithep 
{
  class TPhoton : public TObject
  {
    public:
      TPhoton(){}
      ~TPhoton(){}

      Float_t pt, eta, phi; 	              // kinematics
      Float_t scEt, scEta, scPhi;             // supercluster
      Float_t trkIso04;                       // track isolation (hollow cone)
      Float_t emIso04;                        // ECAL-based isolation
      Float_t hadIso04;                       // HCAL-based isolation
      Float_t HoverE;		              // H/E
      Float_t R9;		              // ratio of energies in 3x3 to SC
      Float_t sigiEtaiEta;                    // eta-width of shower in number of crystals            
      UInt_t  scID;                           // supercluster ID (for matching to electron superclusters)
      Bool_t  hasPixelSeed;                   // supercluster has pixel seed?     
      ULong_t hltMatchBits;  	              // bits from matching with HLT primitives 
    
    ClassDef(TPhoton,1)
  };  
}
#endif
