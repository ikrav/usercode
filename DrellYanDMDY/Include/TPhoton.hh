#ifndef EWKANA_NTUPLER_TPHOTON_HH
#define EWKANA_NTUPLER_TPHOTON_HH

#include <TObject.h>

namespace mithep 
{
  class TPhoton : public TObject
  {
    public:
      TPhoton():
      pt(0), eta(0), phi(0), scEt(0), scEta(0), scPhi(0), trkIso04(0), trkIso04NoPV(0), emIso04(0), hadIso04(0), HoverE(0),
      R9(0), sigiEtaiEta(0), sigiPhiiPhi(0), scID(0), hasPixelSeed(0), passSpikeFilter(0), passEleVetoConvRec(0), 
      passConvId(0), hltMatchBits(0)  
      {}
      ~TPhoton(){}

      Float_t pt, eta, phi; 	              // kinematics
      Float_t scEt, scEta, scPhi;             // supercluster
      Float_t trkIso04;                       // track isolation (hollow cone)
      Float_t trkIso04NoPV;                   // track isolation without PV-constraint
      Float_t emIso04;                        // ECAL-based isolation
      Float_t hadIso04;                       // HCAL-based isolation
      Float_t HoverE;		              // H/E
      Float_t R9;		              // ratio of energies in 3x3 to SC
      Float_t sigiEtaiEta;                    // eta-width of shower in number of crystals
      Float_t sigiPhiiPhi;                    // phi-width of shower in number of crystals           
      UInt_t  scID;                           // supercluster ID (for matching to electron superclusters)
      Bool_t  hasPixelSeed;                   // supercluster has pixel seed?
      Bool_t  passSpikeFilter;                // pass spike filter?
      Bool_t  passEleVetoConvRec;             // 
      Bool_t  passConvId;                     // 
      ULong_t hltMatchBits;  	              // bits from matching with HLT primitives 
    
    ClassDef(TPhoton,2)
  };  
}
#endif
