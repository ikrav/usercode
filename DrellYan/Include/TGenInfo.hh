#ifndef EWKANA_NTUPLER_TGENINFO_HH
#define EWKANA_NTUPLER_TGENINFO_HH

#include <TObject.h>

namespace mithep 
{
  // Generator level info data object
  class TGenInfo : public TObject
  {
    public:
      TGenInfo(){}
      ~TGenInfo(){}
      
      UInt_t  npho;                      // number of FSR photons
      Int_t   id_1, id_2;                // parton ID
      Float_t x_1, x_2;		         // parton momentum fraction
      Float_t weight;		         // event weight
      Float_t vmass, vpt, vy, vphi;      // boson info
      Float_t mass, pt, y, phi;          // dilepton info
      Float_t pt_1, eta_1, phi_1;        // lepton info
      Float_t pt_2, eta_2, phi_2;  
      Float_t phopt, phoeta, phophi;     // leading photon kinematics
      Float_t decx, decy, decz;	         // boson decay vertex

      // SC matched to electrons
      Float_t scEt_1, scEta_1;
      Float_t scEt_2, scEta_2;
      Float_t scMass;
      	  
    ClassDef(TGenInfo,2)
  };
}
#endif
