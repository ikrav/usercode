#ifndef EWKANA_NTUPLER_TGENINFO_HH
#define EWKANA_NTUPLER_TGENINFO_HH

#include <TObject.h>

namespace mithep 
{
  // Generator level info data object
  class TGenInfo : public TObject
  {
    public:
      TGenInfo():
      npho(0), id_1(0), id_2(0), lid_1(0), lid_2(0), x_1(0), x_2(0), weight(0), 
      vmass(0), vpt(0), vy(0), vphi(0), vpt_1(0), veta_1(0), vphi_1(0), vpt_2(0), veta_2(0), vphi_2(0),
      mass(0), pt(0), y(0), phi(0), pt_1(0), eta_1(0), phi_1(0), pt_2(0), eta_2(0), phi_2(0),
      phopt(0), phoeta(0), phophi(0), decx(0), decy(0), decz(0)
      {}
      ~TGenInfo(){}
      
      UInt_t  npho;                      // number of FSR photons
      Int_t   id_1, id_2;                // parton ID
      Int_t   lid_1, lid_2;              // lepton ID
      Float_t x_1, x_2;		         // parton momentum fraction
      Float_t weight;		         // event weight
      Float_t vmass, vpt, vy, vphi;      // boson info
      Float_t vpt_1, veta_1, vphi_1;
      Float_t vpt_2, veta_2, vphi_2;
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
