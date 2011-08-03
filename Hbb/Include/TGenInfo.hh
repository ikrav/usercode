#ifndef MITHTT_NTUPLER_TGENINFO_HH
#define MITHTT_NTUPLER_TGENINFO_HH

#include <TObject.h>

namespace mithep 
{
  // Generator level info data object
  class TGenInfo : public TObject
  {
    public:
      TGenInfo(){}
      ~TGenInfo(){}
      
      Float_t weight;				// event weight
      Int_t   pid_1, pid_2;			// parton ID
      Float_t x_1, x_2;				// parton momentum fraction
      Int_t   id_a, id_b;			// boson IDs
      Float_t vmass_a, vpt_a, vy_a, vphi_a;	// boson A info
      Float_t vmass_b, vpt_b, vy_b, vphi_b;	// boson B info
      Int_t   id_1_a, id_2_a;			// lepton/quark IDs
      Int_t   id_1_b, id_2_b;			// lepton/quark IDs
      Float_t pt_1_a, eta_1_a, phi_1_a;		// lepton info
      Float_t pt_2_a, eta_2_a, phi_2_a;  
      Float_t pt_1_b, eta_1_b, phi_1_b;
      Float_t pt_2_b, eta_2_b, phi_2_b;  
      Float_t decx, decy, decz;			// boson decay vertex
      	  
    ClassDef(TGenInfo,2)
  };
}
#endif
