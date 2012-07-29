#ifndef EWKANA_NTUPLER_TVERTEX_HH
#define EWKANA_NTUPLER_TVERTEX_HH

#include <TObject.h>

namespace mithep 
{
  class TVertex : public TObject
  {
    public:
      TVertex():
      nTracksFit(0), ndof(0), chi2(0), sumPt(0), x(0), y(0), z(0)
      {}
      ~TVertex(){}

      UInt_t  nTracksFit;
      Float_t ndof;      
      Float_t chi2;
      Float_t sumPt;
      Float_t x,y,z;					 

    ClassDef(TVertex,2)
  };
}
#endif
