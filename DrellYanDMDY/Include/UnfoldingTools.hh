#ifndef UnfoldingTools
#define UnfoldingTools

#include <TString.h>
#include <TVectorD.h>

namespace unfolding{

  void  unfold(TVectorD &vin, TVectorD &vout, TString unfoldingConstFileName);

  void  propagateErrorThroughUnfolding(TVectorD &errorIn, 
					TVectorD &errorPropagated,
					TString unfoldingConstFileName);
  
  void calculateTotalUnfoldingSystError(TVectorD &yieldsBeforeUnfolding, 
				   TVectorD &systUnfolding, 
					TString fullUnfoldingConstFileName,
					TString extraUnfoldingErrorsFileName);

  bool checkBinningConsistency(TString fileName);

}

#endif
