#ifndef UnfoldingTools
#define UnfoldingTools

#include <TString.h>
#include <TVectorD.h>

// all functions return three possible values:
//  -1 : the needed file could not be opened
//   1 : all checked ok. Calculation is done
//   0 : inconsistent binning detected

namespace unfolding{

  int  unfold(TVectorD &vin, TVectorD &vout, TString unfoldingConstFileName);

  int  propagateErrorThroughUnfolding(TVectorD &errorIn, 
					TVectorD &errorPropagated,
					TString unfoldingConstFileName);
  
  int calculateTotalUnfoldingSystError(TVectorD &yieldsBeforeUnfolding, 
				   TVectorD &systUnfolding, 
					TString fullUnfoldingConstFileName,
					TString extraUnfoldingErrorsFileName);

  int checkBinningConsistency(TString fileName);

}

#endif
