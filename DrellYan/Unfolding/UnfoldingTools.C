//
// This file contains tools related to unfolding
//
#include <TFile.h>
#include <TMatrixD.h>

#include "../Include/UnfoldingTools.hh"
#include "../Include/DYTools.hh"

namespace unfolding {

//-----------------------------------------------------------------
// Function that executes unfolding
//-----------------------------------------------------------------
  int  unfold(const TVectorD &vin, TVectorD &vout, TString unfoldingConstFileName)
  {
    
    // Read unfolding constants
    printf("unfold: Load constants from <%s>\n",unfoldingConstFileName.Data()); fflush(stdout);
    
    int res=checkBinningConsistency(unfoldingConstFileName);
    if (res!=1) return res;
    
    TFile fileConstants(unfoldingConstFileName); // file had to exist to reach this point
    TMatrixD DetResponse             = *(TMatrixD *)fileConstants.FindObjectAny("DetResponse");
    TMatrixD DetInvertedResponse     = *(TMatrixD *)fileConstants.FindObjectAny("DetInvertedResponse");
    TMatrixD DetInvertedResponseErr  = *(TMatrixD *)fileConstants.FindObjectAny("DetInvertedResponseErr");
    TVectorD BinLimitsArray          = *(TVectorD *)fileConstants.FindObjectAny("BinLimitsArray");

    int nBins = BinLimitsArray.GetNoElements()-1;
    // Apply unfolding matrix
    vout = 0;
    for(int i=0; i<nBins; i++){
      for(int j=0; j<nBins; j++){
	vout[i] += DetInvertedResponse(j,i) * vin[j];
	
      }
    }
    
    fileConstants.Close();
    return 1;
  }

//-----------------------------------------------------------------
// Function that propagates systematic errors through unfolding
//-----------------------------------------------------------------
 int  propagateErrorThroughUnfolding(const TVectorD &errorIn, 
					TVectorD &errorPropagated,
					TString unfoldingConstFileName)
  {

    // Read unfolding constants
    printf("propagateErrorThroughUnfolding: Load constants from <%s>\n",unfoldingConstFileName.Data()); fflush(stdout);
    
    int res=checkBinningConsistency(unfoldingConstFileName);
    if (res!=1) return res;

    TFile fileConstants(unfoldingConstFileName);
    TMatrixD DetResponse             = *(TMatrixD *)fileConstants.FindObjectAny("DetResponse");
    TMatrixD DetInvertedResponse     = *(TMatrixD *)fileConstants.FindObjectAny("DetInvertedResponse");
    TMatrixD DetInvertedResponseErr  = *(TMatrixD *)fileConstants.FindObjectAny("DetInvertedResponseErr");
    TVectorD BinLimitsArray          = *(TVectorD *)fileConstants.FindObjectAny("BinLimitsArray");
    int nBins = BinLimitsArray.GetNoElements()-1;

    // Apply unfolding matrix
    errorPropagated = 0;
    for(int i=0; i<nBins; i++){
      for(int j=0; j<nBins; j++){
	errorPropagated[i] += pow( DetInvertedResponse   (j,i) * errorIn[j], 2);
      }
      errorPropagated[i] = sqrt(errorPropagated[i]);
    }
    
    fileConstants.Close();
    return 1;
  }

  // ------------------------------------------

  // This function adds together all pieces of unfolding systematics
  int calculateTotalUnfoldingSystError(const TVectorD &yieldsBeforeUnfolding, 
				   TVectorD &systUnfolding, 
				   TString fullUnfoldingConstFileName,
				   TString extraUnfoldingErrorsFileName){

    int res=checkBinningConsistency(fullUnfoldingConstFileName);
    if (res!=1) return res;

    TFile fileConstants(fullUnfoldingConstFileName);
    TMatrixD DetResponse             = *(TMatrixD *)fileConstants.FindObjectAny("DetResponse");
    TMatrixD DetInvertedResponse     = *(TMatrixD *)fileConstants.FindObjectAny("DetInvertedResponse");
    TMatrixD DetInvertedResponseErr  = *(TMatrixD *)fileConstants.FindObjectAny("DetInvertedResponseErr");
    TVectorD BinLimitsArray          = *(TVectorD *)fileConstants.FindObjectAny("BinLimitsArray");
    int nBins = BinLimitsArray.GetNoElements()-1;

    // Estimate unfolding error due to uncertainty of unfolding matrix elements
    TVectorD systElementsError(nBins);
    systElementsError = 0;

    for(int i=0; i<nBins; i++){
      for(int j=0; j<nBins; j++){
	systElementsError[i] += pow( DetInvertedResponseErr(j,i) * yieldsBeforeUnfolding[j], 2);
      }
      systElementsError[i] = sqrt(systElementsError[i]);
    }

    // Read relative unfolding systematic error from a file.
    // This error covers various aspects of unfolding systematics
    // not related to MC statistics. It needs to be calculated
    // separately.
    TFile fileExtraUnfoldingErrors(extraUnfoldingErrorsFileName);
    if( ! fileExtraUnfoldingErrors.IsOpen()){
      printf("ERROR: required file with unfolding errors %s is not found!\n",
	     extraUnfoldingErrorsFileName.Data());
      return -1;
    }
    TVectorD systOtherSourcesPercent 
      = *(TVectorD *)fileExtraUnfoldingErrors.FindObjectAny("unfoldingSystPercent");
    if( systOtherSourcesPercent.GetNoElements() != DYTools::nMassBins){
      printf("ERROR: Wrong binning of the unfolding systematics array!\n");
      assert(0);
    }
    // For absolute errors we need to know unfolded yields
    TVectorD yieldsAfterUnfolding(nBins);
    unfold(yieldsBeforeUnfolding, yieldsAfterUnfolding, fullUnfoldingConstFileName);
    // Calculate absolute error from other sources
    TVectorD systOtherSources(nBins);
    for(int i=0; i<nBins; i++){
      systOtherSources[i] = (systOtherSourcesPercent[i]/100.0) * yieldsAfterUnfolding[i];

    }

    // Add all pieces of unfolding systematics together
    for(int i=0; i<nBins; i++){
      systUnfolding[i] = sqrt( systElementsError[i]*systElementsError[i] 
			       + systOtherSources[i] * systOtherSources[i]);
    }
    fileExtraUnfoldingErrors.Close();
    fileConstants.Close();
    return 1;
  }

  // ------------------------------------------

  int checkBinningConsistency(TString fileName){

    int result = 1;

    TFile fileConstants(fileName);
    if (!fileConstants.IsOpen()) return -1;
    TMatrixD DetResponse             = *(TMatrixD *)fileConstants.FindObjectAny("DetResponse");
    TMatrixD DetInvertedResponse     = *(TMatrixD *)fileConstants.FindObjectAny("DetInvertedResponse");
    TMatrixD DetInvertedResponseErr  = *(TMatrixD *)fileConstants.FindObjectAny("DetInvertedResponseErr");
    TVectorD BinLimitsArray          = *(TVectorD *)fileConstants.FindObjectAny("BinLimitsArray");
    
    // Check that the binning is consistent
    bool checkResult = true;
    int nBins = BinLimitsArray.GetNoElements()-1;
    if( DetResponse.GetNrows()    != DYTools::nMassBins ) checkResult = false;
    if( BinLimitsArray.GetNoElements() != DYTools::nMassBins+1) checkResult = false;
    for(int i=0; i<nBins+1; i++){
      if( BinLimitsArray[i] != DYTools::massBinLimits[i] )
	checkResult = false;
    }

    fileConstants.Close();
    if( !checkResult ){
      printf("unfold: ERROR: inconsistent binning in the inputs\n");
      result = 0;
    }else
      printf("unfold: Binning in the inputs is consistent\n");
    
    return result;
  }

} // end of namespace 
