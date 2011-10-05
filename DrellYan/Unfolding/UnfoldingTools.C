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
  void  unfold(TVectorD &vin, TVectorD &vout, TString unfoldingConstFileName)
  {
    
    // Read unfolding constants
    printf("unfold: Load constants\n"); fflush(stdout);
    
    if(!checkBinningConsistency(unfoldingConstFileName))
      assert(0);
    
    TFile fileConstants(unfoldingConstFileName);
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
    return;
  }

//-----------------------------------------------------------------
// Function that propagates systematic errors through unfolding
//-----------------------------------------------------------------
  void  propagateErrorThroughUnfolding(TVectorD &errorIn, 
					TVectorD &errorPropagated,
					TString unfoldingConstFileName)
  {

    // Read unfolding constants
    printf("unfold: Load constants\n"); fflush(stdout);
    
    if(!checkBinningConsistency(unfoldingConstFileName))
      assert(0);

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
    return;
  }

  // This function adds together all pieces of unfolding systematics
  void calculateTotalUnfoldingSystError(TVectorD &yieldsBeforeUnfolding, 
				   TVectorD &systUnfolding, 
				   TString fullUnfoldingConstFileName,
				   TString extraUnfoldingErrorsFileName){

    if(!checkBinningConsistency(fullUnfoldingConstFileName))
      assert(0);

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
      assert(0);
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
  }

  bool checkBinningConsistency(TString fileName){

    bool result = true;

    TFile fileConstants(fileName);
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
      result = false;
    }else
      printf("unfold: Binning in the inputs is consistent\n");
    
    return result;
  }

} // end of namespace 
