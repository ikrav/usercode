//
// This file contains tools related to unfolding
//
#include <TFile.h>
#include <TMatrixD.h>

#include "../Include/UnfoldingTools.hh"

namespace unfolding {

//-----------------------------------------------------------------
// Function that executes unfolding
//-----------------------------------------------------------------
  int  unfold(TVectorD &vin, TVectorD &vout, const TString &unfoldingConstFileName)
  {
    
    // Read unfolding constants
    std::cout << "unfold: Load constants from <" << unfoldingConstFileName 
	      << ">" << std::endl;
    
    int res=checkBinningConsistency(unfoldingConstFileName);
    if (res!=1) return res;
    if (!checkFlatVectorRanges(vin,vout,unfoldingConstFileName," operates")) return 0;
    
    TFile fileConstants(unfoldingConstFileName); // file had to exist to reach this point
    //TMatrixD DetResponse             = *(TMatrixD *)fileConstants.FindObjectAny("DetResponse");
    TMatrixD DetInvertedResponse     = *(TMatrixD *)fileConstants.FindObjectAny("DetInvertedResponse");
    //TMatrixD DetInvertedResponseErr  = *(TMatrixD *)fileConstants.FindObjectAny("DetInvertedResponseErr");

    // Apply unfolding matrix
    vout = 0;
    int nBins = DYTools::getTotalNumberOfBins();
    for(int i=0; i<nBins; i++){
      for(int j=0; j<nBins; j++){
	vout[i] += DetInvertedResponse(j,i) * vin[j];
	
      }
    }
    
    fileConstants.Close();
    return 1;
  }

  //-----------------------------------------------------------------
 
 int  unfoldTrueToReco(TVectorD &vin, TVectorD &vout, const TString &unfoldingConstFileName)
  {

    // Read unfolding constants
    std::cout << "unfold: Load constants from <" << unfoldingConstFileName 
	      << ">" << std::endl;

    int res=checkBinningConsistency(unfoldingConstFileName);
    if (res!=1) return res;
    if (!checkFlatVectorRanges(vin,vout,unfoldingConstFileName," operates")) return 0;

    TFile fileConstants(unfoldingConstFileName); // file had to exist to reach this point
   TMatrixD DetResponse             = *(TMatrixD *)fileConstants.FindObjectAny("DetResponse");
   //TMatrixD DetInvertedResponse     = *(TMatrixD *)fileConstants.FindObjectAny("DetInvertedResponse");
   //TMatrixD DetInvertedResponseErr  = *(TMatrixD *)fileConstants.FindObjectAny("DetInvertedResponseErr");

   // Apply unfolding matrix
   vout = 0;
   int nBins = DYTools::getTotalNumberOfBins();
   for(int i=0; i<nBins; i++){
     for(int j=0; j<nBins; j++){
       vout[i] += DetResponse(j,i) * vin[j];
     }
   }

   fileConstants.Close();
   return 1;
 }


//-----------------------------------------------------------------
// Function that packs matrix to a vector and executes the unfolding
//-----------------------------------------------------------------

  int  unfold(TMatrixD &vinM, TMatrixD &voutM, const TString &unfoldingConstFileName, TVectorD &vin, TVectorD &vout)
  {
    
    // Read unfolding constants
    std::cout << "unfold(M): Load constants from <" << unfoldingConstFileName << ">" << std::endl;
    
    int res=checkBinningConsistency(unfoldingConstFileName);
    if (res!=1) return res;
    if (!checkRangesMY(vinM,"vinM",voutM,"voutM") ||
	!checkFlatVectorRanges(vin,vout,unfoldingConstFileName," derived flat vectors")) return 0;
    
    TFile fileConstants(unfoldingConstFileName); // file had to exist to reach this point
    TMatrixD DetResponse             = *(TMatrixD *)fileConstants.FindObjectAny("DetResponse");
    TMatrixD DetInvertedResponse     = *(TMatrixD *)fileConstants.FindObjectAny("DetInvertedResponse");
    TMatrixD DetInvertedResponseErr  = *(TMatrixD *)fileConstants.FindObjectAny("DetInvertedResponseErr");

    // Pack the matrix
    if (!flattenMatrix(vinM, vin)) return 0;

    // Apply unfolding matrix
    vout = 0;
    int nBins = DYTools::getTotalNumberOfBins();
    for(int i=0; i<nBins; i++){
      for(int j=0; j<nBins; j++){
	vout[i] += DetInvertedResponse(j,i) * vin[j];
      }
    }

    // Unpack the final matrix
    if (!deflattenMatrix(vout, voutM)) return 0;
    
    fileConstants.Close();
    return 1;
  }

//-----------------------------------------------------------------
// Function that propagates systematic errors through unfolding
//-----------------------------------------------------------------
 int  propagateErrorThroughUnfolding(TVectorD &errorIn, 
				     TVectorD &errorPropagated,
				     const TString &unfoldingConstFileName)
  {

    // Read unfolding constants
    std::cout << "propagateErrorThroughUnfolding: Load constants from <" << unfoldingConstFileName << ">" << std::endl;
    
    int res=checkBinningConsistency(unfoldingConstFileName);
    if (res!=1) return res;
    if (!checkFlatVectorRanges(errorIn,errorPropagated,unfoldingConstFileName," operates")) return 0;

    TFile fileConstants(unfoldingConstFileName);
    TMatrixD DetResponse             = *(TMatrixD *)fileConstants.FindObjectAny("DetResponse");
    TMatrixD DetInvertedResponse     = *(TMatrixD *)fileConstants.FindObjectAny("DetInvertedResponse");
    TMatrixD DetInvertedResponseErr  = *(TMatrixD *)fileConstants.FindObjectAny("DetInvertedResponseErr");

    // Apply unfolding matrix
    int nBins = DYTools::getTotalNumberOfBins();
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

//-----------------------------------------------------------------
// Function that propagates flattens the matrix and propages the 
// systematic errors through unfolding
//-----------------------------------------------------------------
  int  propagateErrorThroughUnfolding(TMatrixD &errorInM,
				      TMatrixD &errorPropagatedM,
				      const TString &unfoldingConstFileName,
				      TVectorD &errorIn, 
				      TVectorD &errorPropagated)
  {

    // Read unfolding constants
    std::cout << "propagateErrorThroughUnfolding(M): Load constants from <" << unfoldingConstFileName << ">" << std::endl;
    
    int res=checkBinningConsistency(unfoldingConstFileName);
    if (res!=1) return res;
    if (!checkRangesMY(errorInM,"errorInM",errorPropagatedM,"errorPropagatedM") ||
	!checkFlatVectorRanges(errorIn,errorPropagated,unfoldingConstFileName," operates")) return 0;

    TFile fileConstants(unfoldingConstFileName);
    TMatrixD DetResponse             = *(TMatrixD *)fileConstants.FindObjectAny("DetResponse");
    TMatrixD DetInvertedResponse     = *(TMatrixD *)fileConstants.FindObjectAny("DetInvertedResponse");
    TMatrixD DetInvertedResponseErr  = *(TMatrixD *)fileConstants.FindObjectAny("DetInvertedResponseErr");

    // Pack the matrix
    if (!flattenMatrix(errorInM, errorIn)) return 0;

    // Apply unfolding matrix
    int nBins = DYTools::getTotalNumberOfBins();
    errorPropagated = 0;
    for(int i=0; i<nBins; i++){
      for(int j=0; j<nBins; j++){
	errorPropagated[i] += pow( DetInvertedResponse   (j,i) * errorIn[j], 2);
      }
      errorPropagated[i] = sqrt(errorPropagated[i]);
    }
    
    // Unpack the final matrix
    if (!deflattenMatrix(errorPropagated, errorPropagatedM)) return 0;

    fileConstants.Close();
    return 1;
  }

  // ------------------------------------------
  /*
  // This function adds together all pieces of unfolding systematics
  int calculateTotalUnfoldingSystError(TVectorD &yieldsBeforeUnfolding, 
				   TVectorD &systUnfolding, 
				   const TString &fullUnfoldingConstFileName,
				   const TString &extraUnfoldingErrorsFileName){

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
  */

  // ------------------------------------------

  // This function adds together all pieces of unfolding systematics
  int calculateTotalUnfoldingSystErrorFlat(
	  TVectorD &yieldsBeforeUnfolding, 
	  TVectorD &systUnfolding, 
	  const TString &fullUnfoldingConstFileName,
	  const TString &extraUnfoldingErrorsFileName){

    int res=checkBinningConsistency(fullUnfoldingConstFileName);
    if (res!=1) return res;
    if (!checkFlatVectorRanges(yieldsBeforeUnfolding,systUnfolding,"calculateTotalUnfoldingSystErrorFlat")) return 0;
  
    TFile fileConstants(fullUnfoldingConstFileName);
    TMatrixD DetResponse             = *(TMatrixD *)fileConstants.FindObjectAny("DetResponse");
    TMatrixD DetInvertedResponse     = *(TMatrixD *)fileConstants.FindObjectAny("DetInvertedResponse");
    TMatrixD DetInvertedResponseErr  = *(TMatrixD *)fileConstants.FindObjectAny("DetInvertedResponseErr");

    int nBins = DYTools::getTotalNumberOfBins();
 
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
    TMatrixD* systOtherSourcesPercentPtr
      = (TMatrixD *)fileExtraUnfoldingErrors.FindObjectAny("unfoldingSystPercent");
    assert(systOtherSourcesPercentPtr);
    assert(checkRangesMY(*systOtherSourcesPercentPtr,
			 TString("unfoldingSystPercent from ") + extraUnfoldingErrorsFileName));

    // For absolute errors we need to know unfolded yields
    TVectorD yieldsAfterUnfolding(nBins);
    unfold(yieldsBeforeUnfolding, yieldsAfterUnfolding, fullUnfoldingConstFileName);
    // Calculate absolute error from other sources
    TVectorD systOtherSourcesPercentV(nBins),systOtherSourcesV(nBins);
    flattenMatrix(*systOtherSourcesPercentPtr, systOtherSourcesPercentV);
    for(int i=0; i<nBins; i++){
      systOtherSourcesV[i] = (systOtherSourcesPercentV[i]/100.0) * yieldsAfterUnfolding[i];
    }

    // Add all pieces of unfolding systematics together
    for(int i=0; i<nBins; i++){
      systUnfolding[i] = sqrt( systElementsError[i]*systElementsError[i] 
			       + systOtherSourcesV[i] * systOtherSourcesV[i]);
    }
    fileExtraUnfoldingErrors.Close();
    fileConstants.Close();

    if (systOtherSourcesPercentPtr) delete systOtherSourcesPercentPtr;
    return 1;
  }

  // ------------------------------------------

  int checkBinningConsistency(const TString &fileName){

    int result = 1;

    TFile fileConstants(fileName);
    if (!fileConstants.IsOpen()) return -1;
    result=(checkBinningArrays(fileConstants)) ? 1:0;
    TMatrixD *DetResponsePtr             = (TMatrixD *)fileConstants.FindObjectAny("DetResponse");
    TMatrixD *DetInvertedResponsePtr     = (TMatrixD *)fileConstants.FindObjectAny("DetInvertedResponse");
    TMatrixD *DetInvertedResponseErrPtr  = (TMatrixD *)fileConstants.FindObjectAny("DetInvertedResponseErr");

    if (!DetResponsePtr || !DetInvertedResponsePtr || !DetInvertedResponseErrPtr) {
      std::cout << "unfolding::checkBinningConsistency: failed to locate DetResponse, DetInvertedResponse, or DetInvertedResponseErr on file <" << fileName << ">\n";

      result=-1;
    }
    
    // Check that the binning is consistent
    bool checkResult = (result==1) ? true : false;
    int nBins = DYTools::getTotalNumberOfBins();
    if ( DetResponsePtr && ( DetResponsePtr->GetNrows() != nBins )) checkResult = false;

    fileConstants.Close();
    if( !checkResult ){
      std::cout << "unfold: ERROR: inconsistent binning in the inputs\n";
      if (result==1) result = 0;
    }else {
      //if (result==1) std::cout << "unfold: Binning in the inputs is consistent\n";
      //else std::cout << "unfold: CheckBinningConsistency -- Code ERROR\n";
      if (result!=1) std::cout << "unfold: CheckBinningConsistency -- Code ERROR\n";
    }
    
    if (DetResponsePtr) delete DetResponsePtr;
    if (DetInvertedResponsePtr) delete DetInvertedResponsePtr;
    if (DetInvertedResponseErrPtr) delete DetInvertedResponseErrPtr;

    return result;
  }

// -------------------------------------------------------

  void writeBinningArrays(TFile &fout) {
    fout.cd();
    TVectorD mass(nMassBins+1);
    TVectorD rapidityCounts(nMassBins);
    for (int i=0; i<nMassBins+1; i++) mass[i]=massBinLimits[i];
    for (int i=0; i<nMassBins  ; i++) rapidityCounts[i]=nYBins[i];
    mass.Write("massBinning");
    rapidityCounts.Write("rapidityCounts");
  }

  // -------------------------------------------------------

  bool checkBinningArrays(TFile &fin) {
    const char *fncname="unfolding::checkBinningArrays: ";
    TString fileInfo=TString("on file <") + fin.GetName() + ">";
    fin.cd();
    TVectorD* mass= (TVectorD*)fin.FindObjectAny("massBinning");
    TVectorD* rapidityCounts= (TVectorD*)fin.FindObjectAny("rapidityCounts");
    if (!mass || !rapidityCounts) {
      std::cout << fncname << "file <" << fin.GetName() << "> does not contain keys massBinning and/or rapidityCounts\n";
      return false;
    }
    bool comparisonOk=checkBinningRanges(*mass,*rapidityCounts, fin.GetName());
    if (mass) delete mass;
    if (rapidityCounts) delete rapidityCounts;
    return comparisonOk;
  }

  // -------------------------------------------------------

  bool checkBinningRanges(const TVectorD &mass, const TVectorD &rapidityCounts, const TString &fname) {
    const char *fncname="unfolding::checkBinningRanges: ";
    TString fileInfo=TString("on file <") + fname + ">";

    bool massOk=true, rapidityOk=true;
    if (mass.GetNoElements() != nMassBins+1) {
      std::cout << fncname 
		<< " number of mass bins " << fileInfo
		<< " is " << mass.GetNoElements() 
		<< " while " << (nMassBins+1) << " is expected\n";
      massOk=false;
    }
    
    if (rapidityCounts.GetNoElements() != nMassBins ) {
      std::cout << fncname
		<< "number of mass bins in rapidityCounts " << fileInfo
		<< " is " << rapidityCounts.GetNoElements() 
		<< " while " << nMassBins << " is expected\n";
      rapidityOk=false;
    }

    if (massOk) {
      for(int i=0; i<nMassBins+1; i++){
	if( massBinLimits[i] != mass[i] ) {
	  std::cout << fncname
		    << " mass limit " << fileInfo
		    << " at i=" << i 
		    << " is " << mass[i] 
		    << " instead of expected " << massBinLimits[i] << "\n";
	  massOk=false;
	}
      }
    }
    if (rapidityOk) {
      for (int i=0; i<nMassBins; i++) {
	if ( nYBins[i] != rapidityCounts[i] ) {
	  std::cout << fncname 
		    << "y bin count " << fileInfo
		    << " at i=" << i 
		    << " is " << rapidityCounts[i] 
		    << " instead of expected " << nYBins[i] << "\n";
	  rapidityOk=false;
	}
      }
    }
    return (massOk && rapidityOk);
  }

  // -----------------------------------------

  bool checkRangesMY(const TMatrixD &M, const TString &name, int verb) {
    bool ok=true;
    if (verb) std::cout << "checkRangesMY: check matrix with name=<" << name << ">\n";
    int nMaxYBins= DYTools::findMaxYBins();
    if(( M.GetNrows() != nMassBins ) || ( M.GetNcols() != nMaxYBins )) {
      std::cout << "checkRangesMY(TMatrixD): " << name << " is not good: "
		<< M.GetNrows() << 'x' << M.GetNcols() << " instead of "
		<< nMassBins << 'x' << nMaxYBins << "\n";
      ok=false;
    }
    return ok;
  }

  // -----------------------------------------

  bool checkRangesMY(const TMatrixD &M1, const TString &name1, const TMatrixD &M2, const TString &name2, int verb) {
    if (verb) std::cout << "checkRangesMY(" << name1 << ',' << name2 << ")\n";
    bool ok= (checkRangesMY(M1,name1,verb) && checkRangesMY(M2,name2,verb));
    return ok;
  }

  // -----------------------------------------

  // check that the vectors have totalNumberOfBins elements
  bool checkFlatVectorRanges(const TVectorD &v1, const TVectorD &v2, const TString &info, const TString info2) {
    int nBins = DYTools::getTotalNumberOfBins();
    if ( (v1.GetNoElements() != nBins) ||
	 (v2.GetNoElements() != nBins) ) {
      std::cout << "unfolding::checkFlatVectorRanges(info=<"
		<< info << info2 << ">)"
		<< ": the provided variables have sizes "
		<< v1.GetNoElements() << " and "
		<< v2.GetNoElements() 
		<< " instead of expected getTotalNumberOfBins="
		<< nBins << "\n";
      return false;
    }
    return true;
  }


  // -----------------------------------------

} // end of namespace 


