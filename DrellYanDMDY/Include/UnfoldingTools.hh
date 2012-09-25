#ifndef UnfoldingTools
#define UnfoldingTools

#include <TString.h>
#include <TVectorD.h>
#include <TMatrixD.h>

#include "../Include/DYTools.hh"

// all functions return three possible values:
//  -1 : the needed file could not be opened
//   1 : all checked ok. Calculation is done
//   0 : inconsistent binning detected

namespace unfolding{

  int  unfold(const TVectorD &vin, TVectorD &vout, const TString &unfoldingConstFileName);
  int  unfoldTrueToReco(const TVectorD &vin, TVectorD &vout, const TString &unfoldingConstFileName);
  // pack MatrixD to flat-indexed VectorD and apply the unfolding matrix 
  // from a file
  int  unfold(const TMatrixD &vinM, TMatrixD &voutM, const TString &unfoldingConstFileName,
	      TVectorD &vinFlat, TVectorD &voutFlat);


  int  unfoldFSR(const TVectorD &vin, TVectorD &vout, const TString &unfoldingConstFileName, const TString &correctionsFileName);
  int  unfoldFSRTrueToReco(const TVectorD &vin, TVectorD &vout, const TString &unfoldingConstFileName, const TString &correctionsFileName);

  int  unfoldFSR(const TMatrixD &vinM, TMatrixD &voutM, 
		 const TString &unfoldingConstFileName, 
		 const TString &correctionsFileName,
		 TVectorD &vinFlat, TVectorD &voutFlat);

  int  propagateErrorThroughUnfolding(const TVectorD &errorIn, 
				      TVectorD &errorPropagated,
				      const TString &unfoldingConstFileName);
  

  int  propagateErrorThroughUnfolding(const TMatrixD &errorInMatrix, 
				      TMatrixD &errorPropagatedMatrix,
				      const TString &unfoldingConstFileName,
				      TVectorD &errorIn,
				      TVectorD &errorPropagated);

  int  propagateErrorThroughFsrUnfolding(const TVectorD &errorIn, 
				      TVectorD &errorPropagated,
					 const TString &unfoldingConstFileName,
					 const TString &correctionsFileName
					 );

  int  propagateErrorThroughFsrUnfolding(const TMatrixD &errorInMatrix, 
				      TMatrixD &errorPropagatedMatrix,
				      const TString &unfoldingConstFileName,
				      const TString &correctionsFileName,
				      TVectorD &errorIn,
					 TVectorD &errorPropagated);
  

  /* int calculateTotalUnfoldingSystError(const TVectorD &yieldsBeforeUnfolding, 
				       TVectorD &systUnfolding, 
				       const TString &fullUnfoldingConstFileName,
				       const TString &extraUnfoldingErrorsFileName);*/

  int calculateTotalUnfoldingSystErrorFlat
          (const TVectorD &yieldsBeforeUnfolding, 
	   TVectorD &systUnfolding, 
	   const TString &fullUnfoldingConstFileName,
	   const TString &extraUnfoldingErrorsFileName);
 
  int checkBinningConsistency(const TString &fileName);

  //  convert m[nMassBins][ybins] -> v[flat_idx]
  int flattenMatrix(const TMatrixD &m, TVectorD &v) {
    for (int i=0; i<DYTools::nMassBins; ++i) {
      for (int yi=0; yi<DYTools::nYBins[i]; ++yi) {
	int flatIdx=DYTools::findIndexFlat(i,yi);
	v[flatIdx]=m[i][yi];
      }
    }
    return 1;
  }

  //  convert v[flat_idx] -> m[nMassBins][ybins]
  int deflattenMatrix(const TVectorD &v, TMatrixD &m) {
    for (int i=0; i<DYTools::nMassBins; ++i) {
      for (int yi=0; yi<DYTools::nYBins[i]; ++yi) {
	int flatIdx=DYTools::findIndexFlat(i,yi);
	m[i][yi]=v[flatIdx];
      }
    }
    return 1;
  }

  void writeBinningArrays(TFile &fout);
  bool checkBinningArrays(TFile &fin);

  // assume M[mass][rapidity]
  bool checkRangesMY(const TMatrixD &M, const TString &name, int verb=1);
  bool checkRangesMY(const TMatrixD &M1, const TString &name1, const TMatrixD &M2, const TString &name2, int verb=1);
  // assume Vector[flat_index]
  bool checkRangesFI(const TVectorD &M, const TString &name, int verb=1);

  // check that mass and rapidityCounts match massBinLimits and nYBins
  bool checkBinningRanges(const TVectorD &mass, const TVectorD &rapidityCounts, const TString &fname);

  // check that the vectors have totalNumberOfBins elements
  bool checkFlatVectorRanges(const TVectorD &v1, const TVectorD &v2, const TString &info, const TString info2="");
}


// -------------------------------------------------------

// -------------------------------------------------------

#endif
