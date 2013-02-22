//
// This file computes the effect of pile-up reweighting
// on each analysis step.
//
// Inputs: the three text files containing the yields
// recorded for each analysis step. The files correspond
// to the default pile-up reweighting, and the cases when
// pile-up reweighting is varied.
//
// 
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <assert.h>
#include "TString.h"
#include "TObjString.h"
#include "TObjArray.h"

const int nbins = 40;
const double limits[nbins+1] = 
  {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 
   81, 86, 91, 96, 101, 106, 110, 115, 120, 126, 133, 141, 
   150, 160, 171, 185, 200, 220, 243, 273, 320, 380, 440, 
   510, 600, 1000, 1500}; // 40 bins


void getArrays(TString dirName, TString fileName,
	       double *nObs, double *nUnf, 
	       double *nEffCorr, double *nAccCorr,
	       double *nFsrCorr);

void getYieldsArray(vector<TString> &fileContent, int ncol, double *data);
void checkObservedYieldsMatch(double *nObsDef, double *nObsPlus, double *nObsMinus);
void computeError(double *dataRefDef, double *dataCorrDef,
		  double *dataRefPlus, double *dataCorrPlus,
		  double *dataRefMinus, double *dataCorrMinus,
		  double *absError, double *relError);

void pileupSystematicsBreakdown(){
  
  TString dirNameDef = "/home/hep/ikrav/releases/another_UserCode_v2/UserCode/ikravchenko/DrellYanDMDY/CrossSection/tables1D/";
  TString fileNameDef = "tables1DrawsignalunfoldedeffcorrectedacccorrectedFSRcorrected-2.txt";
  
  TString dirNamePlus = "/home/hep/ikrav/releases/another_UserCode_v2_pileup_syst/UserCode/ikravchenko/DrellYanDMDY/CrossSection/tables1D/";
  TString fileNamePlus = "tables1DrawsignalunfoldedeffcorrectedacccorrectedFSRcorrected-4.txt";
  
  TString dirNameMinus = "/home/hep/ikrav/releases/another_UserCode_v2_pileup_syst_minus/UserCode/ikravchenko/DrellYanDMDY/CrossSection/tables1D/";
  TString fileNameMinus = "tables1DrawsignalunfoldedeffcorrectedacccorrectedFSRcorrected-4.txt";
  
  double nObsDef[nbins], nUnfDef[nbins],nEffCorrDef[nbins], 
    nAccCorrDef[nbins], nFsrCorrDef[nbins];

  double nObsPlus[nbins], nUnfPlus[nbins], nEffCorrPlus[nbins],
    nAccCorrPlus[nbins], nFsrCorrPlus[nbins];

  double nObsMinus[nbins], nUnfMinus[nbins], nEffCorrMinus[nbins], 
    nAccCorrMinus[nbins], nFsrCorrMinus[nbins]; 
  
  getArrays(dirNameDef, fileNameDef, nObsDef, nUnfDef, 
	    nEffCorrDef, nAccCorrDef, nFsrCorrDef);
  getArrays(dirNamePlus, fileNamePlus, nObsPlus, nUnfPlus, 
	    nEffCorrPlus, nAccCorrPlus, nFsrCorrPlus);
  getArrays(dirNameMinus, fileNameMinus, nObsMinus, nUnfMinus, 
	    nEffCorrMinus, nAccCorrMinus, nFsrCorrMinus);
  

  // Check that the observed yields, the starting point,
  // is identical for all variations of pile-up
  checkObservedYieldsMatch(nObsDef, nObsPlus, nObsMinus);

  // Compute errors
  double normUnfError[nbins]={0}, normEffError[nbins]={0}, 
				    normAccError[nbins]={0},
				      normFsrError[nbins]={0};
  double relUnfError[nbins]={0}, relEffError[nbins]={0}, relAccError[nbins]={0},
							   relFsrError[nbins]={0};
											       
  printf("Compute effect of unfolding\n");
  computeError(nObsDef,     nUnfDef,
	       nObsPlus,    nUnfPlus,
	       nObsMinus,   nUnfMinus,
	       normUnfError, relUnfError);

  printf("Compute effect of efficiency correction\n");
  computeError(nUnfDef,     nEffCorrDef,
	       nUnfPlus,    nEffCorrPlus,
	       nUnfMinus,   nEffCorrMinus,
	       normEffError, relEffError);

  printf("Compute effect of acceptance correction\n");
  computeError(nEffCorrDef,     nAccCorrDef,
	       nEffCorrPlus,    nAccCorrPlus,
	       nEffCorrMinus,   nAccCorrMinus,
	       normAccError, relAccError);

  printf("Compute effect of FSR correction\n");
  computeError(nAccCorrDef,     nFsrCorrDef,
	       nAccCorrPlus,    nFsrCorrPlus,
	       nAccCorrMinus,   nFsrCorrMinus,
	       normFsrError, relFsrError);

  printf("Compute total effect\n");
  double relTotError[nbins]={0}, normTotError[nbins]={0};
  computeError(nObsDef,     nFsrCorrDef,
	       nObsPlus,    nFsrCorrPlus,
	       nObsMinus,   nFsrCorrMinus,
	       normTotError, relTotError);
  

  
  // Summary table for relative errors
  printf("\nSummary table: relative errors\n");
  printf(" mass       Unf,%%      Eff, %%       Acc, %%       FSR, %%      Total, %%\n");
  for(int i=0; i<40; i++){
    printf("%4.0f-%4.0f  %6.1f     %6.1f      %6.1f      %6.1f      %6.1f\n", 
	   limits[i], limits[i+1],
	   relUnfError[i]*100, relEffError[i]*100, 
	   relAccError[i]*100, relFsrError[i]*100, relTotError[i]*100);

  }  

  
  // Summary table for normalized relative errors
  printf("\nSummary table: relative errors on the normalized to Z peak yields\n");
  printf(" mass       Unf,%%      Eff, %%       Acc, %%       FSR, %%      Total, %%\n");
  for(int i=0; i<40; i++){
    printf("%4.0f-%4.0f  %6.1f     %6.1f      %6.1f      %6.1f      %6.1f\n", 
	   limits[i], limits[i+1],
	   normUnfError[i]*100, normEffError[i]*100, 
	   normAccError[i]*100, normFsrError[i]*100, normTotError[i]*100);

  }  

  return;
}

void getArrays(TString dirName, TString fileName, double *nObs, double *nUnf,
	       double *nEffCorr, double *nAccCorr, double *nFsrCorr)
{

  // Compose the full file name and open the input file
  TString fullFileName = dirName + fileName;
  ifstream fin(fullFileName);
  printf("Open file %s\n", fullFileName.Data());
  
  // Read the content of the file into a vector of strings
  vector<TString> fileContent;
  string str;
  int nDataLines = 0;
  while( getline(fin,str) ){
    TString thisLine(str);
    // Discard all lines that do not contain the yield values
    // indicated by the \pm that happens only in those lines
    if( ! thisLine.Contains("pm") ) continue;
    nDataLines++;
    // The number of lines should never be above the total number
    // of mass bins.
    if( nDataLines > 40 ) assert(0);
    fileContent.push_back(thisLine);
  }
  
  // The layout of columns is:
  //         mass range | observed yields | unfolded | eff-corrected | acc-corrected | fsr-corrected

  // Get the yields out a column at a time
  getYieldsArray(fileContent, 1, nObs);
  getYieldsArray(fileContent, 2, nUnf);
  getYieldsArray(fileContent, 3, nEffCorr);
  getYieldsArray(fileContent, 4, nAccCorr);
  getYieldsArray(fileContent, 5, nFsrCorr);
  
}


void getYieldsArray(vector<TString> &fileContent, int ncol, double *data){

  // Now parse the strings and get the yields out
  for(int i=0; i<nbins; i++){
    TString thisString = fileContent.at(i);
    
    // First, split into substrings using & as a boundary
    TObjArray *columns = thisString.Tokenize("&");
    TObjString *tmpObjString;
    TObject *tmpObject;

    // Save the observed yields, column 1
    tmpObject = columns->At(ncol);
    tmpObjString = (TObjString*)tmpObject;
    TString obsYieldsString = tmpObjString->GetString();
    obsYieldsString.ReplaceAll("\\pm", "");
    obsYieldsString.ReplaceAll("$", "");
    // Now the string contains two numbers: the yield and the error.
    // Pull out the first number.
    TObjArray *yieldWithError = obsYieldsString.Tokenize(" ");
    tmpObject = yieldWithError->At(0);
    tmpObjString = (TObjString*)tmpObject;
    TString yieldsString = tmpObjString->GetString();
    data[i] = yieldsString.Atof();
  }

}

void checkObservedYieldsMatch(double *nObsDef, double *nObsPlus, double *nObsMinus){

  for(int i=0; i<nbins; i++){
    if( nObsDef[i] != nObsPlus[i])
      assert(0);
    if( nObsDef[i] != nObsMinus[i])
      assert(0);
  }
  printf("The observed yields match, as expected.\n");
  return;
}

void computeError(double *dataRefDef, double *dataCorrDef,
		  double *dataRefPlus, double *dataCorrPlus,
		  double *dataRefMinus, double *dataCorrMinus,
		  double *normError, double *relError)
{

  // Compute normalizations: add bins from 9 to 21
  double normRefDef = 0;
  double normRefPlus = 0;
  double normRefMinus = 0;
  double normCorrDef = 0;
  double normCorrPlus = 0;
  double normCorrMinus = 0;
  for(int i=9; i<=21; i++){
    normRefDef    += dataRefDef[i];
    normRefPlus   += dataRefPlus[i];
    normRefMinus  += dataRefMinus[i];
    normCorrDef   += dataCorrDef[i];
    normCorrPlus  += dataCorrPlus[i];
    normCorrMinus += dataCorrMinus[i];
  }
  
  printf("Effect of variation on normalization:\n");
  printf("   prev: plus = %f  minus = %f       next: plus %f  minus = %f\n",
	 normRefPlus/normRefDef, normRefMinus/normRefDef,
	 normCorrPlus/normCorrDef, normCorrMinus/normCorrDef);
	 

  // Compute all errors in this loop
  double diff1, diff2;
  for(int i=0; i<nbins; i++){

    // Compute relative errors
    double refDef    = dataRefDef[i];
    double corrDef   = dataCorrDef[i];

    if( refDef == 0 || corrDef == 0 ){
      printf("Denominator is zero!\n");
      relError[i] = 0;
      normError[i] = 0;
      continue;
    }

    diff1 = fabs(dataCorrDef[i] - dataCorrPlus[i]);
    diff2 = fabs(dataCorrDef[i] - dataCorrMinus[i]);
    double cumulThis = (diff1>diff2) ?  dataCorrPlus[i] : dataCorrMinus[i];
    double cumulPrev = (diff1>diff2) ?  dataRefPlus[i] : dataRefMinus[i];

    relError[i] = fabs(1.0-(cumulThis/corrDef) / (cumulPrev/refDef));

    // Compute normalized relative errors
    double relNormRefDef   = dataRefDef[i] / normRefDef;
    double relNormRefPlus  = dataRefPlus[i] / normRefPlus;
    double relNormRefMinus = dataRefMinus[i] / normRefMinus;
    double relNormCorrDef   = dataCorrDef[i] / normCorrDef;
    double relNormCorrPlus  = dataCorrPlus[i] / normCorrPlus;
    double relNormCorrMinus = dataCorrMinus[i] / normCorrMinus;

    diff1 = fabs( relNormCorrDef - relNormCorrPlus);
    diff2 = fabs( relNormCorrDef - relNormCorrMinus);
    cumulThis = (diff1>diff2) ?  relNormCorrPlus : relNormCorrMinus;
    cumulPrev = (diff1>diff2) ?  relNormRefPlus : relNormRefMinus;

    normError[i] = fabs(1.0-(cumulThis/relNormCorrDef) / (cumulPrev/relNormRefDef));

//     printf("%.0f-%.0f\n", limits[i], limits[i+1]);
//     printf("    the yields: Reference: def= %.0f   plus= %.0f  minus= %.0f  Corrected: def= %.0f  plus= %.0f minus= %.0f\n",
// 	   dataRefDef[i], dataRefPlus[i], dataRefMinus[i],
// 	   dataCorrDef[i], dataCorrPlus[i], dataCorrMinus[i]);
//     printf("    norms:      Reference: def= %.0f   plus= %.0f  minus= %.0f  Corrected: def= %.0f  plus= %.0f minus= %.0f\n",
// 	   normRefDef , normRefPlus , normRefMinus ,
// 	   normCorrDef , normCorrPlus , normCorrMinus );

  }

}


