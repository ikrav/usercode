#ifndef latexPrintouts_HH
#define latexPrintouts_HH

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system

#include <vector>                   // STL vector class

#include "../Include/DYTools.hh"        // helper class for organizing input ntuple files
     
#endif

void latexPrintoutAcceptance2D(TMatrixD accv, TMatrixD accErrv, TString producedBy);
void latexPrintoutAcceptance1D(TMatrixD accv, TMatrixD accErrv, TString producedBy);
void latexPrintoutEfficiency2D(TMatrixD effv, TMatrixD effErrv, TString producedBy);
void latexPrintoutEfficiency1D(TMatrixD effv, TMatrixD effErrv, TString producedBy);
void latexPrintoutScaleFactors2D(TMatrixD scalev, TMatrixD scaleErrv, TString producedBy);
void latexPrintoutScaleFactors1D(TMatrixD scalev, TMatrixD scaleErrv, TString producedBy);
void latexPrintoutFsr2D(TMatrixD corrv, TMatrixD corrErrv, TString producedBy);
void latexPrintoutFsr1D(TMatrixD corrv, TMatrixD corrErrv, TString producedBy);
void latexPrintoutFsrInAcceptance2D(TMatrixD corrv, TMatrixD corrErrv, TString producedBy);
void latexPrintoutFsrInAcceptance1D(TMatrixD corrv, TMatrixD corrErrv, TString producedBy);

void latexPrintoutOneValue2D(TMatrixD value, TMatrixD valueErr, TString producedBy, TString valueName, TString baseOfReferenceName, TString tableName);
void latexPrintoutOneValue1D(TMatrixD value, TMatrixD valueErr, TString producedBy, TString valueName, TString baseOfReferenceName, TString tableName);

void latexPrintoutTwoColumns2D(const int nValues, int* valuesType, TMatrixD** values, TMatrixD** valuesErr1, TMatrixD** valuesErr2, TString producedBy, TString* valuesName, TString baseOfReferenceName, TString tableName);
void latexPrintoutTwoColumns1D(const int nValues, int* valuesType, TMatrixD** values, TMatrixD** valuesErr1, TMatrixD** valuesErr2,TString producedBy, TString* valuesName, TString baseOfReferenceName, TString tableName);



#endif
