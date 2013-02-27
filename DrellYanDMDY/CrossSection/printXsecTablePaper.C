// This script prints all r-shape cross sections in the exact format
// required for the paper (as of now EWK-11-007).

#include <TFile.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TString.h>


void printXsecTablePaper(){

  TFile f1("../root_files/DY_m10+pr+a05+o03+pr_4839pb_fsrUnfGood/xSec_results_1D.root");
  TMatrixD *xsecMatrixPreFsr = (TMatrixD*)f1.Get("normXSecByBin");
  TMatrixD *xsecErrMatrixPreFsr = (TMatrixD*)f1.Get("normXSecErrByBin");


  TFile f2("../root_files/DY_m10+pr+a05+o03+pr_4839pb_fsrUnfGood/xSecDET_results_1D.root");
  TMatrixD *xsecMatrixPreFsrDET = (TMatrixD*)f2.Get("normXSec");
  TMatrixD *xsecErrMatrixPreFsrDET = (TMatrixD*)f2.Get("normXSecErr");

  TFile f3("../root_files/DY_m10+pr+a05+o03+pr_4839pb_fsrUnfGood/xSecPostFsr_results_1DFull2011_hltEffOld.root");
  TMatrixD *xsecMatrixPostFsr = (TMatrixD*)f3.Get("normXSec");
  TMatrixD *xsecErrMatrixPostFsr = (TMatrixD*)f3.Get("normXSecErr");
  TMatrixD *xsecErrStatMatrixPostFsr = (TMatrixD*)f3.Get("normXSecErrStat");
  TMatrixD *xsecErrSystMatrixPostFsr = (TMatrixD*)f3.Get("normXSecErrSyst");

  TFile f4("../root_files/DY_m10+pr+a05+o03+pr_4839pb_fsrUnfGood/xSecPostFsrDET_results_1DFull2011_hltEffOld.root");
  TMatrixD *xsecMatrixPostFsrDET = (TMatrixD*)f4.Get("normXSec");
  TMatrixD *xsecErrStatMatrixPostFsrDET = (TMatrixD*)f4.Get("normXSecErrStat");
  TMatrixD *xsecErrSystMatrixPostFsrDET = (TMatrixD*)f4.Get("normXSecErrSyst");
  TVectorD *massBinLimits = (TVectorD*)f4.Get("massBinLimits");


  printf("%f\n", (*xsecMatrixPreFsr)(0,0));

  printf("loaded\n"); fflush(stdout);


  printf("Cross sections normalized to Z peak AND bin width\n");
  printf("           PreFsr                PostFsr               PreFsrDET         PostFsrDET\n");
  printf(" Mass bin [GeV] & $r^{i}_{pre FSR}$ & $r^{i}_{post FSR}$ & $r^{i}_{pre FSR, DET}$ & $r^{i}_{post FSR, DET}$ \\\\ \n");
  for(int i=0; i<40; i++){

    double binw = (*massBinLimits)[i+1] - (*massBinLimits)[i];

    double a1 = (*xsecErrStatMatrixPostFsr)(i,0);
    double a2 = (*xsecErrSystMatrixPostFsr)(i,0);
    double errPostFsr = sqrt(a1*a1+a2*a2);

    double e1 = (*xsecErrStatMatrixPostFsrDET)(i,0);
    double e2 = (*xsecErrSystMatrixPostFsrDET)(i,0);
    double errPostFsrDET = sqrt(e1*e1+e2*e2);

    TString format1 = " &$ %8.2f \\pm %8.2f $";
    if( i>=30 && i<35){
      format1 = " &$ %8.2f \\pm %8.2f $";
    }else if( i>=35 && i<38 ){
      format1 = " &$ %8.3f \\pm %8.3f $";
    }else if( i>= 38) {
      format1 = " &$ %8.4f \\pm %8.4f $";
    }
    TString format2 = format1;
    TString format3 = format1;
    TString format4 = format1;
    printf(" %4.0f-%4.0f  ", (*massBinLimits)(i), (*massBinLimits)(i+1));
    printf((format1+format2+format3+format4+TString(" \\\\ \\hline \n")).Data(), 
	   (*xsecMatrixPreFsr)(i,0)*100000, (*xsecErrMatrixPreFsr)(i,0)*100000,
	   (*xsecMatrixPostFsr)(i,0)*100000/binw, errPostFsr*100000/binw,
	   (*xsecMatrixPreFsrDET)(i,0)*100000/binw, (*xsecErrMatrixPreFsrDET)(i,0)*100000/binw,
	   (*xsecMatrixPostFsrDET)(i,0)*100000/binw, errPostFsrDET*100000/binw);
    
  }

  printf("\n\n Printout to provide inputs for combination with muons\n\n");
  printf("r-shape for the electron channel normalized to Z peak and to bin width, pre-FSR in full acceptance\n");
  printf(" Bin  MassMin MassMax         r-value     total error\n");
  for(int i=0; i<40; i++){
    printf(" %3d   %4.0f    %4.0f  ", i, (*massBinLimits)(i), (*massBinLimits)(i+1));
    printf("    %15.6g       %15.6g\n", (*xsecMatrixPreFsr)(i,0), (*xsecErrMatrixPreFsr)(i,0) );

  }


}
