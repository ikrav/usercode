#include "TFile.h"
#include "TVectorD.h"
#include "TH1F.h"
#include "TObjString.h"
#include "TString.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TStyle.h"
#include "../Include/DYTools.hh"
#include "../Include/MyTools.hh"
#include "../Include/plotFunctions.hh"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

Int_t minMutualMultiple();
Int_t minMutualMultipleTwo(Int_t n1, Int_t n2);
Bool_t checkMatrixSize(TMatrixD m);

TString subtractBackground(const TString conf){


  // Read from configuration file only the location of the root files
  TString inputDir;
  Double_t lumi;
  Bool_t doWeight;
  ifstream ifs;
  ifs.open(conf.Data());
  assert(ifs.is_open());
  string line;
  while(getline(ifs,line)) {
    if(line[0]=='#') continue;
    stringstream ss1(line); ss1 >> lumi;
    getline(ifs,line);
    stringstream ss2(line); ss2 >> doWeight;
    getline(ifs,line);
    inputDir = TString(line);
    break;
  }
  ifs.close();
  inputDir.ReplaceAll("selected_events","yields");

  TFile file(inputDir+TString("/yields.root"));

  TVectorD* vec=(TVectorD *)file.Get("dummySampleCount");
  Int_t NSamples=vec->GetNoElements();
  TObjString* sample_names[NSamples];
  TMatrixD* yields[NSamples];
  TMatrixD* yieldsSumw2[NSamples];
  //since we have variable binning in rapidity, we have to determine the maximum number of Y-bins to create matrices
  int nYBinsMax=findMaxYBins();
  //fill out yields matrices from yields.root file

  TString sn[NSamples];

  for (int i=0; i<NSamples; i++)
    {
      TString nameName="sample_name_";
      nameName+=i;
      sample_names[i]=(TObjString *)file.Get(nameName);
      sn[i]=sample_names[i]->String();
      TString yieldsName="yields_";
      yieldsName+=sn[i];
      yields[i]=(TMatrixD *)file.Get(yieldsName);
      TString yieldsSumw2Name="yieldsSumw2_";
      yieldsSumw2Name+=sn[i];
      yieldsSumw2[i]=(TMatrixD *)file.Get(yieldsSumw2Name);

      if (!checkMatrixSize(*yields[i])) return "yields[i]: wrong size of matrix";
      if (!checkMatrixSize(*yieldsSumw2[i])) return "yieldsSumw2[i]: wrong size of matrix";
    }

  file.Close();




  bool useTrue2eBgDataDriven = false;
  bool useFakeBgDataDriven = false;
    
  printf("Calculate total backgrounds\n");
  if(useTrue2eBgDataDriven)
    printf("    use data driven estimates for true dielectron backgrounds\n");
  else
    printf("    use MC estimates for true dielectron backgrounds\n");
  if(useFakeBgDataDriven)
    printf("    use data driven estimates for fake dielectron backgrounds\n");
  else
    printf("    use MC estimates for fake dielectron backgrounds\n");


  TMatrixD observedYields(DYTools::nMassBins2D,nYBinsMax);
  TMatrixD observedYieldsErrorSquared(DYTools::nMassBins2D,nYBinsMax);
  TMatrixD signalYields(DYTools::nMassBins2D,nYBinsMax);
  TMatrixD signalYieldsError(DYTools::nMassBins2D,nYBinsMax);
  TMatrixD signalYieldsErrorSyst(DYTools::nMassBins2D,nYBinsMax);

  // Matrices to store backgrounds
  TMatrixD true2eBackground(DYTools::nMassBins2D,nYBinsMax);
  TMatrixD true2eBackgroundError(DYTools::nMassBins2D,nYBinsMax);
  TMatrixD true2eBackgroundErrorSyst(DYTools::nMassBins2D,nYBinsMax);
  TMatrixD wzzz(DYTools::nMassBins2D,nYBinsMax);
  TMatrixD wzzzError(DYTools::nMassBins2D,nYBinsMax);
  TMatrixD wzzzErrorSyst(DYTools::nMassBins2D,nYBinsMax);
  TMatrixD fakeEleBackground(DYTools::nMassBins2D,nYBinsMax);
  TMatrixD fakeEleBackgroundError(DYTools::nMassBins2D,nYBinsMax);
  TMatrixD fakeEleBackgroundErrorSyst(DYTools::nMassBins2D,nYBinsMax);


  true2eBackground = 0;
  true2eBackgroundError = 0;
  true2eBackgroundErrorSyst = 0;
  wzzz = 0;
  wzzzError = 0;
  wzzzErrorSyst = 0;
  fakeEleBackground = 0;
  fakeEleBackgroundError = 0;
  fakeEleBackgroundErrorSyst = 0;
  // The total background  
  TMatrixD totalBackground(DYTools::nMassBins2D,nYBinsMax);
  TMatrixD totalBackgroundError(DYTools::nMassBins2D,nYBinsMax);
  TMatrixD totalBackgroundErrorSyst(DYTools::nMassBins2D,nYBinsMax);


  // Calculate true dielectron background, which includes
  // WW, ttbar, and DY->tautau. By choice, we do not include WZ and ZZ.
  if(useTrue2eBgDataDriven)
    {
      TFile fTrueDataDriven(inputDir+TString("/true2eBkgDataPoints.root"));
      TMatrixD true2eBackgroundFromData        = *(TMatrixD*)fTrueDataDriven.Get("true2eBackgroundFromData");
      TMatrixD true2eBackgroundFromDataError     = *(TMatrixD*)fTrueDataDriven.Get("true2eBackgroundFromDataError");
      TMatrixD true2eBackgroundFromDataErrorSyst = *(TMatrixD*)fTrueDataDriven.Get("true2eBackgroundFromDataErrorSyst");
      if (!checkMatrixSize(true2eBackgroundFromData)) return "true2eBackground: wrong size of matrix";
      if (!checkMatrixSize(true2eBackgroundFromDataError)) return "true2eBackgroundError: wrong size of matrix";
      if (!checkMatrixSize(true2eBackgroundFromDataErrorSyst)) return "true2eBackgroundErrorSyst: wrong size of matrix";
      true2eBackground          = true2eBackgroundFromData;
      true2eBackgroundError     = true2eBackgroundFromDataError;
      true2eBackgroundErrorSyst = true2eBackgroundFromDataErrorSyst;
    }
  else
    {

       for(int i=0; i<DYTools::nMassBins2D; i++)
         for (int j=0; j<DYTools::nYBins[i]; j++)
           {
              true2eBackground(i,j)=0;
              true2eBackgroundError(i,j)=0;
              true2eBackgroundErrorSyst(i,j)=0;
           } 
         double mSyst=0;
         for (int k=0; k<NSamples; k++)
            {
              TMatrixD aux1(DYTools::nMassBins2D,nYBinsMax);
              TMatrixD aux2(DYTools::nMassBins2D,nYBinsMax);
              aux1=yields[k]->GetSub(0,DYTools::nMassBins2D-1,0,nYBinsMax-1);
              aux2=yieldsSumw2[k]->GetSub(0,DYTools::nMassBins2D-1,0,nYBinsMax-1);
              for (int i=0; i<DYTools::nMassBins2D; i++)
                for (int j=0; j<DYTools::nYBins[i]; j++)
                  {
                    if (sn[k]=="ttbar" || sn[k]=="ztt" || sn[k]=="ww")
                      {
  
                        true2eBackground(i,j)+=aux1(i,j);
                        true2eBackgroundError(i,j)+=aux2(i,j);
                        // Use ballpark numbers: 0% systematics on DY->tautau, 50% on ttbar, 100% on WW
                        if (sample_names[k]==(TObjString*)"ztt") mSyst=0.0;
                        if (sample_names[k]==(TObjString*)"ttbar") mSyst=0.5;
                        if (sample_names[k]==(TObjString*)"ww") mSyst=1.0;
                        true2eBackgroundErrorSyst(i,j)+=mSyst*mSyst*aux1(i,j)*aux1(i,j);

                      }
                   }
            }
  
    }

  // Calculate WZ and ZZ backgrounds
   for(int i=0; i<DYTools::nMassBins2D; i++)
     for (int j=0; j<DYTools::nYBins[i]; j++)
       {
          {
             wzzz(i,j)=0;
             wzzzError(i,j)=0;
             wzzzErrorSyst(i,j)=0;
          }
       } 
    for (int k=0; k<NSamples; k++)
      {
         if (sn[k]=="zz" || sn[k]=="wz")
           {
             TMatrixD aux1(DYTools::nMassBins2D,nYBinsMax);
             TMatrixD aux2(DYTools::nMassBins2D,nYBinsMax);
             aux1=yields[k]->GetSub(0,DYTools::nMassBins2D-1,0,nYBinsMax-1);
             aux2=yieldsSumw2[k]->GetSub(0,DYTools::nMassBins2D-1,0,nYBinsMax-1);
             for (int i=0; i<DYTools::nMassBins2D; i++)
               for (int j=0; j<DYTools::nYBins[i]; j++)
                 {
                     wzzz(i,j)+=aux1(i,j);
                     wzzzError(i,j)+=aux2(i,j);
                     wzzzErrorSyst(i,j)+=aux1(i,j)*aux1(i,j);
                 }
           }
       }


    // Calculate qcd and wjets backgrounds
    if(useFakeBgDataDriven)
      {
        TFile fFakeDataDriven(inputDir+TString("/fakeBkgDataPoints.root"));
        TMatrixD fakeEleBackgroundFromData          = *(TMatrixD*)fFakeDataDriven.Get("fakeEleBackgroundFromData");
        TMatrixD fakeEleBackgroundFromDataError     = *(TMatrixD*)fFakeDataDriven.Get("fakeEleBackgroundFromDataError");
        TMatrixD fakeEleBackgroundFromDataErrorSyst = *(TMatrixD*)fFakeDataDriven.Get("fakeEleBackgroundFromDataErrorSyst");
        if (!checkMatrixSize(fakeEleBackground)) return "fakeEleBackground: wrong size of matrix";
        if (!checkMatrixSize(fakeEleBackgroundError)) return "fakeEleBackgroundError: wrong size of matrix";
        if (!checkMatrixSize(fakeEleBackgroundErrorSyst)) return "fakeEleBackgroundErrorSyst: wrong size of matrix";
        fakeEleBackground          = fakeEleBackgroundFromData;
        fakeEleBackgroundError     = fakeEleBackgroundFromDataError;
        fakeEleBackgroundErrorSyst = fakeEleBackgroundFromDataErrorSyst;

       }
    else
      {
       for(int i=0; i<DYTools::nMassBins2D; i++)
         for (int j=0; j<DYTools::nYBins[i]; j++)
           {
              fakeEleBackground(i,j)=0;
              fakeEleBackgroundError(i,j)=0;
              fakeEleBackgroundErrorSyst(i,j)=0;
           } 

         for (int k=0; k<NSamples; k++)
           {
              if (sn[k]=="qcd" || sn[k]=="wjets")
                 {
                    TMatrixD aux1(DYTools::nMassBins2D,nYBinsMax);
                    TMatrixD aux2(DYTools::nMassBins2D,nYBinsMax);
                    aux1=yields[k]->GetSub(0,DYTools::nMassBins2D-1,0,nYBinsMax-1);
                    aux2=yieldsSumw2[k]->GetSub(0,DYTools::nMassBins2D-1,0,nYBinsMax-1);
                    for (int i=0; i<DYTools::nMassBins2D; i++)
                      for (int j=0; j<DYTools::nYBins[i]; j++)
                        {
                           fakeEleBackground(i,j)+=aux1(i,j);
                           fakeEleBackgroundError(i,j)+=aux2(i,j);
                           fakeEleBackgroundErrorSyst(i,j)+=0.5*0.5*aux1(i,j)*aux1(i,j);
                        }
                 }
            }
      }

    // Calculate the total background
    for (int i=0; i<DYTools::nMassBins2D; i++)
      for (int j=0; j<DYTools::nYBins[i]; j++)
        {
           totalBackground(i,j)=true2eBackground(i,j) + wzzz(i,j) + fakeEleBackground(i,j);
           totalBackgroundError(i,j)=sqrt( true2eBackgroundError(i,j) * true2eBackgroundError(i,j) +
				    wzzzError(i,j) * wzzzError(i,j) +
				    fakeEleBackgroundError(i,j) * fakeEleBackgroundError(i,j) );
           totalBackgroundErrorSyst(i,j)=sqrt( true2eBackgroundErrorSyst(i,j) * true2eBackgroundErrorSyst(i,j) +
				    wzzzErrorSyst(i,j) * wzzzErrorSyst(i,j) +
				    fakeEleBackgroundErrorSyst(i,j) * fakeEleBackgroundErrorSyst(i,j) );
         }



  // Loop over bins and perform background subtraction
  printf("Subtract background from observed yields\n");

  for (int k=0; k<NSamples; k++)
    {
      if (sn[k]=="data")
         {
            observedYields=*yields[k];
            observedYieldsErrorSquared=*yieldsSumw2[k];
              for (int i=0; i<DYTools::nMassBins2D; i++)
                for (int j=0; j<DYTools::nYBins[i]; j++)
                  {
                       signalYields(i,j) = observedYields(i,j) - totalBackground(i,j);
                       signalYieldsError(i,j) = sqrt( observedYieldsErrorSquared(i,j) + 
				 totalBackgroundError(i,j) * totalBackgroundError(i,j) );
                       signalYieldsErrorSyst(i,j) = totalBackgroundErrorSyst(i,j);
                  }
          }
    }

  TMatrixD bkgRatesUsual(DYTools::nMassBins2D,nYBinsMax);
  for (int i=0; i<DYTools::nMassBins2D; i++)
      for (int j=0; j<DYTools::nYBins[i]; j++)
           bkgRatesUsual(i,j)=100.0*totalBackground(i,j)/signalYields(i,j);     


  PlotMatrixVariousBinning(bkgRatesUsual,"bkgRatesPercent","COLZ");
 
  // Save sideband-subtracted signal yields
  TFile fileOut(inputDir+TString("/yields_bg-subtracted.root"),"recreate");
  signalYields         .Write("YieldsSignal");
  signalYieldsError    .Write("YieldsSignalErr");
  signalYieldsErrorSyst.Write("YieldsSignalSystErr");
  fileOut.Close();

  return "Ok";

}

Bool_t checkMatrixSize(TMatrixD m)
{  
  if ((m.GetNrows()==DYTools::nMassBins2D) && (m.GetNcols()==DYTools::findMaxYBins()))
  return 1;
  else return 0;
}
