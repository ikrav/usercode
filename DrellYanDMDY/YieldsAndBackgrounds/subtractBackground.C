#include "TFile.h"
#include "TVectorD.h"
#include "TH1F.h"
#include "TObjString.h"
#include "TString.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TStyle.h"
#include "../Include/DYTools.hh"
#include "../Include/DYToolsUI.hh"
#include "../Include/MyTools.hh"
#include "../Include/CPlot.hh"
#include "../Include/plotFunctions.hh"
#include "../Include/UnfoldingTools.hh"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

Int_t minMutualMultiple();
Int_t minMutualMultipleTwo(Int_t n1, Int_t n2);
Bool_t checkMatrixSize(TMatrixD m);

const int correct_error_code=1;
const int wzzzSystError_is_100percent=1;

TString subtractBackground(const TString conf,
		   DYTools::TSystematicsStudy_t runMode=DYTools::NORMAL,
			   const TString plotsDirExtraTag=""){


  std::cout << "\n\nRun mode: " << SystematicsStudyName(runMode) << "\n";
  switch(runMode) {
  case DYTools::NORMAL:
  case DYTools::ESCALE_STUDY:
  case DYTools::ESCALE_STUDY_RND:
    break;
  default:
    std::cout << "subtractBackground is not ready for runMode=" << SystematicsStudyName(runMode) << "\n";
    throw 2;
  }
 
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================

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

  // sOutDir is a static data member in the CPlot class.
  // There is a strange crash of the whole ROOT session well after
  // this script is executed when one attempts to exit ROOT, with 
  // a dump of memory map. This happens only on UNL Tier3, but
  // there is no sign of a problem on any other computer.
  //   The consequence of this variable is not set is that the plots
  // will be created in the local directory rather than the
  // one configured through sOutDir.
//   CPlot::sOutDir = outputDir + TString("/plots");

  if ((runMode==DYTools::ESCALE_STUDY) || (runMode==DYTools::ESCALE_STUDY_RND)) {
    CPlot::sOutDir = "plots_escale/";
  }
  else CPlot::sOutDir = "plots";
  CPlot::sOutDir += plotsDirExtraTag;


  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //============================================================================================================== 


  // yields.root
  TString yieldsFile=inputDir+TString("/yields") + DYTools::analysisTag + TString(".root");
  std::cout << "yieldsFile=<" << yieldsFile << ">\n";
  TFile file(yieldsFile);
  if (!file.IsOpen()) {
    std::cout << "failed to open a file <" << yieldsFile << ">\n";
    throw 2;
  }

  TVectorD* vec=(TVectorD *)file.Get("dummySampleCount");
  Int_t NSamples=vec->GetNoElements();
  TObjString* sample_names[NSamples];
  TMatrixD* yields[NSamples];
  TMatrixD* yieldsSumw2[NSamples];
  //since we have variable binning in rapidity, we have to determine the maximum number of Y-bins to create matrices
  int nYBinsMax=DYTools::findMaxYBins();
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


  TMatrixD observedYields(DYTools::nMassBins,nYBinsMax);
  TMatrixD observedYieldsErrorSquared(DYTools::nMassBins,nYBinsMax);
  TMatrixD signalYields(DYTools::nMassBins,nYBinsMax);
  TMatrixD signalYieldsError(DYTools::nMassBins,nYBinsMax);
  TMatrixD signalYieldsErrorSyst(DYTools::nMassBins,nYBinsMax);

  // Matrices to store backgrounds
  TMatrixD true2eBackground(DYTools::nMassBins,nYBinsMax);
  TMatrixD true2eBackgroundError(DYTools::nMassBins,nYBinsMax); // squared!
  TMatrixD true2eBackgroundErrorSyst(DYTools::nMassBins,nYBinsMax);
  TMatrixD wzzz(DYTools::nMassBins,nYBinsMax);
  TMatrixD wzzzError(DYTools::nMassBins,nYBinsMax);
  TMatrixD wzzzErrorSyst(DYTools::nMassBins,nYBinsMax);
  TMatrixD fakeEleBackground(DYTools::nMassBins,nYBinsMax);
  TMatrixD fakeEleBackgroundError(DYTools::nMassBins,nYBinsMax);
  TMatrixD fakeEleBackgroundErrorSyst(DYTools::nMassBins,nYBinsMax);


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
  TMatrixD totalBackground(DYTools::nMassBins,nYBinsMax);
  TMatrixD totalBackgroundError(DYTools::nMassBins,nYBinsMax);
  TMatrixD totalBackgroundErrorSyst(DYTools::nMassBins,nYBinsMax);


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

       for(int i=0; i<DYTools::nMassBins; i++)
         for (int j=0; j<DYTools::nYBins[i]; j++)
           {
              true2eBackground(i,j)=0;
              true2eBackgroundError(i,j)=0;
              true2eBackgroundErrorSyst(i,j)=0;
           } 
         double mSyst=0;
         for (int k=0; k<NSamples; k++)
            {
              TMatrixD aux1(DYTools::nMassBins,nYBinsMax);
              TMatrixD aux2(DYTools::nMassBins,nYBinsMax);
              aux1=yields[k]->GetSub(0,DYTools::nMassBins-1,0,nYBinsMax-1);
              aux2=yieldsSumw2[k]->GetSub(0,DYTools::nMassBins-1,0,nYBinsMax-1);
              for (int i=0; i<DYTools::nMassBins; i++)
                for (int j=0; j<DYTools::nYBins[i]; j++)
                  {
                    if (sn[k]=="ttbar" || sn[k]=="ztt" || sn[k]=="ww")
                      {
  
                        true2eBackground(i,j)+=aux1(i,j);
                        true2eBackgroundError(i,j)+=aux2(i,j);
                        // Use ballpark numbers: 0% systematics on DY->tautau, 50% on ttbar, 100% on WW
                        if (sn[k]=="ztt") mSyst=0.0;
                        if (sn[k]=="ttbar") mSyst=0.5;
                        if (sn[k]=="ww") mSyst=1.0;
                        true2eBackgroundErrorSyst(i,j)+=mSyst*mSyst*aux1(i,j)*aux1(i,j);
			//std::cout << sn[k] << ", " << sample_names[k]->String() << ": true2eBackgroundErrorSys(" << i << "," << j << ")+=(" << mSyst << '*' << aux1(i,j) << ")^2=" << (mSyst*mSyst*aux1(i,j)*aux1(i,j)) << "\n";

                      }
                   }
            }
  
    }

  // Calculate WZ and ZZ backgrounds, if true2eBackground is not available
  if (!useTrue2eBgDataDriven) {
    for (int k=0; k<NSamples; k++)
      {
         if (sn[k]=="zz" || sn[k]=="wz")
           {
             TMatrixD aux1(DYTools::nMassBins,nYBinsMax);
             TMatrixD aux2(DYTools::nMassBins,nYBinsMax);
             aux1=yields[k]->GetSub(0,DYTools::nMassBins-1,0,nYBinsMax-1);
             aux2=yieldsSumw2[k]->GetSub(0,DYTools::nMassBins-1,0,nYBinsMax-1);
             for (int i=0; i<DYTools::nMassBins; i++)
               for (int j=0; j<DYTools::nYBins[i]; j++)
                 {
                     wzzz(i,j)+=aux1(i,j);
                     wzzzError(i,j)+=aux2(i,j);
		     if (wzzzSystError_is_100percent)
		       wzzzErrorSyst(i,j)+=aux1(i,j);
		     else 
		       wzzzErrorSyst(i,j)+=aux1(i,j)*aux1(i,j);

		     //std::cout << "wzzzErrorSyst(" << i << "," << j << ")+=" <<  aux1(i,j) << "^2=" << (aux1(i,j)*aux1(i,j)) << "\n";
                 }
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
       for(int i=0; i<DYTools::nMassBins; i++)
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
                    TMatrixD aux1(DYTools::nMassBins,nYBinsMax);
                    TMatrixD aux2(DYTools::nMassBins,nYBinsMax);
                    aux1=yields[k]->GetSub(0,DYTools::nMassBins-1,0,nYBinsMax-1);
                    aux2=yieldsSumw2[k]->GetSub(0,DYTools::nMassBins-1,0,nYBinsMax-1);
                    for (int i=0; i<DYTools::nMassBins; i++)
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
    for (int i=0; i<DYTools::nMassBins; i++)
      for (int j=0; j<DYTools::nYBins[i]; j++)
        {
           totalBackground(i,j)=true2eBackground(i,j) + wzzz(i,j) + fakeEleBackground(i,j);
	   if (correct_error_code) {
           totalBackgroundError(i,j)=sqrt( true2eBackgroundError(i,j) +
					   wzzzError(i,j) +
					   fakeEleBackgroundError(i,j) );
	   if (wzzzSystError_is_100percent) {
           totalBackgroundErrorSyst(i,j)=sqrt( true2eBackgroundErrorSyst(i,j) +
				       wzzzErrorSyst(i,j)*wzzzErrorSyst(i,j) +
					       fakeEleBackgroundErrorSyst(i,j) );
	   }
	   else {
           totalBackgroundErrorSyst(i,j)=sqrt( true2eBackgroundErrorSyst(i,j) +
					       wzzzErrorSyst(i,j) +
					       fakeEleBackgroundErrorSyst(i,j) );
	   }
	   }
	   else {
           totalBackgroundError(i,j)=sqrt( true2eBackgroundError(i,j) * true2eBackgroundError(i,j) +
				    wzzzError(i,j) * wzzzError(i,j) +
				    fakeEleBackgroundError(i,j) * fakeEleBackgroundError(i,j) );
           totalBackgroundErrorSyst(i,j)=sqrt( true2eBackgroundErrorSyst(i,j) * true2eBackgroundErrorSyst(i,j) +
				    wzzzErrorSyst(i,j) * wzzzErrorSyst(i,j) +
				    fakeEleBackgroundErrorSyst(i,j) * fakeEleBackgroundErrorSyst(i,j) );
	   }
         }



  // Loop over bins and perform background subtraction
  printf("Subtract background from observed yields\n");

  for (int k=0; k<NSamples; k++)
    {
      if (sn[k]=="data")
         {
            observedYields=*yields[k];
            observedYieldsErrorSquared=*yieldsSumw2[k];
              for (int i=0; i<DYTools::nMassBins; i++)
                for (int j=0; j<DYTools::nYBins[i]; j++)
                  {
                       signalYields(i,j) = observedYields(i,j) - totalBackground(i,j);
                       signalYieldsError(i,j) = sqrt( observedYieldsErrorSquared(i,j) + 
				 totalBackgroundError(i,j) * totalBackgroundError(i,j) );
                       signalYieldsErrorSyst(i,j) = totalBackgroundErrorSyst(i,j);
                  }
          }
    }

  TMatrixD bkgRatesUsual(DYTools::nMassBins,nYBinsMax);
  for (int i=0; i<DYTools::nMassBins; i++)
      for (int j=0; j<DYTools::nYBins[i]; j++)
           bkgRatesUsual(i,j)=100.0*totalBackground(i,j)/signalYields(i,j);


  //
  // data-MC residual difference: 
  // create weight factors to make MC prediction signal-like
  //
  TMatrixD zeePredictedYield(DYTools::nMassBins,nYBinsMax);
  int zeeFound=0;
  for (int i=0; i<NSamples; ++i) {
    if (sn[i] == "zee") {
      zeePredictedYield= *yields[i];
      zeeFound=1;
      break;
    }
  }
  if (!zeeFound) {
    std::cout << "failed to locate MC zee sample\n";
    throw 2;
  }

  // 1. determine counts in Z-peak area in data and MC
  double countDataNoBkg=0., countMCsignal=0.;
  for (int i=0; i<DYTools::nMassBins; i++) {
    if ((DYTools::massBinLimits[i]>59.999) && (DYTools::massBinLimits[i]<120.0001)) {
      for (int yi=0; yi<DYTools::nYBins[i]; ++yi) {
	countDataNoBkg += signalYields[i][yi];
      }
      for (int yi=0; yi<DYTools::nYBins[i]; ++yi) {
	countMCsignal += zeePredictedYield[i][yi];
      }
    }
  }
  if (fabs(countDataNoBkg)<1e-7) {
    std::cout << "found no data events in Z-peak region\n";
    throw 2;
  }
  if (fabs(countMCsignal)<1e-7) {
    std::cout << "found no data events in Z-peak region\n";
    throw 2;
  }
  // 2. create weight maps
  TMatrixD zeeMCShapeReweight(DYTools::nMassBins,nYBinsMax);
  zeeMCShapeReweight=0;
  double ZpeakFactor = countMCsignal/double(countDataNoBkg);
  for (int i=0; i<DYTools::nMassBins; i++) {
    for (int yi=0; yi<DYTools::nYBins[i]; ++yi) {
      double mc = zeePredictedYield[i][yi];
      double weight=(fabs(mc)<1e-7) ?  
	0. :  ( ZpeakFactor * signalYields[i][yi] / mc );
      zeeMCShapeReweight[i][yi] = weight;
    }
  }
 

 
  // Save sideband-subtracted signal yields
  // inputDir+TString("/yields_bg-subtracted.root")
  TString outFileName=inputDir + TString("/yields_bg-subtracted") + DYTools::analysisTag + TString(".root");
  TFile fileOut(outFileName,"recreate");
  signalYields         .Write("YieldsSignal");
  signalYieldsError    .Write("YieldsSignalErr");   // not squared
  signalYieldsErrorSyst.Write("YieldsSignalSystErr"); // not squared
  zeeMCShapeReweight   .Write("ZeeMCShapeReweight");
  unfolding::writeBinningArrays(fileOut);
  fileOut.Close();

  TString outFileNamePlots=outFileName;
  outFileNamePlots.Replace(outFileNamePlots.Index(".root"),sizeof(".root"),"-plots.root");
  TFile *foutPlots=new TFile(outFileNamePlots,"recreate");
  PlotMatrixVariousBinning(bkgRatesUsual,"bkgRatesPercent","COLZ",foutPlots);
  PlotMatrixVariousBinning(signalYields,"signalYields","COLZ",foutPlots);
  PlotMatrixVariousBinning(signalYieldsError,"signalYieldsError","COLZ",foutPlots);
  foutPlots->Close();

  if (0) {
    int yi=0;
    // Print tables of yields and background
    printf("\n\n\t\tTables for yi=%d\n\n",yi);
    
    // Table 1: split background into true, wz/zz, and qcd
    printf(" Note: stat error in signal yield contain stat error on background,\n");
    printf("   and syst error on signal yield contains syst error on background\n");
    printf("mass range            observed       true2e-bg         wz-zz-bg               fake-bg                 total-bg            signal\n");
    for(int i=0; i<DYTools::nMassBins; i++){
      printf("%5.1f-%5.1f GeV: ", DYTools::massBinLimits[i], DYTools::massBinLimits[i+1]);
      printf(" %7.0f+-%3.0f ", observedYields[i][yi], sqrt(observedYieldsErrorSquared[i][yi]));
      printf(" %5.1f+-%4.1f+-%4.1f ", true2eBackground[i][yi], sqrt(true2eBackgroundError[i][yi]), sqrt(true2eBackgroundErrorSyst[i][yi]));
      printf(" %6.2f+-%4.2f+-%4.2f ", wzzz[i][yi], sqrt(wzzzError[i][yi]), sqrt(wzzzErrorSyst[i][yi]));
      printf(" %5.1f+-%5.1f+-%5.1f ", fakeEleBackground[i][yi], sqrt(fakeEleBackgroundError[i][yi]), sqrt(fakeEleBackgroundErrorSyst[i][yi]));
      printf("    %5.1f+-%4.1f+-%4.1f ", totalBackground[i][yi], totalBackgroundError[i][yi], totalBackgroundErrorSyst[i][yi]);
      printf("    %8.1f+-%5.1f+-%5.1f ", signalYields[i][yi], signalYieldsError[i][yi], signalYieldsErrorSyst[i][yi]);
      printf("\n");
    }

    // Table 2: combined true2e and WZ/ZZ backgrounds only
    printf("\n  only true2e-bg + ww-wz\n");
    printf("mass range      true2e, includingwz/zz\n");
    for(int i=0; i<DYTools::nMassBins; i++){
      printf("%5.1f-%5.1f GeV: ", DYTools::massBinLimits[i], DYTools::massBinLimits[i+1]);
      double val = true2eBackground[i][yi] + wzzz[i][yi];
      double err=0, sys=0;
      if (correct_error_code){
      err = sqrt(true2eBackgroundError[i][yi]
		 + wzzzError[i][yi]);
      if (wzzzSystError_is_100percent) {
      sys = sqrt(true2eBackgroundErrorSyst[i][yi]
		 + wzzzErrorSyst[i][yi] * wzzzErrorSyst[i][yi]);
      }
      else {
      sys = sqrt(true2eBackgroundErrorSyst[i][yi]
		 + wzzzErrorSyst[i][yi]);
      }
      }
      else {
      err = sqrt(true2eBackgroundError[i][yi]*true2eBackgroundError[i][yi]
			+ wzzzError[i][yi]*wzzzError[i][yi]);
      sys = sqrt(true2eBackgroundErrorSyst[i][yi]*true2eBackgroundErrorSyst[i][yi]
			+ wzzzErrorSyst[i][yi]*wzzzErrorSyst[i][yi]);
      }
      printf(" %5.1f+-%4.1f+-%4.1f ", val,err, sys);
      printf("\n");
    }

    // Table 3: Systematic error on signal yields assuming that it includes
    // only the syst. error on the background.
    printf("\n  Systematics, %% relative to background subtracted yields\n");
    printf("mass range            subtr-signal    total-bg      syst-from-bg-frac      syst-from-bg-percent\n");
    for(int i=0; i<DYTools::nMassBins; i++){
      printf("%5.1f-%5.1f GeV: ", DYTools::massBinLimits[i], DYTools::massBinLimits[i+1]);
      printf("    %8.1f+-%5.1f+-%4.1f ", signalYields[i][yi], signalYieldsError[i][yi],signalYieldsErrorSyst[i][yi]);
      printf("    %5.1f+-%4.1f+-%4.1f ", totalBackground[i][yi], totalBackgroundError[i][yi], totalBackgroundErrorSyst[i][yi]);
      printf("    %6.4f ", totalBackgroundErrorSyst[i][yi]/signalYields[i][yi]);
      printf("    %6.1f ", totalBackgroundErrorSyst[i][yi]*100.0/signalYields[i][yi]);
      printf("\n");
    }
  }
 
  return "Ok";

}

Bool_t checkMatrixSize(TMatrixD m)
{  
  if ((m.GetNrows()==DYTools::nMassBins) && (m.GetNcols()==DYTools::findMaxYBins()))
  return 1;
  else return 0;
}
