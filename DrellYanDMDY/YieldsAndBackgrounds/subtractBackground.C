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
#include "../Include/MitStyleRemix.hh"
#include "../Include/ComparisonPlot.hh"
#include "../Include/plotFunctions.hh"
#include "../Include/UnfoldingTools.hh"
#include "../Include/latexPrintouts.hh"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

Int_t minMutualMultiple();
Int_t minMutualMultipleTwo(Int_t n1, Int_t n2);
Bool_t checkMatrixSize(TMatrixD m);

void unsquareMatrixElements(TMatrixD &m);

template<class T>
inline T SQR(const T &x) { return (x)*(x); }

// -----------------------------------------------------------------------------

TString subtractBackground(const TString conf,
		   DYTools::TSystematicsStudy_t runMode=DYTools::NORMAL,
			   const TString plotsDirExtraTag="",
			   int performPUReweight=1){


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
  Bool_t hasData=false;

  ifstream ifs;
  ifs.open(conf.Data());
  assert(ifs.is_open());
  string line;
  int state=0;
  while(getline(ifs,line)) {
    if(line[0]=='#') continue;
    if(line[0]=='%') { 
      state++; 
      continue; 
    }
    if(state==0) {
      stringstream ss1(line); ss1 >> lumi;
      getline(ifs,line);
      stringstream ss2(line); ss2 >> doWeight;
      getline(ifs,line);
      inputDir = TString(line);
      TString escaleTag_loc,format_loc;
      getline(ifs,line);
      stringstream ss3(line); ss3 >> escaleTag_loc;
      getline(ifs,line);
      format_loc = TString(line);
      state++;
    }
    else if (state==1) {
      hasData=true;
    }
    else break;
  }
  ifs.close();
  if (performPUReweight) inputDir.ReplaceAll("selected_events","yields");
  else inputDir.ReplaceAll("selected_events","yields_noPU");

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
  else CPlot::sOutDir = TString("plots") + DYTools::analysisTag;
  if (!performPUReweight) CPlot::sOutDir.Append("_noPU");
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
  TMatrixD observedYieldsErr(DYTools::nMassBins,nYBinsMax);
  TMatrixD observedYieldsErrorSquared(DYTools::nMassBins,nYBinsMax);
  TMatrixD signalYields(DYTools::nMassBins,nYBinsMax);
  TMatrixD signalYieldsError(DYTools::nMassBins,nYBinsMax);
  TMatrixD signalYieldsErrorSyst(DYTools::nMassBins,nYBinsMax);

  // Matrices to store backgrounds
  TMatrixD true2eBackground(DYTools::nMassBins,nYBinsMax);
  TMatrixD true2eBackgroundError(DYTools::nMassBins,nYBinsMax);
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
  // WW, ttbar, Wt, and DY->tautau. By choice, we do not include WZ and ZZ.
  if(useTrue2eBgDataDriven)
    {
      TString true2eFName=inputDir+TString("/true2eBkgDataPoints_");
      if (DYTools::study2D) true2eFName.Append("2D");
      else true2eFName.Append("1D");
      true2eFName.Append(".root");
      TFile fTrueDataDriven(true2eFName);
      assert(unfolding::checkBinningArrays(fTrueDataDriven));
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

      double mSyst=0;
      for (int k=0; k<NSamples; k++) {
	if (sn[k]=="ttbar" || sn[k]=="ztt" || sn[k]=="ww" || sn[k]=="wtop") {
	  TMatrixD aux1= (* yields[k] );
	  TMatrixD aux2= (* yieldsSumw2[k] );
	  // Use ballpark numbers: 0% systematics on DY->tautau, 50% on ttbar, 100% on WW
	  if (sn[k]=="ztt") mSyst=0.0;
	  else if (sn[k]=="ttbar") mSyst=0.5;
	  else if ((sn[k]=="ww") || (sn[k]=="wtop")) mSyst=1.0;
	  for (int i=0; i<DYTools::nMassBins; i++)
	    for (int j=0; j<DYTools::nYBins[i]; j++) {
	      true2eBackground(i,j)+=aux1(i,j);
	      true2eBackgroundError(i,j)+=aux2(i,j);
	      true2eBackgroundErrorSyst(i,j)+= mSyst*mSyst*aux1(i,j)*aux1(i,j);
	      //std::cout << sn[k] << ", " << sample_names[k]->String() << ": true2eBackgroundErrorSys(" << i << "," << j << ")+=(" << mSyst << '*' << aux1(i,j) << ")^2=" << (mSyst*mSyst*aux1(i,j)*aux1(i,j)) << "\n";
	      
	    }
	}
      }
      // extract root to have unsquared error
      unsquareMatrixElements(true2eBackgroundError);
      unsquareMatrixElements(true2eBackgroundErrorSyst);
    }

  // Calculate WZ and ZZ backgrounds, if true2eBackground is not available
  if (!useTrue2eBgDataDriven) {
    for (int k=0; k<NSamples; k++) {
      if (sn[k]=="zz" || sn[k]=="wz") {
	TMatrixD aux1= (* yields[k] );
	TMatrixD aux2= (* yieldsSumw2[k] );
	for (int i=0; i<DYTools::nMassBins; i++)
	  for (int j=0; j<DYTools::nYBins[i]; j++) {
	    wzzz(i,j)+=aux1(i,j);
	    wzzzError(i,j)+=aux2(i,j);
	    wzzzErrorSyst(i,j)+=aux1(i,j);
	    //std::cout << "wzzzErrorSyst(" << i << "," << j << ")+=" <<  aux1(i,j) << "^2=" << (aux1(i,j)*aux1(i,j)) << "\n";
	  }
      }
    }
    // extract root to have unsquared error
    unsquareMatrixElements(wzzzError);
    // make 100% error //unsquareMatrixElements(wzzzErrorSyst);
  }


  // Calculate qcd and wjets backgrounds
  if(useFakeBgDataDriven)
      {
        TString fakeEEFName(inputDir+TString("/fakeBkgDataPoints_"));
	if (DYTools::study2D) fakeEEFName.Append("2D");
	else fakeEEFName.Append("1D");
	fakeEEFName.Append(".root");
        TFile fFakeDataDriven(fakeEEFName);
	assert(unfolding::checkBinningArrays(fFakeDataDriven));
        TMatrixD fakeEleBackgroundFromData          = *(TMatrixD*)fFakeDataDriven.Get("fakeBackgroundFromData");
        TMatrixD fakeEleBackgroundFromDataError     = *(TMatrixD*)fFakeDataDriven.Get("fakeBackgroundFromDataError");
        TMatrixD fakeEleBackgroundFromDataErrorSyst = *(TMatrixD*)fFakeDataDriven.Get("fakeBackgroundFromDataErrorSyst");
        if (!checkMatrixSize(fakeEleBackgroundFromData)) return "fakeEleBackgroundFromData: wrong size of matrix";
        if (!checkMatrixSize(fakeEleBackgroundFromDataError)) return "fakeEleBackgroundFromDataError: wrong size of matrix";
        if (!checkMatrixSize(fakeEleBackgroundFromDataErrorSyst)) return "fakeEleBackgroundFromDataErrorSyst: wrong size of matrix";
        fakeEleBackground          = fakeEleBackgroundFromData;
        fakeEleBackgroundError     = fakeEleBackgroundFromDataError;
        fakeEleBackgroundErrorSyst = fakeEleBackgroundFromDataErrorSyst;

       }
    else
      {

	for (int k=0; k<NSamples; k++) {
	  if (sn[k]=="qcd" || sn[k]=="wjets") {
	    TMatrixD aux1= (* yields[k]);
	    TMatrixD aux2= (* yieldsSumw2[k]);
	    for (int i=0; i<DYTools::nMassBins; i++)
	      for (int j=0; j<DYTools::nYBins[i]; j++) {
		fakeEleBackground(i,j)+=aux1(i,j);
		fakeEleBackgroundError(i,j)+=aux2(i,j);
		fakeEleBackgroundErrorSyst(i,j)+=0.5*0.5*aux1(i,j)*aux1(i,j);
	      }
	  }
	}
	// extract root to have unsquared error
	unsquareMatrixElements(fakeEleBackgroundError);
	unsquareMatrixElements(fakeEleBackgroundErrorSyst);
      }

  // total true2e background
  TMatrixD totalTrue2eBackground=true2eBackground;
  TMatrixD totalTrue2eBackgroundError=true2eBackgroundError;
  TMatrixD totalTrue2eBackgroundErrorSyst=true2eBackgroundErrorSyst;
  if (!useTrue2eBgDataDriven) {
    // add wzzz contribution
    for (int i=0; i<DYTools::nMassBins; ++i) {
      for (int j=0; j<DYTools::nYBins[i]; ++j) {
	totalTrue2eBackground(i,j)= true2eBackground(i,j) + wzzz(i,j);
	totalTrue2eBackgroundError(i,j)= sqrt( SQR(true2eBackgroundError(i,j)) +
					       SQR(wzzzError(i,j)) );
	totalTrue2eBackgroundErrorSyst(i,j)=
	  sqrt( SQR(true2eBackgroundErrorSyst(i,j)) +
		SQR(wzzzErrorSyst(i,j)) );
      }
    }
  }

    // Calculate the total background
    for (int i=0; i<DYTools::nMassBins; i++)
      for (int j=0; j<DYTools::nYBins[i]; j++) 
        {
           totalBackground(i,j)=true2eBackground(i,j) + wzzz(i,j) + fakeEleBackground(i,j);
           totalBackgroundError(i,j)=sqrt( SQR(true2eBackgroundError(i,j)) +
					   SQR(wzzzError(i,j)) +
					   SQR(fakeEleBackgroundError(i,j)) );
	   totalBackgroundErrorSyst(i,j)=sqrt( SQR(true2eBackgroundErrorSyst(i,j)) +
					       SQR(wzzzErrorSyst(i,j)) +
					       SQR(fakeEleBackgroundErrorSyst(i,j)) );
	}



  // Loop over bins and perform background subtraction
  printf("Subtract background from observed yields\n");

  for (int k=0; k<NSamples; k++) {
    if (sn[k]=="data") {
      observedYields=*yields[k];
      observedYieldsErrorSquared=*yieldsSumw2[k];
      observedYieldsErr=0;
      for (int i=0; i<DYTools::nMassBins; i++)
	for (int j=0; j<DYTools::nYBins[i]; j++) {
	  double sy=observedYields(i,j) - totalBackground(i,j);
	  signalYields(i,j) = (sy>0) ? sy : 0.;
	  signalYieldsError(i,j) = 
	    sqrt( observedYieldsErrorSquared(i,j) + 
		  SQR(totalBackgroundError(i,j)) );
	  signalYieldsErrorSyst(i,j) = totalBackgroundErrorSyst(i,j);
	  observedYieldsErr(i,j) = sqrt(observedYieldsErrorSquared(i,j));
	}
    }
  }

  TMatrixD bkgRatesUsual(DYTools::nMassBins,nYBinsMax);
  for (int i=0; i<DYTools::nMassBins; i++) {
    for (int j=0; j<DYTools::nYBins[i]; j++) {
      bkgRatesUsual(i,j)= (signalYields(i,j)>0.) ?
	100.0*totalBackground(i,j)/signalYields(i,j) : 0.;
    }
  }


  //
  // data-MC residual difference: 
  // create weight factors to make MC prediction signal-like
  //
  TMatrixD zeePredictedYield(DYTools::nMassBins,nYBinsMax);
  TMatrixD zeePredictedYieldErr(DYTools::nMassBins,nYBinsMax);
  int zeeFound=0;
  for (int i=0; i<NSamples; ++i) {
    if (sn[i] == "zee") {
      zeePredictedYield= *yields[i];
      zeePredictedYieldErr=0;
      zeePredictedYieldErr= *yieldsSumw2[i];
      unsquareMatrixElements(zeePredictedYieldErr);
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
    if ((DYTools::massBinLimits[i]>59.999) && (DYTools::massBinLimits[i+1]<120.0001)) {
      for (int yi=0; yi<DYTools::nYBins[i]; ++yi) {
	countDataNoBkg += signalYields[i][yi];
      }
      for (int yi=0; yi<DYTools::nYBins[i]; ++yi) {
	countMCsignal += zeePredictedYield[i][yi];
      }
    }
  }
  if (fabs(countDataNoBkg)<1e-7) {
    std::cout << "found no data events in Z-peak region. Aborting\n";
    assert(0);
  }
  if (fabs(countMCsignal)<1e-7) {
    std::cout << "found no mc events in Z-peak region. Aborting\n";
    assert(0);
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
  /*
  zeePredictedYield. Write("mcYieldsSignal");
  zeePredictedYieldErr.Write("mcYieldsSignalErr");  // not squared
  zeePredictedYield *= (1/ZpeakFactor);
  zeePredictedYieldErr *= (1/ZpeakFactor);
  zeePredictedYield.Write("mcYieldsSignalScaled");
  zeePredictedYieldErr.Write("mcYieldsSignalScaledErr");
  observedYields.Write("observedYields");
  observedYieldsErr.Write("observedYieldsErr"); // not squared
  totalBackground.Write("totalBackground");
  totalBackgroundError.Write("totalBackgroundErr");
  totalBackgroundErrorSyst.Write("totalBackgroundErrSyst");
  totalTrue2eBackground.Write("true2eBkgr");
  totalTrue2eBackgroundError.Write("true2eBkgrErr");
  totalTrue2eBackgroundErrorSyst.Write("true2eBkgrSystErr");
  fakeEleBackground.Write("fake2eBkgr");
  fakeEleBackgroundError.Write("fake2eBkgrErr");
  fakeEleBackgroundErrorSyst.Write("fake2eBkgrSystErr");
  */
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
    for(int i=0; i<DYTools::nMassBins; i++){
      // Print tables of yields and background
      if ( (DYTools::study2D==1) ||
	   ((DYTools::study2D==0) && (i==0)) ) {
	   
	printf("\n\n\t\tTables for iMass=%d\n\n",i);
	
	// Table 1: split background into true, wz/zz, and qcd
	printf(" Note: stat error in signal yield contain stat error on background,\n");
	printf("   and syst error on signal yield contains syst error on background\n");
	printf("mass range   rapidity range   observed       true2e-bg         wz-zz-bg               fake-bg                 total-bg            signal\n");
      }
      for (int yi=0; yi<DYTools::nYBins[i]; ++yi) {
	double absYmin=0, absYmax=0;
	DYTools::findAbsYValueRange(i,yi, absYmin,absYmax);
	printf("%5.1f-%5.1f GeV ", DYTools::massBinLimits[i], DYTools::massBinLimits[i+1]);
	printf("%4.2f-%4.2f ", absYmin,absYmax);
	printf(" %7.0f+-%3.0f ", observedYields[i][yi], observedYieldsErr[i][yi]);
	printf(" %5.1f+-%4.1f+-%4.1f ", true2eBackground[i][yi], true2eBackgroundError[i][yi], true2eBackgroundErrorSyst[i][yi]);
	printf(" %6.2f+-%4.2f+-%4.2f ", wzzz[i][yi], wzzzError[i][yi], wzzzErrorSyst[i][yi]);
	printf(" %5.1f+-%5.1f+-%5.1f ", fakeEleBackground[i][yi], fakeEleBackgroundError[i][yi], fakeEleBackgroundErrorSyst[i][yi]);
	printf("    %5.1f+-%4.1f+-%4.1f ", totalBackground[i][yi], totalBackgroundError[i][yi], totalBackgroundErrorSyst[i][yi]);
	printf("    %8.1f+-%5.1f+-%5.1f ", signalYields[i][yi], signalYieldsError[i][yi], signalYieldsErrorSyst[i][yi]);
	printf("\n");
      }
    }
    std::cout << std::endl;
  }

  if (1) {
    for(int i=0; i<DYTools::nMassBins; i++){
      // Print tables of yields and background
      if ( (DYTools::study2D==1) ||
	   ((DYTools::study2D==0) && (i==0)) ) {
	   
	printf("\n\n\t\tTables for iMass=%d\n\n",i);
	
	// Table 1: split background into true, wz/zz, and qcd
	printf(" Note: stat error in signal yield contain stat error on background,\n");
	printf("   and syst error on signal yield contains syst error on background\n");
	printf("mass range   rapidity range   observed       true2e-bg         fake-bg                 total-bg            signal\n");
      }
      for (int yi=0; yi<DYTools::nYBins[i]; ++yi) {
	double absYmin=0, absYmax=0;
	DYTools::findAbsYValueRange(i,yi, absYmin,absYmax);
	printf("%5.1f-%5.1f GeV ", DYTools::massBinLimits[i], DYTools::massBinLimits[i+1]);
	printf("%4.2f-%4.2f ", absYmin,absYmax);
	printf(" %7.0f+-%3.0f ", observedYields[i][yi], observedYieldsErr[i][yi]);
	printf(" %5.1f+-%4.1f+-%4.1f ", true2eBackground[i][yi], true2eBackgroundError[i][yi], true2eBackgroundErrorSyst[i][yi]);
	printf(" %5.1f+-%5.1f+-%5.1f ", fakeEleBackground[i][yi], fakeEleBackgroundError[i][yi], fakeEleBackgroundErrorSyst[i][yi]);
	printf("    %5.1f+-%4.1f+-%4.1f ", totalBackground[i][yi], totalBackgroundError[i][yi], totalBackgroundErrorSyst[i][yi]);
	printf("    %8.1f+-%5.1f+-%5.1f ", signalYields[i][yi], signalYieldsError[i][yi], signalYieldsErrorSyst[i][yi]);
	printf("\n");
      }
    }
    std::cout << std::endl;
  }


  if (0) {
    int yi=0;
    // Table 2: combined true2e and WZ/ZZ backgrounds only
    printf("\n  only true2e-bg + ww-wz\n");
    printf("mass range      true2e, includingwz/zz\n");
    for(int i=0; i<DYTools::nMassBins; i++){
      printf("%5.1f-%5.1f GeV: ", DYTools::massBinLimits[i], DYTools::massBinLimits[i+1]);
      double val = true2eBackground[i][yi] + wzzz[i][yi];
      double err=0, sys=0;
      err = sqrt(SQR(true2eBackgroundError[i][yi])
		 + SQR(wzzzError[i][yi]));
      sys = sqrt(SQR(true2eBackgroundErrorSyst[i][yi])
		 + SQR(wzzzErrorSyst[i][yi]) );
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

//Latex printout
  if (DYTools::study2D==1)
     latexPrintoutBackgroundRates2D(observedYields, observedYieldsErr, 
                                  totalBackground, totalBackgroundError, 
                                  totalBackgroundErrorSyst, bkgRatesUsual, 
                                    "YieldsAndBackgrounds/subtractBackground.C");
   else if (DYTools::study2D==0)
     latexPrintoutBackgroundRates1D(observedYields, observedYieldsErr, 
                                  totalBackground, totalBackgroundError, 
                                  totalBackgroundErrorSyst, bkgRatesUsual, 
                                    "YieldsAndBackgrounds/subtractBackground.C");

 
  //
  // debug plots
  //
  const int makeDebugPlots=1;
  if (makeDebugPlots) {
    if (0) {
      const int plot_idx=0;
      const int fnc_of_rapidity=0;
      TMatrixD tmpErr=zeeMCShapeReweight;
      tmpErr=0;      
      PlotMatrixMYSlices(plot_idx,fnc_of_rapidity,zeeMCShapeReweight,  "zeeMCShapeReweight");
    }

    if (1) {
      std::vector<int> indices;
      indices.push_back(0);
      std::vector<TMatrixD> matrV;
      std::vector<TMatrixD> matrErrV;
      std::vector<TString> labelV;
      TMatrixD tmpErr=zeeMCShapeReweight;  tmpErr=0;
      TMatrixD mcRew=zeePredictedYield;
      for (int i=0; i<zeePredictedYield.GetNrows(); ++i) {
	for (int j=0; j<zeePredictedYield.GetNcols(); ++j) {
	  mcRew(i,j) = zeePredictedYield(i,j) * zeeMCShapeReweight(i,j);
	}
      }
      matrV.push_back(signalYields); matrV.push_back(zeePredictedYield);
      matrV.push_back(mcRew);
      matrErrV.push_back(tmpErr); matrErrV.push_back(tmpErr);
      matrErrV.push_back(tmpErr);
      labelV.reserve(matrV.size());
      labelV.push_back("signalYields (data)");
      labelV.push_back("zeeYields (MC)");
      labelV.push_back("reweighted MC");
      PlotMatrixMYSlices(indices,0,matrV, matrErrV, labelV, "dataVsMC",
		       "hist", NULL, "dataVsMC");
    }

    if (1) {
      const int iYBin=0;
      const int perMassBinWidth=0;
      const int perRapidityBinWidth=0;
      TCanvas *canvZpeak=MakeCanvas("canvZpeak","canvZpeak",800,900);
      TH1F* hDataNoBkg=extractMassDependence("hDataNoBkg","", 
					     signalYields,signalYieldsError,
					     iYBin,
					     perMassBinWidth,perRapidityBinWidth);
      TH1F *hZee=extractMassDependence("hZee","",
				       zeePredictedYield,zeePredictedYieldErr,
				       iYBin,
				       perMassBinWidth,perRapidityBinWidth);
      if (!hZee) {
	std::cout << "\n\n\tError: failed to locate Zee sample\n\n";
      }
      else {
	ComparisonPlot_t cp(ComparisonPlot_t::_ratioPlain,"compPlot","",
			    "mass [GeV]", "counts", "MC/data");
	std::cout << "hZee normalization factor=" << (hDataNoBkg->Integral()/hZee->Integral()) << "\n";
	hZee->Scale(hDataNoBkg->Integral()/hZee->Integral());
	hZee->SetMarkerStyle(24);
	double dy=0.2;
	cp.SetRatioYRange(1-dy,1+dy);
#ifndef _check_Zpeak
	cp.SetLogx();
#endif
	removeError(hZee);
	canvZpeak->Divide(1,2);
	cp.PreparePads(canvZpeak);
	cp.AddHist1D(hDataNoBkg, "data signal", "LPE", kBlack, 1,0,1);
	cp.AddHist1D(hZee, "MC (normalized)", "LP same", kBlue, 1,0,1);
	cp.Draw(canvZpeak,false,"png");
	SaveCanvas(canvZpeak,canvZpeak->GetName());
      }
    }

  }


  return "Ok";

}

Bool_t checkMatrixSize(TMatrixD m)
{  
  if (m.GetNrows()==DYTools::nMassBins) {
    if ((m.GetNcols()==DYTools::findMaxYBins()) ||
	(m.GetNcols()==DYTools::nYBinsMax) ) {
      return 1;
    }
  }

  std::cout << "m.dims (" << m.GetNrows() << " x " << m.GetNcols() << ") instead of expected (" << DYTools::nMassBins << " x " << DYTools::findMaxYBins() << ") or (" << DYTools::nMassBins << " x " << DYTools::nYBinsMax << ")" << std::endl;
  return 0;
}


void unsquareMatrixElements(TMatrixD &m) {
  for (int i=0; i<m.GetNrows(); ++i) {
    for (int j=0; j<m.GetNcols(); ++j) {
      m(i,j) = sqrt(m(i,j));
    }
  }
}
