#include "TFile.h"
#include "TVectorD.h"
#include "TH1F.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

void subtractBackground(const TString conf){


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

  TFile file(inputDir+TString("/massHist.root"));

  TH1F *data  = (TH1F*) file.Get("data"); 
  TH1F *zee   = (TH1F*) file.Get("zee");   bool zeeMc = true;
  TH1F *ztt   = (TH1F*) file.Get("ztt");   bool zttMc = true;
  TH1F *qcd   = (TH1F*) file.Get("qcd");   bool qcdMc = true;
  TH1F *ttbar = (TH1F*) file.Get("ttbar"); bool ttbarMc = true;
  TH1F *wjets = (TH1F*) file.Get("wjets"); bool wjetsMc = true;
  TH1F *ww   = (TH1F*) file.Get("ww");     bool wwMc = true;
  TH1F *wz   = (TH1F*) file.Get("wz");     bool wzMc = true;
  TH1F *zz   = (TH1F*) file.Get("zz");     bool zzMc = true;

  // Make sure that all MC predictions that are expected according
  // to boolean flags above are present. Data has to be present always.
  if( !data) return;
  if( !zee   && zeeMc   ) return;
  if( !ztt   && zttMc   ) return;
  if( !qcd   && qcdMc   ) return;
  if( !ttbar && ttbarMc ) return;
  if( !wjets && wjetsMc ) return;
  if( !ww    && wwMc    ) return;
  if( !wz    && wzMc    ) return;
  if( !zz    && zzMc    ) return;

  // Close the file in such a way that TTrees do not disappear
  data->SetDirectory(0);
  if( zee  ) zee->SetDirectory(0);
  if( ztt  ) ztt->SetDirectory(0);
  if( qcd  ) qcd->SetDirectory(0);
  if( ttbar) ttbar->SetDirectory(0);
  if( wjets) wjets->SetDirectory(0);
  if( ww   ) ww->SetDirectory(0);
  if( wz   ) wz->SetDirectory(0);
  if( zz   ) zz->SetDirectory(0);
  file.Close();

  // Read data driven background estimates
  bool useTrue2eBgDataDriven = false;
  TFile fTrueDataDriven(inputDir+TString("/DY_bg_true_data-driven_tmp20110719.root"));
  TVectorD true2eBackgroundFromData          = *(TVectorD*)fTrueDataDriven.Get("true2eBackgroundFromData");
  TVectorD true2eBackgroundFromDataError     = *(TVectorD*)fTrueDataDriven.Get("true2eBackgroundFromDataError");
  TVectorD true2eBackgroundFromDataErrorSyst = *(TVectorD*)fTrueDataDriven.Get("true2eBackgroundFromDataErrorSyst");
  bool useFakeBgDataDriven = false;
  TFile fFakeDataDriven(inputDir+TString("/DY_bg_fake_data-driven_tmp20110719.root"));
  TVectorD fakeEleBackgroundFromData          = *(TVectorD*)fFakeDataDriven.Get("fakeEleBackgroundFromData");
  TVectorD fakeEleBackgroundFromDataError     = *(TVectorD*)fFakeDataDriven.Get("fakeEleBackgroundFromDataError");
  TVectorD fakeEleBackgroundFromDataErrorSyst = *(TVectorD*)fFakeDataDriven.Get("fakeEleBackgroundFromDataErrorSyst");

  // A few consistency checks
  bool checkResult = true;
  int nMassBins = data->GetNbinsX();
  if( useTrue2eBgDataDriven && nMassBins != true2eBackgroundFromData.GetNoElements()) checkResult = false;
  if( useFakeBgDataDriven   && nMassBins != fakeEleBackgroundFromData.GetNoElements()) checkResult = false;
  if( !checkResult){
    printf("ERROR: inconsistent binning in the inputs\n");
    assert(0);
  }else
    printf("Binning in the inputs is consistent\n");

  // figure out mass binning
  const TArrayD *binsArray = data->GetXaxis()->GetXbins();
  TVectorD binLimits(nMassBins+1);
  for(int i=0; i<=nMassBins; i++)
    binLimits[i] = binsArray->GetAt(i);

  // Define several arrays we will need
  // The primary output: the arrays that will be written for use
  // in the next stages of the analysis
  TVectorD observedYields(nMassBins);
  TVectorD observedYieldsError(nMassBins);
  TVectorD signalYields(nMassBins);
  TVectorD signalYieldsError(nMassBins);
  TVectorD signalYieldsErrorSyst(nMassBins);
  // Arrays to store backgrounds
  TVectorD true2eBackground(nMassBins);
  TVectorD true2eBackgroundError(nMassBins);
  TVectorD true2eBackgroundErrorSyst(nMassBins);
  TVectorD wzzz(nMassBins);
  TVectorD wzzzError(nMassBins);
  TVectorD wzzzErrorSyst(nMassBins);
  TVectorD fakeEleBackground(nMassBins);
  TVectorD fakeEleBackgroundError(nMassBins);
  TVectorD fakeEleBackgroundErrorSyst(nMassBins);
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
  TVectorD totalBackground(nMassBins);
  TVectorD totalBackgroundError(nMassBins);
  TVectorD totalBackgroundErrorSyst(nMassBins);

  // Loop over mass bins and calculate total background and its error
  printf("Calculate total backgrounds\n");
  if(useTrue2eBgDataDriven)
    printf("    use data driven estimates for true dielectron backgrounds\n");
  else
    printf("    use MC estimates for true dielectron backgrounds\n");
  if(useFakeBgDataDriven)
    printf("    use data driven estimates for fake dielectron backgrounds\n");
  else
    printf("    use MC estimates for fake dielectron backgrounds\n");
  for(int i=0; i<nMassBins; i++){
    // Calculate true dielectron background, which includes
    // WW, ttbar, and DY->tautau. By choice, we do not include WZ and ZZ.
    if(useTrue2eBgDataDriven){
      // Use data-driven background measurement
      true2eBackground[i]          = true2eBackgroundFromData[i];
      true2eBackgroundError[i]     = true2eBackgroundFromDataError[i];
      true2eBackgroundErrorSyst[i] = true2eBackgroundFromDataErrorSyst[i];
    }else{
      // Use Monte Carlo
      double zttYield = 0, zttError = 0;
      double ttbarYield = 0, ttbarError = 0;
      double wwYield = 0, wwError = 0;
      if(ztt){
	zttYield = ztt->GetBinContent(i+1);
	zttError = ztt->GetBinError(i+1);
      }
      if(ttbar){
	ttbarYield = ttbar->GetBinContent(i+1);
	ttbarError = ttbar->GetBinError(i+1);
      }
      if(ww){
	wwYield = ww->GetBinContent(i+1);
	wwError = ww->GetBinError(i+1);
      }
      true2eBackground[i]      = zttYield + ttbarYield + wwYield;
      true2eBackgroundError[i] = sqrt(  zttError * zttError
					+ ttbarError * ttbarError
					+ wwError * wwError);
      // Use ballpark numbers: 0% systematics on DY->tautau, 50% on ttbar, 100% on WW
      true2eBackgroundErrorSyst[i] = sqrt(  0*0*zttYield * zttYield
					    + 0.5*0.5*ttbarYield * ttbarYield
					    + 1.0*1.0*wwYield * wwYield);
    }
    // Calculate WZ and ZZ backgrounds
    double wzYield = 0, wzError = 0;
    double zzYield = 0, zzError = 0;
    if( wz ){
      wzYield = wz->GetBinContent(i+1);
      wzError = wz->GetBinError(i+1);
    }
    if( zz ){
      zzYield = zz->GetBinContent(i+1);
      zzError = zz->GetBinError(i+1);
    }
    wzzz[i] = wzYield + zzYield;
    wzzzError[i] = sqrt( wzError*wzError + zzError*zzError);
    // Use 100% error on MC estimate of WZ and ZZ backgrounds
    wzzzErrorSyst[i] = wzzz[i];

    // Calculate fake dielectrion background, which includes
    // QCD and W+jets
    if(useFakeBgDataDriven){
      // Use data-driven background measurement
      fakeEleBackground[i]          = fakeEleBackgroundFromData[i];
      fakeEleBackgroundError[i]     = fakeEleBackgroundFromDataError[i];
      fakeEleBackgroundErrorSyst[i] = fakeEleBackgroundFromDataErrorSyst[i];
    }else{
      // Use Monte Carlo
      double qcdYield = 0, qcdError = 0;
      double wjetsYield = 0, wjetsError = 0;
      if(qcd){
	qcdYield = qcd->GetBinContent(i+1);
	qcdError = qcd->GetBinError(i+1);
      }
      if(wjets){
	wjetsYield = wjets->GetBinContent(i+1);
	wjetsError = wjets->GetBinError(i+1);
      }
      fakeEleBackground[i]      = qcdYield + wjetsYield;
      fakeEleBackgroundError[i] = sqrt(  qcdError * qcdError
					+ wjetsError * wjetsError);
      // Use ballpark numbers: 50% on QCD, 50% on W+jets
      fakeEleBackgroundErrorSyst[i] = sqrt(  0.5*0.5*qcdYield * qcdYield
					    + 0.5*0.5*wjetsYield * wjetsYield);
    }

    // Calculate the total background
    totalBackground[i] = true2eBackground[i] + wzzz[i] + fakeEleBackground[i];
    totalBackgroundError[i] = sqrt( true2eBackgroundError[i] * true2eBackgroundError[i] +
				    wzzzError[i] * wzzzError[i] +
				    fakeEleBackgroundError[i] * fakeEleBackgroundError[i] );
    totalBackgroundErrorSyst[i] = sqrt( true2eBackgroundErrorSyst[i] * true2eBackgroundErrorSyst[i] +
				    wzzzErrorSyst[i] * wzzzErrorSyst[i] +
				    fakeEleBackgroundErrorSyst[i] * fakeEleBackgroundErrorSyst[i] );
  }

  // Loop over bins and perform background subtraction
  printf("Subtract background from observed yields\n");
  for(int i=0; i<nMassBins; i++){
    observedYields[i] = data->GetBinContent(i+1);
    observedYieldsError[i] = data->GetBinError(i+1);
    signalYields[i] = observedYields[i] - totalBackground[i];
    signalYieldsError[i] = sqrt( observedYieldsError[i] * observedYieldsError[i] + 
				 totalBackgroundError[i] * totalBackgroundError[i] );
    signalYieldsErrorSyst[i] = totalBackgroundErrorSyst[i];
  }

  // Print tables of yields and background

  // Table 1: split background into true, wz/zz, and qcd
  printf(" Note: stat error in signal yield contain stat error on background,\n");
  printf("   and syst error on signal yield contains syst error on background\n");
  printf("mass range            observed       true2e-bg         wz-zz-bg               fake-bg                 total-bg            signal\n");
  for(int i=0; i<nMassBins; i++){
    printf("%5.1f-%5.1f GeV: ", binLimits[i], binLimits[i+1]);
    printf(" %7.0f+-%3.0f ", observedYields[i], observedYieldsError[i]);
    printf(" %5.1f+-%4.1f+-%4.1f ", true2eBackground[i], true2eBackgroundError[i], true2eBackgroundErrorSyst[i]);
    printf(" %6.2f+-%4.2f+-%4.2f ", wzzz[i], wzzzError[i], wzzzErrorSyst[i]);
    printf(" %5.1f+-%5.1f+-%5.1f ", fakeEleBackground[i], fakeEleBackgroundError[i], fakeEleBackgroundErrorSyst[i]);
    printf("    %5.1f+-%4.1f+-%4.1f ", totalBackground[i], totalBackgroundError[i], totalBackgroundErrorSyst[i]);
    printf("    %8.1f+-%5.1f+-%5.1f ", signalYields[i], signalYieldsError[i], signalYieldsErrorSyst[i]);
    printf("\n");
  }

  // Table 2: combined true2e and WZ/ZZ backgrounds only
  printf("\n  only true2e-bg + ww-wz\n");
  printf("mass range      true2e, includingwz/zz\n");
  for(int i=0; i<nMassBins; i++){
    printf("%5.1f-%5.1f GeV: ", binLimits[i], binLimits[i+1]);
    double val = true2eBackground[i] + wzzz[i];
    double err = sqrt(true2eBackgroundError[i]*true2eBackgroundError[i]
		      + wzzzError[i]*wzzzError[i]);
    double sys = sqrt(true2eBackgroundErrorSyst[i]*true2eBackgroundErrorSyst[i]
		      + wzzzErrorSyst[i]*wzzzErrorSyst[i]);
    printf(" %5.1f+-%4.1f+-%4.1f ", val,err, sys);
    printf("\n");
  }

  // Table 3: Systematic error on signal yields assuming that it includes
  // only the syst. error on the background.
  printf("\n  Systematics, %% relative to background subtracted yields\n");
  printf("mass range            subtr-signal    total-bg      syst-from-bg-frac      syst-from-bg-percent\n");
  for(int i=0; i<nMassBins; i++){
    printf("%5.1f-%5.1f GeV: ", binLimits[i], binLimits[i+1]);
    printf("    %8.1f+-%5.1f+-%4.1f ", signalYields[i], signalYieldsError[i],signalYieldsErrorSyst[i]);
    printf("    %5.1f+-%4.1f+-%4.1f ", totalBackground[i], totalBackgroundError[i], totalBackgroundErrorSyst[i]);
    printf("    %6.4f ", totalBackgroundErrorSyst[i]/signalYields[i]);
    printf("    %6.1f ", totalBackgroundErrorSyst[i]*100.0/signalYields[i]);
    printf("\n");
  }

  // Save sideband-subtracted signal yields
  TFile fileOut(inputDir+TString("/yields_bg-subtracted.root"),"recreate");
  signalYields         .Write("YieldsSignal");
  signalYieldsError    .Write("YieldsSignalErr");
  signalYieldsErrorSyst.Write("YieldsSignalSystErr");
  binLimits            .Write("BinLimitsForYields");
  fileOut.Close();

}
