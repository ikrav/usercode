#include "TFile.h"
#include "TVectorD.h"
#include "TH1F.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "../Include/CSample.hh"  // helper class for organizing input ntuple files


void subtractBackground(const TString conf){

  std::cout << "conf=" << conf << "\n";

  // Read from configuration file only the location of the root files
  TString inputDir;
  Double_t lumi;
  Bool_t doWeight;
  ifstream ifs;
  ifs.open(conf.Data());
  if (!ifs.is_open()) {
    std::cout << "failed to open file <" << conf << ">\n";
    assert(ifs.is_open());
  }

  string line;
  vector<TString>  snamev;    // sample name (for output file)
  vector<CSample*> samplev;   // data/MC samples
  Int_t state=0;
  while(getline(ifs,line)) {
    if(line[0]=='#') continue;
    if(line[0]=='%') { 
      state++; 
      continue; 
    }
    if(line[0]=='$') {
      samplev.push_back(new CSample());
      stringstream ss(line);
      string chr;
      string sname;
      Int_t color;
      ss >> chr >> sname >> color;
      string label = line.substr(line.find('@')+1);
      snamev.push_back(sname);
      samplev.back()->label = label;
      samplev.back()->color = color;
      continue;
    }
    
    if(state==0) {  // general settings
      stringstream ss1(line); ss1 >> lumi;
      getline(ifs,line);
      stringstream ss2(line); ss2 >> doWeight;
      getline(ifs,line);
      inputDir = TString(line);
      getline(ifs,line);
      // backwards compatibility for the input file
      if (line.size()>3) {  // escale is defined
	TString escaleTag_loc=TString(line);
	getline(ifs,line);
	// check that it was correct to use this work-around
	if (line.find('%')!=std::string::npos) {
	  std::cout << "backwards-compatibility code failure\n";
	  return;
	}
      }
      TString format_loc = TString(line);
      
    } else if(state==1) {  // define data sample
      string fname;
      Double_t xsec;
      string json;
      stringstream ss(line);
      ss >> fname >> xsec >> json;
      samplev.back()->fnamev.push_back(fname);
      samplev.back()->xsecv.push_back(xsec);
      samplev.back()->jsonv.push_back(json);
    
    } else if(state==2) {  // define MC samples
      string fname;
      Double_t xsec;
      stringstream ss(line);
      ss >> fname >> xsec;
      samplev.back()->fnamev.push_back(fname);
      samplev.back()->xsecv.push_back(xsec);
    }
  }

  ifs.close();
  inputDir.ReplaceAll("selected_events","yields");
  std::cout << "inputDir=" << inputDir << "\n";

  // check the requested samples
  //bool zeeMcReq=false; // have it always "true"
  bool zttMcReq=false;
  bool qcdMcReq=false;
  bool ttbarMcReq=false;
  bool wjetsMcReq=false;
  bool wwMcReq=false;
  bool wzMcReq=false;
  bool zzMcReq=false;
  for (unsigned int i=0; i<snamev.size(); ++i) {
    if (snamev[i] == "ttbar") ttbarMcReq=true;
    else if (snamev[i] == "wjets") wjetsMcReq=true;
    else if (snamev[i] == "ww") wwMcReq=true;
    else if (snamev[i] == "wz") wzMcReq=true;
    else if (snamev[i] == "zz") zzMcReq=true;
    else if (snamev[i] == "ztt") zttMcReq=true;
    else if (snamev[i] == "qcd") qcdMcReq=true;
    //else if (snamev[i] == "zee") zeeMcReq=true;
  }

  if (!ttbarMcReq) std::cout << "\n\tWarning: ttbar is not requested\n\n";
  if (!wjetsMcReq) std::cout << "\n\tWarning: wjets is not requested\n\n";
  if (!wwMcReq) std::cout << "\n\tWarning: ww is not requested\n\n";
  if (!wzMcReq) std::cout << "\n\tWarning: wz is not requested\n\n";
  if (!zzMcReq) std::cout << "\n\tWarning: zz is not requested\n\n";
  if (!zttMcReq) std::cout << "\n\tWarning: ztt is not requested\n\n";
  if (!qcdMcReq) std::cout << "\n\tWarning: qcd is not requested\n\n";

  TFile file(inputDir+TString("/massHist.root"));
  if (!file.IsOpen()) {
    std::cout << "failed to open a file <" << inputDir << "/massHist.root>\n";
    throw 2;
  }

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
  if( !data) { std::cout << "data is requested but missing\n"; return; }
  if( !zee   && zeeMc   ) { std::cout << "zeeMc is requested but missing\n"; return; }
  if( zttMcReq   && !ztt   && zttMc   ) { std::cout << "zttMc is requested but missing\n"; return; }
  if( qcdMcReq   && !qcd   && qcdMc   ) { std::cout << "qcd is requested but missing\n"; return; }
  if( ttbarMcReq && !ttbar && ttbarMc ) { std::cout << "ttbar is requested but missing\n"; return; }
  if( wjetsMcReq && !wjets && wjetsMc ) { std::cout << "wjets is requested but missing\n"; return; }
  if( wwMcReq    && !ww    && wwMc    ) { std::cout << "ww is requested but missing\n"; return; }
  if( wzMcReq    && !wz    && wzMc    ) { std::cout << "wz is requested but missing\n"; return; }
  if( zzMcReq    && !zz    && zzMc    ) { std::cout << "zz is requested but missing\n"; return; }

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
  bool useTrue2eBgDataDriven = true;
  TFile fTrueDataDriven(inputDir+TString("/true2eBkgDataPoints.root"));
  TVectorD true2eBackgroundFromData          = *(TVectorD*)fTrueDataDriven.Get("true2eBackgroundFromData");
  TVectorD true2eBackgroundFromDataError     = *(TVectorD*)fTrueDataDriven.Get("true2eBackgroundFromDataError");
  TVectorD true2eBackgroundFromDataErrorSyst = *(TVectorD*)fTrueDataDriven.Get("true2eBackgroundFromDataErrorSyst");
  bool useFakeBgDataDriven = true;
  TFile fFakeDataDriven(inputDir+TString("/fakeBkgDataPoints.root"));
  TVectorD fakeEleBackgroundFromData          = *(TVectorD*)fFakeDataDriven.Get("fakeBackgroundFromData");
  TVectorD fakeEleBackgroundFromDataError     = *(TVectorD*)fFakeDataDriven.Get("fakeBackgroundFromDataError");
  TVectorD fakeEleBackgroundFromDataErrorSyst = *(TVectorD*)fFakeDataDriven.Get("fakeBackgroundFromDataErrorSyst");

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
  if (!fileOut.IsOpen()) {
    std::cout << "failed to create a file <" << inputDir << "/yields_bg-subtracted.root>\n";
    throw 2;
  }
  signalYields         .Write("YieldsSignal");
  signalYieldsError    .Write("YieldsSignalErr");
  signalYieldsErrorSyst.Write("YieldsSignalSystErr");
  binLimits            .Write("BinLimitsForYields");
  fileOut.Close();

}
