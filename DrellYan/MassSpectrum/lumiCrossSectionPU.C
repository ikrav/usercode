#include <TROOT.h>
#include <TBenchmark.h>             // class to track macro running statistics
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TVectorD.h>
#include <TVector.h>
#include <TGaxis.h>
#include <vector>                   // STL vector class
#include <iostream>                 // dard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings
#include "../Include/CSample.hh"
#include "../Include/ZeeData.hh"
#include "../Include/ElectronEnergyScale.hh"
#include "../Include/DYTools.hh"
#include "../MassSpectrum/storeData.hh"
#include "../Include/MyTools.hh"
#include "../Include/MitStyleRemix.hh"
#include "../Include/CPlot.hh"


// Possible various combinations of the switches
// 
// rebinLuminosity=2 -- use 5 data taking periods (built-in)
// rebinLuminosity=0 -- use luminosity as given in the file (experimental!)
// rebinLuminosity=1 -- load luminosity file (e.g. run_luminosity.root) 
//          and bunch the run-ranges according to what is 
//          luminosityBlockSize value
// rebinLuminosity=3 -- load luminosity file, bin according to the run range
//          the width is controlled by runRangeWidth
//

// run_luminosity*root  path is changed to "./"

//const char *luminosityFileName="run_luminosity_ver1.root"; // incorrect, yet working root file
//const char *luminosityFileName="run_luminosity.root"; // path will be prepended
//const char *luminosityFileName="run_luminosity150.root"; // path will be prepended
const double luminosityBlockSize=150;
const int rebinLuminosity=0;  // set to 1 for run_luminosity.root
const char *luminosityFileName="run_luminosity5_ik.root"; // path will be prepended   
//const double luminosityBlockSize=1500;

const int runRangeWidth=2000;


//  A correction to make largest luminosity version
//  to appear first in the list
const int swapLumiMethods=1; 


// Correct for no Zs after run about 179890
const int bruteForceCorrectionForLumi=1;
const UInt_t bruteForceCutOff=179890;

// specify whether rho_vs_pu.root should be taken from a local directory
const int localRhoVsPUFile=0;

const int xAxis_in_RR=0;

// ---------------------------------------------------------------------

void PrepareHistoStyle(TH1F* h, int color, int markerStyle=20) {
  h->GetXaxis()->SetTitleOffset(1.10);
  h->GetYaxis()->SetTitleOffset(1.40);
  h->SetLineColor(color);
  h->SetMarkerColor(color);
  h->SetMarkerStyle(markerStyle);
}


// ---------------------------------------------------------------------

void PrintHisto(TH1F *histo) {
  std::cout << "histo: " << histo->GetName() << "\n";
  std::cout << " center  value   valueerr\n";
  for(int i=1; i<=histo->GetNbinsX(); i++) {
    printf(" %f    %f    %f\n",histo->GetBinCenter(i),histo->GetBinContent(i),histo->GetBinError(i));
  }
}

// ---------------------------------------------------------------------

void PrintHisto(TH1F *histo1, TH1F *histo2) {
  std::cout << "histo1: " << histo1->GetName() << ", histo2: " << histo2->GetName() << "\n";
  std::cout << " center  value1 value1err   value2 value2err\n";
  for(int i=1; i<=histo1->GetNbinsX(); i++) {
    printf(" %f    %f    %f",histo1->GetBinCenter(i),histo1->GetBinContent(i),histo1->GetBinError(i));
    printf("    ");
    printf("   %f    %f\n",histo2->GetBinContent(i),histo2->GetBinError(i));
  }
}

// ---------------------------------------------------------------------

void PrintVecs(const char *msg, const std::vector<TString> &names, const vector<vector<double>*> &numbers, const vector<TString> &format);
//void PrepareLuminosity(int version, std::vector<LumiInfo_t> &luminosity);
int PrepareLuminosity(int version, std::vector<LumiInfo_t> &luminosity, const TString &fname, double dLumi, int &lumiSectionCount, double **lumiSections );

int LoadEffScaleFactors(const TString &fname, TVectorD &rho_barrel, TVectorD &rho_endcap);

void MakePlots(const char *title, const char *name1, std::vector<TH1F*> &data, const char *name2, std::vector<TH1F*> &data2, const std::vector<int> &colors, int plotTwo, int plot_set);

int WriteToDir(TFile *fout, const TString &dirName, const std::vector<std::vector<TH1F*>*> &hVV);


// ---------------------------------------------------------------------

void lumiCrossSectionPU(TString conf = "../config_files/data.conf") {

  gBenchmark->Start("lumiCrossSectionPU");

  //--------------------------------------------------------------------------------------------------------------
  // Global constants and useful variables
  //==============================================================================================================


  char buf[200];
  const double factor_acceptance = 5598660/double(11428900);
  const double factor_fsr= 11428900/double(11840986);
   
  const int applySmearing=1;

  vector<int> colorV;
  vector<int> marker;

  int excludeJuly2011runs=1;
  LumiInfo_t lumiJuly2011(171050,171578,1.0);

  colorV.reserve(5);
  colorV.push_back(kBlack); colorV.push_back(kBlack); 
  colorV.push_back(46); colorV.push_back(kGreen+2); colorV.push_back(kBlue+1);
  marker.reserve(5);
  marker.push_back(20); 
  marker.push_back(20); marker.push_back(24); 
  marker.push_back(22); marker.push_back(26);
  
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================

  TString  outputDir;  // output directory
  Double_t totalLuminosity;              // luminosity (pb^-1)
  Bool_t   doWeight;          // weight events?
  TString  format;            // plot format

  vector<TString>  snamev;    // sample name (for output file)
  vector<CSample*> samplev;   // data/MC samples
    
  //
  // parse .conf file
  //
  ifstream ifs;
  ifs.open(conf.Data());
  assert(ifs.is_open());
  string line;
  Int_t state=0;
  TString escaleTag = "Date20110901_EPS11_default";
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
      stringstream ss1(line); ss1 >> totalLuminosity;
      getline(ifs,line);
      stringstream ss2(line); ss2 >> doWeight;
      getline(ifs,line);
      outputDir = TString(line);
      getline(ifs,line);
      // backwards compatibility for the input file
      if (line.size()>3) {  // escale is defined
	escaleTag=TString(line);
	getline(ifs,line);
	// check that it was correct to use this work-around
	if (line.find('%')!=std::string::npos) {
	  std::cout << "backwards-compatibility code failure\n";
	  return;
	}
      }
      format = TString(line);
      
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
  if (samplev[0]->fnamev.size()==0) {
    std::cout << "no data files specified in the input\n";
    return;
  }
  
  TString yieldsDir=outputDir;
  yieldsDir.Replace(yieldsDir.Index("selected_events"), sizeof("selected_events")-1, "yields");
  TString constantsDir=outputDir;
  constantsDir.Replace(constantsDir.Index("selected_events"), sizeof("selected_events")-1, "constants");
  

  // 
  // Set up energy scale corrections
  //
  ElectronEnergyScale escale(escaleTag);
  assert(escale.isInitialized());
  escale.print();

  //
  // Access samples and fill histograms
  //  
  ZeeData data;
  TFile *infile=0;
  TTree *eventTree=0; 
  vector<TH1F*> hPVv;
  hPVv.reserve(samplev.size());

  /*
  // Open file with number of PV distributions for pile-up reweighting
  const TString fnamePV = outputDir + TString("/npv.root");
  TFile *pvfile = new TFile(fnamePV);
  assert(pvfile);
  for(UInt_t isam=0; isam<samplev.size(); isam++) {
    TString hName = TString("hNGoodPV_")+snamev[isam];
    TH1F *hPVThis = (TH1F*) pvfile->Get(hName); assert(hPVThis);
    if (!hPVThis) { std::cout << "failed to get PV histo <" << hName << ">\n"; return; }
    hPVThis->SetDirectory(0);
    hPVv.push_back(hPVThis);
  }
  pvfile->Close();
  delete pvfile;
  */


  // load efficiency scale factors as a function of PU
  TVectorD rhoBarrel, rhoEndcap;
  TString rhoVsPU_path= (localRhoVsPUFile) ? TString(".") : constantsDir;
  if (!LoadEffScaleFactors(TString(rhoVsPU_path + TString("/rho_vs_pu.root")),rhoBarrel,rhoEndcap)) return;


  // Cache selected events
  vector<SelectedEventData_t> dataV;
  vector<vector<SelectedEventData_t>*> mcV;
  vector<SelectedEventData_t> *selV=NULL;
  vector<SelectedEventData_t> *mcSignal=NULL;

  mcV.reserve(samplev.size());
  for(UInt_t isam=0; isam<samplev.size(); isam++) {
    if (isam==0) {
      selV=&dataV; 
      mcV.push_back(NULL);
    }
    else {
      selV=new vector<SelectedEventData_t>(); 
      mcV.push_back(selV); 
    }
    if ( snamev[isam] == "zee" ) mcSignal=selV;
    
    const TString fname = outputDir + TString("/ntuples/") + snamev[isam] + TString("_select.root");
    cout << "Processing " << fname << "..." << endl;   
    infile = new TFile(fname);
    assert(infile); 

    // Get the TTree and set branch address
    eventTree = (TTree*)infile->Get("Events"); assert(eventTree); 
    eventTree->SetBranchAddress("Events",&data);

    selV->reserve(eventTree->GetEntries());
    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
      eventTree->GetEntry(ientry);
      selV->push_back(SelectedEventData_t(data));
    }
    infile->Close();
    delete infile;
  }

  
  // Load MC efficiency as a function of PU
  TH1F *hMCSignalEff=new TH1F("hMCSignalEff","hMCSignalEff",DYTools::nPVBinCount,DYTools::nPVLimits);
  hMCSignalEff->Sumw2();
  {
    TString effConstFileName=constantsDir + TString("/event_efficiency_constants.root");
    TFile fa(effConstFileName);
    if (!fa.IsOpen()) {
      std::cout << "failed to open the file <" << effConstFileName << ">\n";
      return;
    }
    TVectorD nEventsZPeakPURaw, nPassZPeakPURaw;
    nEventsZPeakPURaw.Read("nEventsZPeakPURawArray");
    nPassZPeakPURaw.Read("nPassZPeakPURawArray");
    if ((nEventsZPeakPURaw.GetNoElements()==0) ||
	(nPassZPeakPURaw.GetNoElements()==0)) {
      std::cout << "nEventsZPeakPURawArray[" << nEventsZPeakPURaw.GetNoElements() << ", nPassZPeakPURawArray[" << nPassZPeakPURaw.GetNoElements() << "]\n";
      return;
    }

    TH1F *hMCSignalTotal=new TH1F("hMCSignalTotal","",DYTools::nPVBinCount,DYTools::nPVLimits);
    hMCSignalTotal->Sumw2();
    for (int i=0; i<nEventsZPeakPURaw.GetNoElements(); ++i) hMCSignalTotal->Fill(i, nEventsZPeakPURaw[i]);
    for (int i=0; i<nPassZPeakPURaw.GetNoElements(); ++i)   hMCSignalEff->Fill(i, nPassZPeakPURaw[i]);
    //PrintHisto(hMCSignalTotal);  PrintHisto(hMCSignalEff);
    if (0) {
      TCanvas *cx= new TCanvas("ctestCount","ctestCount",600,600);
      hMCSignalTotal->GetXaxis()->SetTitle("number of good vertices");
      hMCSignalTotal->GetYaxis()->SetTitle("count");
      hMCSignalTotal->DrawCopy("PE1");
      int ci=kRed+1; hMCSignalEff->SetMarkerColor(ci); hMCSignalEff->SetLineColor(ci);
      hMCSignalEff->DrawCopy("PE1 same");
      cx->Update();
    }

    hMCSignalEff->Divide(hMCSignalTotal);
    for (int i=1; i<hMCSignalEff->GetNbinsX(); ++i) {
      double e=hMCSignalEff->GetBinContent(i);
      hMCSignalEff->SetBinError(i, sqrt(e*(1-e)/hMCSignalTotal->GetBinContent(i)));
    }
    //PrintHisto(hMCSignalEff);
    delete hMCSignalTotal;
    fa.Close();
  }
  if (0) {
    TCanvas *cx= new TCanvas("ctestEff","ctestEff",600,600);
    hMCSignalEff->GetXaxis()->SetTitle("number of good vertices");
    hMCSignalEff->GetYaxis()->SetTitle("efficiency");
    PrepareHistoStyle(hMCSignalEff,kBlack);
    hMCSignalEff->DrawCopy();
    cx->Update();
    cx->SaveAs("canvas_efficiencyPU.pdf");
  }


  // Create PV distributions for MC
  for (unsigned int isam=0; isam<samplev.size(); ++isam) {
    sprintf(buf,"hPV_%s",snamev[isam].Data());
    TH1F *hPVData= new TH1F(buf,"",DYTools::nPVBinCount,DYTools::nPVLimits);
    hPVv.push_back(hPVData);
    const vector<SelectedEventData_t>* dt= (isam==0) ? &dataV : mcV[isam];
    if (!dt) { std::cout << "dt is null for isam=" << isam << "\n"; return ; }
    for (unsigned int i=0; i<dt->size(); ++i) {
      if ((*dt)[i].massInsideRange(60,120,escale,(applySmearing && (isam>0)))) {
	hPVData->Fill((*dt)[i].nGoodPV, (*dt)[i].weight);
      }
    }
    //PrintHisto(hPVData);
  }


  std::vector<std::vector<double>*> mcZContribV;

  // histograms for each lumiCalc version
  std::vector<TH1F*> hLumiSigmaV;
  std::vector<TH1F*> hMCLumiSigmaV;
  //std::vector<TH1F*> hLumiSigmaPerLumiV;
  //std::vector<TH1F*> hMCLumiSigmaPerLumiV;
  std::vector<TH1F*> hLumiZCountV, hLumiSignalCountV, hLumiNbkgrFracV;
  std::vector<TH1F*> hLumiAvgPUV;
  std::vector<TH1F*> hLumiAvgEffV, hLumiAvgRhoV, hLumiAvgRhoEffV;
  std::vector<TH1F*> hLumiLumiV;

  std::vector<std::vector<TH1F*>*> hLumiZCountPUVV;
  std::vector<std::vector<TH1F*>*> hLumiZMCCountPUVV;
  std::vector<std::vector<TH1F*>*> hLumiNbkgrPUVV;
  std::vector<std::vector<TH1F*>*> hLumiSignalPUVV;
  std::vector<std::vector<TH1F*>*> hLumiNbkgrFracPUVV;
  std::vector<std::vector<TH1F*>*> hLumiAvgEffPUVV;
  std::vector<std::vector<TH1F*>*> hLumiAvgRhoPUVV;
  std::vector<std::vector<TH1F*>*> hLumiAvgRhoEffPUVV;
  std::vector<std::vector<TH1F*>*> hLumiSigmaPUVV;
  std::vector<std::vector<TH1F*>*> hLumiSigmaNormPUVV; // sigma * Nobs/Nobs[PU]

  for (int lumiVersion=1; lumiVersion<=4; ++lumiVersion) {
    //if (lumiVersion>2) break;
    int lumiSectionCount=0;
    double *lumiSections=NULL;
    std::vector<LumiInfo_t> lumiInfoV;
    if (lumiVersion==0) { 
      LumiInfo_t lumiTmp(0,200000, totalLuminosity); lumiInfoV.push_back(lumiTmp);  // catch all
    }
    else {
      //TString lumiBlockInfoFile=constantsDir + TString("/") + TString(luminosityFileName);
      TString lumiBlockInfoFile=TString("./") + TString(luminosityFileName);
      if (!PrepareLuminosity(lumiVersion,lumiInfoV,lumiBlockInfoFile,luminosityBlockSize,lumiSectionCount,&lumiSections)) {
	std::cout << "failed to prepare lumiInfoV for version=" << lumiVersion << "\n";
	return;
      }
    }

    sprintf(buf,"hSigma_%d",lumiVersion);
    TH1F *hSigma=new TH1F(buf,"",lumiSectionCount,lumiSections);
    PrepareHistoStyle(hSigma,colorV[lumiVersion],marker[lumiVersion]);
    hLumiSigmaV.push_back(hSigma);

    //sprintf(buf,"hSigmaPerLumi_%d",lumiVersion);
    //TH1F *hSigmaPerLumi=new TH1F(buf,"",lumiSectionCount,lumiSections);
    //PrepareHistoStyle(hSigmaPerLumi,colorV[lumiVersion],marker[lumiVersion]);
    //hLumiSigmaPerLumiV.push_back(hSigmaPerLumi);

    sprintf(buf,"hMCSigma_%d",lumiVersion);
    TH1F *hMCSigma=new TH1F(buf,"",lumiSectionCount,lumiSections);
    PrepareHistoStyle(hMCSigma,colorV[lumiVersion],marker[lumiVersion]);
    hMCLumiSigmaV.push_back(hMCSigma);

    //sprintf(buf,"hMCSigmaPerLumi_%d",lumiVersion);
    //TH1F *hMCSigmaPerLumi=new TH1F(buf,"",lumiSectionCount,lumiSections);
    //PrepareHistoStyle(hMCSigmaPerLumi,colorV[lumiVersion],marker[lumiVersion]);
    //hMCLumiSigmaPerLumiV.push_back(hMCSigmaPerLumi);

    sprintf(buf,"hZCount_%d",lumiVersion);
    TH1F *hZCount=new TH1F(buf,"",lumiSectionCount,lumiSections);
    PrepareHistoStyle(hZCount,colorV[lumiVersion],marker[lumiVersion]);
    hLumiZCountV.push_back(hZCount);

    //sprintf(buf,"hNbkgr_%d",lumiVersion);
    //TH1F *hNbkgr=new TH1F(buf,"",lumiSectionCount,lumiSections);
    //PrepareHistoStyle(hNbkgr,colorV[lumiVersion],marker[lumiVersion]);
    //hLumiNbkgrV.push_back(hNbkgr);

    sprintf(buf,"hSignalCount_%d",lumiVersion);
    TH1F *hSignalCount=new TH1F(buf,"",lumiSectionCount,lumiSections);
    PrepareHistoStyle(hSignalCount,colorV[lumiVersion],marker[lumiVersion]);
    hLumiSignalCountV.push_back(hSignalCount);

    sprintf(buf,"hNbkgrFrac_%d",lumiVersion);
    TH1F *hNbkgrFrac=new TH1F(buf,"",lumiSectionCount,lumiSections);
    PrepareHistoStyle(hNbkgrFrac,colorV[lumiVersion],marker[lumiVersion]);
    hLumiNbkgrFracV.push_back(hNbkgrFrac);

    sprintf(buf,"hAvgPU_%d",lumiVersion);
    TH1F *hAvgPU=new TH1F(buf,"",lumiSectionCount,lumiSections);
    PrepareHistoStyle(hAvgPU,colorV[lumiVersion],marker[lumiVersion]);
    hLumiAvgPUV.push_back(hAvgPU);

    sprintf(buf,"hAvgEff_%d",lumiVersion);
    TH1F *hAvgEff=new TH1F(buf,"",lumiSectionCount,lumiSections);
    PrepareHistoStyle(hAvgEff,colorV[lumiVersion],marker[lumiVersion]);
    hLumiAvgEffV.push_back(hAvgEff);

    sprintf(buf,"hAvgRho_%d",lumiVersion);
    TH1F *hAvgRho=new TH1F(buf,"",lumiSectionCount,lumiSections);
    PrepareHistoStyle(hAvgRho,colorV[lumiVersion],marker[lumiVersion]);
    hLumiAvgRhoV.push_back(hAvgRho);

    sprintf(buf,"hAvgRhoEff_%d",lumiVersion);
    TH1F *hAvgRhoEff=new TH1F(buf,"",lumiSectionCount,lumiSections);
    PrepareHistoStyle(hAvgRhoEff,colorV[lumiVersion],marker[lumiVersion]);
    hLumiAvgRhoEffV.push_back(hAvgRhoEff);

    sprintf(buf,"hLumi_%d",lumiVersion);
    TH1F *hLumi=new TH1F(buf,"",lumiSectionCount,lumiSections);
    PrepareHistoStyle(hLumi,colorV[lumiVersion],marker[lumiVersion]);
    hLumiLumiV.push_back(hLumi);

    std::vector<TH1F*>* hZCountPUV=new std::vector<TH1F*>();
    hLumiZCountPUVV.push_back(hZCountPUV);
    std::vector<TH1F*>* hZMCCountPUV=new std::vector<TH1F*>();
    hLumiZMCCountPUVV.push_back(hZMCCountPUV);
    std::vector<TH1F*>* hNbkgrPUV=new std::vector<TH1F*>();
    hLumiNbkgrPUVV.push_back(hNbkgrPUV);
    std::vector<TH1F*>* hSignalPUV=new std::vector<TH1F*>();
    hLumiSignalPUVV.push_back(hSignalPUV);
    std::vector<TH1F*>* hNbkgrFracPUV=new std::vector<TH1F*>();
    hLumiNbkgrFracPUVV.push_back(hNbkgrFracPUV);
    std::vector<TH1F*>* hAvgEffPUV=new std::vector<TH1F*>();
    hLumiAvgEffPUVV.push_back(hAvgEffPUV);
    std::vector<TH1F*>* hAvgRhoPUV=new std::vector<TH1F*>();
    hLumiAvgRhoPUVV.push_back(hAvgRhoPUV);
    std::vector<TH1F*>* hAvgRhoEffPUV=new std::vector<TH1F*>();
    hLumiAvgRhoEffPUVV.push_back(hAvgRhoEffPUV);
    std::vector<TH1F*>* hSigmaPUV=new std::vector<TH1F*>();
    hLumiSigmaPUVV.push_back(hSigmaPUV);
    std::vector<TH1F*>* hSigmaNormPUV=new std::vector<TH1F*>();
    hLumiSigmaNormPUVV.push_back(hSigmaNormPUV);
 

    std::vector<double> dataZv, mcSignalZv, mcBkgrZv; 
    std::vector<double> sigmaV;
    
    // count per lumi block
    dataZv.reserve(lumiInfoV.size());
    mcSignalZv.reserve(lumiInfoV.size());
    mcBkgrZv.reserve(lumiInfoV.size());
    mcZContribV.reserve(lumiInfoV.size());
    double sumLumi=0; // for checking
    for (unsigned int lumi_i=0; lumi_i<lumiInfoV.size(); ++lumi_i) {
      LumiInfo_t lumi(lumiInfoV[lumi_i]);
      
      // prepare container for sigma[PU] calculations
      sprintf(buf,"hZCountPU_lv%d_%u_%u",lumiVersion,lumi.runNumMin,lumi.runNumMax);
      TH1F* hZCountPU=new TH1F(buf,buf,DYTools::nPVBinCount,DYTools::nPVLimits);
      sprintf(buf,"hZMCCountPU_lv%d_%u_%u",lumiVersion,lumi.runNumMin,lumi.runNumMax);
      TH1F* hZMCCountPU=new TH1F(buf,buf,DYTools::nPVBinCount,DYTools::nPVLimits);
      sprintf(buf,"hNbkgrPU_lv%d_%u_%u",lumiVersion,lumi.runNumMin,lumi.runNumMax);
      TH1F* hNbkgrPU=new TH1F(buf,buf,DYTools::nPVBinCount,DYTools::nPVLimits);
      sprintf(buf,"hSumEffPU_lv%d_%u_%u",lumiVersion,lumi.runNumMin,lumi.runNumMax);
      TH1F* hSumEffPU=new TH1F(buf,buf,DYTools::nPVBinCount,DYTools::nPVLimits);
      sprintf(buf,"hSumRhoPU_lv%d_%u_%u",lumiVersion,lumi.runNumMin,lumi.runNumMax);
      TH1F* hSumRhoPU=new TH1F(buf,buf,DYTools::nPVBinCount,DYTools::nPVLimits);
      sprintf(buf,"hSumRhoEffPU_lv%d_%u_%u",lumiVersion,lumi.runNumMin,lumi.runNumMax);
      TH1F* hSumRhoEffPU=new TH1F(buf,buf,DYTools::nPVBinCount,DYTools::nPVLimits);

      // obtain PV distribution in data
      sprintf(buf,"hPVData_lv%d_%u_%u",lumiVersion,lumi.runNumMin,lumi.runNumMax);
      TH1F *hPVData= new TH1F(buf,buf,DYTools::nPVBinCount,DYTools::nPVLimits);
      double dataZCount=0, bkgrZCount=0;
      double sumEff=0;  // sum of efficiencies for the average
      double sumRho=0; // sum of eff. scale factors for the average
      double sumRhoEff=0; // sum of (eff * rho) for the average
      double sumPU=0; // sum of the number of vertices
      for (unsigned int i=0; i<dataV.size(); ++i) {
	if (excludeJuly2011runs && lumiJuly2011.insideRange(dataV[i].runNum)) continue;
	if (lumi.insideRange(dataV[i].runNum)) {
	  if (dataV[i].nGoodPV>40) std::cout << "\n\n\t\tdata event has " << dataV[i].nGoodPV << "a good primary vertices\n\n\n";
	  if (dataV[i].massInsideRange(60,120) && 
	      dataV[i].nGoodPV)  // at least 1 good vertex
	  {
	    hPVData->Fill(dataV[i].nGoodPV, dataV[i].weight);
	    dataZCount++;                          // count Z candidates in data
	    hZCountPU->Fill(dataV[i].nGoodPV, dataV[i].weight);
	    sumPU += dataV[i].nGoodPV;
	    if (dataV[i].nGoodPV>40) std::cout << "data event has " << dataV[i].nGoodPV << "a good primary vertices\n";
	    int idx=hMCSignalEff->FindBin(dataV[i].nGoodPV);
	    if (idx==0) std::cout << "idx=0 for goodPV=" << dataV[i].nGoodPV << "\n";
	    if ((idx==-1) || (idx>hMCSignalEff->GetNbinsX())) {
	      std::cout << "got bin idx=" << idx << " from hMCSignalEff. dataV[i].nGoodPV=" << dataV[i].nGoodPV << "\n";
	      return;
	    }
	    const double currEff=hMCSignalEff->GetBinContent(idx);
	    sumEff += currEff; // add-up the efficiency for averaging
	    int arr_idx=idx-1; // correction for TVector indexing
	    // first has to be the leading electron!
	    double rho_factor= (dataV[i].firstIsInBarrel()) ? 
	      rhoBarrel[arr_idx] : rhoEndcap[arr_idx];
	    rho_factor *= (dataV[i].secondIsInBarrel()) ? 
	      rhoBarrel[arr_idx] : rhoEndcap[arr_idx];
	    sumRho += rho_factor;
	    sumRhoEff += (rho_factor * currEff);
	    hSumEffPU->Fill(dataV[i].nGoodPV, currEff);
	    hSumRhoPU->Fill(dataV[i].nGoodPV, rho_factor);
	    hSumRhoEffPU->Fill(dataV[i].nGoodPV, currEff * rho_factor);
	  }
	}
      }


      if (dataZCount==0) {
	std::cout << "No valid Z candidates in lumi block " << lumi << "\n";
	continue;
      }

      dataZv.push_back(dataZCount);
      vector<double> *zContrib=new vector<double>();
      zContrib->reserve(samplev.size());
      mcZContribV.push_back(zContrib);
      zContrib->push_back(0);

      // work with backgrounds
      for (UInt_t isam=1; isam<samplev.size(); isam++) {
	if (isam==0) continue;
	
	// prepare PV weights. It's data PV distribution divided by MC PV distribution
	sprintf(buf,"hPVWeights_%s_%u_%u",snamev[isam].Data(),lumi.runNumMin,lumi.runNumMax);
	TH1F* puWeights=(TH1F*)hPVData->Clone(buf);
	puWeights->Scale( hPVv[isam]->GetSumOfWeights() / hPVData->GetSumOfWeights());
	if (0) PrintHisto(puWeights,hPVData);
	puWeights->Divide(hPVv[isam]);
	if (0) PrintHisto(puWeights);
	//
	// add backgrounds
	//
	Double_t Nbkgr=0.;
	const vector<SelectedEventData_t> *bkgr = mcV[isam]; 
	const double lumiReweight=lumi.lumiWeight/totalLuminosity;
	for (vector<SelectedEventData_t>::const_iterator it=bkgr->begin(); it!=bkgr->end(); it++) {
	  if (it->massInsideRange(60,120,escale,applySmearing)) {
	    int pvIdx=puWeights->FindBin( it->nGoodPV );
	    if (pvIdx>0) {
	      const double eventWeight= it->weight;
	      const double puWeight=puWeights->GetBinContent(pvIdx);	    
	      const double full_weight= eventWeight * puWeight * lumiReweight;
	      Nbkgr+=full_weight;
	      if ( snamev[isam] != "zee" ) hNbkgrPU->Fill( it->nGoodPV, full_weight );
	      else hZMCCountPU->Fill( it->nGoodPV, full_weight );
	    }
	  }
	}
	if (1) { std::cout << "Nbkgr=" << Nbkgr << ", hNbkgrPU->Sum=" << hNbkgrPU->Integral() << "\n"; }
	zContrib->push_back(Nbkgr);
	if ( snamev[isam] != "zee" ) bkgrZCount+=Nbkgr;
	else mcSignalZv.push_back(Nbkgr);
      }
      mcBkgrZv.push_back(bkgrZCount);

      double Nsignal=dataZCount-bkgrZCount;
      double avgEff=sumEff/dataZCount;
      double avgRho=sumRho/dataZCount;
      double avgRhoEff=sumRhoEff/dataZCount;
      double sigma=Nsignal/( lumi.lumiWeight * avgRhoEff * factor_acceptance * factor_fsr );
      double sigmaMC= mcSignalZv.back()/ ( lumi.lumiWeight * avgRhoEff * factor_acceptance * factor_fsr );
      double xLumi=sumLumi + 0.5* lumi.lumiWeight;
      if (xAxis_in_RR) xLumi=0.5*(lumi.runNumMin+lumi.runNumMax);
      double avgPU=sumPU/dataZCount;
      hSigma->Fill( xLumi, sigma );
      hMCSigma->Fill( xLumi, sigmaMC );
      //hSigmaPerLumi->Fill( sumLumi + 0.5* lumi.lumiWeight, sigma/lumi.lumiWeight );
      //hMCSigmaPerLumi->Fill( sumLumi + 0.5 * lumi.lumiWeight, sigmaMC/lumi.lumiWeight );
      hZCount->Fill( xLumi, dataZCount );
      hSignalCount->Fill( xLumi, Nsignal );
      hNbkgrFrac->Fill( xLumi, 100.*bkgrZCount/double(dataZCount) );
      hAvgPU->Fill( xLumi, avgPU );
      hAvgEff->Fill( xLumi, avgEff );
      hAvgRho->Fill( xLumi, avgRho );
      hAvgRhoEff->Fill( xLumi, avgRhoEff );
      hLumi->Fill( xLumi, lumi.lumiWeight );
      sumLumi += lumi.lumiWeight;
      std::cout << "sumLumi=" << sumLumi << ", lumiWeight=" << lumi.lumiWeight << ", dataZCount=" << dataZCount << ", bkgrZCount=" << bkgrZCount << ", avgRhoEff=" << avgRhoEff << ", sigma=" << sigma << ", sigmaMC=" << sigmaMC << "\n";


      sprintf(buf,"hNsignalPU_lv%d_%u_%u",lumiVersion,lumi.runNumMin,lumi.runNumMax);
      TH1F* hNSignalPU=new TH1F(buf,buf,DYTools::nPVBinCount,DYTools::nPVLimits);
      for (int i=1; i<hZCountPU->GetNbinsX(); ++i) {
	double nsignal= hZCountPU->GetBinContent(i) - hNbkgrPU->GetBinContent(i);
	hNSignalPU->SetBinContent( i, nsignal );
	double zcountErr= hZCountPU->GetBinError(i);
	double nbkgrErr = hNbkgrPU->GetBinError(i);
	double nsignalErr=sqrt( zcountErr*zcountErr + nbkgrErr*nbkgrErr );
	hNSignalPU->SetBinError( i, nsignalErr );
      }
      PrintHisto(hZCountPU,hNbkgrPU);
      sprintf(buf,"hNbkgrFracPU_lv%d_%u_%u",lumiVersion,lumi.runNumMin,lumi.runNumMax);
      TH1F* hNbkgrFracPU=(TH1F*)hNbkgrPU->Clone(buf);
      hNbkgrFracPU->SetTitle(buf);
      hNbkgrFracPU->Divide(hZCountPU);
      sprintf(buf,"hAvgEffPU_lv%d_%u_%u",lumiVersion,lumi.runNumMin,lumi.runNumMax);
      TH1F *hAvgEffPU= (TH1F*)hSumEffPU->Clone(buf);
      hAvgEffPU->SetTitle(buf);
      hAvgEffPU->Divide(hZCountPU);
      sprintf(buf,"hAvgRhoPU_lv%d_%u_%u",lumiVersion,lumi.runNumMin,lumi.runNumMax);
      TH1F *hAvgRhoPU= (TH1F*)hSumRhoPU->Clone(buf);
      hAvgRhoPU->Divide(hZCountPU);
      sprintf(buf,"hAvgRhoEffPU_lv%d_%u_%u",lumiVersion,lumi.runNumMin,lumi.runNumMax);
      TH1F *hAvgRhoEffPU= (TH1F*)hSumRhoEffPU->Clone(buf);
      hAvgRhoEffPU->SetTitle(buf);
      hAvgRhoEffPU->Divide(hNSignalPU);

      sprintf(buf,"hSigmaPU_lv%d_%u_%u",lumiVersion,lumi.runNumMin,lumi.runNumMax);
      TH1F *hSigmaPU= (TH1F*)hNSignalPU->Clone(buf);
      hSigmaPU->SetTitle(buf);
      hSigmaPU->Scale(1/ ( lumi.lumiWeight * factor_acceptance * factor_fsr ) );
      hSigmaPU->Divide(hAvgRhoEffPU);

      sprintf(buf,"hSigmaNormPU_lv%d_%u_%u",lumiVersion,lumi.runNumMin,lumi.runNumMax);
      TH1F *hSigmaNormPU= (TH1F*)hSigmaPU->Clone(buf);
      hSigmaNormPU->SetTitle(buf);
      hSigmaNormPU->Scale( dataZCount );
      hSigmaNormPU->Divide(hZCountPU);

      hSigmaPUV->push_back(hSigmaPU);
      hSigmaNormPUV->push_back(hSigmaNormPU);
      hZCountPUV->push_back(hZCountPU);
      hZMCCountPUV->push_back(hZMCCountPU);
      hNbkgrPUV->push_back(hNbkgrPU);
      hSignalPUV->push_back(hNSignalPU);
      hNbkgrFracPUV->push_back(hNbkgrFracPU);
      hAvgEffPUV->push_back(hAvgEffPU);
      hAvgRhoPUV->push_back(hAvgRhoPU);
      hAvgRhoEffPUV->push_back(hAvgRhoEffPU);

      delete hSumEffPU;
      delete hSumRhoPU;
      delete hSumRhoEffPU;
    }
  
    vector<vector<double>*> numbersV;
    vector<TString> namesV, formatV;
    namesV.push_back("dataZ"); numbersV.push_back(&dataZv);  formatV.push_back("%9.1lf");
    namesV.push_back("mcSignalZ"); numbersV.push_back(&mcSignalZv); formatV.push_back("%9.1lf");
    namesV.push_back("mcBkgrZ"); numbersV.push_back(&mcBkgrZv); formatV.push_back("%9.1lf");
    PrintVecs("Number of events in the 60-120 mass range",namesV,numbersV,formatV);

    //PrintHisto(hSigma);
    //PrintHisto(hMCSigma);
    //PrintHisto(hSigmaPerLumi);
    //PrintHisto(hMCSigmaPerLumi);
    //std::cout << "Z counts in the luminosity block\n";
    //PrintHisto(hSigma,hMCSigma);
    //PrintHisto(hSigmaPerLumi,hMCSigmaPerLumi);
    //PrintHisto(hZCount,hSignalCount);
    //PrintHisto(hAvgPU);
  }

   int the_set=1;
  MakePlots("Sigma","Sigma", hLumiSigmaV, "avgPU",hLumiLumiV, colorV, 1, the_set);
  //MakePlots("Sigma","Sigma", hLumiSigmaV, "SigmaPerLumi",hLumiSigmaPerLumiV, colorV, 1);
  the_set=2;
  MakePlots("ZCount","ZCount", hLumiZCountV, "NbkgrFrac",hLumiNbkgrFracV, colorV, 1, the_set);
  the_set=3;
  MakePlots("AvgEff","AvgEff", hLumiAvgEffV, "AvgRho",hLumiAvgRhoV, colorV, 1, the_set);
  the_set=4;
  MakePlots("avgPU","avgPU", hLumiAvgPUV,  "dummy", hLumiAvgPUV, colorV, 0, the_set);
  the_set=5;
  MakePlots("avgRhoEff","avgRhoEff", hLumiAvgRhoEffV,  "dummy", hLumiAvgPUV, colorV, 0, the_set);

  TFile *f=new TFile("dataLumiCrossSectionPU.root","recreate");
  WriteToDir(f,"Nobs" ,hLumiZCountPUVV);
  WriteToDir(f,"Nbkgr",hLumiNbkgrPUVV);
  WriteToDir(f,"Nsignal",hLumiSignalPUVV);
  WriteToDir(f,"NMCsignal",hLumiZMCCountPUVV);
  WriteToDir(f,"NbkgrFrac",hLumiNbkgrFracPUVV);
  WriteToDir(f,"AvgEff",hLumiAvgEffPUVV);
  WriteToDir(f,"AvgRho",hLumiAvgRhoPUVV);
  WriteToDir(f,"AvgRhoEff",hLumiAvgRhoEffPUVV);
  WriteToDir(f,"CrossSection",hLumiSigmaPUVV);
  WriteToDir(f,"NormCrossSection",hLumiSigmaNormPUVV);
  f->Close();
  delete f;

  gBenchmark->Show("lumiCrossSectionPU");
  return;
}

// ---------------------------------------------------------------------

void PrintVecs(const char *msg, const std::vector<TString> &names, const vector<vector<double>*> &numbers, const std::vector<TString> &format) {
  std::cout << "\n\n";
  if (msg) std::cout << msg << "\n";
  for (unsigned int i=0; i<names.size(); ++i) {
    printf(" %10s",names[i].Data());
  }
  std::cout << std::endl;
  for (unsigned int j=0; j<numbers[0]->size(); ++j) {
    for (unsigned int i=0; i<numbers.size(); ++i) {
      printf(format[i].Data(),(*numbers[i])[j]);
      printf(" ");
    }
    printf("\n");
  }
  printf("\n");
  return;
}

// ---------------------------------------------------------------------

int PrepareLuminosity(int version, std::vector<LumiInfo_t> &luminosity, const TString &fname, double dLumi, int &lumiSectionCount, double **lumiSections ) {
  if (swapLumiMethods) { 
    if (version==1) version=2;
    else if (version==2) version=1;
  }
  luminosity.clear();
  lumiSectionCount=0;
  if (*lumiSections) delete *lumiSections;

  TVector runNumMin,runNumMax;
  TVectorD lumiWeights;
  if (rebinLuminosity!=2) {
    TFile f(fname);
    if (!f.IsOpen()) {
      std::cout << "failed to open a file <" << fname << ">\n";
      return 0;
  }
    runNumMin.Read("runMin");
    runNumMax.Read("runMax");
    switch(version) {
    case 1: lumiWeights.Read("lumiNoCorr"); break;
    case 2: lumiWeights.Read("lumiCorrV2"); break;
    case 3: lumiWeights.Read("lumiCorrV3"); break;
    case 4: lumiWeights.Read("lumiCorrPix");break;
    default:
      std::cout << "PrepareLuminosity does not know version=" << version << "\n";
      return 0;
    }
    f.Close();
    const int n=runNumMin.GetNoElements();
    if ((n==0) || (n!=runNumMax.GetNoElements()) || (n!=lumiWeights.GetNoElements())) {
      std::cout << "count mismatch : runNumMin[" << runNumMin.GetNoElements() << "], runNumMax[" << runNumMax.GetNoElements() << "], lumiWeight[" << lumiWeights.GetNoElements() << "]\n";
      return 0;
    }

    if (1) {
      std::cout << "\nRaw luminosities :\n";
      for (int i=0; i<runNumMin.GetNoElements(); ++i) {
	std::cout << " run " << runNumMin[i] << " -- " << runNumMax[i] << "    "  << lumiWeights[i] << "\n";
      }
      std::cout << "\n";
    }
  }
  
  const double tolerance=0.1*dLumi;
  double lumiChunk=0;
  double lumi=0;
  LumiInfo_t info;
  int start=1;
  std::vector<double> lumiBins;
  lumiBins.push_back(0.);
  lumiChunk=0;
  if (rebinLuminosity==0) {
    lumi=0.;
    if (xAxis_in_RR) lumiBins[0]=runNumMin[0];
    for (int i=0; i<runNumMin.GetNoElements(); ++i) {
      if (bruteForceCorrectionForLumi && (UInt_t(runNumMin[i])>=bruteForceCutOff)) continue;
      info.assign(UInt_t(runNumMin[i]), UInt_t(runNumMax[i]), lumiWeights[i]);
      luminosity.push_back(info);
      lumi += lumiWeights[i];
      if (xAxis_in_RR) lumiBins.push_back(runNumMax[i]);
      else lumiBins.push_back(lumi);
    }
  }
  else if (rebinLuminosity==2) {
    const char *line1="0 160405 163869       204.124       219.13       217.196       214.995";
    const char *line2="1 165088 167913       884.97       949.633       936.575       953.26";
    const char *line3="2 170249 172619       191.145       203.078       199.382       203.619";
    const char *line4="3 172620 173692       640.639       688.585       674.813       698.526";
    const char *line5="4 175832 180252       2573.89       2730.69       2636.96       2712.5";
    std::vector<std::string> lines;
    lines.push_back(line1); lines.push_back(line2); 
    lines.push_back(line3); lines.push_back(line4);
    lines.push_back(line5);
    int idx,runMin,runMax;
    double lumiW[4];
    lumi=0.;
    for (unsigned int i=0; i<lines.size(); ++i) {
      stringstream ss(lines[i]);
      ss >> idx >> runMin >> runMax >> lumiW[0] >> lumiW[1] >> lumiW[2] >> lumiW[3];
      info.assign(runMin,runMax,lumiW[version-1]);
      luminosity.push_back(info);
      lumi += lumiW[version-1];
      lumiBins.push_back(lumi);
    }
  }
  else if (rebinLuminosity==1) {
    const int n=runNumMin.GetNoElements();
    if (xAxis_in_RR) lumiBins[0]=runNumMin[0];
    for (int i=0; i<runNumMin.GetNoElements(); ++i) {
      if (start) { info.runNumMin=UInt_t(runNumMin[i]); start=0; }
      double lumiChunk1 = lumiChunk + lumiWeights[i];
      int nextChunkIsLarge=((i<n-1) && (lumiWeights[i+1]>dLumi)) ? 1:0;
      int nextChunkIsGood=((i<n-1) && (lumiChunk1 + lumiWeights[i+1] < dLumi)) ? 1:0;
      int done=0;
      if (bruteForceCorrectionForLumi && (i<n-1) && (UInt_t(runNumMin[i+1])>bruteForceCutOff)) { 
	done=1; 
	nextChunkIsLarge=1; 
      }
      if (!rebinLuminosity || ( lumiChunk1 > dLumi) || (( lumiChunk1 > dLumi-tolerance) && !nextChunkIsGood) || nextChunkIsLarge ) {
	info.runNumMax=UInt_t(runNumMax[i]); start=1; lumiChunk=0.;
	info.lumiWeight=lumiChunk1;
	luminosity.push_back(info);
	lumi += lumiChunk1;
	if (xAxis_in_RR) lumiBins.push_back(runNumMax[i]);
	else lumiBins.push_back(lumi);
      }
      else lumiChunk=lumiChunk1;
      if (done) break;
    }
    if (start!=1) {
      info.runNumMax=UInt_t(runNumMax[n-1]);
      info.lumiWeight=lumiChunk;
      luminosity.push_back(info);
      lumi += lumiChunk;
      if (xAxis_in_RR) lumiBins.push_back(runNumMax[n-1]);
      else lumiBins.push_back(lumi);
    }
  }
  else if (rebinLuminosity==3) {
    const int n=runNumMin.GetNoElements();
    UInt_t runMin=UInt_t(runNumMin[0]);
    dLumi=0; lumi=0;
    if (xAxis_in_RR) lumiBins[0]=runNumMin[0];
    for (int i=0; i<n; ++i) {
      UInt_t runMax=UInt_t(runNumMax[i]);
      int done= (i==n-1) ? 1:0;
      if (bruteForceCorrectionForLumi && (UInt_t(runNumMin[i])>bruteForceCutOff)) done=1;
      if (( runMax>= (runMin/runRangeWidth)*runRangeWidth+runRangeWidth ) ||
	  done ) {
	info.assign(runMin, runMax-1, dLumi);
	luminosity.push_back(info);
	lumi+=dLumi;
	if (xAxis_in_RR) lumiBins.push_back(runNumMax[i]);
	else lumiBins.push_back(lumi);
	dLumi=lumiWeights[i];
	runMin=runMax;
      }
      else {
	dLumi+=lumiWeights[i];
      }
      if (done) break;
    }
  }
  else {
    std::cout << "unclear value for rebinLuminosity=" << rebinLuminosity << "\n";
    return 0;
  }

  // create the array
  lumiSectionCount=lumiBins.size()-1;
  *lumiSections=new double[lumiBins.size()];
  for (unsigned int i=0; i<lumiBins.size(); ++i) (*lumiSections)[i]=lumiBins[i];

  // print luminosity vector
  const int print_eol=1;
  PrintVec("luminosityInfo ",luminosity, print_eol);

  return 1;
}

// ---------------------------------------------------------------------

int ConvertVecsToHistos(const std::vector<TVectorD*> &vecs, int binCount, const double *binLimits, std::vector<TH1F*> &histos, const std::vector<TString> &names) {
  histos.reserve(vecs.size());
  for (unsigned int veci=0; veci<vecs.size(); ++veci) {
    const TVectorD *vec=vecs[veci];
    TH1F* h=new TH1F(names[veci],"",binCount,binLimits);
    histos.push_back(h);
    for (int i=0; i<vec->GetNoElements(); ++i) {
      h->SetBinContent(i+1,(*vec)[i]);
      h->SetBinError(i+1, 0.);
    }
  }
  return 1;
}

// ---------------------------------------------------------------------

int ConvertVecsToHistos(const std::vector<TVectorD*> &vecs, const std::vector<TVectorD*> &vecErrs, int binCount, const double *binLimits, std::vector<TH1F*> &histos, const std::vector<TString> &names) {
  histos.reserve(vecs.size());
  for (unsigned int veci=0; veci<vecs.size(); ++veci) {
    const TVectorD *vec=vecs[veci];
    const TVectorD *vecErr=vecErrs[veci];
    TH1F* h=new TH1F(names[veci],"",binCount,binLimits);
    histos.push_back(h);
    for (int i=0; i<vec->GetNoElements(); ++i) {
      h->SetBinContent(i+1,(*vec)[i]);
      h->SetBinError(i+1, (*vecErr)[i]);
    }
  }
  return 1;
}

// ---------------------------------------------------------------------

int LoadEffScaleFactors(const TString &fname, TVectorD &rho_barrel, TVectorD &rho_endcap) {
  //int res=1;
  //int res=1;

  TFile fEffScaleFactors(fname);
  TVector pu_bin_min = *(TVector*)fEffScaleFactors.Get("pu_bin_min");
  TVector pu_bin_max = *(TVector*)fEffScaleFactors.Get("pu_bin_max");
  TVectorD rho_reco_barrel = *(TVectorD*)fEffScaleFactors.Get("rho_reco_barrel");
  TVectorD rho_reco_endcap = *(TVectorD*)fEffScaleFactors.Get("rho_reco_endcap");
  TVectorD rho_id_barrel = *(TVectorD*)fEffScaleFactors.Get("rho_id_barrel");
  TVectorD rho_id_endcap = *(TVectorD*)fEffScaleFactors.Get("rho_id_endcap");
  TVectorD rho_hlt_barrel = *(TVectorD*)fEffScaleFactors.Get("rho_hlt_barrel");
  TVectorD rho_hlt_endcap = *(TVectorD*)fEffScaleFactors.Get("rho_hlt_endcap");

  TVectorD rho_reco_barrel_err = *(TVectorD*)fEffScaleFactors.Get("rho_reco_barrel_Err");
  TVectorD rho_reco_endcap_err = *(TVectorD*)fEffScaleFactors.Get("rho_reco_endcap_Err");
  TVectorD rho_id_barrel_err = *(TVectorD*)fEffScaleFactors.Get("rho_id_barrel_Err");
  TVectorD rho_id_endcap_err = *(TVectorD*)fEffScaleFactors.Get("rho_id_endcap_Err");
  TVectorD rho_hlt_barrel_err = *(TVectorD*)fEffScaleFactors.Get("rho_hlt_barrel_Err");
  TVectorD rho_hlt_endcap_err = *(TVectorD*)fEffScaleFactors.Get("rho_hlt_endcap_Err");

  if (0) {
    std::cout << "rho factors\n";
    std::cout << " pu_bin_min  pu_bin_max  rho_reco_barrel rho_reco_endcap rho_id_barrel rho_id_endcap rho_hlt_barrel rho_hlt_endcap\n";
    for (int i=0; i<rho_reco_barrel.GetNoElements(); ++i) {
      printf("  %2.0lf -- %2.0lf  %6.4lf %6.4lf  %6.4lf %6.4lf  %6.4lf %6.4lf\n",pu_bin_min[i],pu_bin_max[i],rho_reco_barrel[i],rho_reco_endcap[i],rho_id_barrel[i],rho_id_endcap[i],rho_hlt_barrel[i],rho_hlt_endcap[i]);
    }
  }

  if ((rho_reco_barrel.GetNoElements() != DYTools::nPVBinCount) ||
      (rho_reco_endcap.GetNoElements() != DYTools::nPVBinCount) ||
      (rho_id_barrel.GetNoElements() != DYTools::nPVBinCount) ||
      (rho_id_endcap.GetNoElements() != DYTools::nPVBinCount) ||
      (rho_hlt_barrel.GetNoElements() != DYTools::nPVBinCount) ||
      (rho_hlt_endcap.GetNoElements() != DYTools::nPVBinCount)) {
    std::cout << "Number of pv Bin count is different. DYTools::nPVBinCount=" << DYTools::nPVBinCount << ", rho_vs_pu.root has rho_reco_barrel.GetNoElements()=" << rho_reco_barrel.GetNoElements() << "\n";
    return 0;
  }

  if (0) {
    std::vector<TH1F*> histos;
    std::vector<TVectorD*> vecs;
    std::vector<TVectorD*> vecErrs;
    std::vector<TString> names;
    std::vector<TString> labels;
    vecs.reserve(6); vecErrs.reserve(6); names.reserve(6); labels.reserve(6);
    vecs.push_back(&rho_reco_barrel); names.push_back("reco_barrel"); labels.push_back("reco");
    vecs.push_back(&rho_id_barrel); names.push_back("id_barrel"); labels.push_back("id");
    vecs.push_back(&rho_hlt_barrel); names.push_back("hlt_barrel"); labels.push_back("hlt");
    vecs.push_back(&rho_reco_endcap); names.push_back("reco_endcap"); labels.push_back("reco");
    vecs.push_back(&rho_id_endcap); names.push_back("id_endcap"); labels.push_back("id");
    vecs.push_back(&rho_hlt_endcap); names.push_back("hlt_endcap"); labels.push_back("hlt");

    vecErrs.push_back(&rho_reco_barrel_err);
    vecErrs.push_back(&rho_id_barrel_err);
    vecErrs.push_back(&rho_hlt_barrel_err); 
    vecErrs.push_back(&rho_reco_endcap_err);
    vecErrs.push_back(&rho_id_endcap_err); 
    vecErrs.push_back(&rho_hlt_endcap_err); 
    if (1) {
      // no errors
      ConvertVecsToHistos(vecs,DYTools::nPVBinCount,DYTools::nPVLimits, histos, names); 
    }
    else {
      // with errors
      ConvertVecsToHistos(vecs,vecErrs,DYTools::nPVBinCount,DYTools::nPVLimits, histos, names);
    }

    int color[3];
    int markers[3];
    color[0]=kBlack; color[1]=46; color[2]=kBlue+2;
    markers[0]=20; markers[1]=22; markers[2]=24;
    TCanvas *c= MakeCanvas("canvas_effScaleFactors","canvas_effScaleFactors", 1200, 600);
    c->Divide(2,1);
    CPlot cpB("cpBarrel","barrel", "PU","efficiency scale factor");
    CPlot cpE("cpEndcap","endcap", "PU","efficiency scale factor");
    for (unsigned int i=0; i<3; ++i) {
      PrepareHistoStyle(histos[i],color[i],markers[i]);
      cpB.AddHist1D(histos[i],labels[i],"LP",color[i],1,0,1);
    }
    for (unsigned int i=3; i<6; ++i) {
      PrepareHistoStyle(histos[i],color[i-3],markers[i-3]);
      cpE.AddHist1D(histos[i],labels[i],"LP",color[i-3],1,0,1);
    }
    cpB.SetLegend(0.2,0.2,0.45,0.4);
    cpE.SetLegend(0.2,0.2,0.45,0.4);
    cpB.Draw(c,false,".png",1);
    cpE.Draw(c,false,".png",2);
    c->Update();
    c->SaveAs("canvas_effScaleFactors.pdf");
    //return 0;
  }

  rho_barrel.ResizeTo(rho_id_barrel.GetNoElements());
  rho_endcap.ResizeTo(rho_id_endcap.GetNoElements());
  for (int i=0; i<rho_id_barrel.GetNoElements(); ++i) {
    rho_barrel[i] = rho_reco_barrel[i] * rho_id_barrel[i] * rho_hlt_barrel[i];
  }
  for (int i=0; i<rho_id_endcap.GetNoElements(); ++i) {
    rho_endcap[i] = rho_reco_endcap[i] * rho_id_endcap[i] * rho_hlt_endcap[i];
  }
  return 1;
}

// ---------------------------------------------------------------------

void MakePlots(const char *title, const char *name1, std::vector<TH1F*> &data1, const char *name2, std::vector<TH1F*> &data2, const std::vector<int> &color, int twoPlots, int plot_set) {
  TString canvasName=TString("canvas_") + TString(name1);
  TCanvas *c= MakeCanvas(canvasName,title, 600+600*twoPlots,600);
  c->Divide(1+twoPlots,1);
  char xlabel[50], ylabel[50], ylabel2[50];
  double ymin1=0, ymax1=1e5, ymin2=0, ymax2=1e5;

  //if (xAxisLumiIn_fm) sprintf(xlabel, "#int#font[12]{L}dt [fb^{-1}]"); else
  sprintf(xlabel, "#int#font[12]{L}dt [pb^{-1}]");

  //sprintf(xlabel2,"%s",xlabel);
  int scale_is_large=((rebinLuminosity==2) || ( TString(luminosityFileName) == TString("run_luminosity5.root") )) ? 1:0;
  if (rebinLuminosity==3) scale_is_large=2;
  switch(plot_set) {
  case 1:
    sprintf(ylabel, "#sigma");
    sprintf(ylabel2, "#font[12]{L_{block}}");
    ymin1=800; ymax1=1600;
    ymin2=100; ymax2=2*luminosityBlockSize;
    switch (scale_is_large) { 
    case 1: ymax2=3000; break;
    case 2: ymax2=1000; break;
    }
    break;
  case 2:
    sprintf(ylabel , "N_{obs}");
    sprintf(ylabel2, "N_{bkgr}/N_{obs} [%%]");
    ymin1=20000; ymax1=55000;
    ymin2=0; ymax2=1;
    switch (scale_is_large) {
    case 1: ymax1=600000; break;
    case 2: ymax1=200000; break;
    }
    TGaxis::SetMaxDigits(4);
    break;
  case 3:
    sprintf(ylabel , "<#epsilon>");
    sprintf(ylabel2, "<#rho>");
    ymin1=0.4; ymax1=0.7;
    ymin2=0.8; ymax2=1.2;
    //TGaxis::SetMaxDigits(3);
    break;
  case 4:
    sprintf(ylabel, "<PU>");
    //sprintf(ylabel2, "#sigma");
    ymin1=0; ymax1=30;
    ymin2=800; ymax2=1600;
    break;
  case 5:
    sprintf(ylabel, "<#epsilon#rho>");
    //sprintf(ylabel2, "#sigma");
    ymin1=0.3; ymax1=0.7;
    ymin2=800; ymax2=1600;
    break;
  default:
    sprintf(ylabel, "variable 1");
    sprintf(ylabel2, "variable 2");
    ymin1=0; ymax1=160000;
    ymin2=0; ymax2=5000;
  }
  if (xAxis_in_RR) {
    TGaxis::SetMaxDigits(6);
  }


  std::vector<TString> labels;
  if (swapLumiMethods) {
  labels.push_back("lumiCorrV2");
  labels.push_back("lumiNoCorr");
  }
  else {
    labels.push_back("lumiNoCorr");
    labels.push_back("lumiCorrV2");
  }
  labels.push_back("lumiCorrV3");
  labels.push_back("lumiCorrPix");

  //const char *drawOption="E";
  const char *drawOption="LP";

  if (1) {
  CPlot cp1(name1,"",xlabel,ylabel);  
  for (unsigned int i=0; i<data1.size(); ++i) {
    //std::cout << "data1[" << i << "]:\n"; PrintHisto(data1[i]);
    cp1.AddHist1D(data1[i],labels[i],drawOption,color[i+1],1,0,1);
  }
  cp1.SetYRange(ymin1-0.001,ymax1+0.001);
  cp1.Draw(c,false,".png",1);
  }

  if (twoPlots) {
  CPlot cp2(name2,"",xlabel,ylabel2);
  for (unsigned int i=0; i<data2.size(); ++i) {
    //std::cout << "data2[" << i << "]:\n"; PrintHisto(data2[i]);
    cp2.AddHist1D(data2[i],labels[i],drawOption,color[i+1],1,0,1);
  }
  cp2.SetYRange(ymin2-0.001,ymax2+0.001);
  cp2.Draw(c,false,".png",2);
  }
  c->Update();
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

int WriteToDir(TFile *fout, const TString &dirName, const std::vector<std::vector<TH1F*>*> &hVV) {
  if (hVV.size()!=4) {
    std::cout << "WriteToDir is foreseen for vec of 4 vectors\n";
    return 0;
  }

  std::vector<TString> labels;
  if (swapLumiMethods) {
  labels.push_back("lumiCorrV2");
  labels.push_back("lumiNoCorr");
  }
  else {
    labels.push_back("lumiNoCorr");
    labels.push_back("lumiCorrV2");
  }
  labels.push_back("lumiCorrV3");
  labels.push_back("lumiCorrPix");

  // start from the top dir
  fout->cd();
  // and create a subdirectory with a specified name
  TDirectory *dtop= fout->mkdir(dirName);
  dtop->cd();
  // loop over vector items
  for (unsigned int vi=0; vi<hVV.size(); ++vi) {
    // create a subdirectory for each version
    TDirectory *d= dtop->mkdir(labels[vi]);
    d->cd();
    // write the histograms
    const std::vector<TH1F*>* hV=hVV[vi];
    for (unsigned int i=0; i<hV->size(); ++i) {
      (*hV)[i]->Write();
    }
  }
  return 1;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
