#include <TROOT.h>
#include <TBenchmark.h>             // class to track macro running statistics
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TVectorD.h>
#include <TVector.h>
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

int LoadEffScaleFactors(const TString &fname, TVectorD &rho_recoId_barrel, TVectorD &rho_recoId_endcap, TVectorD &rho_hlt_barrel, TVectorD &rho_hlt_endcap);

void MakePlots(const char *title, const char *name1, std::vector<TH1F*> &data, const char *name2, std::vector<TH1F*> &data2, const std::vector<int> &colors);


// ---------------------------------------------------------------------

void lumiCrossSection(TString conf = "../config_files/data.conf") {

  gBenchmark->Start("lumiCrossSection");

  //--------------------------------------------------------------------------------------------------------------
  // Global constants and useful variables
  //==============================================================================================================


  char buf[200];
  const double factor_acceptance = 5598660/double(11428900);
  const double factor_fsr= 11428900/double(11840986);
  const double luminosityBlockSize=150;
  
  const int applySmearing=1;

  vector<int> colorV;
  vector<int> marker;

  colorV.reserve(5);
  colorV.push_back(kBlack); colorV.push_back(kBlack); 
  colorV.push_back(46); colorV.push_back(kGreen+2); colorV.push_back(kBlue+2);
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

  // load efficiency scale factors as a function of PU
  TVectorD rhoRecoIdBarrel, rhoRecoIdEndcap, rhoHltBarrel, rhoHltEndcap;
  if (!LoadEffScaleFactors(TString(constantsDir + TString("/rho_vs_pu.root")),rhoRecoIdBarrel,rhoRecoIdEndcap,rhoHltBarrel,rhoHltEndcap)) return;

  
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

    TH1F *hMCSignalTotal=new TH1F("hMCSignalTotal","",DYTools::nPVBinCount,DYTools::nPVLimits);
    hMCSignalTotal->Sumw2();
    for (int i=0; i<nEventsZPeakPURaw.GetNoElements(); ++i) hMCSignalTotal->Fill(i, nEventsZPeakPURaw[i]);
    for (int i=0; i<nPassZPeakPURaw.GetNoElements(); ++i)   hMCSignalEff->Fill(i, nPassZPeakPURaw[i]);
    //PrintHisto(hMCSignalTotal);  PrintHisto(hMCSignalEff);
    if (0) {
      TCanvas *cx= new TCanvas("ctestCount","ctestCount",600,600);
      hMCSignalTotal->DrawCopy("PE1");
      hMCSignalTotal->GetXaxis()->SetTitle("number of good vertices");
      hMCSignalTotal->GetYaxis()->SetTitle("count");
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
    hMCSignalEff->DrawCopy();
    cx->Update();
  }


  // Create PV distributions for MC
  for (unsigned int isam=0; isam<samplev.size(); ++isam) {
    sprintf(buf,"hPV_%s",snamev[isam].Data());
    TH1F *hPVData= new TH1F(buf,"",DYTools::nPVBinCount,DYTools::nPVLimits);
    hPVv.push_back(hPVData);
    const vector<SelectedEventData_t>* dt= (isam==0) ? &dataV : mcV[isam];
    if (!dt) { std::cout << "dt is null for isam=" << isam << "\n"; return ; }
    for (unsigned int i=0; i<dt->size(); ++i) {
      if ((*dt)[i].massInsideRange(60,120)) {
	hPVData->Fill((*dt)[i].nGoodPV, (*dt)[i].weight);
      }
    }
    //PrintHisto(hPVData);
  }


  std::vector<std::vector<double>*> mcZContribV;

  // histograms for each lumiCalc version
  std::vector<TH1F*> hLumiSigmaV;
  std::vector<TH1F*> hMCLumiSigmaV;
  std::vector<TH1F*> hLumiSigmaPerLumiV;
  std::vector<TH1F*> hMCLumiSigmaPerLumiV;

  for (int lumiVersion=1; lumiVersion<=4; ++lumiVersion) {
    //if (lumiVersion>2) break;
    int lumiSectionCount=0;
    double *lumiSections=NULL;
    std::vector<LumiInfo_t> lumiInfoV;
    if (lumiVersion==0) { 
      LumiInfo_t lumiTmp(0,200000, totalLuminosity); lumiInfoV.push_back(lumiTmp);  // catch all
    }
    else {
      TString lumiBlockInfoFile=constantsDir + TString("/run_luminosity.root");
      if (!PrepareLuminosity(lumiVersion,lumiInfoV,lumiBlockInfoFile,luminosityBlockSize,lumiSectionCount,&lumiSections)) {
	std::cout << "failed to prepare lumiInfoV for version=" << lumiVersion << "\n";
	return;
      }
    }

    sprintf(buf,"hSigma_%d",lumiVersion);
    TH1F *hSigma=new TH1F(buf,"",lumiSectionCount,lumiSections);
    PrepareHistoStyle(hSigma,colorV[lumiVersion],marker[lumiVersion]);
    hLumiSigmaV.push_back(hSigma);

    sprintf(buf,"hSigmaPerLumi_%d",lumiVersion);
    TH1F *hSigmaPerLumi=new TH1F(buf,"",lumiSectionCount,lumiSections);
    PrepareHistoStyle(hSigmaPerLumi,colorV[lumiVersion],marker[lumiVersion]);
    hLumiSigmaPerLumiV.push_back(hSigmaPerLumi);

    sprintf(buf,"hMCSigma_%d",lumiVersion);
    TH1F *hMCSigma=new TH1F(buf,"",lumiSectionCount,lumiSections);
    PrepareHistoStyle(hMCSigma,colorV[lumiVersion],marker[lumiVersion]);
    hMCLumiSigmaV.push_back(hMCSigma);

    sprintf(buf,"hMCSigmaPerLumi_%d",lumiVersion);
    TH1F *hMCSigmaPerLumi=new TH1F(buf,"",lumiSectionCount,lumiSections);
    PrepareHistoStyle(hMCSigmaPerLumi,colorV[lumiVersion],marker[lumiVersion]);
    hMCLumiSigmaPerLumiV.push_back(hMCSigmaPerLumi);

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
      
      // obtain PV distribution in data
      sprintf(buf,"hPVData_%u_%u",lumi.runNumMin,lumi.runNumMax);
      TH1F *hPVData= new TH1F(buf,buf,DYTools::nPVBinCount,DYTools::nPVLimits);
      double dataZCount=0, bkgrZCount=0;
      double sumEff=0;  // sum of efficiencies for the average
      double sumRho=0; // sum of eff. scale factors for the average
      for (unsigned int i=0; i<dataV.size(); ++i) {
	if (lumi.insideRange(dataV[i].runNum)) {
	  hPVData->Fill(dataV[i].nGoodPV, dataV[i].weight);
	  if (dataV[i].massInsideRange(60,120) && 
	      dataV[i].nGoodPV)  // at least 1 good vertex
	    {
	    dataZCount++;                          // count Z candidates in data
	    int idx=hMCSignalEff->FindBin(dataV[i].nGoodPV);
	    if (idx==0) std::cout << "idx=0 for goodPV=" << dataV[i].nGoodPV << "\n";
	    if ((idx==-1) || (idx>hMCSignalEff->GetNbinsX())) {
	      std::cout << "got bin idx=" << idx << " from hMCSignalEff. dataV[i].nGoodPV=" << dataV[i].nGoodPV << "\n";
	      return;
	    }
	    sumEff += hMCSignalEff->GetBinContent(idx); // add-up the efficiency for averaging
	    idx--; // correction for TVector indexing
	    // first has to be the leading electron!
	    double rho_factor= (dataV[i].firstIsInBarrel()) ? 
	      rhoRecoIdBarrel[idx]*rhoHltBarrel[idx] : rhoRecoIdEndcap[idx]*rhoHltEndcap[idx];
	    rho_factor *= (dataV[i].secondIsInBarrel()) ? 
	      rhoRecoIdBarrel[idx]*rhoHltBarrel[idx] : rhoRecoIdEndcap[idx]*rhoHltEndcap[idx];
	    sumRho += rho_factor;
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
	    }
	  }
	}
	zContrib->push_back(Nbkgr);
	if ( snamev[isam] != "zee" ) bkgrZCount+=Nbkgr;
	else mcSignalZv.push_back(Nbkgr);
      }
      mcBkgrZv.push_back(bkgrZCount);

      double Nsignal=dataZCount-bkgrZCount;
      double avgEff=sumEff/Nsignal;
      double avgRho=sumRho/Nsignal;
      double sigma=Nsignal/( lumi.lumiWeight * avgEff * avgRho * factor_acceptance * factor_fsr );
      double sigmaMC= mcSignalZv.back()/ ( lumi.lumiWeight * avgEff * avgRho * factor_acceptance * factor_fsr );
      hSigma->Fill( sumLumi + 0.5* lumi.lumiWeight, sigma );
      hMCSigma->Fill( sumLumi + 0.5 * lumi.lumiWeight, sigmaMC );
      hSigmaPerLumi->Fill( sumLumi + 0.5* lumi.lumiWeight, sigma/lumi.lumiWeight );
      hMCSigmaPerLumi->Fill( sumLumi + 0.5 * lumi.lumiWeight, sigmaMC/lumi.lumiWeight );
      sumLumi += lumi.lumiWeight;
      std::cout << "sumLumi=" << sumLumi << ", dataZCount=" << dataZCount << ", bkgrZCount=" << bkgrZCount << ", avgEff=" << avgEff << ", sigma=" << sigma << ", sigmaMC=" << sigmaMC << "\n";
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
    std::cout << "Z counts per lumi\n";
    PrintHisto(hSigmaPerLumi,hMCSigmaPerLumi);
  }

  MakePlots("Sigma","Sigma", hLumiSigmaV, "SigmaPerLumi",hLumiSigmaPerLumiV, colorV);

  gBenchmark->Show("lumiCrossSection");
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
  luminosity.clear();
  lumiSectionCount=0;
  if (*lumiSections) delete *lumiSections;

  TFile f(fname);
  if (!f.IsOpen()) {
    std::cout << "failed to open a file <" << fname << ">\n";
    return 0;
  }
  TVectorD runNumMin,runNumMax, lumiWeights;
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
  if ((n!=runNumMax.GetNoElements()) || (n!=lumiWeights.GetNoElements())) {
    std::cout << "count mismatch : runNumMin[" << runNumMin.GetNoElements() << "], runNumMax[" << runNumMax.GetNoElements() << "], lumiWeight[" << lumiWeights.GetNoElements() << "]\n";
    return 0;
  }

  if (1) {
    std::cout << "\nRaw luminositied :\n";
    for (int i=0; i<runNumMin.GetNoElements(); ++i) {
      std::cout << " run " << runNumMin[i] << " -- " << runNumMax[i] << "    "  << lumiWeights[i] << "\n";
    }
    std::cout << "\n";
  }
  
  const double tolerance=0.1*dLumi;
  double lumiChunk=0;
  double lumi=0;
  LumiInfo_t info;
  int start=1;
  std::vector<double> lumiBins;
  lumiBins.push_back(0.);
  lumiChunk=0;
  for (int i=0; i<n; ++i) {
    if (start) { info.runNumMin=UInt_t(runNumMin[i]); start=0; }
    double lumiChunk1 = lumiChunk + lumiWeights[i];
    int nextChunkIsLarge=((i<n-1) && (lumiWeights[i+1]>dLumi)) ? 1:0;
    if (( lumiChunk1 > dLumi) || ( lumiChunk1 > dLumi-tolerance) || nextChunkIsLarge ) {
      info.runNumMax=UInt_t(runNumMax[i]); start=1; lumiChunk=0.;
      info.lumiWeight=lumiChunk1;
      luminosity.push_back(info);
      lumi += lumiChunk1;
      lumiBins.push_back(lumi);
    }
    else lumiChunk=lumiChunk1;
  }
  if (start!=1) {
    info.runNumMax=UInt_t(runNumMax[n-1]);
    info.lumiWeight=lumiChunk;
    luminosity.push_back(info);
    lumi += lumiChunk;
    lumiBins.push_back(lumi);
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

int LoadEffScaleFactors(const TString &fname, TVectorD &rho_recoId_barrel, TVectorD &rho_recoId_endcap, TVectorD &rho_hlt_barrel, TVectorD &rho_hlt_endcap) {
  //int res=1;

  TFile fEffScaleFactors(fname);
  TVector pu_bin_min = *(TVector*)fEffScaleFactors.Get("pu_bin_min");
  TVector pu_bin_max = *(TVector*)fEffScaleFactors.Get("pu_bin_max");
  TVectorD rho_reco_barrel = *(TVectorD*)fEffScaleFactors.Get("rho_reco_barrel");
  TVectorD rho_reco_endcap = *(TVectorD*)fEffScaleFactors.Get("rho_reco_endcap");
  TVectorD rho_id_barrel = *(TVectorD*)fEffScaleFactors.Get("rho_id_barrel");
  TVectorD rho_id_endcap = *(TVectorD*)fEffScaleFactors.Get("rho_id_endcap");
  rho_hlt_barrel.ResizeTo(rho_id_barrel.GetNoElements());
  rho_hlt_endcap.ResizeTo(rho_id_endcap.GetNoElements());
  rho_hlt_barrel = *(TVectorD*)fEffScaleFactors.Get("rho_hlt_barrel");
  rho_hlt_endcap = *(TVectorD*)fEffScaleFactors.Get("rho_hlt_endcap");

  //TVectorD rho_reco_barrel_err = *(TVectorD*)fEffScaleFactors.Get("rho_reco_barrel_Err");
  //TVectorD rho_reco_endcap_err = *(TVectorD*)fEffScaleFactors.Get("rho_reco_endcap_Err");
  //TVectorD rho_id_barrel_err = *(TVectorD*)fEffScaleFactors.Get("rho_id_barrel_Err");
  //TVectorD rho_id_endcap_err = *(TVectorD*)fEffScaleFactors.Get("rho_id_endcap_Err");
  //TVectorD rho_hlt_barrel_err = *(TVectorD*)fEffScaleFactors.Get("rho_hlt_barrel_Err");
  //TVectorD rho_hlt_endcap_err = *(TVectorD*)fEffScaleFactors.Get("rho_hlt_endcap_Err");

  if (1) {
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

  rho_recoId_barrel.ResizeTo(rho_reco_barrel);
  for (int i=0; i<rho_id_barrel.GetNoElements(); ++i) {
    rho_recoId_barrel[i] = rho_reco_barrel[i] * rho_id_barrel[i];
  }
  rho_recoId_endcap.ResizeTo(rho_reco_endcap);
  for (int i=0; i<rho_id_endcap.GetNoElements(); ++i) {
    rho_recoId_endcap[i] = rho_reco_endcap[i] * rho_id_endcap[i];
  }
  return 1;
}

// ---------------------------------------------------------------------

void MakePlots(const char *title, const char *name1, std::vector<TH1F*> &data1, const char *name2, std::vector<TH1F*> &data2, const std::vector<int> &color) {
  TString canvasName=TString("canvas_") + TString(name1);
  int twoPlots=0;
  TCanvas *c= MakeCanvas(canvasName,title, 600+600*twoPlots,600);
  c->Divide(1+twoPlots,1);
  char xlabel[50], ylabel[50], ylabel2[50];
  sprintf(xlabel, "#int#font[12]{L}dt [pb^{-1}]");
  sprintf(ylabel, "#sigma");
  sprintf(ylabel2, "#sigma/LumiBin");

  std::vector<TString> labels;
  labels.push_back("lumiNoCorr");
  labels.push_back("lumiCorrV2");
  labels.push_back("lumiCorrV3");
  labels.push_back("lumiCorrPix");

  if (1) {
  CPlot cp1(name1,"",xlabel,ylabel);  
  for (unsigned int i=0; i<data1.size(); ++i) {
    //std::cout << "data1[" << i << "]:\n"; PrintHisto(data1[i]);
    cp1.AddHist1D(data1[i],labels[i],"E",color[i+1]);
  }
  //cp1.SetYRange(800,1500);
  cp1.SetYRange(800,2500.001);
  cp1.Draw(c,false,".png",1);
  }

  if (twoPlots) {
  CPlot cp2(name2,"",xlabel,ylabel2);
  for (unsigned int i=0; i<data2.size(); ++i) {
    //std::cout << "data2[" << i << "]:\n"; PrintHisto(data2[i]);
    cp2.AddHist1D(data2[i],labels[i],"E",color[i+1]);
  }
  cp2.SetYRange(0,30);
  cp2.Draw(c,false,".png",2);
  }
  c->Update();
}

// ---------------------------------------------------------------------
