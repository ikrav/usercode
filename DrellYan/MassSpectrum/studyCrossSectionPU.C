#include <TROOT.h>
#include <TBenchmark.h>             // class to track macro running statistics
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TVectorD.h>
#include <TVector.h>
#include <TGaxis.h>
#include <TKey.h>
#include <vector>                   // STL vector class
#include <iostream>                 // dard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings
#include <assert.h>
using namespace std;
#include "../Include/CSample.hh"
#include "../Include/ZeeData.hh"
#include "../Include/ElectronEnergyScale.hh"
#include "../Include/DYTools.hh"
#include "../MassSpectrum/storeData.hh"
#include "../Include/MyTools.hh"
#include "../Include/MitStyleRemix.hh"
#include "../Include/CPlot.hh"

// ---------------------------------------------------------------------

const double luminosityBlockSize=150;
const int rebinLuminosity=0;  // set to 1 for run_luminosity.root
const char *luminosityFileName="run_luminosity5.root"; // path will be prepended   

//  A correction to make largest luminosity version
//  to appear first in the list
const int swapLumiMethods=1; 


const int xAxis_in_RR=0;

typedef enum { _PS_CrossSection, _PS_NormCrossSection } TPlotSetKind_t;

// ---------------------------------------------------------------------
void PrepareHistoStyle(TH1F* h, int color, int markerStyle=20);
TH1F *CopyHisto(const TString &name, TH1F *h);

void MakePlots(const char *title, const char *name1, std::vector<TH1F*> &data, const char *name2, std::vector<TH1F*> &data2, const std::vector<int> &colors, int plotTwo, int plot_set);
void MakePlots5(const char *title, std::vector<std::vector<TH1F*>*> &data, const std::vector<TString> &labels, const std::vector<int> &colors, const std::vector<int> &markers, TPlotSetKind_t plotSet);

int ReadFromDir(TFile *fout, const TString &dirName, std::vector<std::vector<TH1F*>*> &hVV, std::vector<TString> &names);


// ---------------------------------------------------------------------

void studyCrossSectionPU(TString plotFile = "./dataLumiCrossSectionPU.root") {

  gBenchmark->Start("studyCrossSectionPU");

  vector<int> colors;
  vector<int> markers;
  colors.reserve(5);
  colors.push_back(kBlack); colors.push_back(38); 
  colors.push_back(46); colors.push_back(kGreen+2); colors.push_back(kBlue+1);
  //colors.push_back(k
  markers.reserve(5);
  markers.push_back(27); 
  markers.push_back(20); markers.push_back(24); 
  markers.push_back(22); markers.push_back(26);

  //--------------------------------------------------------------------------------------------------------------
  // Work
  //==============================================================================================================

  TFile *fin=new TFile(plotFile);
  if (!fin) {
    std::cout << "failed to open the file <" << plotFile << ">\n";
    return;
  }

  std::vector<TString> labels;
  std::vector<std::vector<TH1F*>*> hLumiSigmaPUVV, hLumiSigmaNormPUVV;
  ReadFromDir(fin,"CrossSection",hLumiSigmaPUVV,labels);
  ReadFromDir(fin,"NormCrossSection",hLumiSigmaNormPUVV,labels);
  fin->Close();

  MakePlots5("CrossSectionPU",hLumiSigmaPUVV,labels,colors,markers,_PS_CrossSection);
  MakePlots5("NormCrossSectionPU",hLumiSigmaNormPUVV,labels,colors,markers,_PS_NormCrossSection);

  gBenchmark->Show("studyCrossSectionPU");
  return;
}


// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void PrepareHistoStyle(TH1F* h, int color, int markerStyle) {
  h->GetXaxis()->SetTitleOffset(1.10);
  h->GetYaxis()->SetTitleOffset(1.40);
  h->SetLineColor(color);
  h->SetMarkerColor(color);
  h->SetMarkerStyle(markerStyle);
}

// ---------------------------------------------------------------------

TH1F* CopyHisto(const TString &name, TH1F *h) {
  const TArrayD *bins=h->GetXaxis()->GetXbins();
  double *b=new double[bins->GetSize()];
  for (int i=0; i<bins->GetSize(); ++i) {
    b[i]=(*bins)[i];
  }
  TH1F *histo=new TH1F(name,name,h->GetNbinsX(),b);
  for (int i=1; i<h->GetNbinsX(); ++i) {
    histo->SetBinContent(i, h->GetBinContent(i));
    histo->SetBinError(i, h->GetBinError(i));
  }
  delete b;
  return histo;
}

// ---------------------------------------------------------------------

int ReadFromDir(TFile *fout, const TString &dirName, std::vector<std::vector<TH1F*>*> &hVV, std::vector<TString> &names) {
  hVV.clear(); names.clear();

  fout->cd(dirName);
  TDirectory *sourceDir = gDirectory;
  TList *subDirList=sourceDir->GetListOfKeys();
  int count=subDirList->GetEntries();
  std::cout << "there are " << count << " keys\n";
  hVV.reserve(count);
  TIter subDirs(subDirList);
  TKey *key=NULL;
  while ( (key = (TKey*)subDirs()) ) {
    const TString sDirName=key->GetName();
    std::cout << "sDirName=" << sDirName << "\n";
    names.push_back( sDirName );
    sourceDir->cd( sDirName );
    TDirectory *currDir = gDirectory;
    TList *histosList=currDir->GetListOfKeys();
    TIter histos( histosList );
    std::vector<TH1F*> *hV = new std::vector<TH1F*>();
    hVV.push_back(hV);
    hV->reserve(histosList->GetEntries());
    while( (key = (TKey*)histos()) ) {
      std::cout << "key->GetName() = " << key->GetName() << "\n";
      TObject *obj= key->ReadObj();
      if ( obj->IsA()->InheritsFrom( TH1F::Class() ) ) {
	TH1F *h=(TH1F*)obj;
	h->SetDirectory(0);
	hV->push_back(h);
      }
      else {
	delete hV;
      }
    }
  }
  return 1;
}

// ---------------------------------------------------------------------
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

void MakePlots5(const char *title, std::vector<std::vector<TH1F*>*> &data, const std::vector<TString> &labels, const std::vector<int> &colors, const std::vector<int> &markers, TPlotSetKind_t plotSet) {
  
  std::cout << "title=" << title << ", data.size=" << data.size();
  if (data.size()) std::cout << ", data[0]->size()=" << data[0]->size(); 
  std::cout << "\n";

  char xlabel[50], ylabel[50];
  double ymin=0, ymax=1e4;

  //if (xAxisLumiIn_fm) sprintf(xlabel, "#int#font[12]{L}dt [fb^{-1}]"); else
  sprintf(xlabel, "nGoodPV");

  switch(plotSet) {
  case _PS_CrossSection:
    sprintf(ylabel, "#sigma per PU bin");
    ymin=0; ymax=500;
    break;
  case _PS_NormCrossSection:
    sprintf(ylabel, "N_{obs} #sigma / N_{obs,pu}");
    ymin=0; ymax=2000;
    TGaxis::SetMaxDigits(3);
    break;
  default:
    sprintf(ylabel, "variable 1");
  }
  if (xAxis_in_RR) {
    TGaxis::SetMaxDigits(6);
  }

  std::vector<CPlot*> cplots;
  std::vector<TString> names;
  names.push_back("2011A_May10");
  names.push_back("2011A_PromptReco");
  names.push_back("2011A_05Aug");
  names.push_back("2011A_03Oct");
  names.push_back("2011B");

  std::vector<TString> letter;
  letter.push_back("(a)");
  letter.push_back("(b)");
  letter.push_back("(c)");
  letter.push_back("(d)");
  letter.push_back("(e)");
  letter.push_back("(f)");

  TString drawOption="LP";

  TString canvasName=TString("canvas_") + TString(title);
  int baseSize=250;
  TCanvas *canvas= MakeCanvas(canvasName,title, 2*baseSize,3*baseSize);
  canvas->Divide(2,3);

  for (unsigned int ci=0; ci<data[0]->size(); ++ci) {
    TString cName="canvas_" + names[ci];
    CPlot *cp=new CPlot(cName,"",xlabel,ylabel);
    cp->SetYRange(ymin-0.001,ymax+0.001);
    cplots.push_back(cp);
    for (unsigned int i=0; i<data.size(); ++i) {
      TH1F* h=(*data[i])[ci];
      PrepareHistoStyle(h,colors[i],markers[i]);
      cp->AddHist1D(h,labels[i],drawOption,colors[i],1,0,1);
      if (plotSet==_PS_NormCrossSection) cp->TransLegend(0.0,-0.12);
    }
    cp->AddTextBox(names[ci],0.25,0.92,0.75,0.99, 0,kBlack,kWhite);
    cp->AddTextBox(letter[ci],0.18,0.82,0.28,0.88, 0,kBlack,kWhite);
    cp->Draw(canvas,false,".png",ci+1);
  }

  if (1) {
  int idx=3;
  TString cTName="canvas_"+labels[idx];
  drawOption="LP";
  CPlot *cTime=new CPlot(cTName,"",xlabel,ylabel);
  cTime->SetYRange(ymin-0.001,ymax+0.001);
  cplots.push_back(cTime);
  for (unsigned int i=0; i<data[idx]->size(); ++i) {
    //if (( plotSet == _PS_CrossSection ) && (i%2==1)) continue;
    if (i%2==1) continue;
    TH1F* ho=(*data[idx])[i];
    TString tmp_name=ho->GetName() + TString("_tmp");
    //TH1F* h=CopyHisto(tmp_name,ho);
    TH1F* h=(TH1F*)(*data[idx])[i]->Clone(tmp_name);
    PrepareHistoStyle(h,colors[i],markers[i]);
    cTime->AddHist1D(h,names[i],drawOption,colors[i],1,0,1);
    if (plotSet==_PS_NormCrossSection) cTime->TransLegend(0.0,-0.15);
  }
  TString text=TString("method=") + labels[idx];
  cTime->AddTextBox(text,0.25,0.92,0.75,0.99, 0,kBlack,kWhite);
  cTime->AddTextBox(letter[5],0.18,0.82,0.28,0.88, 0,kBlack,kWhite);
  cTime->Draw(canvas,false,".png",6);
  }
  canvas->Update();
  return ;
}


// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
