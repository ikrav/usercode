
//#include "../Include/ElectronEnergyScale.hh"
#include "../ElectronEnergyScaleAdv.C"
#include "../Include/MitStyleRemix.hh"
#include "../Include/CPlot.hh"
//#include "../Include/MyTools.hh"
#include <math.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TText.h>
#include <TGaxis.h>

#include "fileLocation.hh"


void DrawShadedText(TText *label, Double_t x, Double_t y, const char *text, int color=kBlack) {
  label->SetTextColor(kGray); label->DrawTextNDC(x+0.005,y-0.005,text);
  label->SetTextColor(color); label->DrawTextNDC(x,y,text);
}

typedef enum { _EESF_None=0, _EESF_Data=1, _EESF_MC=2 } TEESFactorType_t;
typedef enum { _EBSet_None=0, _EBSet_BB=1, _EBSet_BE=2, _EBSet_EE=4, _EBSet_All=8 } TEndcapBarrelSet_t;
typedef enum { _Landscape=0, _Portrait=1 } TPageOrientation_t;

// -------------------------------------------------------------

struct TH1FQuartet_t {
  TString FName;
  TH1F *FAll, *FBB, *FBE, *FEE;
  TH1FQuartet_t(TString name="") : FName(name), FAll(NULL), FBB(NULL), FBE(NULL), FEE(NULL) {}

  const TString& Name() const { return FName; }
  void Name(const TString &name) { FName=name; }

  // ---------------------------------------------

  void Clear() { 
    FName.Clear();
    if (FAll) { delete FAll; FAll=0; }
    if (FBB)  { delete FBB;  FBB=0; }
    if (FBE)  { delete FBE;  FBE=0; }
    if (FEE)  { delete FEE;  FEE=0; }
  }
  
  // ---------------------------------------------

  int Create(const TString &name, const TString &hnameBase, int numBins, double xmin, double xmax) {
    TString hname;
    Clear();
    FName=name;
    hname=hnameBase; hname.Append("_all");
    FAll=new TH1F(hname.Data(),hname.Data(),numBins,xmin,xmax);
    hname=hnameBase; hname.Append("_BB");
    FBB =new TH1F(hname.Data(),hname.Data(),numBins,xmin,xmax);
    hname=hnameBase; hname.Append("_BE");
    FBE =new TH1F(hname.Data(),hname.Data(),numBins,xmin,xmax);
    hname=hnameBase; hname.Append("_EE");
    FEE =new TH1F(hname.Data(),hname.Data(),numBins,xmin,xmax);

    FAll->Sumw2(); FBB->Sumw2(); FBE->Sumw2(); FEE->Sumw2();
    return 1;
  }

  // ---------------------------------------------

  void SetMarker(int markerStyle, double markerSize=1.) {
    FAll->SetMarkerStyle(markerStyle); FAll->SetMarkerSize(markerSize);
    FBB->SetMarkerStyle(markerStyle); FBB->SetMarkerSize(markerSize);
    FBE->SetMarkerStyle(markerStyle); FBE->SetMarkerSize(markerSize);
    FEE->SetMarkerStyle(markerStyle); FEE->SetMarkerSize(markerSize);
  }

  // ---------------------------------------------

  void SetColor(int color) {
    FAll->SetLineColor(color); FAll->SetMarkerColor(color);
    FBB->SetLineColor(color); FBB->SetMarkerColor(color);
    FBE->SetLineColor(color); FBE->SetMarkerColor(color);
    FEE->SetLineColor(color); FEE->SetMarkerColor(color);
  }

  // ---------------------------------------------

  void SetFillColor(int color) {
    FAll->SetFillColor(color); FAll->SetFillStyle(1001);
    FBB->SetFillColor(color); FBB->SetFillStyle(1001);
    FBE->SetFillColor(color); FBE->SetFillStyle(1001);
    FEE->SetFillColor(color); FEE->SetFillStyle(1001);
  }

  // ---------------------------------------------
  
  void SetXTitleOffset(double offset) {
    FAll->GetXaxis()->SetTitleOffset(offset);
    FBB->GetXaxis()->SetTitleOffset(offset);
    FBE->GetXaxis()->SetTitleOffset(offset);
    FEE->GetXaxis()->SetTitleOffset(offset);
  }
  // ---------------------------------------------
  
  void SetYTitleOffset(double offset, int set=15) {
    if (set&8) FAll->GetYaxis()->SetTitleOffset(offset);
    if (set&1) FBB->GetYaxis()->SetTitleOffset(offset);
    if (set&2) FBE->GetYaxis()->SetTitleOffset(offset);
    if (set&4) FEE->GetYaxis()->SetTitleOffset(offset);
  }

  // ---------------------------------------------

  template<class ESCorr_t>
  int Fill(const EtaEtaMassData_t &eem, const ESCorr_t &esf, TEESFactorType_t applyEES) {
    double mass=0;
    switch(applyEES) {
    case _EESF_None: mass=eem.mass(); break;
    case _EESF_Data: mass=eem.mass()*sqrt(esf.getEnergyScaleCorrection(eem.eta1()) * esf.getEnergyScaleCorrection(eem.eta2())); break;
    case _EESF_MC: mass=eem.mass() + esf.generateMCSmear(eem.eta1(),eem.eta2()); break;
    default:
      std::cout << "unprepared applyEES=" << applyEES << " in TH1FQuartet_t::Fill\n";
      throw 2;
    }
    int isEB1=(fabs(eem.eta1())<1.5) ? 1:0;
    int isEB2=(fabs(eem.eta2())<1.5) ? 1:0;
    FAll->Fill(mass);
    switch(isEB1+isEB2) {
    case 0: FEE->Fill(mass); break;
    case 1: FBE->Fill(mass); break;
    case 2: FBB->Fill(mass); break;
    default:
      std::cout << "incorrect isEB1+isEB2=" << (isEB1+isEB2) << " in TH1FQuartet_t::Fill\n";
      throw 2;
    }
    return 1;
  }

  // ---------------------------------------------

  template<class ESCorr_t>
  int Fill(const std::vector<EtaEtaMassData_t> &eemV, const ESCorr_t &esf, TEESFactorType_t applyEES) {
    int res=1;
    for (unsigned int i=0; res && (i<eemV.size()); ++i) {
      res=this->Fill(eemV[i],esf,applyEES);
    }
    if (!res) std::cout << "error in TH1FQuartet_t::Fill(vector)\n";
    return res;
  }

  // ---------------------------------------------

  int Plot(CPlot &cp, TEndcapBarrelSet_t theSet, const TEESFactorType_t marks) const {
    const TH1F *hs;
    TString explLabel;
    switch(theSet) {
    case _EBSet_BB: hs=FBB; explLabel="barrel-barrel"; break;
    case _EBSet_BE: hs=FBE; explLabel="barrel-endcap"; break;
    case _EBSet_EE: hs=FEE; explLabel="endcap-endcap"; break;
    case _EBSet_All: hs=FAll; explLabel="all regions"; break;
    case _EBSet_None:
    default:
      std::cout << "TH1FQuartet::Plot cannot handle theSet=" << theSet << "\n";
      return 0;
    }
    TString hname=hs->GetName();
    hname.Append("1");
    TH1F *h=(TH1F*) hs->Clone(hname);
    TString label, options;
    int color;
    switch(marks) {
    case _EESF_Data: 
      label="data"; options="E"; color=kBlack; 
      cp.AddHist1D(h, label, options, color);
      break;
    case _EESF_MC: 
      label="MC"; options="hist "; color=426; 
      cp.AddHist1D(h, label, options, color, 1, 1);
      break;
    default:
      std::cout << "TH1FQuartet::Plot cannot handle marks=" << marks << "(TEESFActorType_t)\n";
      return 0;
    }
    cp.AddTextBox(explLabel, 0.65,0.5,0.949,0.55, 0,kBlack,kWhite);
    return 1;
  }

  // ---------------------------------------------

  void Normalize(const TH1FQuartet_t &data) {
    double nAll0=(data.FAll) ? data.FAll->Integral() : 1;
    if (fabs(nAll0)<1e-6) nAll0=1;
    double nBB0=(data.FBB) ? data.FBB->Integral() : 1;
    if (fabs(nBB0) <1e-6) nBB0=1;
    double nBE0=(data.FBE) ? data.FBE->Integral() : 1;
    if (fabs(nBE0) <1e-6) nBE0=1;
    double nEE0=(data.FEE) ? data.FEE->Integral() : 1;
    if (fabs(nEE0) <1e-6) nEE0=1;

    double nAll=(FAll) ? FAll->Integral() : 1;
    if (fabs(nAll)<1e-6) nAll=1;
    double nBB=(FBB) ? FBB->Integral() : 1;
    if (fabs(nBB) <1e-6) nBB=1;
    double nBE=(FBE) ? FBE->Integral() : 1;
    if (fabs(nBE) <1e-6) nBE=1;
    double nEE=(FEE) ? FEE->Integral() : 1;
    if (fabs(nEE) <1e-6) nEE=1;
    
    if (FAll) FAll->Scale(nAll0/nAll);
    if (FBB)  FBB->Scale(nBB0/nBB);
    if (FBE)  FBE->Scale(nBE0/nBE);
    if (FEE)  FEE->Scale(nEE0/nEE);
    return;
  }

  // ---------------------------------------------

  void PrintIntegrals() const {
    std::cout << "TH1FQuartet(" << FName << ") integrals:\n";
    std::cout << "  all           :" << ((FAll) ? FAll->Integral() : 0.) << "\n";
    std::cout << "  barrel-barrel :" << ((FBB) ? FBB->Integral() : 0.) << "\n";
    std::cout << "  barrel-endcap :" << ((FBE) ? FBE->Integral() : 0.) << "\n";
    std::cout << "  endcap-endcap :" << ((FEE) ? FEE->Integral() : 0.) << "\n";
  }

};

// -------------------------------------------------------------
// -------------------------------------------------------------

struct Distributions_t {
  TString FName;
  TH1FQuartet_t FData,FMC;

  Distributions_t(TString name="") : FName(name), FData(),FMC() {}

  const TString& Name() const { return FName; }
  void Name(const TString &name) { FName=name; }
  
  // ---------------------------------------------

  int Create(const TString &name, const TString &distrNameBase, int numBins, double xmin, double xmax) {
    FName=name;
    TString nameData=name;    nameData.Append("_Data");
    TString hnameData=distrNameBase;    hnameData.Insert(1,"Data");
    TString nameMC=name;   nameMC.Append("_MC");
    TString hnameMC=distrNameBase;    hnameMC.Insert(1,"MC");
    int res=
      FData.Create(nameData,hnameData,numBins,xmin,xmax) &&
      FMC.Create(nameMC,hnameMC,numBins,xmin,xmax);
    if (res) {
      FMC.SetMarker(1); // dot
      if (1) {
	FData.SetMarker(24); // empty circle
      }
      else {
	FData.SetMarker(20,0.9);
      }
      FMC.SetColor(426); FMC.SetFillColor(426);
      FData.SetYTitleOffset(1.4);
      FMC.SetYTitleOffset(1.4);
      FData.SetYTitleOffset(1.5,int(_EBSet_All));
      FMC.SetYTitleOffset(1.5,int(_EBSet_All));
    }
    if (!res) std::cout << "error in Distributions_t::Create\n";
    return res;
  }

  // ---------------------------------------------

  template<class ESCorr_t>
  int Fill(const std::vector<EtaEtaMassData_t> &eemV, const ESCorr_t &esf, TEESFactorType_t dataType, TEESFactorType_t applyEES) {
    int res=1;
    if (applyEES!=_EESF_None) applyEES=dataType;
    switch(dataType) {
    case _EESF_Data: res=FData.Fill(eemV,esf,applyEES); break;
    case _EESF_MC:   res=FMC.Fill(eemV,esf,applyEES);   break;
    default:
      std::cout << "wrong dataType=" << dataType << " in Distributions_t::Fill\n";
      res=0;
    }
    return res;
  }

  // ---------------------------------------------

  template<class ESCorr_t>
  int Fill(const std::vector<std::vector<EtaEtaMassData_t>*> &eemV, const ESCorr_t &esf, TEESFactorType_t dataType, TEESFactorType_t applyEES) {
    int res=1;
    for (unsigned int i=0; res && (i<eemV.size()); ++i) {
      res=this->Fill(*eemV[i],esf,dataType,applyEES);
    }
    if (!res) std::cout << "error in Distributions_t::Fill(v<v>)\n";
    return res;
  }

  // ---------------------------------------------

  template<class ESCorr_t>
  int Fill(const std::vector<std::vector<EtaEtaMassData_t>*> &eemDataV, const std::vector<std::vector<EtaEtaMassData_t>*> &eemMCV, const ESCorr_t &esf, TEESFactorType_t applyEES) {
    int res=
      this->Fill(eemDataV,esf,_EESF_Data,applyEES) &&
      this->Fill(eemMCV  ,esf,_EESF_MC  ,applyEES);
    if (!res) std::cout << "error in Distributions_t::Fill(dataV,mcV)\n";
    else {
      this->PrintIntegrals();
    }
    return res;
  }

  // ---------------------------------------------

  int Plot(CPlot &cp, TEndcapBarrelSet_t theSet) const {
    int res=
      FMC.Plot(cp,theSet,_EESF_MC) &&
      FData.Plot(cp,theSet,_EESF_Data);


    //TCanvas c("c","c",300,300);
    //cp.SetYRange(0,1e4);
    //cp.Draw(&c,false,"png");
    //c.Update();
    //char x; std::cin >> x; std::cin.clear(); std::cin.sync();

    //int res=
    // FData.Plot(cp,theSet,_EESF_Data);
    //FMC.Plot(cp,theSet,_EESF_MC);
    if (!res) std::cout << "error in Distributions_t::Plot(cp,theSet=" << theSet << ")\n";
    return res;
  }

  // ---------------------------------------------

  int CreatePlots(std::vector<CPlot*> &cpV, const TString &cplotTitle="", const TString &xlabel="mass [GeV/c^{2}]", const TString &ylabel="counts", const std::vector<double> *ymax=NULL) const {
    ClearVec(cpV);
    cpV.reserve(4);
    int res=1;
    for (unsigned int i=0; res && (i<4); ++i) {
      CPlot *cp=new CPlot(cplotTitle,"",xlabel,ylabel);
      std::cout << "i=" << i << ". cp is null? " << ((cp) ? "no" : "yes") << std::endl;
      //std::cout << "is legend null? " << ((cp->GetLegend()) ? "no" : "yes") << std::endl;
      cpV.push_back(cp);
      TEndcapBarrelSet_t theSet=TEndcapBarrelSet_t(1<<i);
      std::cout << "theSet=" << theSet << "\n";
      res=this->Plot(*cp,theSet);

      if (ymax && (i<ymax->size()))  cp->SetYRange(0,(*ymax)[i]);
      else cp->SetYRange(0,1e5);
    }
    if (!res) std::cout << "error in Distributions_t::Plot(cpV)\n";
    return res;
  }

  // ---------------------------------------------

  void PrintIntegrals() const {
    std::cout << "Distributions(" << FName << ") :\n";
    FData.PrintIntegrals();
    FMC.PrintIntegrals();
  }

  // ---------------------------------------------

  void NormalizeMCtoData() {
    FMC.Normalize(FData);
  }

  // ---------------------------------------------
};


// -------------------------------------------------------------




// --------------------------------------------

void eesCorrectionDemo(const TString &electronEnergyScaleString) {
  std::cout << "\nElectron energy scale correction demo for eesString=<" << electronEnergyScaleString << ">\n";
  TString eesString=electronEnergyScaleString;

  // Some globals
  CPlot::sOutDir = "plots";
  TGaxis::SetMaxDigits(3); // show no more than 3 digits on axes
  TPageOrientation_t pageOrientation= (0) ? _Portrait : _Landscape;
  

  ElectronEnergyScale esc(electronEnergyScaleString);
  //esc.setCalibrationSet(ElectronEnergyScale::UNCORRECTED);
  esc.print();
  if (!esc.isInitialized()) {
    std::cout << "failed to initialize\n";
    return;
  }

  //TString path="../../root_files/selected_events/DY_m10+pr+a05+o03+pr_4680pb/ntuples/";
  //path="/media/spektras/DYee2011/2ndProd-ntuple/";
  //TString eemFileExp= path; eemFileExp.Append("data_bpVer1_EtaEtaM.root");
  //TString eemFileMC=path; eemFileMC.Append("zee_bpVer1_EtaEtaM.root");
  TString eemFileExp, eemFileMC;
  setEEMFileLocation(eemFileExp,eemFileMC);

  double xmin=60, xmax=120;
  std::vector<std::vector<EtaEtaMassData_t>*> mcDataAll,expDataAll;
  if (!esc.loadEEMFile( eemFileMC, mcDataAll,xmin,xmax) ||
      !esc.loadEEMFile(eemFileExp,expDataAll,xmin,xmax)) {
    std::cout << "failed to load data\n";
    return;
  }

  std::cout << "Data loaded" << std::endl;

  Distributions_t hSetCorr, hSetUnCorr;
  if (!hSetCorr.Create("ESCorr","hCorrected_", 60,60.,120.) ||
      !hSetUnCorr.Create("ESUnCorr","hUnCorrected_",60,60.,120.)) {
    std::cout << "failed to create distributions" << std::endl;
    return;
  }

  std::cout << "Distributions created" << std::endl;

  //
  // _EESF_None means uncorrected, otherwise it means that the appropriate
  // correction has to be applied
  //
  if (!hSetCorr.Fill(expDataAll,mcDataAll,esc,_EESF_Data) ||
      !hSetUnCorr.Fill(expDataAll,mcDataAll,esc,_EESF_None)) {
    std::cout << "failed to fill distributions\n";
    return;
  }

  std::cout << "Distributions filled" << std::endl;

  hSetCorr.NormalizeMCtoData();
  hSetUnCorr.NormalizeMCtoData();

  std::cout << "Distributions normalized" << std::endl;

  std::vector<CPlot*> cpQCorr, cpQUnCorr;
  std::vector<double> ymax;
  ymax.reserve(4);
  ymax.push_back(1e5);
  ymax.push_back(5e4);
  ymax.push_back(1e4);
  ymax.push_back(2e5);
  //for (unsigned int i=0; i<ymax.size(); i++) ymax[i]=1e8;
  
  TString xlabel="m(e^{+}e^{-}) [GeV]";
  TString ylabel="Events/1.0 GeV";
  if (!hSetCorr.CreatePlots(cpQCorr, "corrected",xlabel,ylabel,&ymax) ||
      !hSetUnCorr.CreatePlots(cpQUnCorr, "uncorrected",xlabel,ylabel,&ymax)) {
    std::cout << "failed to prepare CPlots\n";
    return;
  }

  std::cout << "Distributions plotted" << std::endl;

  // shift legend
  for (unsigned int i=0; i<4; ++i) {
    cpQCorr[i]->TransLegend(0.1,-0.05);
    cpQUnCorr[i]->TransLegend(0.1,-0.05);
  }

  // add luminosity text
  TString lumiText="#int#font[12]{L}dt = 4.7 fb^{-1}";
  for (unsigned int i=0; i<4; ++i) {
    cpQCorr[i]->AddTextBox(lumiText, 0.2,0.7,0.43,0.9, 0);
    cpQUnCorr[i]->AddTextBox(lumiText, 0.2,0.7,0.43,0.9, 0);
  }

  // add correction text
  TString corrText="corrected";
  TString uncorrText="uncorrected";
  for (unsigned int i=0; i<4; ++i) {
    cpQCorr[i]->AddTextBox(corrText, 0.4,0.92,0.73,0.98, 0);
    cpQUnCorr[i]->AddTextBox(uncorrText, 0.4,0.92,0.73,0.98, 0);
  }


  // Display the distributions
  int cwidth =900 - 300*int(pageOrientation);
  int cheight=600 + 300*int(pageOrientation);
  TCanvas *c = (1) ? MakeCanvas("c","c",cwidth,cheight) : NULL;
  bool save=0;
  TString format="png";
  if (c) {
    if (pageOrientation==_Portrait) c->Divide(2,3); else c->Divide(3,2);
    if (1) for (int i=0; i<3; ++i) {
      int idx=(pageOrientation==_Portrait) ? (2*i+1) : (i+1);
      cpQUnCorr[i]->Draw(c,save,format,idx);;
    }
    if (1) for (int i=0; i<3; ++i) {
      int idx=(pageOrientation==_Portrait) ? (2*i+2) : (i+4);
      cpQCorr[i]->Draw(c,save,format,idx);
    }

    std::cout << "Distributions showed" << std::endl;
    c->Update();
  }


  TCanvas *cAll = (0) ? MakeCanvas("cAll","cAll",1200,600) : NULL;
  if (cAll) {
    cAll->Divide(2);
    cpQUnCorr[3]->Draw(cAll,save,format,1);
    cpQCorr[3]->Draw(cAll,save,format,2);
    cAll->Update();
  }
}
