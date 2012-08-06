#ifndef myFitModels4_H
#define myFitModels4_H

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TCanvas.h>                // class for drawing
#include <TH1F.h>                   // 1D histograms
#include <Rtypes.h>   // colors
#include <iostream>                 // standard I/O
#include <vector>                   // STL vector class
#include <string>                   // C++ string class

#include "RooAbsReal.h"
//#include "RooAbsRealLValue.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooGaussModel.h"
#include "RooGExpModel.h"
//#include "RooChebychev.h"
#include "RooAddPdf.h"
#include <RooMoment.h>

//#include "RooPolyVar.h"
//#include "RooFitResult.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooGlobalFunc.h"

#include "RooExtendPdf.h"
#include "RooBreitWigner.h"
#include "RooVoigtian.h"
#include "RooBifurGauss.h"
#include "RooCBShape.h"
#include "RooFFTConvPdf.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooAddModel.h"
#include "RooMinuit.h"
#include "RooCategory.h"
#include "RooFormulaVar.h"
#include "RooSimultaneous.h"
#include "RooGenericPdf.h"

#include "TAxis.h"
#include "RooPlot.h"
#include "TLatex.h"

#ifndef __myLib__
#include "ProcessArgs.h"
#include "TPlotsAcc.hh"
#endif
#endif


//=== FUNCTION DECLARATIONS ======================================================================================


//=== MAIN MACRO =================================================================================================

extern int cPlotTheFit;
extern int cDumpTheFitContents;
extern int cExtraCanvasIdx;
extern int saveCFormat;
extern std::string saveFigFileTag;

struct DataInfos_t {
  std::vector<TString> *labelv;   // legend label
  std::vector<Int_t>   *colorv;   // color in plots
  std::vector<Int_t>   *linev;    // line style
  std::string FName;

  DataInfos_t() : labelv(0), colorv(0), linev(0), FName() {}
  DataInfos_t(std::vector<TString> *labs, std::vector<Int_t> *cols, std::vector<Int_t> *line_styles) : labelv(labs), colorv(cols), linev(line_styles), FName() {}
  DataInfos_t(DataInfos_t &info) : labelv(info.labelv), colorv(info.colorv), linev(info.linev), FName(info.FName) {}
  const std::vector<TString>* labels() const { return labelv; }
  const std::vector<Int_t>* colors() const { return colorv; }
  const std::vector<Int_t>* linestyles() const { return linev; }
  void Name(std::string name) { FName=name; }
  const std::string Name() const { return FName; }
  const char *getName() const { return FName.c_str(); }
  int NameOk() const { return (FName.size())?1:0; }
  TString label(UInt_t i) const { return (*labelv)[i]; }
  Int_t color(UInt_t i) const { return (*colorv)[i]; }
  Int_t linestyle(UInt_t i) const { return (*linev)[i]; }

  void Assign(std::vector<TString> *labs, std::vector<Int_t> *cols, std::vector<Int_t> *line_styles) {
    labelv=labs; colorv=cols; linev=line_styles;
  }

  void Assign(DataInfos_t &info) {
    labelv=info.labelv; colorv=info.colorv; linev=info.linev;
  }
};

//const int cFit1Gauss=1;
//const int cFit2Gauss=2;
//const int cFit2GaussCheby=3;
//const int cFitTestGaussRF102=4;
//const int cFitTestGaussRF102Demo=41; /// demo of a successful fit and unsuccessful fit
//const int cFitTestGaussRF102M=6;
//const int cFit1GaussTest=5;
//const int cFitTestExplicit1=71;
//const int cFitTestExplicit2=72;
//const int cFitExplicit2=7;
//const int cFitTestExplicit3=81;
//const int cFitExplicit3=8;
//const int cFitBreitWigner=9;
//const int cFitBreitWignerPF=900;  // pass-fail
//const int cFitBreitWignerPFStem=910;  // pass-fail stem model
//const int cFitBreitWignerPFStemT=911;  // pass-fail stem model
//const int cFitBreitWignerPFStemT2=912;  // pass-fail stem model
//const int cFitBreitWignerEff=90;
//const int cFitBreitWignerM1=91;
//const int cFitBreitWignerM1Eff=92;
const int cFitBreitWignerIlyaSpecial=100;
const int cFitBreitWignerIlyaLiteral=101;
//const int cFitBreitWignerIlyaSpecialT=102;
//const int cFitBreitWignerIlyaLiteralT=103;
const int cFitPFStemMCGauss=800;

#ifndef __myLib
int FitModelAllowed(int fit_model, int label_size, int fname_size, int pass_data_size, int fail_data_size);
#endif

//int Fit1GaussModel(const RooDataSet *gr, RooRealVar &x, RooRealVar &y, const DataInfos_t &info);
//int Fit2GaussModel(const RooDataSet *gr, RooRealVar &x, RooRealVar &y, const DataInfos_t &info);
//int Fit2GaussChebyModel(const RooDataSet *gr, RooRealVar &x, RooRealVar &y, const DataInfos_t &info);
//int FitTestGaussRF102(const RooDataSet *gr, RooRealVar &x, RooRealVar &y, const DataInfos_t &info);
//int FitTestGaussRF102Demo(const RooDataSet *gr, RooRealVar &x, RooRealVar &y, const DataInfos_t &info);
//int FitTestGaussRF102M(const RooDataSet *gr, RooRealVar &x, RooRealVar &y, const DataInfos_t &info); /// uses toy model 2 and calls Fit2GaussChebyModel
//int Fit1GaussModelTest(const RooDataSet *gr, RooRealVar &x, RooRealVar &y, const DataInfos_t &info);

//int FitExplicit1(const RooDataSet *gr, RooRealVar &x, RooRealVar &y, const DataInfos_t &info); // try to construct an explicit fitting function
//int FitExplicit2(const RooDataSet *gr, RooRealVar &x, RooRealVar &y, const DataInfos_t &info); // try to construct an explicit fitting function

//int FitTestExplicit1(const RooDataSet *gr, RooRealVar &x, RooRealVar &y, const DataInfos_t &info); /// uses toy model 2 and calls FitExplicit1
//int FitTestExplicit2(const RooDataSet *gr, RooRealVar &x, RooRealVar &y, const DataInfos_t &info); /// uses toy model 2 and calls FitExplicit2

//int FitExplicit2H(const RooDataHist *gr, RooRealVar &x, RooRealVar &y, const DataInfos_t &info); // try to construct an explicit fitting function
//int FitExplicit3H(const RooDataHist *gr, RooRealVar &x, RooRealVar &y, const DataInfos_t &info); // try to construct an explicit fitting function
//int FitTestExplicit2H(const RooDataHist *gr, RooRealVar &x, RooRealVar &y, const DataInfos_t &info); /// uses toy model 2 and calls FitExplicit2H
//int FitTestExplicit3H(const RooDataHist *gr, RooRealVar &x, RooRealVar &y, const DataInfos_t &info); /// uses toy model 2 and calls FitExplicit3H


//template<class rooDataClass_t>
//int FitBreitWigner(const rooDataClass_t *gr, RooRealVar &x, RooRealVar &y, const DataInfos_t &info);
//template<class rooDataClass_t>
//int FitBreitWignerEff(const rooDataClass_t *gr, RooRealVar &x, RooRealVar &y, const DataInfos_t &info);

// multi fit
//template<class rooDataClass_t>
//int FitBreitWignerV(const std::vector<rooDataClass_t*> &ds, RooRealVar &x, RooRealVar &y, const DataInfos_t &info);
//template<class rooDataClass_t>
//int FitBreitWignerVPF(const std::vector<rooDataClass_t*> &ds, RooRealVar &x, RooRealVar &y, const DataInfos_t &info);
//template<class rooDataClass_t>
//int FitBreitWignerVPFStem(const std::vector<rooDataClass_t*> &ds, RooRealVar &x, RooRealVar &y, const DataInfos_t &info);
//template<class rooDataClass_t>
//int FitBreitWignerVM1(const std::vector<rooDataClass_t*> &ds, RooRealVar &x, RooRealVar &y, const DataInfos_t &info);
//template<class rooDataClass_t>
//int FitBreitWignerVM1Eff(const std::vector<rooDataClass_t*> &ds, RooRealVar &x, RooRealVar &y, const DataInfos_t &info);


int FitBreitWignerIlyaSpecial(const std::vector<RooDataSet*> &dataV, RooRealVar &x, RooRealVar &y, const DataInfos_t &info);
int FitBreitWignerIlyaLiteral(const std::vector<RooDataSet*> &dataV, RooRealVar &x, RooRealVar &y, const DataInfos_t &info);

//template<class rooDataClass_t>
//int FitBreitWignerVPFStemT(const std::vector<rooDataClass_t*> &ds, RooRealVar &x, RooRealVar &y, const DataInfos_t &info);
//template<class rooDataClassBinned_t, class rooDataClassUnbinned_t>
//int FitBreitWignerVPFStemT2(const std::vector<rooDataClassBinned_t*> &ds_binned, const std::vector<rooDataClassUnbinned_t*> &ds, RooRealVar &x, RooRealVar &y, const DataInfos_t &info);
//template<class rooDataClass_t>
//int FitBreitWignerIlyaSpecialT(const std::vector<rooDataClass_t*> &dataV, RooRealVar &x, RooRealVar &y, const DataInfos_t &info);
//template<class rooDataClass_t>
//int FitBreitWignerIlyaLiteralT(const std::vector<rooDataClass_t*> &dataV, RooRealVar &x, RooRealVar &y, const DataInfos_t &info);

int FitPFStemMCGauss(const std::vector<RooDataSet*> &ds, RooRealVar &x, RooRealVar &y, const DataInfos_t &info);


// --------------------------------------------------------------------
// --------------------------------------------------------------------

void PrepareNumCharVec(unsigned int count, std::vector<char> &idxV);
void PrepareCharVec(unsigned int count, std::vector<char> &idxV);
void PrepareNameEndings(unsigned int count, std::vector<char> &idxV, std::vector<TString> &nameEndings);

void PrepareNumCharVec(unsigned int count, std::vector<TString> &idxV);
void PrepareCharVec(unsigned int count, std::vector<TString> &idxV);
void PrepareNameEndings(unsigned int count, std::vector<TString> &idxV, std::vector<TString> &nameEndings);

template<class Class_t>
void PrintValues(const char *msg, const std::vector<Class_t*> &vec);
template<class Class_t>
void PrintValues(const char *msg, const std::vector<std::vector<Class_t*>*> &vec);
template<class Class_t>
void PrintPrint(const char *msg, const std::vector<Class_t*> &vec);
template<class Class_t>
void PrintPrint(const char *msg, const std::vector<std::vector<Class_t*>*> &vec);


int makeHTML(const TString outDir, const TString html_fname, const TString description, const std::vector<std::string> &add_html_lines, const std::vector<TString> &nameEndings, const std::vector<TString> &nameEndingsForFile, RooRealVar &x, const std::vector<RooDataSet*> &combData, const std::vector<RooAddPdf*> &simPdf, const std::vector<std::vector<int>*>* combine=NULL, const std::vector<TString> *combination_names=NULL);
int makeHTML(const TString outDir, const TString html_fname, const std::vector<TString> &description, const std::vector<std::string> &add_html_lines, const std::vector<TString> &nameEndings, const std::vector<TString> &nameEndingsForFile, std::vector<RooRealVar*> &x, const std::vector<const std::vector<RooDataSet*>*> &combData, const std::vector<const std::vector<RooAddPdf*>*> &simPdf, const std::vector<std::vector<int>*>* combine=NULL, const std::vector<TString> *combination_names=NULL);


// --------------------------------------------------------------------
// --------------------------------------------------------------------

RooDataSet* ConvertToDataSet(const std::vector<double> *dt, RooRealVar &x, RooRealVar &y);
std::vector<RooDataSet*>* ConvertToDataSet(const std::vector<      std::vector<double>*> &dt, RooRealVar &x, RooRealVar &y, int reverse_order=0);
std::vector<RooDataSet*>* ConvertToDataSet(const std::vector<const std::vector<double>*> &dt, const std::vector<std::string> &names, RooRealVar &x, RooRealVar &y);

#ifndef __myLib__
std::vector<RooDataSet*>* ConvertToDataSet(const std::vector<mithep::AccEffData_t*> &dataP, const std::vector<mithep::AccEffData_t*> &dataF, RooRealVar &x, RooRealVar &y);
std::vector<std::vector<double>*>* ConvertMassToVectorDbl(const std::vector<mithep::AccEffData_t*> &dataP, const std::vector<mithep::AccEffData_t*> &dataF);
#endif

template<class HistoClass_t>
HistoClass_t* ConvertToHisto(const RooDataSet &dt, const TString &name, const HistoClass_t *templateHisto=NULL, const char *variable_name="x");
template<class HistoClass_t>
std::vector<HistoClass_t*>* ConvertToHisto(const std::vector<RooDataSet*> &dt, const std::vector<TString> &name, const HistoClass_t *templateHisto=NULL, const char *variable_name="x");

template<class TRooAbsRealDescendant_t, class HistoClass_t>
  inline HistoClass_t* CreateHisto(const TRooAbsRealDescendant_t& rooVar, const char *name, const RooRealVar *xAxisArg, const RooCmdArg &args1, const HistoClass_t *dummyH=NULL) {
  HistoClass_t *h;
  h=(HistoClass_t*)rooVar.createHistogram(name,*(RooAbsRealLValue*)xAxisArg,args1);
  if (dummyH) std::cout << "dummyH is not null! (ignored)\n";
  return h;
}

template<class TRooAbsRealDescendant_t, class HistoClass_t>
  inline HistoClass_t* CreateHisto(const TRooAbsRealDescendant_t& rooVar, const char *name, const RooRealVar *xAxisArg, int bin_count, const double *bins, const HistoClass_t *dummyH=NULL);

template<class TRooAbsRealDescendant_t, class HistoClass_t>
  inline HistoClass_t* CreateHisto(const TRooAbsRealDescendant_t& rooVar, const char *name, const RooRealVar *xAxisArg, int bin_count, const double *bins, const HistoClass_t *dummyH=NULL) {
  assert(bins);
  HistoClass_t *h=(HistoClass_t*)dummyH->Clone(name);
  if (1) { // both versions give the same central values
    const double dx=0.2;
    std::cout << "CreateHisto(TRooAbsRealDescendant): granularity dx=" << dx << "\n";
    double xmin=bins[0];
    double xmax=bins[bin_count];
    int locBinCount=int((xmax-xmin)/dx+0.001);
    RooCmdArg arg=RooFit::Binning(locBinCount,xmin,xmax);
    //std::cout << "localBinCount=" << locBinCount << ", for x range " << xmin << ".." << xmax << "\n";
    HistoClass_t *htmp=CreateHisto(rooVar,"tmpH",xAxisArg,arg,dummyH);
    //std::cout << "the histogram got created" << std::endl;
    int iloc=0;
    double xc;
    xc=iloc*dx+0.5*dx+bins[0];
    for (int i=0; i<bin_count; ++i) {
      do {
	if ((xc>=bins[i]) && (xc<bins[i+1])) {
	  h->Fill(xc,htmp->GetBinContent(iloc+1));
	  //std::cout << "i=" << i << ", iloc=" << iloc << ", xc=" << xc << ", bins[i]=" << bins[i] << "\n";
 	}
	iloc++;
	xc=iloc*dx+0.5*dx+bins[0];
     } while (xc<bins[i+1]);
    }
    delete htmp;    
  }
  else {
    for (int i=0; i<bin_count; ++i) {
      RooCmdArg arg=RooFit::Binning(1,bins[i],bins[i+1]);
      HistoClass_t *htmp=(HistoClass_t*)rooVar.createHistogram(name,*(RooAbsRealLValue*)xAxisArg,arg);
      if (0) {
	h->SetBinContent(i+1,htmp->GetBinContent(1));
	h->SetBinError(i+1,htmp->GetBinError(1));
      }
      else {
	double center=htmp->GetBinCenter(1);
	h->Fill(center, htmp->GetBinContent(1));
      }
      delete htmp;
    }
  }
  return h;
}

// --------------------------------------------------------------------
// --------------------------------------------------------------------

int CreatePassFailStemFormula(unsigned int count, std::vector<RooRealVar*> &eff, RooRealVar &NEvents, std::vector<RooFormulaVar*> &formula);

int SetInitialPassFailEffValues(unsigned int count, const std::vector<RooDataSet*> &data, std::vector<RooRealVar*> &eff, int &total_event_count);

RooDataSet* CreateCombinedData(TString name,TString descr, const std::vector<RooDataSet*> &data, RooRealVar &x, RooCategory &sample, const std::vector<TString> &categories);
RooSimultaneous* CreateCombinedModel(TString name, TString descr, std::vector<RooAddPdf*> &model, RooCategory &sample, const std::vector<TString> &categories);

// --------------------------------------------------------------------
// --------------------------------------------------------------------

int CreateRooCategories(unsigned int count, TString name_base, const std::vector<TString> &name_endings, RooCategory &sample, std::vector<TString> &categories);
int CreateRooCategories(unsigned int count, TString name, const std::vector<char> &idxV, RooCategory &sample, std::vector<TString> &categories);

template<class THistoClass_t>
int CreateHistos(unsigned int count, TString name, std::vector<THistoClass_t*> &histos, int bin_count, double x_min, double x_max);

template<class THistoClass_t>
int CreateHistos(unsigned int count, TString name, std::vector<THistoClass_t*> &histos, int bin_count, const double *bin_limits);

template<class Class_t>
int CreateRooArgList(RooArgList &list, std::vector<Class_t*> &objs, int append=0);
template<class Class_t>
inline int CreateRooArgList(std::vector<Class_t*> &objs, RooArgList &list, int append=0);
template<class Class_t>
inline int CreateRooArgList(std::vector<Class_t*> &objs, std::vector<RooArgList*> &list, int append=0);

int CreateRooRealVar(unsigned int count, TString name, std::vector<RooRealVar*> &vars, double start_value, double min_value=-9999, double max_value=-9999);
int CreateRooRealVar(unsigned int count, const TString &name_base, const std::vector<TString> &name_endings, std::vector<RooRealVar*> &vars, double start_value, double min_value=-9999, double max_value=-9999);

int CreateRooExponential(unsigned int count, TString name, RooRealVar &x, std::vector<RooRealVar*> &tauV, std::vector<RooExponential*> &expos);

int CreateRooHistPdf(const RooRealVar &x, const std::vector<RooDataSet*> &data, TString name1, std::vector<RooDataHist*> &histos, TString name2, std::vector<RooHistPdf*> &hpdfs, int from_idx=-1, int to_idx=-1);

template<class THistoClass_t>
int CreateRooHistPdf(const RooRealVar &x, const std::vector<THistoClass_t*> &data, TString name1, std::vector<RooDataHist*> &histos, TString name2, std::vector<RooHistPdf*> &hpdfs, int from_idx=-1, int to_idx=-1);

template<class XClass_t, class Gauss_t, class Mean_t, class Width_t>
int CreateRooGaussian(unsigned int count, TString name, std::vector<XClass_t*> &x, std::vector<Mean_t*> &mean, std::vector<Width_t*> &width, std::vector<Gauss_t*> &gauss);
template<class XClass_t, class Gauss_t, class Mean_t, class Width_t>
int CreateRooGaussian(unsigned int count, TString name, std::vector<XClass_t*> &x, Mean_t &mean, std::vector<Width_t*> &width, std::vector<Gauss_t*> &gauss);

template<class Gauss_t, class Mean_t, class Width_t>
int CreateRooGaussian(unsigned int count, TString name, RooRealVar &x, std::vector<Mean_t*> &mean, std::vector<Width_t*> &width, std::vector<Gauss_t*> &gauss);
template<class Gauss_t, class Mean_t, class Width_t>
  int CreateRooGaussian(unsigned int count, TString name, RooRealVar &x, Mean_t &mean, std::vector<Width_t*> &width, std::vector<Gauss_t*> &gauss);
template<class Gauss_t, class Mean_t, class Width_t>
int CreateRooGaussian(unsigned int count, TString name, RooRealVar &x, std::vector<Mean_t*> &mean, Width_t &width, std::vector<Gauss_t*> &gauss);
template<class Gauss_t, class Mean_t, class Width_t>
int CreateRooGaussian(unsigned int count, TString name, RooRealVar &x, Mean_t &mean, Width_t &width, std::vector<Gauss_t*> &gauss);

template<class GaussExp_t, class Mean_t, class Width_t, class Tau_t>
  int CreateRooGaussExp(unsigned int count, TString name, RooRealVar &x, Mean_t &mean, std::vector<Width_t*> &width, std::vector<Tau_t*> &tau, std::vector<GaussExp_t*> &gauss);
template<class GaussExp_t, class Mean_t, class Width_t, class Tau_t>
  int CreateRooGaussExp(unsigned int count, TString name, RooRealVar &x, std::vector<Mean_t*> &mean, std::vector<Width_t*> &width, std::vector<Tau_t*> &tau, std::vector<GaussExp_t*> &gauss);

template<class Mean_t, class Width1_t, class Width2_t>
  int CreateRooBifurGauss(unsigned int count, TString name, RooRealVar &x, Mean_t &mean, std::vector<Width1_t*> &width1, std::vector<Width2_t> &width2, std::vector<RooBifurGauss*> &bifGauss);

template<class Mean_t, class Width_t, class Alpha_t, class CBParN_t>
int CreateRooCBShape(unsigned int count, TString name, RooRealVar &x, Mean_t &mean, std::vector<Width_t*> &width, std::vector<Alpha_t*> &alpha, std::vector<CBParN_t*> &parN, std::vector<RooCBShape*> &cbshape);
template<class Mean_t, class Width_t, class Alpha_t, class CBParN_t>
int CreateRooCBShape(unsigned int count, TString name, RooRealVar &x, std::vector<Mean_t*> &mean, std::vector<Width_t*> &width, std::vector<Alpha_t*> &alpha, std::vector<CBParN_t*> &parN, std::vector<RooCBShape*> &cbshape);

template<class Mean_t, class Width_t>
  int CreateRooBreitWigner(unsigned int count, TString name, RooRealVar &x, Mean_t &mean, std::vector<Width_t*> &width, std::vector<RooBreitWigner*> &bwV);

template<class Mean_t, class Width_t, class Sigma_t>
  int CreateRooVoigtian(unsigned int count, TString name, RooRealVar &x, Mean_t &mean, std::vector<Width_t*> &width, std::vector<Sigma_t*> &sigma, std::vector<RooVoigtian*> &bwV);


template<class FormClass_t, class CountClass_t>
  int CreateRooExtendPdf(unsigned int count, TString name, std::vector<FormClass_t*> &form, std::vector<CountClass_t*> &counts, std::vector<RooExtendPdf*> &extended);
template<class FormClass_t, class CountClass_t>
  int CreateRooExtendPdf(unsigned int count, TString name, std::vector<FormClass_t*> &form, CountClass_t &countVar, std::vector<RooExtendPdf*> &extended);
template<class FormClass_t, class CountClass_t>
  int CreateRooExtendPdf(unsigned int count, TString namebase, const std::vector<TString> &nameEndings, std::vector<FormClass_t*> &form, CountClass_t &countVar, std::vector<RooExtendPdf*> &extended);

template<class TRooAddPdf_t>
int CreateRooAddPdf(unsigned int count, TString name, std::vector<RooArgList*> &lists, std::vector<TRooAddPdf_t*> &model);
template<class TRooAddPdf_t>
int CreateRooAddPdf(unsigned int count, TString name, const std::vector<TString> &endings, std::vector<RooArgList*> &lists, std::vector<TRooAddPdf_t*> &model);

int CreateRooAddPdf(unsigned int count, TString name, std::vector<RooExtendPdf*> &var1, std::vector<RooExtendPdf*> &var2, std::vector<RooAddPdf*> &model);
//int CreateRooAddPdf(unsigned unsigned int count, TString name_base, const std::vector<TString> &endings, std::vector<RooExtendPdf*> &var1, std::vector<RooExtendPdf*> &var2, std::vector<RooAddPdf*> &model);

template<class TRooAddModel_t, class TCoefList_t>
int CreateRooAddModel(unsigned int count, TString name, std::vector<RooArgList*> &lists, std::vector<TCoefList_t*> &coefs, std::vector<TRooAddModel_t*> &model);
template<class TRooAddModel_t, class TCoefList_t>
int CreateRooAddModel(unsigned int count, TString name, const std::vector<TString> &endings, std::vector<RooArgList*> &lists, std::vector<TCoefList_t*> &coefs, std::vector<TRooAddModel_t*> &model);

template<class Class_t>
int CreateRooFormulaVar(unsigned int count, TString name, TString formulaString, RooArgList &defArgs, std::vector<Class_t*> &vars, std::vector<RooFormulaVar*> &formula);
template<class Class1_t, class Class2_t>
int CreateRooFormulaVar(unsigned int count, TString name, TString formulaString, RooArgList &defArgs, std::vector<Class1_t*> &var1, std::vector<Class2_t*> &var2, std::vector<RooFormulaVar*> &formula);

template<class Resol_t, class Shape_t, class Conv_t>  // e.g. RooGaussModel, RooHistoPdf & RooNumConvPdf; RooGaussian, RooFormula & RooFFTConvPdf;   note: FFTConv may need a different ordering than NumConv
int CreateRooConvPdf(unsigned int count, TString name, RooRealVar &x, std::vector<Resol_t*> &resol, std::vector<Shape_t*> &shape, std::vector<Conv_t*> &convolution);



// --------------------------------------------------------------------
// --------------------------------------------------------------------

template<class HistoClass_t>
HistoClass_t* ConvertToHisto(const RooDataSet &dt, const TString &name, const HistoClass_t *templateHisto=NULL, const char *variable_name="x") {
  assert(name);
  HistoClass_t *h=(templateHisto) ? (HistoClass_t*)templateHisto->Clone(name) : new HistoClass_t(name.Data(),name.Data(),101,1,100);
  assert(h);
  for (int ii=0; ii<dt.numEntries(); ii++) {
    const RooArgSet *arg=dt.get(ii);
    h->Fill(arg->getRealValue(variable_name,0,kTRUE));
  }
  return h;
}


template<class HistoClass_t>
std::vector<HistoClass_t*>* ConvertToHisto(const std::vector<RooDataSet*> &dt, const std::vector<TString> &names, const HistoClass_t* templateHisto=NULL, const char *variable_name="x") {
  std::vector<HistoClass_t*>* hh=new std::vector<HistoClass_t*>();
  for (unsigned int i=0; i<dt.size(); ++i) {
    hh->push_back(ConvertToHisto(*dt[i],names[i],templateHisto,variable_name));
  }
  return hh;
}

// --------------------------------------------------------------------
// --------------------------------------------------------------------

template<class THistoClass_t>
int CreateHistos(unsigned int count, TString name, std::vector<THistoClass_t*> &histos, int bin_count, double x_min, double x_max) {
  histos.reserve(count);
  size_t pos=name.Length()-1;
  for (unsigned int i=0; i<count; ++i) {
    name[pos]=char(int('a')+i);
    histos.push_back(new THistoClass_t(name,name,bin_count,x_min,x_max));
    histos.back()->SetDirectory(0);
  }
  return 1;
}

template<class THistoClass_t>
int CreateHistos(unsigned int count, TString name, std::vector<THistoClass_t*> &histos, int bin_count, const double *bin_limits) {
  histos.reserve(count);
  size_t pos=name.Length()-1;
  for (unsigned int i=0; i<count; ++i) {
    name[pos]=char(int('a')+i);
    histos.push_back(new THistoClass_t(name,name,bin_count,bin_limits));
    histos.back()->SetDirectory(0);
  }
  return 1;
}


template<class Class_t>
int CreateRooArgList(RooArgList &list, std::vector<Class_t*> &objs, int append) {
  if (!append) list.removeAll();
  for (unsigned int i=0; i<objs.size(); ++i) {
    list.add(*objs[i]);
  }
  return 1;
}


template<class Class_t>
int CreateRooArgList(std::vector<Class_t*> &objs, std::vector<RooArgList*> &list, int append) {
  if (!append) {
    for (unsigned int i=0; i<list.size(); i++) delete list[i];
    list.clear();
  }  
  else {
    if (objs.size()!=list.size()) {
      std::cout << "CreateRooArgList: append=1, sizes are different, objs.size=" << objs.size() << ", exsisting list.size=" << list.size() << "\n";
      return 0;
    }
  }
  list.reserve(objs.size());
  TString name="args_X"; size_t pos=name.Length()-1;
  for (unsigned int i=0; i<objs.size(); ++i) {
    name[pos]=char(int('a')+i);
    if (!append) {
      list.push_back(new RooArgList(name));
    }
    list[i]->add(*objs[i]);
  }
  return 1;
}


template<class THistoClass_t>
int CreateRooHistPdf(const RooRealVar &x, const std::vector<THistoClass_t*> &data, TString name1, std::vector<RooDataHist*> &histos, TString name2, std::vector<RooHistPdf*> &hpdfs, int from_idx, int to_idx) {
  int N=data.size();
  if (from_idx<0) from_idx=0;
  if (to_idx<0) to_idx=N;
  if ((from_idx>N) || (to_idx>N) || (to_idx-from_idx<=0)) {
    std::cout << "CreateRooHistPdf(HistoClass) for names=" << name1 << "and " << name2 << ": from_idx=" << from_idx << ", to_idx=" << to_idx << "\n";
    return 0;
  }
  N=to_idx-from_idx;
  histos.reserve(N); hpdfs.reserve(N);
  size_t pos1=name1.Length()-1;
  size_t pos2=name2.Length()-1;
  if (1) {
    for (int i=from_idx; i<to_idx; i++) {
      const char c=char(int('a')-from_idx+i);
      name1[pos1]=c; name2[pos2]=c;
      histos.push_back(new RooDataHist(name1,name1,RooArgSet(x),data[i]));
      hpdfs.push_back(new RooHistPdf(name2,name2,RooArgSet(x),*histos.back()));
    }
  }
  return 1;
}


template<class XClass_t, class Gauss_t, class Mean_t, class Width_t>
int CreateRooGaussian(unsigned int count, TString name, std::vector<XClass_t*> &x, std::vector<Mean_t*> &mean, std::vector<Width_t*> &width, std::vector<Gauss_t*> &gauss) {
  if ((count!=x.size()) || (count!=mean.size()) || (count!=width.size()))  {
    std::cout << "CreateRooGaussian<T> for name=" << name << ": count=" << count << ", x.size=" << x.size() << ", mean.size=" << mean.size() << ", width.size=" << width.size() << "\n";
    return 0;
  }
  gauss.reserve(count);
  size_t pos=name.Length()-1;
  for (unsigned int i=0; i<count; ++i) {
    name[pos]=char(int('a')+i);
    gauss.push_back(new Gauss_t(name,name, *x[i],*mean[i],*width[i]));
  }
  return 1;
}

template<class XClass_t, class Gauss_t, class Mean_t, class Width_t>
int CreateRooGaussian(unsigned int count, TString name, std::vector<XClass_t*> &x, Mean_t &mean, std::vector<Width_t*> &width, std::vector<Gauss_t*> &gauss) {
  if ((count!=x.size()) || (count!=width.size()))  {
    std::cout << "CreateRooGaussian<T> for name=" << name << ": count=" << count << ", x.size=" << x.size() << ", width.size=" << width.size() << "\n";
    return 0;
  }
  gauss.reserve(count);
  size_t pos=name.Length()-1;
  for (unsigned int i=0; i<count; ++i) {
    name[pos]=char(int('a')+i);
    gauss.push_back(new Gauss_t(name,name, *x[i],mean,*width[i]));
  }
  return 1;
}


template<class Gauss_t, class Mean_t, class Width_t>
int CreateRooGaussian(unsigned int count, TString name, RooRealVar &x, std::vector<Mean_t*> &mean, std::vector<Width_t*> &width, std::vector<Gauss_t*> &gauss) {
  if ((count!=mean.size()) || (count!=width.size()))  {
    std::cout << "CreateRooGaussian<T> for name=" << name << ": count=" << count << ", mean.size=" << mean.size() << ", width.size=" << width.size() << "\n";
    return 0;
  }
  gauss.reserve(count);
  size_t pos=name.Length()-1;
  for (unsigned int i=0; i<count; ++i) {
    name[pos]=char(int('a')+i);
    gauss.push_back(new Gauss_t(name,name, x,*mean[i],*width[i]));
  }
  return 1;
}

template<class Gauss_t, class Mean_t, class Width_t>
int CreateRooGaussian(unsigned int count, TString name, RooRealVar &x, std::vector<Mean_t*> &mean, Width_t &width, std::vector<Gauss_t*> &gauss) {
  if (count!=mean.size())  {
    std::cout << "CreateRooGaussian<T> for name=" << name << ": count=" << count << ", mean.size=" << mean.size() << "\n";
    return 0;
  }
  gauss.reserve(count);
  size_t pos=name.Length()-1;
  for (unsigned int i=0; i<count; ++i) {
    name[pos]=char(int('a')+i);
    gauss.push_back(new Gauss_t(name,name, x,*mean[i],width));
  }
  return 1;
}

template<class Gauss_t, class Mean_t, class Width_t>
int CreateRooGaussian(unsigned int count, TString name, RooRealVar &x, Mean_t &mean, std::vector<Width_t*> &width, std::vector<Gauss_t*> &gauss) {
  if (count!=width.size())  {
    std::cout << "CreateRooGaussian<T> for name=" << name << ": count=" << count << ", width.size=" << width.size() << "\n";
    return 0;
  }
  gauss.reserve(count);
  size_t pos=name.Length()-1;
  for (unsigned int i=0; i<count; ++i) {
    name[pos]=char(int('a')+i);
    gauss.push_back(new Gauss_t(name,name, x,mean,*width[i]));
  }
  return 1;
}


template<class Gauss_t, class Mean_t, class Width_t>
int CreateRooGaussian(unsigned int count, TString name, RooRealVar &x, Mean_t &mean, Width_t &width, std::vector<Gauss_t*> &gauss) {
  gauss.reserve(count);
  size_t pos=name.Length()-1;
  for (unsigned int i=0; i<count; ++i) {
    name[pos]=char(int('a')+i);
    gauss.push_back(new Gauss_t(name,name, x,mean,width));
  }
  return 1;
}


template<class GaussExp_t, class Mean_t, class Width_t, class Tau_t>
int CreateRooGaussExp(unsigned int count, TString name, RooRealVar &x, Mean_t &mean, std::vector<Width_t*> &width, std::vector<Tau_t*> &tau, std::vector<GaussExp_t*> &gauss) {
  if ((count!=width.size()) || (count!=tau.size()))  {
    std::cout << "CreateRooGaussExp<T> for name=" << name << ": count=" << count << ", width.size=" << width.size() << ", tau.size=" << tau.size() << "\n";
    return 0;
  }
  gauss.reserve(count);
  size_t pos=name.Length()-1;
  for (unsigned int i=0; i<count; ++i) {
    name[pos]=char(int('a')+i);
    gauss.push_back(new GaussExp_t(name,name, x,mean,*width[i],*tau[i]));
  }
  return 1;
}

template<class GaussExp_t, class Mean_t, class Width_t, class Tau_t>
int CreateRooGaussExp(unsigned int count, TString name, RooRealVar &x, std::vector<Mean_t*> &mean, std::vector<Width_t*> &width, std::vector<Tau_t*> &tau, std::vector<GaussExp_t*> &gauss) {
  if ((count!=mean.size()) || (count!=width.size()) || (count!=tau.size()))  {
    std::cout << "CreateRooGaussian<T> for name=" << name << ": count=" << count << ", mean.size=" << mean.size() << ", width.size=" << width.size() << ", tau.size=" << tau.size() << "\n";
    return 0;
  }
  gauss.reserve(count);
  size_t pos=name.Length()-1;
  for (unsigned int i=0; i<count; ++i) {
    name[pos]=char(int('a')+i);
    gauss.push_back(new GaussExp_t(name,name, x,*mean[i],*width[i],*tau[i]));
  }
  return 1;
}

template<class Mean_t, class Width1_t, class Width2_t>
  int CreateRooBifurGauss(unsigned int count, TString name, RooRealVar &x, Mean_t &mean, std::vector<Width1_t*> &width1, std::vector<Width2_t> &width2, std::vector<RooBifurGauss*> &bifGauss) {
  if ((count!=width1.size()) || (count!=width2.size())) {
    std::cout << "CreateRooBifurGauss<T> for name=" << name << ": count=" << count << ", width1.size=" << width1.size() << ", width2.size=" << width2.size() << "\n";
    return 0;
  }
  bifGauss.reserve(count);
  size_t pos=name.Length()-1;
  for (unsigned int i=0; i<count; ++i) {
    name[pos]=char(int('a')+i);
    bifGauss.push_back(new RooBifurGauss(name,name, x,mean,*width1[i],*width2[i]));
  }
  return 1;
}

template<class Mean_t, class Width_t>
int CreateRooBreitWigner(unsigned int count, TString name, RooRealVar &x, Mean_t &mean, std::vector<Width_t*> &width, std::vector<RooBreitWigner*> &bwV) {
  if ((count!=width.size()))  {
    std::cout << "CreateRooBreightWigner<T> for name=" << name << ": count=" << count << ", width.size=" << width.size() << "\n";
    return 0;
  }
  bwV.reserve(count);
  size_t pos=name.Length()-1;
  for (unsigned int i=0; i<count; ++i) {
    name[pos]=char(int('a')+i);
    bwV.push_back(new RooBreitWigner(name,name, x,mean,*width[i]));
  }
  return 1;
}

template<class Mean_t, class Width_t, class Sigma_t>
  int CreateRooVoigtian(unsigned int count, TString name, RooRealVar &x, Mean_t &mean, std::vector<Width_t*> &width, std::vector<Sigma_t*> &sigma, std::vector<RooVoigtian*> &voigtianV) {
  if ((count!=width.size()) || (count!=sigma.size()))  {
    std::cout << "CreateRooBreightWigner<T> for name=" << name << ": count=" << count << ", width.size=" << width.size() << ", sigma.size=" << sigma.size() << "\n";
    return 0;
  }
  voigtianV.reserve(count);
  size_t pos=name.Length()-1;
  for (unsigned int i=0; i<count; ++i) {
    name[pos]=char(int('a')+i);
    voigtianV.push_back(new RooVoigtian(name,name, x,mean,*width[i],*sigma[i]));
  }
  return 1;
}


template<class Mean_t, class Width_t, class Alpha_t, class CBParN_t>
int CreateRooCBShape(unsigned int count, TString name, RooRealVar &x, std::vector<Mean_t*> &mean, std::vector<Width_t*> &width, std::vector<Alpha_t*> &alpha, std::vector<CBParN_t*> &parN, std::vector<RooCBShape*> &cbshape) {
  if ((count!=mean.size()) || (count!=width.size()) || (count!=alpha.size()) || (count!=parN.size()))  {
    std::cout << "CreateRooCBShape<T> for name=" << name << ": count=" << count << ", mean.size=" << mean.size() << ", width.size=" << width.size() << ", alpha.size=" << alpha.size() << ", parN.size=" << parN.size() << "\n";
    return 0;
  }
  cbshape.reserve(count);
  size_t pos=name.Length()-1;
  for (unsigned int i=0; i<count; ++i) {
    name[pos]=char(int('a')+i);
    cbshape.push_back(new RooCBShape(name,name, x,*mean[i],*width[i],*alpha[i],*parN[i]));
  }
  return 1;
}

template<class Mean_t, class Width_t, class Alpha_t, class CBParN_t>
int CreateRooCBShape(unsigned int count, TString name, RooRealVar &x, Mean_t &mean, std::vector<Width_t*> &width, std::vector<Alpha_t*> &alpha, std::vector<CBParN_t*> &parN, std::vector<RooCBShape*> &cbshape) {
  if ((count!=width.size()) || (count!=alpha.size()) || (count!=parN.size()))  {
    std::cout << "CreateRooCBShape<T> for name=" << name << ": count=" << count << ", width.size=" << width.size() << ", alpha.size=" << alpha.size() << ", parN.size=" << parN.size() << "\n";
    return 0;
  }
  cbshape.reserve(count);
  size_t pos=name.Length()-1;
  for (unsigned int i=0; i<count; ++i) {
    name[pos]=char(int('a')+i);
    cbshape.push_back(new RooCBShape(name,name, x,mean,*width[i],*alpha[i],*parN[i]));
  }
  return 1;
}




template<class FormClass_t, class CountClass_t>
 int CreateRooExtendPdf(unsigned int count, TString name, std::vector<FormClass_t*> &form, std::vector<CountClass_t*> &counts, std::vector<RooExtendPdf*> &extended) {
  if ((count==0) || (count!=form.size()) || (count!=counts.size())) {
    std::cout << "CreateRooExtendPdf for name=" << name << ": count=" << count << ", form.size=" << form.size() << ", counts.size=" << counts.size() << "\n";
    return 0;
  }
  extended.reserve(count);
  size_t pos=name.Length()-1;
  for (unsigned int i=0; i<count; ++i) {
    const int span=(int('z')-int('a')+1);
    int rem=i%span;
    char ch=char(int('a')+rem);
    if (i/span>0) {
      name[pos-1]=char('a'+(i/span-1));
    }
    name[pos]=ch;
    extended.push_back(new RooExtendPdf(name,name, *form[i], *counts[i]));
  }
  return 1;
}

template<class FormClass_t, class CountClass_t>
 int CreateRooExtendPdf(unsigned int count, TString name, std::vector<FormClass_t*> &form, CountClass_t &countVar, std::vector<RooExtendPdf*> &extended) {
  if ((count==0) || (count!=form.size())) {
    std::cout << "CreateRooExtendPdf for name=" << name << ": count=" << count << ", form.size=" << form.size() << "\n";
    return 0;
  }
  extended.reserve(count);
  size_t pos=name.Length()-1;
  for (unsigned int i=0; i<count; ++i) {
    name[pos]=char(int('a')+i);
    const char *cn=name.Data();
    extended.push_back(new RooExtendPdf(cn,cn, *form[i], countVar));
  }
  return 1;
}

template<class FormClass_t, class CountClass_t>
int CreateRooExtendPdf(unsigned int count, TString namebase, const std::vector<TString> &nameEndings, std::vector<FormClass_t*> &form, CountClass_t &countVar, std::vector<RooExtendPdf*> &extended) {
  if ((count==0) || (count!=form.size())) {
    std::cout << "CreateRooExtendPdf for namebase=" << namebase << ": count=" << count << ", form.size=" << form.size() << "\n";
    return 0;
  }
  if (count!=nameEndings.size()) {
    std::cout << "CreateRooExtendPdf for namebase=" << namebase << ": count=" << count << ", nameEndings.size=" << nameEndings.size() << "\n";
  }
  extended.reserve(count);
  TString name;
  for (unsigned int i=0; i<count; ++i) {
    name=namebase+nameEndings[i];
    const char *cn=name.Data();
    extended.push_back(new RooExtendPdf(cn,cn, *form[i], countVar));
  }
  return 1;
}


// ---------------------------------------------------------------------

template<class TRooAddPdf_t>
int CreateRooAddPdf(unsigned int count, TString name, std::vector<RooArgList*> &lists, std::vector<TRooAddPdf_t*> &model) {
  if ((count==0) || (count!=lists.size())) {
    std::cout << "CreateRooAddPdf<T> for name=" << name << ": count=" << count << ", list.size=" << lists.size() << "\n";
    return 0;
  }
  model.reserve(count);
  size_t pos=name.Length()-1;
  for (unsigned int i=0; i<count; ++i) {
    name[pos]=char(int('a')+i);
    model.push_back(new TRooAddPdf_t(name,name,*lists[i]));
  }
  return 1;
}
// ---------------------------------------------------------------------

template<class TRooAddPdf_t>
int CreateRooAddPdf(unsigned int count, TString name_base, const std::vector<TString> &endings, std::vector<RooArgList*> &lists, std::vector<TRooAddPdf_t*> &model) {
  if ((count==0) || (count!=lists.size()) || (count>endings.size())) {
    std::cout << "CreateRooAddPdf<T> for name_base=" << name_base << ": count=" << count << ", list.size=" << lists.size() << ", endings.size=" << endings.size() << "\n";
    return 0;
  }
  model.reserve(count);
  for (unsigned int i=0; i<count; ++i) {
    TString name=name_base+endings[i];
    model.push_back(new TRooAddPdf_t(name,name,*lists[i]));
  }
  return 1;
}

// ---------------------------------------------------------------------

template<class TRooAddModel_t, class TCoefList_t>
  int CreateRooAddModel(unsigned int count, TString name, std::vector<RooArgList*> &lists, std::vector<TCoefList_t*> &coefs, std::vector<TRooAddModel_t*> &model) {
  if ((count==0) || (count!=lists.size()) || (count!=coefs.size())) {
    std::cout << "CreateRooAddModel<T> for name=" << name << ": count=" << count << ", list.size=" << lists.size() << ", coefs.size=" << coefs.size() << "\n";
    return 0;
  }
  model.reserve(count);
  size_t pos=name.Length()-1;
  for (unsigned int i=0; i<count; ++i) {
    name[pos]=char(int('a')+i);
    model.push_back(new TRooAddModel_t(name,name,*lists[i],*coefs[i]));
  }
  return 1;
}
// ---------------------------------------------------------------------

template<class TRooAddModel_t, class TCoefList_t>
int CreateRooAddModel(unsigned int count, TString name_base, const std::vector<TString> &endings, std::vector<RooArgList*> &lists, std::vector<TCoefList_t*> &coefs, std::vector<TRooAddModel_t*> &model) {
  if ((count==0) || (count!=lists.size()) || (count>endings.size()) || (count!=coefs.size())) {
	std::cout << "CreateRooAddModel<T> for name_base=" << name_base << ": count=" << count << ", list.size=" << lists.size() << ", coefs_list.size=" << coefs.size() << ", endings.size=" << endings.size() << "\n";
    return 0;
  }
  model.reserve(count);
  for (unsigned int i=0; i<count; ++i) {
    TString name=name_base+endings[i];
    model.push_back(new TRooAddModel_t(name,name,*lists[i],*coefs[i]));
  }
  return 1;
}




template<class Class_t>
int CreateRooFormulaVar(unsigned int count, TString name, TString formulaString, RooArgList &defArgs, std::vector<Class_t*> &vars, std::vector<RooFormulaVar*> &formula) {
  if (count!=vars.size()) {
    std::cout << "CreateRooFormulaVar<T> for name=" << name << " count=" << count << ", vars.size=" << vars.size() << "\n";
    return 0;
  }
  formula.reserve(count);
  size_t pos=name.Length()-1;
  char buf[30];
  for (unsigned int i=0; i<count; ++i) {
    sprintf(buf,"%s_%d",name.Data(),i);    
    RooArgList args(defArgs);
    args.add(*vars[i]);
    formula.push_back(new RooFormulaVar(buf,formulaString,args));
  }
  return 1;
}

template<class Class1_t, class Class2_t>
int CreateRooFormulaVar(unsigned int count, TString name, TString formulaString, RooArgList &defArgs, std::vector<Class1_t*> &var1, std::vector<Class2_t*> &var2, std::vector<RooFormulaVar*> &formula) {
  if ((count!=var1.size()) || (count!=var2.size())) {
    std::cout << "CreateRooFormulaVar<T,T> for name=" << name << " count=" << count << ", var1.size=" << var1.size() << ", var2.size=" << var2.size() << "\n";
    return 0;
  }
  formula.reserve(count);
  size_t pos=name.Length()-1;
  char buf[30];
  for (unsigned int i=0; i<count; ++i) {
    sprintf(buf,"%s_%d",name.Data(),i);    
    RooArgList args(defArgs);
    args.add(*var1[i]); args.add(*var2[i]);
    formula.push_back(new RooFormulaVar(buf,formulaString,args));
  }
  return 1;
}

template<class Resol_t, class Shape_t, class Conv_t>  // e.g. RooGaussModel, RooHistoPdf & RooNumConvPdf; RooGaussian, RooFormula & RooFFTConvPdf
int CreateRooConvPdf(unsigned int count, TString name, RooRealVar &x, std::vector<Resol_t*> &resol, std::vector<Shape_t*> &shape, std::vector<Conv_t*> &convolution) {
  if ((count!=resol.size()) || (count!=shape.size()))  {
    std::cout << "CreateConvPdf<T> for name=" << name << ": count=" << count << ", resol.size=" << resol.size() << ", shape.size=" << shape.size() << "\n";
    return 0;
  }
  convolution.reserve(count);
  size_t pos=name.Length()-1;
  for (unsigned int i=0; i<count; ++i) {
    name[pos]=char(int('a')+i);
    convolution.push_back(new Conv_t(name,name, x,*resol[i],*shape[i]));
  }
  return 1;
}

// --------------------------------------------------------------------
// --------------------------------------------------------------------

/*
inline
int FitMyModel(int model, const std::vector<double> *data, const DataInfos_t &info) {
  int res=0;
  //RooRealVar x("x","x",0);
  //RooRealVar y("y","y",0);
  //RooDataSet *ds=ConvertToDataSet(data,x,y);

  switch(model) {
    //case cFit1Gauss: res=Fit1GaussModel(ds,x,y,info); break;
    //case cFit2Gauss: res=Fit2GaussModel(ds,x,y,info); break;
    //case cFit2GaussCheby: res=Fit2GaussChebyModel(ds,x,y,info); break;
    //case cFitTestGaussRF102: res=FitTestGaussRF102(ds,x,y,info); break;
    //case cFitTestGaussRF102Demo: res=FitTestGaussRF102Demo(ds,x,y,info); break;
    //case cFit1GaussTest: res=Fit1GaussModelTest(ds,x,y,info); break;
    //case cFitTestGaussRF102M: res=FitTestGaussRF102M(ds,x,y,info); break;
    //case cFitTestExplicit1: res=FitTestExplicit1(ds,x,y,info); break;
    //case cFitTestExplicit2: res=FitTestExplicit2(ds,x,y,info); break;
    //case cFitExplicit2: res=FitExplicit2(ds,x,y,info); break;
    //case cFitBreitWigner: res=FitBreitWigner(ds,x,y,info); break;
    //case cFitBreitWignerEff: res=FitBreitWignerEff(ds,x,y,info); break;
  default:
    std::cout << "FitMyModel is not ready for model=" << model << "\n";
  }
  if (!res) std::cout << " FitMyModel encountered a problem\n";
  return (res) ? 1:0;
}
*/
// --------------------------------------------------------------------
/*
inline
int FitMyMultiModel(int model, const std::vector<std::vector<double>*> &data, const DataInfos_t &info) {
  int res=0;
  //RooRealVar x("x","x",0);
  //RooRealVar y("y","y",0);
  std::cout << "FitMyMultiModel(" << model << ", Vector)" << std::endl;
  //std::vector<RooDataSet*> *ds=ConvertToDataSet(data,x,y);
  //std::cout << "converted to vector<RooDataSet>" << std::endl;

  switch(model) {
    //case cFitBreitWigner: res=FitBreitWignerV(*ds,x,y,info); break;
    //case cFitBreitWignerPF: res=FitBreitWignerVPF(*ds,x,y,info); break;
    //case cFitBreitWignerPFStem: res=FitBreitWignerVPFStem(*ds,x,y,info); break;
#ifdef enableFitRebin
    //case cFitBreitWignerPFStemT: res=FitBreitWignerVPFStemT(*ds,x,y,info); break;
#endif
    //case cFitBreitWignerPFStemT2: res=FitBreitWignerVPFStemT2(*ds,*ds,x,y,info); break;
    //case cFitBreitWignerM1: res=FitBreitWignerVM1(*ds,x,y,info); break;
    //case cFitBreitWignerM1Eff: res=FitBreitWignerVM1Eff(*ds,x,y,info); break;
    //case cFitBreitWignerIlyaSpecial: res=FitBreitWignerIlyaSpecial(*ds,x,y,info); break;
    //case cFitBreitWignerIlyaLiteral: res=FitBreitWignerIlyaLiteral(*ds,x,y,info); break;
    //case cFitPFStemMCGauss: res=FitPFStemMCGauss(*ds,x,y,info); break;
  default:
    std::cout << "FitMyMultiModel is not ready for model=" << model << "\n";
  }
  if (!res) std::cout << " FitMyMultiModel encountered a problem\n";
  return (res) ? 1:0;
 }
*/
// --------------------------------------------------------------------

#ifndef __myLib__
// current version!!
inline
int FitMyMultiModel(int model, const std::vector<mithep::AccEffData_t*> &dataP, const std::vector<mithep::AccEffData_t*> &dataF, const DataInfos_t &info) {
  int res=0;
  RooRealVar x("x","x",0);
  RooRealVar y("y","y",0);
  std::cout << "FitMyMultiModel(" << model << ", Vector)" << std::endl;
  std::vector<RooDataSet*> *ds=ConvertToDataSet(dataP,dataF,x,y);
  std::cout << "converted to vector<RooDataSet>" << std::endl;

  switch(model) {
    //case cFitBreitWigner: res=FitBreitWignerV(*ds,x,y,info); break;
    //case cFitBreitWignerPF: res=FitBreitWignerVPF(*ds,x,y,info); break;
    //case cFitBreitWignerPFStem: res=FitBreitWignerVPFStem(*ds,x,y,info); break;
#ifdef enableFitRebin
    //case cFitBreitWignerPFStemT: res=FitBreitWignerVPFStemT(*ds,x,y,info); break;
#endif
    //case cFitBreitWignerPFStemT2: res=FitBreitWignerVPFStemT2(*ds,*ds,x,y,info); break;
    //case cFitBreitWignerM1: res=FitBreitWignerVM1(*ds,x,y,info); break;
    //case cFitBreitWignerM1Eff: res=FitBreitWignerVM1Eff(*ds,x,y,info); break;
    case cFitBreitWignerIlyaSpecial: res=FitBreitWignerIlyaSpecial(*ds,x,y,info); break;
    case cFitBreitWignerIlyaLiteral: res=FitBreitWignerIlyaLiteral(*ds,x,y,info); break;
 #ifdef enableFitRebin
    //case cFitBreitWignerIlyaSpecialT: res=FitBreitWignerIlyaSpecialT(*ds,x,y,info); break;
    //case cFitBreitWignerIlyaLiteralT: res=FitBreitWignerIlyaLiteralT(*ds,x,y,info); break;
#endif
    case cFitPFStemMCGauss: res=FitPFStemMCGauss(*ds,x,y,info); break;
  default:
    std::cout << "FitMyMultiModel is not ready for model=" << model << "\n";
  }
  if (!res) std::cout << " FitMyMultiModel(PF) encountered a problem\n";
  return (res) ? 1:0;
}
#endif

// --------------------------------------------------------------------

//TH1F* toyModel1(int nbins, double xmin, double xmax, double gauss_avg, double gauss_stddev, double distortion_amplitude, double distortion_spread, double distortion_center, int count=1000);
//TH1F* toyModel2(int nbins, double xmin, double xmax, double gauss_avg, double gauss_stddev, double distortion_amplitude, double distortion_spread, double distortion_center, int count=1000);


//TH1* mytoyModel1(int nbins, double xmin, double xmax);


//RooDataSet* toyDSModel1(RooRealVar &x, double gauss_avg, double gauss_stddev, double distortion_amplitude, double distortion_spread, double distortion_center, int count);
//RooDataSet* toyDSModel2(RooRealVar &x, double gauss_avg, double gauss_stddev, double distortion_amplitude, double distortion_spread, double distortion_center, int count);


//RooDataSet* mytoyDSModel1(RooRealVar &x, int count);

// --------------------------------------------------------------------
// --------------------------------------------------------------------

//#include "MyFitModels5_inc.h"

//#ifdef enableFitRebin
//#include "MyFitModels5_inc2.h"
//#endif

// ----------------------------------------------------------
// ----------------------------------------------------------

template<class Class_t>
void PrintValues(const char *msg, const std::vector<std::vector<Class_t*>*> &vec) {
  std::cout << msg << "[" << vec.size() << "]: ";
  char buf[30];
  for (unsigned int i=0; i<vec.size(); ++i) {
    sprintf(buf,"%s_%d",msg,i);
    PrintValues(buf,*vec[i]);
  }
  std::cout << "\n";
}

template<class Class_t>
void PrintValues(const char *msg, const std::vector<Class_t*> &vec) {
  std::cout << msg << "[" << vec.size() << "]: ";
  for (unsigned int i=0; i<vec.size(); ++i) {
    std::cout << " " << vec[i]->getVal();
  }
  std::cout << "\n";
}

template<class Class_t>
void PrintPrint(const char *msg, const std::vector<std::vector<Class_t*>*> &vec) {
  std::cout << msg << "[" << vec.size() << "]: ";
  char buf[30];
  for (unsigned int i=0; i<vec.size(); ++i) {
    sprintf(buf,"%s_%d",msg,i);
    PrintPrint(buf,*vec[i]);
  }
  std::cout << "\n";
}

template<class Class_t>
void PrintPrint(const char *msg, const std::vector<Class_t*> &vec) {
  std::cout << msg << "[" << vec.size() << "]: ";
  for (unsigned int i=0; i<vec.size(); ++i) {
    std::cout << " "; vec[i]->Print(); std::cout << "\n";
  }
  //std::cout << "\n";
}

//------------------------------------------------------------------------------------------------------------------------

inline void PrintRooDataSetCounts(const char *msg, const std::vector<RooDataSet*> &arr) {
  if (msg) std::cout << msg;
  std::cout << " vec[" << arr.size() << "] of RooDataSet sizes: ";
  for (unsigned int i=0; i<arr.size(); ++i) {
    std::cout << " " << arr[i]->numEntries();
  }
  std::cout << "\n";
  return;
}


// ----------------------------------------------------------
// ----------------------------------------------------------

// ======================================================================

inline TLatex* CreateChiSquareText(double chi2, int is_ndof=0, int color=kRed+3) {
  char buf[60];
  char format[50];
  char name[20];
  if (is_ndof) sprintf(name,"#chi^{2}/N_{dof}"); else sprintf(name,"#chi^{2}");
  if ((chi2>1e-1) && (chi2<1e5)) sprintf(format,"%s=%c5.2lf",name,'%');
  else if ((chi2<2e-1) && (chi2>1e-4)) sprintf(format,"%s=%c8.5lf",name,'%');
  else sprintf(format,"%s=%c6.2e",name,'%');
  sprintf(buf,format,chi2);
  TLatex *txt=new TLatex(0.2,0.82,buf); 
  txt->SetNDC(kTRUE); txt->SetTextSize(0.04); txt->SetTextColor(color);
  return txt;
}

inline TLatex* CreateChiSquareText(RooPlot* frame, int color=kRed+3) {
  assert(frame);
  TLatex* txt=CreateChiSquareText(frame->chiSquare());
  txt->SetTextColor(color);
  frame->addObject(txt);
  return txt;
}

// ======================================================================


#endif

