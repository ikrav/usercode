#ifndef fitMassScaleFactorsE_C
#define fitMassScaleFactorsE_C


#define __myLib__
#include <TROOT.h>
#include <TH1D.h>
#include <vector>
#include <iostream>
#include <assert.h>

namespace mithep {
  typedef std::vector<double> AccEffData_t;
  typedef TH1D myHistoClass_t;
};

namespace DYTools {
  int massRangeSet=13;
  const int nMassBinsSet11 = 40; // 2011 binning
  const double massBinLimits11[nMassBinsSet11 + 1] = 
    {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 
     81, 86, 91, 96, 101, 106, 110, 115, 120, 126, 133, 141, 
     150, 160, 171, 185, 200, 220, 243, 273, 320, 380, 440, 
     510, 600, 1000, 1500}; // 40 bins
  const double ZMassMin11= 15;
  const double ZMassMax11= 1500;

  // 2011 mass binning Z zone part
  const int nMassBinsSet12 = 13;
  const double massBinLimits12[nMassBinsSet12 + 1] = 
    {60, 64, 68, 72, 76, 
     81, 86, 91, 96, 101, 106, 110, 115, 120}; // 40-27 bins
  const double ZMassMin12= 60;
  const double ZMassMax12= 120;

  // for comparison with RooFit auto-binning
  const int nMassBinsSet13 = int((120-60)/0.6 + 1e-3);
  const double ZMassMin13= 60;
  const double ZMassMax13=120;
};

namespace DYTools {

  //TEnergyScaleFactorsWorkCase_t g_esfWorkCase;
  //TEnergyScaleFactorsFitModel_t g_esfFitModel;
  int scEtaRangeSet;
  int nScEtaBins3;
  double* scEtaBinLimits3;

  inline int NumberOfMassBins() {
    int result =-1;
    switch(massRangeSet) {
    case 11: result=nMassBinsSet11; break;
    case 12: result=nMassBinsSet12; break;
    case 13: result=nMassBinsSet13; break;
     default:
      std::cout << "\n\nDYTools2011::NumberOfMassBins requested for unprepared massRangeSet=" << massRangeSet << "\n";
    }    
    return result;
  }

  inline
  double MassMin() {
    double result =-1;
    switch(massRangeSet) {
    case 11: result=ZMassMin11; break;
    case 12: result=ZMassMin12; break;
    case 13: result=ZMassMin13; break;
    default:
      std::cout << "\n\nDYTools2011::MassMin requested for unprepared massRangeSet=" << massRangeSet << "\n";
    }    
    return result;
  }
  
  inline
  double MassMax() {
    double result =-1;
    switch(massRangeSet) {
    case 11: result=ZMassMax11; break;
    case 12: result=ZMassMax12; break;
    case 13: result=ZMassMax13; break;
    default:
      std::cout << "\n\nDYTools2011::MassMax requested for unprepared massRangeSet=" << massRangeSet << "\n";
    }    
    return result;
  }

  inline void copyMassBins(int count, double *dest, const double *source) {
    assert(dest); assert(source);
    for (int i=0; i<count; ++i) dest[i]=source[i];
  }

  inline 
  double* getMassBins() {
    int n=NumberOfMassBins();
    double *bins=new double[n+1];
    switch(massRangeSet) {
    case 11: copyMassBins(n+1,bins,massBinLimits11); break;
    case 12: copyMassBins(n+1,bins,massBinLimits12); break;
    default:
      double zmin=MassMin();
      double zmax=MassMax();
      double dh=(zmax-zmin)/n;
      for (int i=0; i<n+1; i++) { bins[i]=zmin+i*dh; }
    }
    return bins;
  }
};


#include "ElectronEnergyScaleAdv.C"
#include "MyFitModels.hh"
#include "MyFitModels.cc"

double fit_range_min=60;
double fit_range_max=120;
int no_plot=0;

#include "fitMassScaleFactorsE.hh"
#include "fitMassScaleFactorsE.cc"

// ----------------------------------------------------------------

template<class HistoClass_t>
HistoClass_t* combineEtaEtaHistos(const std::vector<HistoClass_t*> &source, const std::vector<int> &combine, const TString &combination_name) {
  if (combine.size()==0) {
    std::cout << "combineEtaEtaHistos: combine.size=0. Returning NULL\n";
    return NULL;
  }
  HistoClass_t *h=(HistoClass_t*)source[combine[0]]->Clone(combination_name);
  for (unsigned int i=1; i<combine.size(); ++i) {
    h->Add(source[combine[i]]);
  }
  return h;
}

// ----------------------------------------------------------------

template<class HistoClass_t>
int combineEtaEtaHistos(const std::vector<HistoClass_t*> &source, const std::vector<std::vector<int>*> &combine, const std::vector<TString> &combination_names, std::vector<HistoClass_t*> &result) {
  result.clear(); 
  if (combine.size()!=combination_names.size()) {
    std::cout << "combineEtaEtaHistos(V): combine.size=" << combine.size() << ", combination_names.size=" << combination_names.size() << "\n";
    return 0;
  }
  if ((source.size()==0) || (combine.size()==0)) {
    std::cout << "combineEtaEtaHistos(V): source.size=" << source.size() << ", combine.size=" << combine.size() << " (neither can be 0)\n";
    return 0;
  }
  result.reserve(source.size());
  TString xaxis_name=source[0]->GetXaxis()->GetTitle();
  for (unsigned int i=0; i<combine.size(); ++i) {
    HistoClass_t *tmpH=combineEtaEtaHistos(source, *combine[i], combination_names[i]);
    if (!tmpH) {
      std::cout << "failed to create a combination in combineEtaEtaHistos(V)\n";
      return 0;
    }
    tmpH->GetXaxis()->SetTitle(xaxis_name.Data());
    result.push_back(tmpH);
  }
  return 1;
}

// ----------------------------------------------------------------

#include "Test.hh"


int LoadFile(const char *filename, std::vector<std::vector<double>*> &massMC, std::vector<std::vector<double>*> &et1_mc, std::vector<std::vector<double>*> &et2_mc, std::vector<std::vector<double>*> &massExp, std::vector<std::vector<double>*> &et1_exp, std::vector<std::vector<double>*> &et2_exp) {
  std::cout << "loading <" << filename << ">" << std::endl;
  FILE *fin=fopen(filename,"r");
  if (!fin) {
    std::cout << "failed to open the file with masses\n";
    return 0;
  }
  const int bufsize=100;
  char buf[bufsize];
  fgets(buf,bufsize,fin);
  fgets(buf,bufsize,fin);
  unsigned int mcCount=0;
  fscanf(fin,"%u",&mcCount); fgets(buf,bufsize,fin);
  massMC.reserve(mcCount);
  std::vector<double> *vec=NULL;
  for (unsigned int mci=0; mci<mcCount; ++mci) {
    unsigned int t=0;
    fscanf(fin,"%u",&t);
    //std::cout << "loading mci=" << mci << ", expects " << t << " counts\n";
    vec=new std::vector<double>();
    massMC.push_back(vec);
    vec->reserve(t);
    double x;
    for (unsigned int i=0; i<t; ++i) {
      fscanf(fin,"%lf",&x);
      vec->push_back(x);
    }
    //std::cout << "last mass value=" << x << "\n";
  }
  fgets(buf,bufsize,fin);
  //std::cout << "buf=" << buf << "\n";
  fgets(buf,bufsize,fin);
  //std::cout << "buf=" << buf << "\n";
  unsigned int expCount=0;
  fscanf(fin,"%u",&expCount); fgets(buf,bufsize,fin);
  //std::cout << "expecting " << expCount << " exp data vecs\n";
  massExp.reserve(expCount);
  for (unsigned int expi=0; expi<expCount; ++expi) {
    unsigned int t=0;
    fscanf(fin,"%u",&t);
    vec=new std::vector<double>();
    massExp.push_back(vec);
    vec->reserve(t);
    double x;
    for (unsigned int i=0; i<t; ++i) {
      fscanf(fin,"%lf",&x);
      vec->push_back(x);
    }
  }
  fclose(fin);
  std::cout << "loaded " << massMC.size() << " MC and " << massExp.size() << " exp mass vectors\n";
  if (massMC.size()<10) {
    const unsigned int bin_count=(unsigned int)(trunc(sqrt(2*massMC.size())));
    std::cout << "massMC sizes:\n";
    unsigned int k=0;
    for (unsigned int i=0; i<bin_count; ++i) {
      for (unsigned int j=0; j<i; j++) printf("          ");
      for (unsigned int j=i; j<bin_count; ++j, ++k) printf(" %9d",int(massMC[k]->size()));
      printf("\n");
    }
  }
  if (massExp.size()<10) {
    const unsigned int bin_count=(unsigned int)(trunc(sqrt(2*massExp.size())));
    std::cout << "massExp sizes:\n";
    unsigned int k=0;
    for (unsigned int i=0; i<bin_count; ++i) {
      for (unsigned int j=0; j<i; j++) printf("          ");
      for (unsigned int j=i; j<bin_count; ++j, ++k) printf(" %9d",int(massExp[k]->size()));
      printf("\n");
    }
  }
  if (0) { PrintVVecCounts("et1_mc ",et1_mc); PrintVVecCounts("et2_mc ",et2_mc); }
  if (0) { PrintVVecCounts("et1_exp ",et1_exp); PrintVVecCounts("et2_exp ",et2_exp); }
  
  return ((massMC.size()==massExp.size()) && (massMC.size()>0)) ? 1:0;
}


// ---------------------------------------------------

int CreateComparisonPlots(ElectronEnergyScaleAdv_t &sf, const std::vector<std::vector<double>*> &massMC, const std::vector<std::vector<double>*> &massExp, std::vector<std::string> &html_lines, RooRealVar &mass, ScaledSmearedDataAux_t &aux, std::vector<RooAddPdf*> &mcModelV, std::vector<RooDataSet*> &expModelV, const std::vector<std::vector<int>*> &combine, const std::vector<TString> &combination_names) {
  int expect_count=sf.ExpectCount();
  int massBinCount=DYTools::NumberOfMassBins();
  double *massBins=DYTools::getMassBins();

  std::vector<mithep::myHistoClass_t*> unsmearedMCHV,unscaledExpHV;
  std::vector<mithep::myHistoClass_t*> smearedMCHV,scaledExpHV;

  std::cout << "aux.fitModel=" << aux.fitModel << "\n";
  std::cout << "calling sf.PrepareScaledSmearedData(0,0)" << std::endl;
  if (!sf.PrepareScaledSmearedData(massMC,massExp,massBinCount,massBins,unsmearedMCHV,unscaledExpHV,"NoSF",0,0)) {
    std::cout << "failed with unprocessed data collections\n";
    return 0;
  }
  std::cout << "passed" << std::endl;
  if (0 && (sf.FitModelKind()==_ESFModel_gauss)) {
    html_lines.push_back("<br />Using per-event smearing<br>\n");
    if (sf.FitModelKind()==_ESFModel_gauss) {
      if (!sf.PrepareScaledSmearedData(massMC,massExp,massBinCount,massBins,smearedMCHV,scaledExpHV,"SF",1,1)) {
	std::cout << "failed with basic data collections\n";
	return 0;
      }
    }
  }
  else {
    html_lines.push_back("<br />Using template smearing<br>\n");
    char buf[50];
    mithep::myHistoClass_t *hbase=new mithep::myHistoClass_t("hbase","hbase",massBinCount,massBins);
    hbase->GetXaxis()->SetTitle("mass");
    hbase->GetYaxis()->SetTitle("count");
    hbase->GetYaxis()->SetTitleOffset(1.6);
    smearedMCHV.reserve(expect_count); scaledExpHV.reserve(expect_count);
    for (int i=0; i<expect_count; ++i) {
      if (i<10) sprintf(buf,"smearedMC_0%d",i); else sprintf(buf,"smearedMC_%d",i);
      mithep::myHistoClass_t *hmc=NULL;
      if (massMC[i]->size()>0) {
	hmc=CreateHisto(*mcModelV[i],buf,&mass,massBinCount,massBins,hbase);
	double scale=unsmearedMCHV[i]->Integral()/hmc->Integral();
	hmc->Scale(scale);
      }
      else {
	std::cout << "ScanAndScale: data contains no points. Creating an empty histo for " << buf << "\n";
	hmc=(mithep::myHistoClass_t*)hbase->Clone(buf);
      }
      smearedMCHV.push_back(hmc);
    }
    for (int i=0; i<expect_count; ++i) {
      if (i<10) sprintf(buf,"smearedExp_0%d",i); else sprintf(buf,"smearedExp_%d",i);
      mithep::myHistoClass_t *hexp=ConvertToHisto(*expModelV[i],buf,hbase);
      scaledExpHV.push_back(hexp);
    }
    delete hbase;
  }

  std::string file_tag="cmp";
  const int print_discrepancy_table=0;
  const int print_descriptions_separately=0;
  if (!PrepareComparisonTablesAndPlots(file_tag.c_str(),unsmearedMCHV,unscaledExpHV,smearedMCHV,scaledExpHV, combine,combination_names, print_discrepancy_table, print_descriptions_separately, histoDir_USER.c_str(), html_lines)) {
    std::cout << "failed in PrepareComparisonTablesAndPlots\n";
    return 0;
  }
  return 1;
}

#include "CPlotMdf.cc"

// ---------------------------------------------------


#endif
