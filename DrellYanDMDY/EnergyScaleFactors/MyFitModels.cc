#include "MyFitModels.hh"
#include <map>

#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooFitResult.h"
#include "RooDataHist.h"
#include "RooCategory.h"
#include "RooMinuit.h"
#include "RooGenericPdf.h"

#include "RooChebychev.h"
#include "RooPolyVar.h"
#include <RooNumConvPdf.h>
#include "TPaveText.h"

#include "TAxis.h"
#include "RooPlot.h"
#include "TText.h"
#include "TLatex.h"

#include "TFile.h"
#include "HelpingTools.hh"


using RooFit::LineColor;
using RooFit::LineStyle;
using RooFit::LineWidth;

int saveCFormat=1;
std::string saveFigFileTag="";

#ifdef __myLib__
extern int no_plot;
extern double fit_range_min;
extern double fit_range_max;
const int fit_difference=0;
const std::vector<int> fit_gated;
const std::vector<double> fit_gated_bin;
#endif

//=== FUNCTION DECLARATIONS ======================================================================================

#ifndef __myLib__
int FitModelAllowed(int fit_model, int label_size, int fname_size, int pass_data_size, int fail_data_size) {
  return 1;
}
#endif



//====================================================================================================

int cPlotTheFit=1;
int cDumpTheFitContents=0;
int cSaveTheFit=0;
int cExtraCanvasIdx=0;



// ======================================================================
// ======================================================================

double gaussRandom(double mean, double stddev) {
  double accum=0.;
  for (int i=0; i<12; ++i) accum+=rand();
  return (accum/double(RAND_MAX)-6.)*stddev+mean;
}

// ======================================================================
// ======================================================================

void PrepareNumCharVec(unsigned int count, std::vector<char> &idxV) {
  idxV.clear();
  idxV.reserve(count);
  for (unsigned int i=0; i<count; i++) {
    char c=char('0'+i+1);
    if (i>8) c=char('A'+i-9);
    if (c>'Z') c=char('a'+i-9-'Z'+'A'-1);
    idxV.push_back(c);
  }  
}

void PrepareCharVec(unsigned int count, std::vector<char> &idxV) {
  idxV.clear();
  idxV.reserve(count);
  for (unsigned int i=0; i<count; i++) {
    char c=char('A'+i);
    if (c>'Z') c=char('a'+i-'Z'+'A'-1);
    if (c>'z') {
      std::cout << "i=" << i << ". switched\n";
      c=char(192+4+i-'z');
    }
    std::cout << " (" << i << ")" << c; // << "\n";
    idxV.push_back(c);
  }  
  std::cout << "\n";
}

void PrepareNameEndings(unsigned int count, std::vector<char> &idxV, std::vector<TString> &nameEndings) {
  unsigned int N=count*(count+1)/2;
  nameEndings.reserve(N);
  if (idxV.size()<N) {
    if (N<9) PrepareNumCharVec(N,idxV);
    else PrepareCharVec(N,idxV);
  }
  unsigned int k=0;
  char buff[5];
  buff[2]='_'; buff[4]='\0';
  for (unsigned int i=0; i<count; ++i) {
    for (unsigned int j=i; j<count; ++j, ++k) {
      buff[0]=idxV[i]; buff[1]=idxV[j]; buff[3]=idxV[k];
      nameEndings.push_back(TString(buff));
    }
  }
}

void PrepareNumCharVec(unsigned int count, std::vector<TString> &idxV) {
  idxV.clear();
  idxV.reserve(count);
  int shift=0;
  char buf[5];
  buf[0]='0'; buf[1]='0'; buf[2]='0'; buf[3]='0'; buf[4]='\0';
  if (count>=9) shift=1;
  else if (count>=99) shift=2;
  else if (count>=999) shift=3;
  for (unsigned int i=0; i<count; i++) {
    if (shift) {
      if ((i==9) || (i==99) || (i==999)) shift--;
    }
    sprintf(buf+shift,"%d",i+1);
    idxV.push_back(TString(buf));
  }  
}

void PrepareCharVec(unsigned int count, std::vector<TString> &idxV) {
  idxV.clear();
  idxV.reserve(count);
  char buf[3];
  buf[2]='\0';
  for (unsigned int i=0; i<count; i++) {
    char c1=char('a'+(i/24));
    char c2=char('a'+(i%24));
    buf[0]=c1; buf[1]=c2;
    idxV.push_back(TString(buf));
  }
  std::cout << "\n";
}

void PrepareNameEndings(unsigned int count, std::vector<TString> &idxV, std::vector<TString> &nameEndings) {
  int N=count*(count+1)/2;
  nameEndings.reserve(N);
  PrepareCharVec(N,idxV);
  unsigned int k=0;
  char buff[8];
  buff[4]='_'; buff[7]='\0';
  for (unsigned int i=0; i<count; ++i) {
    const char *b1=idxV[i].Data();
    for (unsigned int j=i; j<count; ++j, ++k) {
      const char *b2=idxV[j].Data();
      const char *b3=idxV[k].Data();
      buff[0]=b1[0]; buff[1]=b1[1]; 
      buff[2]=b2[0]; buff[3]=b2[1];
      buff[5]=b3[0]; buff[6]=b3[1];
      nameEndings.push_back(TString(buff));
    }
  }
}


// ======================================================================
// ======================================================================

RooDataSet* ConvertToDataSet(const std::vector<double> *dt, RooRealVar &x, RooRealVar &y) {
  std::cout << "Converting vector of doubles to a RooDataSet\n";
  Double_t xmin=1e9, xmax=-1e9;
  Double_t ymin=1e9, ymax=-1e9;
  for (std::vector<double>::const_iterator it=dt->begin(); it!=dt->end(); it++) {
    if (*it < xmin) xmin=*it;
    if (*it > xmax) xmax=*it;
  }
  std::cout << "ConvertToDataSet deduced: xmin=" << xmin << ", xmax=" << xmax << "\n";
  if (fit_range_min>=0) {
    xmin=fit_range_min;
    std::cout << " forcing xmin to fit_range_min=" << fit_range_min << "\n";
  }
  if (fit_range_max>=0) {
    xmax=fit_range_max;
    std::cout << " forcing xmax to fit_range_max=" << fit_range_max << "\n";
  }
  ymin=xmin; ymax=xmax;

  x.setRange(xmin,xmax);
  //y.setRange(ymin,ymax);
  RooArgSet value(x);
  RooDataSet *d = new RooDataSet("d","d",RooArgSet(x));
  for (std::vector<double>::const_iterator it=dt->begin(); it!=dt->end(); ++it) {
    if ((*it>=xmin) && (*it<=xmax)) {
      x.setVal(*it);
      d->add(value);
    }
  }
  if (0) y.setRange(ymin,ymax);
  return d;
}



// ----------------------------------------------------------------------


std::vector<RooDataSet*>* ConvertToDataSet(const std::vector<std::vector<double>*> &dtV, RooRealVar &x, RooRealVar &y, int reverse_order) {
  std::cout << "Converting vector<vector of doubles> to a vector<RooDataSet*> (fit_difference=" << fit_difference << ")\n";

  Double_t xmin=1e9, xmax=-1e9;
  Double_t ymin=1e9, ymax=-1e9;
  for (unsigned int i=0; i<dtV.size(); ++i) {
    const std::vector<double> *dt= dtV[i];
    for (std::vector<double>::const_iterator it=dt->begin(); it!=dt->end(); it++) {
      if (*it < xmin) xmin=*it;
      if (*it > xmax) xmax=*it;
    }
  }
  std::cout << "ConvertToDataSet(vector) deduced: xmin=" << xmin << ", xmax=" << xmax << "\n";
  if (fit_range_min>=0) {
    xmin=fit_range_min;
    std::cout << " forcing xmin to fit_range_min=" << fit_range_min << "\n";
  }
  if (fit_range_max>=0) {
    xmax=fit_range_max;
    std::cout << " forcing xmax to fit_range_max=" << fit_range_max << "\n";
  }
  ymin=xmin; ymax=xmax;

  x.setRange(xmin,xmax);
  //y.setRange(ymin,ymax);
  RooArgSet value(x);
  
  std::vector<RooDataSet*>* dd=new std::vector<RooDataSet*>();
  assert(dd);
  dd->reserve(dtV.size());
  char name[5];
  std::vector<unsigned int> counts;
  counts.reserve(dtV.size());
  for (unsigned int ii=0; ii<dtV.size(); ++ii) {
    unsigned int i=(reverse_order) ? dtV.size()-1-ii : ii;
    std::cout << "i=" << i << ", dtV[i].size=" << dtV[i]->size() << "\n";
    unsigned int count=0;
    sprintf(name,"d%u",i);
    RooDataSet *d = new RooDataSet(name,name,RooArgSet(x));
    const std::vector<double> *dt=dtV[i];
    std::cout << "i=" << i << "\n";
    for (std::vector<double>::const_iterator it=dt->begin(); it!=dt->end(); ++it) {
      if ((*it>=xmin) && (*it<=xmax)) {
	x.setVal(*it);
	d->add(value);
	count++;
	//std::cout << " " << *it;
      }
    }
    //std::cout << "\n";
    dd->push_back(d);
    counts.push_back(count);
  }

  std::cout << "counts : "; for (unsigned int i=0; i<counts.size(); ++i) { std::cout << " " << counts[i]; } std::cout << "\n";
  if (0) y.setRange(ymin,ymax);
  return dd;
}


// ======================================================================
// ----------------------------------------------------------------------
// ======================================================================

std::vector<RooDataSet*>* ConvertToDataSet(const std::vector<const std::vector<double>*> &dtV, std::vector<std::string> &names, RooRealVar &x, RooRealVar &y) {
  std::cout << "Converting vector<const vector of doubles> to a vector<RooDataSet*>\n";
  if (dtV.size()!=names.size()) {
    std::cout << " dtV.size=" << dtV.size() << ", names.size=" << names.size() << std::endl;
    return NULL;
  }
  if (!dtV.size()) {
    std::cout << "ConvertToDataSet(vector of vectors): no data provided (dtV.size()=0)\n";
    return NULL;
  }
  std::cout << "dtV.size=" << dtV.size() << std::endl;
  Double_t xmin=1e9, xmax=-1e9;
  Double_t ymin=1e9, ymax=-1e9;
  for (unsigned int i=0; i<dtV.size(); ++i) {
    const std::vector<double> *dt= dtV[i];
    for (std::vector<double>::const_iterator it=dt->begin(); it!=dt->end(); it++) {
      if (*it < xmin) xmin=*it;
      if (*it > xmax) xmax=*it;
    }
  }
  std::cout << "ConvertToDataSet(vector) deduced: xmin=" << xmin << ", xmax=" << xmax << "\n";
  if (fit_range_min>=0) {
    xmin=fit_range_min;
    std::cout << " forcing xmin to fit_range_min=" << fit_range_min << "\n";
  }
  if (fit_range_max>=0) {
    xmax=fit_range_max;
    std::cout << " forcing xmax to fit_range_max=" << fit_range_max << "\n";
  }
  //std::cout << "here" << std::endl;
  ymin=xmin; ymax=xmax;

  //std::cout << "setting the range of x to " << xmin << " .. " << xmax << std::endl;
  x.setRange(xmin,xmax);
  //y.setRange(ymin,ymax);
  std::cout << "range set" << std::endl;

  RooArgSet value(x);
  
  //std::cout << "here" << std::endl;

  std::vector<RooDataSet*>* dd=new std::vector<RooDataSet*>();
  assert(dd);
  dd->reserve(dtV.size());
  std::vector<unsigned int> counts;
  counts.reserve(dtV.size());
  for (unsigned int ii=0; ii<dtV.size(); ++ii) {
    //std::cout << "ii=" << ii << std::endl;
    //unsigned int i=dtV.size()-1-ii;
    unsigned int i=ii;
    if (!dtV[i]) std::cout << "null dtV[" << i << "]" << std::endl;
    //std::cout << "i=" << i << ", dtV[i].size=" << dtV[i]->size() << std::endl;
    unsigned int count=0;
    const char *name=names[i].c_str();
    RooDataSet *d = new RooDataSet(name,name,RooArgSet(x));
    const std::vector<double> *dt=dtV[i];
    //std::cout << "i=" << i << "\n";
    for (std::vector<double>::const_iterator it=dt->begin(); it!=dt->end(); ++it) {
      if ((*it>=xmin) && (*it<=xmax)) {
	x.setVal(*it);
	d->add(value);
	count++;
	//std::cout << " " << *it;
      }
    }
    //std::cout << "\n";
    dd->push_back(d);
    counts.push_back(count);
  }

  std::cout << "counts : "; for (unsigned int i=0; i<counts.size(); ++i) { std::cout << " " << counts[i]; } std::cout << "\n";
  if (0) y.setRange(ymin,ymax);
  return dd;
}

// ======================================================================
// ----------------------------------------------------------------------
// ======================================================================

#ifndef __myLib__
std::vector<RooDataSet*>* ConvertToDataSet(const std::vector<mithep::AccEffData_t*> &dataP, const std::vector<mithep::AccEffData_t*> &dataF, RooRealVar &x, RooRealVar &y) {
  std::cout << "Converting vector<AccEffData> to a vector<RooDataSet*>\n";
  std::vector<const std::vector<double>*> sel;
  std::vector<std::string> names;
  std::vector<RooDataSet*> *rset, *rsetP, *rsetF;
  RooRealVar x1,y1;
  int gated=(fit_gated.size()) ? 1:0;
  std::cout << "gated=" << gated << "\n";
  std::vector<mithep::AccEffData_t*> gatedP;
  std::vector<mithep::AccEffData_t*> gatedF;
  if (gated) {
    gatedP.reserve(dataP.size());
    gatedF.reserve(dataF.size());
    std::cout << "Gating is on:\n";
    for (unsigned int i=0; i<fit_gated.size(); ++i) {
      switch(fit_gated[i]) {
      case 1: std::cout << "eta gate " << DYTools::scEtaBinMin(fit_gated_bin[i]) << "..." << DYTools::scEtaBinMax(fit_gated_bin[i]) << "\n"; break;
      case 2: std::cout << "Et gate " << DYTools::scEtBinMin(fit_gated_bin[i]) << "..." << DYTools::scEtBinMax(fit_gated_bin[i]) << "\n"; break;
      default:
	std::cout << "ConvertToDataSet(PF) is not ready for fit_gated[i]=" << fit_gated[i] << std::endl;
	return NULL;
      }
    }
    for (unsigned int i=0; i<dataP.size(); ++i) {
      mithep::AccEffData_t *dt=new mithep::AccEffData_t();
      dt->FillGated(*dataP[i],fit_gated,fit_gated_bin);
      gatedP.push_back(dt);
    }
    for (unsigned int i=0; i<dataF.size(); ++i) {
      mithep::AccEffData_t *dt=new mithep::AccEffData_t();
      dt->FillGated(*dataF[i],fit_gated,fit_gated_bin);
      gatedF.push_back(dt);
    }
    if ((gatedP.size()==0) || (gatedF.size()==0)) {
      std::cout << "gated data size was 0 (pass count=" << gatedP.size() << ", fail count=" << gatedF.size() << ")\n";
      return NULL;
    }
  }
  const std::vector<mithep::AccEffData_t*> *passP=(gated) ? &gatedP : &dataP;
  const std::vector<mithep::AccEffData_t*> *passF=(gated) ? &gatedF : &dataF;
  for (int pass=0; pass<2; pass++) {
    const std::vector<mithep::AccEffData_t*> *data=(pass) ? passP : passF;
    std::cout << "pass=" << pass << ", data->size=" << data->size() << std::endl;
    names.clear();
    sel.clear();
    for (unsigned int i=0; i<data->size(); ++i) {
      names.push_back(std::string((*data)[i]->HRefMZName));
      sel.push_back(& (*data)[i]->HRefMZ);
    }
    if (pass) rsetP=ConvertToDataSet(sel,names,x,y);
    else rsetF=ConvertToDataSet(sel,names,x,y);
  }
  std::cout << "loop ended" << std::endl;
  if (rsetF) {
    rset=rsetF;
    if (rsetP) {
      rset->reserve(rsetF->size() + rsetP->size());
      rset->insert(rset->end(), rsetP->begin(), rsetP->end());
      delete rsetP;
    }
  }
  else {
    rset=rsetP;
  }
  if (gatedP.size()) for (unsigned int i=0; i<gatedP.size(); ++i) delete gatedP[i];
  if (gatedF.size()) for (unsigned int i=0; i<gatedF.size(); ++i) delete gatedF[i];
  std::cout << "counts : "; for (unsigned int i=0; i<rset->size(); ++i) { std::cout << " " << (*rset)[i]->numEntries(); } std::cout << "\n";
  return rset;
}
#endif

// ======================================================================

#ifndef __myLib__
std::vector<std::vector<double>*>* ConvertMassToVectorDbl(const std::vector<mithep::AccEffData_t*> &dataP, const std::vector<mithep::AccEffData_t*> &dataF) {
  std::cout << "Converting vector<AccEffData> to a vector<vector<dbl>*>\n";
  std::vector<std::vector<double>*> selP;
  std::vector<std::vector<double>*> selF;
  int gated=(fit_gated.size()) ? 1:0;
  std::cout << "gated=" << gated << "\n";
  std::vector<mithep::AccEffData_t*> gatedP;
  std::vector<mithep::AccEffData_t*> gatedF;
  if (gated) {
    gatedP.reserve(dataP.size());
    gatedF.reserve(dataF.size());
    std::cout << "Gating is on:\n";
    for (unsigned int i=0; i<fit_gated.size(); ++i) {
      switch(fit_gated[i]) {
      case 1: std::cout << "eta gate " << DYTools::scEtaBinMin(fit_gated_bin[i]) << "..." << DYTools::scEtaBinMax(fit_gated_bin[i]) << "\n"; break;
      case 2: std::cout << "Et gate " << DYTools::scEtBinMin(fit_gated_bin[i]) << "..." << DYTools::scEtBinMax(fit_gated_bin[i]) << "\n"; break;
      default:
	std::cout << "ConvertMassToVectorDbl(PF) is not ready for fit_gated[i]=" << fit_gated[i] << std::endl;
	return NULL;
      }
    }
    for (unsigned int i=0; i<dataP.size(); ++i) {
      mithep::AccEffData_t *dt=new mithep::AccEffData_t();
      dt->FillGated(*dataP[i],fit_gated,fit_gated_bin);
      gatedP.push_back(dt);
    }
    for (unsigned int i=0; i<dataF.size(); ++i) {
      mithep::AccEffData_t *dt=new mithep::AccEffData_t();
      dt->FillGated(*dataF[i],fit_gated,fit_gated_bin);
      gatedF.push_back(dt);
    }
    if ((gatedP.size()==0) || (gatedF.size()==0)) {
      std::cout << "gated data size was 0 (pass count=" << gatedP.size() << ", fail count=" << gatedF.size() << ")\n";
      return NULL;
    }
  }
  const std::vector<mithep::AccEffData_t*> *passP=(gated) ? &gatedP : &dataP;
  const std::vector<mithep::AccEffData_t*> *passF=(gated) ? &gatedF : &dataF;
  for (int pass=0; pass<2; pass++) {
    const std::vector<mithep::AccEffData_t*> *data=(pass) ? passP : passF;
    std::vector<std::vector<double>*> *sel=(pass) ? &selP : &selF;
    std::cout << "pass=" << pass << ", data->size=" << data->size() << std::endl;
    for (unsigned int i=0; i<data->size(); ++i) {
      sel->push_back(new std::vector<double>((*data)[i]->HRefMZ));
    }
  }
  std::cout << "loop ended" << std::endl;
  std::vector<std::vector<double>*> *sel=new std::vector<std::vector<double>*>();
  if (selP.size()) sel->insert(sel->end(),selP.begin(),selP.end());
  if (selF.size()) sel->insert(sel->end(),selF.begin(),selF.end());
  if (gatedP.size()) for (unsigned int i=0; i<gatedP.size(); ++i) delete gatedP[i];
  if (gatedF.size()) for (unsigned int i=0; i<gatedF.size(); ++i) delete gatedF[i];
  return sel;
}
#endif

// ======================================================================
// ======================================================================

#ifndef __myLib__

int CreatePassFailStemFormula(int count, std::vector<RooRealVar*> &eff, RooRealVar &nSigEv, std::vector<RooFormulaVar*> &passf) {
  if ((count==0) || (count!=eff.size()+1)) {
    std::cout << "CreatePassFailStemFormula: count=" << count << ", eff.size=" << eff.size() << "\n";
    return 0;
  }
  passf.reserve(count);
  switch(count) {
  case 1:
    passf.push_back(new RooFormulaVar("pass","pass","nSigEv",RooArgList(nSigEv)));
    break;
  case 2: {
    passf.push_back(new RooFormulaVar("fail_1","fail_1","(1-eff_a) * nSigEv",RooArgList(*eff[0],nSigEv)));
    passf.push_back(new RooFormulaVar("pass_1","pass_1","eff_a * nSigEv",RooArgList(*eff[0],nSigEv)));
  }
    break;
  case 3: {
    passf.push_back(new RooFormulaVar("fail_1","fail_1","(1-eff_a) * nSigEv",RooArgList(*eff[0],nSigEv)));
    passf.push_back(new RooFormulaVar("fail_2","fail_2","eff_a * (1-eff_b) * nSigEv",RooArgList(*eff[0],*eff[1],nSigEv)));
    passf.push_back(new RooFormulaVar("pass_2","pass_2","eff_a * eff_b * nSigEv",RooArgList(*eff[0],*eff[1],nSigEv)));
  }
    break;
  case 4: {
    passf.push_back(new RooFormulaVar("fail_1","fail_1","(1-eff_a) * nSigEv",RooArgList(*eff[0],nSigEv)));
    passf.push_back(new RooFormulaVar("fail_2","fail_2","eff_a * (1-eff_b) * nSigEv",RooArgList(*eff[0],*eff[1],nSigEv)));
    passf.push_back(new RooFormulaVar("fail_3","fail_3","eff_a * eff_b * (1-eff_c) * nSigEv",RooArgList(*eff[0],*eff[1],*eff[2],nSigEv)));
    passf.push_back(new RooFormulaVar("pass_3","pass_3","eff_a * eff_b * eff_c * nSigEv",RooArgList(*eff[0],*eff[1],*eff[2],nSigEv)));
  }
    break;
  default:
    std::cout << "CreatePassFailStemFormula preparation is not ready for count=" << count << "\n";
    return 0;
  }
  return 1;
}

// --------------------------------------------------------------------

int SetInitialPassFailEffValues(int count, const std::vector<RooDataSet*> &data, std::vector<RooRealVar*> &eff, int &total_event_count) {
  if ((count!=data.size()) || (count!=eff.size()+1)) {
    std::cout << "SetInitialPassFaillEffValues: count=" << count << ", data.size()=" << data.size() << ", eff.size=" << eff.size() << "\n";
    return 0;
  }
  // Signal events are bound to 1 number
  std::vector<int> cntEvt(count+1);
  cntEvt.back() = 0;
  for (unsigned int i=0; i<data.size(); ++i) {
    unsigned int idx=data.size()-1-i;
    cntEvt[idx]= cntEvt[idx+1] + data[idx]->numEntries();
  }
  if (1) {
    std::cout << "got event counts:\n";
    for (unsigned int i=0; i<data.size(); ++i) {
      std::cout << " sample i=" << i << " has " << data[i]->numEntries() << " entries\n";
    }
    for (unsigned int i=0; i<cntEvt.size(); ++i) {
      std::cout << " sum until sample i=" << i << ": " << cntEvt[i] << "\n";
    }
    //std::cout << "\n";
  }

  total_event_count=cntEvt[0]; 
  for (unsigned int i=0; i<eff.size(); ++i) {
    eff[i]->setVal(1-0.9*cntEvt[i]/double(cntEvt[0]));
  }

  std::cout << "total_event_count=" << total_event_count << "\n";
  std::cout << "eff="; for (unsigned int i=0; i<eff.size(); i++) std::cout << " " << eff[i]->getVal(); std::cout << "\n";
  return 1;
}
#endif

// --------------------------------------------------------------------

RooDataSet* CreateCombinedData(TString name,TString descr, const std::vector<RooDataSet*> &data, RooRealVar &x, RooCategory &sample, const std::vector<TString> &categories) {
  unsigned int count=data.size();
  if ((count!=categories.size())) {
    std::cout << "CreateCombinedData for " << name << ": data.size=" << data.size() << ", categories.size=" << categories.size() << "\n";
    return NULL;
  }
  std::map<std::string,RooDataSet*> sets_by_caths;
  for (unsigned int i=0; i<count; ++i) {
    sets_by_caths.insert(std::pair<std::string,RooDataSet*>(std::string(categories[i]),data[i]));
  }
  RooDataSet *combData=new RooDataSet(name,descr, x, RooFit::Index(sample),RooFit::Import(sets_by_caths));
  return combData;
}

// --------------------------------------------------------------------

RooSimultaneous* CreateCombinedModel(TString name, TString descr, std::vector<RooAddPdf*> &model, RooCategory &sample, const std::vector<TString> &categories) {
  unsigned int count=model.size();
  if ((count!=categories.size())) {
    std::cout << "CreateCombinedModel for " << name << ": model.size=" << model.size() << ", categories.size=" << categories.size() << "\n";
    return 0;
  }
  
  RooSimultaneous *simPdf=new RooSimultaneous(name,descr,sample);
  for (unsigned int i=0; i<count; ++i) {
    simPdf->addPdf(*model[i],categories[i]);
  }
  return simPdf;
}

// ======================================================================
// ======================================================================

// --------------------------------------------------------------------

int CreateRooCategories(unsigned int count, TString name_base, const std::vector<TString> &name_endings, RooCategory &sample, std::vector<TString> &categories) {
  if (count<name_endings.size()) {
    std::cout << "CreateRooCategories: got fewer endings (" << name_endings.size() << ") than needed (" << count << ")\n";
    return 0;
  }
  categories.clear();
  TString name;
  for (unsigned int i=0; i<count; ++i) {
    name=name_base+name_endings[i];
    sample.defineType(name);
    categories.push_back(name);
  }
  return 1;
}

// --------------------------------------------------------------------

int CreateRooCategories(unsigned int count, TString name, const std::vector<char> &idxV, RooCategory &sample, std::vector<TString> &categories) {
  if (count<idxV.size()) {
    std::cout << "CreateRooCategories: got fewer endings (" << idxV.size() << ") than needed (" << count << ")\n";
    return 0;
  }
  categories.clear();
  size_t pos=name.Length()-1;
  for (unsigned int i=0; i<count; ++i) {
    name[pos]=idxV[i];
    sample.defineType(name);
    categories.push_back(name);
  }
  return 1;
}


// ----------------------------------------------------------------------

//#error x
int CreateRooRealVar(unsigned int count, TString name, std::vector<RooRealVar*> &vars, double start_value, double min_value, double max_value) {
  vars.reserve(count);
  size_t pos=name.Length()-1;
  if (min_value<-9998) {
    for (unsigned int i=0; i<count; ++i) {
      name[pos]=char(int('a')+i);
      vars.push_back(new RooRealVar(name,name,start_value));
    }
  }
  else if (max_value<-9998) {
    for (unsigned int i=0; i<count; ++i) {
      name[pos]=char(int('a')+i);
      vars.push_back(new RooRealVar(name,name,start_value,min_value));
    }
  }
  else {
    for (unsigned int i=0; i<count; ++i) {
      name[pos]=char(int('a')+i);
      vars.push_back(new RooRealVar(name,name,start_value,min_value,max_value));
    }
  }
  return 1;
}

// ======================================================================

int CreateRooRealVar(unsigned int count, const TString &name_base, const std::vector<TString> &name_endings, std::vector<RooRealVar*> &vars, double start_value, double min_value, double max_value) {
  if (name_endings.size()!=count) {
    std::cout << "CreateRooRealVar(count,name,name_endings): got count=" << count << " and " << name_endings.size() << " name endings\n";
    return 0;
  }
  vars.reserve(count);
  TString name;
  if (min_value<-9998) {
    for (unsigned int i=0; i<count; ++i) {
      name=name_base+name_endings[i];
      vars.push_back(new RooRealVar(name,name,start_value));
    }
  }
  else if (max_value<-9998) {
    for (unsigned int i=0; i<count; ++i) {
      name=name_base+name_endings[i];
      vars.push_back(new RooRealVar(name,name,start_value,min_value));
    }
  }
  else {
    for (unsigned int i=0; i<count; ++i) {
      name=name_base+name_endings[i];
      vars.push_back(new RooRealVar(name,name,start_value,min_value,max_value));
    }
  }
  return 1;
}

// ---------------------------------------------------------------------

int CreateRooExponential(unsigned int count, TString name, RooRealVar &x, std::vector<RooRealVar*> &tau, std::vector<RooExponential*> &expos) {
  if (count!=tau.size()) {
    std::cout << "CreateRooExponential for name=" << name << ": count=" << count << ", tau.size=" << tau.size() << "\n";
    return 0;
  }
  expos.reserve(tau.size());
  size_t pos=name.Length()-1;
  for (unsigned int i=0; i<tau.size(); ++i) {
    name[pos]=char(int('a')+i);
    expos.push_back(new RooExponential(name,name, x,*tau[i]));
  }
  return 1;
}

// ---------------------------------------------------------------------


int CreateRooHistPdf(const RooRealVar &x, const std::vector<RooDataSet*> &data, TString name1, std::vector<RooDataHist*> &histos, TString name2, std::vector<RooHistPdf*> &hpdfs, int from_idx, int to_idx) {
  int N=data.size();
  if (from_idx<0) from_idx=0;
  if (to_idx<0) to_idx=N;
  if ((from_idx>N) || (to_idx>N) || (to_idx-from_idx<=0)) {
    std::cout << "CreateRooHistPdf(RooDataSet) for names=" << name1 << "and " << name2 << ": from_idx=" << from_idx << ", to_idx=" << to_idx << "\n";
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
      histos.push_back(data[i]->binnedClone(name1,name1));
      hpdfs.push_back(new RooHistPdf(name2,name2,RooArgSet(x),*histos.back()));
    }
  }
  /*
  else {
    for (unsigned int i=from_idx; i<to_idx; i++) {
      RooDataHist *dh=data[i]->binnedClone(name1,name1);
	use double *massBins_loc=DYTools::getMassBins();
	HRefMZ = new TH1F(HRefMZName,HRefMZName, DYTools::NumberOfMassBins(), massBins_loc);
	delete [] massBins_loc;
      TH1F *h=dh->createHistogram("h",DYTools::NumberOfMassBins(),DYTools::MassMin(),DYTools::MassMax());
      const char c=char(int('a')-from_idx+i);
      name1[pos1]=c; name2[pos2]=c;
      histos.push_back(new RooDataHist(name1,name1,RooArgSet(x),h));
      hpdfs.push_back(new RooHistPdf(name2,name2,RooArgSet(x),*histos.back()));
    }
  }
  */
  return 1;
}

// ---------------------------------------------------------------------

int CreateRooAddPdf(unsigned int count, TString name, std::vector<RooExtendPdf*> &var1, std::vector<RooExtendPdf*> &var2, std::vector<RooAddPdf*> &model) {
  if ((count==0) || (count!=var1.size()) || (count!=var2.size())) {
    std::cout << "CreateRooAddPdf for name=" << name << ": count=" << count << ", var1.size=" << var1.size() << ", var2.size=" << var2.size() << "\n";
    return 0;
  }
  model.reserve(count);
  size_t pos=name.Length()-1;
  for (unsigned int i=0; i<count; ++i) {
    name[pos]=char(int('a')+i);
    model.push_back(new RooAddPdf(name,name,RooArgList(*var1[i],*var2[i])));
  }
  return 1;
}

// ======================================================================
// ======================================================================


/*  
int Fit1GaussModel(const RooDataSet *ds, RooRealVar &x, RooRealVar &y, const DataInfos_t &info) {
  std::cout << "Doing Fit1GaussModel\n";
  Double_t  xavg = 0.5 * (x.getMin() + x.getMax());

  RooDataSet *data=new RooDataSet(*ds);

  x.Print();
  data->Print();

  // Create one Gaussian PDFs g1(x,mean1,sigma) and their parameters
  Double_t Mee_center=90;
  RooRealVar mean1("mean1","mean of gaussian 1",Mee_center,20,420); //0.8*Mee_center) ;
  RooRealVar sigma1("sigma1","width of gaussian 1",3,0.2,80); //2);

  RooGaussian sig("sig","Signal",x,mean1,sigma1) ;

  RooRealVar a("a","a",-1e9,1e9);
  a.setVal(0.5);
  
  TCanvas* c = (cPlotTheFit) ? (new TCanvas("the_fit","the_fit",800,400)) : NULL ;
  RooPlot *xframe1 = (cPlotTheFit) ? x.frame(RooFit::Title("Initial frame")) : NULL ;
  RooPlot *xframe2 = (cPlotTheFit) ? x.frame(RooFit::Title("Fitted frame")) : NULL ;

  if (cPlotTheFit) {
    data->plotOn(xframe1,RooFit::LineColor(info.color(0)),RooFit::LineStyle(info.linestyle(0))) ;
    sig.plotOn(xframe1,RooFit::LineColor(info.color(1)),RooFit::LineStyle(info.linestyle(1)));
    c->Divide(2);
    c->cd(1) ; gPad->SetLeftMargin(0.15) ; xframe1->GetYaxis()->SetTitleOffset(1.6) ; xframe1->Draw() ;
    c->Update();
  }

  //sig.fitTo(data);

  // Sum the signal and background 
  //RooAddPdf  model("model","g*a",RooArgList(sig),RooArgList(a)) ;
  //RooGaussian model(sig);
  
  // Perform fit and save result
  RooFitResult* fit_result = (fit_dont_fit) ? NULL : sig.fitTo(*data) ;
  //RooFitResult* fit_result = model.chi2FitTo(data) ;
  //RooFitResult* fit_result = model.chi2FitTo(data,Save()) ; // rf602
  // Summary printing: Basic info plus final values of floating fit parameters
  
  if (cPlotTheFit && 1) {
    data->plotOn(xframe2,RooFit::LineColor(info.color(0)),RooFit::LineStyle(info.linestyle(0))) ;
    sig.plotOn(xframe2,RooFit::LineColor(info.color(1)),RooFit::LineStyle(info.linestyle(1)));
    //bkg.plotOn(xframe2,RooFit::LineColor(info.color(2)),RooFit::LineStyle(info.linestyle(2)));

    //model.plotOn(xframe2,RooFit::LineColor(info.color(3)),RooFit::LineStyle(info.linestyle(3))) ;

    // Overlay the background component of model with a dashed line
    //model.plotOn(xframe2,RooFit::Components(sig),RooFit::LineColor(kGreen+2),RooFit::LineStyle(kDashed)) ;
    
    c->cd();
    c->cd(2) ; gPad->SetLeftMargin(0.15) ; xframe2->GetYaxis()->SetTitleOffset(1.6) ; xframe2->Draw() ;
    c->Update();
    if (cSaveTheFit) c->SaveAs("the_fit_1G.png");
    std::cout << "updated c: enter a char...\n"; char xxx; std::cin >> xxx;
  }
  
  if (cDumpTheFitContents) {
    fit_result->Print() ;
	  
    // Verbose printing: Basic info, values of constant paramaters, initial and
    // final values of floating parameters, global correlations
    fit_result->Print("v") ;
    
    // Open new ROOT file save save result 
    TFile f("rf607_fitresult.root","RECREATE") ;
    fit_result->Write("rf607") ;
    f.Close() ;
  }
    
  if (0 && cPlotTheFit) {
    // Plot data and PDF overlaid
    RooPlot* xframe = x.frame(RooFit::Title("Example of composite pdf=(sig1+sig2)+bkg")) ;
    data->plotOn(xframe) ;
    //model.plotOn(xframe) ;
    
    c->cd();
    gPad->SetLeftMargin(0.15) ; xframe->GetYaxis()->SetTitleOffset(1.6) ; xframe->Draw() ;

    c->Update();
    std::cout << "updated c 2nd: enter a char...\n"; char xxx; std::cin >> xxx;
    }

  if (cDumpTheFitContents) {
    // Print structure of composite p.d.f.
    //model.Print("t") ;
  }
  return 1;
}
*/
// ======================================================================
/*  
int Fit2GaussModel(const RooDataSet *ds, RooRealVar &x, RooRealVar &y, const DataInfos_t &info) {
  std::cout << "Doing Fit2GaussModel\n";
  RooDataSet data(*ds);

  Double_t  xavg = 0.5 * (x.getMin() + x.getMax());

  // Create a binned dataset that imports contents of TH1 and associates its contents to observable 'x'
  const char *data_name="my_histo";
  RooCategory cat("cat","cat");
  cat.defineType("SampleA");

  // Create one Gaussian PDFs g1(x,mean1,sigma) and their parameters
  Double_t Mee_center=90;
  RooRealVar mean1("mean1","mean of gaussian 1",Mee_center) ;
  RooRealVar sigma1("sigma1","width of gaussian 1",2);
  RooRealVar mean2("mean2","mean of gaussian 1",0.95*xavg) ;
  RooRealVar sigma2("sigma2","width of gaussian 1",40);

  RooGaussian sig("sig","Signal",x,mean1,sigma1) ;
  RooGaussian bkg("bkg","Background",x,mean2,sigma2) ;
  
  TCanvas* c = (cPlotTheFit) ? (new TCanvas("the_fit","the_fit",800,400)) : NULL ;
  RooPlot *xframe1 = (cPlotTheFit) ? x.frame(RooFit::Title("Initial frame")) : NULL ;
  RooPlot *xframe2 = (cPlotTheFit) ? x.frame(RooFit::Title("Fitted frame")) : NULL ;

  if (cPlotTheFit) {
    data.plotOn(xframe1,RooFit::LineColor(info.color(0)),RooFit::LineStyle(info.linestyle(0))) ;
    sig.plotOn(xframe1,RooFit::LineColor(info.color(1)),RooFit::LineStyle(info.linestyle(1)));
    bkg.plotOn(xframe1,RooFit::LineColor(info.color(2)),RooFit::LineStyle(info.linestyle(2)));
    c->Divide(2);
    c->cd(1) ; gPad->SetLeftMargin(0.15) ; xframe1->GetYaxis()->SetTitleOffset(1.6) ; xframe1->Draw() ;
    c->Update();
  }

  // Sum the signal and background 
  RooRealVar bkgfrac("bkgfrac","fraction of background",0.1,0.,.2) ;
  RooAddPdf  model("model","g+a",RooArgList(bkg,sig),bkgfrac) ;
  
  // Perform fit and save result
  RooFitResult* fit_result = (fit_dont_fit) ? NULL : model.fitTo(data) ;
  //RooFitResult* fit_result = model.chi2FitTo(data) ;
  //RooFitResult* fit_result = model.chi2FitTo(data,Save()) ; // rf602
  // Summary printing: Basic info plus final values of floating fit parameters
  
  if (cPlotTheFit && 1) {
    data.plotOn(xframe2,RooFit::LineColor(info.color(0)),RooFit::LineStyle(info.linestyle(0))) ;
    //sig.plotOn(xframe2,RooFit::LineColor(info.color(1)),RooFit::LineStyle(info.linestyle(1)));
    //bkg.plotOn(xframe2,RooFit::LineColor(info.color(2)),RooFit::LineStyle(info.linestyle(2)));

    //model.plotOn(xframe2,RooFit::LineColor(info.color(3)),RooFit::LineStyle(info.linestyle(3))) ;

    // Overlay the background component of model with a dashed line
    model.plotOn(xframe2,RooFit::Components(sig),RooFit::LineColor(kGreen+2),RooFit::LineStyle(kDashed)) ;
    model.plotOn(xframe2,RooFit::Components(bkg),RooFit::LineColor(kBlue),RooFit::LineStyle(kDashed)) ;
    
    // Overlay the background+sig components of model with a dotted line
    //model.plotOn(xframe2,RooFit::Components(RooArgSet(bkg,sig)),RooFit::LineColor(info.color(4)),RooFit::LineStyle(kDotted)) ;
    model.plotOn(xframe2,RooFit::Components(RooArgSet(bkg,sig)),RooFit::LineColor(kRed),RooFit::LineStyle(kDotted)) ;

    c->cd();
    c->cd(2) ; gPad->SetLeftMargin(0.15) ; xframe2->GetYaxis()->SetTitleOffset(1.6) ; xframe2->Draw() ;
    c->Update();
    if (cSaveTheFit) c->SaveAs("the_fit.png");
    std::cout << "updated c: enter a char...\n"; char xxx; std::cin >> xxx;
  }
  
  if (cDumpTheFitContents) {
    fit_result->Print() ;
	  
    // Verbose printing: Basic info, values of constant paramaters, initial and
    // final values of floating parameters, global correlations
    fit_result->Print("v") ;
    
    // Open new ROOT file save save result 
    TFile f("rf607_fitresult.root","RECREATE") ;
    fit_result->Write("rf607") ;
    f.Close() ;
  }
    
  if (0 && cPlotTheFit) {
    // Plot data and PDF overlaid
    RooPlot* xframe = x.frame(RooFit::Title("Example of composite pdf=(sig1+sig2)+bkg")) ;
    data.plotOn(xframe) ;
    model.plotOn(xframe) ;
    
    // Overlay the background component of model with a dashed line
    model.plotOn(xframe,RooFit::Components(bkg),RooFit::LineStyle(kDashed)) ;
    
    // Overlay the background+sig components of model with a dotted line
    model.plotOn(xframe,RooFit::Components(RooArgSet(bkg,sig)),RooFit::LineStyle(kDotted)) ;
    

    TCanvas* c = new TCanvas("rf101_basics","rf101_basics",800,400) ;
    c->cd();
    gPad->SetLeftMargin(0.15) ; xframe->GetYaxis()->SetTitleOffset(1.6) ; xframe->Draw() ;

    c->Update();
    std::cout << "updated c 2nd: enter a char...\n"; char xxx; std::cin >> xxx;
  }

  if (cDumpTheFitContents) {
    // Print structure of composite p.d.f.
    model.Print("t") ;
  }

  return 1;
}
*/
// ======================================================================
/*  
int Fit2GaussChebyModel(const RooDataSet *ds, RooRealVar &x, RooRealVar &y, const DataInfos_t &info) {
  std::cout << "Doing Fit2GaussChebyModel\n";

  Double_t  xavg = 0.5 * (x.getMin() + x.getMax());

  // Create one Gaussian PDFs g1(x,mean1,sigma) and their parameters
  Double_t Mee_center=90;
  RooRealVar mean1("mean1","mean of gaussian 1",Mee_center,20.,500.) ;
  RooRealVar sigma1("sigma1","width of gaussian 1",2,0.,100.);
  RooRealVar mean2("mean2","mean of gaussian 1",0.7*xavg,20.,500.) ;
  RooRealVar sigma2("sigma2","width of gaussian 1",5,0.,100.);
  RooGaussian sig1("sig1","Signal component 1",x,mean1,sigma1) ;
  RooGaussian sig2("sig2","Signal component 2",x,mean2,sigma2) ;
  
  // Build Chebychev polynomial p.d.f.  
  RooRealVar a0("a0","a0",0.1,-1.,1.) ;
  RooRealVar a1("a1","a1",-0.1,-1.,1.) ;
  RooChebychev bkg("bkg","Background",x,RooArgSet(a0,a1)) ;
  RooPolyVar bkg_poly("bkg_poly","Background(poly)",x,RooArgSet(a0,a1));
  // A d d   s i g n a l   c o m p o n e n t s 
  // ------------------------------------------
  
  // Sum the signal components into a composite signal p.d.f.
  RooRealVar sig1frac("sig1frac","fraction of component 1 in signal",0.8,0.7,1.) ;
  RooAddPdf sig("sig","Signal",RooArgList(sig1,sig2),sig1frac) ;
  // Sum the composite signal and background 
  RooRealVar bkgfrac("bkgfrac","fraction of background",0.05,0.,0.1) ;
  RooAddPdf  model("model","g1+g2+a",RooArgList(bkg,sig),bkgfrac) ;
  
  // Perform fit and save result
  RooDataSet* data=new RooDataSet(*ds);
  RooFitResult* fit_result = NULL; 
  //RooFitResult* fit_result = model.chi2FitTo(data,RooFit::Save()) ; // rf602
  // Summary printing: Basic info plus final values of floating fit parameters
  
  if (cPlotTheFit) {
    RooPlot *xframe = x.frame(RooFit::Title("Gaussian p.d.f.")) ;
    data->plotOn(xframe);
    sig1.plotOn(xframe,LineColor(kRed));
    sig2.plotOn(xframe,LineColor(kGreen+1));
    bkg.plotOn(xframe);
    fit_result=(fit_dont_fit) ? NULL : model.fitTo(*data) ;
    RooPlot *xframe2 = x.frame(RooFit::Title("Gaussian p.d.f. with data"));
    model.plotOn(xframe2,LineStyle(kDashed),LineColor(kBlack),LineWidth(1));
    data->plotOn(xframe2);
    sig1.plotOn(xframe2,LineColor(kRed));
    sig2.plotOn(xframe2,LineColor(kGreen+1));
    bkg.plotOn(xframe2);

    TCanvas* c = new TCanvas("rf101_basics","rf101_basics",800,400) ;
    c->Divide(2) ;
    c->cd(1) ; gPad->SetLeftMargin(0.15) ; xframe->GetYaxis()->SetTitleOffset(1.6) ; xframe->Draw() ;
    c->cd(2) ; gPad->SetLeftMargin(0.15) ; xframe2->GetYaxis()->SetTitleOffset(1.6) ; xframe2->Draw() ;
    c->Update();
    c->SaveAs("fit.png"); std::cout << "saved as fit.png\n";
    std::cout << "updated c: enter a char...\n"; char xxx; std::cin >> xxx;
  }
  else { fit_result=(fit_dont_fit) ? NULL : model.fitTo(*data); }
  
  if (cDumpTheFitContents && fit_result) {
    fit_result->Print() ;
	  
    // Verbose printing: Basic info, values of constant paramaters, initial and
    // final values of floating parameters, global correlations
    fit_result->Print("v") ;
    
    // Open new ROOT file save save result 
    TFile f("rf607_fitresult.root","RECREATE") ;
    fit_result->Write("rf607") ;
    f.Close() ;
  }
    
  if (0 && cPlotTheFit) {
    // Plot data and PDF overlaid
    RooPlot* xframe = x.frame(RooFit::Title("Example of composite pdf=(sig1+sig2)+bkg")) ;
    data->plotOn(xframe) ;
    model.plotOn(xframe) ;
    
    // Overlay the background component of model with a dashed line
    model.plotOn(xframe,RooFit::Components(bkg),RooFit::LineStyle(kDashed),RooFit::LineColor(kGray)) ;
    model.plotOn(xframe,RooFit::Components(sig1),RooFit::LineStyle(kDashed),RooFit::LineColor(kRed)) ;
    
    // Overlay the background+sig2 components of model with a dotted line
    model.plotOn(xframe,RooFit::Components(RooArgSet(bkg,sig2)),RooFit::LineStyle(kDotted),RooFit::LineColor(kGreen+4)) ;
    

    TCanvas* c = new TCanvas("rf101_basics","rf101_basics",800,400) ;
    c->cd();
    gPad->SetLeftMargin(0.15) ; xframe->GetYaxis()->SetTitleOffset(1.6) ; xframe->Draw() ;

    c->Update();
    c->SaveAs("fit.png"); std::cout << "saved as fit.png\n";
    std::cout << "updated c 2nd: enter a char...\n"; char xxx; std::cin >> xxx;
  }

  if (cDumpTheFitContents) {
    // Print structure of composite p.d.f.
    model.Print("t") ;
  }

  return 1;
}
*/
// ======================================================================

/*
int FitTestGaussRF102(const RooDataSet *dummy, RooRealVar &x, RooRealVar &y, const DataInfos_t &info) {
  std::cout << "FitTestGaussRF102\n";
  double avg=120.;
  double stddev=10.5;
  double alpha=0.005;
  double ampl=0.2;
  double center=avg*0.8;
  int count=1000;

  RooDataSet* dh=toyDSModel2(x,avg,stddev,ampl,alpha,center,count);
  RooDataSet* dh1=toyDSModel1(x,avg,stddev,0,alpha,center,count);

  RooPlot* frame = x.frame(RooFit::Title("toy model")) ;

  dh->plotOn(frame,LineColor(info.color(600)),LineStyle(1)) ; 
  //dh1->plotOn(frame,LineColor(info.color(400)),LineStyle(2)) ; 

  // Fit a Gaussian p.d.f to the data
  RooRealVar mean("mean","mean",100,10,1000) ;
  RooRealVar sigma("sigma","sigma",3,0.1,10000) ;
  RooGaussian gauss("gauss","gauss",x,mean,sigma) ;
  gauss.fitTo(*dh) ;
  gauss.plotOn(frame) ;
  std::cout << "gauss: "; gauss.Print();
  std::cout << "  - mean : "; mean.Print();
  std::cout << "  - sigma: "; sigma.Print();
  
    TCanvas* c4 = new TCanvas("rf101_basics","rf101_basics",800,400) ;
    c4->cd();
    gPad->SetLeftMargin(0.15) ; frame->GetYaxis()->SetTitleOffset(1.6) ; frame->Draw() ;

    c4->Update();
    c4->SaveAs("fit.png");
    //std::cout << "updated c4 : enter a char...\n"; char xxx; std::cin >> xxx;

  return 1;
}
*/

// ======================================================================

/*
int FitTestGaussRF102Demo(const RooDataSet *dummy, RooRealVar &x, RooRealVar &y, const DataInfos_t &info) {
  std::cout << "FitTestGaussRF102Demo\n";
  double avg=120.;
  double stddev=10.5;
  double alpha=0.005;
  double ampl=0.008;
  double center=avg*0.8;
  int count=1000;

  RooDataSet* dh1=toyDSModel1(x,avg,stddev,0,alpha,center,count);
  RooDataSet* dh2=toyDSModel2(x,avg,stddev,ampl,alpha,center,count);
  RooDataSet* dh3=toyDSModel2(x,avg,stddev,0.015,alpha,center,count);

  RooPlot* frame1 = x.frame(RooFit::Title("gauss toy model")) ;
  RooPlot* frame2 = x.frame(RooFit::Title("a bit distorted gauss t.m."));
  RooPlot* frame3 = x.frame(RooFit::Title("distorted gauss toy model"));

  dh1->plotOn(frame1,LineColor(info.color(600)),LineStyle(1)) ; 
  dh2->plotOn(frame2,LineColor(info.color(400)),LineStyle(2)) ; 
  dh3->plotOn(frame3,LineColor(info.color(400)),LineStyle(2)) ; 

  // Fit a Gaussian p.d.f to the data
  RooRealVar mean1("mean1","mean1",100,10,1000) ;
  RooRealVar sigma1("sigma1","sigma1",3,0.1,10000) ;
  RooGaussian gauss1("gauss1","gauss1",x,mean1,sigma1) ;
  gauss1.fitTo(*dh1) ;
  gauss1.plotOn(frame1) ;
  std::cout << "gauss1: "; gauss1.Print();
  std::cout << "  - mean1 : "; mean1.Print();
  std::cout << "  - sigma1: "; sigma1.Print();
  
  RooRealVar mean2("mean2","mean2",100,10,1000) ;
  RooRealVar sigma2("sigma2","sigma2",3,0.1,10000) ;
  RooGaussian gauss2("gauss2","gauss2",x,mean2,sigma2) ;
  gauss2.fitTo(*dh2) ;
  gauss2.plotOn(frame2) ;
  std::cout << "gauss2: "; gauss2.Print();
  std::cout << "  - mean2 : "; mean2.Print();
  std::cout << "  - sigma2: "; sigma2.Print();
  
  RooRealVar mean3("mean3","mean3",100,10,1000) ;
  RooRealVar sigma3("sigma3","sigma3",3,0.1,10000) ;
  RooGaussian gauss3("gauss3","gauss3",x,mean3,sigma3);
  gauss3.fitTo(*dh3) ;
  gauss3.plotOn(frame3) ;
  std::cout << "gauss3: "; gauss3.Print();
  std::cout << "  - mean3 : "; mean3.Print();
  std::cout << "  - sigma3: "; sigma3.Print();
  
    TCanvas* c4 = new TCanvas("rf101_basics","rf101_basics",1200,400) ;
    c4->Divide(3);
    c4->cd(1); gPad->SetLeftMargin(0.15) ; frame1->GetYaxis()->SetTitleOffset(1.6) ; frame1->Draw() ;
    c4->cd(2); gPad->SetLeftMargin(0.15) ; frame2->GetYaxis()->SetTitleOffset(1.6) ; frame2->Draw() ;
    c4->cd(3); gPad->SetLeftMargin(0.15) ; frame3->GetYaxis()->SetTitleOffset(1.6) ; frame3->Draw() ;

    c4->Update();
    c4->SaveAs("fit.png");
    //std::cout << "updated c4 : enter a char...\n"; char xxx; std::cin >> xxx;

  return 1;
}
*/

// ======================================================================

/*
int FitTestGaussRF102M(const RooDataSet *dummy, RooRealVar &x, RooRealVar &y, const DataInfos_t &info) {
  std::cout << "FitTestGaussRF102M\n";
  double avg=120.;
  double stddev=10.5;
  double alpha=0.005;
  double ampl=0.008;
  double center=avg*0.8;

  RooDataSet* hd=toyDSModel2(x,avg,stddev,ampl,alpha,center, 1000);
  std::cout << "calling to Fit2GaussChebyModel\n";
  return Fit2GaussChebyModel(hd,x,y,info);
}
*/
  
 // ======================================================================
/*
int Fit1GaussModelTest(const RooDataSet *ds, RooRealVar &x, RooRealVar &y, const DataInfos_t &info) {
  std::cout << "Fit1GaussModelTest\n";

  /
  RooDataHist dh("dh","dh",x,RooFit::Import(*hd)) ;
  // Make plot of binned dataset showing Poisson error bars (RooFit default)
  RooPlot* frame = x.frame(RooFit::Title("Imported TH1 with Poisson error bars")) ;
  dh.plotOn(frame) ; 

  // Fit a Gaussian p.d.f to the data
  RooRealVar mean("mean","mean",100,-10,10000) ;
  RooRealVar sigma("sigma","sigma",3,0.1,10000) ;
  RooGaussian gauss("gauss","gauss",x,mean,sigma) ;
  gauss.fitTo(dh) ;
  gauss.plotOn(frame) ;
  std::cout << "gauss: "; gauss.Print();
  std::cout << "  - mean : "; mean.Print();
  std::cout << "  - sigma: "; sigma.Print();
  
    TCanvas* c4 = new TCanvas("rf101_basics","rf101_basics",800,400) ;
    c4->cd();
    gPad->SetLeftMargin(0.15) ; frame->GetYaxis()->SetTitleOffset(1.6) ; frame->Draw() ;

    c4->Update();
    c4->SaveAs("fit.png");
    //std::cout << "updated c4 : enter a char...\n"; char xxx; std::cin >> xxx;
    * /
  return 1;
}
*/

// ----------------------------------------------------------
/*
int FitExplicit1(const RooDataSet *ds, RooRealVar &x, RooRealVar &y, const DataInfos_t &info) {
  std::cout << "entered FitExplicit1. Based on rf601_intminuit.C example. Uses extended fit\n";
  RooDataSet data(*ds);

  // S e t u p   p d f   a n d   l i k e l i h o o d 
  // -----------------------------------------------

  // Observable
  //RooRealVar x("x","x",-20,20) ;

  // Model (intentional strong correlations)
  double ini_mean1=70, ini_sigma1=3;
  double ini_mean2=100, ini_sigma2=4;
  RooRealVar mean1("mean1","mean of g1",ini_mean1, 50,140) ;
  RooRealVar sigma_g1("sigma_g1","width of g1",ini_sigma1,0.5,100) ;
  RooGaussian g1("g1","g1",x,mean1,sigma_g1) ;
  RooRealVar mean2("mean2","mean of g2",ini_mean2, 50,200) ;
  RooRealVar sigma_g2("sigma_g2","width of g2",ini_sigma2,3.0,6.0) ;
  RooGaussian g2("g2","g2",x,mean2,sigma_g2) ;

  RooRealVar frac("frac","frac",0.5,0.0,1.0) ;
  RooAddPdf model("model","model",RooArgList(g1,g2),frac) ;

  RooRealVar nsig1("nsig1","nsignal1",500,0,1e7);
  RooExtendPdf esig1("esig1","esig1",g1,nsig1);
  RooRealVar nsig2("nsig2","nsignal2",500,0,1e7);
  RooExtendPdf esig2("esig2","esig2",g2,nsig2);
  RooAddPdf emodel("emodel","emodel",RooArgList(esig1,esig2));


  // Generate 1000 events
  //RooDataSet* data = model.generate(x,1000) ;
  
  // Construct unbinned likelihood of model w.r.t. data
  RooAbsReal* nll = model.createNLL(data) ;

  // I n t e r a c t i v e   m i n i m i z a t i o n ,   e r r o r   a n a l y s i s
  // -------------------------------------------------------------------------------

  // Create MINUIT interface object
  RooMinuit m(*nll) ;

  // Activate verbose logging of MINUIT parameter space stepping
  m.setVerbose(kTRUE) ;

  // Call MIGRAD to minimize the likelihood
  m.migrad() ;

  // Print values of all paramaters, that reflect values (and error estimates)
  // that are back propagated from MINUIT
  //model.getParameters(x)->Print("s") ;

  // Disable verbose logging
  m.setVerbose(kFALSE) ;

  // Run HESSE to calculate errors from d2L/dp2
  m.hesse() ;

  // Print value (and error) of sigma_g2 parameter, that reflects
  // value and error back propagated from MINUIT
  //sigma_g2.Print() ;

  // Run MINOS on sigma_g2 parameter only
  m.minos(); //sigma_g2) ;

  // Print value (and error) of sigma_g2 parameter, that reflects
  // value and error back propagated from MINUIT
  //sigma_g2.Print() ;

  // make a plot
  RooPlot* frame = x.frame(RooFit::Title("explicit non-extended fit") );
  data.plotOn(frame) ; 
  model.plotOn(frame);
  RooPlot *frame2 = x.frame(RooFit::Title("implicit extended fit") );
  data.plotOn(frame2);
  mean1.setVal(ini_mean1); mean1.setRange(50,140);
  sigma_g1.setVal(ini_sigma1); sigma_g1.setRange(0.5,100);
  mean2.setVal(ini_mean2); mean2.setRange(50,200);
  sigma_g2.setVal(ini_sigma2); sigma_g2.setRange(3,6);
  emodel.fitTo(data,RooFit::Extended(true));
  emodel.plotOn(frame2);

    TCanvas* c4 = new TCanvas("rf101_basics","rf101_basics",800,400) ;
    c4->Divide(2);
    c4->cd(1);
    gPad->SetLeftMargin(0.15) ; frame->GetYaxis()->SetTitleOffset(1.6) ; frame->Draw() ;
    c4->cd(2);
    gPad->SetLeftMargin(0.15) ; frame2->GetYaxis()->SetTitleOffset(1.6) ; frame2->Draw() ;

    c4->Update();
    c4->SaveAs("fit.png");
    //std::cout << "updated c4 : enter a char...\n"; char xxx; std::cin >> xxx;
  return 1;
}
*/


// ----------------------------------------------------------
/*
int FitExplicit2(const RooDataSet *ds, RooRealVar &x, RooRealVar &y, const DataInfos_t &info) {
  std::cout << "entered FitExplicit2. Simplified FitExplicit1. Two Gaussian fit\n";
  RooDataSet data(*ds);

  // S e t u p   p d f   a n d   l i k e l i h o o d 
  // -----------------------------------------------

  // Observable

  RooRealVar mean1("mean1","mean of g1",70, 0,140) ;
  RooRealVar sigma_g1("sigma_g1","width of g1",3,0.5,20) ;
  RooGaussian g1("g1","g1",x,mean1,sigma_g1) ;
  RooRealVar mean2("mean2","mean of g2",140, 50,200) ;
  RooRealVar sigma_g2("sigma_g2","width of g2",4,0.5,600.0) ;
  RooGaussian g2("g2","g2",x,mean2,sigma_g2) ;

  RooRealVar nsig1("nsig1","nsignal1",500,0,1e7);
  RooExtendPdf esig1("esig1","esig1",g1,nsig1);
  RooRealVar nsig2("nsig2","nsignal2",500,0,1e7);
  RooExtendPdf esig2("esig2","esig2",g2,nsig2);
  RooAddPdf emodel("emodel","emodel",RooArgList(esig1,esig2));

  if (1) {
    mean1.setVal(90); mean1.setRange(80,100);
    mean2.setVal(40); //mean2.setRange(40,40); //mean2.setRange(0,500);    
    sigma_g2.setVal(50); sigma_g2.setRange(50,50);
  }

  // Generate 1000 events
  //RooDataSet* data = model.generate(x,1000) ;
  
  RooPlot* frame = x.frame(RooFit::Title("starting frame") );
  data.plotOn(frame) ; 
  //emodel.plotOn(frame);
  RooPlot *frame2 = x.frame(RooFit::Title("extended fit") );
  data.plotOn(frame2);

  if (1) {
    const int print_level=0; // -1 is minimal
    std::cout << "default call to fitTo\n";
    if (1) {
      emodel.fitTo(data,RooFit::Extended(true)); //,RooFit::PrintLevel(print_level));
      std::cout << " fit completed\n";
    }
    else std::cout << "... no fit!\n";
  }
  else {
    std::cout << "\"interactive\" fit to the data\n";
    RooAbsReal *nll=emodel.createNLL(data);
    RooMinuit m(*nll);
    m.migrad();
    m.hesse();
    m.minos();
  }
  emodel.plotOn(frame2);
  emodel.plotOn(frame2,RooFit::Components(RooArgSet(esig1)),RooFit::LineStyle(kDashed),RooFit::LineColor(kRed));
  emodel.plotOn(frame2,RooFit::Components(RooArgSet(esig2)),RooFit::LineStyle(kDotted),RooFit::LineColor(kGreen+2));

  char name[30];
  sprintf(name,"FitExplicit2_%d",cExtraCanvasIdx);
    TCanvas* c4 = new TCanvas(name,name,800,400) ;
    c4->Divide(2);
    c4->cd(1);
    gPad->SetLeftMargin(0.15) ; frame->GetYaxis()->SetTitleOffset(1.6) ; frame->Draw() ;
    c4->cd(2);
    gPad->SetLeftMargin(0.15) ; frame2->GetYaxis()->SetTitleOffset(1.6) ; frame2->Draw() ;

    c4->Update();
    sprintf(name,"fit%d.png",cExtraCanvasIdx);
    c4->SaveAs(name);
    //std::cout << "updated c4 : enter a char...\n"; char xxx; std::cin >> xxx;
  return 1;
}
*/

// ----------------------------------------------------------

/*
int FitTestExplicit1(const RooDataSet *dummy, RooRealVar &x, RooRealVar &y, const DataInfos_t &info) {
  std::cout << "FitTestExplicit1\n";
  double avg=120;
  double stddev=10.5;
  double alpha=0.005;
  double ampl=0.008;
  double center=avg*0.8;

  RooDataSet* hd=toyDSModel2(x,avg,stddev,ampl,alpha,center, 1000);
  std::cout << "calling to FitExplicit1\n";
  return FitExplicit1(hd,x,y,info);
}
*/
  
// ----------------------------------------------------------

/*
int FitTestExplicit2(const RooDataSet *dummy, RooRealVar &x, RooRealVar &y, const DataInfos_t &info) {
  std::cout << "FitTestExplicit2\n";
  double avg=120;
  double stddev=10.5;
  double alpha=0.005;
  double ampl=0.008;
  double center=avg*0.8;

  RooDataSet* hd=toyDSModel2(x,avg,stddev,ampl,alpha,center, 1000);
  std::cout << "calling to FitExplicit2\n";
  return FitExplicit2(hd,x,y,info);
}
*/  

// ----------------------------------------------------------
// ======================================================================
// ======================================================================
// ----------------------------------------------------------

/*
int FitExplicit2H(const RooDataHist *ds, RooRealVar &x, RooRealVar &y, const DataInfos_t &info) {
  std::cout << "entered FitExplicit2H. Simplified FitExplicit1. Two Gaussian fit\n";
  RooDataHist data(*ds);

  // S e t u p   p d f   a n d   l i k e l i h o o d 
  // -----------------------------------------------

  // Observable

  RooRealVar mean1("mean1","mean of g1",70, 0,140) ;
  RooRealVar sigma_g1("sigma_g1","width of g1",3.5,3.5,20) ;
  RooGaussian g1("g1","g1",x,mean1,sigma_g1) ;
  RooRealVar mean2("mean2","mean of g2",140, 50,200) ;
  RooRealVar sigma_g2("sigma_g2","width of g2",4,0.5,600.0) ;
  RooGaussian g2("g2","g2",x,mean2,sigma_g2) ;

  RooRealVar nsig1("nsig1","nsignal1",500,0,1e7);
  RooExtendPdf esig1("esig1","esig1",g1,nsig1);
  RooRealVar nsig2("nsig2","nsignal2",500,0,1e7);
  RooExtendPdf esig2("esig2","esig2",g2,nsig2);
  RooAddPdf emodel("emodel","emodel",RooArgList(esig1,esig2));

  if (1) {
    mean1.setVal(90); mean1.setRange(80,100);
    //sigma_g1.setRange(3,3);
    mean2.setVal(40); mean2.setRange(30,140); //mean2.setRange(0,500);    
    sigma_g2.setVal(50); //sigma_g2.setRange(50,50);
  }

  // Generate 1000 events
  //RooDataSet* data = model.generate(x,1000) ;
  
  RooPlot* frame = x.frame(RooFit::Title("starting frame") );
  data.plotOn(frame) ; 
  emodel.plotOn(frame);
  RooPlot *frame2 = x.frame(RooFit::Title("extended fit") );
  data.plotOn(frame2);

  if (0) {
    const int print_level=0; // -1 is minimal
    std::cout << "default call to fitTo\n";
    if (1) {
      emodel.fitTo(data,RooFit::Extended(true)); //,RooFit::PrintLevel(print_level));
      std::cout << " fit completed\n";
    }
    else std::cout << "... no fit!\n";
  }
  else {
    std::cout << "\"interactive\" fit to the data\n";
    RooAbsReal *nll=emodel.createNLL(data);
    RooMinuit m(*nll);
    m.migrad();
    m.hesse();
    m.minos();
  }
  emodel.plotOn(frame2);
  emodel.plotOn(frame2,RooFit::Components(RooArgSet(esig1)),RooFit::LineStyle(kDashed),RooFit::LineColor(kRed));
  emodel.plotOn(frame2,RooFit::Components(RooArgSet(esig2)),RooFit::LineStyle(kDotted),RooFit::LineColor(kGreen+2));

    TCanvas* c4 = new TCanvas("FitExplicit2","FitExplicit2",800,400) ;
    c4->Divide(2);
    c4->cd(1);
    gPad->SetLeftMargin(0.15) ; frame->GetYaxis()->SetTitleOffset(1.6) ; frame->Draw() ;
    c4->cd(2);
    gPad->SetLeftMargin(0.15) ; frame2->GetYaxis()->SetTitleOffset(1.6) ; frame2->Draw() ;

    c4->Update();
    c4->SaveAs("fit.png");
    if (!no_wait) { std::cout << "updated c4 : enter a char...\n"; char xxx; std::cin >> xxx; }
  return 1;
}
*/

// ----------------------------------------------------------

/*
int FitTestExplicit2H(const RooDataHist *dummy, RooRealVar &x, RooRealVar &y, const DataInfos_t &info) {
  std::cout << "FitTestExplicit2H\n";
  double avg=120;
  double stddev=10.5;
  double alpha=0.005;
  double ampl=0.008;
  double center=avg*0.8;
  int nbins=100;

  TH1F* hist=toyModel2(nbins,x.getMin(),x.getMax(),avg,stddev,ampl,alpha,center, 1000);
  RooDataHist* hd=ConvertToDataHist(hist);
  std::cout << "calling to FitExplicit2H\n";
  return FitExplicit2H(hd,x,y,info);
}
*/

// ----------------------------------------------------------

/*
int FitExplicit3H(const RooDataHist *ds, RooRealVar &x, RooRealVar &y, const DataInfos_t &info) {
  std::cout << "entered FitExplicit3H. One Gaussian and Chebyshev polynomial fit\n";
  RooDataHist data(*ds);

  // S e t u p   p d f   a n d   l i k e l i h o o d 
  // -----------------------------------------------

  // Observable

  RooRealVar mean1("mean1","mean of g1",70, 0,140) ;
  RooRealVar sigma_g1("sigma_g1","width of g1",3.5,3.5,20) ;
  RooGaussian g1("g1","g1",x,mean1,sigma_g1) ;
  //RooRealVar mean2("mean2","mean of g2",140, 50,200) ;
  //RooRealVar sigma_g2("sigma_g2","width of g2",4,0.5,600.0) ;
  //RooGaussian g2("g2","g2",x,mean2,sigma_g2) ;

  // Build Chebychev polynomial p.d.f.  
  RooRealVar a0("a0","a0",0.1,-2.,2.) ;
  RooRealVar a1("a1","a1",-0.1,-2.,2.) ;
  RooRealVar a2("a2","a2",-0.1,-2.,2.) ;
  RooChebychev bkg("bkg","Background",x,RooArgSet(a0,a1,a2)) ;
  //RooPolyVar bkg_poly("bkg_poly","Background(poly)",x,RooArgSet(a0,a1,a2));

  RooRealVar nsig1("nsig1","nsignal1",500,0,1e7);
  RooExtendPdf esig1("esig1","esig1",g1,nsig1);
  //RooRealVar nsig2("nsig2","nsignal2",500,0,1e7);
  //RooExtendPdf esig2("esig2","esig2",g2,nsig2);
  //RooAddPdf emodel("emodel","emodel",RooArgList(esig1,esig2));

  RooRealVar nbkg("nbkg","nbkg",500,0,1e7);
  RooExtendPdf ebkg("ebkg","ebkg",bkg,nbkg);
  RooAddPdf emodel("emodel","emodel",RooArgList(esig1,ebkg));

  if (1) {
    mean1.setVal(90); mean1.setRange(80,100);
    //sigma_g1.setRange(3,3);
    //mean2.setVal(40); mean2.setRange(30,140); //mean2.setRange(0,500);    
    //sigma_g2.setVal(50); //sigma_g2.setRange(50,50);
  }

  // Generate 1000 events
  //RooDataSet* data = model.generate(x,1000) ;
  
  RooPlot* frame = x.frame(RooFit::Title("starting frame") );
  data.plotOn(frame) ; 
  emodel.plotOn(frame);
  RooPlot *frame2 = x.frame(RooFit::Title("extended fit") );
  data.plotOn(frame2);

  if (0) {
    const int print_level=0; // -1 is minimal
    std::cout << "default call to fitTo\n";
    if (1) {
      emodel.fitTo(data,RooFit::Extended(true)); //,RooFit::PrintLevel(print_level));
      std::cout << " fit completed\n";
    }
    else std::cout << "... no fit!\n";
  }
  else {
    std::cout << "\"interactive\" fit to the data\n";
    RooAbsReal *nll=emodel.createNLL(data);
    RooMinuit m(*nll);
    m.migrad();
    m.hesse();
    m.minos();
  }
  emodel.plotOn(frame2);
  emodel.plotOn(frame2,RooFit::Components(RooArgSet(esig1)),RooFit::LineStyle(kDashed),RooFit::LineColor(kRed));
  emodel.plotOn(frame2,RooFit::Components(RooArgSet(ebkg)),RooFit::LineStyle(kDotted),RooFit::LineColor(kGreen+2));

    TCanvas* c4 = new TCanvas("FitExplicit3H","FitExplicit3H",800,400) ;
    c4->Divide(2);
    c4->cd(1);
    gPad->SetLeftMargin(0.15) ; frame->GetYaxis()->SetTitleOffset(1.6) ; frame->Draw() ;
    c4->cd(2);
    gPad->SetLeftMargin(0.15) ; frame2->GetYaxis()->SetTitleOffset(1.6) ; frame2->Draw() ;

    c4->Update();
    c4->SaveAs("fit.png");
    if (!no_wait) { std::cout << "updated c4 : enter a char...\n"; char xxx; std::cin >> xxx; }
  return 1;
}
*/

// ----------------------------------------------------------

/*
int FitTestExplicit3H(const RooDataHist *dummy, RooRealVar &x, RooRealVar &y, const DataInfos_t &info) {
  std::cout << "FitTestExplicit3H\n";
  double avg=120;
  double stddev=10.5;
  double alpha=0.005;
  double ampl=0.008;
  double center=avg*0.8;
  int nbins=100;

  TH1F* hist=toyModel2(nbins,x.getMin(),x.getMax(),avg,stddev,ampl,alpha,center, 1000);
  RooDataHist* hd=ConvertToDataHist(hist);
  std::cout << "calling to FitExplicit3H\n";
  return FitExplicit3H(hd,x,y,info);
}
*/

// ======================================================================
// ======================================================================

/*
TH1F* toyModel1(int nbins, double xmin, double xmax, double gauss_avg, double gauss_stddev, double distortion_amplitude, double distortion_spread, double distortion_center, int count) {
  double avg=gauss_avg;
  double stddev=gauss_stddev;
  double alpha=distortion_spread;
  double ampl=distortion_amplitude;
  double center=distortion_center;
  TH1F* hh = new TH1F("hh","hh",nbins,xmin,xmax);
  hh->SetDirectory(0);
  for (int i=0 ; i<count ; i++) {
    double r=gaussRandom(avg,stddev);
    double v= ampl*r*r*exp(-alpha*(r-center)*(r-center));
    hh->Fill(v+r) ;
  }
  //TH1* hd= hh;
  //return hd;
  return hh;
}
*/

// ======================================================================

/*
TH1F* toyModel2(int nbins, double xmin, double xmax, double gauss_avg, double gauss_stddev, double distortion_amplitude, double distortion_spread, double distortion_center, int count) {
  double avg=gauss_avg;
  double stddev=gauss_stddev;
  double alpha=distortion_spread;
  double ampl=distortion_amplitude;
  double center=distortion_center;
  TH1F* hh = new TH1F("hh","hh",nbins,xmin,xmax);
  hh->SetDirectory(0);
  for (int i=0 ; i<count ; i++) {
    double r=gaussRandom(avg,stddev);
    double v= ampl*r*r*exp(-alpha*(r-center)*(r-center));
    hh->Fill(v+r) ;
    hh->Fill(r);  // difference from toy model 1
  }
  //TH1* hd= hh;
  //return hd;
  return hh;
}
*/

// --------------------------------------------------------------------
// --------------------------------------------------------------------
/*
TH1* mytoyModel1(int nbins, double xmin, double xmax) {
  double avg=120.;
  double stddev=10.5;
  double alpha=0.005;
  double ampl=0.02;
  double center=avg*0.8;
  return toyModel1(nbins,xmin,xmax, avg,stddev, 0,alpha,center);
}
*/

// ======================================================================
// ======================================================================

// ======================================================================
/*
RooDataSet* toyDSModel1(RooRealVar &x, double gauss_avg, double gauss_stddev, double distortion_amplitude, double distortion_spread, double distortion_center, int count) {
  double avg=gauss_avg;
  double stddev=gauss_stddev;
  double alpha=distortion_spread;
  double ampl=distortion_amplitude;
  double center=distortion_center;
  
  RooDataSet *d = new RooDataSet("d","d",RooArgSet(x));
  RooRealVar z(x);
  RooArgSet value(z);

  for (int i=0 ; i<count ; i++) {
    double r=gaussRandom(avg,stddev);
    double v= ampl*r*r*exp(-alpha*(r-center)*(r-center));
    if ((r+v>=x.getMin()) && (r+v<=x.getMax())) {
      z.setVal(r+v); d->add(value);
    }
  }
  return d;
}
*/

// ======================================================================

/*
RooDataSet* toyDSModel2(RooRealVar &x, double gauss_avg, double gauss_stddev, double distortion_amplitude, double distortion_spread, double distortion_center, int count) {
  double avg=gauss_avg;
  double stddev=gauss_stddev;
  double alpha=distortion_spread;
  double ampl=distortion_amplitude;
  double center=distortion_center;
  
  RooDataSet *d = new RooDataSet("d","d",RooArgSet(x));
  RooRealVar z(x);
  RooArgSet value(z);

  for (int i=0 ; i<count ; i++) {
    double r=gaussRandom(avg,stddev);
    double v= ampl*r*r*exp(-alpha*(r-center)*(r-center));
    if ((r+v>=x.getMin()) && (r+v<=x.getMax())) {
      z.setVal(r+v); d->add(value);
    }
    if ((r>=x.getMin()) && (r<=x.getMax())) {  // addition to the toyDSModel1
      z.setVal(r);   d->add(value);
    }
  }
  return d;
}
*/

// --------------------------------------------------------------------
// --------------------------------------------------------------------

/*
RooDataSet* mytoyDSModel1(RooRealVar &x, int count) {
  double avg=120.;
  double stddev=10.5;
  double alpha=0.005;
  double ampl=0.01;
  double center=avg*0.8;
  return toyDSModel1(x, avg,stddev, 0,alpha,center, count);
}
*/

// ======================================================================
#ifndef __myLib__

int FitBreitWignerIlyaSpecial(const std::vector<RooDataSet*> &ds, RooRealVar &x, RooRealVar &y, const DataInfos_t &info) {
  std::cout << "entered FitBreitWignerIlyaSpecial. Uses extended fit\n";
  //std::cout << "First dataset is \"total\", the second is assumed to be \"pass\"\n";
  if (ds.size()!=2) {
    std::cout << "method is foreseen for 2 datasets only\n";
    return 0;
  }

  const unsigned int N=ds.size();
  std::vector<RooDataSet*> dataV;
  dataV.reserve(N);
  for (unsigned int i=0; i<ds.size(); ++i) {
    dataV.push_back(new RooDataSet(*ds[i]));
  }
  
  RooDataSet *ds_fail=dataV[0];
  RooDataSet *ds_pass=dataV[1];

  // Define category to distinguish physics and control sample events
  std::vector<TString> categories;
  categories.reserve(2);
  RooCategory sample("sample","sample");
  categories.push_back("pass");  sample.defineType("pass");
  categories.push_back("fail");  sample.defineType("fail");

  // a combined dataSet
  RooDataSet *combData=new RooDataSet("combData","combined data",x,RooFit::Index(sample),RooFit::Import("pass",*ds_pass),RooFit::Import("fail",*ds_fail));


  // Breit-Wigner
  const int allow_to_vary=1;
  const double bw_mean_0=91.188;
  const double bw_mean_min=(allow_to_vary) ? 80 : bw_mean_0;
  const double bw_mean_max=(allow_to_vary) ? 100 : bw_mean_0;
  const double bw_width_0=2.495;
  const double bw_width_min=(allow_to_vary) ? 1 : bw_width_0;
  const double bw_width_max=(allow_to_vary) ? 5 : bw_width_0;
  RooRealVar bw_mean("bw_mean","bw_mean",bw_mean_0,bw_mean_min,bw_mean_max);
  RooRealVar bw_width("bw_width","bw_width",bw_width_0,bw_width_min,bw_width_max);
  RooBreitWigner bw("bw","Breit Wigner", x, bw_mean, bw_width);
  
  // Crystal ball
  RooRealVar cb_m0("cb_m0","cb_m0",-3,3);
  RooRealVar cb_sigma("cb_sigma","cb_sigma",0,7);
  RooRealVar cb_alpha("cb_alpha","cb_alpha",0.1,2);
  RooRealVar cb_n("cb_n","cb_n",0.1,20);
  RooCBShape cb("CrystalBall","Crystal ball", x, cb_m0, cb_sigma, cb_alpha, cb_n);

  // convolute
  RooFFTConvPdf failing("failing","failing: bw x cb",x,bw,cb);
  
  // passing sample
  RooRealVar nevents("nevents","n_events",500,0,1e7);
  RooRealVar eff("eff","pass fail ratio",0.5,0,1);  
  RooFormulaVar passWeight("passWeight","pass weight","nevents*eff",RooArgList(nevents,eff));
  RooGenericPdf simpleSignal("simpleSignal","simpleSignal","1.0",RooArgList());
  RooExtendPdf ePass("ePass","ePass",simpleSignal,passWeight);

  // failing sample
  RooFormulaVar failWeight("failWeight","fail weight","nevents*(1-eff)",RooArgList(nevents,eff));
  RooExtendPdf eFail("eFail","eFail",failing,failWeight);

  // exponential background
  RooRealVar bkg_tau("bkg_tau","bkg_tau",-0.1,0.1);
  RooExponential bkg("background","expo decay",x,bkg_tau);
  RooRealVar bkgEvents("bkgEvents","bkgEvents",200,0,1e7);
  RooExtendPdf ebkg("ebkg","ebkg",bkg,bkgEvents);
  
  // Construct a simultaneous pdf
  RooAddPdf modelPass("modelPass","modelPass",RooArgList(ePass));
  RooAddPdf modelFail("modelFail","modelFail",RooArgList(eFail,ebkg));
  RooSimultaneous *simPdf=new RooSimultaneous("simPdf","simultaneous pdf",sample);
  simPdf->addPdf(modelPass,"pass");
  simPdf->addPdf(modelFail,"fail");


  // Generate 1000 events
  //RooDataSet* data = model.generate(x,1000) ;

  // make a plot
  //TString frame1Name="initial frame";
  RooPlot* frame1pass = x.frame(RooFit::Title("ini.frame pass") );
  RooPlot* frame1fail = x.frame(RooFit::Title("ini.frame fail") );
  combData->plotOn(frame1pass, RooFit::Cut("sample==sample::pass"), RooFit::LineColor(info.color(0))) ; 
  combData->plotOn(frame1fail, RooFit::Cut("sample==sample::fail"), RooFit::LineColor(info.color(1))) ; 
  simPdf->plotOn(frame1pass,RooFit::Slice(sample,"pass"),RooFit::ProjWData(sample,*combData),RooFit::LineColor(info.color(0)));
  simPdf->plotOn(frame1fail,RooFit::Slice(sample,"fail"),RooFit::ProjWData(sample,*combData),RooFit::LineColor(info.color(1)));
  //RooPlot *frame2 = NULL;
  TString frame2Name;

  if (!fit_dont_fit) simPdf->fitTo(*combData, RooFit::Extended(true), RooFit::NumCPU(2,true), RooFit::Timer(true));
  frame2Name = "implicit fit";

  if (info.NameOk()) {
    frame2Name.Append(" for ");
    frame2Name.Append(info.getName());
  }
  RooPlot *frame2pass = x.frame(RooFit::Title("fit.frame pass"));
  RooPlot *frame2fail = x.frame(RooFit::Title("fit.frame fail"));
  combData->plotOn(frame2pass, RooFit::Cut("sample==sample::pass"), RooFit::LineColor(info.color(0))) ; 
  combData->plotOn(frame2fail, RooFit::Cut("sample==sample::fail"), RooFit::LineColor(info.color(1))) ; 
  simPdf->plotOn(frame2pass,RooFit::Slice(sample,"pass"),RooFit::ProjWData(sample,*combData),RooFit::LineColor(info.color(0)));
  simPdf->plotOn(frame2fail,RooFit::Slice(sample,"fail"),RooFit::ProjWData(sample,*combData),RooFit::LineColor(info.color(1)));
    
   TString canvas_name="BWcb";
    TCanvas* c4 = new TCanvas("BWcb",canvas_name,800,800) ;
    c4->Divide(2,2);
    c4->cd(1); gPad->SetLeftMargin(0.15) ; frame1pass->GetYaxis()->SetTitleOffset(1.6) ; frame1pass->Draw() ;
    c4->cd(2); gPad->SetLeftMargin(0.15) ; frame2pass->GetYaxis()->SetTitleOffset(1.6) ; frame2pass->Draw() ;
    c4->cd(3); gPad->SetLeftMargin(0.15) ; frame1fail->GetYaxis()->SetTitleOffset(1.6) ; frame1fail->Draw() ;
    c4->cd(4); gPad->SetLeftMargin(0.15) ; frame2fail->GetYaxis()->SetTitleOffset(1.6) ; frame2fail->Draw() ;

    c4->Update();
    const char *fname="fit_frame100.png";
    c4->SaveAs(fname);
    std::cout << "saved canvas as " << fname << "\n";
    if (!no_wait) { std::cout << "updated c4 : enter a char...\n"; char xxx; std::cin >> xxx; }
  return 1;
}

// ======================================================================

int FitBreitWignerIlyaLiteral(const std::vector<RooDataSet*> &ds, RooRealVar &x, RooRealVar &y, const DataInfos_t &info) {
  std::cout << "entered FitBreitWignerIlyaLiteral. Uses extended fit\n";
  //std::cout << "First dataset is \"total\", the second is assumed to be \"pass\"\n";
  if (ds.size()!=2) {
    std::cout << "method is foreseen for 2 datasets only\n";
    return 0;
  }

  const unsigned int N=ds.size();
  std::vector<RooDataSet*> dataV;
  dataV.reserve(N);
  for (unsigned int i=0; i<ds.size(); ++i) {
    dataV.push_back(new RooDataSet(*ds[i]));
  }
  
  RooRealVar mass(x);
  RooDataSet *dataUnbinnedPass=dataV[1];
  RooDataSet *dataUnbinnedFail=dataV[0];
  RooDataHist *dataBinnedPass   = NULL;
  RooDataHist *dataBinnedFail   = NULL;

  RooCategory probeType("probeType","probeType");
  probeType.defineType("pass");
  probeType.defineType("fail");
  RooAbsData *data;
 
  // If needed do binned fit
  bool unbinnedFit = true;
  if(unbinnedFit){
    data = new RooDataSet("data","data",mass,RooFit::Index(probeType),
			  RooFit::Import("pass",*dataUnbinnedPass),
			  RooFit::Import("fail",*dataUnbinnedFail));
    cout << endl << "Setting up UNBINNED fit" << endl << endl;
  }else{
    dataBinnedPass   = dataUnbinnedPass->binnedClone("dataBinnedPass","dataBinnedPass");
    dataBinnedFail   = dataUnbinnedFail->binnedClone("dataBinnedFail","dataBinnedFail");
    data = new RooDataHist("data","data",mass,RooFit::Index(probeType),
			   RooFit::Import("pass",*dataBinnedPass),
			   RooFit::Import("fail",*dataBinnedFail));
    cout << endl << "Setting up BINNED fit" << endl << endl;
  }
  // Define the PDFs
  //
  // Common pieces
  //
  // True signal model
  RooRealVar zMass ("zMass" ,"zMass" ,91.188);
  RooRealVar zWidth("zWidth","zWidth",2.495);
  RooBreitWigner bwPdf("bwPdf","bwPdf",mass,zMass,zWidth);
  // Total signal events and efficiency
  RooRealVar nsignal("nsignal","nsignal",1000,0.0,1.0e7);
  RooRealVar eff    ("eff"    ,"eff"    ,0.7,0.0,1.0);
  //  PDF for the PASS sample
  RooFormulaVar nsigPass("nsigPass","nsigPass","@0*@1",RooArgList(nsignal,eff));
  
  //
  // Background
  RooRealVar lambdaBgPass("lambdaBgPass","lambdaBgPass",-0.1, -0.5, 0.0);
  RooExponential bgPassPdf("bgPassPdf","bgPassPdf",mass,lambdaBgPass);
  // Signal
  //     - resolution function
  RooRealVar cbMeanPass("cbMeanPass","cbMeanPass"   ,0.0, -10.0,10.0);
  RooRealVar cbWidthPass("cbWidthPass","cbWidthPass",1.0,  0.1, 5.0);
  RooRealVar cbAlphaPass("cbAlphaPass","cbAlphaPass",5.0,  0.0,20.0);
  RooRealVar cbNPass("cbNPass","cbNPass"            ,1.0,  0.0,10.0);
  RooCBShape cbPassPdf("cbPassPdf","cbPassPdf",mass,cbMeanPass,cbWidthPass,cbAlphaPass,cbNPass);
  //     - realistic model
  mass.setBins(10000,"cache");
  RooFFTConvPdf signalPassPdf("signalPassPdf","signalPassPdf",mass,bwPdf, cbPassPdf);
   // Combine signal and background
   RooRealVar nbgPass ("nbgPass" ,"nbgPass" ,1,0.0,1.0e5);
   RooAbsPdf *passPdf;
   const int COUNTnFIT=1;
   const int FITnFIT=2;
   int mode=(fit_Ilya_version==1) ?  1 : (fit_Ilya_version/10);
   if (mode==0) mode=COUNTnFIT;
   if( mode == COUNTnFIT ){
     std::cout << "COUNTnFIT branch\n";
     RooGenericPdf * simpleSignal = new RooGenericPdf("simpleSignal","simpleSignal", "1.0",RooArgList());
     RooExtendPdf * simpleSignalExtended = new RooExtendPdf("passPdf","passPdf", *simpleSignal, nsigPass);
     passPdf = simpleSignalExtended;
   }
   else if( mode == FITnFIT ){
     std::cout << "FITnFIT branch\n";
     passPdf = new RooAddPdf("passPdf","passPdf",RooArgList(signalPassPdf,bgPassPdf), RooArgList(nsigPass,nbgPass));
   }
   else{
     printf("ERROR: inappropriate mode requested\n");
     return 0;
   }
   //
   //  PDF for the FAIL sample
   //
   // Background
   RooRealVar lambdaBgFail("lambdaBgFail","lambdaBgFail",-0.1, -0.5, 0.0);
   RooExponential bgFailPdf("bgFailPdf","bgFailPdf",mass,lambdaBgFail);
   // Signal
   //     - resolution function
   RooRealVar cbMeanFail("cbMeanFail","cbMeanFail"   ,0.0, -10.0,10.0);
   RooRealVar cbWidthFail("cbWidthFail","cbWidthFail",1.0,   0.1, 5.0);
   RooRealVar cbAlphaFail("cbAlphaFail","cbAlphaFail",5.0,   0.0,20.0);
   RooRealVar cbNFail("cbNFail","cbNFail"            ,1.0,   0.0,10.0);
   RooCBShape  cbFailPdf ("cbFailPdf","cbFailPdf",mass,cbMeanFail,cbWidthFail,cbAlphaFail,cbNFail);
   //     - realistic model
   RooFFTConvPdf signalFailPdf("signalFailPdf","signalFailPdf",mass, bwPdf, cbFailPdf);
   // Combine signal and background
   RooFormulaVar nsigFail("nsigFail","nsigFail","@0*(1.0-@1)",RooArgList(nsignal,eff));
   RooRealVar nbgFail ("nbgFail" ,"nbgFail" ,1,0.0,1.0e5);
   RooAddPdf failPdf("failPdf","failPdf",RooArgList(signalFailPdf,bgFailPdf), RooArgList(nsigFail,nbgFail));
 
   // Combine pass and fail
   RooSimultaneous fullPdf("fullPdf","fullPdf",probeType);
   fullPdf.addPdf(*passPdf,"pass");
   fullPdf.addPdf(failPdf,"fail");
 
 
   // Do the fit
   // Start with a reasonable point and do rough approximation first
   double total = dataUnbinnedPass->numEntries() + dataUnbinnedFail->numEntries();
   nsignal.setVal( 0.99*total);
   eff.setVal(0.90);
   if( mode == FITnFIT ){
     nbgPass.setVal(0.01*total);
     cbAlphaPass.setVal(1.0);
     cbNPass    .setVal(5.0);
     cbAlphaPass.setConstant(kTRUE);
     cbNPass    .setConstant(kTRUE);
   }
   nbgFail.setVal(0.01*total);
   cbAlphaFail.setVal(0.5);
   cbNFail    .setVal(5.0);

   RooFitResult *result=NULL;
   if (!fit_dont_fit) {
   if (fit_Ilya_version%10) {
     std::cout << "fit_Ilya_version=1 (with intermediate fits)\n";
     cbAlphaFail.setConstant(kTRUE);
     cbNFail    .setConstant(kTRUE);
     result = fullPdf.fitTo(*data,
			    RooFit::Extended(kTRUE),
			    RooFit::Save());
     
     // Release shape parameters and refine the fit
     if( mode == FITnFIT ){
       cbAlphaPass.setConstant(kFALSE);
       cbNPass    .setConstant(kFALSE);
     }
     cbAlphaFail.setConstant(kFALSE);
     cbNFail    .setConstant(kFALSE);
     result = fullPdf.fitTo(*data,
			    RooFit::Extended(kTRUE),
			    RooFit::Minos(RooArgSet(eff)),
			    RooFit::Save());
     // If minos fails, refit without minos
     if((fabs(eff.getErrorLo())<5e-5) || (eff.getErrorHi()<5e-5))
       std::cout << "minos failed\n";
       result = fullPdf.fitTo(*data,
			      RooFit::Extended(kTRUE),
			      RooFit::Save());
   }
   else {
     std::cout << "fit_Ilya_version=0 (no intermediate fits)\n";
     cbAlphaFail.setConstant(kFALSE);
     cbNFail    .setConstant(kFALSE);
     cbAlphaPass.setConstant(kFALSE);
     cbNPass    .setConstant(kFALSE);
     result = fullPdf.fitTo(*data,
			    RooFit::Extended(kTRUE),
			    RooFit::Minos(RooArgSet(eff)),
			    RooFit::Save());
   }
   }

   double efficiency       = eff.getVal();
   double efficiencyErrHi  = eff.getErrorHi();
   double efficiencyErrLo  = fabs(eff.getErrorLo());
   std::cout << "eff= " << efficiency << ", errHi=" << efficiencyErrHi << ", errLo=" << efficiencyErrLo << "\n";
 
  if (dump_fit) {
    std::cout << "--------------------------------------------------\n";
    fullPdf.Print();
    std::cout << "Dumping fit parameters\n";
    std::cout << " 1 bw_mean1 " << zMass.getVal() << "\n";
    std::cout << " 2 bw_mean2 " << zMass.getVal() << "\n";
    std::cout << " 3 bw_width1 " << zWidth.getVal() << "\n";
    std::cout << " 4 bw_width2 " << zWidth.getVal() << "\n";
    std::cout << " 5 nSigEv    " << nsignal.getVal() << "\n";
    std::cout << " 6 eff_1     " << eff.getVal() << "\n";
    std::cout << " 7 bkg_tau_1 " << lambdaBgPass.getVal() << "\n";
    std::cout << " 8 bkg_tau_2 " << lambdaBgFail.getVal() << "\n";
    std::cout << " 9 cb_m0_1   " << cbMeanPass.getVal() << "\n";
    std::cout << "10 cb_sigma_1 " << cbWidthPass.getVal() << "\n";
    std::cout << "11 cb_alpha_1 " << cbAlphaPass.getVal() << "\n";
    std::cout << "12 cb_n_1     " << cbNPass.getVal() << "\n";
    std::cout << "13 nbkgr_1    " << nbgPass.getVal() << "\n";
    std::cout << "14 cb_m0_2    " << cbMeanFail.getVal() << "\n";
    std::cout << "15 cb_sigma_2 " << cbWidthFail.getVal() << "\n";
    std::cout << "16 cb_alpha_2 " << cbAlphaFail.getVal() << "\n";
    std::cout << "17 cb_n_2     " << cbNFail.getVal() << "\n";
    std::cout << "\n";
    std::cout << "18 nsignal_1  " << nsigPass.getVal() << "\n";
    std::cout << "19 nsignal_2  " << nsigFail.getVal() << "\n";
  }


   // Draw fit results
   TCanvas *c4=new TCanvas("c4lit","c4lit",800,400);
   c4->Divide(2);
   c4->cd(1);
   //passPad->cd();
   //passPad->Clear();
   RooPlot *framePass = mass.frame();
   dataUnbinnedPass->plotOn(framePass);
   //if(mode == FITnFIT){
     passPdf->plotOn(framePass);
     passPdf->plotOn(framePass,RooFit::Components("bgPassPdf"),LineStyle(kDashed));
     //}
   framePass->Draw();
   //passPad->Update();
 
   c4->cd(2);
   //failPad->cd();
   //failPad->Clear();
   RooPlot *frameFail = mass.frame();
   dataUnbinnedFail->plotOn(frameFail);
   failPdf.plotOn(frameFail);
   failPdf.plotOn(frameFail,RooFit::Components("bgFailPdf"),LineStyle(kDashed));
   frameFail->Draw();
   //failPad->Update();
   if (no_plot) c4->Update();
 
   // Print fit outcome into fit log
   //result->printStream(fitLog,RooPrintable::kValue,RooPrintable::kVerbose);
   //fitLog << endl;
   //printCorrelations(fitLog, result);
   const char *fname_base="fit_frame101";
   std::string fname;
   if (fit_Ilya_version) fname="IK";
   if (!unbinnedFit) fname+="binned";
   if (mode==COUNTnFIT) fname+="-cnf.png";
   else fname+="-fnf.png";
   fname=std::string(fname_base)+fname;
   c4->SaveAs(fname.c_str());
   std::cout << "saved canvas as " << fname.c_str() << "\n";
    //c4->SaveAs("fit_frames.C");
    if (!no_wait) { std::cout << "updated c4 : enter a char...\n"; char xxx; std::cin >> xxx; }
 
  return 1;
}

// ======================================================================


int FitPFStemMCGauss(const std::vector<RooDataSet*> &ds, RooRealVar &x, RooRealVar &y, const DataInfos_t &info) {
  const int N2=ds.size();
  std::cout << "* entered FitBreitWignerVPFStemMCGauss(pass-fail).\n";
  std::cout << " * First entries describe the shape to be smeared by a Gaussian.\n";
  std::cout << " * Pass(Exp),Fail(Exp),Pass(MC),Fail(MC)\n";
  std::cout << "* Signal events are bound by 1 number\n";
  std::cout << "* Signal is (MC shape)@Gaussian,  background is an exponential\n";
  if (N2==0) { std::cout<<"no data provided. exciting\n"; return 0; }
  //if (!fit_pass_fail) { std::cout << "the flag -fit-pass-fail has to be on\n"; return 0; }
  if (N2%2!=0) { std::cout <<"to fit pass-fail an even number of datasets is needed. " << ds.size() << " were provided\n"; return 0; }
  const int N=ds.size()/2;

  // Monte-Carlo shape
  RooRealVar mass(x);
  mass.setBins(10000,"cache");
  mass.setBins(60);
  // following the example of rf901_numintconfig.C
  RooAbsReal::defaultIntegratorConfig()->getConfigSection("RooAdaptiveGaussKronrodIntegrator1D").setRealValue("maxSeg",1000) ;

  // Data signal
  std::vector<RooDataSet*> dataMCV,dataExpV;
  dataMCV.reserve(N); dataExpV.reserve(N);
  // separate and convert to fail-pass sequence
  for (unsigned int i=0; i<N; ++i) {
    int idx=N-1-int(i);
    dataExpV.push_back(new RooDataSet(*ds[idx]));
  }
  for (unsigned int i=N; i<2*N; ++i) {
    int idx=2*N-1-int(i);
    dataMCV.push_back(new RooDataSet(*ds[N+idx]));
  }
  std::cout << "sizes were "; for (unsigned int i=0; i<ds.size(); ++i) std::cout << " " << ds[i]->numEntries(); std::cout << "\n";
  std::cout << " MC counts "; for (unsigned int i=0; i<dataMCV.size(); ++i) std::cout << " " << dataMCV[i]->numEntries(); std::cout << "\n";
  std::cout << " Exp counts "; for (unsigned int i=0; i<dataExpV.size(); ++i) std::cout << " " << dataExpV[i]->numEntries(); std::cout << "\n";

  // create templatess
  std::vector<RooDataHist*> mcHistV;
  std::vector<RooHistPdf*> mcV;
  CreateRooHistPdf(mass,dataMCV,"MCHist_X",mcHistV,"MCShape_X",mcV);

  std::cout << "mcV created" << std::endl;

  // Gaussians
  std::vector<RooRealVar*> g_biasV; 
  std::vector<RooRealVar*> g_widthV;
  CreateRooRealVar(N,"g_biasX",g_biasV,0.2);
  CreateRooRealVar(N,"g_widthX",g_widthV,1,0.1,5);
  std::vector<RooGaussModel*> gaussV; 
  CreateRooGaussian(N,"gaussX",mass,g_biasV,g_widthV,gaussV);

  // convolute
  //std::vector<RooNumConvPdf*> signalV;
  //CreateConvPdf(N,"signalX",mass,gaussV,mcV,signalV);
  std::vector<RooFFTConvPdf*> signalV;
  mass.setBins(10000,"cache");
  CreateRooConvPdf(N,"signalX",mass,mcV,gaussV,signalV);
  
  // exponential background
  std::vector<RooRealVar*> bkg_tauV; //("bkg_tau","bkg_tau",-0.1,0.1);
  std::vector<RooExponential*> bkgV; //("background","expo decay",x,bkg_tau);
  CreateRooRealVar(N,"bkg_tau_X",bkg_tauV,-0.3,0.3);
  CreateRooExponential(N,"bkground_X",mass,bkg_tauV,bkgV);


  std::vector<RooRealVar*> eff;
  CreateRooRealVar(N-1,"eff_X",eff,0,1);
  int total_event_count=0;
  SetInitialPassFailEffValues(N,dataExpV,eff,total_event_count);
  RooArgSet calcAsymmErr;
  for (unsigned int i=0; i<eff.size(); ++i) {
    calcAsymmErr.add(*eff[i]);
  }
  
  RooRealVar nSigEv("nSigEv","nSigEv",total_event_count,0,1e9);
  // eff[0]  - related to the number of events that have passed [0]
  // (1-eff[0]) - related to the number of events that have failed [0]
  // the provided pass-fail data is (fail, (fail, (fail,.., (fail_last,passall))))

  std::vector<RooFormulaVar*> passf;
  CreatePassFailStemFormula(N,eff,nSigEv,passf);
  if (1) {
    std::cout << "initial passf values:\n";
    for (UInt_t i=0; i<passf.size(); ++i) {
      //std::cout << " #" << i << " " << passf.getValue() << "\n";
      std::cout << " #" << i << "  "; passf[i]->Print(); std::cout << "\n";
    }
  }


  // Extend variables
  std::vector<RooExtendPdf*> esigV; //("esig","esig",signal,nsignal);
  std::vector<RooRealVar*> nbkgV; //("nbkgr","nbkgr",500,0,1e7);
  std::vector<RooExtendPdf*> ebkgV; //("ebkg","ebkg",bkg,nbkg);
  CreateRooExtendPdf(N,"esig_X",signalV,passf,esigV);
  CreateRooRealVar(N,"nbkgr_X",nbkgV,500,0,1e8);
  CreateRooExtendPdf(N,"ebkgr_X",bkgV,nbkgV,ebkgV);


  // Combine the models
  std::vector<RooAddPdf*> model;
  CreateRooAddPdf(N,"modelX",esigV,ebkgV,model);

  // Define category to distinguish physics and control sample events
  TString name1,name2;
  size_t pos1,pos2;
  const char idxV[]={'1','2','3','4','5','6','7','8','9','A'};
  std::vector<TString> categories;
  categories.reserve(N);
  RooCategory sample("sample","sample");
  name1="failX"; pos1=name1.Length()-1;
  for (unsigned int i=0; i<N; ++i) {
    const char c=idxV[i]; name1[pos1]=c;
    if (i==N-1) { name1="passX"; name1[pos1]=idxV[i-1]; }
    categories.push_back(name1);
    sample.defineType(name1);
  }
  std::cout << "categories: "; for (unsigned int i=0; i<N; ++i) std::cout << " " << categories[i]; std::cout << "\n";


  // a Combined data set 
  RooDataSet *combData=CreateCombinedData("combData","combined data", dataExpV, mass,sample,categories);
  // Construct a simultaneous pdf
  RooSimultaneous *simPdf=CreateCombinedModel("simPdf","simultaneous pdf",model,sample,categories);

  // fitting

  TString frame1Name="initial frame";
  if (info.NameOk()) {
    frame1Name.Append(" for ");
    frame1Name.Append(info.getName());
  }

  vector<RooPlot*> frames;
  RooPlot* frame = x.frame(RooFit::Title(frame1Name) );
  // Plot all data according to the category
  for (unsigned int i=0; i<N; i++) {
    name1="sample==sample::" + categories[i];
    std::cout << "plot combData with name=" << name1 << "\n";
    combData->plotOn(frame,RooFit::Cut(name1), RooFit::LineColor(info.color(i)),RooFit::MarkerColor(info.color(i)));
  }
  //dataV[0]->plotOn(frame,RooFit::LineColor(kGreen+2),RooFit::LineStyle(kDashed));
  //if (N>1) dataV[1]->plotOn(frame,RooFit::LineColor(kRed),RooFit::LineStyle(kDashed));

  // plot initial trials
  for (unsigned int i=0; i<N; i++) {
    simPdf->plotOn(frame,RooFit::Slice(sample,categories[i]),RooFit::ProjWData(sample,*combData), RooFit::LineColor(info.color(i)));
  }
  frames.push_back(frame);
  RooPlot *frame2 = NULL;
  TString frame2Name;

  if (fit_explicit) std::cout << "the request to fit explicitly is ignored!\n";

  RooFitResult *fitRes=NULL;
  if (!fit_dont_fit) {
    //const int MIGRAD_max_iters=20;
    //std::cout << " Setting max iterations in MIGRAD to " << MIGRAD_max_iters << "\n";
    //TVirtualFitter::SetMaxIterations(MIGRAD_max_iters);
    for (unsigned int iter=0; iter<1; ++iter) { // no improvement due to more iters
      if (1) {
	fitRes=simPdf->fitTo(*combData,RooFit::Extended(true),RooFit::NumCPU(2,true), RooFit::Timer(true),RooFit::Save(),RooFit::Minos(calcAsymmErr));
      }
      else {
	std::cout << "\n ** Warning ** skewed error\n";
	fitRes=simPdf->fitTo(*combData,RooFit::Extended(true),RooFit::NumCPU(2,true), RooFit::Timer(true),RooFit::Save());
	fitRes=simPdf->fitTo(*combData,RooFit::Extended(true),RooFit::NumCPU(2,true), RooFit::Timer(true),RooFit::Save(),RooFit::Minos(calcAsymmErr));
      }
    }
  }
  frame2Name = "implicit fit";

  if (info.NameOk()) {
    frame2Name.Append(" for ");
    frame2Name.Append(info.getName());
  }
  frame2 = x.frame(RooFit::Title(frame2Name));
    
  // Plot all data according to the category
  if (0) {
    for (unsigned int i=0; i<N; i++) {
      name1="sample==sample::" + categories[i];
      combData->plotOn(frame2,RooFit::Cut(name1), RooFit::LineColor(info.color(i)));
    }
    // plot initial trials
    for (unsigned int i=0; i<N; i++) {
      simPdf->plotOn(frame2,RooFit::Slice(sample,categories[i]),RooFit::ProjWData(sample,*combData), RooFit::LineColor(info.color(i)),RooFit::MarkerColor(info.color(i)));
    }
  }
  else {
    for (unsigned int i=0; i<N; i++) {
      frame2Name=categories[i];
      RooPlot *framex= x.frame(RooFit::Title(frame2Name));
      name1="sample==sample::" + categories[i];
      std::cout << "\n plot " << name1 << "\n";
      combData->plotOn(framex,RooFit::Cut(name1), RooFit::LineColor(info.color(i)));
      //combData->plotOn(framex,RooFit::Cut(name1), RooFit::LineColor(info.color(i)),RooFit::MarkerColor(info.color(i)));
      simPdf->plotOn(framex,RooFit::Slice(sample,categories[i]),RooFit::ProjWData(sample,*combData), RooFit::LineColor(info.color(i)));
      //ebkgV[i]->plotOn(framex,RooFit::LineColor(kGreen+2), RooFit::LineStyle(kDashed));
      //for (unsigned int mi=0; (mi<i) && (mi<i); ++mi) {
      unsigned int mi=i;
	name1="esig_X"; pos1=name1.Length()-1;
	name2="ebkgr_X"; pos2=name2.Length()-1;
	const char c=idxV[i]; name1[pos1]=c; name2[pos2]=c;
	model[i]->plotOn(framex,RooFit::Components(name1.Data()),RooFit::LineStyle(kDashed),RooFit::Normalization(1.,RooAbsReal::RelativeExpected));
	model[i]->plotOn(framex,RooFit::Components(name2.Data()),RooFit::LineStyle(kDotted),RooFit::Normalization(1.,RooAbsReal::RelativeExpected));
	//}
      frames.push_back(framex);
    }
  }

   //model.plotOn(frame2,RooFit::Components(esig),RooFit::LineColor(kGreen+2),RooFit::LineStyle(kDashed));
   //model.plotOn(frame2,RooFit::Components(ebkg),RooFit::LineColor(kRed+2),RooFit::LineStyle(kDotted));
   //esig.plotOn(frame2,RooFit::LineColor(kGreen+2));
   //ebkg.plotOn(frame2,RooFit::LineColor(kRed+2));

   TString canvas_name="BWcb";
   TCanvas* c4 = NULL;
   if (0) {
     c4=new TCanvas("BWcb",canvas_name,800,400) ;
     c4->Divide(2);
     c4->cd(1);
     gPad->SetLeftMargin(0.15) ; frame->GetYaxis()->SetTitleOffset(1.6) ; frame->Draw() ;
     c4->cd(2);
     gPad->SetLeftMargin(0.15) ; frame2->GetYaxis()->SetTitleOffset(1.6) ; frame2->Draw() ;
   }
   else {
     int dx=(frames.size()>3) ? 2:frames.size();
     int dy=(frames.size()>3) ? 2:1;
     c4=new TCanvas("BWcb",canvas_name,400*dx,400*dy);
     c4->Divide(dx,dy);
     for (unsigned int i=0; i<frames.size(); ++i) {
       c4->cd(i+1);
       gPad->SetLeftMargin(0.15) ; frames[i]->GetYaxis()->SetTitleOffset(1.6) ; frames[i]->Draw() ;
     }
   }
   const char *outfile_name=(fit_rebinned) ? "fit_frame800R.png" : "fit_frame800.png";

    c4->Update();
    c4->SaveAs(outfile_name);
    //c4->SaveAs("fit_frames.C");

  if (dump_fit) {
    std::cout << "--------------------------------------------------\n";
    //simPdf->Print("t");
    simPdf->Print();
    std::cout << "Dumping fit parameters\n";
    for (unsigned int i=0; i<model.size(); ++i) {
      std::cout << "model #" << i+1 << " : ";// << model[i]->getVal() << "\n";
      model[i]->Print(); //std::cout << "\n";
      //std::cout << "error: " << model[i]->getError() << "\n";
    }
    for (unsigned int i=0; i<esigV.size(); ++i) {
      std::cout << "esigV #" << i+1 << " : ";
      esigV[i]->Print(); //std::cout << "\n";
    }
    for (unsigned int i=0; i<passf.size(); ++i) {
      std::cout << "passeff(passf) #" << i+1 << " : ";
      passf[i]->Print(); //std::cout << "\n";
    }
    std::cout << "derived signal event counts:\n";
    for (unsigned int i=0; i<passf.size(); ++i) {
      std::cout << " " << (i+1) << " nsignal_" << (i+1) << "   " << passf[i]->getVal() << " ";
      if (fitRes) std::cout << passf[i]->getPropagatedError(*fitRes);
      std::cout << "\n";
    }
    for (unsigned int i=0; i<ebkgV.size(); ++i) {
      std::cout << "ebkgV #" << i+1 << " : ";
      ebkgV[i]->Print(); //std::cout << "\n";
    }
    for (unsigned int i=0; i<nbkgV.size(); ++i) {
      std::cout << "nbkgV #" << i+1 << " : ";
      nbkgV[i]->Print(); //std::cout << "\n";
    }
    if (0) {
      std::cout << "another dump\n";
      for (unsigned int i=0; i<model.size(); ++i) {
	std::cout << "model #" << i+1 << " : ";// << model[i]->getVal() << "\n";
	model[i]->printCompactTree(); //std::cout << "\n";
	//std::cout << "error: " << model[i]->getError() << "\n";
      }
    }
    std::cout << "Efficiency errors :\n";
    for (unsigned int i=0; i<eff.size(); ++i) {
      std::cout << eff[i]->GetName() << " = " << eff[i]->getVal() << " +" << eff[i]->getErrorHi() << " " << eff[i]->getErrorLo() << "\n";
    }
    std::cout << "Event number errors :\n";
    for (unsigned int i=0; i<passf.size(); ++i) {
      std::cout << passf[i]->GetName() << " = " << passf[i]->getVal() << "  +/-" << ((fitRes) ? passf[i]->getPropagatedError(*fitRes) : 0.0) << "\n";
    }
    std::cout << "--------------------------------------------------\n";
  }

  std::cout << "plot saved as " << outfile_name << "\n";
    if (!no_wait) { std::cout << "updated c4 : enter a char...\n"; char xxx; std::cin >> xxx; }
  return 1;
}

#endif // __myLib__

// ======================================================================


// ---------------------------------------------------------

int makeHTML(const TString outDir, const TString html_fname, const TString description, const std::vector<std::string> &add_html_lines, const std::vector<TString> &nameEndings, const std::vector<TString> &nameEndingsForFile, RooRealVar &x, const std::vector<RooDataSet*> &combData, const std::vector<RooAddPdf*> &simPdf, const std::vector<std::vector<int>*>* combine, const std::vector<TString>* combination_names) {
  std::cout << "Entered makeHTML\n";
  // new version
  std::vector<TString> descriptionV;
  descriptionV.push_back(description);
  std::vector<const std::vector<RooDataSet*>*> combDataV;
  std::vector<const std::vector<RooAddPdf*>*> simPdfV;
  std::vector<RooRealVar*> xV;
  combDataV.push_back(&combData);
  simPdfV.push_back(&simPdf);
  xV.push_back(&x);
  return makeHTML(outDir,html_fname,descriptionV,add_html_lines,nameEndings,nameEndingsForFile,xV,combDataV,simPdfV,combine,combination_names);
}

// ---------------------------------------------------------

int makeHTML(const TString outDir, const TString html_fname, const std::vector<TString> &descriptionV, const std::vector<std::string> &add_html_lines, const std::vector<TString> &nameEndings, const std::vector<TString> &nameEndingsForFile, std::vector<RooRealVar*> &x, const std::vector<const std::vector<RooDataSet*>*> &combDataV, const std::vector<const std::vector<RooAddPdf*>*> &simPdfV, const std::vector<std::vector<int>*>* combine, const std::vector<TString>* combination_names) {
  std::cout << "Entered makeHTML(comparison)" << std::endl;
  const char *extension=".png";
  if (no_plot) std::cout << "WARNING: no_plot=" << no_plot << "\n";
  const unsigned int N=combDataV.size();
  if ((N!=simPdfV.size()) || (N!=descriptionV.size()) || (N!=x.size())) {
    std::cout << " combDataV.size=" << N << ", simPdfV.size=" << simPdfV.size() << ", descriptionV.size=" << descriptionV.size() << ", x.size=" << x.size() << ". bailing out\n";
    return 0;
  }
  for (unsigned int i=0; i<N; ++i) {
    if (combDataV[i]->size() != simPdfV[i]->size()) {
      std::cout << " combDataV[" << i << "]->size=" << combDataV[i]->size() << " versus simPdfV[" << i << "]->size=" << simPdfV[i]->size() << "\n";
      return 0;
    }
  }

  char buf[100];
  std::vector<TString> figfnames;
  std::vector<TText*> chi2Values; // only a container
  if (saveFigFileTag.size()) sprintf(buf,"cx_%s",saveFigFileTag.c_str()); else sprintf(buf,"cx");
  TCanvas *cx=(no_plot) ? NULL : new TCanvas(buf,buf,600,600);
  //std::vector<int> insertDiv;

  if (combine && combination_names) {
    if (combine->size() != combination_names->size()) {
      std::cout << "makeHTML got " << combine->size() << " combination instructions and " << combination_names->size() << " combination names\n";
      return 0;
    }
    chi2Values.reserve(simPdfV.size()*combine->size());
    for (unsigned int i=0; i<combine->size(); ++i) {
      TString name=(*combination_names)[i];
      TString figNameBase=(saveFigFileTag.size()) ? (saveFigFileTag+name) : name;
      for (int t=0; t<figNameBase.Length(); t++) {
	if ((figNameBase[t]==' ') ||
	    (figNameBase[t]=='-')) figNameBase[t]='_';
      }
      figNameBase.Append(extension);
      
      const std::vector<int>* comb=(*combine)[i];
      std::cout << "combine : "; for (unsigned int ii=0; ii<comb->size(); ++ii) std::cout << " " << (*comb)[ii]; std::cout << "\n";
      //insertDiv.push_back(0);
      //if ((i>0) && (comb->size()!=(*combine)[i-1]->size())) *insertDiv.back()=1;

      for (unsigned int part=0; part<N; ++part) {
	const std::vector<RooDataSet*> *combData= combDataV[part];
	const std::vector<RooAddPdf*> *simPdf= simPdfV[part];
	TString figName="partX"; 
	figName[figName.Length()-1]=char(int('1')+part);
	if (N>1) figName.Append(figNameBase); else figName=figNameBase;
	figfnames.push_back(figName);
	RooPlot* frame = x[part]->frame(RooFit::Title(name));


	std::vector<double> nEvtRel;
	nEvtRel.reserve((*comb).size());

	RooDataSet *rawData0=(*combData)[(*comb)[0]];
	RooDataHist *data=rawData0->binnedClone("data");
	RooDataSet *collectData=(RooDataSet*)rawData0->Clone("collectData");
	//data->plotOn(frame,RooFit::LineColor(kGreen),RooFit::MarkerColor(kGreen));
	nEvtRel.push_back(rawData0->numEntries());
	for (unsigned int ci=1; ci<comb->size(); ++ci) {
	  RooDataSet *rawData=(*combData)[(*comb)[ci]];
	  nEvtRel.push_back(rawData->numEntries());
	  RooDataHist *dt=rawData->binnedClone();
	  collectData->append(*rawData);
	  //dt->plotOn(frame,RooFit::LineColor(kGreen+2*ci),RooFit::MarkerColor(kGreen+2*ci));
	  //data->add(*combData[(*comb)[ci]]);
	  data->add(*dt);
	}
	data->plotOn(frame,RooFit::Name("data"),RooFit::LineColor(kBlack),RooFit::MarkerColor(kBlack));
	/*
	  if (0) {
	  // not working
	  RooAbsReal* chi=model.createChi2(*data,RooFit::Extended(kTRUE),RooFit::Range(60,120),RooFit::DataError(RooAbsData::Poisson));
	  sprintf(buf,"chi^2=%6.2e",chi->getVal());
	  TText *txt=new TText(0.2,0.8,buf);
	  chi2Values.push_back(txt);
	  txt->SetNDC(kTRUE); txt->SetTextSize(0.04); txt->SetTextColor(kBlue+2);
	  frame->addObject(txt);
	  }
	*/
	//! delete data; // postponed


	double nEvtTot=0;
	for (unsigned int ci=0; ci<nEvtRel.size(); ++ci) nEvtTot+=nEvtRel[ci];
	if (nEvtTot==double(0.)) {
	  std::cout << "Error: nEvtTot=0!!\n"; return 0;
	}
	std::vector<RooRealVar*> weights;
	CreateRooRealVar(comb->size(),"weights_X",weights,1.);
	for (unsigned int ci=0; ci<weights.size(); ++ci) {
	  weights[ci]->setVal(nEvtRel[ci]/nEvtTot);
	}
	std::vector<RooAddPdf*> includedModel;
	includedModel.reserve(comb->size());
	std::cout << "combination name=" << name << "\n";
	for (unsigned int ci=0; ci<comb->size(); ++ci) {
	  //std::cout << "including model " << (*simPdf)[(*comb)[ci]]->GetName() << " with norm=" << (*simPdf)[(*comb)[ci]]->getNorm() << "\n";
	  includedModel.push_back((*simPdf)[(*comb)[ci]]);
	}
	std::vector<RooExtendPdf*> weightedModel;
	CreateRooExtendPdf(comb->size(),"weightedModel_X",includedModel,weights, weightedModel);
      
	RooArgList args;
	for (unsigned int ci=0; ci<weightedModel.size(); ++ci) {
	  args.add(*weightedModel[ci]);
	}
	RooAddPdf model("model","model",args);

	//model.plotOn(frame,RooFit::LineColor(kRed+2),RooFit::Normalization(1,RooAbsReal::RelativeExpected));
	model.plotOn(frame,RooFit::LineColor(kRed+2));
	//model.plotOn(frame,RooFit::LineColor(kRed+4),RooFit::LineStyle(kDashed),RooFit::Normalization(1,RooAbsReal::RelativeExpected));
	//model.plotOn(frame,RooFit::LineColor(kGreen+4),RooFit::LineStyle(kDashed),RooFit::Normalization(RooAbsReal::RelativeExpected));
	//model.forceNumInt(kTRUE);
	//model.plotOn(frame,RooFit::LineColor(kGreen+2),RooFit::Normalization(1,RooAbsReal::RelativeExpected));
	//model.plotOn(frame,RooFit::LineColor(kGreen+4),RooFit::LineStyle(kDashed),RooFit::Normalization(1,RooAbsReal::RelativeExpected));
	chi2Values.push_back(CreateChiSquareText(frame)); // add info about chi2
	//std::cout << "last step reached" << std::endl;
	if (!no_plot) {
	  cx->cd();
	  gPad->SetLeftMargin(0.15);
	  frame->GetYaxis()->SetTitleOffset(1.6);
	  frame->Draw();	  
	  TString figSaveName=outDir+figName;
	  cx->SaveAs(figSaveName);
	  std::cout << "model norm = " << model.getNorm() << std::endl;
	  if (saveCFormat) {
	    sprintf(buf,"%s",figSaveName.Data());
	    char *p=strstr(buf,extension);
	    if (p) {
	      sprintf(p,".C");
	      cx->SaveAs(buf);
	      RemoveWrongStringFromFile(buf,"[x]");
	    }
	    if (1) {
	      if (saveFigFileTag.size()) sprintf(buf,"ct_%s",saveFigFileTag.c_str()); else sprintf(buf,"ct");
	      TCanvas *ct=new TCanvas(buf,buf,600,600);
	      double *DYMassBins=DYTools::getMassBins();
	      int bincount_tmp=DYTools::NumberOfMassBins();
	      bincount_tmp=60;
	      TH1D histoBase("histoBase","histoBase",bincount_tmp,fit_range_min,fit_range_max);
	      histoBase.GetXaxis()->SetTitle("mass");
	      histoBase.GetYaxis()->SetTitle("events");
	      TH1D *histoExp=ConvertToHisto(*collectData,"histoExp",&histoBase);
	      RooCmdArg args1=RooFit::Binning(bincount_tmp,fit_range_min,fit_range_max);
	      TH1D *histoMC=CreateHisto(model,"histoMC",x[part],args1,&histoBase);

	      //TH1D* histoMC=(TH1D*)model.createHistogram("histoMC",*x[part],);
	      histoMC->SetLineColor(kGreen+1);
	      //TH1D *histoExp=new TH1D("histoExp","histoExp",bincount_tmp,fit_range_min,fit_range_max);
	      histoExp->SetLineColor(kRed+2);
	      //if (collectData && (collectData->numEntries()>0)) {
	      //std::cout << "collecting Data\n";
	      //for (int ii=0; ii<collectData->numEntries(); ii++) {
	      //  const RooArgSet *arg=collectData->get(ii);
	      //  histoExp->Fill(arg->getRealValue("x",0,kTRUE));
		//}
		//histoMC->Scale(->Integral();
		//histoExp->Scale(1/collectData->numEntries());
	      //}
	      //if (histoExp->Integral()!=0) histoExp->Scale(1/histoExp->Integral());
	      if (histoMC->Integral()) histoMC->Scale(histoExp->Integral()/histoMC->Integral());
	      double ymax1=histoMC->GetMaximum();
	      double ymax2=histoExp->GetMaximum();
	      double ymax=(ymax1>ymax2) ? ymax1:ymax2;
	      histoBase.SetMaximum(1.05*ymax);
	      histoBase.Draw();
	      histoMC->Draw("AL* same");
	      histoExp->Draw("histo same");
	      std::cout << "chi2=" << histoMC->Chi2Test(histoExp,"WW CHI2") << "\n"; //std::cout << " chi2/ndof=" << histoMC->Chi2Test(histoExp,"WW CHI2/NDOF") << "\n";
	      TLatex *chi2Val=CreateChiSquareText(histoMC->Chi2Test(histoExp,"WW CHI2"),0);
	      //TLatex *chi2ValNdof=CreateChiSquareText(histoMC->Chi2Test(histoExp,"WW CHI2/NDOF"),1);
	      chi2Val->Draw();
	      //chi2ValNdof->SetY(0.75);
	      //chi2ValNdof->Draw();
	      
	      p=buf+(strlen(buf)-1);
	      ct->Update();
	      sprintf(p,"_test_histo.C");
	      ct->SaveAs(buf);
	      sprintf(p,"_test_histo.png");
	      ct->SaveAs(buf);
	      delete ct;
	      delete DYMassBins;
	      delete chi2Val; 
	      //delete chi2ValNdof;
	    }
	  }
	}
	if (data) delete data;
	if (collectData) delete collectData;
      }
    }
    std::cout << "combined plots created\n";
  }

  TPaveText noDataPointsTxt(0.2,0.45,0.8,0.6,"NDC br");
  noDataPointsTxt.SetFillStyle(0);
  noDataPointsTxt.AddText("no data points");

  for (unsigned int part=0; part<N; ++part) {
    const std::vector<RooDataSet*> *combData= combDataV[part];
    const std::vector<RooAddPdf*> *simPdf= simPdfV[part];    
    for (unsigned int i=0; i<combData->size(); ++i) {
      TString extra=saveFigFileTag;
      if (N>1) { 
	extra.Append("partX "); extra[extra.Length()-2]=char('1'+part); 
      }
      TString name=extra+"set " + nameEndings[i];
      if (extra.Length()>0) extra[extra.Length()-1]='_';
      TString figName=extra+"set"+nameEndingsForFile[i]+extension;
      figfnames.push_back(figName);
      RooPlot* frame = x[part]->frame(RooFit::Title(name));
      
      if ((*combData)[i]->numEntries()==0) {
	cx->cd();
	noDataPointsTxt.Draw();
      }
      else {
	(*combData)[i] ->plotOn(frame,RooFit::LineColor(kBlack),RooFit::MarkerColor(kBlack));
	(*simPdf)[i]->plotOn(frame,RooFit::LineColor(kRed+2)); //,RooFit::MarkerColor(kRed+2));
      }
      chi2Values.push_back(CreateChiSquareText(frame,kBlack)); // add info about chi2
      if (!no_plot) {
	cx->cd();
	gPad->SetLeftMargin(0.15);
	frame->GetYaxis()->SetTitleOffset(1.6);
	frame->Draw();
	TString figSaveName=outDir+figName;
	std::cout << "HERE figSaveName=" << figSaveName << std::endl;
	cx->SaveAs(figSaveName);
	if (saveCFormat) {
	  sprintf(buf,"%s",figSaveName.Data());
	  char *p=strstr(buf,extension);
	  if (p) {
	    sprintf(p,".C");
	    cx->SaveAs(buf);
	  }
	}
	std::cout << "saved" << std::endl;
      }
    }
  }
  std::cout << "plots created" << std::endl;


  std::ofstream htmlfile;
  TString fname=outDir+html_fname;
  std::cout << "saving a file " << fname << "\n";
  htmlfile.open(fname.Data());

  htmlfile << "<!DOCTYPE html" << endl;
  htmlfile << "    PUBLIC \"-//W3C//DTD HTML 3.2//EN\">" << endl;
  htmlfile << "<html>" << endl;
  TString title = (N==1) ? descriptionV[0] : "Comparison";
  htmlfile << "<head><title>" << title << "</title></head>" << endl;
  htmlfile << "<body bgcolor=\"EEEEEE\">" << endl;
  htmlfile << "<h3 style=\"text-align:left; color:DD6600;\">";
  if (N>1) {
    htmlfile << "Comparison of &lt;" << descriptionV[0];
    for (unsigned int i=1; i<descriptionV.size(); ++i) {
      if (i+1 == descriptionV.size()) htmlfile << "&gt; and &lt;"; else htmlfile << "&gt;, &lt;";
      htmlfile << descriptionV[i];
    }
    htmlfile << "&gt;</h3>" << endl;
  }
  else { htmlfile << descriptionV[0] << "</h3>" << endl; }


  for (unsigned int i=0; i<add_html_lines.size(); i++) {
    htmlfile << add_html_lines[i];
  }

  htmlfile << "<br><small>&chi;<sup>2</sup> values in plots are obtained from RooPlot::chiSquare()</small><br><br>" << endl;
  htmlfile << "<table border=\"0\" cellspacing=\"2\" width=\"100%\">" << endl;  
  if (N>1) {
    //htmlfile << "<hr />" << endl;
    htmlfile << "<tr>" << endl;
    for (unsigned int i=0; i<N; ++i) {
      htmlfile << "<td width=\"25%\">" << descriptionV[i] << "</td>" << endl;
    }
    htmlfile << "</tr>" << endl;
    //htmlfile << "<hr />" << endl;
  }

  unsigned int plot_i=0;
  if (combine && combination_names) {
    for (unsigned int ii=0; ii<(*combine).size(); ++ii) {
      htmlfile << "<tr>" << endl;
      for (unsigned int i=0; i<N; i++, plot_i++) {
	//std::cout << "printing (built) ii=" << ii << ", plot_i=" << plot_i << ", figfnames[plot_i]=" << figfnames[plot_i] << "\n";
	htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"" << figfnames[plot_i] << "\"><img src=\"" << figfnames[plot_i] << "\" alt=\"" << figfnames[plot_i] << "\" width=\"100%\"></a></td>" << endl;
      }
      htmlfile << "<td width=\"25%\">" << (*combination_names)[ii] << "</td>" << endl;
    }
    htmlfile << "</tr>" << endl;
    //htmlfile << "<hr />" << endl;
  }


  for (unsigned int i=0; i<N; ++i) {
    if (N>1) htmlfile << " <tr></tr><tr><td width=\"25%\"><h4 style=\"text-align:left; color:DD6600;\">" << descriptionV[i] << " plots</h4></td></tr>" << endl;

    for (unsigned int ii=0; ii<combDataV[i]->size(); ii+=3) {
      htmlfile << "<tr>" << endl;
      for (unsigned int t=0; (t<3) && (ii+t<combDataV[i]->size()); ++t, ++plot_i) {
	//std::cout << "printing (comb) ii=" << ii << ", plot_i=" << plot_i << ", figfnames[plot_i]=" << figfnames[plot_i] << "\n";
	htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"" << figfnames[plot_i] << "\"><img src=\"" << figfnames[plot_i] << "\" alt=\"" << figfnames[plot_i] << "\" width=\"100%\"></a></td>" << endl;
      }
      htmlfile << "</tr>" << endl;  
    }
  }
    
  htmlfile << "</table>" << endl;
  htmlfile << "<hr />" << endl;

  htmlfile << "</body>" << endl;
  htmlfile << "</html>" << endl;
  htmlfile.close();

  std::cout << "HTML file saved" << std::endl;
  return 1;
}

// ======================================================================
