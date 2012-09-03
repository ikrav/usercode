#ifndef MYTOOLS_HH
#define MYTOOLS_HH

#include <TGraphAsymmErrors.h>
#include <TH1F.h>
#include "TCanvas.h"
#include "TString.h"
#include "TSystem.h"
#include "TMatrixD.h"
#include <iostream>
#include <math.h>
#include "../Include/CPlot.hh"
#include "../Include/DYTools.hh"

//#include "RooStats/FeldmanCousins.h"

//namespace toolbox 
//{
// Double_t calcEff(const Int_t    pass, const Int_t    total, Double_t *errl=0, Double_t *errh=0, Int_t method=0);
// Double_t calcEff(const Double_t pass, const Double_t total, Double_t *errl=0, Double_t *errh=0, Int_t method=0);

//Double_t deltaR(const Double_t eta1, const Double_t phi1, const Double_t eta2, const Double_t phi2);

//Double_t deltaPhi(const Double_t phi1, const Double_t phi2);

//Int_t roundToInt(const Double_t x);
//}

//------------------------------------------------------------------------------------------------------------------------
// Double_t toolbox::calcEff(Int_t pass, Int_t total, Double_t *errl, Double_t *errh, Int_t method)
// {
//   // method: 0 -> Bayes Divide
//   //         1 -> Feldman-Cousins 
//   //         2 -> Clopper-Pearson
  
//   Double_t r = (total>0) ? (Double_t)pass/(Double_t)total : 0;
//   if(errl) *errl = 0;
//   if(errh) *errh = 0;
    
//   const Double_t conf = 0.68269;
  
//   if(method==0) {    
//     TGraphAsymmErrors g;
//     Double_t mode, low, high;
//     g.Efficiency(pass,total,conf,mode,low,high);
//     if(errl) *errl = mode - low;
//     if(errh) *errh = high - mode;
//   }

//   if(method==1) {
// //    FeldmanCousins fc;
// //    fc.SetConfidenceLevel(conf);
//   }

//   if(method==2) {
//   }
 
//   return r;
// }

// Double_t toolbox::calcEff(Double_t pass, Double_t total, Double_t *errl, Double_t *errh, Int_t method) 
// {
//   // Round values to whole numbers first
//   return calcEff(roundToInt(pass),roundToInt(total),errl,errh,method);
// }

//------------------------------------------------------------------------------------------------------------------------
namespace toolbox {

inline
Double_t deltaR(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2)
{
  const Double_t pi = 3.14159265358979;
  Double_t dphi = fabs(phi1-phi2);
  while (dphi>pi)
    dphi = fabs(dphi - 2.0*pi);
    
  Double_t deta = eta1-eta2;
  
  return sqrt(dphi*dphi + deta*deta);
}

}
//------------------------------------------------------------------------------------------------------------------------
namespace toolbox {

inline
Double_t deltaPhi(Double_t phi1, Double_t phi2) 
{
  // Compute dPhi between two given angles. Results is in [0,pi].
  const Double_t pi = 3.14159265358979;
  Double_t dphi = fabs(phi1-phi2);
  while (dphi>pi)
    dphi = fabs(dphi - 2.0*pi);

  return dphi;
}

}
//------------------------------------------------------------------------------------------------------------------------
namespace toolbox {

inline
Int_t roundToInt(Double_t x)
{
  if(x>0)
    return ((x-floor(x)) < (ceil(x)-x)) ? (Int_t)floor(x) : (Int_t)ceil(x);
  else
    return ((x-floor(x)) < (ceil(x)-x)) ? (Int_t)ceil(x) : (Int_t)floor(x);
}

}
//------------------------------------------------------------------------------------------------------------------------

template<class T>
inline
void PrintVec(const char *msg, const std::vector<T>& vec, int prneol=0) {
  if (msg) std::cout << msg;
  std::cout << "vec[" << vec.size() << "]: ";
  for (unsigned int i=0; i<vec.size(); ++i) {
    if (prneol) std::cout << "\n" << i << ") ";
    std::cout << " " << vec[i];
  }
  if (prneol) std::cout << "\n";
}

//------------------------------------------------------------------------------------------------------------------------

inline
void SaveCanvas(TCanvas* canv, const TString &canvName, TString destDir=CPlot::sOutDir) 
{
  gSystem->mkdir(destDir,kTRUE);
  gSystem->mkdir(destDir+TString("/png"),kTRUE);
  gSystem->mkdir(destDir+TString("/pdf"),kTRUE);
  gSystem->mkdir(destDir+TString("/root"),kTRUE);

  TString saveName=destDir+TString("/png/");
  saveName+=canvName;
  saveName+=".png";
  canv->SaveAs(saveName);
  saveName.ReplaceAll("png","pdf");
  canv->SaveAs(saveName);
  saveName.ReplaceAll("pdf","root");
  canv->SaveAs(saveName);
  return;
}
// ----------------------------------------------------------

template<class T>
inline
void ClearVec(std::vector<T*> &vec) {
  for (unsigned int i=0; i<vec.size(); ++i) if (vec[i]) delete vec[i];
  vec.clear();
}

//------------------------------------------------------------------------------------------------------------------------

inline
void HERE(const char *msg) {
  std::cout << ((msg) ? msg : "HERE") << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------

inline bool PosOk(size_t pos) { return (pos==std::string::npos) ? 0 : 1; }

//----------------------------------------------------------------------

template<class T>
inline bool PosOk(const std::string &s, const T& substr) { 
  return (s.find(substr)==std::string::npos) ? 0 : 1; 
}

//------------------------------------------------------------------------------------------------------------------------

inline
int printHisto(std::ostream& out, const TH1F* histo) {
  if (!histo) {
    out << "printHisto: histo is null\n";
    return 0;
  }
  char buf[100];
  out << "values of " << histo->GetName() << "\n";
  for(int i=1; i<=histo->GetNbinsX(); i++) {
    double x=histo->GetBinLowEdge(i);
    double w=histo->GetBinWidth(i);
    sprintf(buf," %5.2f-%5.2f    %f    %f\n",
	    x,x+w,histo->GetBinContent(i),histo->GetBinError(i));
    out << buf;
  }
  return 1;
}

//------------------------------------------------------------------------------------------------------------------------

inline int printHisto(const TH1F* histo) { return printHisto(std::cout, histo); }

//------------------------------------------------------------------------------------------------------------------------

inline
TH1F *extractRapidityDependence(const TString &name, const TString &title,
				const TMatrixD &m, const TMatrixD &mErr,
				int iMassBin, int perMassBinWidth=0) {
  TString hName= name + Form("_massBin%d",iMassBin);
  TString hTitle= title + Form("_massBin%d",iMassBin);
  TH1F *h=new TH1F(name,title,DYTools::nYBins[iMassBin],DYTools::yRangeMin,DYTools::yRangeMax);
  h->SetDirectory(0);
  h->Sumw2();
  for (int iY=0; iY<DYTools::nYBins[iMassBin]; ++iY) {
    double factor= (perMassBinWidth==0) ? 
      1 : 1/(DYTools::massBinLimits[iMassBin+1] - DYTools::massBinLimits[iMassBin]);
    h->SetBinContent(iY+1, m[iMassBin][iY]*factor);
    h->SetBinError(iY+1, fabs(mErr[iMassBin][iY])*factor);
  }
  return h;
}

// -----------------------------------------------------------------------------

inline
TH1F *extractMassDependence(const TString &name, const TString &title,
			    const TMatrixD &m, const TMatrixD &mErr,
			    int iYBin, 
			    int perMassBinWidth=1, int perRapidityBinWidth=0) {
  if (perRapidityBinWidth) std::cout << "\n\tWARNING: extractMassDependence: perRapidityBinWidth=1 is experimental\n\n";
  TString hName= name + Form("_yBin%d",iYBin);
  TString hTitle= title + Form("_yBin%d",iYBin);
  TH1F *h=new TH1F(name,title,DYTools::nMassBins,DYTools::massBinLimits);
  h->SetDirectory(0);
  h->Sumw2();
  for (int iM=0; iM<DYTools::nMassBins; ++iM) {
    bool lastBin= (iM==DYTools::nMassBins-1) ? true : false;  
    if (lastBin && (iYBin!=0)) {
      double yc=DYTools::findAbsYValue(0,iYBin);
      iYBin=DYTools::findAbsYBin(iM,yc);
    }
    double val=m[iM][iYBin];
    if (iM==DYTools::nMassBins-1) {
      val *= DYTools::nYBins[iM]/double(DYTools::nYBins[0]);
    }
    if (perRapidityBinWidth) {
      val *= (DYTools::yRangeMax-DYTools::yRangeMin)/double(DYTools::nYBins[iM]);
    }
    if (perMassBinWidth) val/=(DYTools::massBinLimits[iM+1]-DYTools::massBinLimits[iM]);
    h->SetBinContent(iM+1, val);
    h->SetBinError(iM+1, fabs(mErr[iM][iYBin]));
  }
  return h;
}



//-----------------------------------------------------------------

inline
void printYields(const TString &name, const TMatrixD &cs, const TMatrixD &csErr, const TMatrixD &csErrSyst, int printSystError=1) {
  std::cout << "\nprintYields for name=<" << name << ">\n";
  if ((cs.GetNrows() != csErr.GetNrows()) ||
      (cs.GetNrows() != csErrSyst.GetNrows())) {
    printf(" -- numbers of rows is different: %d, %d, %d\n",cs.GetNrows(),csErr.GetNrows(),csErrSyst.GetNrows());
    assert(0);
  }
  if ((cs.GetNcols() != csErr.GetNcols()) ||
      (cs.GetNcols() != csErrSyst.GetNcols())) {
    printf(" -- numbers of rows is different: %d, %d, %d\n",cs.GetNcols(),csErr.GetNcols(),csErrSyst.GetNcols());
    assert(0);
  }
  std::cout << "cs.GetNcols()=" << cs.GetNcols() << ", cs.GetNrows="<< cs.GetNrows() << "\n";
  if (printSystError) {
    for (int ir=0; ir<cs.GetNrows(); ++ir) {
      for (int ic=0; ic<cs.GetNcols(); ++ic) {
	printf(" %6.4lf %8.4e %8.4e\n",cs[ir][ic],csErr[ir][ic],csErrSyst[ir][ic]);
      }
      printf("\n");
    }
  }
  else {
    for (int ir=0; ir<cs.GetNrows(); ++ir) {
      for (int ic=0; ic<cs.GetNcols(); ++ic) {
	printf(" %6.4lf %8.4e\n",cs[ir][ic],csErr[ir][ic]);
      }
      printf("\n");
    }
  }
  std::cout << "done" << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------

inline
void printSanityCheck(const TMatrixD &val, const TMatrixD &err, const TString &name)
{
  using namespace DYTools;
  std::cout<<"Sanity check printout"<<std::endl;
  if (val.GetNrows()!=nMassBins || val.GetNcols()!=findMaxYBins())
    {
      std::cout<<name<<" matrix has wrong size"<<std::cout;
      return;
    }
  if (err.GetNrows()!=nMassBins || err.GetNcols()!=findMaxYBins())
    {
      std::cout<<name<<"Err matrix has wrong size"<<std::cout;
      return;
    }
  std::cout<<"Nan values of "<<name<<" or/and "<<name<<"Err:"<<std::endl;
  for(int i=0; i<DYTools::nMassBins; i++)
    for (int yi=0; yi<nYBins[i]; ++yi) 
      {
        if ( (val(i,yi)!=val(i,yi)) || (err(i,yi)!=err(i,yi)) )
           std::cout<<name<<"("<<i<<","<<yi<<")="<<val(i,yi)<<", "<<name<<"Err("<<i<<","<<yi<<")="<<err(i,yi)<<std::endl;
      }
  std::cout<<"Large errors ("<<name<<"Errv>0.1*"<<name<<" or "<<name<<"Err>0.1) :"<<std::endl;
  for(int i=0; i<DYTools::nMassBins; i++)
    for (int yi=0; yi<nYBins[i]; ++yi) 
      {
        if ( fabs(err(i,yi))>0.1*fabs(val(i,yi)) || fabs(err(i,yi))>0.1)
           std::cout<<name<<"("<<i<<","<<yi<<")="<<val(i,yi)<<", "<<name<<"Err("<<i<<","<<yi<<")="<<err(i,yi)<<std::endl;
      }
}

//------------------------------------------------------------------------------------------------------------------------

#endif
