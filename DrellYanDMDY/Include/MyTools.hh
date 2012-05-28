#ifndef MYTOOLS_HH
#define MYTOOLS_HH

#include <TGraphAsymmErrors.h>
#include <TH1F.h>
#include "TCanvas.h"
#include "TString.h"
#include "TSystem.h" 
#include <iostream>

//#include "RooStats/FeldmanCousins.h"

namespace toolbox 
{
// Double_t calcEff(const Int_t    pass, const Int_t    total, Double_t *errl=0, Double_t *errh=0, Int_t method=0);
// Double_t calcEff(const Double_t pass, const Double_t total, Double_t *errl=0, Double_t *errh=0, Int_t method=0);

Double_t deltaR(const Double_t eta1, const Double_t phi1, const Double_t eta2, const Double_t phi2);

Double_t deltaPhi(const Double_t phi1, const Double_t phi2);

Int_t roundToInt(const Double_t x);
}

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
Double_t toolbox::deltaR(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2)
{
  const Double_t pi = 3.14159265358979;
  Double_t dphi = fabs(phi1-phi2);
  while (dphi>pi)
    dphi = fabs(dphi - 2.0*pi);
    
  Double_t deta = eta1-eta2;
  
  return sqrt(dphi*dphi + deta*deta);
}

//------------------------------------------------------------------------------------------------------------------------
Double_t toolbox::deltaPhi(Double_t phi1, Double_t phi2) 
{
  // Compute dPhi between two given angles. Results is in [0,pi].
  const Double_t pi = 3.14159265358979;
  Double_t dphi = fabs(phi1-phi2);
  while (dphi>pi)
    dphi = fabs(dphi - 2.0*pi);

  return dphi;
}

//------------------------------------------------------------------------------------------------------------------------

Int_t toolbox::roundToInt(Double_t x)
{
  if(x>0)
    return ((x-floor(x)) < (ceil(x)-x)) ? (Int_t)floor(x) : (Int_t)ceil(x);
  else
    return ((x-floor(x)) < (ceil(x)-x)) ? (Int_t)ceil(x) : (Int_t)floor(x);
}

//------------------------------------------------------------------------------------------------------------------------

template<class T>
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

void SaveCanvas(TCanvas* canv, const TString &canvName, const TString &destDir="plots")
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
void ClearVec(std::vector<T*> &vec) {
  for (unsigned int i=0; i<vec.size(); ++i) if (vec[i]) delete vec[i];
  vec.clear();
}

//------------------------------------------------------------------------------------------------------------------------

void HERE(const char *msg) {
  std::cout << ((msg) ? msg : "HERE") << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------

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

#endif
