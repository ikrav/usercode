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

void SaveCanvas(TCanvas* canv, TString canvName)
{        
  gSystem->mkdir("plots",kTRUE);
  gSystem->mkdir("plots/png",kTRUE);
  gSystem->mkdir("plots/pdf",kTRUE);
  gSystem->mkdir("plots/root",kTRUE);

  TString saveName="plots/png/";
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

#endif
