#ifndef MYTOOLS_HH
#define MYTOOLS_HH

#include <TGraphAsymmErrors.h>
#include <TH1F.h>
#include "TCanvas.h"
#include "TString.h"
#include "TSystem.h"
#include "TMatrixD.h"
#include "TVectorD.h"
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

//--------------------------------------------------------------

template<class Type_t>
inline
void HERE(const char *format, const Type_t &a) {
  std::cout << Form(format,a) << std::endl;
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
int printHisto(std::ostream& out, const TH1F* histo, int exponent=0) {
  if (!histo) {
    out << "printHisto: histo is null\n";
    return 0;
  }
  char buf[100];
  const char *format= (exponent) ? 
    " %5.2f-%5.2f    %e    %e\n" : 
    " %5.2f-%5.2f    %f    %f\n";

  out << "values of " << histo->GetName() << "\n";
  for(int i=1; i<=histo->GetNbinsX(); i++) {
    double x=histo->GetBinLowEdge(i);
    double w=histo->GetBinWidth(i);
    sprintf(buf,format,
	    x,x+w,histo->GetBinContent(i),histo->GetBinError(i));
    out << buf;
  }
  return 1;
}

//------------------------------------------------------------------------------------------------------------------------

inline int printHisto(const TH1F* histo, int exponent=0) { return printHisto(std::cout, histo, exponent); }

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
  std::cout << "dims  m[" << m.GetNrows() << ',' << m.GetNcols() << "], "
	    << " mErr[" << mErr.GetNrows() << ',' << mErr.GetNcols() << "]"
	    << std::endl;
  if ((m.GetNrows()!=DYTools::nMassBins) ||
      (m.GetNcols()!=DYTools::nYBinsMax) ||
      (mErr.GetNrows()!=DYTools::nMassBins) ||
      (mErr.GetNcols()!=DYTools::nYBinsMax)) {
    std::cout << "extractMassDependence: expecting matrices "
	      << DYTools::nMassBins << "x" << DYTools::nYBinsMax 
	      << "  instead of  " 
	      << m.GetNrows() << "x" << m.GetNcols() << "\n";
    assert(0);
  }
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
    double factor=1.;
    if (iM==DYTools::nMassBins-1) {
     factor *= DYTools::nYBins[iM]/double(DYTools::nYBins[0]);
    }
    if (perRapidityBinWidth) {
      factor *= (DYTools::yRangeMax-DYTools::yRangeMin)/double(DYTools::nYBins[iM]);
    }
    if (perMassBinWidth) factor/=(DYTools::massBinLimits[iM+1]-DYTools::massBinLimits[iM]);
    h->SetBinContent(iM+1, val*factor);
    h->SetBinError(iM+1, fabs(mErr[iM][iYBin])*factor);
  }
  return h;
}

//-----------------------------------------------------------------

inline
TH1F *extractMassDependenceSpec(const TString &name, const TString &title,
			    const TMatrixD &m, const TMatrixD &mErr,
			    int iYBin,
			const TVectorD &massGrid, const TVectorD &yBinCount,
			    int perMassBinWidth=1, int perRapidityBinWidth=0) {
  if (perRapidityBinWidth) std::cout << "\n\tWARNING: extractMassDependenceSpec: perRapidityBinWidth=1 is experimental\n\n";
  std::cout << "dims  m[" << m.GetNrows() << ',' << m.GetNcols() << "], "
	    << " mErr[" << mErr.GetNrows() << ',' << mErr.GetNcols() << "]"
	    << ", massGrid[" << massGrid.GetNoElements() << "]"
	    << ", yBinCount[" << yBinCount.GetNoElements() << "]"
	    << std::endl;
  TVectorD yBinCountCheck(massGrid.GetNoElements()-1);
  yBinCountCheck=1;
  if (iYBin!=0) {
    std::cout << "extractMassDependenceSpec: iYBin should be 0\n";
    assert(0);
  }
  if ((m.GetNrows()!=massGrid.GetNoElements()-1) ||
      (m.GetNcols()!=1) ||
      (mErr.GetNrows()!=massGrid.GetNoElements()-1) ||
      (mErr.GetNcols()!=1)) {
    std::cout << "extractMassDependenceSpec: expecting matrices "
	      << massGrid.GetNoElements() << "x" << 1
	      << "  instead of  " 
	      << m.GetNrows() << "x" << m.GetNcols() << "\n";
    assert(0);
  }
  TString hName= name + Form("_yBin%d",iYBin);
  TString hTitle= title + Form("_yBin%d",iYBin);
  double *massBins=new double[massGrid.GetNoElements()];
  for (int i=0; i<massGrid.GetNoElements(); ++i) {
    massBins[i]=massGrid[i];
    //std::cout << "massBins[i=" << i << "]=" << massBins[i] << "\n";
  }
  TH1F *h=new TH1F(name,title,massGrid.GetNoElements()-1,massBins);
  delete massBins;

  h->SetDirectory(0);
  h->Sumw2();
  for (int iM=0; iM<m.GetNrows(); ++iM) {
    double val=m[iM][iYBin];
    double factor=1.;
    if (perRapidityBinWidth) {
      factor *= (DYTools::yRangeMax-DYTools::yRangeMin)/double(yBinCount[iM]);
    }
    if (perMassBinWidth) factor/=(massGrid[iM+1]-massGrid[iM]);
    h->SetBinContent(iM+1, val*factor);
    h->SetBinError(iM+1, fabs(mErr[iM][iYBin])*factor);
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

inline void printCSMatrixValues(const TString &name, const TMatrixD &cs, const TMatrixD &csErr, const TMatrixD &csSystErr, int printSystError=1) {
  printYields(name,cs,csErr,csSystErr,printSystError);
}


//------------------------------------------------------------------------------------------------------------------------

inline void printProgress(const char *msg, int idx, int idxMax) {
  double r=trunc(idx/double(idxMax)*1000)*0.1;
  std::cout << msg << idx << "/" << idxMax << " (" << r << "%)\n";
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

template<class Histo_t>
inline
int AppendToHistoName(Histo_t* h, const TString &add) {
  if (!h) {
    std::cout << "AppendToHistoName got null histo ptr\n";
    return 0;
  }
  h->SetName(h->GetName() + add);
  return 1;
}

//------------------------------------------------------------------------------------------------------------------------

inline
void removeError(TH1F* h) {
  std::cout << "nulifying error for " << h->GetName() << "\n";
  for (int ibin=0; ibin<=h->GetNbinsX(); ++ibin) {
    h->SetBinError(ibin,0);
  }
}

//------------------------------------------------------------------------------------------------------------------------

inline 
void subdivideBinWeightByLinearApprox(
   double m_1, double wsigma_1to0, 
   double m0, double wsigma0to1, double wsigma0to1Err,
   double m1, 
   double wsigma1to2, // needed if m2>0
   double m2,
   double mStar,
   double &wSigmaStar1, double &wSigmaStarErr1,
   double &wSigmaStar2, double &wSigmaStarErr2
   ) {
  // linear approximation
  // the cross section is per bin!
  // (m_1,m0, wsigma_1to0); (m0,m1, wsigma0to1) ; (m1,m2, wsigma1to2)
  // We want to divide (m0,m1) to (m0,mStar), (mStar,m2)
  // The cross section at the center of a bin 
  // ( mc_05, wsigma_1to0/w0 ); ( mc05, wsigma0to1/w1 ) ; ( mc15, wsigma1to2/w2 )
  // or
  // ( mc_05, sigma_1to0 ) ; ( mc05, sigma0to1 ); ( mc15, sigma1to2 )
  // defines linear dependencies
  //  sigma= (sigma0to1 - sigma_1to0)/(mc05-mc_05) * ( m - mc_05 ) + sigma_1to0
  // and
  //  sigma= (sigma1to2 - sigma0to1)/(mc15-mc05) * ( m - mc05 ) + sigma0to1

  const int debug=0;
  const int warn_on_branches=0;
  const int pure_linear_branch=0; // recommended is 0 -- preserves total count in all cases
  const double thr=-9e8;  // if m2<thr, ignore sigma1to2

  if (debug) {
    std::cout << "subdivideBinWeightByLinearApprox: \n";
    std::cout << "  m_1=" << m_1 << ", m0=" << m0 << ", wsigma_1to0=" << wsigma_1to0 << "\n";
    std::cout << "  m0=" << m0 << ", m1=" << m1 << ", wsigma0to1=( " << wsigma0to1 << " pm " << wsigma0to1Err << " )\n";
    if (m2<thr) std::cout << " m2, wsigma1to2 are ignored\n";
    else std::cout << "  m1=" << m1 << ", m2=" << m2 << ", wsigma1to2=" << wsigma1to2 << "\n";
    std::cout << "  mStar=" << mStar << "\n";
    std::cout << "\n";
  }

  double mc_05=0.5*(m_1 + m0);
  double w0=m0-m_1;
  double w1=m1-m0;
  double sigma_1to0 = wsigma_1to0 / w0;
  double sigma0to1 =  wsigma0to1 / w1;
  double slope=(sigma0to1 - sigma_1to0)*2./(m1 - m_1);
  if (debug) std::cout << "slope=" << slope << " (1)\n";
  if ((m2>thr) && pure_linear_branch) {
    // Try to improve using slope with respect to the next bin
    double w2=m2-m1;
    double sigma1to2 =  wsigma1to2 / w2;
    double slopePlus= (sigma1to2 - sigma0to1) * 2./(m2 - m0);
    slope = 0.5*( slope + slopePlus );
    if (debug) std::cout << "slope=" << slope << " (2)\n";
  }

  if ((mStar>=m0) && (mStar<=m1)) {     // interpolation

    // center of the bin to consider, and the width
    double mcStar1=0.5*(m0 + mStar);
    double wStar1= mStar-m0;
    if (debug) std::cout << " mcStar1=" << mcStar1 << ", wStar1=" << wStar1 << "\n";
    // the cross section there
    double sigma1_div_width= slope*( mcStar1 - mc_05 ) + sigma_1to0;
    if (debug) std::cout << " sigma1_div_width=" << sigma1_div_width << "\n";
    // 1st answer
    wSigmaStar1 = sigma1_div_width * wStar1;

    double wStar2= m1-mStar;
    if (pure_linear_branch) {  // Pure linear approximation - use formula
      if (warn_on_branches) std::cout << " branch: pure linear approximation\n";
      // center to the other part of the bin, and the width
      double mcStar2=0.5*(mStar + m1);
      double sigma2_div_width= slope*( mcStar2 - mc_05 ) + sigma_1to0;
      // 2nd answer
      wSigmaStar2 = sigma2_div_width* wStar2;
    }
    else {  // Take the remaining counts
      if (warn_on_branches) std::cout << " branch: taking remaining area\n";
      wSigmaStar2= wsigma0to1 - wSigmaStar1;
    }
    
    // Error calculation
    double err=wsigma0to1Err/w1;
    double frac1=  wStar1 / sqrt( wStar1*wStar1 + wStar2*wStar2 );
    double frac2=  wStar2 / sqrt( wStar1*wStar1 + wStar2*wStar2 );
    wSigmaStarErr1= frac1*err *wStar1;
    wSigmaStarErr2= frac2*err *wStar2;

  }
  else { // extrapolation
    wSigmaStar1=0; wSigmaStarErr1=0;
    wSigmaStar2=0; wSigmaStarErr2=0;
    if (mStar>m1) {
      if (warn_on_branches) std::cout << " branch:  Extrapolation to the right\n";
      double mc01=0.5*(m0 + m1);
      double sigmaAtM1= slope*( m1 - mc01 ) + sigma0to1;
      double sigmaAtMStar= slope*( mStar - mc01 ) + sigma0to1;
      double w= mStar-m1;
      if (debug) std::cout << "sigmaAtM1=" << sigmaAtM1 << ", sigmaAtMStar=" << sigmaAtMStar << ", w=" << w << "\n";
      wSigmaStar1 = 0.5*(sigmaAtM1 + sigmaAtMStar) * w;
      wSigmaStarErr1= (wsigma0to1Err/w1) *w;
    }
    else {
      // this is not a branch -- this is an important warning
      std::cout << "\n\tExtrapolation to the left is not possible by construction\n";
    }
  }

  return;
}

// --------------------------------
// --------------------------------

inline
TH1F* createZpeakHisto(const char *hname="hZpeak_mass", const char *htitle="Z-peak region") {
  const int mbCount=13;
  const double mbins[mbCount+1]={  60,  64,  68,  72, 76, 
				   81,  86,  91,  96,101,
				  106, 110, 115, 120 };
  TH1F *h= new TH1F(hname,htitle,mbCount,mbins);
  h->Sumw2();
  h->SetDirectory(0);
  std::cout << "mcCount=13\n";
  return h;
}

// --------------------------------

inline
TH1F* createZpeakHisto1(const char *hname="hZpeak_mass", const char *htitle="Z-peak region", double mass_min=60., double mass_max=120., int binCount=60) {
  TH1F *h= new TH1F(hname,htitle,binCount,mass_min,mass_max);
  h->Sumw2();
  h->SetDirectory(0);
  return h;
}

// --------------------------------

template<class TH1F_t>
inline
int printZpeakInfo(TH1F_t *h) {
  if (!h) return 0;
  std::cout << Form("(mean,mean_err)=(%6.4lf,%6.4lf)", h->GetMean(),h->GetMeanError())
            << " in " << h->GetName() << "\n";
  return 1;
}

// --------------------------------

inline
TString DayAndTimeTag()
{
   time_t ltime;
   ltime=time(NULL);
   TString str = TString(asctime(localtime(&ltime)));
   str.ReplaceAll(" ","_");
   str.ReplaceAll(":","");
   return str;
}

//------------------------------------------------------------------------------------------------------------------------

#endif
