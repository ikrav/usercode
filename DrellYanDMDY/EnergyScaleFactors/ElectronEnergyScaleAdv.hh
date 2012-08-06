#ifndef ElectronEnergyScaleAdv_H
#define ElectronEnergyScaleAdv_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#ifndef __myLib__ 
#define UsePlotsAcc
#endif

#define DoSmearing

#ifdef UsePlotsAcc
#include "TPlotsAcc.hh"
#include "RooDataSet.h"
#endif

#ifdef DoSmearing
#include <Math/ProbFuncMathCore.h>
#ifndef UseTH1D
#define UseTH1D
#endif
#endif

#ifdef UseTH1D
#include <TH1D.h>
namespace mithep {
  typedef TH1D myHistoClass_t;
}
#endif
#ifdef UseTH1F
#include <TH1F.h>
namespace mithep {
  typedef TH1F myHistoClass_t;
}
#endif

#include <RooRealVar.h>
#include <TApplication.h>


#include "../Include/EtaEtaMass.hh"
#include "../Include/ElectronEnergyScale.hh" // class is used indirectly

typedef enum { _ESFWC_none=0, _ESFWC_1bin=1, _ESFWC_2bin, _ESFWC_3bin, _ESFWC_4bin, _ESFWC_5bin,_ESFWC_6bin, _ESFWC_3EB3EE, _ESFWC_3EB3EEa, _ESFWC_4EB3EE, _ESFWC_BothSided, _ESFWC_2binNegs, _ESFWC_3binNegs, _ESFWC_4binNegs, _ESFWC_5binNegs, _ESFWC_6binNegs, _ESFWC_3EB3EENegs, _ESFWC_4EB3EENegs, _ESFWC_5EB3EENegs, _ESFWC_Separate, _ESFWC_5EBNegs, _ESFWC_6EBNegs, _ESFWC_4EENegs, _ESFWC_6EB3EENegs, _ESFWC_7EB3EENegs, _ESFWC_8EB3EENegs  } TEnergyScaleFactorsWorkCase_t;
typedef enum { _ESFWC_BEall, _ESFWC_EB=1, _ESFWC_EE, _ESFWC_EBNeg, _ESFWC_EENeg } TEnergyScaleFactorsBECase_t;
typedef enum { _ESFWC_any, _ESFWC_EBEB, _ESFWC_EBEE, _ESFWC_EEEE /*, _ESFWC_EBNegEBNeg, _ESFWC_EBNegEENeg, _ESFWC_EENegEENeg*/ , _ESFWC_all } TEnergyScaleFactorsBEDoubleCase_t;

typedef enum { _ESFModel_none=0, _ESFModel_gauss, _ESFModel_gaussExp, _ESFModel_doubleGauss, _ESFModel_bifurGauss, _ESFModel_cball, _ESFModel_breitWigner, _ESFModel_voigtian, _ESFModel_has_shift, _ESFModel_gaussShift, _ESFModel_cballShift } TEnergyScaleFactorsFitModel_t;

#ifndef __myLib__
TEnergyScaleFactorsWorkCase_t IdentifyESFWorkCase(const char *s);
TEnergyScaleFactorsFitModel_t IdentifyESFFitModel(const char *s, int verb=0);
#endif

struct ElectronEnergyScaleAdv_t;

int DetermineESFWorkCase(const char *filename, ElectronEnergyScaleAdv_t &esf);
TEnergyScaleFactorsWorkCase_t DetermineESFWorkCase(const TApplication &theApp);
TEnergyScaleFactorsFitModel_t DetermineESFFitModel(const TApplication &theApp);

// -----------------------------------------------------------------------
// -----------------------------------------------------------------------


struct ElectronEnergyScaleAdv_t {
  int FWorkCase;
  int FFitModel;
  int FScFCount, FSmFCount; // scale factors count, smear factors count
  int FEtaDivisionCount;
  int FInvertedSF; // inverted scale factors
  double FMCtoData;
  double *FEtaBinLimits;
  double *FScaleFactors;
  double *FSmearings;
  double *FScaleFactorErrLo;
  double *FScaleFactorErrHi;
  double *FSmearingErrLo;
  double *FSmearingErrHi;
  ElectronEnergyScaleAdv_t() : 
    FWorkCase(0), FFitModel(0), FScFCount(0), FSmFCount(0),
    FEtaDivisionCount(0), FInvertedSF(0), FMCtoData(0),
    FEtaBinLimits(0), 
    FScaleFactors(0), FSmearings(0), 
    FScaleFactorErrLo(0), FScaleFactorErrHi(0),
    FSmearingErrLo(0), FSmearingErrHi(0) {}

  ElectronEnergyScaleAdv_t(const ElectronEnergyScaleAdv_t &sf) : 
    FWorkCase(0), FFitModel(0), FScFCount(0), FSmFCount(0),
    FEtaDivisionCount(0), FInvertedSF(0), FMCtoData(0),
    FEtaBinLimits(0), 
    FScaleFactors(0), FSmearings(0),
    FScaleFactorErrLo(0), FScaleFactorErrHi(0),
    FSmearingErrLo(0), FSmearingErrHi(0)
       { this->Assign(sf); }

  //ElectronEnergyScaleAdv_t(const TString &escaleTagName) : ElectronEnergyScale(escaleTagName) {}

  ~ElectronEnergyScaleAdv_t() { Clear(); }

  void Clear() {
    FWorkCase=0; FFitModel=0; FScFCount=0; FSmFCount=0;
    FMCtoData=0; FInvertedSF=0; FEtaDivisionCount=0;
    if (FEtaBinLimits) delete [] FEtaBinLimits;
    if (FScaleFactors) delete [] FScaleFactors;
    if (FScaleFactorErrLo) delete [] FScaleFactorErrLo;
    if (FScaleFactorErrHi) delete [] FScaleFactorErrHi;
    if (FSmearings) delete [] FSmearings;
    if (FSmearingErrLo) delete [] FSmearingErrLo;
    if (FSmearingErrHi) delete [] FSmearingErrHi;
    FEtaBinLimits=0; FScaleFactors=0; FSmearings=0;
    FScaleFactorErrLo=0; FScaleFactorErrHi=0;
    FSmearingErrLo=0; FSmearingErrHi=0;
  }

  int WorkCase() const { return FWorkCase; }
  TEnergyScaleFactorsWorkCase_t WorkCaseKind() const { return TEnergyScaleFactorsWorkCase_t(FWorkCase); }
  int WorkCase(int wcase, int fitModel);  // The first method to call!!
  int FitModel() const { return FFitModel; }
  TEnergyScaleFactorsFitModel_t FitModelKind() const { return TEnergyScaleFactorsFitModel_t(FFitModel); }
  int ScFCount() const { return FScFCount; }
  int SmFCount() const { return FSmFCount; }
  double MCtoData() const { return FMCtoData; }
  void MCtoData(double x) { FMCtoData=x; }
  int EtaDivisionCount() const { return FEtaDivisionCount; }
  int BinCount() const { return FEtaDivisionCount; }
  int ExpectCount() const { return FEtaDivisionCount*(FEtaDivisionCount+1)/2; }
  //int Count() const { return FEtaDivisionCount; }
  //int Valid(int i) const { return ((i>=0) && (i<FEtaDivisionCount)) ? 1:0; }
  int ValidScFIdx(int i) const { return ((i>=0) && (i<FScFCount)) ? 1:0; }
  int ValidSmFIdx(int i) const { return ((i>=0) && (i<FSmFCount)) ? 1:0; }
  int ValidBinIdx(int i) const { return ((i>=0) && (i<=FEtaDivisionCount)) ? 1:0; }
  const double* EtaBinLimits() const { return FEtaBinLimits; }
  double EtaBinLimit(int i) const { return (ValidBinIdx(i)) ? FEtaBinLimits[i] : 9999; }
  const double* ScaleFactors() const { return FScaleFactors; }
  double *EditScaleFactors() { return FScaleFactors; }
  int InvertedSF() const { return FInvertedSF; }
  void InvertedSF(int inv) { FInvertedSF=inv; }
  double ScaleFactor(int i) const { return (ValidScFIdx(i)) ? FScaleFactors[i] : 9999.; }

  double ProperScaleFactorValue(int idx, int shift_i=0) const {
    int i=idx+shift_i*FEtaDivisionCount;
    if (!ValidScFIdx(i)) return 9999.;
    double x=9999.;
    int selection=(shift_i==0) ? FInvertedSF : 0;
    switch(selection) {
    case 0: x=FScaleFactors[i]; break;
    case 1: x=1/FScaleFactors[i]; break;
    case 2: x=sqrt(FScaleFactors[i]); break;
    case 3: x=1/sqrt(FScaleFactors[i]); break;
    case 4: x=FScaleFactors[i]*FScaleFactors[i]; break;
    default:
      std::cout << "ProperScaleFactorValue: not ready for FInvertedSF=" << FInvertedSF << "\n";
    }
    return x;
  }

  double ProperScaleFactorValueForPrint(int idx, int shift_i=0) const {
    int i=idx+shift_i*FEtaDivisionCount;
    if (!ValidScFIdx(i)) return 9999.;
    double x=9999.;
    int selection=(shift_i==0) ? FInvertedSF : 0;
    switch(selection) {
    case 0: x=FScaleFactors[i]; break;
    case 1: x=1/FScaleFactors[i]; break;
    case 2: x=FScaleFactors[i]; break;
    case 3: x=FScaleFactors[i]; break;
    default:
      std::cout << "ProperScaleFactorValue: not ready for FInvertedSF=" << FInvertedSF << "\n";
    }
    return x;
  }

  double ProperScaleFactorErrValueForPrint(int i) const {
    if (!ValidScFIdx(i)) return 9999.;
    double x=9999.;
    switch(FInvertedSF) {
    case 2: case 3: x=ScaleFactorErr(i); break;
    case 1: case 0: break;
    default:
      std::cout << "ProperScaleFactorErrValue: not ready for FInvertedSF=" << FInvertedSF << "\n";
    }
    return x;
  }

  std::string ProperScaleFactorValueWithErrForPrint(int i) const;
  std::string ProperSmearFactorValueWithErrForPrint(int i) const;

  double ProperSmearFactorValue(int idx, int shift_i=0) const {
    int i=idx+shift_i*FEtaDivisionCount;
    if (!ValidSmFIdx(i)) return 9999.0;
    return FSmearings[i];
  }

  int EtaBinIdx(double eta) const { 
    int idx=-1;
    if (this->WorkCaseKind()<_ESFWC_BothSided) { eta=fabs(eta); }
    for (int i=0; (idx==-1) && (i<FEtaDivisionCount); ++i) {
      if ((eta>=FEtaBinLimits[i]) && (eta<FEtaBinLimits[i+1])) {
	//std::cout << "located eta in bin i=" << i << "\n";
	idx=i;
      }
    }
    return idx;
  }

  int PrepareEtaEtaBinIdx(int etaBin1, int etaBin2) const {
    int bin=-1;
    if (ValidBinIdx(etaBin1) && ValidBinIdx(etaBin2)) {
      // i<=j, idx = [ sum_{k=0}^{k < i-1} (N-k) ] + j-i = (2N-1-i)i/2 + j
      if (etaBin1>etaBin2) { bin=etaBin1; etaBin1=etaBin2; etaBin2=bin; }
      bin=(2*FEtaDivisionCount-1-etaBin1)*etaBin1/2 + etaBin2;
    }
    return bin;
  }

  int PrepareEtaEtaBinIdx(double eta1, double eta2) const {
    int etaBin1=this->EtaBinIdx(eta1);    
    int etaBin2=this->EtaBinIdx(eta2);
    int bin= PrepareEtaEtaBinIdx(etaBin1,etaBin2);
    //std::cout << "
    return bin;
  }

  int PrepareEtaBinIdx(double eta1, double eta2, int &etaBin1, int &etaBin2) const {
    etaBin1=this->EtaBinIdx(eta1);    etaBin2=this->EtaBinIdx(eta2);
    if (!ValidBinIdx(etaBin1) || !ValidBinIdx(etaBin2)) {
      std::cout << "PrepareEtaBinIdx: Failed to prepare eta bin indices for eta=(" << eta1 << ',' << eta2 << "). Deduced bin idx=(" << etaBin1 << ',' << etaBin2 << ")\n";
      return 0;
    }
    return 1;
  }

  double GetScaleFactorValue(int eta_idx_i, int eta_idx_j,int shift=0) const {
    double val1=ProperScaleFactorValue(eta_idx_i,shift);
    double val2=ProperScaleFactorValue(eta_idx_j,shift);
    return (shift==0) ? (val1*val2) : (val1+val2);
  }
  double GetScaleFactorValue(double eta_i, double eta_j, int shift=0) const {
    int etaBin1,etaBin2;
    if (this->PrepareEtaBinIdx(eta_i,eta_j, etaBin1,etaBin2)) {
      return GetScaleFactorValue(etaBin1,etaBin2,shift);
    }
    return 0;
  }

  double GetSmearFactorValue(int eta_idx_i, int eta_idx_j,int shift=0) const {
    double f1=ProperSmearFactorValue(eta_idx_i,shift);
    double f2=ProperSmearFactorValue(eta_idx_j,shift);
    //std::cout << "eta_idx_i=" << eta_idx_i << ", eta_idx_j=" << eta_idx_j << ", f1=" << f1 << ", f2=" << f2 << "\n";
    return sqrt(f1*f1+f2*f2);
  }
  double GetSmearFactorValue(double eta_i, double eta_j, int shift=0) const {
    int etaBin1,etaBin2;
    if (this->PrepareEtaBinIdx(eta_i,eta_j, etaBin1,etaBin2)) {
      return GetSmearFactorValue(etaBin1,etaBin2,shift);
    }
    return 0;
  }


  static int GetEtaBinLimitCount(TEnergyScaleFactorsWorkCase_t theCase);
  static TEnergyScaleFactorsWorkCase_t GetNegsWorkCaseKind(TEnergyScaleFactorsWorkCase_t wcase);
  static int GetScaleFactorsCount(TEnergyScaleFactorsWorkCase_t wcase, TEnergyScaleFactorsFitModel_t fitModel);
  static int GetSmearFactorsCount(TEnergyScaleFactorsWorkCase_t wcase, TEnergyScaleFactorsFitModel_t fitModel);
  static int FillEtaBinLimits(TEnergyScaleFactorsWorkCase_t theCase, double *b);

#ifdef __hiLevelService__
  int GetSmearIdx_gamma(int &idx_min, int &idx_max) const; // Breit-Wigner, Voigtian
  int GetSmearIdx_sigma(int &idx_min, int &idx_max) const; // Gaussian, GaussExp, Crystal ball
  int GetSmearIdx_smearSigma(int &idx_min, int &idx_max) const; // Voigtian
  int GetSmearIdx_sigmaL(int &idx_min, int &idx_max) const; // BifurGauss
  int GetSmearIdx_sigmaR(int &idx_min, int &idx_max) const; // BifurGauss
  int GetSmearIdx_tau(int &idx_min, int &idx_max) const; // GaussExp
  int GetSmearIdx_parN(int &idx_min, int &idx_max) const; // Crystal ball
  int GetSmearIdx_alpha(int &idx_min, int &idx_max) const; // Crystal ball

  int GetSmearIdxShift_gamma() const { return 0; } // Breit-Wigner, Voigtian
  int GetSmearIdxShift_sigma() const { return 0; } // Gaussian, GaussExp, Crystal ball
  int GetSmearIdxShift_smearSigma() const { return 1; } // Voigtian
  int GetSmearIdxShift_sigmaL() const { return 0; } // BifurGauss
  int GetSmearIdxShift_sigmaR() const { return 1; } // BifurGauss
  int GetSmearIdxShift_tau() const { return 0; } // GaussExp
  int GetSmearIdxShift_parN() const { return 1; } // Crystal ball
  int GetSmearIdxShift_alpha() const { return 1; } // Crystal ball
#endif


  // prepare arrays of indices to collate the data
  static std::vector<int>* BinIndices(TEnergyScaleFactorsWorkCase_t theCase, TEnergyScaleFactorsBECase_t eb);
  static std::vector<int>* DoubleBinIndices(TEnergyScaleFactorsWorkCase_t theCase, TEnergyScaleFactorsBEDoubleCase_t ebeb);
  static std::string WorkCaseName(TEnergyScaleFactorsWorkCase_t wcase);
  static std::string WorkCaseShortName(TEnergyScaleFactorsWorkCase_t wcase);
  static std::string FitModelName(TEnergyScaleFactorsFitModel_t fmodel);
  static std::string FitModelShortName(TEnergyScaleFactorsFitModel_t fmodel);
  static std::string BEDoubleCaseName(TEnergyScaleFactorsBEDoubleCase_t beCase);

  std::string WorkCaseName() const { return ElectronEnergyScaleAdv_t::WorkCaseName(this->WorkCaseKind()); }
  std::string WorkCaseShortName() const { return ElectronEnergyScaleAdv_t::WorkCaseShortName(this->WorkCaseKind()); }
  std::string FitModelName() const { return ElectronEnergyScaleAdv_t::FitModelName(this->FitModelKind()); }
  std::string FitModelShortName() const { return ElectronEnergyScaleAdv_t::FitModelShortName(this->FitModelKind()); }

  const double* Smearings() const { return FSmearings; }
  double *EditSmearings() { return FSmearings; }
  double Smearing(int i) const { if (ValidSmFIdx(i)) return FSmearings[i]; else return 0; }
  double SmearingFactor(int i) const { if (ValidSmFIdx(i)) return FSmearings[i]; else return 0.; }
  
  const double* ScaleFactorErr() const { return FScaleFactorErrLo; }
  double *EditScaleFactorErr() { return FScaleFactorErrLo; }
  double ScaleFactorErr(int i) const { return (ValidScFIdx(i)) ? FScaleFactorErrLo[i] : 9999; }
  const double* ScaleFactorErrLo() const { return FScaleFactorErrLo; }
  double *EditScaleFactorErrLo() { return FScaleFactorErrLo; }
  double ScaleFactorErrLo(int i) const { return (ValidScFIdx(i)) ? FScaleFactorErrLo[i] : 9999; }
  const double* ScaleFactorErrHi() const { return FScaleFactorErrHi; }
  double *EditScaleFactorErrHi() { return FScaleFactorErrHi; }
  double ScaleFactorErrHi(int i) const { return (ValidScFIdx(i)) ? FScaleFactorErrHi[i] : 9999; }
  const double* SmearingErr() const { return FSmearingErrLo; }
  double *EditSmearingErr() { return FSmearingErrLo; }
  double SmearingErr(int i) const { return (ValidSmFIdx(i)) ? FSmearingErrLo[i] : 9999; }
  double SmearingFactorErr(int i) const { return (ValidSmFIdx(i)) ? FSmearingErrLo[i] : 9999; }
  double SmearFactorErr(int i) const { return (ValidSmFIdx(i)) ? FSmearingErrLo[i] : 9999; }
  const double* SmearingErrLo() const { return FSmearingErrLo; }
  double *EditSmearingErrLo() { return FSmearingErrLo; }
  double SmearingErrLo(int i) const { return (ValidSmFIdx(i)) ? FSmearingErrLo[i] : 9999; }
  const double* SmearingErrHi() const { return FSmearingErrHi; }
  double *EditSmearingErrHi() { return FSmearingErrHi; }
  double SmearingErrHi(int i) const { return (ValidSmFIdx(i)) ? FSmearingErrHi[i] : 9999; }

  int ChangeToNegs(const ElectronEnergyScaleAdv_t &sf);

  int MultiplyScaleFactor(int i, double multiplier) {
    if (ValidScFIdx(i)) FScaleFactors[i] *= multiplier; else return 0;
    return 1;
  }
  void MultiplyScaleFactor(double multiplier) {
    int imin=0; 
    int imax=FEtaDivisionCount;
    //if (this->FitModelKind()==_ESFModel_voigtian) {
    //  imin=FEtaDivisionCount; imax=2*FEtaDivisionCount;
    //}
    for (int i=imin; i<imax; ++i) FScaleFactors[i]*=multiplier;
  }

  int MultiplySmearFactor(int i, double multiplier) {
    if (ValidSmFIdx(i)) FSmearings[i] *= multiplier; else return 0;
    return 1;
  }
  void MultiplySmearFactor(double multiplier) {
    int imin=0; 
    int imax=FEtaDivisionCount;
    if (this->FitModelKind()==_ESFModel_voigtian) {
      imin=FEtaDivisionCount; imax=2*FEtaDivisionCount;
    }
    for (int i=imin; i<imax; ++i) FSmearings[i]*=multiplier;
  }

  mithep::myHistoClass_t *CreateHisto(TString namebase, double xmin=-3, double xmax=3, int plot_scale=1, int use_negative_eta_values=0) const;

  mithep::myHistoClass_t *CreateHisto_OnlyScale(TString namebase, double xmin=-3, double xmax=3, int use_negative_eta_values=0) const {
    return CreateHisto(namebase,xmin,xmax,1,use_negative_eta_values);
  }
  mithep::myHistoClass_t *CreateHisto_OnlySmear(TString namebase, double xmin=-3, double xmax=3, int use_negative_eta_values=0) const {
    return CreateHisto(namebase,xmin,xmax,0,use_negative_eta_values);
  }

  int Assign(const ElectronEnergyScaleAdv_t &sf);
  int Assign(const RooRealVar &mc_to_data, const std::vector<RooRealVar*> &scaling, int scaling_formula, const std::vector<RooRealVar*> &smearing);
  int Assign(TEnergyScaleFactorsFitModel_t model, const RooRealVar &mc_to_data,  const std::vector<RooRealVar*> &scaling, int scaling_formula, const std::vector<RooRealVar*> &smearing, const std::vector<RooRealVar*> *sm_pars1=NULL, const std::vector<RooRealVar*> *sm_pars2=NULL, const std::vector<RooRealVar*> *sc_pars3=NULL);

  int Load(const char *fname);
  int Load(const std::string &fname) { return Load(fname.c_str()); }
  int SaveASCIIFile(const char *filename) const; 
  int SaveASCIIFile(const char *mainfilename, const char *options) const;
  int LoadASCIIFile(const char *filename) { return DetermineESFWorkCase(filename,*this); }
  int SaveRootFile(const char *filename) const; // not ready
  int LoadRootFile(const char *filename); // not ready
  int LoadParameterSet(const char *eta_ranges_name, const char *fit_model_name, const char *file_name_base="dir-ESFPars20110930/testESF");

#ifdef UsePlotsAcc
  int ScaleValues(std::vector<mithep::AccEffData_t*> &data) const; // scale Proper(si)*Proper(sj)
  int ScaleValues(std::vector<RooDataSet*> &data, RooRealVar &x, const char *varName="x") const; // scale si*sj
  int InvScaleValues(std::vector<mithep::AccEffData_t*> &data) const; // scale 1/(si*sj)
#endif
  int ScaleValues(std::vector<std::vector<double>*> &data) const; // scale Proper(si)*Proper(sj)
  int InvScaleValues(std::vector<std::vector<double>*> &data) const; // scale 1/(Proper(si)*Proper(sj))

#ifdef DoSmearing
  template<class HistoClass_t>
  int ScaleMass(double mass, int etaBin1, int etaBin2, double &scaled_mass, double *weight, HistoClass_t *h, double *appliedScaleFactor=NULL, int do_scale=1) const;
  template<class HistoClass_t>
  int ScaleMCMass(double mass, int etaBin1, int etaBin2, double &scaled_mass, double *weight, HistoClass_t *h, double *appliedScaleFactor=NULL) const;
  template<class HistoClass_t>
  int SmearMass(double mass, int etaBin1, int etaBin2, int massBinCount, const double *massBins, double *result_distribution, double *weight, HistoClass_t *h, int do_smear=1) const;

  template<class HistoClass_t>
  int ScaleMass(double mass, double eta1, double eta2, double &scaled_mass, double *weight, HistoClass_t *h, double *appliedScaleFactor=NULL) const {
    int etaBin1,etaBin2;
    int res=this->PrepareEtaBinIdx(eta1,eta2,etaBin1,etaBin2);
    if (res) res=this->ScaleMass(mass,etaBin1,etaBin2,scaled_mass,weight,h,appliedScaleFactor);
    if (!res) std::cout << "error in ElectronEnergyScaleAdv_t::ScaleMass(dbl,dbl,dbl)\n";
    return res;
  }

  template<class HistoClass_t>
  int ScaleMCMass(double mass, double eta1, double eta2, double &scaled_mass, double *weight, HistoClass_t *h, double *appliedScaleFactor=NULL) const {
    int etaBin1,etaBin2;
    int res=this->PrepareEtaBinIdx(eta1,eta2,etaBin1,etaBin2);
    if (res) res=this->ScaleMCMass(mass,etaBin1,etaBin2,scaled_mass,weight,h,appliedScaleFactor);
    if (!res) std::cout << "error in ElectronEnergyScaleAdv_t::ScaleMCMass(dbl,dbl,dbl)\n";
    return res;
  }

  template<class HistoClass_t>
  int SmearMass(double mass, double eta1, double eta2, int massBinCount, const double *massBins, double *result_distribution, double *weight, HistoClass_t *h) const {
    int etaBin1,etaBin2;
    int res=this->PrepareEtaBinIdx(eta1,eta2,etaBin1,etaBin2);
    if (res) res=this->SmearMass(mass,etaBin1,etaBin2,massBinCount,massBins,result_distribution,weight,h);
    if (!res) std::cout << "error in ElectronEnergyScaleAdv_t::SmearMass(dbl,dbl,dbl)\n";
    return res;
  }


  template<class HistoClass_t>
  int ScaleAndSmearMass(double mass, double eta1, double eta2, int binCount, const double *massBins, double *result_distribution, double *weight, HistoClass_t *h, double *appliedScaleFactor=NULL) const;
  template<class HistoClass_t>
  int ScaleAndSmearMCMass(double mass, double eta1, double eta2, int binCount, const double *massBins, double *result_distribution, double *weight, HistoClass_t *h, double *appliedScaleFactor=NULL) const;
  
#endif

#ifdef DoSmearing
#ifdef UsePlotsAcc
  template<class HistoClass_t>
  int PrepareScaledSmearedData(const std::vector<mithep::AccEffData_t*> &dataMC, const std::vector<mithep::AccEffData_t*> &dataExp, int massBinCount, const double *massBins, std::vector<HistoClass_t*> &smearedMCV, std::vector<HistoClass_t*> &scaledExpV, const TString &extra, int do_smear=1, int do_scale=1);
#endif
#endif

  template<class HistoClass_t, class doubleVPtr_t>
  int PrepareScaledSmearedData(const std::vector<doubleVPtr_t> &dataMC, const std::vector<doubleVPtr_t> &dataExp, int massBinCount, const double *massBins, std::vector<HistoClass_t*> &smearedMCV, std::vector<HistoClass_t*> &scaledExpV, const TString &extra, int do_smear=1, int do_scale=1);
  
  
  int Verify_EnoughScalingVars(const char *calling_function_name, int needs_positions, const std::vector<RooRealVar*> &vars, const std::vector<RooRealVar*> *scPars2) const;
  int Verify_EnoughSmearingVars(const char *calling_function_name, int needs_positions, const std::vector<RooRealVar*> &vars, const std::vector<RooRealVar*> *smPars2, const std::vector<RooRealVar*> *smPars3) const;

  int CopyScalingFactors(std::vector<RooRealVar*> &vars, std::vector<RooRealVar*> *scPars2) const; // N
  int TransferScalingFactors(std::vector<RooRealVar*> &vars, std::vector<RooRealVar*> *scPars2) const; // N*(N+1)/2
  int CopySmearingFactors(std::vector<RooRealVar*> &vars, std::vector<RooRealVar*> *smPars2, std::vector<RooRealVar*> *smPars3) const; // N
  int TransferSmearingFactors(std::vector<RooRealVar*> &vars, std::vector<RooRealVar*> *smPars2, std::vector<RooRealVar*> *smPars3) const;  // N*(N+1)/2
  //int CopyValues(std::vector<RooRealVar*> &scale, std::vector<RooRealVar*> &smear) const;

  int PrintValuesHTMLTable(std::ostream &out, int append=0, const std::vector<const ElectronEnergyScaleAdv_t*> *others=NULL, const std::vector<std::string> *names=NULL) const;
  int PrintValuesHTMLTable(const TString *fout_name=NULL, const std::vector<const ElectronEnergyScaleAdv_t*> *others=NULL, const std::vector<std::string> *names=NULL) const;
  int PrintValuesHTMLTable(std::vector<std::string> &html_lines, const std::vector<const ElectronEnergyScaleAdv_t*> *others=NULL, const std::vector<std::string> *names=NULL) const;
  template<class Class_t>
    int PrintValuesHTMLTable(const Class_t *var, const std::vector<ElectronEnergyScaleAdv_t*> *others, const std::vector<std::string> *names=NULL) const;

  void PrepareHRNameEndings(std::vector<TString> &hrEnds) const;

#ifdef EtaEtaMass_H
  int ProcessEEMFile(const char *mc_file_name, const char *data_file_name, std::vector<std::vector<double>*> &mcData, std::vector<std::vector<double>*> &expData);
  int ProcessEEMFileApproximateMCWeight(const char *mc_file_name, const char *data_file_name, std::vector<std::vector<double>*> &mcData, std::vector<std::vector<double>*> &expData, double unitWeight=0.);
    //int ProcessEEMFileApproximateMCWeights(const char *mc_file_name, const char *data_file_name, std::vector<std::vector<double>*> &mcData, std::vector<std::vector<double>*> &expData);
  //int ProcessEEMFile(const char *mc_file_name, const char *data_file_name, std::vector<std::vector<double>*> &mcData, std::vector<std::vector<double>*> &mcWeights, std::vector<std::vector<double>*> &expData);
#endif

  void Print() const { std::cout << (*this); }

  friend std::ostream& operator<<(std::ostream &out, const ElectronEnergyScaleAdv_t &sf) {
    out << "EnergyScaleFactors{ WorkCase" << sf.WorkCase();
    out << '(' << ElectronEnergyScaleAdv_t::WorkCaseName(sf.WorkCaseKind());
    out << "), FitModel=" << sf.FitModel();
    out << '(' << ElectronEnergyScaleAdv_t::FitModelName(sf.FitModelKind());
    out << "), scaleFactorsCount=" << sf.ScFCount() << ", smearFactorsCount=" << sf.SmFCount();
    out << ", MCtoData=" << sf.MCtoData() << ", EtaDivisionCount=" << sf.EtaDivisionCount() << ", InvertedSF=" << sf.InvertedSF() << ": "; 
    for (int i=0; i<=sf.EtaDivisionCount(); ++i) {
      out << " " << sf.EtaBinLimit(i);
    }
    out << "; scalings: ";
    for (int i=0; i<sf.ScFCount(); ++i) {
      out << " " << sf.ScaleFactor(i);
    }
    out << "; smearings: ";
    for (int i=0; i<sf.SmFCount(); ++i) {
      out << " " << sf.Smearing(i);
    }
    out << "}";
    return out;
  }
};

// -----------------------------------------------------------------------
// -----------------------------------------------------------------------

void PrepareCombinationInfos(TEnergyScaleFactorsWorkCase_t workCase, TEnergyScaleFactorsBEDoubleCase_t theCase, std::vector<std::vector<int>*> &combine, std::vector<TString> &combination_names);

inline void PrepareCombinationInfos(const ElectronEnergyScaleAdv_t &sf, TEnergyScaleFactorsBEDoubleCase_t theCase, std::vector<std::vector<int>*> &combine, std::vector<TString> &combination_names) {
  PrepareCombinationInfos(sf.WorkCaseKind(),theCase,combine,combination_names);
}

// -----------------------------------------------------------------------
// -----------------------------------------------------------------------

#ifdef DoSmearing
template<class HistoClass_t>
int ElectronEnergyScaleAdv_t::ScaleMass(double mass, int etaBin1, int etaBin2, double &scaled_mass, double *weight, HistoClass_t *h, double *appliedScaleFactor, int do_scale) const {
  double factor=(do_scale) ? this->ProperScaleFactorValue(etaBin1)*this->ProperScaleFactorValue(etaBin2) : 1.0;
  scaled_mass=mass*factor;
  if (appliedScaleFactor) *appliedScaleFactor=factor;
  if (weight && h) {
    h->Fill(scaled_mass,*weight);
  }
  return 1;
}
#endif

// -----------------------------------------------------------------------

#ifdef DoSmearing
template<class HistoClass_t>
int ElectronEnergyScaleAdv_t::ScaleMCMass(double mass, int etaBin1, int etaBin2, double &scaled_mass, double *weight, HistoClass_t *h, double *appliedScaleFactor) const {
  double factor=1/(this->ProperScaleFactorValue(etaBin1)*this->ProperScaleFactorValue(etaBin2));
  if (appliedScaleFactor) *appliedScaleFactor=factor;
  scaled_mass=mass*factor;
  if (weight && h) {
    h->Fill(scaled_mass,*weight);
  }
  return 1;
}
#endif

// -----------------------------------------------------------------------

#ifdef DoSmearing
template<class HistoClass_t>
int ElectronEnergyScaleAdv_t::SmearMass(double mass, int etaBin1, int etaBin2, int massBinCount, const double *massBins, double *result_distribution, double *weight, HistoClass_t *h, int do_smear) const {
  //std::cout << "entered SmearMass" << std::endl;
  assert(massBins); assert(result_distribution);
  if (do_smear) {
    double s1=this->SmearingFactor(etaBin1); 
    double s2=this->SmearingFactor(etaBin2);
    double sigma=sqrt(s1*s1+s2*s2);
    //RooGaussian g("gauss","gauss",mass,sigma);
    double cdf_lower=ROOT::Math::gaussian_cdf(massBins[0],sigma,mass);
    double cdf_upper=0;
    for (int i=0; i<massBinCount; ++i) {
      cdf_upper=ROOT::Math::gaussian_cdf(massBins[i+1],sigma,mass);
      result_distribution[i]=cdf_upper-cdf_lower;
      cdf_lower=cdf_upper;
    }
  }
  else {
    for (int i=0; i<massBinCount; ++i) {
      //std::cout << "i=" << i << "/" << massBinCount << std::endl;
      if ((mass>=massBins[i]) && (mass<massBins[i+1])) {
	result_distribution[i]=1.;
      }
      else result_distribution[i]=0;
    }
    //std::cout << "test passed" << std::endl;
  }
  if (weight && h) {
    for (int i=0; i<massBinCount; ++i) {
      double mcenter=0.5*(massBins[i]+massBins[i+1]);
      h->Fill(mcenter,(*weight)*result_distribution[i]);
    }
  }
  //std::cout << "leaving SmearMass" << std::endl;
  return 1;
}
#endif

// -----------------------------------------------------------------------

#ifdef DoSmearing
template<class HistoClass_t>
int ElectronEnergyScaleAdv_t::ScaleAndSmearMass(double mass, double eta1, double eta2, int massBinCount, const double *massBins, double *result_distribution, double *weight, HistoClass_t *h, double *appliedScaleFactor) const {
  assert(massBins); assert(result_distribution);
  double scaled_mass;
  int etaBin1,etaBin2;
  int res=this->PrepareEtaBinIdx(eta1,eta2,etaBin1,etaBin2);
  if (res) res=this->ScaleMass<HistoClass_t>(mass,etaBin1,etaBin2,scaled_mass,NULL,NULL,appliedScaleFactor);
  if (res) res=this->SmearMass(scaled_mass,etaBin1,etaBin2,massBinCount,massBins,result_distribution,weight,h);
  if (!res) std::cout << "error in ElectronEnergyScaleAdv_t::ScaleAndSmearMass\n";
  return res;
}
#endif
// -----------------------------------------------------------------------

#ifdef DoSmearing
template<class HistoClass_t>
int ElectronEnergyScaleAdv_t::ScaleAndSmearMCMass(double mass, double eta1, double eta2, int massBinCount, const double *massBins, double *result_distribution, double *weight, HistoClass_t *h, double *appliedScaleFactor) const {
  assert(massBins); assert(result_distribution);
  double scaled_mass;
  int etaBin1,etaBin2;
  int res=this->PrepareEtaBinIdx(eta1,eta2,etaBin1,etaBin2);
  if (res) res=this->ScaleMCMass<HistoClass_t>(mass,etaBin1,etaBin2,scaled_mass,NULL,NULL,appliedScaleFactor);
  if (res) res=this->SmearMass(scaled_mass,etaBin1,etaBin2,massBinCount,massBins,result_distribution,weight,h);
  if (!res) std::cout << "error in ElectronEnergyScaleAdv_t::ScaleAndSmearMCMass\n";
  return res;
}
#endif

// -----------------------------------------------------------------------

#ifdef DoSmearing
#ifdef UsePlotsAcc
template<class HistoClass_t>
int ElectronEnergyScaleAdv_t::PrepareScaledSmearedData(const std::vector<mithep::AccEffData_t*> &dataMC, const std::vector<mithep::AccEffData_t*> &dataExp, int massBinCount, const double *massBins, std::vector<HistoClass_t*> &smearedMCV, std::vector<HistoClass_t*> &scaledExpV, const TString &extra, int do_smear, int do_scale) {
  const char *fncname="PrepareScaledSmearedData";
  const int expect_size=this->FEtaDivisionCount*(this->FEtaDivisionCount+1)/2;
  if ((expect_size!=dataMC.size()) || (expect_size!=dataExp.size())) {
    std::cout << fncname << ": etaDivisionCount=" << this->FEtaDivisionCount << ", expect_size=" << expect_size << ", dataMC.size=" << dataMC.size() << ", dataExp.size=" << dataExp.size() << "\n";
    return 0;
  }
  // prepare histograms
  std::vector<TString> hrEnds;
  this->PrepareHRNameEndings(hrEnds);
  char buf[30];
  int k=0;
  for (int ib=0; ib<this->FEtaDivisionCount; ++ib) {
    for (int jb=ib; jb<this->FEtaDivisionCount; ++jb, ++k) {
      sprintf(buf,"smeared%sMC_%d_%d__%d",extra.Data(),ib,jb,k);
      HistoClass_t *hmcnew=new HistoClass_t(buf,hrEnds[k].Data(),massBinCount,massBins);
      hmcnew->GetXaxis()->SetTitle("mass");
      hmcnew->GetYaxis()->SetTitle("count");
      smearedMCV.push_back(hmcnew);
      //smearedMCV.back()->Set
      sprintf(buf,"scaled%sExp_%d_%d__%d",extra.Data(),ib,jb,k);
      HistoClass_t* hexpnew=new HistoClass_t(buf,hrEnds[k].Data(),massBinCount,massBins);
      hexpnew->GetXaxis()->SetTitle("mass");
      hexpnew->GetYaxis()->SetTitle("count");
      scaledExpV.push_back(hexpnew);
    }
  }
  //std::cout << "prepare smeared mc" << std::endl;
  // prepare smeared mc
  k=0;
  double result_distribution[massBinCount];
  double weight=1;
  for (int ib=0; ib<this->FEtaDivisionCount; ++ib) {
    for (int jb=ib; jb<this->FEtaDivisionCount; ++jb, ++k) {
      const std::vector<double>* masses=dataMC[k]->MassConstPtr();
      HistoClass_t *histo=smearedMCV[k];
      for (unsigned int i=0; i<masses->size(); ++i) {
	//std::cout << "smearing mass k=" << k << ", element=" << i << std::endl;
	this->SmearMass((*masses)[i],ib,jb,massBinCount,massBins,result_distribution,&weight,histo, do_smear);
      }
    }
  }
  //std::cout << "prepare scaled exp" << std::endl;
  // prepare scaled exp
  k=0;
  double scaled_mass=0;
  //double weight=1;
  for (int ib=0; ib<this->FEtaDivisionCount; ++ib) {
    for (int jb=ib; jb<this->FEtaDivisionCount; ++jb, ++k) {
      const std::vector<double>* masses=dataExp[k]->MassConstPtr();
      HistoClass_t *histo=scaledExpV[k];
      for (unsigned int i=0; i<masses->size(); ++i) {
	this->ScaleMass((*masses)[i],ib,jb,scaled_mass,&weight,histo,NULL,do_scale);
      }
    }
  }
  //std::cout << "done" << std::endl;
  return 1;
}
#endif
#endif

// -----------------------------------------------------------------------

template<class HistoClass_t, class doubleVPtr_t>
int ElectronEnergyScaleAdv_t::PrepareScaledSmearedData(const std::vector<doubleVPtr_t> &dataMC, const std::vector<doubleVPtr_t> &dataExp, int massBinCount, const double *massBins, std::vector<HistoClass_t*> &smearedMCV, std::vector<HistoClass_t*> &scaledExpV, const TString &extra, int do_smear, int do_scale) {
  const char *fncname="PrepareScaledSmearedData";
  const unsigned int expect_size=this->FEtaDivisionCount*(this->FEtaDivisionCount+1)/2;
  if ((expect_size!=dataMC.size()) || (expect_size!=dataExp.size())) {
    std::cout << fncname << ": etaDivisionCount=" << this->FEtaDivisionCount << ", expect_size=" << expect_size << ", dataMC.size=" << dataMC.size() << ", dataExp.size=" << dataExp.size() << "\n";
    return 0;
  }
  // prepare histograms
  std::vector<TString> hrEnds;
  this->PrepareHRNameEndings(hrEnds);
  char buf[30];
  int k=0;
  for (int ib=0; ib<this->FEtaDivisionCount; ++ib) {
    for (int jb=ib; jb<this->FEtaDivisionCount; ++jb, ++k) {
      sprintf(buf,"smeared%sMC_%d_%d__%d",extra.Data(),ib,jb,k);
      HistoClass_t *hmcnew=new HistoClass_t(buf,hrEnds[k].Data(),massBinCount,massBins);
      hmcnew->GetXaxis()->SetTitle("mass");
      hmcnew->GetYaxis()->SetTitle("count");
      smearedMCV.push_back(hmcnew);
      //smearedMCV.back()->Set
      sprintf(buf,"scaled%sExp_%d_%d__%d",extra.Data(),ib,jb,k);
      HistoClass_t* hexpnew=new HistoClass_t(buf,hrEnds[k].Data(),massBinCount,massBins);
      hexpnew->GetXaxis()->SetTitle("mass");
      hexpnew->GetYaxis()->SetTitle("count");
      scaledExpV.push_back(hexpnew);
    }
  }
  //std::cout << "prepare smeared mc" << std::endl;
  // prepare smeared mc
  k=0;
  double result_distribution[massBinCount];
  double weight=1;
  for (int ib=0; ib<this->FEtaDivisionCount; ++ib) {
    for (int jb=ib; jb<this->FEtaDivisionCount; ++jb, ++k) {
      const std::vector<double>* masses=(const std::vector<double>*)(dataMC[k]);
      HistoClass_t *histo=smearedMCV[k];
      for (unsigned int i=0; i<masses->size(); ++i) {
	//std::cout << "smearing mass k=" << k << ", element=" << i << std::endl;
	this->SmearMass((*masses)[i],ib,jb,massBinCount,massBins,result_distribution,&weight,histo, do_smear);
      }
    }
  }
  //std::cout << "prepare scaled exp" << std::endl;
  // prepare scaled exp
  k=0;
  double scaled_mass=0;
  //double weight=1;
  for (int ib=0; ib<this->FEtaDivisionCount; ++ib) {
    for (int jb=ib; jb<this->FEtaDivisionCount; ++jb, ++k) {
      const std::vector<double>* masses=(const std::vector<double>*)(dataExp[k]);
      HistoClass_t *histo=scaledExpV[k];
      for (unsigned int i=0; i<masses->size(); ++i) {
	this->ScaleMass((*masses)[i],ib,jb,scaled_mass,&weight,histo,NULL,do_scale);
      }
    }
  }
  //std::cout << "done" << std::endl;
  return 1;
}


// -----------------------------------------------------------------------

template<class Class_t>
int ElectronEnergyScaleAdv_t::PrintValuesHTMLTable(const Class_t *var, const std::vector<ElectronEnergyScaleAdv_t*> *others, const std::vector<std::string> *names) const {
  std::vector<const ElectronEnergyScaleAdv_t*> tmp;
  if (!others) {
    std::cout << "ElectronEnergyScaleAdv_t::PrintValuesHTMLTable<T> got null others\n";
    return 0;
  }
  for (unsigned int i=0; i<others->size(); ++i) {
    tmp.push_back((const ElectronEnergyScaleAdv_t*)(*others)[i]);
  }
  int res=this->PrintValuesHTMLTable(var,&tmp,names);
  return res;
}


// -----------------------------------------------------------------------
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
// -----------------------------------------------------------------------

#endif
