#ifndef fitMassScaleFactorsE_H
#define fitMassScaleFactorsE_H

#include "MyFitModels.hh"
#include "ElectronEnergyScaleAdv.hh"

extern double massResAccuracy; // for fitting

// ------------------------------------------------------------------

// fit with a Gaussian model
//int FitMassSmearingIIc(int bin_count, const std::vector<RooDataSet*> &mc_mass, const std::vector<RooDataSet*> &mc_pt, const std::vector<RooDataSet*> &exp, RooRealVar &x, RooRealVar &y, const DataInfos_t &info);

//#ifndef __myLib__
// fit with a model
int FitMassSmearingIIe(std::ostream &out, unsigned int bin_count, const std::vector<RooDataSet*> &mc_mass, const std::vector<RooDataSet*> &mc_Et1, const std::vector<RooDataSet*> &mc_Et2, const std::vector<RooDataSet*> &exp, RooRealVar &x, RooRealVar &y, const DataInfos_t &info);

//#endif

inline 
int FitMassSmearingIIe(std::ostream &out, unsigned int bin_count, const std::vector<std::vector<double>*> &mc_mass, const std::vector<std::vector<double>*> &mc_Et1, const std::vector<std::vector<double>*> &mc_Et2, const std::vector<std::vector<double>*> &exp_mass) {
  RooRealVar x("x","x",0,2000);
  RooRealVar y("y","y",-1e10,1e10);
  std::vector<RooDataSet*> *dsMC,*dsExp;
  std::vector<RooDataSet*> dsMC_Et1,dsMC_Et2;
  if (0) { PrintVVecCounts("mc_Et1 ",mc_Et1); PrintVVecCounts("mc_Et2 ",mc_Et2); }
  out << "mc_mass.size()=" << mc_mass.size() << ", exp_mass.size()=" << exp_mass.size() << "\n";
  dsMC=ConvertToDataSet(mc_mass,x,y);
  out << "dsMC->size= " << dsMC->size() << "\n";
  dsExp=ConvertToDataSet(exp_mass,x,y);
  out << "chk dsMC->size= " << dsMC->size() << "\n";
  out << "dsExp->size= " << dsExp->size() << "\n";
  std::vector<TString> labelV;
  labelV.push_back(TString("MC")); labelV.push_back(TString("Exp"));
  std::vector<Int_t> mcolorsInt,lineV;
  mcolorsInt.push_back(kBlue+1); mcolorsInt.push_back(kRed+1);
  lineV.push_back(kSolid); lineV.push_back(kSolid);
  DataInfos_t info;
  info.Assign(&labelV,&mcolorsInt,&lineV);
  int fit_res=FitMassSmearingIIe(out,bin_count,*dsMC,dsMC_Et1,dsMC_Et2,*dsExp,x,y,info);
  if (dsMC) delete dsMC;
  if (dsExp) delete dsExp;
  return fit_res;
}
//#endif


// ------------------------------------------------------------------

namespace DYTools {
  extern TEnergyScaleFactorsWorkCase_t g_esfWorkCase;
  extern TEnergyScaleFactorsFitModel_t g_esfFitModel;
};

namespace DYTools {
  inline int PrepareScEtaBinInfo(TEnergyScaleFactorsWorkCase_t esfWorkCase) {
    DYTools::g_esfWorkCase=esfWorkCase;
    DYTools::scEtaRangeSet=(esfWorkCase<_ESFWC_BothSided) ? 3 : -3;  // put -3 if negative eta value bins are used
    int n=ElectronEnergyScaleAdv_t::GetEtaBinLimitCount(esfWorkCase);
    DYTools::nScEtaBins3=n;
    DYTools::scEtaBinLimits3 = new double[n+1];
    ElectronEnergyScaleAdv_t::FillEtaBinLimits(esfWorkCase, DYTools::scEtaBinLimits3);
    return 1;
  }
};


// ------------------------------------------------------------------

struct ScaledSmearedDataAux_t {
  int fitModel;
  std::vector<mithep::myHistoClass_t*> massH;
  std::vector<RooDataHist*> mcHistV;
  std::vector<RooHistPdf*> mcV;
  RooRealVar bias;
  std::vector<RooRealVar*> gWidthV;
  std::vector<RooGaussModel*> gaussV;
  std::vector<RooGaussModel*> gauss2V;
  std::vector<RooAddModel*> mixedV;
  std::vector<RooBifurGauss*> bifurGaussV;
  std::vector<RooRealVar*> smPars2V, smPars3V;
  std::vector<RooGExpModel*> gaussExpV;
  std::vector<RooCBShape*> cballV;
  std::vector<RooBreitWigner*> bwShapeV;
  std::vector<RooVoigtian*> voigtianV;
  std::vector<RooFFTConvPdf*> mcSignalV;
  std::vector<RooExtendPdf*> mcESigV;
  // cbShift
  std::vector<RooMoment*> calcMeanV;
  std::vector<RooFormulaVar*> shiftFromMeanV;
  std::vector<RooCBShape*> cballShiftedV;

  RooRealVar MCOverData;

  ScaledSmearedDataAux_t(const ElectronEnergyScaleAdv_t &sf) : 
    fitModel(sf.FitModel()),
    bias("bias","bias",0),
    MCOverData("MCOverData","MC_over_data",sf.MCtoData())
  {}
};

// ------------------------------------------------------------------

#ifndef __myLib__
int PrepareScaledSmearedData(const ElectronEnergyScaleAdv_t &sf, const std::vector<std::vector<double>*> &massMC, const std::vector<std::vector<double>*> &massExp, const std::vector<std::vector<double>*> &Et1_exp, const std::vector<std::vector<double>*> &Et2_exp, RooRealVar &mass, ScaledSmearedDataAux_t &aux, std::vector<RooAddPdf*> &mcModelV, std::vector<RooDataSet*> &expModelV, int scaleMC=0);


inline int PrepareScaledSmearedData(const ElectronEnergyScaleAdv_t &sf, const std::vector<mithep::AccEffData_t*> &dataMC, const std::vector<mithep::AccEffData_t*> &dataExp, const std::vector<mithep::AccEffData_t*> &dataExp_Et2, RooRealVar &mass, ScaledSmearedDataAux_t &aux, std::vector<RooAddPdf*> &mcModelV, std::vector<RooDataSet*> &expModelV, int scaleMC=0) {
  std::vector<std::vector<double>*> massMC,massExp,Et1_exp,Et2_exp;
  massMC.reserve(dataMC.size()); 
  massExp.reserve(dataExp.size()); Et1_exp.reserve(dataExp.size());
  for (unsigned int i=0; i<dataMC.size(); ++i) {
    massMC.push_back(dataMC[i]->GetClonePtr(0));
  }
  for (unsigned int i=0; i<dataExp.size(); ++i) {
    massExp.push_back(dataExp[i]->GetClonePtr(0));
  }
  //for (unsigned int i=0; i<dataExp.size(); ++i) {
  //Et1_exp.push_back(dataExp[i]->GetClonePtr(1));
  //}
  if (0) {
    std::cout << "dumping the data\n";
    FILE *fout=fopen("fit_test_masses_2bins_20111014.dat","w");
    fprintf(fout,"! file created for test purposes\n");
    fprintf(fout,"! MC data\n");
    fprintf(fout,"%u ! MC mass vec.size\n",massMC.size());
    for (unsigned int i=0; i<massMC.size(); ++i) {
      fprintf(fout,"  %u \n", massMC[i]->size());
      for (unsigned int ii=0; ii<massMC[i]->size(); ++ii) {
	fprintf(fout," %6.4lf",(*massMC[i])[ii]);
      }
      fprintf(fout,"\n");
    }
    fprintf(fout,"! Exp data\n");
    fprintf(fout,"%u ! Exp mass vec.size\n",massExp.size());
    for (unsigned int i=0; i<massExp.size(); ++i) {
      fprintf(fout,"  %u \n", massExp[i]->size());
      for (unsigned int ii=0; ii<massExp[i]->size(); ++ii) {
	fprintf(fout," %6.4lf",(*massExp[i])[ii]);
      }
      fprintf(fout,"\n");
    }
    fclose(fout);
    return 0;
  }
  return PrepareScaledSmearedData(sf,massMC,massExp,Et1_exp,Et2_exp,mass,aux,mcModelV,expModelV,scaleMC);
}
#endif

// ------------------------------------------------------------------

//int DetermineESFWorkCase(const char *filename);
//int DetermineESFFitModel(const char *filename);
//int DetermineESFWorkCaseAndFitModel(const char *filename, int &work_case, int &fit_model);

// ------------------------------------------------------------------


#endif
