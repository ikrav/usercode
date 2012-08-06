#define __myLib__
#include <TROOT.h>
//#include <TH1D.h>
//#include <vector>

//namespace mithep {
//  typedef std::vector<double> AccEffData_t;
//  typedef TH1D myHistoClass_t;
//};

//#include "ElectronEnergyScaleAdv.hh"
#include "fitMassScaleFactorsE.C"
#include "HelpingTools.cc"
#include "../Include/EtaEtaMass.hh"

void example_rec_run(const char *workCaseStr, const char *fitModelStr)
{
  saveFigFileTag="rec_";

  std::vector<std::vector<double>*> massMC,massExp;
  std::vector<std::vector<double>*> Et1_mc,Et2_mc; // dummies
  std::vector<std::vector<double>*> Et1_exp,Et2_exp; // dummies

  TEnergyScaleFactorsWorkCase_t work_case=_ESFWC_6bin;
  work_case=_ESFWC_2bin;
  work_case=IdentifyESFWorkCase(workCaseStr);
  if (work_case==_ESFWC_none) { 
    std::cout << "failed to determine work case\n"; 
    return;
  }
  TEnergyScaleFactorsFitModel_t fit_model=_ESFModel_gauss;
  fit_model=_ESFModel_breitWigner;
  fit_model=_ESFModel_voigtian;
  fit_model=IdentifyESFFitModel(fitModelStr,1);
  if (fit_model==_ESFModel_none) {
    std::cout << "failed to determine fit model\n";
    return;
  }

  std::string wcName=ElectronEnergyScaleAdv_t::WorkCaseShortName(work_case);
  std::string fmName=ElectronEnergyScaleAdv_t::FitModelShortName(fit_model);
  std::string taskNameD=wcName+std::string("-")+fmName;
  std::string taskNameU=wcName+std::string("_")+fmName;
  const int bufsize=100;
  char buf[bufsize];

  std::string eem_path="/media/spektras/DYee2011/2ndProd-ntuple/";
  std::string mcEEMFile =eem_path + std::string("zee_bpVer1_EtaEtaM.root");
  std::string expEEMFile=eem_path + std::string("data_bpVer1_EtaEtaM.root");

  sprintf(buf,"dir-Work-%s-%s-rec/",wcName.c_str(),fmName.c_str());
  histoDir_USER=buf;
  gSystem->mkdir(histoDir_USER.c_str(),true);

  ElectronEnergyScaleAdv_t esf;
  if (!esf.WorkCase(work_case,fit_model)) {
    return ;
  }
  if (!esf.ProcessEEMFile(mcEEMFile.c_str(),expEEMFile.c_str(),massMC,massExp)) 
    //if (!esf.ProcessEEMFileApproximateMCWeight(mcEEMFile.c_str(),expEEMFile.c_str(),massMC,massExp)) 
    {
    std::cout << "failed to load data\n";
    return;
  }

  PrintVVecCounts("massMC ",massMC);
  PrintVVecCounts("massExp",massExp);

  //return;

  sprintf(buf,"dir-Work-%s-fit/testESF_%s.inp",taskNameD.c_str(),taskNameU.c_str());
  if (!esf.LoadASCIIFile(buf)) {
    std::cout << "failed to load ESF from <" << buf << ">\n";
    return;
  }
  RooRealVar mass("x","mass",fit_range_min,fit_range_max); // call it "x"!
  ScaledSmearedDataAux_t aux(esf);
  std::vector<RooAddPdf*> mcModelV;
  std::vector<RooDataSet*> expModelV;
  mass.setBins(10000,"cache");

  int res=1;

  if (res) res=PrepareScaledSmearedData(esf,massMC,massExp,Et1_exp,Et2_exp,mass,aux,mcModelV,expModelV);

  std::vector<std::vector<int>*> combine;
  std::vector<TString> combination_names;
  std::vector<std::string> html_lines;   
  if (res) PrepareCombinationInfos(esf,_ESFWC_all, combine, combination_names);

  if (0 && res) { // create comparison plots (optional - many plots)
    res=CreateComparisonPlots(esf,massMC,massExp,html_lines,mass,aux,mcModelV,expModelV,combine,combination_names);
  }

  if (res) {
    int bin_count=esf.EtaDivisionCount();
    //int expect_count=bin_count*(bin_count+1)/2;
    std::vector<TString> idxV;
    std::vector<TString> nameEndings;
    PrepareNameEndings(bin_count,idxV,nameEndings);

    std::vector<TString> hr_nameEndings;
    if (!esf.PrintValuesHTMLTable(html_lines)) {
      std::cout << "failed in sf.PrintValuesHTMLTable\n";
      return;
    }
    esf.PrepareHRNameEndings(hr_nameEndings);

    res=makeHTML(histoDir_USER.c_str(),"0scanAndScale.html", "scanAndScale inspection plot", html_lines, hr_nameEndings,nameEndings, mass,expModelV,mcModelV, &combine, &combination_names);
  }
}


