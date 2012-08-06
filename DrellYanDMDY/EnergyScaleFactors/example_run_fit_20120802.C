#define __myLib__
#include <TROOT.h>
#include <TSystem.h>
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


void example_fit_run(const char *workCaseStr, const char *fitModelStr)
{
  saveFigFileTag="fit_";

  std::vector<std::vector<double>*> massMC,massExp;
  std::vector<std::vector<double>*> Et1_mc,Et2_mc; // dummies
  std::vector<std::vector<double>*> Et1_exp,Et2_exp; // dummies

  TEnergyScaleFactorsWorkCase_t work_case=_ESFWC_2bin;
  work_case=IdentifyESFWorkCase(workCaseStr);
  if (work_case==_ESFWC_none) { 
    std::cout << "failed to determine work case\n"; 
    return;
  }
  TEnergyScaleFactorsFitModel_t fit_model=_ESFModel_cballShift;
  fit_model=_ESFModel_gauss;
  fit_model=IdentifyESFFitModel(fitModelStr,1);
  if (fit_model==_ESFModel_none) {
    std::cout << "failed to determine fit model\n";
    return;
  }

  std::string wcName=ElectronEnergyScaleAdv_t::WorkCaseShortName(work_case);
  std::string fmName=ElectronEnergyScaleAdv_t::FitModelShortName(fit_model);
  const int bufsize=100;
  char buf[bufsize];

  std::string eem_path="/media/spektras/DYee2011/2ndProd-ntuple/";
  std::string mcEEMFile =eem_path + std::string("zee_bpVer1_EtaEtaM.root");
  std::string expEEMFile=eem_path + std::string("data_bpVer1_EtaEtaM.root");

  sprintf(buf,"dir-Work-%s-%s-fit/",wcName.c_str(),fmName.c_str());
  //sprintf(buf,"dir-WorkFEWZ-%s-%s-fit/",wcName.c_str(),fmName.c_str());

  histoDir_USER=buf;
  std::cout << "histoDir_USER=" << histoDir_USER << "\n";
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

  //PrintVVecCounts("massMC ",massMC);
  //PrintVVecCounts("massExp",massExp);

  //return ;

  // Adjust integrator
  RooAbsReal::defaultIntegratorConfig()->setEpsAbs(1e-5) ;
  RooAbsReal::defaultIntegratorConfig()->setEpsRel(1e-5) ;
  RooAbsReal::defaultIntegratorConfig()->getConfigSection("RooIntegrator1D").setRealValue("maxSteps",50) ;
  //RooAbsReal::defaultIntegratorConfig()->Print("v") ;


  DYTools::PrepareScEtaBinInfo(work_case);
  DYTools::g_esfFitModel=fit_model;
  std::string foutName;
  foutName=histoDir_USER + std::string("fit.log");
  fpos_t pos;
  int fd=0;
  if (foutName.size()) {
    fgetpos(stdout, &pos);
    fd = dup(fileno(stdout));
    freopen(foutName.c_str(), "w", stdout);
    if (!fd) {
      std::cout << " failed to open the file <" << foutName << ">\n";
      return;
    }
  }
  int res=FitMassSmearingIIe(std::cout,DYTools::nScEtaBins3,massMC,Et1_mc,Et2_mc,massExp);
  if (foutName.size()) {
    std::cout << "fit returned code " << res << "\n";
    fflush(stdout);
    dup2(fd, fileno(stdout));
    close(fd);
    clearerr(stdout);
    fsetpos(stdout,&pos);
  }
  std::cout << "fit returned code " << res << "\n";
  return;
}


