#include "ElectronEnergyScaleAdv.hh"
#include <assert.h>
#include <TVectorT.h>
#include <TVector.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <algorithm>
#include <TCanvas.h>
#include "HelpingTools.hh"


inline int PosOk(size_t n) { return (n!=std::string::npos) ? 1:0; }

// ------------------------------------------------------------
// ------------------------------------------------------------

TEnergyScaleFactorsWorkCase_t IdentifyESFWorkCase(const char *s) {
  TEnergyScaleFactorsWorkCase_t wc=_ESFWC_none;
  if (strstr(s,"2binNegs")) { wc=_ESFWC_2binNegs; }
  else if (strstr(s,"3binNegs")) { wc=_ESFWC_3binNegs; }
  else if (strstr(s,"4binNegs")) { wc=_ESFWC_4binNegs; }
  else if (strstr(s,"5binNegs")) { wc=_ESFWC_5binNegs; }
  else if (strstr(s,"6binNegs")) { wc=_ESFWC_6binNegs; }
  else if (strstr(s,"5EBNegs")) wc=_ESFWC_5EBNegs;
  else if (strstr(s,"6EBNegs")) wc=_ESFWC_6EBNegs;
  else if (strstr(s,"4EENegs")) wc=_ESFWC_4EENegs;
  else if (strstr(s,"3EB3EEbinNegs")) { wc=_ESFWC_3EB3EENegs; }
  else if (strstr(s,"3EB3EENegs")) { wc=_ESFWC_3EB3EENegs; }
  else if (strstr(s,"4EB3EEbinNegs")) { wc=_ESFWC_4EB3EENegs; }
  else if (strstr(s,"4EB3EENegs")) { wc=_ESFWC_4EB3EENegs; }
  else if (strstr(s,"5EB3EEbinNegs")) { wc=_ESFWC_5EB3EENegs; }
  else if (strstr(s,"5EB3EENegs")) { wc=_ESFWC_5EB3EENegs; }
  else if (strstr(s,"6EB3EEbinNegs")) { wc=_ESFWC_6EB3EENegs; }
  else if (strstr(s,"6EB3EENegs")) { wc=_ESFWC_6EB3EENegs; }
  else if (strstr(s,"7EB3EEbinNegs")) { wc=_ESFWC_7EB3EENegs; }
  else if (strstr(s,"7EB3EENegs")) { wc=_ESFWC_7EB3EENegs; }
  else if (strstr(s,"8EB3EEbinNegs")) { wc=_ESFWC_8EB3EENegs; }
  else if (strstr(s,"8EB3EENegs")) { wc=_ESFWC_8EB3EENegs; }
  else if (strstr(s,"3EB3EEa")) { wc=_ESFWC_3EB3EEa; }
  else if (strstr(s,"3EB3EEabins")) { wc=_ESFWC_3EB3EEa; }
  else if (strstr(s,"3EB3EE")) { wc=_ESFWC_3EB3EE; }
  else if (strstr(s,"3EB3EEbins")) { wc=_ESFWC_3EB3EE; }
  else if (strstr(s,"4EB3EE")) { wc=_ESFWC_4EB3EE; }
  else if (strstr(s,"4EB3EEbins")) { wc=_ESFWC_4EB3EE; }
  else if (strstr(s,"2bin")) { wc=_ESFWC_2bin; }
  else if (strstr(s,"3bin")) { wc=_ESFWC_3bin; }
  else if (strstr(s,"4bin")) { wc=_ESFWC_4bin; }
  else if (strstr(s,"5bin")) { wc=_ESFWC_5bin; }
  else if (strstr(s,"6bin")) { wc=_ESFWC_6bin; }
  else {
    std::cout << "IdentifyWorkCase(" << s << "): failed to recognize the work case\n";
   }
   return wc;
 }

// ------------------------------------------------------------

std::string ElectronEnergyScaleAdv_t::WorkCaseName(TEnergyScaleFactorsWorkCase_t wcase) {
  std::string name;
  switch(wcase) {
  case _ESFWC_1bin: name="1 bin"; break;
  case _ESFWC_2bin: name="2 bins"; break;
  case _ESFWC_3bin: name="3 bins"; break;
  case _ESFWC_4bin: name="4 bins"; break;
  case _ESFWC_5bin: name="5 bins"; break;
  case _ESFWC_6bin: name="6 bins"; break;
  case _ESFWC_3EB3EE: name="3 EB and 3 EE bins"; break;
  case _ESFWC_3EB3EEa: name="3 EB and 3 EE bins (ver.a)"; break;
  case _ESFWC_4EB3EE: name="4 EB and 3 EE bins"; break;
  case _ESFWC_2binNegs: name="2 bins on each eta side"; break;
  case _ESFWC_3binNegs: name="3 bins on each eta side"; break;
  case _ESFWC_4binNegs: name="4 bins on each eta side"; break;
  case _ESFWC_5binNegs: name="5 bins on each eta side"; break;
  case _ESFWC_6binNegs: name="6 bins on each eta side"; break;
  case _ESFWC_3EB3EENegs: name="3 EB and 3 EE bins for each eta side"; break;
  case _ESFWC_4EB3EENegs: name="4 EB and 3 EE bins for each eta side"; break;
  case _ESFWC_5EB3EENegs: name="5 EB and 3 EE bins for each eta side"; break;
  case _ESFWC_6EB3EENegs: name="6 EB and 3 EE bins for each eta side"; break;
  case _ESFWC_7EB3EENegs: name="7 EB and 3 EE bins for each eta side"; break;
  case _ESFWC_8EB3EENegs: name="8 EB and 3 EE bins for each eta side"; break;
  case _ESFWC_5EBNegs: name="5 EB bins for each eta side"; break;
  case _ESFWC_6EBNegs: name="6 EB bins for each eta side"; break;
  case _ESFWC_4EENegs: name="4 EE bins for each eta side"; break;
  default:
    std::cout << "ElectronEnergyScaleAdv_t:: WorkCaseName is not ready for wcase=" << wcase << "\n";
    char buf[30];
    sprintf(buf,"unknown wcase=%d",wcase);
    name=buf;
  }
  return name;
}
// ------------------------------------------------------------

std::string ElectronEnergyScaleAdv_t::WorkCaseShortName(TEnergyScaleFactorsWorkCase_t wcase) {
  std::string name;
  switch(wcase) {
  case _ESFWC_1bin: name="1bin"; break;
  case _ESFWC_2bin: name="2bins"; break;
  case _ESFWC_3bin: name="3bins"; break;
  case _ESFWC_4bin: name="4bins"; break;
  case _ESFWC_5bin: name="5bins"; break;
  case _ESFWC_6bin: name="6bins"; break;
  case _ESFWC_3EB3EE: name="3EB3EEbins"; break;
  case _ESFWC_3EB3EEa: name="3EB3EEabins"; break;
  case _ESFWC_4EB3EE: name="4EB3EEbins"; break;
  case _ESFWC_2binNegs: name="2binNegs"; break;
  case _ESFWC_3binNegs: name="3binNegs"; break;
  case _ESFWC_4binNegs: name="4binNegs"; break;
  case _ESFWC_5binNegs: name="5binNegs"; break;
  case _ESFWC_6binNegs: name="6binNegs"; break;
  case _ESFWC_3EB3EENegs: name="3EB3EENegs"; break;
  case _ESFWC_4EB3EENegs: name="4EB3EENegs"; break;
  case _ESFWC_5EB3EENegs: name="5EB3EENegs"; break;
  case _ESFWC_6EB3EENegs: name="6EB3EENegs"; break;
  case _ESFWC_7EB3EENegs: name="7EB3EENegs"; break;
  case _ESFWC_8EB3EENegs: name="8EB3EENegs"; break;
  case _ESFWC_5EBNegs: name="5EBNegs"; break;
  case _ESFWC_6EBNegs: name="6EBNegs"; break;
  case _ESFWC_4EENegs: name="4EENegs"; break;
  default:
    std::cout << "ElectronEnergyScaleAdv_t:: WorkCaseShortName is not ready for wcase=" << wcase << "\n";
    char buf[30];
    sprintf(buf,"unknown wcase=%d",wcase);
    name=buf;
  }
  return name;
}

 // ------------------------------------------------------------

 TEnergyScaleFactorsFitModel_t IdentifyESFFitModel(const char *sOrig, int verb) {
   TEnergyScaleFactorsFitModel_t fm=_ESFModel_none;
   if (!strstr(sOrig,"fit model")) {
     if (verb) std::cout << "failed to locate \"fit model\" in <" << sOrig << ">\n";
     return fm;
   }
   if (strlen(sOrig)>299) return fm;
   char s[300];
   sprintf(s,sOrig);
   for (unsigned int i=0; i<strlen(sOrig); ++i) {
     if ((s[i]>='A') && (s[i]<='Z')) s[i]=char(int(s[i])-int('A')+'a');
   }
   if (verb) std::cout << "studying <" << s << ">\n";
   if (strstr(s,"fit model none")) { fm=_ESFModel_none; }
   else if (strstr(s,"gaussexp") || strstr(s,"gaussExp") || strstr(s,"GaussExp")) { fm=_ESFModel_gaussExp; }
   else if (strstr(s,"doublegauss") || strstr(s,"doubleGauss") || strstr(s,"DoubleGauss")) { fm=_ESFModel_doubleGauss; }
   else if (strstr(s,"bifurgauss") || strstr(s,"bifurGauss") || strstr(s,"BifurGauss")) { fm=_ESFModel_bifurGauss; }
   else if (strstr(s,"gausshift")) fm=_ESFModel_gaussShift;
   else if (strstr(s,"gauss") || strstr(s,"Gauss")) { fm=_ESFModel_gauss; }
   else if (strstr(s,"cballshift")) { fm=_ESFModel_cballShift; }
   else if (strstr(s,"cball") || strstr(s,"CBall") || strstr(s,"crystal ball")) { fm=_ESFModel_cball; }
   else if (strstr(s,"breitwigner")) fm=_ESFModel_breitWigner;
   else if (strstr(s,"voigtian")) fm=_ESFModel_voigtian;
   else {
     std::cout << "IdentifyESFFitModel(" << s << "): failed to recognize the fit model\n";
   }
   if (verb) std::cout << "identified model=" << fm << "\n";
   return fm;
 }

 // ------------------------------------------------------------
 // ------------------------------------------------------------

 int ElectronEnergyScaleAdv_t::GetEtaBinLimitCount(TEnergyScaleFactorsWorkCase_t theCase) {
   int count=0;
   //std::cout << "ElectronEnergyScaleAdv_t::GetEtaBinLimitCount(theCase=" << theCase << ")\n";
   switch(theCase) {
   case _ESFWC_1bin:
   case _ESFWC_2bin:
   case _ESFWC_3bin:
   case _ESFWC_4bin:
   case _ESFWC_5bin:
   case _ESFWC_6bin: 
     count=theCase; 
     //std::cout << "here!\n";
     break;
   case _ESFWC_3EB3EE: count=6; break;
   case _ESFWC_3EB3EEa: count=6; break;
   case _ESFWC_4EB3EE: count=7; break;
   case _ESFWC_2binNegs: count=4; break;
   case _ESFWC_3binNegs: count=6; break;
   case _ESFWC_4binNegs: count=8; break;
   case _ESFWC_5binNegs: count=10; break;
   case _ESFWC_6binNegs: count=12; break;
   case _ESFWC_3EB3EENegs: count=12; break;
   case _ESFWC_4EB3EENegs: count=14; break;
   case _ESFWC_5EB3EENegs: count=16; break;
   case _ESFWC_6EB3EENegs: count=18; break;
   case _ESFWC_7EB3EENegs: count=20; break;
   case _ESFWC_8EB3EENegs: count=22; break;
   case _ESFWC_5EBNegs: count=10; break;
   case _ESFWC_6EBNegs: count=12; break;
   case _ESFWC_4EENegs: count=8+1; break;
   default:
     std::cout << "EnergyScaleFactors::GetEtaBinLimitCount(" << theCase << "): code not ready\n";
   }
   //std::cout << "count=" << count << "\n";
   return count;
 }

 // ------------------------------------------------------------

 int ElectronEnergyScaleAdv_t::GetScaleFactorsCount(TEnergyScaleFactorsWorkCase_t wcase, TEnergyScaleFactorsFitModel_t fitModel) {
   int binCount=ElectronEnergyScaleAdv_t::GetEtaBinLimitCount(wcase);
   int n=0;
   switch(fitModel) {
   case _ESFModel_none: break;
   case _ESFModel_gauss: 
   case _ESFModel_gaussExp:
   case _ESFModel_doubleGauss:
   case _ESFModel_bifurGauss:
   case _ESFModel_cball:
   case _ESFModel_breitWigner:
   case _ESFModel_voigtian:
     n=binCount; break;
   case _ESFModel_gaussShift:
   case _ESFModel_cballShift:
     n=2*binCount; break;
   default:
     std::cout << "ESFAdv::GetScaleFactorsCount: unindentified model=" << fitModel << "\n";
   }
   return n;
 }

 // ------------------------------------------------------------

 int ElectronEnergyScaleAdv_t::GetSmearFactorsCount(TEnergyScaleFactorsWorkCase_t wcase, TEnergyScaleFactorsFitModel_t fitModel) {
   int binCount=ElectronEnergyScaleAdv_t::GetEtaBinLimitCount(wcase);
   int n=0;
   switch(fitModel) {
   case _ESFModel_none: break;
   case _ESFModel_gauss: n=binCount; break;
   case _ESFModel_gaussExp: n=2*binCount; break;
   case _ESFModel_doubleGauss: n=3*binCount; break;
   case _ESFModel_bifurGauss: n=2*binCount; break;
   case _ESFModel_cball: n=3*binCount; break;
   case _ESFModel_breitWigner: n=binCount; break;
   case _ESFModel_voigtian: n=2*binCount; break;
   case _ESFModel_gaussShift: n=binCount; break;
   case _ESFModel_cballShift: n=3*binCount; break;
   default:
     std::cout << "ESFAdv::GetSmearFactorsCount: unindentified model=" << fitModel << "\n";
   }
   return n;
 }

 // ------------------------------------------------------------

 int ElectronEnergyScaleAdv_t::FillEtaBinLimits(TEnergyScaleFactorsWorkCase_t theCase, double *b) {
   switch(theCase) {
   case _ESFWC_1bin: b[0]=0.0; b[1]=2.5; break;
   case _ESFWC_2bin: b[0]=0.0; b[1]=1.5; b[2]=2.5; break;
   case _ESFWC_3bin: b[0]=0.0; b[1]=0.8; b[2]=1.5; b[3]=2.5; break;
   case _ESFWC_4bin: b[0]=0.0; b[1]=1.0; b[2]=1.5; b[3]=2.; b[4]=2.5; break;
   case _ESFWC_5bin: b[0]=0.0; b[1]=0.5; b[2]=1.0; b[3]=1.5; b[4]=2.; b[5]=2.5; break;
   case _ESFWC_6bin:  b[0]=0.0; b[1]=0.4; b[2]=0.8; b[3]=1.2; b[4]=1.5; b[5]=2.0; b[6]=2.5; break;
   case _ESFWC_3EB3EE: b[0]=0.0; b[1]=0.6; b[2]=1.2; b[3]=1.5; b[4]=1.9; b[5]=2.2; b[6]=2.5; break;
   case _ESFWC_3EB3EEa: b[0]=0.0; b[1]=0.6; b[2]=1.2; b[3]=1.5; b[4]=2.1; b[5]=2.3; b[6]=2.5; break;
   case _ESFWC_4EB3EE: b[0]=0.0; b[1]=0.4; b[2]=0.8; b[3]=1.2; b[4]=1.5; b[5]=1.9; b[6]=2.2; b[7]=2.5; break;
    case _ESFWC_2binNegs: b[0]=-2.5; b[1]=-1.5; b[2]=0; b[3]=1.5; b[4]=2.5; break;
   case _ESFWC_3binNegs: b[0]=-2.5; b[1]=-1.5; b[2]=-0.8; b[3]=0.; b[4]=0.8; b[5]=1.5; b[6]=2.5; break;
   case _ESFWC_4binNegs: b[0]=-2.5; b[1]=-2.; b[2]=-1.5; b[3]=-1.0; b[4]=0.0; b[5]=1.0; b[6]=1.5; b[7]=2.; b[8]=2.5; break;
   case _ESFWC_5binNegs: b[0]=-2.5; b[1]=-2.; b[2]=-1.5; b[3]=-1.0; b[4]=-0.5; b[5]=0.0; b[6]=0.5; b[7]=1.0; b[8]=1.5; b[9]=2.; b[10]=2.5; break;
   case _ESFWC_6binNegs: b[0]=-2.5; b[1]=-2.; b[2]=-1.5; b[3]=-1.2; b[4]=-0.8; b[5]=-0.4; b[6]=0.0; b[7]=0.4; b[8]=0.8; b[9]=1.2; b[10]=1.5; b[11]=2.; b[12]=2.5; break;
   case _ESFWC_3EB3EENegs: b[0]=-2.5; b[1]=-2.2; b[2]=-1.9; b[3]=-1.5; b[4]=-1.2; b[5]=-0.6; b[6]=0.0; b[7]=0.6; b[8]=1.2; b[9]=1.5; b[10]=1.9; b[11]=2.2; b[12]=2.5; break;
   case _ESFWC_4EB3EENegs: b[0]=-2.5; b[1]=-2.2; b[2]=-1.9; b[3]=-1.5; b[4]=-1.2; b[5]=-0.8; b[6]=-0.4; b[7]=0.0; b[8]=0.4; b[9]=0.8; b[10]=1.2; b[11]=1.5; b[12]=1.9; b[13]=2.2; b[14]=2.5; break;
   case _ESFWC_5EB3EENegs: b[0]=-2.5; b[1]=-2.2; b[2]=-1.9; b[3]=-1.5; b[4]=-1.2; b[5]=-0.9; b[6]=-0.6; b[7]=-0.3; b[8]=0.0; b[9]=0.3; b[10]=0.6; b[11]=0.9; b[12]=1.2; b[13]=1.5; b[14]=1.9; b[15]=2.2; b[16]=2.5; break;
   case _ESFWC_6EB3EENegs: b[0]=-2.5; b[1]=-2.2; b[2]=-1.9; b[3]=-1.5; b[4]=-1.25; b[5]=-1.; b[6]=-0.75; b[7]=-0.5; b[8]=-0.25; b[9]=0.; b[10]=0.25; b[11]=0.5; b[12]=0.75; b[13]=1.0; b[14]=1.25; b[15]=1.5; b[16]=1.9; b[17]=2.2; b[18]=2.5; break;
   case _ESFWC_7EB3EENegs: b[0]=-2.5; b[1]=-2.2; b[2]=-1.9; b[3]=-1.5; b[4]=-1.2; b[5]=-1.; b[6]=-0.8; b[7]=-0.6; b[8]=-0.4; b[9]=-0.2; b[10]=0.; b[11]=0.2; b[12]=0.4; b[13]=0.6; b[14]=0.8; b[15]=1.0; b[16]=1.2; b[17]=1.5; b[18]=1.9; b[19]=2.2; b[20]=2.5; break;
   case _ESFWC_8EB3EENegs: b[0]=-2.5; b[1]=-2.2; b[2]=-1.9; b[3]=-1.5; b[4]=-1.33; b[5]=-1.14; b[6]=-0.95; b[7]=-0.76; b[8]=-0.57; b[9]=-0.38; b[10]=-0.19; b[11]=0.; b[12]=0.19; b[13]=0.38; b[14]=0.57; b[15]=0.76; b[16]=0.95; b[17]=1.14; b[18]=1.33; b[19]=1.5; b[20]=1.9; b[21]=2.2; b[22]=2.5; break;
   case _ESFWC_5EBNegs: b[0]=-1.5; b[1]=-1.2; b[2]=-0.9; b[3]=-0.6; b[4]=-0.3; b[5]=0.; b[6]=0.3; b[7]=0.6; b[8]=0.9; b[9]=1.2; b[10]=1.5; break;
   case _ESFWC_6EBNegs: b[0]=-1.5; b[1]=-1.25; b[2]=-1; b[3]=-0.75; b[4]=-0.5; b[5]=-0.25; b[6]=0.; b[7]=0.25; b[8]=0.5; b[9]=0.75; b[10]=1.; b[11]=1.25; b[12]=1.5; break;
   case _ESFWC_4EENegs: b[0]=-2.5; b[1]=-2.25; b[2]=-2.0; b[3]=-1.75; b[4]=-1.5; b[5]=1.5; b[6]=1.75; b[7]=2.0; b[8]=2.25; b[9]=2.5; break;
   default:
     std::cout << "EnergyScaleFactors::FillEtaBinLimits(" << theCase << "): code not ready /eta values/\n";
     return 0;
   }
   return 1;
 }

// ------------------------------------------------------------

TEnergyScaleFactorsWorkCase_t ElectronEnergyScaleAdv_t::GetNegsWorkCaseKind(TEnergyScaleFactorsWorkCase_t wcase1) {
  std::string name1=ElectronEnergyScaleAdv_t::WorkCaseShortName(wcase1);
  std::string name2=name1.substr(0,name1.size()-1)+std::string("Negs");
  TEnergyScaleFactorsWorkCase_t wcase2=IdentifyESFWorkCase(name2.c_str());
  std::cout << " GenNegsWorkCaseKind: transformed " << name1 << " to " << ElectronEnergyScaleAdv_t::WorkCaseShortName(wcase2) << "\n";
  return wcase2;
}

// ------------------------------------------------------------
#ifdef __hiLevelService__
int ElectronEnergyScaleAdv_t::GetSmearIdx_gamma(int &idx_min, int &idx_max) const { // Breit-Wigner, Voigtian

  idx_min=0*FEtaDivisionCount; idx_max=1*FEtaDivisionCount;
  return 1;
}
#endif

// ------------------------------------------------------------
#ifdef __hiLevelService__
int ElectronEnergyScaleAdv_t::GetSmearIdx_sigma(int &idx_min, int &idx_max) const { // Gaussian, GaussExp, Crystal ball
  idx_min=0*FEtaDivisionCount; idx_max=1*FEtaDivisionCount;
  return 1;
}
#endif

 // ------------------------------------------------------------
#ifdef __hiLevelService__
int ElectronEnergyScaleAdv_t::GetSmearIdx_smearSigma(int &idx_min, int &idx_max) const { // Voigtian
  idx_min=1*FEtaDivisionCount; idx_max=2*FEtaDivisionCount;
  return 1;
}
#endif

// ------------------------------------------------------------
#ifdef __hiLevelService__
int ElectronEnergyScaleAdv_t::GetSmearIdx_sigmaL(int &idx_min, int &idx_max) const  {// BifurGauss
  idx_min=0*FEtaDivisionCount; idx_max=1*FEtaDivisionCount; 
  return 1;
}
#endif

 // ------------------------------------------------------------
#ifdef __hiLevelService__
int ElectronEnergyScaleAdv_t::GetSmearIdx_sigmaR(int &idx_min, int &idx_max) const { // BifurGauss
  idx_min=1*FEtaDivisionCount; idx_max=2*FEtaDivisionCount;
  return 1;
}
#endif

 // ------------------------------------------------------------
#ifdef __hiLevelService__
int ElectronEnergyScaleAdv_t::GetSmearIdx_tau(int &idx_min, int &idx_max) const { // GaussExp
  idx_min=1*FEtaDivisionCount; idx_max=2*FEtaDivisionCount;
  return 1;
}
#endif

 // ------------------------------------------------------------
#ifdef __hiLevelService__
int ElectronEnergyScaleAdv_t::GetSmearIdx_parN(int &idx_min, int &idx_max) const { // Crystal ball
  idx_min=1*FEtaDivisionCount; idx_max=2*FEtaDivisionCount;
  return 1;
}
#endif

// ------------------------------------------------------------
#ifdef __hiLevelService__
int ElectronEnergyScaleAdv_t::GetSmearIdx_alpha(int &idx_min, int &idx_max) const { // Crystal ball
  idx_min=2*FEtaDivisionCount; idx_max=3*FEtaDivisionCount;
  return 1;
}
#endif

// ------------------------------------------------------------

 std::vector<int>* ElectronEnergyScaleAdv_t::BinIndices(TEnergyScaleFactorsWorkCase_t theCase, TEnergyScaleFactorsBECase_t eb) {
   std::vector<int>* idx=new std::vector<int>();
   int err=0;
   int division=0;
   int divisionEBneg=-1, divisionEEpos=-1;
   switch(theCase) {
   case _ESFWC_1bin: division=0; break;
   case _ESFWC_2bin: division=1; break;
   case _ESFWC_3bin: division=2; break;
   case _ESFWC_4bin: division=2; break;
   case _ESFWC_5bin: division=3; break;
   case _ESFWC_6bin: division=4; break;
   case _ESFWC_3EB3EE: division=3; break;
   case _ESFWC_3EB3EEa: division=3; break;
   case _ESFWC_4EB3EE: division=4; break;
   case _ESFWC_2binNegs: divisionEBneg=1; division=2; divisionEEpos=3; break;
   case _ESFWC_3binNegs: divisionEBneg=1; division=3; divisionEEpos=5; break;
   case _ESFWC_4binNegs: divisionEBneg=2; division=4; divisionEEpos=6; break;
   case _ESFWC_5binNegs: divisionEBneg=2; division=5; divisionEEpos=8; break;
   case _ESFWC_6binNegs: divisionEBneg=2; division=6; divisionEEpos=10; break;
   case _ESFWC_3EB3EENegs: divisionEBneg=3; division=6; divisionEEpos=9; break;
   case _ESFWC_4EB3EENegs: divisionEBneg=3; division=7; divisionEEpos=11; break;
   case _ESFWC_5EB3EENegs: divisionEBneg=3; division=8; divisionEEpos=13; break;
   case _ESFWC_6EB3EENegs: divisionEBneg=3; division=9; divisionEEpos=15; break;
   case _ESFWC_7EB3EENegs: divisionEBneg=3; division=10; divisionEEpos=17; break;
   case _ESFWC_8EB3EENegs: divisionEBneg=3; division=11; divisionEEpos=19; break;
   case _ESFWC_5EBNegs: divisionEBneg=0; division=5; divisionEEpos=10; break;
   case _ESFWC_6EBNegs: divisionEBneg=0; division=6; divisionEEpos=12; break;
   case _ESFWC_4EENegs: divisionEBneg=4; division=4; divisionEEpos=5; break;
  default:
    std::cout << "BinIndices: 'division' not ready for theCase=" << theCase << "\n";
    err=1;
  }
  const int imax=ElectronEnergyScaleAdv_t::GetEtaBinLimitCount(theCase);
  if (theCase<_ESFWC_BothSided) {
    switch(eb) {
    case _ESFWC_BEall: for (int i=0; i<imax; i++) idx->push_back(i); break;
    case _ESFWC_EB: for (int i=0; i<division; ++i) idx->push_back(i); break;
    case _ESFWC_EE: for (int i=division; i<imax; ++i) idx->push_back(i); break;
    default:
      std::cout << "BinIndices: filling is not ready for eb=" << eb << "\n";
      std::cout << "  note: work case is only absolute range\n";
      err=1;
    }
  }
  else {
    switch(eb) {
    case _ESFWC_BEall: for (int i=0; i<imax; i++) idx->push_back(i); break;
    case _ESFWC_EB: 
      if (theCase!=_ESFWC_4EENegs) {
	for (int i=division; i<divisionEEpos; ++i) idx->push_back(i); 
      }
      break;
    case _ESFWC_EE: 
      if ((theCase!=_ESFWC_5EBNegs) && (theCase!=_ESFWC_6EBNegs)) {
	for (int i=divisionEEpos; i<imax; ++i) idx->push_back(i); 
      }
      break;
    case _ESFWC_EBNeg: 
      if (theCase!=_ESFWC_4EENegs) {
	for (int i=divisionEBneg; i<division; ++i) idx->push_back(i); 
      }
      break;
    case _ESFWC_EENeg: 
      if ((theCase!=_ESFWC_5EBNegs) && (theCase!=_ESFWC_6EBNegs)) {
	for (int i=0; i<divisionEBneg; ++i) idx->push_back(i);
      }
      break;
    default:
      std::cout << "BinIndices: filling is not ready for eb=" << eb << " (both sides)\n";
      err=1;
    }
  }
  if (err) { delete idx; return NULL; }
  return idx;
}


// ------------------------------------------------------------

std::vector<int>* ElectronEnergyScaleAdv_t::DoubleBinIndices(TEnergyScaleFactorsWorkCase_t theCase, TEnergyScaleFactorsBEDoubleCase_t ebeb) {
  //std::cout << "ESFAdv::DoubleBinIndices for workCase=" << ElectronEnergyScaleAdv_t::WorkCaseName(theCase) << ", and ebeb=" << ElectronEnergyScaleAdv_t::BEDoubleCaseName(ebeb) << "\n";
  int imax=ElectronEnergyScaleAdv_t::GetEtaBinLimitCount(theCase);
  TEnergyScaleFactorsBECase_t eb1,eb2;
  switch(ebeb) {
  case _ESFWC_any: 
  case _ESFWC_all:
    eb1=_ESFWC_BEall; eb2=eb1; break;
  case _ESFWC_EBEB: eb1=_ESFWC_EB; eb2=eb1; break;
  case _ESFWC_EBEE: eb1=_ESFWC_EB; eb2=_ESFWC_EE; break;
  case _ESFWC_EEEE: eb1=_ESFWC_EE; eb2=eb1; break;
    //case _ESFWC_EBNegEBNeg: eb1=_ESFWC_EBNeg; eb2=eb1; break;
    //case _ESFWC_EBNegEENeg: eb1=_ESFWC_EENeg; eb2=_ESFWC_EBNeg; break;
    //case _ESFWC_EENegEENeg: eb1=_ESFWC_EENeg; eb2=eb1; break;
  default:
    std::cout << "DoubleBinIndices: not ready for ebeb=" << ebeb << "\n";
    return NULL;
  }

  int done=0;
  int iter=0;
  std::vector<int>* idx=new std::vector<int>();  
  do {
    //std::cout << "iter=" << iter << ", eb1=" << eb1 << ", eb2=" << eb2 << "\n";
    std::vector<int>* idx1=ElectronEnergyScaleAdv_t::BinIndices(theCase,eb1);
    std::vector<int>* idx2=ElectronEnergyScaleAdv_t::BinIndices(theCase,eb2);
    if (!idx1 || !idx2) {
      std::cout << "DoubleBinIndices: Failed to obtain bin indices for ebeb=" << ebeb << "\n";
      return NULL;
    }
    //int imax=ElectronEnergyScaleAdv_t::GetEtaBinLimitCount(theCase);
    int k=0;
    for (int i=0; i<imax; i++) {
      for (int j=i; j<imax; ++j, ++k) {
	if ((std::find(idx1->begin(),idx1->end(),i)!=idx1->end()) &&
	    (std::find(idx2->begin(),idx2->end(),j)!=idx2->end())) {
	  idx->push_back(k);
	}
      }
    }
    delete idx1;
    delete idx2;
    if (theCase>_ESFWC_BothSided) {
      switch(ebeb) {
      case _ESFWC_EBEB: 
	if (iter==0) { eb1=_ESFWC_EBNeg; eb2=_ESFWC_EB; }
	else if (iter==1) { eb1=_ESFWC_EBNeg; eb2=eb1; }
	else done=1;
	break;
      case _ESFWC_EBEE:
	if (iter==0) { eb1=_ESFWC_EBNeg; eb2=_ESFWC_EE; }
	else if (iter==1) { eb1=_ESFWC_EENeg; eb2=_ESFWC_EB; }
	else if (iter==2) { eb1=_ESFWC_EENeg; eb2=_ESFWC_EBNeg; }
	else done=1;
	break;
      case _ESFWC_EEEE:
	if (iter==0) { eb1=_ESFWC_EENeg; eb2=_ESFWC_EE; }
	else if (iter==1) { eb1=_ESFWC_EENeg; eb2=_ESFWC_EENeg; }
	else done=1;
	break;
      default:
	done=1;
      }
    }
    else done=1;
    iter++;
  } while(done!=1);
  std::sort(idx->begin(),idx->end());
  return idx;
}


// ------------------------------------------------------------

std::string ElectronEnergyScaleAdv_t::FitModelName(TEnergyScaleFactorsFitModel_t fm) {
  std::string name;
  switch(fm) {
  case _ESFModel_none: name="fit model none"; break;
  case _ESFModel_gauss: name="fit model Gauss"; break;
  case _ESFModel_gaussExp: name="fit model GaussExp"; break;
  case _ESFModel_doubleGauss: name="fit model doubleGauss"; break;
  case _ESFModel_bifurGauss: name="fit model bifurGauss"; break;
  case _ESFModel_cball: name="fit model CBall"; break;
  case _ESFModel_breitWigner: name="fit model breitWigner"; break;
  case _ESFModel_voigtian: name="fit model Voigtian"; break;
  case _ESFModel_gaussShift: name="fit model GaussShift (Gauss with a bias)"; break;
  case _ESFModel_cballShift: name="fit model CBallShift (CBall with a bias)"; break;
  default:
    char buf[30];
    sprintf(buf,"unknown fit model=%d",fm);
    name=buf;
  }
  return name;
}
// ------------------------------------------------------------

std::string ElectronEnergyScaleAdv_t::FitModelShortName(TEnergyScaleFactorsFitModel_t fm) {
  std::string name;
  switch(fm) {
  case _ESFModel_none: name="none"; break;
  case _ESFModel_gauss: name="Gauss"; break;
  case _ESFModel_gaussExp: name="GaussExp"; break;
  case _ESFModel_doubleGauss: name="doubleGauss"; break;
  case _ESFModel_bifurGauss: name="bifurGauss"; break;
  case _ESFModel_cball: name="CBall"; break;
  case _ESFModel_breitWigner: name="BreitWigner"; break;
  case _ESFModel_voigtian: name="Voigtian"; break;
  case _ESFModel_gaussShift: name="gaussShift"; break;
  case _ESFModel_cballShift: name="CBallShift"; break;
  default:
    char buf[30];
    sprintf(buf,"unknown fit model=%d",fm);
    name=buf;
  }
  return name;
}

// ------------------------------------------------------------

std::string ElectronEnergyScaleAdv_t::BEDoubleCaseName(TEnergyScaleFactorsBEDoubleCase_t beCase) {
  std::string s;
  switch(beCase) {
  case _ESFWC_any: s="any EB/EE"; break;
  case _ESFWC_EBEB: s="EB-EB"; break;
  case _ESFWC_EBEE: s="EB-EE"; break;
  case _ESFWC_EEEE: s="EE-EE"; break;
  default:
    std::cout << "ESFAdv::BEDoubleCaseName is not ready for beCase=" << beCase << "\n";
    s="unknown";
  }
  return s;
}

// ------------------------------------------------------------
// ------------------------------------------------------------

int ElectronEnergyScaleAdv_t::SaveASCIIFile(const char *filename) const {
  FILE *fout=fopen(filename,"w");
  fprintf(fout,"! g_esfWorkCase=%d (%s)\n",FWorkCase,this->WorkCaseName().c_str());
  fprintf(fout,"! g_esfWorkCaseShortName= %s\n",this->WorkCaseShortName().c_str());
  fprintf(fout,"! g_esfFitModel=%d (%s)\n",FFitModel,this->FitModelName().c_str());
  char buf[100];
  switch(FInvertedSF) {
  case 1: sprintf(buf,"scaling inverted\n"); break;
  case 2: sprintf(buf,"scaling sqrt\n"); break;
  case 3: sprintf(buf,"scaling invertedsqrt\n"); break;
  case 4: sprintf(buf,"scaling squared\n"); break;
  default:
    std::cout << "SaveASCIIFile: not ready for FInvertedSF=" << FInvertedSF << "\n";
    return 0;
  }
  fprintf(fout,buf);
  fprintf(fout,"! bins ");
  for (int i=0; i<=FEtaDivisionCount; ++i) { 
    fprintf(fout," %3.2lf",FEtaBinLimits[i]);
  }
  fprintf(fout,"\n");
  fprintf(fout,"MCOverData=%3.6lf\n", FMCtoData);
  fprintf(fout,"EtaDivisionCount=%d\n",FEtaDivisionCount);
  fprintf(fout,"ScalingFactorsCount=%d\n",FScFCount);
  fprintf(fout,"SmearingFactorsCount=%d\n",FSmFCount);
  for (int i=0; i<FScFCount; i++) {
    std::cout << "saving scF_" << i << std::endl;
    fprintf(fout,"scale_%d %12.6g %12.6g\n",i,FScaleFactors[i],FScaleFactorErrLo[i]); //,FScaleFactorErrHi[i]);
  }
  switch(FFitModel) {
  case _ESFModel_cball: {
    std::string name="smear";
    int idx0=0;
    for (int i=0; i<FSmFCount; ++i) {
      if (i==FEtaDivisionCount) { name="parN"; idx0=i; }
      else if (i==2*FEtaDivisionCount) { name="alpha"; idx0=i; }
      fprintf(fout,"%s_%d %12.6g %12.6g\n",name.c_str(),i-idx0,FSmearings[i],FSmearingErrLo[i]);
    }
  }
    break;
  default:
    for (int i=0; i<FSmFCount; ++i) {
      fprintf(fout,"smear_%d %12.6g %12.6g\n",i,FSmearings[i],FSmearingErrLo[i]); //,FSmearingErrHi[i]);
    }
  }
  fclose(fout);
  return 1;
}

// ------------------------------------------------------------

int ElectronEnergyScaleAdv_t::SaveASCIIFile(const char *mainfilename, const char *options) const {
  if (strstr(options,"<default>")==NULL) {
    std::cout << "ESFAdv::SaveASCIIFile(mainfilename): cannot understand option=" << options << "\n";
    return 0;
  }
  char buf[300];
  sprintf(buf,mainfilename);
  char *p=strrchr(buf,'/');
  if (!p) { sprintf(buf,"./"); p=buf+1; }
  sprintf(p+1,"testESF_%s_%s.inp",this->WorkCaseShortName().c_str(),this->FitModelShortName().c_str());
  int res=this->SaveASCIIFile(buf);
  if (!res) std::cout << " error in ESFAdv::SaveASCIIFile(mainfilename=" << mainfilename << ")\n";
  else {
    std::cout << "saved parameters to file <" << buf << ">\n";
  }
  return res;
}

// ------------------------------------------------------------

int ElectronEnergyScaleAdv_t::SaveRootFile(const char *filename) const {
  std::cout << "SaveRootFile: not ready\n";
  return 0;
  TFile outF(filename,"RECREATE");
  if (!outF.IsOpen()) {
    std::cout << "SaveRootFile: failed to create the parameters file <" << filename << ">\n";
    return 0;
  }
  //char buf[100];
  
  /*
  TString version="1.0";
  sprintf(buf,"energy scale factors for %s, %s",this->FitModelName().c_str(),this->WorkCaseName().c_str());
  TTree tree("energyScaleFactors",buf);
  TBranch *verBr=tree.Branch("fileVersion","TString",&version,32000,0);
  TVectorT<Int_t> info;
  info.ResizeTo(6);
  info[0]=FWorkCase; info[1]=FFitModel; info[2]=FScFCount; info[3]=FSmFCount;
  info[4]=FEtaDivisionCount; info[5]=FInvertedSF;
  TBranch *infoBr=tree.Branch("info","TVectorT<Int_t>",&info,32000,0);
  TBranch *mcToDataBr=tree.Branch("MCtoData","Double_t",&FMCtoData,32000,0);
  TBranch *etaBinLimitsBr=tree.Branch("etaBinLimits",TVectorT<Double_t>,
  outF.Close();
  return 1;
  */  
  //std::cout << "SaveRootFile is not ready\n";
  std::string version="1.0";
  outF.WriteObject(&version,"fileVersion");
  TString workCaseName=this->WorkCaseName();
  TString fitModelName=this->FitModelName();
  Int_t workCase=FWorkCase;
  outF.WriteObject(&workCase,"workCase");
  outF.WriteObject(&workCaseName,"workCaseName");
  outF.WriteObject(&FFitModel,"fitModel");
  outF.WriteObject(&fitModelName,"fitModelName");
  outF.WriteObject(&FScFCount,"scalingFactorCount");
  outF.WriteObject(&FSmFCount,"smearingFactorCount");
  outF.WriteObject(&FEtaDivisionCount,"etaDivisionCount");
  outF.WriteObject(&FInvertedSF,"invertedSFKind");
  outF.WriteObject(&FMCtoData,"mcToData");
  outF.WriteObject(FEtaBinLimits,"etaBinLimitsV");
  outF.WriteObject(FScaleFactors,"scaleFactorsV");
  outF.WriteObject(FSmearings,"smearingFactorsV");
  outF.WriteObject(FScaleFactorErrLo,"scaleFactorErrLoV");
  outF.WriteObject(FScaleFactorErrHi,"scaleFactorErrHiV");
  outF.WriteObject(FSmearingErrLo,"smearingFactorErrLoV");
  outF.WriteObject(FSmearingErrHi,"smearingFactorErrHiV");
  outF.Close();
  return 1;
}

// ------------------------------------------------------------

int ElectronEnergyScaleAdv_t::LoadRootFile(const char *filename) {
  std::cout << "LoadRootFile: not ready\n";
  return 0;
  TFile outF(filename);
  if (!outF.IsOpen()) {
    std::cout << "LoadRootFile: failed to create the parameters file <" << filename << ">\n";
    return 0;
  }
  //char buf[100];
  int err=0;
  std::string *pVersion=new std::string();
  std::string *pWorkCaseName=new std::string();
  std::string *pFitModelName=new std::string();
  int *pWorkCase=new int();
  int *pFitModel=new int();
  int *pScFCount=new int;
  int *pSmFCount=new int;
  int *pEtaDivisionCount=new int;
  int *pInvertedSF=new int;
  double *pMCtoData=new double;
  do {
    outF.GetObject("fileVersion",pVersion);
    if (pVersion && (*pVersion!="1.0")) {
      err=1; std::cout << "unknown file version (" << *pVersion << ")\n";
      break;
    }
    if (pVersion) std::cout << "the value obtained is <" << *pVersion << ">\n";
    outF.GetObject("workCase",pWorkCase);
    outF.GetObject("workCaseName",pWorkCaseName);
    outF.GetObject("fitModel",pFitModel);
    outF.GetObject("fitModelName",pFitModelName);
    if (pWorkCase && pFitModel) {
      this->WorkCase(*pWorkCase,*pFitModel);
      if (pWorkCaseName && (this->WorkCaseName()!=(*pWorkCaseName))) {
        std::cout << "version mismatch? workCaseName=" << this->WorkCaseName() << ", file gave " << (*pWorkCaseName) << "\n";
	err=1; break;
      }
    }
    else { 
      std::cout << "failed to obtain pWorkCase or pFitModel\n";
      if (pWorkCase) std::cout << " - pWorkCase=" << *pWorkCase << "\n";
      if (pFitModel) std::cout << " - pFitModel=" << *pFitModel << "\n";
      err=1; break;
    }
    outF.GetObject("scalingFactorCount",pScFCount);
    outF.GetObject("smearingFactorCount",pSmFCount);
    outF.GetObject("etaDivisionCount",pEtaDivisionCount);
    if (pScFCount && pSmFCount && pEtaDivisionCount) {
      if ((FScFCount!=*pScFCount) || (FSmFCount!=*pSmFCount) || (FEtaDivisionCount!=*pEtaDivisionCount)) {
	std::cout << " count value mismatch: " << FScFCount << " vs " << *pScFCount << ", " << FSmFCount << " vs " << *pSmFCount << ", " << FEtaDivisionCount << " vs " << *pEtaDivisionCount << "\n";
	err=1; break;
      }
    }
    else {
      std::cout << "failed to obtain pScFCount, pSmFCount or pEtaDivisionCount\n";
      err=1; break;
    }
    outF.GetObject("invertedSFKind",pInvertedSF);
    if (pInvertedSF) this->InvertedSF(*pInvertedSF); else { err=1; break; }
    outF.GetObject("mcToData",pMCtoData);
    if (pMCtoData) this->MCtoData(*pMCtoData); else { err=1; break; }
    outF.GetObject("etaBinLimitsV",this->FEtaBinLimits);
    outF.GetObject("scaleFactorsV",this->FScaleFactors);
    outF.GetObject("smearingFactorsV",this->FSmearings);
    outF.GetObject("scaleFactorErrLoV",this->FScaleFactorErrLo);
    outF.GetObject("scaleFactorErrHiV",this->FScaleFactorErrHi);
    outF.GetObject("smearingFactorErrLoV",this->FSmearingErrLo);
    outF.GetObject("smearingFactorErrHiV",this->FSmearingErrHi);
  } while (0);
  outF.Close();
  if (err) {
    std::cout << "LoadRootFile(" << filename << "): error was encountered\n";
    return 0;
  }
  return 1;
}

// ------------------------------------------------------------

int ElectronEnergyScaleAdv_t::LoadParameterSet(const char *eta_ranges_name, const char *fit_model_name, const char *file_name_base) {
  char buf[300];
  char temp[100];
  sprintf(temp,"fit model %s",fit_model_name);
  sprintf(buf,"%s_%s_%s.inp",file_name_base,this->WorkCaseShortName(IdentifyESFWorkCase(eta_ranges_name)).c_str(),this->FitModelShortName(IdentifyESFFitModel(temp,1)).c_str());
  int res=this->LoadASCIIFile(buf);
  if (!res) std::cout << " in method LoadParameterSet(" << eta_ranges_name << ", " << fit_model_name  << "," << file_name_base << "\n";
  return res;
}

// ------------------------------------------------------------

std::string ElectronEnergyScaleAdv_t::ProperScaleFactorValueWithErrForPrint(int i) const {
  double val=this->ProperScaleFactorValueForPrint(i);
  double valErr=this->ProperScaleFactorErrValueForPrint(i);
  char buf[20];
  if (valErr<9000) sprintf(buf,"%6.4lf +/- %4.1e",val,valErr);
  else sprintf(buf,"%6.4lf",val);
  return std::string(buf);
}

// ------------------------------------------------------------

std::string ElectronEnergyScaleAdv_t::ProperSmearFactorValueWithErrForPrint(int i) const {
  double val=this->SmearingFactor(i);
  double valErr=this->SmearingFactorErr(i);
  char buf[20];
  if (valErr<9000) sprintf(buf,"%6.4lf +/- %4.1e",val,valErr);
  else sprintf(buf,"%6.4lf",val);
  return std::string(buf);
}

// ------------------------------------------------------------

int ElectronEnergyScaleAdv_t::WorkCase(int wcase, int fitModel) {
  Clear();
  FWorkCase=wcase; FFitModel=fitModel;
  FMCtoData=1.0;
  TEnergyScaleFactorsWorkCase_t wcaseKind=TEnergyScaleFactorsWorkCase_t(wcase);
  FEtaDivisionCount=ElectronEnergyScaleAdv_t::GetEtaBinLimitCount(wcaseKind);
  FScFCount=ElectronEnergyScaleAdv_t::GetScaleFactorsCount(wcaseKind,TEnergyScaleFactorsFitModel_t(fitModel));
  FSmFCount=ElectronEnergyScaleAdv_t::GetSmearFactorsCount(wcaseKind,TEnergyScaleFactorsFitModel_t(fitModel));

  //std::cout << "got FEtaDivisionCount=" << FEtaDivisionCount << "\n";
  if ((FEtaDivisionCount<=0) || (FScFCount<=0) || (FSmFCount<=0)) {
    std::cout << "error in ElectronEnergyScaleAdv_t::WorkCase(work_case=" << wcase << ", fit_model=" << fitModel << ")\n";
    std::cout << "  derived EtaDivisionCount=" << FEtaDivisionCount << ", ScaleFactorCount=" << FScFCount << ", SmearFactorCount=" << FSmFCount << "\n";
    return 0;
  }

  FEtaBinLimits=new double[FEtaDivisionCount+1];
  assert(FEtaBinLimits);
  if (!ElectronEnergyScaleAdv_t::FillEtaBinLimits(wcaseKind,FEtaBinLimits)) {
    std::cout << "error in WorkCase=" << wcase << "\n";
    return 0;
  }
  
  FScaleFactors=new double [FScFCount];
  FSmearings=new double[FSmFCount];
  FScaleFactorErrLo=new double [FScFCount];
  FScaleFactorErrHi=new double [FScFCount];
  FSmearingErrLo=new double[FSmFCount];
  FSmearingErrHi=new double[FSmFCount];
  assert(FScaleFactors); assert(FSmearings);
  for (int i=0; i<FScFCount; i++) FScaleFactors[i]=9999.0;
  for (int i=0; i<FSmFCount; i++) FSmearings[i]=9999.0;
  for (int i=0; i<FScFCount; i++) FScaleFactorErrLo[i]=0.0;
  for (int i=0; i<FScFCount; i++) FScaleFactorErrHi[i]=0.0;
  for (int i=0; i<FSmFCount; i++) FSmearingErrLo[i]=0.0;
  for (int i=0; i<FSmFCount; i++) FSmearingErrHi[i]=0.0;
  return 1;
}

// ------------------------------------------------------------

int ElectronEnergyScaleAdv_t::Assign(const ElectronEnergyScaleAdv_t &sf) {
  this->Clear();
  if (!this->WorkCase(sf.WorkCase(),sf.FitModel())) { std::cout << "ElectronEnergyScaleAdv_t::Assign failed\n"; return 0; }
  this->FMCtoData=sf.MCtoData();
  this->FInvertedSF=sf.InvertedSF();
  for (int i=0; i<=this->BinCount(); ++i) {
    this->FEtaBinLimits[i]=sf.EtaBinLimit(i);
  }
  for (int i=0; i<this->ScFCount(); ++i) {
    this->FScaleFactors[i]=sf.ScaleFactor(i);
  }
  for (int i=0; i<this->SmFCount(); ++i) {
    this->FSmearings[i]=sf.Smearing(i);
  }
  for (int i=0; i<this->ScFCount(); ++i) {
    this->FScaleFactorErrLo[i]=sf.ScaleFactorErrLo(i);
  }
  for (int i=0; i<this->ScFCount(); ++i) {
    this->FScaleFactorErrHi[i]=sf.ScaleFactorErrHi(i);
  }
  for (int i=0; i<this->SmFCount(); ++i) {
    this->FSmearingErrLo[i]=sf.SmearingErrLo(i);
  }
  for (int i=0; i<this->SmFCount(); ++i) {
    this->FSmearingErrHi[i]=sf.SmearingErrHi(i);
  }
  return 1;
}

// ------------------------------------------------------------

int ElectronEnergyScaleAdv_t::Assign(const RooRealVar &mc_to_data, const std::vector<RooRealVar*> &scaling, int scaling_formula, const std::vector<RooRealVar*> &smearing) {
  if ((FScFCount!=int(scaling.size())) ||
      (FSmFCount!=int(smearing.size()))) {
    std::cout << "ElectronEnergyScaleAdv_t::Assign(RooRealVar): counts are inconsistent: FEtaDivisionCount=" << FEtaDivisionCount << ", FScFCount=" << FScFCount << ", FSmFCount=" << FSmFCount << ", scaling.size=" << scaling.size() << ", smearing.size=" << smearing.size() << "\n";
    return 0;
  }
  FInvertedSF=(scaling_formula)? 2:1;
  FMCtoData=mc_to_data.getVal();
  for (unsigned int i=0; i<scaling.size(); ++i) {
    FScaleFactors[i]=scaling[i]->getVal();
  }
  for (unsigned int i=0; i<smearing.size(); ++i) {
    FSmearings[i]=smearing[i]->getVal();
  }
  for (unsigned int i=0; i<scaling.size(); ++i) {
    FScaleFactorErrLo[i]=scaling[i]->getErrorLo();
  }
  for (unsigned int i=0; i<scaling.size(); ++i) {
    FScaleFactorErrHi[i]=scaling[i]->getErrorHi();
  }
  for (unsigned int i=0; i<smearing.size(); ++i) {
    FSmearingErrLo[i]=smearing[i]->getErrorLo();
  }
  for (unsigned int i=0; i<smearing.size(); ++i) {
    FSmearingErrHi[i]=smearing[i]->getErrorHi();
  }
  return 1;
}

// ------------------------------------------------------------

int ElectronEnergyScaleAdv_t::Assign(TEnergyScaleFactorsFitModel_t model, const RooRealVar &mc_to_data, const std::vector<RooRealVar*> &scaling, int scaling_formula, const std::vector<RooRealVar*> &smearing, const std::vector<RooRealVar*>* pars1, const std::vector<RooRealVar*>* pars2, const std::vector<RooRealVar*>* pars3) {
  const char *fncname="ESFAdv::Assign(model)";
  if (model!=this->FitModelKind()) {
    std::cout << "ESFAdv::Assign(model): models are different. Object fitModel=(" << ElectronEnergyScaleAdv_t::FitModelName(this->FitModelKind()) << "), suppliedModel=(" << ElectronEnergyScaleAdv_t::FitModelName(model) << "\n";
    return 0;
  }
  FInvertedSF=(scaling_formula)? 2:1;
  FMCtoData=mc_to_data.getVal();

  std::cout<< fncname << ". Counts scaling[" << scaling.size() << "], smearing[" << smearing.size() << "]";
  if (pars1) std::cout << ", pars1[" << pars1->size() << "]";
  if (pars2) std::cout << ", pars2[" << pars2->size() << "]";
  if (pars3) std::cout << ", pars3[" << pars3->size() << "]";
  std::cout << "\n";
  //std::cout << "DEBUG: this object: " << *this << std::endl;


  int res=1;
  int needsPars1=0;
  int needsPars2=0;
  int needsPars3=0;

  switch(model) {
  case _ESFModel_gauss: break;
  case _ESFModel_gaussExp: needsPars1=1; break;
  case _ESFModel_doubleGauss: needsPars1=1; needsPars2=1; break;
  case _ESFModel_bifurGauss: needsPars1=1; break;
  case _ESFModel_cball: needsPars1=1; needsPars2=1; break;
  case _ESFModel_breitWigner: break;
  case _ESFModel_voigtian: needsPars1=1; break;
  case _ESFModel_gaussShift: needsPars1=1; needsPars3=1; break;
  case _ESFModel_cballShift: needsPars1=1; needsPars2=1; needsPars3=1; break;
  default:
    std::cout << " ! unindentified model=" << model << "\n";
    res=0;
  }

  if (needsPars1 && !pars1) {
    std::cout << " ! the supplied pars1 is null\n";
    res=0;
  }
  if (needsPars2 && !pars2) {
    std::cout << " ! the supplied pars2 is null\n";
    res=0;
  }
  if (needsPars3 && !pars3) {
    std::cout << " ! the supplied pars3 is null\n";
    res=0;
  }
  std::cout << "here res=" << res << std::endl;
  //return 0;
  if (res) {
    if (needsPars1 && (smearing.size() != pars1->size())) {
      std::cout << " incorrect extra smearing parameters1 size (" << pars1->size() << "), expected=" << smearing.size() << "\n";
    }
    if (needsPars2 && (smearing.size() != pars2->size())) {
      std::cout << " incorrect extra smearing parameters2 size (" << pars2->size() << "), expected=" << smearing.size() << "\n";
    }
    if (needsPars3 && (scaling.size() != pars3->size())) {
      std::cout << " incorrect extra scaling parameters1 size (" << pars3->size() << "), expected=" << scaling.size() << "\n";
    }
  }

  if (!res) {
    std::cout << fncname << " model=" << model << "(" << ElectronEnergyScaleAdv_t::FitModelName(model) << ")\n";
    return res;
  }

  for (unsigned int i=0; i<scaling.size(); ++i) {
    FScaleFactors[i]=scaling[i]->getVal();
  }
  for (unsigned int i=0; i<scaling.size(); ++i) {
    FScaleFactorErrLo[i]=scaling[i]->getErrorLo();
  }
  for (unsigned int i=0; i<scaling.size(); ++i) {
    FScaleFactorErrHi[i]=scaling[i]->getErrorHi();
  }
  if (needsPars3) {
    for (unsigned int i=scaling.size(); i<scaling.size()+pars3->size(); ++i) {
      unsigned int ii=i-scaling.size();
      FScaleFactors[i]=(*pars3)[ii]->getVal();
    }
    for (unsigned int i=scaling.size(); i<scaling.size()+pars3->size(); ++i) {
      unsigned int ii=i-scaling.size();
      FScaleFactorErrLo[i]=(*pars3)[ii]->getErrorLo();
    }
    for (unsigned int i=scaling.size(); i<scaling.size()+pars3->size(); ++i) {
      unsigned int ii=i-scaling.size();
      FScaleFactorErrHi[i]=(*pars3)[ii]->getErrorHi();
    }
  }
  
  for (unsigned int i=0; i<smearing.size(); ++i) {
    FSmearings[i]=smearing[i]->getVal();
  }
  for (unsigned int i=0; i<smearing.size(); ++i) {
    FSmearingErrLo[i]=smearing[i]->getErrorLo();
  }
  for (unsigned int i=0; i<smearing.size(); ++i) {
    FSmearingErrHi[i]=smearing[i]->getErrorHi();
  }
  if (needsPars1) {
    for (unsigned int i=smearing.size(); i<smearing.size()+pars1->size(); ++i){
      unsigned int ii=i-smearing.size();
      FSmearings[i]=(*pars1)[ii]->getVal();
    }
    for (unsigned int i=smearing.size(); i<smearing.size()+pars1->size(); ++i){
      unsigned int ii=i-smearing.size();
      FSmearingErrLo[i]=(*pars1)[ii]->getErrorLo();
    }
    for (unsigned int i=smearing.size(); i<smearing.size()+pars1->size(); ++i){
      unsigned int ii=i-smearing.size();
      FSmearingErrHi[i]=(*pars1)[ii]->getErrorHi();
    }
  }
  if (needsPars2) {
    for (unsigned int i=2*smearing.size(); i<2*smearing.size()+pars2->size(); ++i){
      unsigned int ii=i-2*smearing.size();
      FSmearings[i]=(*pars2)[ii]->getVal();
    }
    for (unsigned int i=2*smearing.size(); i<2*smearing.size()+pars2->size(); ++i){
      unsigned int ii=i-2*smearing.size();
      FSmearingErrLo[i]=(*pars2)[ii]->getErrorLo();
    }
    for (unsigned int i=2*smearing.size(); i<2*smearing.size()+pars2->size(); ++i){
      unsigned int ii=i-2*smearing.size();
      FSmearingErrHi[i]=(*pars2)[ii]->getErrorHi();
    }
  }
  std::cout << "DEBUG: assigned obj: " << *this << std::endl;
  //std::cout << "DEBUG : leaving Assign" << std::endl;
  return res;
}

// ------------------------------------------------------------

int ElectronEnergyScaleAdv_t::Load(const char *fname) {
  std::cout << "Entered ElectronEnergyScaleAdv_t::Load" << std::endl;
  if (FEtaDivisionCount==0) {
    std::cout << "ElectronEnergyScaleAdv_t::Load the WorkCase had to be set\n";
    return 0;
  }
  std::string s;
  size_t p;
  unsigned int scaling_idx_max=0;
  int cycle=1;
  int all_parameters_found=0;
  while (cycle) {
    std::cout << "cycle=" << cycle << std::endl;
    int not_ready=0;
    std::vector<std::string> targets;
    targets.reserve(2);
    //targets.push_back("MCOverData");
    targets.push_back("scale_");

    std::ifstream fin;
    fin.open(fname);
    if (cycle==2) targets.push_back("smear_");
    else {
      if (this->FitModelKind()==_ESFModel_gauss) {
	targets.push_back("smear_");
      }
      else if (this->FitModelKind()==_ESFModel_gaussExp) {
	targets.push_back("expTau_");
      }
      else if (this->FitModelKind()==_ESFModel_doubleGauss) {
	targets.push_back("sigma1_");
	targets.push_back("sigma2_");
	targets.push_back("relfrac_");
      }
      else if (this->FitModelKind()==_ESFModel_bifurGauss) {
	targets.push_back("sigmaL_");
	targets.push_back("sigmaR_");
      }
      else if (this->FitModelKind()==_ESFModel_cball) {
	targets.push_back("smear_");
	targets.push_back("parN_");
	targets.push_back("alpha_");
      }
      else if (this->FitModelKind()==_ESFModel_breitWigner) {
	targets.push_back("gamma_");
      }
      else if (this->FitModelKind()==_ESFModel_voigtian) {
	targets.push_back("gamma_");
	targets.push_back("sigma_");
      }
      else { not_ready=1; }
    }
    if (this->FitModelKind()==_ESFModel_gaussShift) {
      targets.clear();
      targets.push_back("scale_");
      targets.push_back("bias_");
      scaling_idx_max=targets.size();
      targets.push_back("smear_");
      not_ready=0;
    }
    else if (this->FitModelKind()==_ESFModel_cballShift) {
      targets.clear();
      targets.push_back("scale_");
      targets.push_back("bias_");
      scaling_idx_max=targets.size();
      targets.push_back("smear_");
      if (cycle==1) {
	targets.push_back("parN_");
	targets.push_back("alpha_");
      }
      not_ready=0;
    }
    if (not_ready) {
      std::cout << "WARNING ESFAdv::Load is not adapted\n";
    }

    if (1) {
      std::cout << "search for targets[" << targets.size() << "]: ";
      for (unsigned int i=0; i<targets.size(); ++i) {
	std::cout << " {" << targets[i] << "}";
	if (i==scaling_idx_max) std::cout << " /end-of-scalings/ ";
      }
      std::cout << "\n";
    }

    std::vector<int> parsFound(targets.size());
    
    while (!fin.eof() && getline(fin,s)) {
      //std::cout << "studying <" << s << "\n";
      p=s.find("MCOverData");
      if (PosOk(p)) {
	//std::cout << "studying <" << s << "\n";
	FMCtoData=atof(s.c_str()+p+11);
	std::cout << "got MCtoData=" << FMCtoData << "\n";
      }
      else if (PosOk(s.find("scaling "))) {
	if (PosOk(s.find("scaling invertedsqrt"))) FInvertedSF=3;
	else if (PosOk(s.find("scaling inverted"))) FInvertedSF=1;
	else if (PosOk(s.find("scaling sqrt"))) { 
	  FInvertedSF=2;
	  //std::cout << "WARNING: forcing 'scaling invertedsqrt'\n"; FInvertedSF=3;
	  //std::cout << "WARNING: forcing 'scaling not inverted'\n"; FInvertedSF=0;
	}
	else if (PosOk(s.find("scaling squared"))) FInvertedSF=4;
	else FInvertedSF=0;
	std::cout << "got InvertedSF=" << FInvertedSF << "\n";
      }
      else {
	for (unsigned int ti=0; ti<targets.size(); ++ti) {
	  const std::string *t= & targets[ti];
	  const int tlen=(*t).size();
	  const int idxMax=(ti<=scaling_idx_max) ? FScFCount : FSmFCount;
	  double *dest=(ti<=scaling_idx_max) ? FScaleFactors : FSmearings;
	  double *destErr=(ti<=scaling_idx_max) ? FScaleFactorErrLo : FSmearingErrLo;
	  p=s.find(*t);
	  if (PosOk(p)) {
	    parsFound[ti]=1;
	    //std::cout << "studying <" << s << "\n";
	    int idx=0;
	    char ct=s[p+tlen];
	    if ((ct>='0') && (ct<='9')) { idx=atoi(s.c_str()+p+tlen); }
	    else {
	      idx=int(s[p+tlen]-'a');
	      if (idx<0) idx+=int('a'-'0')-1;
	    }
	    if (ti>1) idx += FEtaDivisionCount*(ti-1);
	    std::cout << (*t) << " idx=" << idx << "\n";
	    if (idx>idxMax) {
	      std::cout << "idx>idxMax: idx=" << idx << ", idxMax=" << idxMax << "\n";
	      return 0;
	    }
	    p=s.find_first_not_of(' ',p+tlen+2);
	    dest[idx]=atof(s.c_str()+p);
	    p=s.find_first_of(' ',p);
	    destErr[idx]=atof(s.c_str()+p);
	    std::cout << "got " << (*t) << " value=" << dest[idx] << ", err=" << destErr[idx] << "\n";
	  }
	}
      }
    }
    fin.close();

    all_parameters_found=1;
    for (unsigned int i=0; i<parsFound.size(); ++i) {
      if (parsFound[i]==0) {
	all_parameters_found=0;
	std::cout << "ElectronEnergyScaleAdv::Load: failed to locate all parameters\n";
	if (cycle==1) {
	  cycle=3;
	  std::cout << "... Will try to use another set of names\n";
	}
      }
    }
    cycle--;
    if (cycle==1) cycle=0;
  }

  std::cout << "Leaving ElectronEnergyScaleAdv_t::Load" << std::endl;
  return all_parameters_found;
}

// ------------------------------------------------------------

#ifdef UsePlotsAcc
int ElectronEnergyScaleAdv_t::ScaleValues(std::vector<mithep::AccEffData_t*> &data) const {
  int expect_count=FEtaDivisionCount*(FEtaDivisionCount+1)/2;
  if (data.size()!=expect_count) {
    std::cout << "ScaleValues: expected " << expect_count << " datasets, got=" << data.size() << "\n";
    return 0;
  }
  int k=0;
  for (unsigned int i=0; i<FEtaDivisionCount; ++i) {
    for (unsigned int j=i; j<FEtaDivisionCount; ++j, ++k) {
      double xfactor=ProperScaleFactorValue(i)*ProperScaleFactorValue(j);
      if (k<=20) std::cout << "scaling by " << ProperScaleFactorValue(i) << " * " << ProperScaleFactorValue(j) << "\n";
      if (xfactor>9000) {
	std::cout << " ScaleValues: unprepared  FInvertedSF=" << FInvertedSF << " case\n";
	return 0;
      }
      mithep::AccEffData_t* d=data[k];
      d->ScaleEnergy(xfactor); 
      if (this->FitModelKind()>_ESFModel_has_shift) {
	double shift=GetScaleFactorValue(int(i),int(j),1);
	d->ShiftEnergy(shift);
      }
    }
  }
  return 1;
}
#endif

// ------------------------------------------------------------

#ifdef UsePlotsAcc
int ElectronEnergyScaleAdv_t::ScaleValues(std::vector<RooDataSet*> &data, RooRealVar &x, const char *variable_name) const {
  int expect_count=FEtaDivisionCount*(FEtaDivisionCount+1)/2;
  if (data.size()!=expect_count) {
    std::cout << "ScaleValues(RooDataSet): expected " << expect_count << " datasets, got=" << data.size() << "\n";
    return 0;
  }

  std::vector<RooDataSet*> vec=data;
  data.clear();
  int k=0;
  for (unsigned int i=0; i<FEtaDivisionCount; ++i) {
    for (unsigned int j=i; j<FEtaDivisionCount; ++j, ++k) {
      double xfactor=ProperScaleFactorValue(i)*ProperScaleFactorValue(j);
      if (k<20) std::cout << "scaling by " << ProperScaleFactorValue(i) << " * " << ProperScaleFactorValue(j) << "\n";
      if (xfactor>9000) {
	std::cout << " ScaleValues: unprepared  FInvertedSF=" << FInvertedSF << " case\n";
	return 0;
      }
      
      RooArgSet value(x);
      RooDataSet *dt=new RooDataSet("dt","dt",RooArgSet(x));
      data.push_back(dt);
      RooDataSet *v=vec[k];
      for (unsigned int ii=0; ii<v->numEntries(); ++ii) {
	const RooArgSet *arg=v->get(ii);
	if (this->FitModelKind()>_ESFModel_has_shift) {
	  double shift=GetScaleFactorValue(int(i),int(j),1);
	  x.setVal(xfactor * arg->getRealValue(variable_name,0,kTRUE) + shift);
	}
	else {
	  x.setVal(xfactor * arg->getRealValue(variable_name,0,kTRUE));
	}
	dt->add(value);
      }
    }
  }
  return 1;
}
#endif

// ------------------------------------------------------------

int ElectronEnergyScaleAdv_t::ScaleValues(std::vector<std::vector<double>*> &data) const {
  unsigned int expect_count=FEtaDivisionCount*(FEtaDivisionCount+1)/2;
  if (data.size()!=expect_count) {
    std::cout << "ScaleValues(doubleVV): expected " << expect_count << " datasets, got=" << data.size() << "\n";
    return 0;
  }
  const unsigned int etaDivCount=(unsigned int)(FEtaDivisionCount);
  int k=0;
  for (unsigned int i=0; i<etaDivCount; ++i) {
    for (unsigned int j=i; j<etaDivCount; ++j, ++k) {
      double xfactor=ProperScaleFactorValue(i)*ProperScaleFactorValue(j);
      if (k<20) std::cout << "scaling by " << ProperScaleFactorValue(i) << " * " << ProperScaleFactorValue(j) << "\n";
      if (xfactor>9000) {
	std::cout << " ScaleValues: unprepared  FInvertedSF=" << FInvertedSF << " case\n";
	return 0;
      }
      std::vector<double>* d=data[k];
      double shift=0.;
      if (this->FitModelKind()>_ESFModel_has_shift) {
	shift=GetScaleFactorValue(int(i),int(j),1);
      }
      for (unsigned int ii=0; ii<d->size(); ++ii) {
	(*d)[ii]=((*d)[ii]*xfactor + shift);
      }
    }
  }
  return 1;
}

// ------------------------------------------------------------

#ifdef UsePlotsAcc
int ElectronEnergyScaleAdv_t::InvScaleValues(std::vector<mithep::AccEffData_t*> &data) const {
  int expect_count=FEtaDivisionCount*(FEtaDivisionCount+1)/2;
  if (data.size()!=expect_count) {
    std::cout << "ScaleValues: expected " << expect_count << " datasets, got=" << data.size() << "\n";
    return 0;
  }
  int k=0;
  for (unsigned int i=0; i<FEtaDivisionCount; ++i) {
    for (unsigned int j=i; j<FEtaDivisionCount; ++j, ++k) {
      double xfactor=1/(ProperScaleFactorValue(i)*ProperScaleFactorValue(j));
      if (k<20) std::cout << "scaling by 1/(" << ProperScaleFactorValue(i) << " * " << ProperScaleFactorValue(j) << ")\n";
      if (xfactor>9000) {
	std::cout << " ScaleValues: unprepared  FInvertedSF=" << FInvertedSF << " case\n";
	return 0;
      }
      mithep::AccEffData_t* d=data[k];
      d->ScaleEnergy(xfactor); 
      if (this->FitModelKind()>_ESFModel_has_shift) {
	double shift=GetScaleFactorValue(int(i),int(j),1);
	d->ShiftEnergy(-shift);
      }
    }
  }
  return 1;
}
#endif


// ------------------------------------------------------------

int ElectronEnergyScaleAdv_t::InvScaleValues(std::vector<std::vector<double>*> &data) const {
  unsigned int expect_count=FEtaDivisionCount*(FEtaDivisionCount+1)/2;
  if (data.size()!=expect_count) {
    std::cout << "InvScaleValues(doubleVV): expected " << expect_count << " datasets, got=" << data.size() << "\n";
    return 0;
  }
  const unsigned int etaDivCount=(unsigned int)(FEtaDivisionCount);
  int k=0;
  for (unsigned int i=0; i<etaDivCount; ++i) {
    for (unsigned int j=i; j<etaDivCount; ++j, ++k) {
      double xfactor=1/(ProperScaleFactorValue(i)*ProperScaleFactorValue(j));
      if (k<20) std::cout << "scaling by 1/(" << ProperScaleFactorValue(i) << " * " << ProperScaleFactorValue(j) << ")\n";
      if (xfactor>9000) {
	std::cout << " ScaleValues: unprepared  FInvertedSF=" << FInvertedSF << " case\n";
	return 0;
      }
      std::vector<double>* d=data[k];
      double shift=0.;
      if (this->FitModelKind()>_ESFModel_has_shift) {
	shift=GetScaleFactorValue(int(i),int(j),1);
      }
      for (unsigned int ii=0; ii<d->size(); ++ii) {
	(*d)[ii]=((*d)[ii]*xfactor + shift);
      }
    }
  }
  return 1;
}

// ------------------------------------------------------------

int ElectronEnergyScaleAdv_t::Verify_EnoughSmearingVars(const char *testee, int needs_positions, const std::vector<RooRealVar*> &vars, const std::vector<RooRealVar*> *smPars2, const std::vector<RooRealVar*> *smPars3) const {
  const char *fncname="Verify_EnoughSmearingVars";
  std::string msgStart=fncname; msgStart+=" for "; msgStart+=testee;
  int needsSMParCount=ElectronEnergyScaleAdv_t::GetSmearFactorsCount(this->WorkCaseKind(),this->FitModelKind());
  int needsSM2=(needsSMParCount>needs_positions) ? 1:0;
  int needsSM3=(needsSMParCount>2*needs_positions) ? 1:0;
  int err=0;
  if (needsSM2 && !smPars2) { err=1; }
  if (needsSM3 && !smPars3) { err=2; }
  if (err) {
    std::cout << msgStart << ": needsSm2=" << needsSM2 << " smPars2=NULL? " << ((smPars2==NULL) ? "yes":"no") << "\n";
    std::cout << msgStart << ": needsSm3=" << needsSM3 << " smPars3=NULL? " << ((smPars3==NULL) ? "yes":"no") << "\n";
    return 0;
  }
  unsigned int needs_positions_ui=(unsigned int)(needs_positions);

  if (needs_positions_ui!=vars.size()) {
    std::cout << msgStart << ": expected " << needs_positions << " main variables, got " << vars.size() << "\n";
    return 0;
  }
  if (needsSM2 && (needs_positions_ui!=smPars2->size())) {
    std::cout << msgStart << ": expected " << needs_positions << " smPars2 variables, got " << smPars2->size() << "\n";
    return 0;
  }
  if (needsSM3 && (needs_positions_ui!=smPars3->size())) {
    std::cout << msgStart << ": expected " << needs_positions << " smPars3 variables, got " << smPars3->size() << "\n";
    return 0;
  }
  return 1;
}

// ------------------------------------------------------------

int ElectronEnergyScaleAdv_t::Verify_EnoughScalingVars(const char *testee, int needs_positions, const std::vector<RooRealVar*> &vars, const std::vector<RooRealVar*> *scPars2) const {
  const char *fncname="Verify_EnoughScalingVars";
  std::string msgStart=fncname; msgStart+=" for "; msgStart+=testee;
  int needsSCParCount=ElectronEnergyScaleAdv_t::GetScaleFactorsCount(this->WorkCaseKind(),this->FitModelKind());
  int needsSC2=(needsSCParCount>needs_positions) ? 1:0;
  int err=0;
  if (needsSC2 && !scPars2) { err=1; }
  if (err) {
    std::cout << msgStart << ": needsSc2=" << needsSC2 << " scPars2=NULL? " << ((scPars2==NULL) ? "yes":"no") << "\n";
    std::cout << "needs_positions=" << needs_positions << ", needsSCParCount=" << needsSCParCount << "\n";
    return 0;
  }
  unsigned int needs_positions_ui=(unsigned int)(needs_positions);

  if (needs_positions_ui!=vars.size()) {
    std::cout << msgStart << ": expected " << needs_positions << " main variables, got " << vars.size() << "\n";
    return 0;
  }
  if (needsSC2 && (needs_positions_ui!=scPars2->size())) {
    std::cout << msgStart << ": expected " << needs_positions << " scPars2 variables, got " << scPars2->size() << "\n";
    return 0;
  }
  return 1;
}

// ------------------------------------------------------------

int ElectronEnergyScaleAdv_t::CopyScalingFactors(std::vector<RooRealVar*> &vars, std::vector<RooRealVar*> *scPars2) const {
  if (!Verify_EnoughScalingVars("CopyScalingFactors",FEtaDivisionCount,vars,scPars2)) return 0;

  for (unsigned int i=0; i<vars.size(); ++i) vars[i]->setVal(FScaleFactors[i]);

  if (scPars2) { 
    for (unsigned int i=vars.size(); i<2*vars.size(); ++i) {
      (*scPars2)[i-vars.size()]->setVal(FScaleFactors[i]);
    }
  }
  return 1;
}

// ------------------------------------------------------------

int ElectronEnergyScaleAdv_t::TransferScalingFactors(std::vector<RooRealVar*> &vars, std::vector<RooRealVar*> *scPars2) const {
  unsigned int expected_count=FEtaDivisionCount*(FEtaDivisionCount+1)/2;
  if (!Verify_EnoughScalingVars("TransferScalingFactors",expected_count,vars,scPars2)) return 0;
  std::cout << "TransferScalingFactors is not handling the sqrt/invsqrt nature of the factors\n"; return 0;
  unsigned int k=0;
  unsigned int etaDivisionCount=(unsigned int)(FEtaDivisionCount);
  for (unsigned int i=0; i<etaDivisionCount; ++i) {
    for (unsigned int j=i; j<etaDivisionCount; ++j, ++k) {
      vars[k]->setVal(sqrt(FScaleFactors[i]*FScaleFactors[i]));
    }
  }
  if (scPars2) {
    k=0;
    for (unsigned int ii=0; ii<etaDivisionCount; ++ii) {
      unsigned int i=ii+etaDivisionCount;
      for (unsigned int jj=ii; jj<etaDivisionCount; ++jj, ++k) {
	unsigned int j=jj+etaDivisionCount;
	(*scPars2)[k]->setVal(sqrt(FScaleFactors[i]*FScaleFactors[i]+FScaleFactors[j]*FScaleFactors[j]));
      }
    }
  }
  return 1;
}

// ------------------------------------------------------------

int ElectronEnergyScaleAdv_t::CopySmearingFactors(std::vector<RooRealVar*> &vars, std::vector<RooRealVar*> *smPars2, std::vector<RooRealVar*> *smPars3) const {
  if (!Verify_EnoughSmearingVars("CopySmearingFactors",FEtaDivisionCount,vars,smPars2,smPars3)) return 0;

  for (unsigned int i=0; i<vars.size(); ++i) vars[i]->setVal(FSmearings[i]);

  if (smPars2 && (smPars2->size()==vars.size())) { 
    for (unsigned int i=vars.size(); i<2*vars.size(); ++i) {
      (*smPars2)[i-vars.size()]->setVal(FSmearings[i]);
    }
  }
  if (smPars3 && (smPars3->size()==vars.size())) { 
    for (unsigned int i=2*vars.size(); i<3*vars.size(); ++i) {
      (*smPars3)[i-2*vars.size()]->setVal(FSmearings[i]);
    }
  }
  return 1;
}

// ------------------------------------------------------------

int ElectronEnergyScaleAdv_t::TransferSmearingFactors(std::vector<RooRealVar*> &vars, std::vector<RooRealVar*> *smPars2, std::vector<RooRealVar*> *smPars3) const {
  unsigned int expected_count=FEtaDivisionCount*(FEtaDivisionCount+1)/2;
  if (!Verify_EnoughSmearingVars("TransferSmearingFactors",expected_count,vars,smPars2,smPars3)) return 0;

  unsigned int k=0;
  unsigned int etaDivisionCount=(unsigned int)(FEtaDivisionCount);
  for (unsigned int i=0; i<etaDivisionCount; ++i) {
    for (unsigned int j=i; j<etaDivisionCount; ++j, ++k) {
      vars[k]->setVal(sqrt(FSmearings[i]*FSmearings[i]+FSmearings[j]*FSmearings[j]));
    }
  }
  if (smPars2) {
    k=0;
    for (unsigned int ii=0; ii<etaDivisionCount; ++ii) {
      unsigned int i=ii+etaDivisionCount;
      for (unsigned int jj=ii; jj<etaDivisionCount; ++jj, ++k) {
	unsigned int j=jj+etaDivisionCount;
	(*smPars2)[k]->setVal(sqrt(FSmearings[i]*FSmearings[i]+FSmearings[j]*FSmearings[j]));
      }
    }
  }
  if (smPars3) {
    k=0;
    for (unsigned int ii=0; ii<etaDivisionCount; ++ii) {
      unsigned int i=ii+2*etaDivisionCount;
      for (unsigned int jj=ii; jj<etaDivisionCount; ++jj, ++k) {
	unsigned int j=jj+2*etaDivisionCount;
	(*smPars3)[k]->setVal(sqrt(FSmearings[i]*FSmearings[i]+FSmearings[j]*FSmearings[j]));
      }
    }
  }
  return 1;
}

// ------------------------------------------------------------

//int ElectronEnergyScaleAdv_t::CopyValues(std::vector<RooRealVar*> &scale, std::vector<RooRealVar*> &smear) const {
//  if ((scale.size()!=FScCount) || (smear.size()!=FSmCount)) {
//  }
//  return 1;
//}

// ------------------------------------------------------------

int ElectronEnergyScaleAdv_t::PrintValuesHTMLTable(const TString *foutname, const std::vector<const ElectronEnergyScaleAdv_t*> *others, const std::vector<std::string> *names) const {
  int res=1;
  if (foutname) {
    std::ofstream htmlfile;
    htmlfile.open(foutname->Data());
    res=this->PrintValuesHTMLTable(htmlfile,0,others,names);
    htmlfile.close();
  }
  else {
    res=this->PrintValuesHTMLTable(std::cout,1,others,names);
  }
  if (!res) std::cout << "error in PrintValuesHTMLTable(foutname)\n";
  return res;
}

// ------------------------------------------------------------

int ElectronEnergyScaleAdv_t::PrintValuesHTMLTable(std::ostream& out, int append, const std::vector<const ElectronEnergyScaleAdv_t*> *others, const std::vector<std::string> *names) const {
  const char *fncname="PrintValuesHTMLTable";
  //const char *tdwidth15="<td width=\"25%\">";
  const char *tdwidth25="<td width=\"25%\">";
  int ok=1;
  if (others) {
    for (unsigned int i=0; ok && (i<others->size()); ++i) {
      if ((*others)[i]->WorkCase()!=FWorkCase) {
	std::cout << "WorkCase=" << FWorkCase << ", other[" << i << "]->WorkCase=" << (*others)[i]->WorkCase() << "\n";
	ok=0;
      }
    }
    if (!ok) {
      std::cout << fncname << ": the supplied work cases are different\n";
      return 0;
    }
  }
  if (others && names) {
    if (others->size()+1 != names->size()) {
      std::cout << " other->size=" << others->size() << ", names->size=" << names->size() << " (expected 1 extra name)\n";
      return 0;
    }
  }

  unsigned int N=(unsigned int)(FEtaDivisionCount);

  // header
  if (!append) {
    out << "<!DOCTYPE html" << endl;
    out << "    PUBLIC \"-//W3C//DTD HTML 3.2//EN\">" << endl;
    out << "<html>" << endl;
    out << "<head><title>Mass scale factors</title></head>" << endl;
    out << "<body bgcolor=\"EEEEEE\">" << endl;
  }
  else out << "<br>" << endl;

  out << "<h3 style=\"text-align:left; color:DD6600;\">Mass scale factors</h3>" << endl;

  out << "<table border=\"0\" cellspacing=\"2\" width=\"100%\">" << endl;  
  out << "<tr>" << endl;
  out << "<td width=\"90%\">" << "FitModel: <strong>" << ElectronEnergyScaleAdv_t::FitModelName(this->FitModelKind()) << "</strong>; etaRangeSet: " << ElectronEnergyScaleAdv_t::WorkCaseName(this->WorkCaseKind()) << ", there are " << FEtaDivisionCount << " eta ranges</td>" << endl;
  out << "</tr>" << endl;
  // (
  out << "<tr><td width=\"90%\"> Eta division points: "; 
  for (unsigned int i=0; i<=N; ++i) {
    //out << " " << (i+1) << ") " << FEtaBinLimits[i] << ";"; 
    out << "  " << FEtaBinLimits[i]; 
  }
  out << "</td></tr>" << endl;
  out << "</table>" << endl;
  // values
  out << "<table border=\"0\" cellspacing=\"2\" width=\"100%\">" << endl;  
  out << "<tr>" << endl;
  out << tdwidth25 << "Parameter</td>";
  out << tdwidth25 << "Eta range</td>";
  if (names) {
    for (unsigned int i=0; i<names->size(); ++i) {
      out << tdwidth25 << (*names)[i] << "</td>";
    }
    out << endl;
  }
  else {
    out << tdwidth25 << "values</td>" << endl;
  }
  out << "</tr>" << endl;
  out << "<tr>" << endl << tdwidth25 << "MCtoData</td>" << tdwidth25 << "</td>" << tdwidth25 << FMCtoData << "</td>"; 
  if (others) {
    for (unsigned int i=0; i<others->size(); ++i) {
      out << tdwidth25 << (*others)[i]->MCtoData() << "</td>";
    }
  }
  out << "</tr>" << endl;
  for (unsigned int pi=0; pi<(unsigned int)(FScFCount); ++pi) {
    unsigned int ibin=(pi%N);
    std::string scaleParName=(pi<N) ? "scale " : "scalePar2 ";
    out << "<tr>" << endl;
    out << tdwidth25 << scaleParName << (ibin+1) << "</td>" << tdwidth25;
    out << " (" << FEtaBinLimits[ibin] << ", " << FEtaBinLimits[ibin+1] << ')';
    out << "</td>";
    out << tdwidth25 << ProperScaleFactorValueWithErrForPrint(pi) << "</td>"; 
    if (others) {
      for (unsigned int i=0; i<others->size(); ++i) {
	out << tdwidth25 << (*others)[i]->ProperScaleFactorValueWithErrForPrint(pi) << "</td>";
      }
    }
  }
  out << "</tr>" << endl;
  for (unsigned int pi=0; pi<(unsigned int)(this->FSmFCount); ++pi) {
    unsigned int ibin=(pi%N);
    std::string smearParName="smear(sigma) ";
    switch (this->FitModelKind()) {
    case _ESFModel_gauss: break;
    case _ESFModel_gaussExp: 
      if (pi>=N) smearParName="smear(tau) "; break;
    case _ESFModel_doubleGauss:
      if (pi>=2*N) smearParName="smear(rel.fract.) ";
      else if (pi>=N) smearParName="smear(sigma2) ";
      else smearParName="smear(sigma) ";
      break;
    case _ESFModel_bifurGauss:
      if (pi>=N) smearParName="smear(sigmaR) ";
      else smearParName="smear(sigmaL) ";
      break;
    case _ESFModel_cball: 
    case _ESFModel_cballShift: 
      if (pi>=2*N) smearParName="smear(alpha) ";
      else if (pi>=N) smearParName="smear(parN) ";
      break;
    case _ESFModel_breitWigner:
      smearParName="smear(gamma) ";
      break;
    case _ESFModel_voigtian:
      if (pi>=N) smearParName="smear(sigma) ";
      else smearParName="smear(gamma) ";
      break;
    default:
      std::cout << "WARNING: unknown smear parameter name for this model\n";
    }
    out << "<tr>" << endl << tdwidth25 << smearParName << (ibin+1) << "</td>" << tdwidth25;
    out << " (" << FEtaBinLimits[ibin] << ',' << FEtaBinLimits[ibin+1] << ')';    
    out << "</td>";
    out << tdwidth25 << ProperSmearFactorValueWithErrForPrint(pi) << "</td>"; 
    if (others) {
      for (unsigned int i=0; i<others->size(); ++i) {
	out << tdwidth25 << (*others)[i]->ProperSmearFactorValueWithErrForPrint(pi) << "</td>";
      }
    }
    out << "</tr>" << endl;
  }
  out << "</table>" << endl;
  if (!append) {
    out << "</body>" << endl;
    out << "</html>" << endl;
  }
  else out << "<br>" << endl;
  return 1;
}

// ------------------------------------------------------------
// ------------------------------------------------------------

void PrepareCombinationInfos(TEnergyScaleFactorsWorkCase_t workCase, TEnergyScaleFactorsBEDoubleCase_t theCase, std::vector<std::vector<int>*> &combine, std::vector<TString> &combination_names) {
  std::vector<int> *indices=NULL;
  if ((theCase==_ESFWC_all) || (theCase==_ESFWC_any)) {
    indices=ElectronEnergyScaleAdv_t::DoubleBinIndices(workCase,_ESFWC_any);
    if (indices && (indices->size()>0)) {
      combine.push_back(indices);
      combination_names.push_back("all EB-EE");
    }
  }
  if ((theCase==_ESFWC_all) || (theCase==_ESFWC_EBEB)) {
    indices=ElectronEnergyScaleAdv_t::DoubleBinIndices(workCase,_ESFWC_EBEB);
    if (indices && (indices->size()>0)) {
      combine.push_back(indices);
      combination_names.push_back("EB-EB");
    }
  }
  if ((theCase==_ESFWC_all) || (theCase==_ESFWC_EBEE)) {
    indices=ElectronEnergyScaleAdv_t::DoubleBinIndices(workCase,_ESFWC_EBEE);
    if (indices && (indices->size()>0)) {
      combine.push_back(indices);
      combination_names.push_back("EB-EE");
    }
  }
  if ((theCase==_ESFWC_all) || (theCase==_ESFWC_EEEE)) {
    indices=ElectronEnergyScaleAdv_t::DoubleBinIndices(workCase,_ESFWC_EEEE);
    if (indices && (indices->size()>0)) {
      combine.push_back(indices);
      combination_names.push_back("EE-EE");
    }
  }
  return;
}

// ------------------------------------------------------------

int ElectronEnergyScaleAdv_t::PrintValuesHTMLTable(std::vector<std::string> &html_lines, const std::vector<const ElectronEnergyScaleAdv_t*> *others, const std::vector<std::string> *names) const {
  const char* destFile="dummy.html";
  std::ofstream fout;
  fout.open(destFile);
  this->PrintValuesHTMLTable(fout,1,others,names);
  fout.close();
  std::ifstream fin;
  fin.open(destFile);
  std::string s;
  while (!fin.eof() && getline(fin,s)) {
    html_lines.push_back(s);
  }
  fin.close();
  return (html_lines.size())?1:0;
}

// ------------------------------------------------------------

void ElectronEnergyScaleAdv_t::PrepareHRNameEndings(std::vector<TString> &hrEnds) const {
  hrEnds.clear();
  //std::cout << "entered PrepareHRNameEndings" << std::endl;
  int n=FEtaDivisionCount;
  int k=0;
  hrEnds.reserve(n*(n+1)/2);
  TString name;
  char buf[100];  
  for (int i=0; i<n; ++i) {
    double ib0=FEtaBinLimits[i];
    double ib1=FEtaBinLimits[i+1];
    for (int j=i; j<n; ++j, ++k) {
      double jb0=FEtaBinLimits[j];
      double jb1=FEtaBinLimits[j+1];
      sprintf(buf,"%d,%d (%3.1lf..%3.1lf) x (%3.1lf..%3.1lf)",i+1,j+1,ib0,ib1,jb0,jb1);
      //std::cout << "prepared buf=" << buf << std::endl;
      name=buf;
      //std::cout << "assigned" << std::endl;
      hrEnds.push_back(name);
      //std::cout << "inserted" << std::endl;
    }
  }
  //std::cout << "leaving ESF_t::PrepareHRNameEndings" << std::endl;
  return;
}

// ------------------------------------------------------------

mithep::myHistoClass_t* ElectronEnergyScaleAdv_t::CreateHisto(TString name, double xmin, double xmax, int plot_scale, int use_negative_eta_values) const {
  if (plot_scale) name.Append(" scale"); else name.Append(" smear");
  int binCount=0;
  int includeGap=1;
  const double gapHi=1.566;
  const double gapLo=1.4442;
  if ((xmin<0) && (FEtaBinLimits[0]<0)) { binCount=FEtaDivisionCount+3+2*includeGap; }
  else if ((xmin<0) && (FEtaBinLimits[0]==0.)) { binCount=2*FEtaDivisionCount+3+2*includeGap; }
  else if ((xmin==0) && (FEtaBinLimits[0]<0)) { binCount=FEtaDivisionCount/2+2+includeGap; }
  else if ((xmin==0) && (FEtaBinLimits[0]==0.)) { binCount=FEtaDivisionCount+2+includeGap; }
  double bins[binCount];
  double vals[binCount];
  double valErrs[binCount];
  for (int i=0; i<binCount; ++i) vals[i]=0;
  for (int i=0; i<binCount; ++i) valErrs[i]=0;
  std::cout << "xmin=" << xmin << ", xmax=" << xmax << ", FEtaBinLimits[0]=" << FEtaBinLimits[0] << "\n";
  if ((xmin<0) && (FEtaBinLimits[0]<0)) {
    int i=0;
    bins[i]=xmin; i++;
    vals[0]=0;
    valErrs[0]=9999;
    int shift=0;
    for (int j=0; j<FEtaDivisionCount+1; ++j) {
      int idx=i+shift+j;
      bins[idx]=FEtaBinLimits[j];
      vals[idx]=(plot_scale) ? FScaleFactors[j] : FSmearings[j];
      valErrs[idx]=(plot_scale) ? ScaleFactorErr(j) : SmearFactorErr(j);
      if (includeGap) {
	if (FEtaBinLimits[j]==-1.5) {
	  bins[idx]=-gapHi;
	  bins[idx+1]=-gapLo;
	  vals[idx+1]=vals[idx]; valErrs[idx+1]=valErrs[idx]; 
	  vals[idx]=0; valErrs[idx]=0;
	  shift++;
	}
	else if (FEtaBinLimits[j]==1.5) {
	  bins[idx]=gapLo;
	  bins[idx+1]=gapHi;
	  vals[idx+1]=vals[idx]; valErrs[idx+1]=valErrs[idx];
	  vals[idx]=0; valErrs[idx]=0;
	  shift++;
	}
      }
    }
    vals[i+shift+FEtaDivisionCount]=0;
    valErrs[i+shift+FEtaDivisionCount]=0;
    bins[i+shift+FEtaDivisionCount+1]=xmax;
  }
  else if ((xmin<0) && (FEtaBinLimits[0]==0.)) {
    int i=0;
    int shift=0;
    bins[i]=xmin; i++;
    for (int j=0; j<FEtaDivisionCount+1; ++j) {
      int ridx=FEtaDivisionCount-j;
      int idx=i+shift+j;
      bins[idx]=-FEtaBinLimits[ridx];
      vals[idx-1]=(plot_scale) ? FScaleFactors[ridx] : FSmearings[ridx];
      valErrs[idx-1]=(plot_scale) ? ScaleFactorErr(ridx) : SmearFactorErr(ridx);
      if (includeGap && (bins[idx]==-1.5)) {
	bins[idx]=-gapHi; bins[idx+1]=-gapLo;
	vals[idx]=0; //vals[idx-1]=0;
	valErrs[idx]=0; //valErrs[idx-1]=0;
	shift++; 
      }
    }
    i+=(FEtaDivisionCount+1-1); // exclude double 0
    vals[i+shift]=(plot_scale) ? FScaleFactors[0] : FSmearings[0];
    valErrs[i+shift]=(plot_scale) ? ScaleFactorErr(0) : SmearFactorErr(0);
    for (int j=0; j<FEtaDivisionCount+1; ++j) {
      int idx=i+j+shift;
      bins[idx]=FEtaBinLimits[j];
      vals[idx]=(plot_scale) ? FScaleFactors[j] : FSmearings[j];
      valErrs[idx]=(plot_scale) ? ScaleFactorErr(j) : SmearFactorErr(j);
      if (1 && includeGap && (bins[idx]==1.5)) {
	bins[idx]=gapLo; bins[idx+1]=gapHi;
	vals[idx+1]= vals[idx]; vals[idx]=0;
	valErrs[idx+1]=valErrs[idx]; valErrs[idx]=0;
	shift++;
      }
    }
    i+=(FEtaDivisionCount+1);
    bins[i+shift]=xmax;
  }
  else if ((xmin==0) && (FEtaBinLimits[0]<0)) {
    int i=0, shift=0;
    if (!use_negative_eta_values) {
      for (int j=0; j<FEtaDivisionCount+1; ++j) {
	if (FEtaBinLimits[j]>=0) { 
	  bins[i]=FEtaBinLimits[j]; 
	  vals[i]=(plot_scale) ? FScaleFactors[j] : FSmearings[j];
	  valErrs[i]=(plot_scale) ? ScaleFactorErr(j) : SmearFactorErr(j);
	  i++;
	  if (includeGap && (bins[i-1]==1.5)) {
	    bins[i-1]=gapLo;
	    bins[i]=gapHi;
	    vals[i]=vals[i-1]; vals[i-1]=0;
	    valErrs[i]=valErrs[i-1]; valErrs[i-1]=0;
	    i++;
	  }
	}
      }
      bins[i]=xmax;
    }
    else {
      shift=includeGap;
      bins[binCount-1]=xmax;
      for (int j=0; j<FEtaDivisionCount+1; ++j) {
	int idx=FEtaDivisionCount/2-j-1;
	if (FEtaBinLimits[j]<0) { 
	  std::cout << "j=" << j << ", idx=" << idx << "\n";
	  bins[idx+shift+1]=-FEtaBinLimits[j]; 
	  vals[idx+shift]=(plot_scale) ? FScaleFactors[j] : FSmearings[j];
	  valErrs[idx+shift]=(plot_scale) ? ScaleFactorErr(j) : SmearFactorErr(j);
	  if (includeGap && (bins[idx+shift+1]==1.5)) {
	    bins[idx+shift+1]=gapHi;
	    bins[idx+shift]=gapLo;
	    vals[idx+shift-1]=vals[idx+shift]; vals[idx+shift]=0.;
	    valErrs[idx+shift-1]=valErrs[idx+shift]; valErrs[idx+shift]=0.;
	    shift--;
	  }
	}       
      }
      //bins[i]=xmax;
    }
  }
  else if ((xmin==0) && (FEtaBinLimits[0]==0)) {
    int shift=0;
    for (int i=0; i<FEtaDivisionCount+1; ++i) {
      bins[i+shift]=FEtaBinLimits[i];
      vals[i+shift]=(plot_scale) ? FScaleFactors[i] : FSmearings[i];
      valErrs[i+shift]=(plot_scale) ? ScaleFactorErr(i) : SmearFactorErr(i);
      if (includeGap && (bins[i+shift]==1.5)) {
	bins[i+shift]=gapLo;
	bins[i+shift+1]=gapHi;
	vals[i+shift+1]=vals[i+shift]; vals[i+shift]=0;
	valErrs[i+shift+1]=valErrs[i+shift]; valErrs[i+shift]=0;
	shift++;
      }
    }
    bins[FEtaDivisionCount+shift+1]=xmax;
  }
  std::cout << "binCount=" << binCount << "\n";
  std::cout << "bins: "; for (int i=0; i<binCount; ++i) std::cout << " " << bins[i]; std::cout << "\n";
  std::cout << "vals: "; for (int i=0; i<binCount-1; ++i) std::cout << " " << vals[i]; std::cout << "\n";
  std::cout << "valErrs: "; for (int i=0; i<binCount-1; ++i) std::cout << " " << valErrs[i]; std::cout << "\n";
  mithep::myHistoClass_t *h=new mithep::myHistoClass_t(name,name,binCount-1,bins);
  h->SetStats(0);
  int idx=(xmin<0) ? 1:0;
  for (int i=0; i<binCount-2; ++i) {
    h->SetBinContent(i+1+idx,vals[i+idx]);
    h->SetBinError(i+1+idx,valErrs[i+idx]);
  }
  return h;
}

// ------------------------------------------------------------

int ElectronEnergyScaleAdv_t::ChangeToNegs(const ElectronEnergyScaleAdv_t &sf) {
  const char *fncname="ESFAdv::ChangeToNegs: ";
  this->Clear();
  if (sf.WorkCaseKind()>_ESFWC_BothSided) {
    std::cout << fncname << "warning the provided object is already 'negs' (" << sf.WorkCaseShortName() << ")\n";
    int res=this->Assign(sf);
    if (!res) std::cout << " in method " << fncname << "\n";
    return res;
  }
  TEnergyScaleFactorsWorkCase_t wcase=ElectronEnergyScaleAdv_t::GetNegsWorkCaseKind(sf.WorkCaseKind());
  if (wcase!=_ESFWC_none) {
    if (!this->WorkCase(wcase,sf.FitModel())) {
      std::cout << "Error in " << fncname << "\n";
      return 0;
    }
    FMCtoData=0.5*sf.FMCtoData;
    FInvertedSF=sf.FInvertedSF;
    double *source=NULL, *dest=NULL;
    for (int loop=0; loop<6; ++loop) {
      int count=(loop<3) ? sf.FScFCount : sf.FSmFCount;
      switch(loop) {
      case 0: source=sf.FScaleFactors; dest=this->FScaleFactors; break;
      case 1: source=sf.FScaleFactorErrLo; dest=this->FScaleFactorErrLo; break;
      case 2: source=sf.FScaleFactorErrHi; dest=this->FScaleFactorErrHi; break;
      case 3: source=sf.FSmearings; dest=this->FSmearings; break;
      case 4: source=sf.FSmearingErrLo; dest=this->FSmearingErrLo; break;
      case 5: source=sf.FSmearingErrHi; dest=this->FSmearingErrHi; break;
      default:
	std::cout << fncname << "not ready for loop=" << loop << "\n";
	return 0;
      }
      for (int i=0; i<count; ++i) {
	int idx=count-1-i;
	dest[i]=source[idx];
      }
      for (int i=0; i<count; ++i) {
	int idx=count+i;
	dest[idx]=source[i];
      }
    }
  }
  else {
    std::cout << fncname << "failed to convert the work case\n";
    return 0;
  }
  return 1;
}

// ------------------------------------------------------------
// ------------------------------------------------------------

int DetermineESFWorkCase(const char *filename, ElectronEnergyScaleAdv_t &esf) {
  std::cout << "entered DetermineESFWorkCase" << std::endl;
  std::ifstream fin;
  fin.open(filename);
  if (!fin) {
    std::cout << "DetermineESFWorkCase: failed to open a file <" << filename << ">\n";
    return 0;
  }
  std::string s;
  size_t pos;
  int res=0; //int(_ESFWC_6bin);
  TEnergyScaleFactorsFitModel_t fitModel=_ESFModel_none;
  while (!fin.eof() && getline(fin,s)) {
    pos=s.find("g_esfWorkCaseShortName=");
    if (PosOk(pos)) {
      TEnergyScaleFactorsWorkCase_t wcase=IdentifyESFWorkCase(s.c_str());
      res=int(wcase);
      
      std::cout << "got g_esfWorkCase=" << res << " (" << ElectronEnergyScaleAdv_t::WorkCaseName(wcase) << ")\n";
      std::cout << std::endl;
    }
    else if (PosOk(s.find("fit model"))) {
      //std::cout << "trying to identify fit model from " << s << std::endl;
      TEnergyScaleFactorsFitModel_t fitModelLoc=_ESFModel_none;      
      fitModelLoc=IdentifyESFFitModel(s.c_str(),0);
      if (fitModelLoc!=_ESFModel_none) fitModel=fitModelLoc;
    }
  }
  fin.close();
  std::cout << "determined fitModel=" << fitModel << "\n";
  if (res && (fitModel==_ESFModel_none)) fitModel=_ESFModel_gauss;
  if (res && (fitModel!=_ESFModel_none)) {
    std::cout << "Loading the file" << std::endl;
    esf.WorkCase(res,int(fitModel));
    res=esf.Load(filename);
  }
  return res;
}

// ------------------------------------------------------------

TEnergyScaleFactorsWorkCase_t DetermineESFWorkCase(const TApplication &theApp) {
  TEnergyScaleFactorsWorkCase_t key=_ESFWC_none;
  for (int i=1; i<theApp.Argc(); ++i) {
    if (strstr(theApp.Argv(i),"esf") &&
	!strstr(theApp.Argv(i),"esfModel")) {
      key=IdentifyESFWorkCase(theApp.Argv(i));
      if (key!=_ESFWC_none) { 
	int bin_count=ElectronEnergyScaleAdv_t::GetEtaBinLimitCount(key);
	double *bins=new double[bin_count+1];
	ElectronEnergyScaleAdv_t::FillEtaBinLimits(key,bins);
	std::cout << " ** esfWorkCase identified " << ElectronEnergyScaleAdv_t::WorkCaseName(key) << ". Bins[" << bin_count << "]=";
	for (int ib=0; ib<bin_count+1; ib++) std::cout << " " << bins[ib];
	std::cout << "\n";
	delete bins;
      }
      else {
	std::cout << " *!! failed to identify esfWorkCase from <" << theApp.Argv(i) << ">\n";
      }
    }
  }
  return key;
}

// ------------------------------------------------------------

TEnergyScaleFactorsFitModel_t DetermineESFFitModel(const TApplication &theApp) {
  TEnergyScaleFactorsFitModel_t key=_ESFModel_none;
  for (int i=1; i<theApp.Argc(); ++i) {
    if (strstr(theApp.Argv(i),"esfModel")) {
      std::string s="fit model ";
      s+=theApp.Argv(i);
      key=IdentifyESFFitModel(s.c_str(),0);
      if (key!=_ESFModel_none) { 
	std::cout << " ** esfFitModel identified " << ElectronEnergyScaleAdv_t::FitModelName(key) << "\n";
      }
      else {
	std::cout << " *!! failed to identify esfFitModel from <" << theApp.Argv(i) << ">\n";
      }
    }
  }
  return key;
}

// ------------------------------------------------------------
// ------------------------------------------------------------

#ifdef EtaEtaMass_H
int ElectronEnergyScaleAdv_t::ProcessEEMFile(const char *mc_file_name, const char *data_file_name, std::vector<std::vector<double>*> &mcData, std::vector<std::vector<double>*> &expData) {
  EtaEtaMassData_t *eem = new EtaEtaMassData_t();
  int res=1;
  const int etaEtaCount=this->EtaDivisionCount()*(this->EtaDivisionCount()+1)/2;
  mcData.clear(); mcData.reserve(etaEtaCount+1);
  expData.clear(); expData.reserve(etaEtaCount+1);
  // Read each file twice. First time we will determine the number of masses
  // in each (eta1,eta2) bin and allocate memory. The second time we will
  // store the mass value
  const int optimize=1; // optimize memory
  for (int loop=0; res && (loop<2+2*optimize); ++loop) {
    //std::cout << "loop=" << loop << std::endl;
    std::vector<int> counts(etaEtaCount);
    const char *fname=(loop%2==0) ? mc_file_name : data_file_name;
    if (!fname) continue;
    std::vector<std::vector<double>*> *data= (loop%2==0) ? &mcData : &expData;
    TFile *fin = new TFile(fname);
    assert(fin);
    TTree *tree = (TTree*)fin->Get("Data"); assert(tree);
    tree->SetBranchAddress("Data",&eem);
    TBranch *branch = tree->GetBranch("Data");
    UInt_t entryMax = tree->GetEntries();
    for (UInt_t ientry=0; ientry < entryMax; ++ientry) {
      branch->GetEntry(ientry);
      int idx=this->PrepareEtaEtaBinIdx(eem->eta1(),eem->eta2());
      //std::cout << loop << " got " << (*eem) << ", idx=" << idx << "\n";
      if ((idx>=0) && (idx < etaEtaCount)) {
	if (optimize && (loop/2==0)) {
	  counts[idx]++;
	}
	else {
	  std::vector<double>* store=(*data)[idx];
	  //if (store->capacity()<store->size()+5) 
	  //std::cout << "adding to idx=" << idx << ", count/capacity=" << store->size() << "/" << store->capacity() << std::endl;
	  store->push_back(eem->mass());
	}
      }
    }
    //delete branch;
    //delete tree;
    delete fin;
    if (optimize && (loop/2==0)) {
      // allocate memory
      unsigned int k=0;
      for (int i=0; i<this->EtaDivisionCount(); ++i) {
	for (int j=i; j<this->EtaDivisionCount(); ++j, ++k) {
	  std::vector<double>* tmp=new std::vector<double>();
	  assert(tmp);
	  //std::cout << "counts[k]=" << counts[k] << "\n";
	  tmp->reserve(counts[k]+1);
	  data->push_back(tmp);
	}
      }
    }
  }
  return 1;
}
#endif

// ------------------------------------------------------------

#ifdef EtaEtaMass_H
int ElectronEnergyScaleAdv_t::ProcessEEMFileApproximateMCWeight(const char *mc_file_name, const char *data_file_name, std::vector<std::vector<double>*> &mcData, std::vector<std::vector<double>*> &expData, double unitWeight) {
  std::cout << "entered ProcessEEMFileApproximateMCWeight" << std::endl;
  EtaEtaMassData_t *eem = new EtaEtaMassData_t();
  int res=1;
  int debug=0;
  TH1F *hraw=(debug==0) ? 0 : new TH1F("hraw","hraw",60,60.,120.);
  TH1F *h1=(debug==0) ? 0 : new TH1F("h1","h1",60,60.,120.);
  TH1F *h2=(debug==0) ? 0 : new TH1F("h2","h2",60,60.,120.);
  if (h1) { h1->SetLineColor(kGreen); h1->SetMarkerColor(kGreen); }
  if (h2) { h2->SetLineColor(kRed);  h2->SetMarkerColor(kRed) ; }
  if (debug) { assert(hraw); assert(h1); assert(h2); }

  const int etaEtaCount=this->EtaDivisionCount()*(this->EtaDivisionCount()+1)/2;
  mcData.clear(); mcData.reserve(etaEtaCount+1);
  expData.clear(); expData.reserve(etaEtaCount+1);
  std::vector<std::vector<double>*> remMCMasses,remMCWeight;
  remMCMasses.reserve(etaEtaCount+1);
  remMCWeight.reserve(etaEtaCount+1);
  double wsum=0., wcount=0.;

  // Read each file twice. First time we will determine the number of masses
  // in each (eta1,eta2) bin and allocate memory. The second time we will
  // store the mass value
  const int optimize=1; // optimize memory
  for (int loop=0; res && (loop<2+2*optimize); ++loop) {
    //std::cout << "loop=" << loop << std::endl;
    std::vector<double> counts(etaEtaCount);
    const char *fname=(loop%2==0) ? mc_file_name : data_file_name;
    std::vector<std::vector<double>*> *data= (loop%2==0) ? &mcData : &expData;
    TFile *fin = new TFile(fname);
    assert(fin);
    TTree *tree = (TTree*)fin->Get("Data"); assert(tree);
    tree->SetBranchAddress("Data",&eem);
    TBranch *branch = tree->GetBranch("Data"); assert(branch);
    UInt_t entryMax = tree->GetEntries();
    for (UInt_t ientry=0; ientry < entryMax; ++ientry) {
      branch->GetEntry(ientry);
      if (fabs(eem->mass()-90)>35) continue;
      int idx=this->PrepareEtaEtaBinIdx(eem->eta1(),eem->eta2());
      //std::cout << loop << " got " << (*eem) << ", idx=" << idx << "\n";
      if ((idx>=0) && (idx < etaEtaCount)) {
	if (optimize && (loop/2==0)) {
	  counts[idx]+=eem->weight();
	  if ((loop==0) && (fabs(eem->mass()-90)<=30)) {
	    wsum+=eem->weight(); wcount+=1.;
	    if (fabs(unitWeight)>1e-3) { wcount=1; wsum=unitWeight; }
	  }
	}
	else {
	  std::vector<double>* store=(*data)[idx];
	  //if (store->capacity()<store->size()+5) 
	  //std::cout << "adding to idx=" << idx << ", count/capacity=" << store->size() << "/" << store->capacity() << std::endl;
	  if (loop%2==1) { 
	    // data - no change
	    store->push_back(eem->mass());
	  }
	  else {
	    double w=eem->weight()*(wcount/wsum);
	    if (debug && (ientry<100) && (eem->mass()>=60) && (eem->mass()<=120)) std::cout << "mass=" << eem->mass() << ", w_raw=" << eem->weight() << ", w=" << w << "\n";
	    if (hraw) hraw->Fill(eem->mass(),w);
	    while (w>1.0) { 
	      store->push_back(eem->mass()); 
	      w-=1.0; 
	      if (h1) h1->Fill(eem->mass(),1.);
	      if (h2) h2->Fill(eem->mass(),1.);
	    }
	    if (fabs(w)>1e-6) {
	      if (!remMCMasses[idx] || !remMCWeight[idx]) { std::cout << "Null ptrs" << std::endl; assert(0); }
	      remMCMasses[idx]->push_back(eem->mass());
	      remMCWeight[idx]->push_back(w);
	    }
	  }
	}
      }
    }
   
    //delete branch;
    //delete tree;
    delete fin;
    if (optimize && (loop/2==0)) {
      // allocate memory
      unsigned int k=0;
      for (int i=0; i<this->EtaDivisionCount(); ++i) {
	for (int j=i; j<this->EtaDivisionCount(); ++j, ++k) {
	  std::vector<double>* tmp=new std::vector<double>();
	  assert(tmp);
	  //std::cout << "counts[k]=" << counts[k] << "\n";
	  tmp->reserve(int(counts[k])+1);
	  data->push_back(tmp);
	}
      }
      if (loop==0) {
	k=0;
	for (int i=0; i<this->EtaDivisionCount(); ++i) {
	  for (int j=i; j<this->EtaDivisionCount(); ++j, ++k) {
	    std::vector<double>* tmp=new std::vector<double>();	assert(tmp);
	    tmp->reserve(int(counts[k])+1);
	    remMCMasses.push_back(tmp);
	    tmp=new std::vector<double>();  	assert(tmp);
	    tmp->reserve(int(counts[k])+1);
	    remMCWeight.push_back(tmp);
	  }
	}
      }
    }
  }

  // Add up the remaining MC weights
  int idx=0;
  std::cout << "adding up the remaining MC weights" << std::endl;
  for (int i=0; i<this->EtaDivisionCount(); ++i) {
    for (int j=i; j<this->EtaDivisionCount(); ++j, ++idx) {
      if (remMCWeight[idx]->size()==0) continue;
      std::vector<double>* remMasses=remMCMasses[idx];
      std::vector<double>* remWeight=remMCWeight[idx];
      int do_sort=1;
      if (do_sort && !QuickSort(*remMasses,*remWeight)) {
	std::cout << "QuickSort failed\n";
	return 0;
      }
      int cx=1;
      //int ccc=100;
      while (remMasses->size()>2) {
	unsigned int ci=0, cj=cx;
	double m1=(*remMasses)[ci], m2=(*remMasses)[cj];
	for (unsigned int ii=cj+1; ii<remMasses->size(); ++ii) {
	  std::cout << "ci=" << ci << ", m1=" << m1 << "; cj=" << cj << ", m2=" << m2 << "; ii=" << ii << "; mass=" << (*remMasses)[ii] << "\n";
	  if (fabs(m1-m2)>fabs(m1-(*remMasses)[ii])) { cj=ii; m2=(*remMasses)[ii]; }
	  if (fabs(m1-m2)<0.2) { std::cout << "break\n"; break; }
	  if (do_sort && (fabs(m1-(*remMasses)[ii])>0.5)) { std::cout << "sorted. break\n"; break; }
	}
	if (fabs(m1-m2)>0.2) {
	  if ((*remWeight)[ci]>0.6) {
	    mcData[idx]->push_back(m1);
	    if (h2) h2->Fill(m1,1.);
	  }
	  cx=1;
	  //std::cout << "erasing ci=" << ci << ", " << (*remMasses)[ci] << "\n";
	  remMasses->erase(remMasses->begin()+ci);
	  remWeight->erase(remWeight->begin()+ci);
	  std::cout << "modified ci=" << ci << ", " << (*remMasses)[ci] << "\n";
	  //return 0;
	}
	else {
	  double w=(*remWeight)[ci]+(*remWeight)[cj];
	  if (w>=1.) {
	    mcData[idx]->push_back(m1);
	    w-=1.;
	    if (h2) h2->Fill(m1,1.);
	  }
	  (*remWeight)[ci]=w;
	  cx=cj;
	  std::cout << "branch mi-m2<=0.2 : m1=" << m1 << ", m2=" << m2 << "\n";
	  remMasses->erase(remMasses->begin()+cj);
	  remWeight->erase(remWeight->begin()+cj);
	  //ccc--; if (ccc<=0) return 0;
	}
      }
    }
  }
  if (debug) {
    TCanvas *ctest = new TCanvas("ctest","ctest",600,600);
    if (hraw) hraw->Draw();
    if (h1) h1->Draw("LP same");
    if (h2) h2->Draw("LP same");
    ctest->Update();
    std::cout << "enter a symbol...";
    char x; cin >> x;
  }
  return 1;
}
#endif


// ------------------------------------------------------------
// ------------------------------------------------------------
