#include <RooGaussian.h>
#include <RooGaussModel.h>
#include <RooGExpModel.h>
#include <RooNumConvPdf.h>
//#include "RooNumIntConfig.h" // adjust numerical integration
#include "fitMassScaleFactorsE.hh"
//#include "ElectronEnergyScaleAdv.hh"

//#ifndef __myLib__
#include "HelpingTools.hh"
//#endif

using RooFit::LineColor;
using RooFit::LineStyle;
using RooFit::LineWidth;

TEnergyScaleFactorsWorkCase_t DYTools::g_esfWorkCase=_ESFWC_1bin;
//TEnergyScaleFactorsFitModel_t DYTools::g_esfFitModel=_ESFModel_doubleGauss;
TEnergyScaleFactorsFitModel_t DYTools::g_esfFitModel=_ESFModel_gaussExp;
double massResAccuracy=0.25;

const int scaling=0;
const TString scaling_formula=(scaling==0) ? "sqrt(@0*@1)*@2" : "@2/sqrt(@0*@1)";
const int testDYBins=0;
int plotInitialTrials=0;

#ifdef __myLib__
const int debugFriendlyWP80=0;
const int fit_extra_flag=0;
const int fit_explicit=0;
int fit_dont_fit=0;
const int dump_fit=0;
const int no_wait=1;
std::string histoDir_USER="dir-Fit/";
#endif

// ======================================================================
// ======================================================================

//#define IncludeMSF_C
//#include "fitMassScaleFactorsC.cc"
//#endif

// ======================================================================

//#ifndef __myLib__
int FitMassSmearingIIe(std::ostream &out, unsigned int bin_count, const std::vector<RooDataSet*> &mc_mass, const std::vector<RooDataSet*> &mc_Et1,const std::vector<RooDataSet*> &mc_Et2, const std::vector<RooDataSet*> &exp_orig, RooRealVar &x, RooRealVar &y, const DataInfos_t &info) {
  const char *fncname="FitMassSmearingIIe";
  y.setVal(0.);

  std::cout << "entered FitMassSmearingIIe\n";
  std::cout << "g_esfWorkCase=" << DYTools::g_esfWorkCase << " (" << ElectronEnergyScaleAdv_t::WorkCaseName(DYTools::g_esfWorkCase) << ")\n";
  std::cout << "g_esfFitModel=" << DYTools::g_esfFitModel << " (" << ElectronEnergyScaleAdv_t::FitModelName(DYTools::g_esfFitModel) << ")\n";
  if (out != std::cout) {
    out << "entered FitMassSmearingIIe\n";
    out << "g_esfWorkCase=" << DYTools::g_esfWorkCase << " (" << ElectronEnergyScaleAdv_t::WorkCaseName(DYTools::g_esfWorkCase) << ")\n";
    out << "g_esfFitModel=" << DYTools::g_esfFitModel << " (" << ElectronEnergyScaleAdv_t::FitModelName(DYTools::g_esfFitModel) << ")\n";
  }
  if ((mc_Et1.size()!=0) || (mc_Et2.size()!=0)) { 
    std::cout << "Et-dependence not implemented: mc_Et1.size=" << mc_Et1.size() << ", mc_Et2.size=" << mc_Et2.size() << "\n";
    return 0; 
  }
 
    
  const unsigned int N=mc_mass.size();
  //int usePtDependence=((N==mc_Et1.size() && (N==mc_Et2.size())) ? 1:0;
  if (N!=exp_orig.size()) {
    out << " MC and exp sizes do not match (" << mc_mass.size() << " vs " << exp_orig.size() << ")\n";
    return 0;
  }
  //if (usePtDependence) {
  //  out << "using pt dependence\n";
  //}

  std::vector<TString> idxV;
  std::vector<TString> nameEndings;
  PrepareNameEndings(bin_count,idxV,nameEndings);

  const int reorderVars=1;
  const int excludeSomeSets=0;
  std::vector<int> includeSets;
  includeSets.reserve(exp_orig.size());
  for (unsigned int i=0; i<exp_orig.size(); ++i) includeSets.push_back(1);

  if (excludeSomeSets) {
    std::vector<unsigned int> badI;
    std::vector<unsigned int> badJ;
    badI.push_back(1); badJ.push_back(7);
    badI.push_back(3); badJ.push_back(9);
    unsigned int k=0;
    for (unsigned int i=0; i<bin_count; ++i) {
      for (unsigned int j=i; j<bin_count; ++j, ++k) {
	//out << "k=" << k << ", exp_orig.size()=" << exp_orig.size() << std::endl;
	if (!includeSets[k]) continue;
	for (unsigned int ii=0; ii<badI.size(); ++ii) {
	  if ((i+1==badI[ii]) && (j+1==badJ[ii])) {
	    includeSets[k]=0;
	    out << "forced exclusion: " << nameEndings[k] << " (" << exp_orig[k]->numEntries() << ")\n";
	  }
	}
      }
    }
    out << "1st step passed" << std::endl;
  }

  std::vector<RooDataSet*> exp;
  exp.reserve(exp_orig.size());
  unsigned int k=0;
  for (unsigned int i=0; i<bin_count; ++i) {
    for (unsigned int j=i; j<bin_count; ++j, ++k) {
      if (!includeSets[k]) continue;
      if ((exp_orig[k]->numEntries()<200) && !debugFriendlyWP80 && (i!=j)) includeSets[k]=0;
      if (!includeSets[k]) {
	out << "set " << nameEndings[k] << " has only " << exp_orig[k]->numEntries() << " entries -> will be excluded\n";
      }
      if (includeSets[k] || reorderVars) {
	exp.push_back(exp_orig[k]);
      }
      //out << "i=" << i << ", j=" << j << ", k=" << k << ", set " << nameEndings[k] << " includeSets[k]=" << includeSets[k] << std::endl;
    }
  }


  TString name1;
  TString name2;
  //size_t pos1,pos2;
  char buf[70];

  // Monte-Carlo shape
  RooRealVar mass(x);
  //mass.setBins(10000,"fft");
  mass.setBins(20000,"cache");
  //mass.setRange("myFitRange",60.,120.);
 // Adjust maximum number of steps of RooIntegrator1D in the global default configuration
 
  TCanvas *tstc=NULL;//new TCanvas("tstc","tstc",600,600);
  RooPlot *tstf=(tstc) ? mass.frame("tstf") : NULL;

  // Create a coarser mesh by binning in bins of massResAccuracy
  //const double massResAccuracy=0.25;
  const int massResCount=(testDYBins) ? DYTools::NumberOfMassBins() : int((x.getMax()-x.getMin())/massResAccuracy+1.+1e-3);
  std::cout << "x.range = " << x.getMin() << " .. " << x.getMax() << "\n";
  const double *DYMassBins=DYTools::getMassBins();
  out << " ! massResCount=" << massResCount << "\n";

  TString par1Name="smear_X";
  TString par2Name="par2_X";
  TString par3Name="par3_X";
  switch(DYTools::g_esfFitModel) {
  case _ESFModel_gauss:
    break;
  case _ESFModel_gaussExp:
    par2Name="expTau_X";
    break;
  case _ESFModel_doubleGauss: 
    par1Name="sigma1_X"; par2Name="sigma2_X"; par3Name="relfrac_X";
    break;
  case _ESFModel_bifurGauss:
    par1Name="sigmaL_X"; par2Name="sigmaR_X";
    break;
  case _ESFModel_cball: 
  case _ESFModel_cballShift:
    par2Name="parN_X"; par3Name="alpha_X";
    break;
  case _ESFModel_breitWigner:
    par1Name="gamma_X";
    break;
  case _ESFModel_voigtian:
    par1Name="gamma_X"; par2Name="sigma_X";
    break;
  default:
    out << fncname << ": no parameter name assignment is defined for model " << ElectronEnergyScaleAdv_t::FitModelName(DYTools::g_esfFitModel) << "\n";
    return 0;
  }

  // the scaling factors
  std::vector<RooRealVar*> scaleV;   // scale the template
  { 
    int spec=0;
      if (DYTools::g_esfWorkCase==_ESFWC_6binNegs) {
    //if (DYTools::g_esfFitModel==_ESFModel_gauss) {
	spec=1;
      //}
    }
    if (spec) CreateRooRealVar(bin_count,"scale_X",scaleV,0.99,0.8,1.2);  
    else CreateRooRealVar(bin_count,"scale_X",scaleV,0.99,0.9,1.1);
  }

  // the smearing factors
  RooRealVar bias("bias","bias",0.);
  std::vector<RooRealVar*> smearV;   // smear the template
  double smear_min=(DYTools::g_esfFitModel==_ESFModel_voigtian) ? 0.05 : 0.2;
  if (DYTools::g_esfFitModel==_ESFModel_breitWigner) smear_min=0.1;
  CreateRooRealVar(bin_count,par1Name.Data(),smearV,1.,smear_min,3.0);
  std::vector<RooRealVar*> expTauV; 
  CreateRooRealVar(bin_count,par2Name.Data(),expTauV,1.,0.1,10);

  std::vector<RooRealVar*> alphaV; // CB alpha
  if ((DYTools::g_esfFitModel==_ESFModel_cball) ||
      (DYTools::g_esfFitModel==_ESFModel_cballShift)) {
      CreateRooRealVar(bin_count,par3Name.Data(),alphaV,1.,-10,10);
  }
  else if (DYTools::g_esfFitModel==_ESFModel_doubleGauss) {
    CreateRooRealVar(bin_count,par3Name.Data(),alphaV,0.1,-1,1);
  }


  //smearV[2]->setVal(0.25);
  //smearV[2]->setRange(0.2,0.35);

  if (1) {
    if (DYTools::g_esfWorkCase==_ESFWC_6binNegs) {
      const char *modelName=ElectronEnergyScaleAdv_t::FitModelShortName(DYTools::g_esfFitModel).c_str();
      sprintf(buf,"dir-Work-6bins-%s-fit/testESF_6bins_%s.inp",modelName,modelName);
      ElectronEnergyScaleAdv_t sfx,sfy;
      if (!sfx.LoadASCIIFile(buf) ||
	  !sfy.ChangeToNegs(sfx)) {
	std::cout << "failed to obtain the factors or to change them\n";
	return 0;
      }
      std::cout << "Loaded file " << buf << "\n";
      std::cout << sfx << "\n\n";
      std::cout << " copied to " << sfy << std::endl;
      if (!sfy.CopyScalingFactors(scaleV,NULL) ||
	  !sfy.CopySmearingFactors(smearV,&expTauV,&alphaV)) {
	out << "failed to transfer parameters\n";
	return 0;
      }
      PrintRooVec("set scaleV initial values to ",scaleV);
      PrintRooVec("set smearV initial values to ",smearV);
      // state current settings again (needed for the processing code)
      out << "g_esfWorkCase=" << DYTools::g_esfWorkCase << " (" << ElectronEnergyScaleAdv_t::WorkCaseName(DYTools::g_esfWorkCase) << ")\n";
      out << "g_esfFitModel=" << DYTools::g_esfFitModel << " (" << ElectronEnergyScaleAdv_t::FitModelName(DYTools::g_esfFitModel) << ")" << std::endl;
    }
  }
  else {
    // old style corrections
    if (DYTools::g_esfFitModel==_ESFModel_gauss) {
      if (DYTools::g_esfWorkCase==_ESFWC_5binNegs) {
	scaleV[0]->setVal(1.04); scaleV[1]->setVal(0.999); 
	scaleV[2]->setVal(1.01); scaleV[3]->setVal(1.000);
	scaleV[4]->setVal(1.01); scaleV[5]->setVal(1.000);
	scaleV[6]->setVal(1.01); scaleV[7]->setVal(1.014);
	scaleV[8]->setVal(1.00); scaleV[9]->setVal(1.040);
	smearV[0]->setVal(2.00); smearV[1]->setVal(1.44);
	smearV[2]->setVal(1.11); smearV[3]->setVal(0.61);
	smearV[4]->setVal(0.56); smearV[5]->setVal(0.60);
	smearV[6]->setVal(0.71); smearV[7]->setVal(1.10);
	smearV[8]->setVal(1.50); smearV[9]->setVal(2.00);
	PrintVec("set scaleV initial values to ",scaleV);
	PrintVec("set smearV initial values to ",smearV);
      }
      else if (DYTools::g_esfWorkCase==_ESFWC_6binNegs) {
	if (1) {
	  ElectronEnergyScaleAdv_t sfx,sfy;
	  if (!sfx.LoadASCIIFile("dir-TestESF/testESF_6bins_Gauss.inp") ||
	      !sfy.ChangeToNegs(sfx) ||
	      !sfy.CopyScalingFactors(scaleV,NULL) ||
	      !sfy.CopySmearingFactors(smearV,NULL,NULL)) {
	    out << "failed to transfer parameters\n";
	    return 0;
	  }
	  // state current settings again
	  out << "g_esfWorkCase=" << DYTools::g_esfWorkCase << " (" << ElectronEnergyScaleAdv_t::WorkCaseName(DYTools::g_esfWorkCase) << ")\n";
	  out << "g_esfFitModel=" << DYTools::g_esfFitModel << " (" << ElectronEnergyScaleAdv_t::FitModelName(DYTools::g_esfFitModel) << ")\n";
	}
	else {
	  scaleV[0]->setVal(1.04); scaleV[11]->setVal(1.04);
	  smearV[0]->setVal(2.00); smearV[11]->setVal(2.00);
	  smearV[1]->setVal(1.44); smearV[10]->setVal(1.44);
	  smearV[2]->setVal(1.18); smearV[ 9]->setVal(1.18);
	  smearV[3]->setVal(0.80); smearV[ 8]->setVal(0.80);
	  smearV[4]->setVal(0.60); smearV[ 7]->setVal(0.60);
	  smearV[5]->setVal(0.60); smearV[ 6]->setVal(0.60);
	}
	PrintRooVec("set scaleV initial values to ",scaleV);
	PrintRooVec("set smearV initial values to ",smearV);
      }
    }
  }
 
  RooArgList calcAsymmErr;
  if (fit_extra_flag==1) {
    CreateRooArgList(calcAsymmErr,scaleV,0);
    CreateRooArgList(calcAsymmErr,smearV,1);
  }

 // the scaled masses and the smeared events
  std::vector<RooFormulaVar*> scaledMassV;
  std::vector<RooFormulaVar*> mcSmearV;
  std::vector<RooFormulaVar*> mcExpTauV;
  std::vector<RooFormulaVar*> mcCBAlphaV; // CBall model
  switch(scaling) {
  case 0: out << "scaling sqrt (scale MC mass): " << scaling_formula << "\n"; break;
  case 1: out << "scaling invertedsqrt (scale Exp Et or pT)" << scaling_formula << "\n"; break;
  default:
    out << "scaling=" << scaling << ". Scaling formula description not coded\n";
    return 0;
  }
  k=0;
  for (unsigned int i=0; i<bin_count; ++i) {
    for (unsigned int j=i; j<bin_count; ++j, k++) {
      sprintf(buf,"scaledMass%s",nameEndings[k].Data());
      RooArgList args;
      args.Clear(); args.add(*scaleV[i]); args.add(*scaleV[j]); args.add(mass);
      scaledMassV.push_back(new RooFormulaVar(buf,buf,scaling_formula.Data(),args));

      // smearing parameter (sigma)
      sprintf(buf,"mcSmear%s",nameEndings[k].Data());
      RooFormulaVar *smearingIJ=new RooFormulaVar(buf,"sqrt(@0*@0+@1*@1)",RooArgList(*smearV[i],*smearV[j]));
      mcSmearV.push_back(smearingIJ);

      // decay parameter (tau)
      sprintf(buf,"mcDecay%s",nameEndings[k].Data());
      RooFormulaVar *mctauIJ=new RooFormulaVar(buf,"sqrt(@0*@0+@1*@1)",RooArgList(*expTauV[i],*expTauV[j]));
      mcExpTauV.push_back(mctauIJ);

      if ((DYTools::g_esfFitModel==_ESFModel_cball) ||
	  (DYTools::g_esfFitModel==_ESFModel_cballShift) ||
	  (DYTools::g_esfFitModel==_ESFModel_doubleGauss)) {
	  sprintf(buf,"mcCBAlpha%s",nameEndings[k].Data());
	  RooFormulaVar *alphaIJ=new RooFormulaVar(buf,"sqrt(@0*@0+@1*@1)",RooArgList(*alphaV[i],*alphaV[j]));
	  mcCBAlphaV.push_back(alphaIJ);
      }
    }
  }


  // prepare MC template
  std::vector<mithep::myHistoClass_t*> rawMassH;
  std::vector<RooDataHist*> mcMassH;
  std::vector<RooHistPdf*> mcScaledPdfV;
  if (testDYBins) CreateHistos(N,"rawMassH_X",rawMassH,massResCount,DYMassBins);
  else {
    if (massResCount>1) CreateHistos(N,"rawMassH_X",rawMassH,massResCount,x.getMin(),x.getMax());
  }
  mcMassH.reserve(N);
  for (unsigned int i=0; i<N; ++i) {
    sprintf(buf,"mcMassHisto%s",nameEndings[i].Data());
    if (massResCount>1) {
      mc_mass[i]->fillHistogram(rawMassH[i],x);
      mcMassH.push_back(new RooDataHist(buf,buf,mass,rawMassH[i]));
    }
    else {
      mcMassH.push_back(mc_mass[i]->binnedClone(buf));
    }
    sprintf(buf,"mcScaledPdf%s",nameEndings[i].Data());
    mcScaledPdfV.push_back(new RooHistPdf(buf,buf,*scaledMassV[i],mass,*mcMassH.back()));
  }

  // One variable to govern MC/data ratio
  RooRealVar MCOverData("MCOverData","MC_over_data",0.1,1e-6,1e10);

  // the resolution model
  std::vector<RooGExpModel*> gaussExpV;
  std::vector<RooGaussModel*> gaussV;
  std::vector<RooGaussModel*> gauss2V;
  std::vector<RooAddModel*> gaussSumV;
  std::vector<RooBifurGauss*> bifurGaussV;
  std::vector<RooCBShape*> cbShapeV;
  std::vector<RooCBShape*> cbShapeShiftedV;
  std::vector<RooFormulaVar*> cbssVarsV;
  //std::vector<RooMoment*> calcMeanV;
  std::vector<RooFormulaVar*> shiftFromMeanV;
  std::vector<RooBreitWigner*> bwShapeV;
  std::vector<RooVoigtian*> voigtianV;
  //std::vector<RooRealVar*> 
  switch (DYTools::g_esfFitModel) {
  case _ESFModel_gauss:
    //CreateRooGaussian(N,"gauss_X",mass,bias,mcSmearV, gaussV);
    CreateRooGaussian(N,"gauss_X",mass,bias,mcSmearV, gaussV);
    break;
  case _ESFModel_gaussExp:
    CreateRooGaussExp(N,"gaussExp_X",mass,bias,mcSmearV,mcExpTauV, gaussExpV);
    break;
  case _ESFModel_doubleGauss: {
    CreateRooGaussian(N,"gauss1_X",mass,bias,mcSmearV, gaussV);
    CreateRooGaussian(N,"gauss2_X",mass,bias,mcExpTauV, gauss2V);
    std::vector<RooArgList*> list;
    CreateRooArgList(gaussV,list);
    CreateRooArgList(gauss2V,list,1);
    if (0) {
      std::vector<RooArgList*> coefList;
      CreateRooArgList(mcCBAlphaV,coefList);
      CreateRooAddModel(N,"gaussSum_X",nameEndings,list,coefList, gaussSumV);
    }
    else CreateRooAddModel(N,"gaussSum_X",nameEndings,list,mcCBAlphaV, gaussSumV);
  }
    break;
  case _ESFModel_bifurGauss: {
    CreateRooBifurGauss(N,"bifurGauss_X",mass,bias,mcSmearV,mcExpTauV, bifurGaussV);
    break;
  }
  case _ESFModel_cball: {
    CreateRooCBShape(N,"cbShape_X",mass,bias,mcSmearV,mcCBAlphaV,mcExpTauV, cbShapeV);
  }
    break;
  case _ESFModel_cballShift: {
    char formula[100];
    //CreateRooCBShape(N,"cbShapeShifted_X",mass,bias,mcSmearV,mcCBAlphaV,mcExpTauV, cbShapeShiftedV);
    cbssVarsV.reserve(12*N);
    for (unsigned int i=0; i<N; ++i) {
      sprintf(buf,"cbvar_nOverAlpha_%s",nameEndings[i].Data());
      RooFormulaVar *vR=new RooFormulaVar(buf,"@0/@1",RooArgList(*mcExpTauV[i],*mcCBAlphaV[i]));
      sprintf(buf,"cbvar_B_%s",nameEndings[i].Data());
      RooFormulaVar *vB=new RooFormulaVar(buf,"TMath::Power(@0/@1,@1)-@0",RooArgList(*mcExpTauV[i],*mcCBAlphaV[i]));
      const double infty=200.;
      const double initialMean=5.00203;
      sprintf(buf,"cbvar_z20_%s",nameEndings[i].Data());
      sprintf(formula,"@0/(@1+%3.0lf/@2)",infty);
      RooFormulaVar *vz20=new RooFormulaVar(buf,formula,RooArgList(*vR,*vB,*mcSmearV[i]));
      sprintf(buf,"cbvar_z21_%s",nameEndings[i].Data());
      sprintf(formula,"@0/(@1+(%3.0lf+%7.5lf)/@2)",infty,initialMean);
      RooFormulaVar *vz21=new RooFormulaVar(buf,formula,RooArgList(*vR,*vB,*mcSmearV[i]));
      
      sprintf(buf,"cbvar_z20n2_%s",nameEndings[i].Data());
      RooFormulaVar *vz20n2=new RooFormulaVar(buf,"TMath::Power(@0,@1-2)",RooArgList(*vz20,*mcExpTauV[i]));
      sprintf(buf,"cbvar_z20n1_%s",nameEndings[i].Data());
      RooFormulaVar *vz20n1=new RooFormulaVar(buf,"TMath::Power(@0,@1-1)",RooArgList(*vz20,*mcExpTauV[i]));
      sprintf(buf,"cbvar_z21n2_%s",nameEndings[i].Data());
      RooFormulaVar *vz21n2=new RooFormulaVar(buf,"TMath::Power(@0,@1-2)",RooArgList(*vz21,*mcExpTauV[i]));
      sprintf(buf,"cbvar_z21n1_%s",nameEndings[i].Data());
      RooFormulaVar *vz21n1=new RooFormulaVar(buf,"TMath::Power(@0,@1-1)",RooArgList(*vz21,*mcExpTauV[i]));
      sprintf(buf,"cbvar_nom_%s",nameEndings[i].Data());
      sprintf(formula,"1 + @0*@1/(@2-1)*(1-@3) - @1*@1/(@2-2)*(1-@4)");
      RooFormulaVar *vNom=new RooFormulaVar(buf,formula,RooArgList(*vB,*vR,*mcExpTauV[i],*vz20n1,*vz20n2));
      sprintf(buf,"cbvar_denomExp_%s",nameEndings[i].Data());
      sprintf(formula,"%12.9e*exp(0.5*@0*@0)*(RooMath::erf(%12.9e*@0)-1)",sqrt(0.5*4*atan(1)),-1/sqrt(2.));
      RooFormulaVar *vDenomExp = new RooFormulaVar(buf,formula,RooArgList(*mcCBAlphaV[i]));
      sprintf(buf,"cbvar_denomPow1_%s",nameEndings[i].Data());
      RooFormulaVar *vDenomPow1=new RooFormulaVar(buf,"@0/(1-@1)*(1-@2)",RooArgList(*vR,*mcExpTauV[i],*vz20n1));
      sprintf(buf,"cbvar_denomPow2_%s",nameEndings[i].Data());
      sprintf(formula,"%12.9e*@0*(@1*@2/(1-@3)*(@4-@5) - @2*@2/(2-@3)*(@6-@7))",1/initialMean);
      RooFormulaVar *vDenomPow2=new RooFormulaVar(buf,formula,RooArgList(*mcSmearV[i],*vB,*vR,*mcExpTauV[i],*vz20n1,*vz21n1,*vz20n2,*vz21n2));
      sprintf(buf,"shiftFromMean%s",nameEndings[i].Data());
      RooFormulaVar *var=new RooFormulaVar(buf,"@0*@1/(@2+@3+@4)",RooArgList(*mcSmearV[i],*vNom,*vDenomExp,*vDenomPow1,*vDenomPow2));
      shiftFromMeanV.push_back(var);

      unsigned int idx0=cbssVarsV.size();
      cbssVarsV.push_back(vR);
      cbssVarsV.push_back(vB);
      cbssVarsV.push_back(vz20);
      cbssVarsV.push_back(vz21);
      cbssVarsV.push_back(vz20n2);
      cbssVarsV.push_back(vz20n1);
      cbssVarsV.push_back(vz21n2);
      cbssVarsV.push_back(vz21n1);
      cbssVarsV.push_back(vNom);
      cbssVarsV.push_back(vDenomExp);
      cbssVarsV.push_back(vDenomPow1);
      cbssVarsV.push_back(vDenomPow2);
      if (1) {
	for (unsigned int idx=idx0; idx<cbssVarsV.size(); ++idx) {
	  cbssVarsV[idx]->Print();
	}
	var->Print();
      }
    }
    CreateRooCBShape(N,"cbShape_X",mass,shiftFromMeanV,mcSmearV,mcCBAlphaV,mcExpTauV, cbShapeV);
    //for (unsigned int i=0; i<cbShapeV.size(); ++i) {
    //  std::cout << "cbShape[" << i << "] mean = " << cbShapeV[i]->mean(mass)->getVal() << "\n";
    //}
  }
    break;
  case _ESFModel_breitWigner:
    CreateRooBreitWigner(N,"bw_X",mass,bias,mcSmearV, bwShapeV);
    break;
  case _ESFModel_voigtian:
    CreateRooVoigtian(N,"voigtian_X",mass,bias,mcSmearV,mcExpTauV, voigtianV);
    break;
  default:
    out << fncname << ": no resolution model construction defined for model " << ElectronEnergyScaleAdv_t::FitModelName(DYTools::g_esfFitModel) << "\n";
    return 0;
  }

  //Pause("before conv");
  std::vector<RooFFTConvPdf*> mcSignalV;
  switch(DYTools::g_esfFitModel) {
  case _ESFModel_gauss:
    CreateRooConvPdf(N,"mcSignal_X",mass,mcScaledPdfV,gaussV, mcSignalV);
    break;
  case _ESFModel_gaussExp:
    CreateRooConvPdf(N,"mcSignal_X",mass,mcScaledPdfV,gaussExpV, mcSignalV);
    break;
  case _ESFModel_doubleGauss:
    CreateRooConvPdf(N,"mcSignal_X",mass,mcScaledPdfV,gaussSumV, mcSignalV);
    break;   
  case _ESFModel_bifurGauss:
    CreateRooConvPdf(N,"mcSignal_X",mass,mcScaledPdfV,bifurGaussV, mcSignalV);
    break;
  case _ESFModel_cball:
  case _ESFModel_cballShift:
    CreateRooConvPdf(N,"mcSignal_X",mass,mcScaledPdfV,cbShapeV, mcSignalV);
    break;
  case _ESFModel_breitWigner:
    CreateRooConvPdf(N,"mcSignal_X",mass,mcScaledPdfV,bwShapeV, mcSignalV);
    break;
  case _ESFModel_voigtian:
    CreateRooConvPdf(N,"mcSignal_X",mass,mcScaledPdfV,voigtianV, mcSignalV);
    break;
  default:
    out << fncname << ": convolution not defined for model " << ElectronEnergyScaleAdv_t::FitModelName(DYTools::g_esfFitModel) << "\n";
    return 0;
  }
  

  if (1 && tstc) {
    out << "plot it\n";
    mcSignalV[0]->plotOn(tstf,RooFit::LineColor(kRed+2));
    out << "draw test frame\n";
    tstf->Draw();
    tstc->Update();
    Pause("chk");
  }

  if (0) {
    //PrintValues("mcMassV",mcMassV);
    PrintValues("scaleV",scaleV);
    PrintValues("smearV",smearV);
    PrintValues("mcSmearV",mcSmearV);
    //PrintValues("mcEvPosV",mass,mcEvPosV,30,60,120);
    //PrintValues("mcEventsV",mcEventsV);
    //PrintPrint("mcEventsV",mcEventsV);
  }

  // Create the extended contributions to be able to add the events
  std::vector<RooExtendPdf*> mcExtEventsV;
  std::vector<RooArgList*> evListV;
  for (unsigned int i=0; i<mcSignalV.size(); ++i) {
    sprintf(buf,"mcExtEvents%s",nameEndings[i].Data());
    RooExtendPdf *ext=new RooExtendPdf(buf,buf,*mcSignalV[i],MCOverData);
    mcExtEventsV.push_back(ext);
    sprintf(buf,"evList_%s",nameEndings[i].Data());
    RooArgList *arg=new RooArgList(buf);
    arg->add(*ext);
    evListV.push_back(arg);
  }
  
  out << "template extended\n";
  out << "evListV.size=" << evListV.size() << "\n";

  // Combine the models
  std::vector<RooAddPdf*> model;
  CreateRooAddPdf(N,"model",nameEndings,evListV,model);
  std::vector<RooAddPdf*> modelSave; // needed for modified selections
  modelSave=model;
 
  // Define category to distinguish physics and control sample events
  std::vector<TString> categories;
  categories.reserve(N);
  RooCategory sample("sample","sample");
  for (unsigned int i=0; i<N; ++i) {
    if (includeSets[i] || reorderVars) {
      name1="set"; name1.Append(idxV[i]);
      categories.push_back(name1);
      sample.defineType(name1);
    }
  }
  //out << "categories: "; for (unsigned int i=0; i<categories.size(); ++i) out << " " << categories[i]; out << "\n";
  //out << "counts: exp.size=" << exp.size() << ", categories=" << categories.size() << std::endl; //, sample.size=" << sample.numEntries() << "\n";
  //out << "N=" << N << ", includeSets.size=" << includeSets.size() << ", model.size=" << model.size() << std::endl;

  if (reorderVars) {
    std::vector<RooDataSet*> eTmp;
    std::vector<TString> neTmp;
    std::vector<TString> caTmp;
    std::vector<int> includeSetsNew;
    std::vector<RooAddPdf*> modTmp;
    neTmp=nameEndings; nameEndings.clear();
    eTmp=exp; exp.clear();
    caTmp=categories; categories.clear();
    sample.clearTypes();
    modTmp=model; model.clear(); modelSave.clear();
    for (unsigned int dij=0; dij<bin_count; ++dij) {
      k=0;
      for (unsigned int i=0; i<bin_count; ++i) {
	for (unsigned int j=i; j<bin_count; ++j, ++k) {
	  if (j-i==dij) {
	    //out << "i=" << i << ", j=" << j << ", k=" << k << std::endl;
	    nameEndings.push_back(neTmp[k]);
	    includeSetsNew.push_back(includeSets[k]);
	    model.push_back(modTmp[k]);
	    if (includeSets[k]) {
	      exp.push_back(eTmp[k]);
	      categories.push_back(caTmp[k]);
	      sample.defineType(caTmp[k]);
	      modelSave.push_back(modTmp[k]);
	    }
	    //out << "passed" << std::endl;
	  }
	}
      }
    }
    includeSets=includeSetsNew;

    out << "reordered: \n";
    //PrintVec("nameEndings",nameEndings);
    //PrintVec("includeSets",includeSets);
    PrintVec2("nameEndings & includeSets",nameEndings,includeSets);
    //return 0;
  }

  out << "categories: "; for (unsigned int i=0; i<categories.size(); ++i) out << " " << categories[i]; out << "\n";
  out << "counts: exp.size=" << exp.size() << ", categories=" << categories.size() << std::endl; //, sample.size=" << sample.numEntries() << "\n";

  /*
  std::map<std::string,RooDataSet*> sets_by_caths;
  //sets_by_caths.reserve(N);
  for (unsigned int i=0; i<N; ++i) {
    sets_by_caths.insert(std::pair<std::string,RooDataSet*>(std::string(categories[i]),exp[i]));
  }
  RooDataSet *combData=new RooDataSet("combData","combined data", mass, RooFit::Index(sample),RooFit::Import(sets_by_caths));
  */

  RooDataSet *combData=CreateCombinedData("combData","combined data", exp, mass,sample,categories);
  out << "combData created" << std::endl;

  // Construct a simultaneous pdf
  RooSimultaneous *simPdf=new RooSimultaneous("simPdf","simultaneous pdf",sample);
  k=0;
  for (unsigned int i=0; i<N; ++i) {
    if (includeSets[i]) {
      simPdf->addPdf(*model[i],categories[k]);
      k++;
    }
  }
  out << "simPdf has " << k << " models" << std::endl;


  // fitting

  TString frame1Name="initial frame";
  if (info.NameOk()) {
    frame1Name.Append(" for ");
    frame1Name.Append(info.getName());
  }

  vector<RooPlot*> frames;
  RooPlot* frame = (no_plot) ? NULL : x.frame(RooFit::Title(frame1Name) );
  if (frame) {
    // Plot all data according to the category
    for (unsigned int i=0; i<categories.size(); i++) {
      name1="sample==sample::" + categories[i];
      out << "plot combData with name=" << name1 << "\n";
      combData->plotOn(frame,RooFit::Cut(name1), RooFit::LineColor(info.color(i)),RooFit::MarkerColor(info.color(i)));
    }
    //dataV[0]->plotOn(frame,RooFit::LineColor(kGreen+2),RooFit::LineStyle(kDashed));
    //if (N>1) dataV[1]->plotOn(frame,RooFit::LineColor(kRed),RooFit::LineStyle(kDashed));
  }

  if (frame && plotInitialTrials) {
    // plot initial trials
    for (unsigned int i=0; i<categories.size(); i++) {
      simPdf->plotOn(frame,RooFit::Slice(sample,categories[i]),RooFit::ProjWData(sample,*combData), RooFit::LineColor(info.color(i)));
    }
    frames.push_back(frame);
  }
  RooPlot *frame2 = NULL;
  TString frame2Name;

  if (fit_explicit) out << "the request to fit explicitly is ignored!\n";

  if (0) {
    TCanvas *tstc2=new TCanvas("tstc2","tstc2",600,600);
    RooPlot *tstf2=(tstc2) ? mass.frame(RooFit::Title("tstf2")) : NULL;
    mcMassH.back()->plotOn(tstf2,RooFit::LineColor(kBlack));
    mcScaledPdfV.back()->plotOn(tstf2,RooFit::LineColor(kRed+2));
    mcSignalV.back()->plotOn(tstf2,RooFit::LineColor(kGreen+2));
    exp.back()->plotOn(tstf2,RooFit::LineColor(kBlue+2));
    tstf2->Draw();
    tstc2->Update();
    tstc2->SaveAs("tstc2.gif");
    Pause();
    delete tstc2;
    delete tstf2;
    //return 0;
  }
 
  RooFitResult *fitRes=NULL;
  if (!fit_dont_fit) {
    fitRes=simPdf->fitTo(*combData,RooFit::Extended(true),RooFit::NumCPU(2,true), RooFit::Timer(true),RooFit::Save()); //,RooFit::Range("myFitRange"));
    if (fit_extra_flag==1) {
      out << " ** fit_extra_flag=1. Calling the fit again, to get asymm.errors\n";
      fitRes=simPdf->fitTo(*combData,RooFit::Extended(true),RooFit::NumCPU(2,true), RooFit::Timer(true),RooFit::Save(),RooFit::Minos(calcAsymmErr)); //,RooFit::Range("myFitRange"));
    }
  }

  if (0) {
    TCanvas *tstc3=new TCanvas("tstc3","tstc3",600,600);
    RooPlot *tstf3=(tstc3) ? mass.frame(RooFit::Title("tstf3")) : NULL;
    mcMassH.back()->plotOn(tstf3,RooFit::LineColor(kBlack));
    mcScaledPdfV.back()->plotOn(tstf3,RooFit::LineColor(kRed+2));
    mcSignalV.back()->plotOn(tstf3,RooFit::LineColor(kGreen+2));
    exp.back()->plotOn(tstf3,RooFit::LineColor(kBlue+2));
    tstf3->Draw();
    tstc3->Update();
    tstc3->SaveAs("tstc3.gif");
    Pause();
    delete tstc3;
    delete tstf3;
    return 0;
  }

  frame2Name = "implicit fit";

  if (info.NameOk()) {
    frame2Name.Append(" for ");
    frame2Name.Append(info.getName());
  }
  frame2 = (no_plot) ? NULL : x.frame(RooFit::Title(frame2Name));
    
  // Plot all data according to the category
  if (!no_plot) {
    if (0) {
      for (unsigned int i=0; i<categories.size(); i++) {
	name1="sample==sample::" + categories[i];
	combData->plotOn(frame2,RooFit::Cut(name1), RooFit::LineColor(info.color(i)));
      }
      // plot initial trials
      for (unsigned int i=0; i<categories.size(); i++) {
	simPdf->plotOn(frame2,RooFit::Slice(sample,categories[i]),RooFit::ProjWData(sample,*combData), RooFit::LineColor(info.color(i)),RooFit::MarkerColor(info.color(i)));
      }
    }
    else {
      for (unsigned int i=0; i<categories.size(); i++) {
	frame2Name=categories[i];
	RooPlot *framex= x.frame(RooFit::Title(frame2Name));
	name1="sample==sample::" + categories[i];
	out << "\n plot " << name1 << "\n";
	combData->plotOn(framex,RooFit::Cut(name1), RooFit::LineColor(info.color(i)));
	//combData->plotOn(framex,RooFit::Cut(name1), RooFit::LineColor(info.color(i)),RooFit::MarkerColor(info.color(i)));
	simPdf->plotOn(framex,RooFit::Slice(sample,categories[i]),RooFit::ProjWData(sample,*combData), RooFit::LineColor(info.color(i)));
	//ebkgV[i]->plotOn(framex,RooFit::LineColor(kGreen+2), RooFit::LineStyle(kDashed));
	//for (unsigned int mi=0; (mi<i) && (mi<i); ++mi) {
	if (0) {
	  //unsigned int mi=i;
	  //name1="esig_X"; pos1=name1.Length()-1;
	  //name2="ebkgr_X"; pos2=name2.Length()-1;
	  name1="esig_"; name1.Append(idxV[i]);
	  //name2="ebkgr_"; name2.Append(idxV[i]);
	  model[i]->plotOn(framex,RooFit::Components(name1.Data()),RooFit::LineStyle(kDashed),RooFit::Normalization(1.,RooAbsReal::RelativeExpected));
	  //model[i]->plotOn(framex,RooFit::Components(name2.Data()),RooFit::LineStyle(kDotted),RooFit::Normalization(1.,RooAbsReal::RelativeExpected));
	  //}
	}
	frames.push_back(framex);
      }
    }
  }

  if (0 && frame2) {
    //model.plotOn(frame2,RooFit::Components(esig),RooFit::LineColor(kGreen+2),RooFit::LineStyle(kDashed));
    //model.plotOn(frame2,RooFit::Components(ebkg),RooFit::LineColor(kRed+2),RooFit::LineStyle(kDotted));
    //esig.plotOn(frame2,RooFit::LineColor(kGreen+2));
    //ebkg.plotOn(frame2,RooFit::LineColor(kRed+2));
  }
  
  TString canvas_name="BWcb";
  TCanvas* c4 = NULL;
  if (!no_plot) {
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
  }
  const char *outfile_name="fit_mscf.png";
  //const char *outfile_name="fit_mscf.C";

  if (0 && c4) {
    c4->Update();
    c4->SaveAs(outfile_name);
    //c4->SaveAs("fit_frames.C");
  }

  if (histoDir_USER.size()) {
    out << "PREPARING HTML\n" << std::endl;
    ElectronEnergyScaleAdv_t sf;
    sf.WorkCase(DYTools::g_esfWorkCase,DYTools::g_esfFitModel);
    std::vector<RooRealVar*> scPars3;
    if (DYTools::g_esfFitModel==_ESFModel_cballShift) {
      for (unsigned int i=0; i<shiftFromMeanV.size(); ++i) {
	sprintf(buf,"shift%s",nameEndings[i].Data());
	RooRealVar *v=new RooRealVar(buf,buf,shiftFromMeanV[i]->getVal());
	v->setError(0); 
	scPars3.push_back(v);
      }
    }
    if (!sf.Assign(DYTools::g_esfFitModel,MCOverData,scaleV,scaling,smearV,&expTauV,&alphaV,&scPars3)) {
      out << "WARNING: failed to assign sf properly" << std::endl;
    }
    out << "Assigned sf : " << sf << "\n";
    if (scaling==0) { sf.InvertedSF(3); out << "scaling sqrt\n"; }
    //return 1;

    if (1) {
      out << "saving fit parameters to a separate file\n";
      char path[200];
      char locFName[300];
      sprintf(path,"%s",histoDir_USER.c_str());
      char *p=strchr(path,'/');
      *(p+1)='\0';
      sprintf(locFName,"%stestESF_%s_%s.inp",path,sf.WorkCaseShortName().c_str(),sf.FitModelShortName().c_str());
      if (scaling==0) sf.InvertedSF(2);
      sf.SaveASCIIFile(locFName);
      out << "fit parameters saved to a file: " << locFName << "\n";
    }


    std::vector<std::string> html_lines;   
    std::vector<TString> hr_nameEndings;
    sf.PrintValuesHTMLTable(html_lines);
    char buf2[250];
    sprintf(buf2,"<tr><td width=""25%%"">mass bin size (massResAccuracy) in fit=%5.2lf</td></tr>\n",massResAccuracy);
    html_lines.push_back(buf2);
    if (reorderVars) {
      if (exp.size()!=exp_orig.size()) {
	sprintf(buf2,"<tr><td width=""25%%"">WARNING: The plot combinations<br> may be incorrect,<br> since some ranges were omitted.<br> Fit was done to %d cases<br> out of %d</td></tr>\n",int(exp.size()),int(exp_orig.size()));
	html_lines.push_back(buf2);
      }
    }
    sf.PrepareHRNameEndings(hr_nameEndings);
    std::vector<std::vector<int>*> combine;
    std::vector<TString> combination_names;
    PrepareCombinationInfos(sf,_ESFWC_all, combine, combination_names);

    //out << "MAKE HTML\n";
    if (reorderVars) {
      makeHTML(histoDir_USER.c_str(),"fitMassScaleFactors.html", "fitMassScaleFactors5e inspection plot", html_lines, hr_nameEndings,nameEndings, mass,exp,modelSave, &combine, &combination_names);    
    }
    else {
      makeHTML(histoDir_USER.c_str(),"fitMassScaleFactors.html", "fitMassScaleFactors5e inspection plot", html_lines, hr_nameEndings,nameEndings, mass,exp,model, &combine, &combination_names);    
    }
    //out << "MAKE HTML returned\n";
  }

  if (dump_fit) {
    out << "--------------------------------------------------\n";
    //simPdf->Print("t");
    simPdf->Print();
    out << "Dumping fit parameters\n";
    for (unsigned int i=0; i<model.size(); ++i) {
      out << "model #" << i+1 << " : ";// << model[i]->getVal() << "\n";
      model[i]->Print(); //out << "\n";
      //out << "error: " << model[i]->getError() << "\n";
    }
    //for (unsigned int i=0; i<mcEventsV.size(); ++i) {
    //  out << "esigV #" << i+1 << " : ";
    //  esigV[i]->Print(); //out << "\n";
    //}
    /*
    out << "derived signal event counts:\n";
    for (unsigned int i=0; i<passf.size(); ++i) {
      out << " " << (i+1) << " nsignal_" << (i+1) << "   " << passf[i]->getVal() << " ";
      if (fitRes) out << passf[i]->getPropagatedError(*fitRes);
      out << "\n";
    }
    if (0) {
      out << "another dump\n";
      for (unsigned int i=0; i<model.size(); ++i) {
	out << "model #" << i+1 << " : ";// << model[i]->getVal() << "\n";
	model[i]->printCompactTree(); //out << "\n";
	//out << "error: " << model[i]->getError() << "\n";
      }
    }
    out << "Efficiency errors :\n";
    for (unsigned int i=0; i<eff.size(); ++i) {
      out << eff[i]->GetName() << " = " << eff[i]->getVal() << " +" << eff[i]->getErrorHi() << " " << eff[i]->getErrorLo() << "\n";
    }
    out << "Event number errors :\n";
    for (unsigned int i=0; i<passf.size(); ++i) {
      out << passf[i]->GetName() << " = " << passf[i]->getVal() << "  +/-" << ((fitRes) ? passf[i]->getPropagatedError(*fitRes) : 0.0) << "\n";
      }
    */
    out << "--------------------------------------------------\n";
  }

  //if (ctest) { ctest->Update(); no_wait=0; }
  if (!no_plot) out << "plot saved as " << outfile_name << "\n";
    if (!no_wait) { out << "updated c4 : enter a char...\n"; char xxx; std::cin >> xxx; }
    /*  */
  return 1;
}
//#endif

// ======================================================================

#ifndef __myLib__
int PrepareScaledSmearedData(const ElectronEnergyScaleAdv_t &sf, const std::vector<std::vector<double>*> &dataMC, const std::vector<std::vector<double>*> &dataExp, const std::vector<std::vector<double>*> &expEt1V, const std::vector<std::vector<double>*> &expEt2V,  RooRealVar &mass, ScaledSmearedDataAux_t &aux, std::vector<RooAddPdf*> &mcModelV, std::vector<RooDataSet*> &expModelV, int scaleMC) 
#else
int PrepareScaledSmearedData(const ElectronEnergyScaleAdv_t &sf, const std::vector<std::vector<double>*> &dataMC, const std::vector<std::vector<double>*> &dataExp, const std::vector<std::vector<double>*> &expEt1V, const std::vector<std::vector<double>*> &expEt2V,  RooRealVar &mass, ScaledSmearedDataAux_t &aux, std::vector<RooAddPdf*> &mcModelV, std::vector<RooDataSet*> &expModelV) 
#endif
{
  const char *fncname="PrepareScaledSmearedData";
 std::cout << "Entered PrepareScaledSmearedData\n";

 PrintVVecCounts("dataMC ",dataMC);
 PrintVVecCounts("dataExp ",dataExp);

  if (0) { PrintVVecCounts("expEt1V ",expEt1V); PrintVVecCounts("expEt2V ",expEt2V); }

 if (0) {  // for sample code
   // Prepare MC data
   int bin_count=sf.EtaDivisionCount();
   //int expect_count=bin_count*(bin_count+1)/2;
   int largest=0;
   unsigned int largest_count=0;
   unsigned int k=0;
   char signature[50];
   for (int i=0; i<bin_count; ++i) {
     for (int j=i; j<bin_count; ++j, ++k) {
       if (dataExp[k]->size()>largest_count) {
	 largest_count=dataExp[k]->size();
	 largest=k;
	 sprintf(signature,"%d_%d__%d",i,j,k);
	 std::cout << "assigning largest is (i,j)=(" << i << ',' << j << "), k=" << k << "\n";
       }
     }
   }
   char buf[200];
   std::cout << "signature=" << signature << "\n";
   sprintf(buf,"sampleData_Exp_%s.dat",signature);
   FILE *fout=fopen(buf,"w");
   const std::vector<double> *d=dataExp[largest];
   fprintf(fout,"%u   ",(unsigned int)(d->size()));
   for (unsigned int i=0; i<d->size(); ++i) {
     fprintf(fout," %12.8lf",(*d)[i]);
   }
   fclose(fout);
   std::cout << "file <" << buf << "> saved\n";
   sprintf(buf,"sampleData_MC_%d.dat",largest);
   FILE *fout2=fopen(buf,"w");
   const std::vector<double> *d2=dataMC[largest];
   fprintf(fout2,"%u  ",(unsigned int)(d2->size()));
   for (unsigned int i=0; i<d2->size(); ++i) {
     fprintf(fout2," %12.8lf",(*d2)[i]);
   }
   fclose(fout2);
   std::cout << "file <" << buf << "> saved\n";
   std::cout << "forcing a return 0\n";
   return 0;
 }
 
  // scale experimental values
  std::cout << "\nwill use sf=" << sf << "\n";
  std::vector<std::vector<double>*> *dataExp_local=new std::vector<std::vector<double>*>();
  for (unsigned int i=0; i<dataExp.size(); ++i) {
    dataExp_local -> push_back(new std::vector<double>(*dataExp[i]));
  }
#ifndef __myLib__
  if (!scaleMC) 
#endif
    sf.ScaleValues(*dataExp_local);

  std::vector<std::vector<double>*> *dataMC_local=new std::vector<std::vector<double>*>();
  for (unsigned int i=0; i<dataMC.size(); ++i) {
    dataMC_local -> push_back(new std::vector<double>(*dataMC[i]));
  }
#ifndef __myLib__
  if (scaleMC) sf.InvScaleValues(*dataMC_local);
#endif

  {
    //std::vector<std::vector<double>*> emptyVec;
    RooRealVar y;
    RooRealVar x(mass);
    std::vector<RooDataSet*> *tmpExpModelV=ConvertToDataSet(*dataExp_local,x,y,0);
    expModelV.insert(expModelV.end(),tmpExpModelV->begin(),tmpExpModelV->end());
    delete tmpExpModelV;
  }
  delete dataExp_local;

  std::cout << "exp data is ready (" << expModelV.size() << " entries)\n";
  PrintRooDataSetCounts("dataExp ",expModelV);


  // Prepare MC data
  int bin_count=sf.EtaDivisionCount();
  int expect_count=bin_count*(bin_count+1)/2;

  //RooRealVar mass(x);
  //mass.setBins(10000,"fft");

  vector<RooPlot*> frames;
  RooPlot* frameX=mass.frame(RooFit::Title("frameX"));
  //frames.push_back(frameX);
  
  std::vector<char> idxV;
  std::vector<TString> nameEndings;
  PrepareNameEndings(bin_count,idxV,nameEndings);

  // Create a coarser mesh by binning in bins of massResAccuracy
  //use global value // const double massResAccuracy=0.2; 
  const int massResCount=int((mass.getMax()-mass.getMin())/massResAccuracy+1.+1e-3);
  std::cout << " ! massResCount=" << massResCount << "\n";
  //std::vector<TH1F*> massH;
  //#ifdef __myLib__
  //aux.massH.reserve(dataMC.size());
  //char buf[200];
  //for (unsigned int i=0; i<dataMC->size(); ++i) {
  //  sprintf(buf,"mcMassH_%s",nameEndings[i].Data());
  //  aux.massH[i] = (mithep::myHistoClass_t*)dataMC[i]->Clone(buf);
  //}
  //#else
  CreateHistos(expect_count,"mcMassH_X",aux.massH,massResCount,mass.getMin(),mass.getMax());
  //#ifndef __myLib__
  //for (unsigned int i=0; i<dataMC_local->size(); ++i) {
  //  (*dataMC_local)[i]->FillHistogram(0,aux.massH[i]);
  //}
  //#else
  for (unsigned int i=0; i<dataMC_local->size(); ++i) {
    const std::vector<double> *d=(*dataMC_local)[i];
    for (unsigned int ii=0; ii<d->size(); ++ii) {
      aux.massH[i]->Fill((*d)[ii]);
    }
  }
  //#endif
  //#endif

  // create templates
  //std::vector<RooDataHist*> mcHistV;
  //std::vector<RooHistPdf*> mcV;
  CreateRooHistPdf(mass,aux.massH,"MCHist_X",aux.mcHistV,"MCShape_X",aux.mcV);
  std::cout << "mcV created (" << aux.mcV.size() << " entries)" << std::endl;
  //mcHistV[0]->plotOn(frameX);
  aux.mcV[0]->plotOn(frameX);

  // One variable to govern MC/data ratio
  //RooRealVar MCOverData("MCOverData","MC_over_data",sf.MCtoData());
  
  // the bias and the smearing factors
  //RooRealVar bias("bias","bias",0);
  //std::vector<RooRealVar*> gWidthV;   // smear the template
  CreateRooRealVar(expect_count,"gWidth",nameEndings,aux.gWidthV,0);

  switch(sf.FitModel()) {
  case _ESFModel_gauss:
    sf.TransferSmearingFactors(aux.gWidthV,NULL,NULL);
    //std::cout << "got smearing widths: "; PrintRooVec("aux.gWidthV=",aux.gWidthV);
    // smearing Gaussians
    //std::vector<RooGaussModel*> gaussV;
    CreateRooGaussian(expect_count,"gauss_X",mass,aux.bias,aux.gWidthV, aux.gaussV);
    //std::vector<RooFFTConvPdf*> mcSignalV;
    CreateRooConvPdf(expect_count,"mcSignal_X",mass,aux.mcV,aux.gaussV, aux.mcSignalV);
    break;
  case _ESFModel_gaussExp: 
    CreateRooRealVar(expect_count,"expTau_X",nameEndings,aux.smPars2V,0);
    sf.TransferSmearingFactors(aux.gWidthV, &aux.smPars2V,NULL);
    CreateRooGaussExp(expect_count,"gaussExp_X",mass,aux.bias,aux.gWidthV,aux.smPars2V, aux.gaussExpV);
    CreateRooConvPdf(expect_count,"mcSignal_X",mass,aux.mcV,aux.gaussExpV, aux.mcSignalV);
    break;
  case _ESFModel_doubleGauss: {
    CreateRooRealVar(expect_count,"sigma2_X",nameEndings,aux.smPars2V,0);
    CreateRooRealVar(expect_count,"relfrac_X",nameEndings,aux.smPars3V,0);
    sf.TransferSmearingFactors(aux.gWidthV, &aux.smPars2V, &aux.smPars3V);
    CreateRooGaussian(expect_count,"gauss1_X",mass,aux.bias,aux.gWidthV, aux.gaussV);
    CreateRooGaussian(expect_count,"gauss2_X",mass,aux.bias,aux.smPars2V, aux.gauss2V);
    std::vector<RooArgList*> list, coefList;
    CreateRooArgList(aux.gaussV,list); CreateRooArgList(aux.gauss2V,list,1);
    CreateRooArgList(aux.smPars3V,coefList);
    CreateRooAddModel(expect_count,"gaussSum_X",list,coefList, aux.mixedV);
    CreateRooConvPdf(expect_count,"mcSignal_X",mass,aux.mcV,aux.mixedV, aux.mcSignalV);
  }
    break; 
  case _ESFModel_bifurGauss:
    CreateRooRealVar(expect_count,"sigmaR_X",nameEndings,aux.smPars2V,0);
    sf.TransferSmearingFactors(aux.gWidthV,&aux.smPars2V,NULL);
    CreateRooBifurGauss(expect_count,"bifurGauss_X",mass,aux.bias,aux.gWidthV,aux.smPars2V,  aux.bifurGaussV);
    CreateRooConvPdf(expect_count,"mcSignal_X",mass,aux.mcV,aux.bifurGaussV, aux.mcSignalV);
    break;
  case _ESFModel_cball: 
    CreateRooRealVar(expect_count,"parN_X",nameEndings,aux.smPars2V,0);
    CreateRooRealVar(expect_count,"alpha_X",nameEndings,aux.smPars3V,0);
    sf.TransferSmearingFactors(aux.gWidthV, &aux.smPars2V, &aux.smPars3V);
    CreateRooCBShape(expect_count,"cbShape_X",mass,aux.bias,aux.gWidthV,aux.smPars3V,aux.smPars2V, aux.cballV);
    CreateRooConvPdf(expect_count,"mcSignal_X",mass,aux.mcV,aux.cballV, aux.mcSignalV);
    break;
  case _ESFModel_cballShift: {
    CreateRooRealVar(expect_count,"parN_X",nameEndings,aux.smPars2V,0);
    CreateRooRealVar(expect_count,"alpha_X",nameEndings,aux.smPars3V,0);
    sf.TransferSmearingFactors(aux.gWidthV, &aux.smPars2V, &aux.smPars3V);
    CreateRooCBShape(expect_count,"cbShape_X",mass,aux.bias,aux.gWidthV,aux.smPars3V,aux.smPars2V, aux.cballShiftedV);
    char buf[30];
    for (int i=0; i<expect_count; ++i) {
      RooMoment *mean=aux.cballShiftedV[i]->mean(mass);
      aux.calcMeanV.push_back(mean);
      sprintf(buf,"shiftFromMean%s",nameEndings[i].Data());
      RooFormulaVar *var=new RooFormulaVar(buf,"-@0",RooArgList(*mean));
      aux.shiftFromMeanV.push_back(var);
    }
    CreateRooCBShape(expect_count,"cbShape_X",mass,aux.shiftFromMeanV,aux.gWidthV,aux.smPars3V,aux.smPars2V, aux.cballV);
    CreateRooConvPdf(expect_count,"mcSignal_X",mass,aux.mcV,aux.cballV, aux.mcSignalV);
  }
    break;
  case _ESFModel_breitWigner:
    sf.TransferSmearingFactors(aux.gWidthV,NULL,NULL);
    CreateRooBreitWigner(expect_count,"bw_X",mass,aux.bias,aux.gWidthV, aux.bwShapeV);
    CreateRooConvPdf(expect_count,"mcSignal_X",mass,aux.mcV,aux.bwShapeV, aux.mcSignalV);
    break;
  case _ESFModel_voigtian:
    CreateRooRealVar(expect_count,"sigma_X",nameEndings,aux.smPars2V,0);
    sf.TransferSmearingFactors(aux.gWidthV,&aux.smPars2V,NULL);
    CreateRooVoigtian(expect_count,"voigtian_X",mass,aux.bias,aux.gWidthV,aux.smPars2V, aux.voigtianV);
    CreateRooConvPdf(expect_count,"mcSignal_X",mass,aux.mcV,aux.voigtianV, aux.mcSignalV);
    break;
  default:
    std::cout << "not prepared sf.FitModel=" << ElectronEnergyScaleAdv_t::FitModelName(sf.FitModelKind()) << " in " << fncname << "\n";
    return 0;
  }

  //aux.mcSignalV[0]->plotOn(frameX);
  //PrintValues("gWidthV",aux.gWidthV);
  //mass.setVal(80);
  //PrintValues("gaussV",aux.gaussV);
  //PrintValues("mcV",aux.mcV);

  //std::vector<RooExtendPdf*> mcESigV;
  CreateRooExtendPdf(expect_count,"mcSig_X",aux.mcSignalV,aux.MCOverData,aux.mcESigV);

  std::vector<RooArgList*> mcEList;
  CreateRooArgList(aux.mcESigV,mcEList);

  //std::vector<RooAddPdf*> mcModelV;
  CreateRooAddPdf(expect_count,"mcModel",nameEndings,mcEList,mcModelV);

  /*
  RooCategory sample("sample","sample");
  std::vector<TString> categories;
  CreateRooCategories(expect_count,"set",nameEndings,sample,categories);

  RooDataSet *combData=CreateCombinedData("combData","combined data", *expDataSet, mass,sample,categories);
  RooSimultaneous *simPdf=new RooSimultaneous("simPdf","simultaneous pdf",sample);
  for (unsigned int i=0; i<expect_count; i++) {
    simPdf->addPdf(*mcModelV[i],categories[i]);
  }

  int plot_idx=0;
  TString name1="set " + nameEndings[plot_idx];
  RooPlot* frame = x.frame(RooFit::Title(name1));
  frames.push_back(frame);
  */

  /*
  name1="sample==sample::" + categories[plot_idx];
  combData->plotOn(frame,RooFit::Cut(name1),RooFit::LineColor(kRed),RooFit::MarkerColor(kRed));
  simPdf->plotOn(frame,RooFit::Slice(sample,categories[plot_idx]),RooFit::ProjWData(sample,*combData),RooFit::LineColor(kBlue));
  TCanvas *c4=new TCanvas("c4","c4",600,600);
  gPad->SetLeftMargin(0.15) ; 
  frames[0]->GetYaxis()->SetTitleOffset(1.6) ; 
  frames[0]->Draw() ;
  c4->Update();

  const char *outfile_name="scanscale.png";
  c4->SaveAs(outfile_name);
  */
  std::cout << "mcModelVcreated (" << mcModelV.size() << " entries)" << std::endl;
  std::cout << "returning 1" << std::endl;
  return 1;
}

// ======================================================================

/*
int DetermineESFWorkCase(const char *filename) {
  std::ifstream fin;
  fin.open(filename);
  if (!fin) {
    std::cout << "DetermineESFWorkCase: failed to open a file <" << filename << ">\n";
    return 0;
  }
  std::string s;
  size_t pos;
  int res=int(_ESFWC_6bin);
  while (!fin.eof() && getline(fin,s)) {
    pos=s.find("g_esfWorkCase=");
    if (PosOk(pos)) {
      res=atoi(s.c_str()+pos+strlen("g_esfWorkCase="));
      std::cout << "got g_esfWorkCase=" << res << "\n";
    }
  }
  fin.close();
  return res;
}
*/

// ======================================================================
