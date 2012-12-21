#ifndef classXSect_H
#define classXSect_H


#include "../Include/DYTools.hh"
#include "../Include/TriggerSelection.hh"


// -------------------------------------------------------------

typedef enum { _th2011_default, _th2011_mb2011, _th2011_spec2011, _th2011_nnlo } TTheory2011Set_t;
typedef enum { _yields_none=0, _yields_gen, _yields_reco, 
	       _yields_data_observed, _yields_data_signal, 
	       _yields_mc_signal, _yields_mc_signal_scaled,
	       _yields_background_tot, _yields_background_true2e, _yields_background_fake2e, 
	       _yields_efficiency, _yields_efficiencyBB, _yields_efficiencyBE, _yields_efficiencyEE,
	       _yields_eff_NPass, _yields_eff_NTotal, 
	       _yields_eff_NTotalBB, _yields_eff_NTotalBE, _yields_eff_NTotalEE,
	       _yields_eff_NPassBB, _yields_eff_NPassBE, _yields_eff_NPassEE,
	       _yields_eff_scaleFactors, _yields_eff_scaleFactorsFI,
	       _yields_effScale, // eff*ScaleFactors
	       _yields_acc_NTotal, _yields_acc_NPass,
	       _yields_XSec, _yields_XSecNorm,
	       _yields_thXSec, _yields_thXSecDET,
	       _yields_thXSecNorm, _yields_thXSecDETNorm,
	       _yields_fsrCorrFactorsDET_Gen, _yields_fsrCorrFactorsDET_Reco,
	       _yields_fsrCorrFactors_Gen, _yields_fsrCorrFactors_Reco,
	       _yields_Last
} TYieldsKind_t;

const int nStandardColors=5;
//const int standardColors[nStandardColors] = { 814, 894, 906, 855, kMagenta+2 };
const int standardColors[nStandardColors] = { 814, kOrange+8, 855, kMagenta+2, kPink-2 };

const int extraFlag_none=0;
const int extraFlag_fineGrid=1;
const int extraFlag_gridSummer2011=2;
const int extraFlag_gridSummer2011spec=3;
const int extraFlag_2MCfiles=4;
const int extraFlag_2MCfilesFineGrid=5;
const int extraFlag_2MCfiles_debug=6;

const int extraFlagData_oldStyle=7;


// -------------------------------------------------------------

inline
int getStandardColor(int i) { return standardColors[i%nStandardColors]; }

// -------------------------------------------------------------

inline
int csInDET(DYTools::TCrossSectionKind_t kind) {
  int yes=-1;
  switch(kind) {
  case DYTools::_cs_None:
    yes=-1;
    break;
  case DYTools::_cs_preFsr: 
  case DYTools::_cs_preFsrNorm:
  case DYTools::_cs_postFsr:
  case DYTools::_cs_postFsrNorm:
    yes=0;
    break;
  case DYTools::_cs_preFsrDet: 
  case DYTools::_cs_preFsrDetNorm:
  case DYTools::_cs_postFsrDet:
  case DYTools::_cs_postFsrDetNorm:
    yes=1;
    break;
  default:
    std::cout << "csInDET cannot handle this kind=<" << CrossSectionKindName(kind) << ">\n";
    assert(0);
  }
  return yes;
}
// -------------------------------------------------------------

inline
int csPreFsr(DYTools::TCrossSectionKind_t kind) {
  int yes=-1;
  switch(kind) {
  case DYTools::_cs_None:
    yes=-1;
    break;
  case DYTools::_cs_preFsr: 
  case DYTools::_cs_preFsrNorm:
  case DYTools::_cs_preFsrDet: 
  case DYTools::_cs_preFsrDetNorm:
    yes=1;
    break;
  case DYTools::_cs_postFsr:
  case DYTools::_cs_postFsrNorm:
  case DYTools::_cs_postFsrDet:
  case DYTools::_cs_postFsrDetNorm:
    yes=0;
    break;
  default:
    std::cout << "csPreFsr cannot handle this kind=<" << CrossSectionKindName(kind) << ">\n";
    assert(0);
  }
  return yes;
}

// -------------------------------------------------------------

inline
TH1F* readTh(TVectorD &v, TVectorD &vErr, const TriggerSelection &triggers, DYTools::TCrossSectionKind_t theKind, int iMassBin, int useFEWZ, int extraFlag){

  v=0; vErr=0;
  assert(triggers.isDefined()); // get rid of compiler complaint

  //printf("Load data yields\n"); fflush(stdout);
  //TFile fileXsecTh   (TString("../root_files/xSecTh_results.root"));

  TString xSecThResultFileName(TString("../root_files/xSecThExt_") +
			       DYTools::analysisTag + TString("_tmp") +
			       //triggers.triggerConditionsName() +
			       TString(".root"));
  std::cout << "readTh: Load theory predictions for <" << CrossSectionKindName(theKind) << ">, extraFlag=" << extraFlag << "\n";
  if (theKind==DYTools::_cs_None) { std::cout << ".. requested _cs_None\n"; return NULL; }
  TString extra;
  int fineGrid=0;
  switch(extraFlag) {
  case extraFlag_fineGrid: extra= "fineGrid_"; fineGrid=1; break;
  case extraFlag_gridSummer2011: extra="summer2011Grid_"; fineGrid=1; break;
  case extraFlag_gridSummer2011spec: extra="summer2011specGrid_"; fineGrid=1; break;
  case extraFlag_2MCfiles: extra= "2MCfiles_"; break;
  case extraFlag_2MCfilesFineGrid: extra="2MCfiles_fineGrid_"; fineGrid=1; break;
  case extraFlag_2MCfiles_debug: extra= "2MCfiles_debug_"; break;
  default:
    // nothing
    ;
  }
  if (useFEWZ) {
    xSecThResultFileName=TString("../root_files/xSecThExt_") + extra +
      DYTools::analysisTag + TString("_tmp.root");
  }
  else {
    xSecThResultFileName=TString("../root_files/xSecThExt_noFEWZ_") + extra +
      DYTools::analysisTag + TString("_tmp.root");
  }

  //if (extraFlag!=1) xSecThResultFileName="xSecTh_2D_tmp_noWeight.root";
  //if (extraFlag==2) xSecThResultFileName="../root_files/xSecTh_noFEWZ_2D_tmp.root";
  std::cout << "xSecThResultFileName=" << xSecThResultFileName << std::endl;

  TFile fileXsecTh   (xSecThResultFileName);
  if (!fileXsecTh.IsOpen()) {
    std::cout << "readTh: failed to open the file\n";
    return NULL;
  }

  TMatrixD *xSecTh=NULL, *xSecThErr=NULL;
  TH1F *histo=NULL;
  switch(theKind) {
  case DYTools::_cs_preFsr:
    xSecTh= (TMatrixD*)fileXsecTh.FindObjectAny("nGenEvents");
    xSecThErr= (TMatrixD*)fileXsecTh.FindObjectAny("nGenEventsErr");
    break;
  case DYTools::_cs_preFsrNorm:
    xSecTh= (TMatrixD*)fileXsecTh.FindObjectAny("nGenEventsNorm");
    xSecThErr= (TMatrixD*)fileXsecTh.FindObjectAny("nGenEventsNormErr");
    break;
  case DYTools::_cs_preFsrDet:
    xSecTh= (TMatrixD*)fileXsecTh.FindObjectAny("nGenEventsDET");
    xSecThErr= (TMatrixD*)fileXsecTh.FindObjectAny("nGenEventsDETErr");
    break;
  case DYTools::_cs_preFsrDetNorm:
    xSecTh= (TMatrixD*)fileXsecTh.FindObjectAny("nGenEventsDETNorm");
    xSecThErr= (TMatrixD*)fileXsecTh.FindObjectAny("nGenEventsDETNormErr");
    break;
  case DYTools::_cs_postFsrDet:
    xSecTh= (TMatrixD*)fileXsecTh.FindObjectAny("nGenEventsDETrecoPostIdx");
    xSecThErr= (TMatrixD*)fileXsecTh.FindObjectAny("nGenEventsDETrecoPostIdxErr");
    break;
  case DYTools::_cs_postFsrDetNorm:
    xSecTh= (TMatrixD*)fileXsecTh.FindObjectAny("nGenEventsDETrecoPostIdxNorm");
    xSecThErr= (TMatrixD*)fileXsecTh.FindObjectAny("nGenEventsDETrecoPostIdxNormErr");
    break;
  default:
    xSecTh=NULL; xSecThErr=NULL;
  }

  if (!xSecTh || !xSecThErr) {
    std::cout << "readTh: failed to load the arrays\n";
    return NULL;
  }

  //std::cout << "xSecTh="; xSecTh->Print(); std::cout << "\n";
  //std::cout << "xSecThErr="; xSecThErr->Print(); std::cout << std::endl;
  std::cout << "nullifying xSecThErr\n";
  *xSecThErr=0;
    
  // Prepare output yields and errors
  bool ok=true;
  if ((DYTools::study2D==0) && (v.GetNoElements()!=DYTools::nMassBins)) {
    std::cout << "readTh: for 1D case v[DYTools::nMassBins]=" << DYTools::nMassBins << " is expected\n";
    ok=false;
  }
  else if ((DYTools::study2D==1) && (v.GetNoElements()!=DYTools::nYBins[iMassBin])) {
    std::cout << "readTh: for 2D case v[DYTools::nYBins[iMassBin=" << iMassBin << "]=" << DYTools::nYBins[iMassBin] << " is expected, not " << v.GetNoElements() << "\n";
    ok=false;
  }

  TString name=Form("xsecTh_%s_%s",DYTools::analysisTag.Data(),CrossSectionKindName(theKind).Data());
  if (ok && (DYTools::study2D==0)) {
    const int perMassBinWidth=1;
    const int perRapidityBinWidth=0;
    if (fineGrid) {
      TVectorD *massBinEdges= (TVectorD*)fileXsecTh.FindObjectAny("massBinEdges");
      TVectorD *rapidityBinCounts= (TVectorD*)fileXsecTh.FindObjectAny("rapidityBinCount");
      if (!massBinEdges || !rapidityBinCounts) {
	std::cout << "readTh: for extraFlag it is expected that the file <" 
		  << xSecThResultFileName << "> will contain vectors " 
		  << "'massBinEdges' and 'rapidityBinCount'\n";
	assert(0);
      }
      histo= extractMassDependenceSpec(name,name, *xSecTh,*xSecThErr, 0,
				       *massBinEdges, *rapidityBinCounts,
				       perMassBinWidth,perRapidityBinWidth);
      delete massBinEdges;
      delete rapidityBinCounts;
      v=0; vErr=0;
    }
    else {
      histo= extractMassDependence(name,name, *xSecTh, *xSecThErr, 0, 
				   perMassBinWidth,perRapidityBinWidth);
      for(int i=0; i<DYTools::nMassBins; i++){
	v[i] = (*xSecTh)[i][0];
	vErr[i] = (*xSecThErr)[i][0];
      }
    }
  }
  else if (ok && (DYTools::study2D==1)) {
    histo= extractRapidityDependence(name,name, *xSecTh, *xSecThErr, iMassBin, 0);
    //printHisto(std::cout,histo);
    for(int iY=0; iY<DYTools::nYBins[iMassBin]; iY++){
      v[iY] = (*xSecTh)[iMassBin][iY];
      vErr[iY] = (*xSecThErr)[iMassBin][iY];
    }
  }
  
  delete xSecTh;
  delete xSecThErr;

  fileXsecTh.Close();
  if (extraFlag) histo->GetXaxis()->SetRangeUser(15.,1500.);
  return histo;
}

// -------------------------------------------------------------

inline
TH1F* readTheory2D(TVectorD &v, TVectorD &vErr, DYTools::TCrossSectionKind_t theKind, int iMassBin, TString xSecThResultFileName, const char *histName) {
  std::cout << "Load 2D theory predictions\n";
  std::cout << "xSecThResultFileName=" << xSecThResultFileName << std::endl;
  if (theKind==DYTools::_cs_None) { std::cout << ".. requested _cs_None\n"; return NULL; }

  TFile fileXsecTh   (xSecThResultFileName);
  TMatrixD *xSecTh=NULL, *xSecThErr=NULL;
  TH1F *histo=NULL;
  switch(theKind) {
  case DYTools::_cs_preFsr:
    //xSecTh= (TMatrixD*)fileXsecTh.FindObjectAny("nGenEvents");
    //xSecThErr= (TMatrixD*)fileXsecTh.FindObjectAny("nGenEventsErr");
    break;
  case DYTools::_cs_preFsrNorm:
    //xSecTh= (TMatrixD*)fileXsecTh.FindObjectAny("nGenEventsNorm");
    //xSecThErr= (TMatrixD*)fileXsecTh.FindObjectAny("nGenEventsNormErr");
    break;
  case DYTools::_cs_preFsrDet:
    xSecTh= (TMatrixD*)fileXsecTh.FindObjectAny("xSectThDET");
    xSecThErr= (TMatrixD*)fileXsecTh.FindObjectAny("xSectThDETErr");
    break;
  case DYTools::_cs_preFsrDetErr:
    xSecTh= (TMatrixD*)fileXsecTh.FindObjectAny("xSectThDETErr");
    xSecThErr= (TMatrixD*)fileXsecTh.FindObjectAny("xSectThDETSystErrPos");
    break;
  case DYTools::_cs_preFsrDetSystErr:
    xSecTh= (TMatrixD*)fileXsecTh.FindObjectAny("xSectThDETSystErrPos");
    xSecThErr= (TMatrixD*)fileXsecTh.FindObjectAny("xSectThDETErr");
    break;
  case DYTools::_cs_preFsrDetNorm:
    xSecTh= (TMatrixD*)fileXsecTh.FindObjectAny("xSectThDETNorm");
    xSecThErr= (TMatrixD*)fileXsecTh.FindObjectAny("xSectThDETNormErr");
    break;
  case DYTools::_cs_preFsrDetNormErr:
    xSecTh= (TMatrixD*)fileXsecTh.FindObjectAny("xSectThDETNormErr");
    xSecThErr= (TMatrixD*)fileXsecTh.FindObjectAny("xSectThDETNormSystErrPos");
    break;
  default:
    xSecTh=NULL; xSecThErr=NULL;
  }
  fileXsecTh.Close();


  if (!xSecTh || !xSecThErr) {
    std::cout << "readTh: failed to load the arrays\n";
    return NULL;
  }

  // Check that the binning is consistent
  // Prepare output yields and errors
  bool ok=true;
  if ((DYTools::study2D==0) && (v.GetNoElements()!=DYTools::nMassBins)) {
    std::cout << "readTh: for 1D case v[DYTools::nMassBins] is expected\n";
    ok=false;
  }
  else if ((DYTools::study2D==1) && (v.GetNoElements()!=DYTools::nYBins[iMassBin])) {
    std::cout << "readTh: for 2D case v[DYTools::nYBins[iMassBin=" << iMassBin << "] is expected\n";
    ok=false;
  }

  TString name=Form("%s_%s_%s",histName,DYTools::analysisTag.Data(),CrossSectionKindName(theKind).Data());
  if (ok && (DYTools::study2D==0)) {
    const int perMassBinWidth=1;
    const int perRapidityBinWidth=0;
    histo= extractMassDependence(name,name, *xSecTh, *xSecThErr, 0, 
				 perMassBinWidth,perRapidityBinWidth);
    for(int i=0; i<DYTools::nMassBins; i++){
      v[i] = (*xSecTh)[i][0];
      vErr[i] = (*xSecThErr)[i][0];
    }
  }
  else if (ok && (DYTools::study2D==1)) {
    histo= extractRapidityDependence(name,name, *xSecTh, *xSecThErr, iMassBin, 0);
    //printHisto(std::cout,histo);
    for(int iY=0; iY<DYTools::nYBins[iMassBin]; iY++){
      v[iY] = (*xSecTh)[iMassBin][iY];
      vErr[iY] = (*xSecThErr)[iMassBin][iY];
    }
  }
  
  delete xSecTh;
  delete xSecThErr;

  return histo;
}

// -------------------------------------------------------------


inline
TH1F* readData(TVectorD &v, TVectorD &vErr1, TVectorD &vErr2, const TriggerSelection &triggers, DYTools::TCrossSectionKind_t theKind, int iMassBin, TString fname, int extraFlag=0){

  //printf("Load data yields\n"); fflush(stdout);
  std::cout << "Load data yields for iMassBin=" << iMassBin << ", extraFlag=" << extraFlag << std::endl;
  if (theKind==DYTools::_cs_None) { std::cout << ".. requested _cs_None\n"; return NULL; }
  if ((iMassBin>=0) && (DYTools::study2D==0)) {
    std::cout << "readData. iMassBin should be -1, if study2D=0\n";
    assert(0);
  }

  TString xSecResultFileName(TString("../root_files/xSec_results_") + 
			     DYTools::analysisTag + TString("_") +
		   triggers.triggerConditionsName() + TString(".root"));

  xSecResultFileName="dir-LatestResults/xSec_results_2D_Full2011_hltEffOld.root";
  //xSecResultFileName="dir-20120629-results/xSec_results_2D_Full2011_hltEffOld.root";
  if (fname.Length()) xSecResultFileName=fname;
  std::cout << "xSecResultFileName= " << xSecResultFileName << std::endl;

  TFile fileXsecResult   (xSecResultFileName);
  TString xSecName,xSecErrName,xSecSystErrName;
  switch(theKind) {
  case DYTools::_cs_preFsrDetNorm:
  case DYTools::_cs_preFsrNorm:
  case DYTools::_cs_postFsrDetNorm:
  case DYTools::_cs_postFsrNorm:
    xSecName="normXSec"; xSecErrName="normXSecErr"; xSecSystErrName="normXSecErrSyst";
    if (extraFlag==extraFlagData_oldStyle) {
      xSecName="normXSecByBin"; xSecErrName="normXSecErrByBin";
      xSecSystErrName="normXSecErrByBinSyst";
    }
    break;
  case DYTools::_cs_preFsrDet:
  case DYTools::_cs_preFsr:
  case DYTools::_cs_postFsrDet:
  case DYTools::_cs_postFsr:
    xSecName="XSec"; xSecErrName="XSecErr"; xSecSystErrName="XSecSystErr";
    break;
  default:
    std::cout << "readData is not ready for theKind=<" << CrossSectionKindName(theKind) << ">\n";
    assert(0);
  }

  std::cout << "loading " << xSecName << ", " << xSecErrName << ", " << xSecSystErrName << "\n";

  TMatrixD xSecM          = *(TMatrixD *)fileXsecResult.FindObjectAny(xSecName);
  TMatrixD xSecErr1M      = *(TMatrixD *)fileXsecResult.FindObjectAny(xSecErrName);
  TMatrixD xSecErr2M      = *(TMatrixD *)fileXsecResult.FindObjectAny(xSecSystErrName);

  // Check that the binning is consistent
  bool checkResult = true;
  if ((iMassBin>=0) && ( v.GetNoElements() != DYTools::nYBins[iMassBin] )) {
    checkResult = false;
  }
  if( !checkResult ){
    printf("ERROR: Data inconsistent binning in the inputs\n");
    assert(0);
  }else
    printf("readData: Binning in the inputs is consistent\n");

  //printYields(Form("from file <%s>",xSecResultFileName.Data()), xSecM,xSecErr1M,xSecErr2M);

  TH1F* histo=NULL;
  bool ok=true;
  TString name=Form("xsec_%s_%s_massBin%d",DYTools::analysisTag.Data(),CrossSectionKindName(theKind).Data(),iMassBin);
  if (ok && (DYTools::study2D==0)) {
    int perMassBinWidth=(extraFlag == extraFlagData_oldStyle) ? 0:1;
    if (DYTools::study2D) perMassBinWidth=0;
    int perRapidityBinWidth=0;
    std::cout << "calling extracMassDependence with perMassBinWidth=" << perMassBinWidth << "\n";
    histo= extractMassDependence(name,name, xSecM, xSecErr1M, 0, perMassBinWidth,perRapidityBinWidth);
    for(int i=0; i<DYTools::nMassBins; i++){
      v[i] = xSecM[i][0];
      vErr1[i] = xSecErr1M[i][0];
      vErr2[i] = xSecErr2M[i][0];
    }
  }
  else if (ok && (DYTools::study2D==1)) {
    int perMassBinWidth=0;
    histo= extractRapidityDependence(name,name, xSecM, xSecErr1M, iMassBin, perMassBinWidth);
    for(int iY=0; iY<DYTools::nYBins[iMassBin]; iY++){
      v[iY] = xSecM[iMassBin][iY];
      vErr1[iY] = xSecErr1M[iMassBin][iY];
      vErr2[iY] = xSecErr2M[iMassBin][iY];
    }
    //printHisto(std::cout,histo);
  }

  fileXsecResult.Close();
  return histo;
}

// -------------------------------------------------------------

inline
TH1F* readTh1D_MSTW2008(DYTools::TCrossSectionKind_t theKind, const char *histName="th2011",const char *fname="default", TTheory2011Set_t rebin=_th2011_default ) {
  const char *fnameDefault="../root_files/theory/xSectTheory1D_MSTW2008.root";
  TString objName,objErrName;
  switch(theKind) {
  case DYTools::_cs_None: 
    if (theKind==DYTools::_cs_None) { std::cout << ".. requested _cs_None\n"; return NULL; }
    break;
  case DYTools::_cs_preFsrNorm: 
    objName="xSecThNorm"; objErrName="xSecThNormErr";
    switch (rebin) {
    case _th2011_default: break;
    case _th2011_mb2011: objName="xSecThNorm_mb2011"; objErrName="xSecThNormErr_mb2011"; break;
    case _th2011_spec2011: objName="xSecThNorm_spec2011"; objErrName="xSecThNormErr_spec2011"; break;
    case _th2011_nnlo: objName="xSecThNorm_NNLO"; objErrName="xSecThNormErr_NNLO"; break;
    default:
      std::cout << "not ready (cs_preFsrNorm)\n";
      assert(0);
    }
    break;
  case DYTools::_cs_preFsr:
    objName="xSecTh"; objErrName="xSecThErr";
    switch (rebin) {
    case _th2011_default: break;
    case _th2011_mb2011: objName="xSecTh_mb2011"; objErrName="xSecThErr_mb2011"; break;
    case _th2011_spec2011: objName="xSecTh_spec2011"; objErrName="xSecThErr_spec2011"; break;
    case _th2011_nnlo: // not available
    default:
      std::cout << "not ready (cs_preFsr)\n";
      assert(0);
    }
    break;
  default:
    std::cout << "readTh1D_MSTW2008 is not ready for theKind=<" << CrossSectionKindName(theKind) << ">\n";
    assert(0);
  }

  TString fileName=fname;
  if (fileName == "default") fileName=fnameDefault;
  std::cout << "readTh1D_MSTW2008: fileName=<" << fileName << ">" << std::endl;

  TFile f(fileName);
  if (!f.IsOpen()) {
    std::cout << "readTh1D_MSTW2008: failed to open a file <" << fname << ">\n";
    assert(0);
  }
  TString massBinStr="massBinsTh";
  switch (rebin) {
  case _th2011_default: break;
  case _th2011_mb2011: massBinStr="massBinsTh_mb2011"; break;
  case _th2011_spec2011: massBinStr="massBinsTh_spec2011"; break;
  case _th2011_nnlo: massBinStr="massBinsTh_NNLO"; break;
  default:
    std::cout << "not ready\n";
    assert(0);
  }

  TVectorD massBins= *(TVectorD*)f.FindObjectAny(massBinStr);
  int numBins=massBins.GetNoElements()-1;
  TVectorD xsec(numBins), xsecErr(numBins);
  xsec.Read(objName);
  xsecErr.Read(objErrName);
  f.Close();

  double *massRange= new double[numBins+1];
  for (int i=0; i<=numBins; ++i) massRange[i]=massBins[i];
  TH1F *hist=new TH1F(histName,histName, numBins,massRange);
  hist->SetDirectory(0);
  delete [] massRange;
  for (int i=0; i<numBins; ++i) {
    hist->SetBinContent(i+1, xsec[i]);
    hist->SetBinError(i+1, xsecErr[i]);
  }
  return hist;
}

// -------------------------------------------------------------

inline
TH1F* readTh(const TriggerSelection &triggers, DYTools::TCrossSectionKind_t theKind, int iMassBin, int useFEWZ, int extraFlag) {
  int nBins=(DYTools::study2D) ? DYTools::nYBins[iMassBin] : DYTools::nMassBins;
  TVectorD v(nBins);
  TVectorD vErr(nBins);
  return readTh(v,vErr,triggers,theKind,iMassBin,useFEWZ,extraFlag);
}

// -------------------------------------------------------------

inline
TH1F* readTh2D_CT10(DYTools::TCrossSectionKind_t theKind, int iMassBin) {
  int nBins=(DYTools::study2D) ? DYTools::nYBins[iMassBin] : DYTools::nMassBins;
  TVectorD v(nBins);
  TVectorD vErr(nBins);
  return readTheory2D(v,vErr,theKind,iMassBin,"../root_files/theory/xSectTheory_CTEQ10W.root","hCT10");
}

// -------------------------------------------------------------

inline
TH1F* readTh2D_MSTW2008(DYTools::TCrossSectionKind_t theKind, int iMassBin) {
  int nBins=(DYTools::study2D) ? DYTools::nYBins[iMassBin] : DYTools::nMassBins;
  TVectorD v(nBins);
  TVectorD vErr(nBins);
  return readTheory2D(v,vErr,theKind,iMassBin,"../root_files/theory/xSectTheory2D_MSTW2008.root","hMSTW2008");
}

// -------------------------------------------------------------

inline
TH1F* readTh2D_HERAPDF15(DYTools::TCrossSectionKind_t theKind, int iMassBin) {
  int nBins=(DYTools::study2D) ? DYTools::nYBins[iMassBin] : DYTools::nMassBins;
  TVectorD v(nBins);
  TVectorD vErr(nBins);
  return readTheory2D(v,vErr,theKind,iMassBin,"../root_files/theory/xSectTheory2D_HERAPDF15.root","hHERAPDF15");
}

// -------------------------------------------------------------

inline
TH1F* readTh2D_ABKM09(DYTools::TCrossSectionKind_t theKind, int iMassBin) {
  int nBins=(DYTools::study2D) ? DYTools::nYBins[iMassBin] : DYTools::nMassBins;
  TVectorD v(nBins);
  TVectorD vErr(nBins);
  return readTheory2D(v,vErr,theKind,iMassBin,"../root_files/theory/xSectTheory2D_ABKM09.root","hABKM09");
}

// -------------------------------------------------------------

inline
TH1F* readData(const TriggerSelection &triggers, DYTools::TCrossSectionKind_t theKind, int iMassBin, TString fname, int extraFlag=0){
  std::cout << "readData <" << fname << ">" << std::endl;
  int nBins=(DYTools::study2D) ? DYTools::nYBins[iMassBin] : DYTools::nMassBins;
  TVectorD v(nBins);
  TVectorD vErr1(nBins);
  TVectorD vErr2(nBins);
  return readData(v,vErr1,vErr2,triggers,theKind,iMassBin,fname,extraFlag);
}

// -------------------------------------------------------------

inline
TH1F* readYields(const TYieldsKind_t kind, int iMassBin, TString fname, int hasError=0){
  TString name, nameErr;
  switch (kind) {
  case _yields_gen: name="yieldsMcPostFsrGen"; hasError=0; break;
  case _yields_reco: name="yieldsMcPostFsrRec"; hasError=0; break;
  case _yields_data_observed: name="observedYields"; nameErr="observedYieldsErr"; break;
  case _yields_data_signal: name="YieldsSignal"; nameErr="YieldsSignalErr"; break;
  case _yields_mc_signal: name="mcYieldsSignal"; nameErr="mcYieldsSignalErr"; break;
  case _yields_mc_signal_scaled: name="mcYieldsSignalScaled"; nameErr="mcYieldsSignalScaledErr"; break;
  case _yields_background_tot: name="totalBackground"; nameErr="totalBackgroundErr"; break;
  case _yields_background_true2e: name="true2eBkgr"; nameErr="true2eBkgrErr"; break;
  case _yields_background_fake2e: name="fake2eBkgr"; nameErr="fake2eBkgrErr"; break;
  case _yields_efficiency: name="efficiencyArray"; name="efficiencyErrArray"; break;
  case _yields_efficiencyBB: name="efficiencyBB"; name="efficiencyBBErr"; break;
  case _yields_efficiencyBE: name="efficiencyBE"; name="efficiencyBEErr"; break;
  case _yields_efficiencyEE: name="efficiencyEE"; name="efficiencyEEErr"; break;
  case _yields_eff_NPass: name="effEval_nPass"; nameErr="effEval_nPassErr"; break;
  case _yields_eff_NPassBB: name="effEval_nPassBB"; nameErr="effEval_nPassBBErr"; break;
  case _yields_eff_NPassBE: name="effEval_nPassBE"; nameErr="effEval_nPassBEErr"; break;
  case _yields_eff_NPassEE: name="effEval_nPassEE"; nameErr="effEval_nPassEEErr"; break;
  case _yields_eff_NTotal: name="effEval_nTotal"; nameErr="effEval_nTotalErr"; break;
  case _yields_eff_NTotalBB: name="effEval_nTotalBB"; nameErr="effEval_nTotalBBErr"; break;
  case _yields_eff_NTotalBE: name="effEval_nTotalBE"; nameErr="effEval_nTotalBEErr"; break;
  case _yields_eff_NTotalEE: name="effEval_nTotalEE"; nameErr="effEval_nTotalEEErr"; break;
  case _yields_eff_scaleFactors: name="scaleFactor"; nameErr="scaleFactorErr"; break;
  case _yields_eff_scaleFactorsFI: name="scaleFactorFlatIdxArray"; nameErr="scaleFactorErrFlatIdxArray"; break;
  case _yields_effScale: name="effScaleF"; nameErr="effScaleFErr"; break;
  case _yields_acc_NTotal: name="accEval_nTotal"; nameErr="accEval_nTotalErr"; break;
  case _yields_acc_NPass: name="accEval_nPass"; nameErr="accEval_nPassErr"; break;
  case _yields_XSec: name="XSec"; nameErr="XSecErr"; break;
  case _yields_XSecNorm: name="normXSec"; nameErr="normXSecErr"; break;
  case _yields_thXSec: name="nGenEvents"; nameErr="nGenEventsErr"; break;
  case _yields_thXSecDET: name="nGenEventsDET"; nameErr="nGenEventsDETErr"; break;
  case _yields_thXSecNorm: name="nGenEventsNorm"; nameErr="nGenEventsNormErr"; break;
  case _yields_thXSecDETNorm: name="nGenEventsDETNorm"; nameErr="nGenEventsDETNormErr"; break;
  case _yields_fsrCorrFactorsDET_Gen: name="fsrDETcorrFactorsGen"; hasError=0; break;
  case _yields_fsrCorrFactorsDET_Reco: name="fsrDETcorrFactorsReco"; hasError=0; break;
  case _yields_fsrCorrFactors_Gen: name="fsrCorrFactorsGen"; hasError=0; break;
  case _yields_fsrCorrFactors_Reco: name="fsrCorrFactorsReco"; hasError=0; break;
  default: 
    std::cout << "readYields: is not ready for this kind=" << int(kind) << "\n";
    assert(0);
  }

  std::cout << "from file <" << fname << "> get <" << name << "> and <" << nameErr << ">\n";

  TFile file  (fname);

  TMatrixD *yieldsPtr = NULL;
  TMatrixD *yieldsErrPtr= NULL;
  if (kind==_yields_eff_scaleFactorsFI) {
    TVectorD *yieldTmp= (TVectorD*)file.FindObjectAny(name);
    TVectorD *yieldTmpErr= (!hasError && (nameErr.Length()==0)) ? NULL : (TVectorD*)file.FindObjectAny(nameErr);
    if (yieldTmp) {
      yieldsPtr=new TMatrixD(DYTools::nMassBins,DYTools::nYBinsMax);
      unfolding::deflattenMatrix(*yieldTmp,*yieldsPtr);
    }
    if (yieldTmpErr) {
      yieldsErrPtr=new TMatrixD(DYTools::nMassBins,DYTools::nYBinsMax);
      unfolding::deflattenMatrix(*yieldTmpErr,*yieldsErrPtr);
    }
  }
  else {
    yieldsPtr = (TMatrixD *)file.FindObjectAny(name);
    yieldsErrPtr= (!hasError && (nameErr.Length()==0)) ? NULL : (TMatrixD*)file.FindObjectAny(nameErr);
  }
  if (!yieldsPtr) { 
    std::cout << "failed to get <" << name << "> from <" << fname << "\n";
    assert(0);
  }
  if (hasError && !yieldsErrPtr) {
    std::cout << "failed to get <" << nameErr << "> from <" << fname << "\n";
    assert(0);
  }

  TMatrixD yields=*yieldsPtr;
  TMatrixD yieldsErr=(hasError) ? *yieldsErrPtr : yields;
  if (!hasError) yieldsErr=0;

  TH1F* histo=NULL;
  bool ok=true;
  TString histoName=Form("yields_%s_%s_massBin%d",DYTools::analysisTag.Data(),fname.Data(),iMassBin);
  histoName.ReplaceAll("/","_");
  if (ok && (DYTools::study2D==0)) {
    int perMassBinWidth=0;
    int perRapidityBinWidth=0;
    histo= extractMassDependence(histoName,histoName, yields, yieldsErr, 0, perMassBinWidth,perRapidityBinWidth);
  }
  else if (ok && (DYTools::study2D==1)) {
    int perMassBinWidth=0;
    //int perRapidityBinWidth=0;
    histo= extractRapidityDependence(histoName,histoName, yields, yieldsErr, iMassBin, perMassBinWidth);
  }

  //printHisto(std::cout,histo, 1);
  return histo;

}

// -------------------------------------------------------------

inline
TH1F* readYields(int kindId, int iMassBin, TString fname, int hasError=0){
  TYieldsKind_t kind=TYieldsKind_t(kindId);
  if (int(kind) != kindId) {
    std::cout << "failed to identify yieldsKindId=" << kindId << "\n";
    assert(0);
  }
  return readYields(kind,iMassBin,fname,hasError);
}

// -------------------------------------------------------------

#endif
