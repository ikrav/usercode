#ifndef classXSect_H
#define classXSect_H


#include "../Include/DYTools.hh"
#include "../Include/TriggerSelection.hh"


// -------------------------------------------------------------

typedef enum { _th2011_default, _th2011_mb2011, _th2011_spec2011, _th2011_nnlo } TTheory2011Set_t;

const int nStandardColors=4;
const int standardColors[nStandardColors] = { 814, 894, 906, 855 };

const int extraFlag_none=0;
const int extraFlag_fineGrid=1;
const int extraFlag_gridSummer2011=2;
const int extraFlag_gridSummer2011spec=3;
const int extraFlag_2MCfiles=4;
const int extraFlag_2MCfilesFineGrid=5;
const int extraFlag_2MCfiles_debug=6;


const int extraFlagData_oldStyle=0;


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
  case DYTools::_cs_preFsrDetNorm:
    xSecTh= (TMatrixD*)fileXsecTh.FindObjectAny("xSectThDETNorm");
    xSecThErr= (TMatrixD*)fileXsecTh.FindObjectAny("xSectThDETNormErr");
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
    xSecName="normXSec"; xSecErrName="normXSecErr"; xSecSystErrName="normXSecErrSyst";
    if (extraFlagData_oldStyle) {
      xSecName="normXSecByBin"; xSecErrName="normXSecErrByBin";
      xSecSystErrName="normXSecErrByBinSyst";
    }
    break;
  case DYTools::_cs_preFsrDet:
  case DYTools::_cs_preFsr:
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
    int perMassBinWidth=(extraFlagData_oldStyle) ? 0:1;
    if (DYTools::study2D) perMassBinWidth=0;
    int perRapidityBinWidth=0;
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
  case DYTools::_cs_preFsrNorm: 
    objName="xSecThNorm"; objErrName="xSecThNormErr";
    switch (rebin) {
    case _th2011_default: break;
    case _th2011_mb2011: objName="xSecThNorm_mb2011"; objErrName="xSecThNormErr_mb2011"; break;
    case _th2011_spec2011: objName="xSecThNorm_spec2011"; objErrName="xSecThNormErr_spec2011"; break;
    case _th2011_nnlo: objName="xSecThNorm_NNLO"; objErrName="xSecThNormErr_NNLO"; break;
    default:
      std::cout << "not ready\n";
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
      std::cout << "not ready\n";
      assert(0);
    }
    break;
  default:
    std::cout << "readTh1D_MSTW2008 is not ready for theKind=<" << CrossSectionKindName(theKind) << ">\n";
    assert(0);
  }

  TString fileName=fname;
  if (fileName == "default") fileName=fnameDefault;

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
TH1F* readData(const TriggerSelection &triggers, DYTools::TCrossSectionKind_t theKind, int iMassBin, TString fname){
  int nBins=(DYTools::study2D) ? DYTools::nYBins[iMassBin] : DYTools::nMassBins;
  TVectorD v(nBins);
  TVectorD vErr1(nBins);
  TVectorD vErr2(nBins);
  return readData(v,vErr1,vErr2,triggers,theKind,iMassBin,fname);
}

// -------------------------------------------------------------


#endif
