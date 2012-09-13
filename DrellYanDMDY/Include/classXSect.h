#ifndef classXSect_H
#define classXSect_H


#include "../Include/DYTools.hh"
#include "../Include/TriggerSelection.hh"


// -------------------------------------------------------------

const int nStandardColors=4;
const int standardColors[nStandardColors] = { 814, 894, 906, 855 };

const int extraFlag_fineGrid=1;
const int extraFlag_2MCfiles=2;
const int extraFlag_2MCfilesFineGrid=3;
const int extraFlag_2MCfiles_debug=4;



// -------------------------------------------------------------

inline
int getStandardColor(int i) { return standardColors[i%nStandardColors]; }

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
    std::cout << "readTh: for 1D case v[DYTools::nMassBins] is expected\n";
    ok=false;
  }
  else if ((DYTools::study2D==1) && (v.GetNoElements()!=DYTools::nYBins[iMassBin])) {
    std::cout << "readTh: for 2D case v[DYTools::nYBins[iMassBin=" << iMassBin << "] is expected\n";
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
  HERE("leaving readTh");
  return histo;
}

// -------------------------------------------------------------

inline
TH1F* readThCT10(TVectorD &v, TVectorD &vErr, DYTools::TCrossSectionKind_t theKind, int iMassBin) {

  TString xSecThResultFileName="xSectTheory_CTEQ10W.root";
  std::cout << "Load CT10 theory predictions\n";
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

  TString name=Form("xsecThCT10_%s_%s",DYTools::analysisTag.Data(),CrossSectionKindName(theKind).Data());
  HERE("a");
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
  
  HERE("b");
  delete xSecTh;
  delete xSecThErr;

  return histo;
}

// -------------------------------------------------------------


inline
TH1F* readData(TVectorD &v, TVectorD &vErr1, TVectorD &vErr2, const TriggerSelection &triggers, DYTools::TCrossSectionKind_t theKind, int iMassBin, TString fname){

  //printf("Load data yields\n"); fflush(stdout);
  std::cout << "Load data yields for iMassBin=" << iMassBin << std::endl;
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
  std::cout << "xSecResultFileName= " << xSecResultFileName << "\n";

  TFile fileXsecResult   (xSecResultFileName);
  TString xSecName,xSecErrName,xSecSystErrName;
  switch(theKind) {
  case DYTools::_cs_preFsrDetNorm:
  case DYTools::_cs_preFsrNorm:
    xSecName="normXSec"; xSecErrName="normXSecErr"; xSecSystErrName="normXSecErrSyst";
    break;
  case DYTools::_cs_preFsrDet:
  case DYTools::_cs_preFsr:
    xSecName="XSec"; xSecErrName="XSecErr"; xSecSystErrName="XSecSystErr";
    break;
  default:
    std::cout << "readData is not ready for theKind=<" << CrossSectionKindName(theKind) << ">\n";
    assert(0);
  }

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
    histo= extractMassDependence(name,name, xSecM, xSecErr1M, 0, 0,0);
    HERE("loop");
    for(int i=0; i<DYTools::nMassBins; i++){
      HERE("i=%d",i);
      v[i] = xSecM[i][0];
      vErr1[i] = xSecErr1M[i][0];
      vErr2[i] = xSecErr2M[i][0];
    }
  }
  else if (ok && (DYTools::study2D==1)) {
    histo= extractRapidityDependence(name,name, xSecM, xSecErr1M, iMassBin, 0);
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
TH1F* readTh2011(DYTools::TCrossSectionKind_t theKind, const char *histName="th2011",const char *fname="../root_files/xSecTh2011_1D.root") {
  TString objName,objErrName;
  switch(theKind) {
  case DYTools::_cs_preFsrNorm: 
    objName="xSecThNorm"; objErrName="xSecThNormErr";
    break;
  case DYTools::_cs_preFsr:
    objName="xSecTh"; objErrName="xSecThErr";
    break;
  default:
    std::cout << "readTh2011 is not ready for theKind=<" << CrossSectionKindName(theKind) << ">\n";
    assert(0);
  }
  TFile f(fname);
  if (!f.IsOpen()) {
    std::cout << "readTh2011: failed to open a file <" << fname << ">\n";
    assert(0);
  }
  TVectorD massBins= *(TVectorD*)f.FindObjectAny("massBinsTh");
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
  int nUnfoldingBins=DYTools::nUnfoldingBinsMax;
  TVectorD v(nUnfoldingBins);
  TVectorD vErr(nUnfoldingBins);
  return readTh(v,vErr,triggers,theKind,iMassBin,useFEWZ,extraFlag);
}

// -------------------------------------------------------------

inline
TH1F* readThCT10(DYTools::TCrossSectionKind_t theKind, int iMassBin) {
  int nUnfoldingBins=DYTools::nUnfoldingBinsMax;
  TVectorD v(nUnfoldingBins);
  TVectorD vErr(nUnfoldingBins);
  return readThCT10(v,vErr,theKind,iMassBin);
}

// -------------------------------------------------------------

inline
TH1F* readData(const TriggerSelection &triggers, DYTools::TCrossSectionKind_t theKind, int iMassBin, TString fname){

  int nUnfoldingBins=DYTools::nUnfoldingBinsMax;
  TVectorD v(nUnfoldingBins);
  TVectorD vErr1(nUnfoldingBins);
  TVectorD vErr2(nUnfoldingBins);
  return readData(v,vErr1,vErr2,triggers,theKind,iMassBin,fname);
}

// -------------------------------------------------------------

inline
void removeError(TH1F* h) {
  std::cout << "nulifying error for " << h->GetName() << "\n";
  for (int ibin=0; ibin<=h->GetNbinsX(); ++ibin) {
    h->SetBinError(ibin,0);
  }
}

// -------------------------------------------------------------



#endif
