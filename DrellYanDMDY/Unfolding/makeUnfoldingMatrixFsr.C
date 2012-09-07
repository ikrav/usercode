#if !defined(__CINT__) || defined(__MAKECINT__)

// TODO:
//
// Review how systematics is done
//
// Switch to using EventSelection class for dielectron selection
//

#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TClonesArray.h>           // ROOT array class
#include "TMatrixD.h"
#include "TVectorD.h"
#include <TCanvas.h>                // class for drawing
#include <TH1F.h>                   // 1D histograms
#include <TCanvas.h>                // class for drawing
#include <TH1F.h>                   // 1D histograms
#include <TGraphErrors.h>           // graphs
#include <TProfile.h>               // profile histograms
#include <TBenchmark.h>             // class to track macro running statistics
#include <TLorentzVector.h>         // class for Lorentz vector computations
#include <TVectorD.h>
#include <TStyle.h>
#include <TRandom.h>
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <vector>                   // STL vector class
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings

#include "../Include/CPlot.hh"          // helper class for plots
#include "../Include/MitStyleRemix.hh"  // style settings for drawing
#include "../Include/MyTools.hh"        // miscellaneous helper functions
#include "../Include/DYTools.hh"
#include "../Include/plotFunctions.hh"

// define classes and constants to read in ntuple
#include "../Include/EWKAnaDefs.hh"
#include "../Include/TEventInfo.hh"
#include "../Include/TGenInfo.hh"
#include "../Include/TDielectron.hh"   

// Helper functions for Electron ID selection
#include "../Include/EleIDCuts.hh"

#include "../Include/ElectronEnergyScale.hh" //extra smearing
#include "../Include/UnfoldingTools.hh"

  // Trigger info
#include "../Include/TriggerSelection.hh"

#include "../Include/EventSelector.hh"
#include "../Include/FEWZ.hh"
#include "../Include/InputFileMgr.hh"

//for getting matrix condition number
#include <TDecompLU.h>

#endif


//=== FUNCTION DECLARATIONS ======================================================================================

int nUnfoldingBins = DYTools::getTotalNumberOfBins();

void computeNormalizedBinContent(double subset, double subsetErr,
				 double total, double totalErr,
				 double& ratio, double& ratioErr);
void calculateInvertedMatrixErrors(const TMatrixD &T, 
          const TMatrixD &TErrPos, const TMatrixD &TErrNeg,
				   TMatrixD &TinvErr);

// ---------------------------------------------------------------------

inline bool validFlatIndex(int idx) {
  return ( idx != -1 && idx < nUnfoldingBins ) ? true : false;
}

inline bool validFlatIndices(int idx1, int idx2) {
  return (validFlatIndex(idx1) && validFlatIndex(idx2)) ? true : false;
}


//=== Class DECLARATIONS =================================================================================================

class UnfoldingMatrix_t {
public:
  typedef enum { _cDET_Response, _cFSR, _cFSR_DET } TUnfoldingMatrixType_t;
  // DET_Response: ini --> PostFsrGen, fin --> PostFsrReco
  // FSR, FSR_DET: ini --> PreFsrGen,  fin -->PostFsrGen
public:
  TUnfoldingMatrixType_t kind;
  TString name, iniYieldsName, finYieldsName;
  TMatrixD *yieldsIni; //(DYTools::nMassBins,DYTools::nYBinsMax);
  TMatrixD *yieldsFin; //(DYTools::nMassBins,DYTools::nYBinsMax);

  // Matrices for unfolding
  TMatrixD *DetMigration; //(nUnfoldingBins, nUnfoldingBins);
  TMatrixD *DetMigrationErr; //(nUnfoldingBins, nUnfoldingBins);
  TMatrixD *DetResponse; //(nUnfoldingBins, nUnfoldingBins);
  TMatrixD *DetResponseErrPos; //(nUnfoldingBins, nUnfoldingBins);
  TMatrixD *DetResponseErrNeg; //(nUnfoldingBins, nUnfoldingBins);
  TMatrixD *DetInvertedResponse; //(nUnfoldingBins, nUnfoldingBins);
  TMatrixD *DetInvertedResponseErr; //(nUnfoldingBins, nUnfoldingBins);

  TVectorD *DetResponseArr; //(nUnfoldingBins);
  TVectorD *DetInvertedResponseArr; //(nUnfoldingBins);
  TVectorD *DetInvertedResponseErrArr; //(nUnfoldingBins);
  TVectorD *yieldsIniArr; //(nUnfoldingBins)
  TVectorD *yieldsFinArr; //(nUnfoldingBins);

public:
  static TString kindName(TUnfoldingMatrixType_t aKind) {
    TString s="unknown";
    switch(aKind) {
    case _cDET_Response: s="DetResponse"; break;
    case _cFSR: s="FSR"; break;
    case _cFSR_DET: s="FSR_DET"; break;
    }
    return s;
  }

  TString ourKindName() const { return UnfoldingMatrix_t::kindName(this->kind); }

public:
  UnfoldingMatrix_t(TUnfoldingMatrixType_t set_kind, const TString &set_name) : 
    kind(set_kind),
    name(set_name), iniYieldsName(), finYieldsName(),
    yieldsIni(0), yieldsFin(0),
    DetMigration(0), DetMigrationErr(0),
    DetResponse(0), DetResponseErrPos(0), DetResponseErrNeg(0),
    DetResponseArr(0), DetInvertedResponseArr(0), DetInvertedResponseErrArr(0),
    yieldsIniArr(0), yieldsFinArr(0)
  {
    TMatrixD my(DYTools::nMassBins,DYTools::nYBinsMax);
    TMatrixD unf(nUnfoldingBins, nUnfoldingBins);
    TVectorD arr(nUnfoldingBins);
    my=0; unf=0; arr=0;
    yieldsIni=new TMatrixD(my); yieldsFin=new TMatrixD(my);
    DetMigration=new TMatrixD(unf); DetMigrationErr=new TMatrixD(unf);
    DetResponse=new TMatrixD(unf); 
    DetResponseErrPos=new TMatrixD(unf); DetResponseErrNeg=new TMatrixD(unf);
    DetInvertedResponse=new TMatrixD(unf); DetInvertedResponseErr=new TMatrixD(unf);
    DetResponseArr=new TVectorD(arr);
    DetInvertedResponseArr=new TVectorD(arr); DetInvertedResponseErrArr=new TVectorD(arr);
    yieldsIniArr=new TVectorD(arr); yieldsFinArr=new TVectorD(arr);
    UnfoldingMatrix_t::getYieldNames(set_kind, iniYieldsName, finYieldsName);
  }

  void fillIni(int iMassBinGen, int iYBinGen, double fullWeight) {
    using namespace DYTools;
    if ((iMassBinGen==-1) || (iYBinGen==-1)) {
      return;
    }
    if ((iMassBinGen >= nMassBins) ||
	(iYBinGen >= nYBins[iMassBinGen])) {
      std::cout << "UnfoldingMatrix_t::fillIni(" << iMassBinGen << "," << iYBinGen << "): incorrect indices. Max values (" << nMassBins << "," << (((iMassBinGen>=0) && (iMassBinGen<nMassBins)) ? nYBins[iMassBinGen] : nYBinsMax) << "; matrixName=<" << name << ">\n";
      assert(0);
    }
    (*yieldsIni)(iMassBinGen,iYBinGen) += fullWeight;
  }

  void fillFin(int iMassBinReco, int iYBinReco, double fullWeight) {
    using namespace DYTools;
    if ((iMassBinReco==-1) || (iYBinReco==-1)) {
      return;
    }
    if ((iMassBinReco >= nMassBins) ||
	(iYBinReco >= nYBins[iMassBinReco])) {
      std::cout << "UnfoldingMatrix_t::fillPostFin(" << iMassBinReco << "," << iYBinReco << "): incorrect indices. Max values (" << nMassBins << "," << (((iMassBinReco>=0) && (iMassBinReco<nMassBins)) ? nYBins[iMassBinReco] : nYBinsMax) << "; matrixName=<" << name << ">\n";
      assert(0);
    }
    (*yieldsFin)(iMassBinReco,iYBinReco) += fullWeight;
  }

  void fillMigration(int idx1, int idx2, double weight) {
    if ( !validFlatIndices(idx1,idx2) ) {
      std::cout << "UnfoldingMatrix_t::fillMigration: idx1=" << idx1 << ", idx2=" << idx2 << ". Max allowed values: " << nUnfoldingBins << "; matrixName=<" << name << ">\n";
      return;
    }
    (*DetMigration)(idx1,idx2) += weight;
    (*DetMigrationErr)(idx1,idx2) += weight * weight;
  }

  void finalizeDetMigrationErr() {
    for(int i=0; i < (*DetMigration).GetNrows(); i++)
      for(int j=0; j < (*DetMigration).GetNcols(); j++)
	if( (*DetMigrationErr)(i,j) >=0 )
	  (*DetMigrationErr)(i,j) = sqrt( (*DetMigrationErr)(i,j) );
	else {
	  printf("UnfoldingMatrix_t::finalizeDetMigrationErr Error: negative weights in DetMigrationErr\n");
	  std::cout << "matrixName=<" << name << ">\n";
	  assert(0);
	}
  }

  void computeResponseMatrix() {
    double tCentral, tErr;
    for(int igen = 0; igen < (*DetMigration).GetNrows(); igen++){
      // First find the normalization for the given generator level slice
      double nEventsInGenBin = 0;
      double nEventsInGenBinErr = 0;
      for(int ireco = 0; ireco < (*DetMigration).GetNcols(); ireco++){
	nEventsInGenBin += (*DetMigration)(igen,ireco);
	nEventsInGenBinErr += ((*DetMigrationErr)(igen,ireco)*
			       (*DetMigrationErr)(igen,ireco));
      }
      nEventsInGenBinErr = sqrt(nEventsInGenBinErr);
      
      // Now normalize each element and find errors
      for(int ireco = 0; ireco < (*DetMigration).GetNcols(); ireco++){
	tCentral = 0;
	tErr     = 0;
	computeNormalizedBinContent((*DetMigration)(igen,ireco),
				    (*DetMigrationErr)(igen,ireco),
				    nEventsInGenBin,
				    nEventsInGenBinErr,
				    tCentral, tErr);
	(*DetResponse)      (igen,ireco) = tCentral;
	(*DetResponseErrPos)(igen,ireco) = tErr;
	(*DetResponseErrNeg)(igen,ireco) = tErr;
      }
    }
  }

  void invertResponseMatrix() {
  // Find inverted response matrix
    (*DetInvertedResponse) = (*DetResponse);
    Double_t det;
    (*DetInvertedResponse).Invert(&det);
    calculateInvertedMatrixErrors(*DetResponse, *DetResponseErrPos, *DetResponseErrNeg, *DetInvertedResponseErr);
  }

  void prepareFIArrays() {
    int resFlatten=
      (unfolding::flattenMatrix(*DetResponse, *DetResponseArr) == 1) &&
      (unfolding::flattenMatrix(*DetInvertedResponse, *DetInvertedResponseArr) == 1) &&
      (unfolding::flattenMatrix(*DetInvertedResponseErr, *DetInvertedResponseErrArr) == 1) &&
      (unfolding::flattenMatrix(*yieldsIni, *yieldsIniArr) == 1) &&
      (unfolding::flattenMatrix(*yieldsFin, *yieldsFinArr) == 1);
    if (!resFlatten) {
      std::cout << "Error : failed to flatten the arrays\n";
      assert(0);
    }
  }


  static void getYieldNames(TUnfoldingMatrixType_t theKind, TString &iniName, TString &finName) {
    switch(theKind) {
    case _cDET_Response: 
      iniName="yieldsMcPostFsrGen";
      finName="yieldsMcPostFsrRec";
      break;
    case _cFSR:
      iniName="yieldsMcPreFsrGen";
      finName="yieldsMcPostFsrGen";
      break;
    case _cFSR_DET:
      iniName="yieldsMcPreFsrGenDET";
      finName="yieldsMcPostFsrGenDET";
      break;
    default:
      std::cout << "getYieldNames cannot handle this 'kind'=" << theKind << "\n";
      assert(0);
    }
  }

  void getFileNames(const TString &outputDir,
		    const TString &fileTag,
		    TString &matrixFName, TString &yieldsFName) const {
    matrixFName=outputDir + this->name + TString("_unfolding_constants") + 
      fileTag + TString(".root");
    yieldsFName=outputDir + TString("yields_") + this->name + 
      fileTag + TString(".root");
  }

  void autoSaveToFile(const TString &outputDir, const TString &fileTag) const {
    TString matrixFName,yieldsFName;
    this->getFileNames(outputDir,fileTag, matrixFName,yieldsFName);
    std::cout << "saving to files <" << matrixFName << "> and <" << yieldsFName << ">\n";
    this->saveToFile(matrixFName,yieldsFName);
  }
  
  int autoLoadFromFile(const TString &outputDir, const TString &fileTag) {
    TString matrixFName,yieldsFName;
    this->getFileNames(outputDir,fileTag, matrixFName,yieldsFName);
    std::cout << "loading from files <" << matrixFName << "> and <" << yieldsFName << ">\n";
    return this->loadFromFile(matrixFName,yieldsFName);
  }
  

  void saveToFile(const TString &fileName, const TString &refFileName) const {
    std::cout << "UnfoldingMatrix_t::saveToFile(\n  <" << fileName << ">\n  <" << refFileName << ">) for name=" << this->name << "\n";
    TFile fConst(fileName, "recreate" );
    //name.Write("matrixName");
    (*DetMigration)            .Write("DetMigration");
    (*DetMigrationErr)         .Write("DetMigrationErr");
    (*DetResponse)             .Write("DetResponse");
    (*DetResponseErrPos)       .Write("DetResponseErrPos");
    (*DetResponseErrNeg)       .Write("DetResponseErrNeg");
    (*DetInvertedResponse)     .Write("DetInvertedResponse");
    (*DetInvertedResponseErr)  .Write("DetInvertedResponseErr");
    (*DetResponseArr)          .Write("DetResponseFIArray");
    (*DetInvertedResponseArr)  .Write("DetInvertedResponseFIArray");
    (*DetInvertedResponseErrArr).Write("DetInvertedResponseErrFIArray");
    (*yieldsIni).Write(iniYieldsName);
    (*yieldsFin).Write(finYieldsName);
    (*yieldsIniArr).Write(iniYieldsName + TString("FIArray"));
    (*yieldsFinArr).Write(finYieldsName + TString("FIArray"));
    unfolding::writeBinningArrays(fConst);
    fConst.Close();

    // Store reference MC arrays in a file
    TFile fRef(refFileName, "recreate" );
    (*yieldsIni).Write(iniYieldsName);
    (*yieldsFin).Write(finYieldsName);
    (*yieldsIniArr).Write(iniYieldsName + TString("FIArray"));
    (*yieldsFinArr).Write(finYieldsName + TString("FIArray"));
    unfolding::writeBinningArrays(fRef);
    fRef.Close();
  }

  int loadFromFile(const TString &fileName, const TString &refFileName) {
    std::cout << "UnfoldingMatrix_t::loadFromFile(\n  <" << fileName << ">\n  <" << refFileName << ">) for name=" << this->name << "\n";
    TFile fConst(fileName);
    if (!fConst.IsOpen()) {
      std::cout << "failed to open the file <" << fileName << ">\n";
      return 0;
    }
    if (!unfolding::checkBinningArrays(fConst)) {
      fConst.Close();
      return 0;
    }
    (*DetMigration)            .Read("DetMigration");
    (*DetMigrationErr)         .Read("DetMigrationErr");
    (*DetResponse)             .Read("DetResponse");
    (*DetResponseErrPos)       .Read("DetResponseErrPos");
    (*DetResponseErrNeg)       .Read("DetResponseErrNeg");
    (*DetInvertedResponse)     .Read("DetInvertedResponse");
    (*DetInvertedResponseErr)  .Read("DetInvertedResponseErr");
    (*DetResponseArr)          .Read("DetResponseFIArray");
    (*DetInvertedResponseArr)  .Read("DetInvertedResponseFIArray");
    (*DetInvertedResponseErrArr).Read("DetInvertedResponseErrFIArray");
    (*yieldsIni).Read(iniYieldsName);
    (*yieldsFin).Read(finYieldsName);
    (*yieldsIniArr).Read(iniYieldsName + TString("FIArray"));
    (*yieldsFinArr).Read(finYieldsName + TString("FIArray"));
    fConst.Close();

    // Store reference MC arrays in a file
    TFile fRef(refFileName);
    if (!fRef.IsOpen()) {
      std::cout << "failed to open the file <" << refFileName << ">\n";
      return 0;
    }
    if (!unfolding::checkBinningArrays(fRef)) {
      fRef.Close();
      return 0;
    }
    //(*yieldsIni).Read(iniYieldsName);
    //(*yieldsFin).Read(finYieldsName);
    //(*yieldsIniArr).Read(iniYieldsName + TString("FIArray"));
    //(*yieldsFinArr).Read(finYieldsName + TString("FIArray"));
    fRef.Close();
    return 1;
  }


  void printYields() const {
    std::cout << "Yields of matrix=" << name << " (" 
	      << iniYieldsName << " and " << finYieldsName << ")\n";
    for (int ir=0; ir<(*yieldsIni).GetNrows(); ++ir) {
      std::cout << "ir=" << ir << "\n";
      for (int ic=0; ic<(*yieldsIni).GetNcols(); ++ic) {
	printf(" % 9.6lf  % 9.6lf\n",(*yieldsIni)[ir][ic],(*yieldsFin)[ir][ic]);
      }
      printf("\n");
    }
  }

  void printMigration() const {
    std::cout << "DetMigration of <" << name << ">:\n";
    int printSystErr=0;
    TMatrixD zeroErr=*DetMigrationErr;
    zeroErr=0;
    printCSMatrixValues("DetMigration",*DetMigration,*DetMigrationErr,zeroErr,printSystErr);
  }

  void printResponse() const {
    std::cout << "DetResponse,ErrPos,ErrNeg of <" << name << ">:\n";
    int printSystErr=1;
    printCSMatrixValues("DetResponse",*DetResponse,*DetResponseErrPos,*DetResponseErrNeg,printSystErr);
  }

  void printInvResponse() const {
    std::cout << "DetInvertedResponse of <" << name << ">:\n";
    int printSystErr=0;
    TMatrixD zeroErr=*DetInvertedResponseErr;
    zeroErr=0;
    printCSMatrixValues("DetInvertedResponse",*DetInvertedResponse,*DetInvertedResponseErr,zeroErr,printSystErr);
  }

  void printMatrices() const {
    std::string line(80,'-');
    std::cout << "\n" << line << "\n";
    printMigration();
    printResponse();
    printInvResponse();
    std::cout << line << "\n";
  }

  void printConditionNumber() const {
    //matrix condition number
    TDecompLU lu(*DetResponse);
    double condLU=lu.Condition();
    std::cout << "Matrix=" << name << "\n";
    std::cout << " condition number from TDecompLU condLU= " << condLU << std::endl;
    std::cout << " condition number ||DetResponse||*||DetResponseInv||=" << DetResponse->Norm1()*DetInvertedResponse->Norm1() << std::endl;
    std::cout << " chk ROOT bug: -condLU*||DetResponse||=" << (-condLU*DetResponse->Norm1()) << "\n" << std::endl;
  }

  void prepareHResponse(TH2F **hResponse_out=NULL,
			TH2F **hInvResponse_out=NULL,
			TCanvas **canv=NULL,
			CPlot **plotResponse_out=NULL,
			CPlot **plotInvResponse_out=NULL
			) {
    // Plot response and inverted response matrices
    TString kName=this->ourKindName();
    TH2F *hResponse = new TH2F(TString("hResponse_") + kName,"",
			       nUnfoldingBins, -0.5, nUnfoldingBins-0.5,
			       nUnfoldingBins, -0.5, nUnfoldingBins-0.5);
    TH2F *hInvResponse = new TH2F(TString("hInvResponse") + kName,"",
				  nUnfoldingBins, -0.5, nUnfoldingBins-0.5,
				  nUnfoldingBins, -0.5, nUnfoldingBins-0.5);
    hResponse->SetDirectory(0);
    hInvResponse->SetDirectory(0);
    for(int i=0; i<(*DetResponse).GetNrows(); i++){
      for(int j=0; j<(*DetResponse).GetNcols(); j++){
	hResponse->SetBinContent(i,j, (*DetResponse)(i,j));
	hInvResponse->SetBinContent(i,j, (*DetInvertedResponse)(i,j));
      }
    }
    hResponse->GetYaxis()->SetTitleOffset(1.1);
    hInvResponse->GetYaxis()->SetTitleOffset(1.1);


    TString canvName=TString("canvResponse") + kName;
    TCanvas *e1 = MakeCanvas(canvName,canvName,1200,600);
    e1->Divide(2,1);
    AdjustFor2DplotWithHeight(e1);
    CPlot *plotResponse=
      new CPlot(TString("response") + kName,"",
		"flat index gen",
		"flat index reco");
    plotResponse->AddHist2D(hResponse,"COLZ");
    plotResponse->Draw(e1,false,"png",1);

    CPlot *plotInvResponse=
      new CPlot(TString("invResponse") + kName,"",
		"flat index reco",
		"flat index gen");
    plotInvResponse->AddHist2D(hInvResponse,"COLZ");
    plotInvResponse->Draw(e1,false,"png",2);
    e1->Update();
    SaveCanvas(e1,Form("hResponse_%s_",DYTools::analysisTag.Data()) + kName);
  
    if (hResponse_out) *hResponse_out=hResponse;
    if (hInvResponse_out) *hInvResponse_out=hInvResponse;
    if (canv) *canv=e1;
    if (plotResponse_out) *plotResponse_out=plotResponse;
    if (plotInvResponse_out) *plotInvResponse_out=plotInvResponse;
  }

  // ------------------------------------------------------

  TMatrixD* getReconstructionEffect(const UnfoldingMatrix_t &inexact) const {
    TMatrixD *res=new TMatrixD(*yieldsIni);
    *res=0;
    for(int i=0; i < res->GetNrows(); i++){
      for(int j=0; j < res->GetNcols(); j++){
	double nexact = (*yieldsIni)(i,j);
	double nactual = (*inexact.yieldsIni)(i,j);
	if( nexact != 0 )
	  (*res)(i,j) = (nexact-nactual)/nexact;
      }
    }
    return res;
  }

  // ------------------------------------------------------

  // ------------------------------------------------------

};

//=== MAIN MACRO =================================================================================================

void makeUnfoldingMatrixFsr(const TString input, 
			 const TString triggerSetString="Full2011DatasetTriggers",
			 int systematicsMode = DYTools::NORMAL, 
			 int randomSeed = 1, double reweightFsr = 1.0, 
			 double massLimit = -1.0, int debugMode=0)
//systematicsMode 0 (NORMAL) - no systematic calc
//1 (RESOLUTION_STUDY) - systematic due to smearing, 2 (FSR_STUDY) - systematics due to FSR, reweighting
//check mass spectra with reweightFsr = 0.95; 1.00; 1.05  
//mass value until which do reweighting
{

  // check whether it is a calculation
  if (input.Contains("_DebugRun_")) {
    std::cout << "plotDYUnfoldingMatrix: _DebugRun_ detected. Terminating the script\n";
    return;
  }

  // normal calculation
  gBenchmark->Start("makeUnfoldingMatrix");

  if (systematicsMode==DYTools::NORMAL)
    std::cout<<"Running script in the NORMAL mode"<<std::endl;
  else if (systematicsMode==DYTools::RESOLUTION_STUDY)
    std::cout<<"Running script in the RESOLUTION_STUDY mode"<<std::endl;
  else if (systematicsMode==DYTools::FSR_STUDY)
    std::cout<<"Running script in the FSR_STUDY mode"<<std::endl;
  else if (systematicsMode==DYTools::ESCALE_RESIDUAL)
    std::cout << "Running script in the ESCALE_RESIDUAL mode\n";
  else { 
    std::cout<<"requested mode not recognized"<<std::endl;
    assert(0);
  }

  if (debugMode==1) std::cout << "\n\n\tDEBUG MODE is ON\n\n";
  else if (debugMode==-1) std::cout << "\n\n\tLOADING MODE is ON\n\n";
  
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================
  
//   Bool_t doSave  = false;    // save plots?
  TString format = "png";   // output file format

  vector<TString> fnamev;   // file names   
  vector<TString> labelv;   // legend label
  vector<Int_t>   colorv;   // color in plots
  vector<Int_t>   linev;    // line style
  vector<Double_t> xsecv;
  vector<Double_t> lumiv;
  TString          dirTag;
  TString          escaleTag; // Energy scale calibrations tag

  if (1) {
    MCInputFileMgr_t mcInp; // avoid errors from empty lines
    if (!mcInp.Load(input)) {
      std::cout << "Failed to load mc input file <" << input << ">\n";
      return;
    }
    fnamev=mcInp.fileNames();
    labelv=mcInp.labels();
    colorv=mcInp.colors();
    linev=mcInp.lineStyles();
    xsecv=mcInp.xsecs();
    lumiv=mcInp.lumis();
    dirTag=mcInp.dirTag();
    escaleTag=mcInp.escaleTag();
  }
  else {
  ifstream ifs;
  ifs.open(input.Data());
  assert(ifs.is_open());
  string line;
  Int_t state=0;
  while(getline(ifs,line)) {
    if(line[0]=='#') continue;
    if(state == 0){
      dirTag = TString(line);
      getline(ifs,line);
      stringstream ss3(line); ss3 >> escaleTag;
      state++;
      continue;
    }else{
      string fname;
      Int_t color, linesty;
      stringstream ss(line);
      Double_t xsec;
      ss >> fname >> xsec >> color >> linesty;
      string label = line.substr(line.find('@')+1);
      fnamev.push_back(fname);
      labelv.push_back(label);
      colorv.push_back(color);
      linev.push_back(linesty);
      xsecv.push_back(xsec);
      lumiv.push_back(0);
    }
  }
  ifs.close();
  }
  
  // 
  // Set up energy scale corrections
  //
  ElectronEnergyScale escale(escaleTag);
  escale.print();

  if( !escale.isInitialized()) {
    printf("Failed to match escale calibration. Tag: >>%s<<\n", escaleTag.Data());
    assert(0);
  }

  TriggerConstantSet constantsSet=DetermineTriggerSet(triggerSetString);  
  assert ( constantsSet != TrigSet_UNDEFINED );

  // For MC the trigger does not depend on run number
  const bool isData=kFALSE;
  TriggerSelection requiredTriggers(constantsSet, isData, 0);

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================

  //for the FSR case
  const bool useFewzWeights = false;
  const bool cutZPT100 = true;
  FEWZ_t fewz(useFewzWeights,cutZPT100);
  if (useFewzWeights && !fewz.isInitialized()) {
    std::cout << "failed to prepare FEWZ correction\n";
    throw 2;
  }


  TRandom random;
  // The random seeds are needed only if we are running this script in systematics mode
  int seed = randomSeed;
  random.SetSeed(seed);
  gRandom->SetSeed(seed);
  if(systematicsMode==DYTools::RESOLUTION_STUDY) {
    escale.randomizeSmearingWidth(seed);
  }

  // prepare tools for ESCALE_RESIDUAL
  TMatrixD *shapeWeights=NULL;
  if (systematicsMode==DYTools::ESCALE_RESIDUAL) {
    TString shapeFName=TString("../root_files/yields/") + dirTag + 
      TString("/yields_bg-subtracted") + DYTools::analysisTag + TString(".root");
    std::cout << "Obtaining shape weights from <" << shapeFName << ">\n";
    TFile fshape(shapeFName);
    if (!fshape.IsOpen()) {
      std::cout << "failed to open a file <" << shapeFName << ">\n";
      throw 2;
    }
    shapeWeights = (TMatrixD*)fshape.Get("ZeeMCShapeReweight");
    if (!shapeWeights) {
      std::cout << "failed to find object \"ZeeMCShapeReweight\"\n";
      throw 2;
    }
    dirTag += TString("_escale_residual");
    std::cout << "changing dirTag to <" << dirTag << ">\n";
    (*shapeWeights)(0,0)=1; (*shapeWeights)(1,0)=1; (*shapeWeights)(2,0)=1;
    std::cout << "shapeWeights:\n"; shapeWeights->Print(); // return;
  }

  //  
  // Set up histograms
  //
  vector<TH1F*> hZMassv;//, hZMass2v, hZPtv, hZPt2v, hZyv, hZPhiv;  
  
  char hname[100];
  for(UInt_t ifile = 0; ifile<fnamev.size(); ifile++) {
    sprintf(hname,"hZMass_%i",ifile); hZMassv.push_back(new TH1F(hname,"",500,0,1500)); hZMassv[ifile]->Sumw2();
  }

  TH1F *hMassDiff   = new TH1F("hMassDiff","", 100, -30, 30);
  TH1F *hMassDiffBB = new TH1F("hMassDiffBB","", 100, -30, 30);
  TH1F *hMassDiffEB = new TH1F("hMassDiffEB","", 100, -30, 30);
  TH1F *hMassDiffEE = new TH1F("hMassDiffEE","", 100, -30, 30);

  // These histograms will contain (gen-reco) difference 
  // for each (mass, Y) bin in a flattened format
  TH2F *hMassDiffV = new TH2F("hMassDiffV","",
			      nUnfoldingBins, -0.5, nUnfoldingBins-0.5,
			      100, -50.0, 50.0);
  TH2F *hYDiffV = new TH2F("hYDiffV","",
			   nUnfoldingBins, -0.5, nUnfoldingBins-0.5,
			   100, -5.0, 5.0);

//   TH1F *hMassDiffV[nUnfoldingBins];
//   for(int i=0; i<nUnfoldingBins; i++){
//     sprintf(hname,"hMassDiffV_%d",i);
//     hMassDiffV[i] = new TH1F(hname,"",100,-50,50);
//   }

  UnfoldingMatrix_t detResponse(UnfoldingMatrix_t::_cDET_Response,"detResponse");
  UnfoldingMatrix_t detResponseExact(UnfoldingMatrix_t::_cDET_Response,"detResponseExact");
  UnfoldingMatrix_t fsr(UnfoldingMatrix_t::_cFSR, "fsr");
  UnfoldingMatrix_t fsrDET(UnfoldingMatrix_t::_cFSR_DET,"fsrDET"); // only relevant indices are checked for ini,fin
  UnfoldingMatrix_t fsrDETexact(UnfoldingMatrix_t::_cFSR_DET,"fsrDETexact"); // all indices are checked

  //
  // Access samples and fill histograms
  //  
  TFile *infile=0;
  TTree *eventTree=0;  
    
  // Data structures to store info from TTrees
  mithep::TEventInfo    *info = new mithep::TEventInfo();
  mithep::TGenInfo *gen  = new mithep::TGenInfo();
  TClonesArray *dielectronArr = new TClonesArray("mithep::TDielectron");
  
  // loop over samples  
  if (debugMode!=-1) {
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
  
    // Read input file
    cout << "Processing " << fnamev[ifile] << "..." << endl;
    infile = new TFile(fnamev[ifile]); 
    assert(infile);

    // Get the TTrees
    eventTree = (TTree*)infile->Get("Events"); assert(eventTree);

    // Find weight for events for this file
    // The first file in the list comes with weight 1,
    // all subsequent ones are normalized to xsection and luminosity
    double xsec=xsecv[ifile];
    AdjustXSectionForSkim(infile,xsec,eventTree->GetEntries(),1);
    lumiv[ifile] = eventTree->GetEntries()/xsec;
    double scale = lumiv[0]/lumiv[ifile];
    cout << "       -> sample weight is " << scale << endl;

    // Set branch address to structures that will store the info  
    eventTree->SetBranchAddress("Info",&info);                TBranch *infoBr       = eventTree->GetBranch("Info");
    eventTree->SetBranchAddress("Gen",&gen);                  TBranch *genBr = eventTree->GetBranch("Gen");
    eventTree->SetBranchAddress("Dielectron",&dielectronArr); TBranch *dielectronBr = eventTree->GetBranch("Dielectron");
  
    // loop over events    
    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
      if (debugMode && (ientry>1000000)) break;
      if (ientry%1000000==0) { printProgress("ientry=",ientry,eventTree->GetEntriesFast()); }

      genBr->GetEntry(ientry);
      infoBr->GetEntry(ientry);

      double reweight=1.;
      if (systematicsMode!=DYTools::FSR_STUDY) reweight=1.0;
      else if (((gen->mass)-(gen->vmass))>massLimit) reweight=1.0;
      else reweight=reweightFsr;

      double fewz_weight = 1.0;
      if (useFewzWeights) fewz_weight=fewz.getWeight(gen->vmass,gen->vpt,gen->vy);
 
      if (ientry<20) {
	printf("reweight=%4.2lf, fewz_weight=%4.2lf,dE_fsr=%+6.4lf\n",reweight,fewz_weight,(gen->mass-gen->vmass));
      }

      int iMassBinGenPreFsr = DYTools::findMassBin(gen->vmass);
      int iYBinGenPreFsr = DYTools::findAbsYBin(iMassBinGenPreFsr, gen->vy);
      int iMassBinGenPostFsr = DYTools::findMassBin(gen->mass);
      int iYBinGenPostFsr = DYTools::findAbsYBin(iMassBinGenPostFsr, gen->y);
      int idxGenPreFsr = DYTools::findIndexFlat(iMassBinGenPreFsr, iYBinGenPreFsr);
      int idxGenPostFsr = DYTools::findIndexFlat(iMassBinGenPostFsr, iYBinGenPostFsr);

      double fullGenWeight = reweight * scale * gen->weight * fewz_weight;

      fsr.fillIni(iMassBinGenPreFsr,iYBinGenPreFsr, fullGenWeight);
      fsr.fillFin(iMassBinGenPostFsr,iYBinGenPostFsr, fullGenWeight);
      if (validFlatIndices(idxGenPreFsr, idxGenPostFsr)) {
	fsr.fillMigration(idxGenPreFsr,idxGenPostFsr, fullGenWeight);
      }
 
      if ( DYTools::goodEtaPair(gen->veta_1, gen->veta_2) &&
	   DYTools::goodEtPair(gen->vpt_1, gen->vpt_2) ) {
	if (validFlatIndex(idxGenPreFsr)) {
	  fsrDET.fillIni(iMassBinGenPreFsr,iYBinGenPreFsr,fullGenWeight);
	}

	if ( DYTools::goodEtaPair(gen->eta_1, gen->eta_2) &&
	     DYTools::goodEtPair(gen->pt_1, gen->pt_2) ) {
	  if (validFlatIndex(idxGenPostFsr)) {
	    fsrDET.fillFin(iMassBinGenPostFsr,iYBinGenPostFsr,fullGenWeight);
	  }
	  /*
	  if ((iMassBinGenPostFsr>=0) && (iYBinGenPostFsr>=0)) {
	    if ((iMassBinGenPostFsr==1) && (iYBinGenPostFsr==0)) {
	      printf("adding %6.4lf for ientry=%d\n",fullGenWeight,ientry);
	    }
	  }
	  */
	  if (validFlatIndices(idxGenPreFsr,idxGenPostFsr)) {
	    fsrDET.fillMigration(idxGenPreFsr,idxGenPostFsr, fullGenWeight);

	    fsrDETexact.fillIni(iMassBinGenPreFsr,iYBinGenPreFsr,fullGenWeight);
	    fsrDETexact.fillFin(iMassBinGenPostFsr,iYBinGenPostFsr,fullGenWeight);
	    fsrDETexact.fillMigration(idxGenPreFsr,idxGenPostFsr, fullGenWeight);
	  }
	}
      }
	
      if( !(requiredTriggers.matchEventTriggerBit(info->triggerBits, 
						  info->runNum))) 
	continue;

      // loop through dielectrons
      dielectronArr->Clear();
      dielectronBr->GetEntry(ientry);    
      for(Int_t i=0; i<dielectronArr->GetEntriesFast(); i++) {

        const mithep::TDielectron *dielectron = (mithep::TDielectron*)((*dielectronArr)[i]);
	
	// Apply selection
	// Eta cuts
        if((fabs(dielectron->scEta_1)>DYTools::kECAL_GAP_LOW) && (fabs(dielectron->scEta_1)<DYTools::kECAL_GAP_HIGH)) continue;
        if((fabs(dielectron->scEta_2)>DYTools::kECAL_GAP_LOW) && (fabs(dielectron->scEta_2)<DYTools::kECAL_GAP_HIGH)) continue;
	if((fabs(dielectron->scEta_1) > 2.5)       || (fabs(dielectron->scEta_2) > 2.5))       continue;  // outside eta range? Skip to next event...
	
	// Asymmetric SC Et cuts
	if( ! ( ( dielectron->scEt_1 > DYTools::etMinLead  && dielectron->scEt_2 > DYTools::etMinTrail)
		|| ( dielectron->scEt_1 > DYTools::etMinTrail  && dielectron->scEt_2 > DYTools::etMinLead) )) continue;
    	
	// Both electrons must match trigger objects. At least one ordering
	// must match
	if( ! requiredTriggers.matchTwoTriggerObjectsAnyOrder( dielectron->hltMatchBits_1,
							       dielectron->hltMatchBits_2,
							       info->runNum) ) continue;
	
	// *** Smurf ID is superseeded by new selection ***
// 	// The Smurf electron ID package is the same as used in HWW analysis
// 	// and contains cuts like VBTF WP80 for pt>20, VBTF WP70 for pt<10
// 	// with some customization, plus impact parameter cuts dz and dxy
// 	if(!passSmurf(dielectron)) continue;  

	// The selection below is for the EGM working points from spring 2012
	// recommended for both 2011 and 2012 data
	if(!passEGM2011(dielectron, WP_MEDIUM, info->rhoLowEta)) continue;  

        // We have a Z candidate! HURRAY! 

// 	// Apply extra smearing to MC reconstructed dielectron mass
// 	// to better resemble the data
// 	// In systematics mode, use randomized MC smear factors
	double smearingCorrection = (systematicsMode == DYTools::RESOLUTION_STUDY) ?
          escale.generateMCSmearRandomized(dielectron->scEta_1,dielectron->scEta_2) :
          escale.generateMCSmear(dielectron->scEta_1,dielectron->scEta_2);
	double massResmeared = dielectron->mass + smearingCorrection;

	hZMassv[ifile]->Fill(massResmeared,scale * gen->weight);

	//
	// Fill structures for response matrix and bin by bin corrections
	// Note: there is no handling of overflow, underflow at present,
	// those entries are just dropped. This can be improved.
	// The only possible cases are: underflow in mass and overflow in Y.

	// Fill the matrix of post-FSR generator level invariant mass and rapidity
	detResponse.fillIni( iMassBinGenPostFsr, iYBinGenPostFsr, fullGenWeight );

	// Fill the matrix of the reconstruction level mass and rapidity
	int iMassReco = DYTools::findMassBin(massResmeared);
	int iYReco = DYTools::findAbsYBin(iMassReco, dielectron->y);
	detResponse.fillFin( iMassReco, iYReco, fullGenWeight );

	double shape_weight = 1.0;
	if( shapeWeights && iMassReco != -1 && iYReco != -1) {
	    shape_weight = (*shapeWeights)[iMassReco][iYReco];
	    //std::cout << "massResmeared=" << massResmeared << ", iMassReco=" << iMassReco << ", shapeWeight=" << shape_weight << "\n";
	}

	
        // Unlike the mass vs Y reference yields matrices, to prepare the
	// migration matrix we flatten (mass,Y) into a 1D array, and then
	// store (mass,Y in 1D)_gen vs (mass,Y in 1D)_rec
	int iIndexFlatGen  = DYTools::findIndexFlat(iMassBinGenPostFsr, iYBinGenPostFsr);
 	int iIndexFlatReco = DYTools::findIndexFlat(iMassReco, iYReco);
	if ( validFlatIndices(iIndexFlatGen, iIndexFlatReco) ) {
	  double fullWeight = reweight * scale * gen->weight * shape_weight;
	  //std::cout << "adding DetMig(" << iIndexFlatGen << "," << iIndexFlatReco << ") = " << reweight << "*" << scale << "*" << gen->weight << "*" << shape_weight << " = "  << (reweight * scale * gen->weight * shape_weight) << "\n";
	  detResponse.fillMigration(iIndexFlatGen, iIndexFlatReco, fullWeight );
	  detResponseExact.fillIni( iMassBinGenPostFsr, iYBinGenPostFsr, fullGenWeight );
	  detResponseExact.fillFin( iMassReco, iYReco, fullGenWeight );
	  detResponseExact.fillMigration(iIndexFlatGen, iIndexFlatReco, fullWeight );
	}

        Bool_t isB1 = DYTools::isBarrel(dielectron->scEta_1);
        Bool_t isB2 = DYTools::isBarrel(dielectron->scEta_2);

	hMassDiff->Fill(massResmeared - gen->mass);
	if( isB1 && isB2 )
	  hMassDiffBB->Fill(massResmeared - gen->mass);
	if( (isB1 && !isB2) || (!isB1 && isB2) )
	  hMassDiffEB->Fill(massResmeared - gen->mass);
	if( !isB1 && !isB2 )
	  hMassDiffEE->Fill(massResmeared - gen->mass);
	
	hMassDiffV->Fill(iIndexFlatGen, massResmeared - gen->mass);
	hYDiffV   ->Fill(iIndexFlatGen, dielectron->y - gen->y);
// 	if(iIndexFlatGen != -1){
// 	  hMassDiffV[iIndexFlatGen]->Fill(massResmeared - gen->mass);
// 	}

      } // end loop over dielectrons

    } // end loop over events 
    delete infile;
    infile=0, eventTree=0;
  } // end loop over files
  } 
  delete gen;

  //return;

  if (debugMode!=-1) {
  // Compute the errors on the elements of migration matrix
  detResponse.finalizeDetMigrationErr();
  detResponseExact.finalizeDetMigrationErr();
  fsr.finalizeDetMigrationErr();
  fsrDET.finalizeDetMigrationErr();
  fsrDETexact.finalizeDetMigrationErr();

  // Find response matrix, which is simply the normalized migration matrix
  std::cout << "find response matrix" << std::endl;
  detResponse.computeResponseMatrix();
  detResponseExact.computeResponseMatrix();
  fsr.computeResponseMatrix();
  fsrDET.computeResponseMatrix();
  fsrDETexact.computeResponseMatrix();

  std::cout << "find inverted response matrix" << std::endl;
  detResponse.invertResponseMatrix();
  detResponseExact.invertResponseMatrix();
  fsr.invertResponseMatrix();
  fsrDET.invertResponseMatrix();
  fsrDETexact.invertResponseMatrix();

  std::cout << "prepare flat-index arrays" << std::endl;
  detResponse.prepareFIArrays();
  detResponseExact.prepareFIArrays();
  fsr.prepareFIArrays();
  fsrDET.prepareFIArrays();
  fsrDETexact.prepareFIArrays();
  }

  //
  // Store constants and reference arrays in files
  //
  if (debugMode!=-1) std::cout << "store constants in a file" << std::endl;

  TString outputDir(TString("../root_files/constants/")+dirTag);
  if((systematicsMode==DYTools::RESOLUTION_STUDY) || (systematicsMode==DYTools::FSR_STUDY))
    outputDir = TString("../root_files/systematics/")+dirTag;
  gSystem->mkdir(outputDir,kTRUE);
  outputDir.Append("/");

  TString fnameTag="";
  {
    TString u="_";
    switch(systematicsMode) {
    case DYTools::NORMAL: 
      fnameTag=DYTools::analysisTag; 
      break;
    case DYTools::RESOLUTION_STUDY: 
      fnameTag=TString("_seed_") + DYTools::analysisTag + u;
      fnameTag+=seed;
      break;
    case DYTools::FSR_STUDY:
      fnameTag=TString("_reweight_") + DYTools::analysisTag + u;
      fnameTag+= int(100*reweightFsr);
      break;
    case DYTools::ESCALE_RESIDUAL:
      fnameTag=DYTools::analysisTag+TString("_escaleResidual");
      break;
    default:
      std::cout<<"requested mode not recognized when determining fnameTag"<<std::endl;
      assert(0);
    }
  }
  std::cout << "fnameTag=<" << fnameTag << ">\n";
  CPlot::sOutDir=TString("plots") + fnameTag;

  if (debugMode!=-1) {
    detResponse.autoSaveToFile(outputDir,fnameTag);  // detResonse, reference mc arrays
    detResponseExact.autoSaveToFile(outputDir,fnameTag);
    fsr.autoSaveToFile(outputDir,fnameTag);
    fsrDET.autoSaveToFile(outputDir,fnameTag);
    fsrDETexact.autoSaveToFile(outputDir,fnameTag);
  }
  else {
    if (!detResponse.autoLoadFromFile(outputDir,fnameTag) ||
	!detResponseExact.autoLoadFromFile(outputDir,fnameTag) ||
	!fsr.autoLoadFromFile(outputDir,fnameTag) ||
	!fsrDET.autoLoadFromFile(outputDir,fnameTag) ||
	!fsrDETexact.autoLoadFromFile(outputDir,fnameTag)) {
      std::cout << "loading failed\n";
      return;
    }
  }



  //--------------------------------------------------------------------------------------------------------------
  // Make plots 
  //==============================================================================================================  

  std::cout << "making plots" << std::endl;

  TString unfoldingConstFileName, yieldsFileName;
  detResponse.getFileNames(outputDir,fnameTag, unfoldingConstFileName, yieldsFileName);
  TString unfoldingConstantsPlotFName=unfoldingConstFileName;
  unfoldingConstantsPlotFName.Replace(unfoldingConstantsPlotFName.Index(".root"),
				      sizeof(".root"),
				      "_plots.root");
  TFile *fPlots=new TFile(unfoldingConstantsPlotFName,"recreate");
  if (!fPlots) {
    std::cout << "failed to create a file <" << unfoldingConstantsPlotFName << ">\n";
  }
 

  TCanvas *c = MakeCanvas("canvZmass1","canvZmass1",800,600);

  // string buffers
  char ylabel[50];   // y-axis label

  // 
  // Draw DY candidate mass at the reconstruction level. Extra
  // smearing is applied. This figure allows one to judge the 
  // correctness of the weights aplied to different samples from the
  // smoothness of the combined result.
  //
  sprintf(ylabel,"a.u. / %.1f GeV/c^{2}",hZMassv[0]->GetBinWidth(1));
  CPlot plotZMass1("zmass1","","m(ee) [GeV/c^{2}]",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    plotZMass1.AddHist1D(hZMassv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotZMass1.SetLogy();
  plotZMass1.Draw(c);
  SaveCanvas(c,"zmass1");
//   plotZMass1.Draw(c,doSave,format);
//   if (fPlots) { fPlots->cd(); c->Write(); }

  //
  // Draw a plot that illustrates the detector resolution effects.
  // We plot (gen-rec)/gen as a function of mass and rapidity.
  //
  TMatrixD resolutionEffect(DYTools::nMassBins,DYTools::nYBinsMax);
  resolutionEffect = 0;
  for(int i=0; i < resolutionEffect.GetNrows(); i++){
    for(int j=0; j < resolutionEffect.GetNcols(); j++){
      double ngen = (*detResponse.yieldsIni)(i,j);
      double nrec = (*detResponse.yieldsFin)(i,j);
      if( ngen != 0 )
	resolutionEffect(i,j) = (ngen-nrec)/ngen;
    }
  }
  resolutionEffect.Print();
  PlotMatrixVariousBinning(resolutionEffect, "resolution_effect", "LEGO2", NULL);

  //
  // Draw a plot that illustrates the losses due to reconstruction
  // We plot (preFsrExact-preFsr)/preFsrExact as a 
  // function of mass and rapidity.
  //
  TMatrixD *unfRecoEffect=detResponseExact.getReconstructionEffect(detResponse);
  unfRecoEffect->Print();
  PlotMatrixVariousBinning(*unfRecoEffect, "reconstruction_effect", "LEGO2", NULL);
  delete unfRecoEffect;

  TMatrixD *unfFsrDETRecoEffect=fsrDETexact.getReconstructionEffect(fsrDET);
  
  PlotMatrixVariousBinning(*unfFsrDETRecoEffect, "reconstruction_effect_fsrDET", "LEGO2", NULL);
  delete unfFsrDETRecoEffect;

  // Plot response and inverted response matrices
  //std::vector<TH2F*> hResponseV, hInvResponseV;
  //std::vector<TCanvas*> canvV;
  //std::vector<CPlot*> cpResponseV;

  //TH2F *hR, *hIR;
  //TCanvas *e2;
  //CPlot *cpR, *cpIR;
  detResponse.prepareHResponse();
  fsr.prepareHResponse();
  fsrDET.prepareHResponse();
  fsrDETexact.prepareHResponse();
   
  // Create a plot of detector resolution without mass binning
  TCanvas *g = MakeCanvas("canvMassDiff","canvMassDiff",600,600);
  CPlot plotMassDiff("massDiff","","reco mass - gen post-FSR mass [GeV/c^{2}]","a.u.");
  hMassDiffBB->Scale(1.0/hMassDiffBB->GetSumOfWeights());
  hMassDiffEB->Scale(1.0/hMassDiffEB->GetSumOfWeights());
  hMassDiffEE->Scale(1.0/hMassDiffEE->GetSumOfWeights());
  plotMassDiff.AddHist1D(hMassDiffBB,"EB-EB","hist",kBlack);
  plotMassDiff.AddHist1D(hMassDiffEB,"EE-EB","hist",kBlue);
  plotMassDiff.AddHist1D(hMassDiffEE,"EE-EE","hist",kRed);
  plotMassDiff.Draw(g);
  SaveCanvas(g,"massDiff");
//   if (fPlots) g->Write();

  // Create a plot of reco - gen post-FSR mass and rapidity difference 
  TCanvas *h1 = MakeCanvas("canvMassDiffV","canvMassDiffV",600,600);
  CPlot plotMassDiffV("massDiffV","",
		      "flat index",
		      "reco mass - gen post-FSR mass [GeV/c^{2}]");
  plotMassDiffV.AddHist2D(hMassDiffV,"LEGO");
  plotMassDiffV.Draw(h1);
  SaveCanvas(h1,"hMassDiffV");

  // Create a plot of reco - gen post-FSR mass and rapidity difference 
  TCanvas *h2 = MakeCanvas("canvYDiffV","canvYDiffV",600,600);
  CPlot plotYDiffV("massDiffV","",
		      "flat index",
		      "reco Y - gen post-FSR Y");
  plotYDiffV.AddHist2D(hYDiffV,"LEGO");
  plotYDiffV.Draw(h2);
  SaveCanvas(h2,"hYDiffV");

  if (fPlots) {
    fPlots->Close();
    delete fPlots;
    std::cout << "plots saved to a file <" << unfoldingConstantsPlotFName << ">\n";
  }

  //draw errors of Unfolding matrix
  TCanvas *cErrorsResp = MakeCanvas("cErrorsResp","detResponse.DetInvertedResponseErr", 600,600);
  detResponse.DetInvertedResponseErr->Draw("LEGO2");
  cErrorsResp->Update();
  SaveCanvas(cErrorsResp,"cErrorsResp");

  TCanvas *cFsrErrorsResp = MakeCanvas("cErrorsFsr","fsr__.DetInvertedResponseErr", 1200, 600);
  cFsrErrorsResp->Divide(2,1);
  cFsrErrorsResp->cd(1);
  fsr.DetInvertedResponseErr->Draw("LEGO2");
  cFsrErrorsResp->cd(2);
  fsrDET.DetInvertedResponseErr->Draw("LEGO2");
  SaveCanvas(cFsrErrorsResp,"cErrorsFsr");

  

  //--------------------------------------------------------------------------------------------------------------
  // Summary print out
  //==============================================================================================================
  cout << endl;
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  cout << endl; 

  detResponse.printConditionNumber();
  fsr.printConditionNumber();
  fsrDET.printConditionNumber();
  fsrDETexact.printConditionNumber();

  if (0) {
    //detResponse.printMatrices();
    //fsr.printMatrices();
    fsrDET.printMatrices();
  }

  //Print errors of the Unfolding matrix when they exceed 0.1
  /*
  for (int iM=0; iM<DYTools::nMassBins; iM++)
    for (int iY=0; iY<DYTools::nYBins[iM]; iY++)
      for (int jM=0; jM<DYTools::nMassBins; jM++)
        for (int jY=0; jY<DYTools::nYBins[jM]; jY++)
          {
	    int i=DYTools::findIndexFlat(iM,iY);
	    int j=DYTools::findIndexFlat(jM,jY);           
             if (DetInvertedResponseErr(i,j)>0.1)
                {
                   std::cout<<"DetInvertedResponseErr("<<i<<","<<j<<")="<<DetInvertedResponseErr(i,j);
                   std::cout<<", DetInvertedResponse("<<i<<","<<j<<")="<<DetInvertedResponse(i,j)<<std::endl;
                   std::cout<<"(iM="<<iM<<", iY="<<iY<<", jM="<<jM<<", jY="<<jY<<")"<<std::endl<<std::endl;
                }
             if (DetInvertedResponseErr2(i,j)>0.1)
                {
                   std::cout<<"DetInvertedResponseErr2("<<i<<","<<j<<")="<<DetInvertedResponseErr2(i,j);
                   std::cout<<", DetInvertedResponse("<<i<<","<<j<<")="<<DetInvertedResponse(i,j)<<std::endl;
                   std::cout<<"(iM="<<iM<<", iY="<<iY<<", jM="<<jM<<", jY="<<jY<<")"<<std::endl<<std::endl;
                }
          }
  */

  /*
  if (0) {
    // Printout of all constants, uncomment if needed
    //printf("DetCorrFactor:\n"); DetCorrFactor.Print();
    printf("DetMigration:\n"); DetMigration.Print();
    printf("DetResponse:\n"); DetResponse.Print();

    printf("DetInvertedResponse:\n"); DetInvertedResponse.Print();
    //printf("DetInvertedResponseErr:\n"); DetInvertedResponseErr.Print();
    //printf("DetResponseArr:\n"); DetResponseArr.Print();
    //printf("DetInvertedResponseArr:\n"); DetInvertedResponseArr.Print();
    //printf("DetInvertedResonseErrArr:\n"); DetInvertedResponseErrArr.Print();

    //   printf("Detector corr factor numerator:\n");
    //   DetCorrFactorNumerator.Print();

    printf("yieldsMcPostFsrGen:\n");
    yieldsMcPostFsrGen.Print();
    
    printf("yieldsMcPostFsrRec:\n");
    yieldsMcPostFsrRec.Print();


    //   printf("Detector corr factor denominator:\n");
    //   DetCorrFactorDenominator.Print();
    //   printf("yieldsMcPostFsrRecArr:\n");
    //   yieldsMcPostFsrRecArr.Print();
    
    //printf("yieldsMcGen:\n");
    //yieldsMcGen.Print();
  }
  */

  detResponse.printYields();
  fsr.printYields();
  fsrDET.printYields();

  gBenchmark->Show("makeUnfoldingMatrix");
}


//=== FUNCTION DEFINITIONS ======================================================================================
void computeNormalizedBinContent(double subset, double subsetErr,
				 double total, double totalErr,
				 double& ratio, double& ratioErr){
  
  if(total == 0) {
    printf("makeUnfoldingMatrix::Possible problem\n");
    printf("     empty column in the response matrix\n");
    return;
  }
  
  ratio = subset/total;

  // The formula for the ratio = subset/total is obtained by error
  // propagation. The subsetErr and totalErr are NOT assumed ot be
  // the sqrt(subset) and sqrt(total). (If one assume that, the formula
  // below reduces to the familiar sqrt(ratio*(1-ratio)/total) ).
  // The subset and subsetErr are part of the total and totalErr.
  // The formula is easiest to derive if we take "A +- dA" and
  // "B +- dB" as independent numbers, with total = A+B and
  // totalErr^2 = dA^2 + dB^2. One then does error propagation of
  // the expression ratio = A/(A+B).
  // The outcome of it is found below (the absolute error on the ratio)
  ratioErr = (1/total)*sqrt( subsetErr*subsetErr*(1-2*ratio)
			     + totalErr*totalErr*ratio*ratio );

  return;
}

void calculateInvertedMatrixErrors(const TMatrixD &T, 
	  const TMatrixD &TErrPos, const TMatrixD &TErrNeg,
				   TMatrixD &TinvErr){

  // Calculate errors on the inverted matrix by the Monte Carlo
  // method

  Double_t det;
  int nRow = T.GetNrows();
  int nCol = T.GetNcols();
  TMatrixD TinvSum(nRow,nCol);
  TMatrixD TinvSumSquares(nRow,nCol);

  // Reset Matrix where we will be accumulating RMS/sigma:
  TinvSum        = 0;
  TinvSumSquares = 0;
  TinvErr        = 0;

  // Do many tries, accumulate RMS
  int N = 10000;
  for(int iTry = 0; iTry<N; iTry++){
    // Find the smeared matrix
    TMatrixD Tsmeared = T;
    for(int i = 0; i<nRow; i++){
      for(int j = 0; j<nCol; j++){
	double central = T(i,j);
	double sigPos = TErrPos(i,j);
	double sigNeg = TErrNeg(i,j);
 	// Switch to symmetric errors: approximation, but much simpler
	double sig = (sigPos+sigNeg)/2.0;
	Tsmeared(i,j) = gRandom->Gaus(central,sig);
      }
    }
    // Find the inverted to smeared matrix
    TMatrixD TinvSmeared = Tsmeared;
    TinvSmeared.Invert(&det);
    // Accumulate sum and sum of squares for each element
    for(int i2 = 0; i2<nRow; i2++){
      for(int j2 = 0; j2<nCol; j2++){
	TinvSum       (i2,j2) += TinvSmeared(i2,j2);
	TinvSumSquares(i2,j2) += TinvSmeared(i2,j2)*TinvSmeared(i2,j2);
      }
    }
  }

  // Calculate the error matrix
  TMatrixD TinvAverage = TinvSum;
  for(int i = 0; i<nRow; i++){
    for(int j = 0; j<nCol; j++){
      TinvErr(i,j) = sqrt( TinvSumSquares(i,j)/double(N) 
			   - (TinvSum(i,j)/double(N))*(TinvSum(i,j)/double(N)) );
      TinvAverage(i,j) = TinvSum(i,j)/double(N);
    }
  }

  return;
}
