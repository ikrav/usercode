#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TStyle.h>
#include <TCanvas.h>                // class for drawing
#include <TH1F.h>                   // 1D histograms
#include <TGraphErrors.h>           // graphs
#include <TGraph2DErrors.h>           // graphs
#include <TProfile.h>               // profile histograms
#include <TBenchmark.h>             // class to track macro running statistics
#include <TLorentzVector.h>         // class for Lorentz vector computations
#include <TVectorD.h>
#include <TGraphErrors.h>
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
#include "../Include/TGenInfo.hh"

#include "../Include/EventSelector.hh"
#include "../Include/FEWZ.hh"
#include "../Include/UnfoldingTools.hh"
#include "../Include/InputFileMgr.hh"
#endif

//using std::cout;
//using std::endl;

template <class T>
inline T SQR(T x) { return (x)*(x); }

//=== FUNCTION DECLARATIONS ======================================================================================

double getErrorOnRatio(double nTotal, double nTotalErr, double nPass, double nPassErr) {
  double nFail = fabs(nTotal - nPass);
  //if (nTotalErr<nPassErr) {
  //  std::cout << "getErrorOnRatio: nTotalErr=" << nTotalErr << ", nPassErr=" << nPassErr << "\n";
  //}
  double nFailErr = sqrt( fabs(nTotalErr*nTotalErr - nPassErr*nPassErr ) );
  double err = sqrt( (nFail*nFail * nPassErr*nPassErr + 
		      nPass*nPass * nFailErr*nFailErr)
                     / (nTotal*nTotal*nTotal*nTotal) );
  return err;
}



// -----------------------------------------------------------------------

TH1F* plotXsec1D(const char *hname, const TMatrixD &nEvents, const TMatrixD &nEventsErr, CPlot *cp=NULL, int color=kBlack, int fineGrid=0) {
  TH1F *h=NULL;
  if (fineGrid) {
    int count=nEvents.GetNrows();
    if (nEvents.GetNcols()!=1) {
      std::cout << "plotXsec1D: fineGrid=1, but nEvents.GetNcols=" << nEvents.GetNcols() << "\n";
      return NULL;
    }
    h=new TH1F(hname,hname, count,0,count);
    h->SetDirectory(0);
    for (int iM=0; iM<count; ++iM) {
      double sum= nEvents[iM][0];
      double sumW2= nEventsErr[iM][0];
      h->SetBinContent(iM+1, sum);
      h->SetBinError(iM+1, sqrt(sumW2));
    }
  }
  else {
    h=new TH1F(hname,hname, DYTools::nMassBins, DYTools::massBinLimits);
    h->SetDirectory(0);
    for (int iM=0; iM<DYTools::nMassBins; ++iM) {
      double sum=0., sumW2=0.;
      for (int iY=0; iY<DYTools::nYBins[iM]; ++iY) {
	sum+= nEvents[iM][iY];
	sumW2+= nEventsErr[iM][iY];
      }
      h->SetBinContent(iM+1, sum);
      h->SetBinError(iM+1, sqrt(sumW2));
    }
  }

  if (cp) {
    cp->AddHist1D(h,TString(hname),"LP",color,1,0,1);
  }
  return h;
}

// -----------------------------------------------------------------------

TH1F* plotXsec2D(const char *name_base, int iM, const TMatrixD &nEvents, const TMatrixD &nEventsErr, CPlot *cp=NULL, int color=kBlack, int addText=0) {
  double massMin=DYTools::massBinLimits[iM];
  double massMax=DYTools::massBinLimits[iM+1];
  TString hname=Form("%s_mass_%2.0lf_%2.0lf",name_base,massMin,massMax);
  TH1F *h=new TH1F(hname,hname, DYTools::nYBins[iM], DYTools::yRangeMin,DYTools::yRangeMax);
  h->SetDirectory(0);
  for (int iY=0; iY<DYTools::nYBins[iM]; ++iY) {
    h->SetBinContent(iY+1, nEvents[iM][iY]);
    h->SetBinError(iY+1, sqrt(nEventsErr[iM][iY]));
  }
  if (cp) {
    cp->AddHist1D(h,TString(name_base),"L",color,1,0,1);
    if (addText) {
      TString mass_range=Form("%2.0lf<m<%2.0lf GeV",massMin,massMax);
      cp->AddTextBox(mass_range, 0.7,0.70,0.92,0.83, 0,kBlack,kWhite);
    }
  }
  return h;
}

// -----------------------------------------------------------------------


// -----------------------------------------------------------------------


// -----------------------------------------------------------------------

//=== MAIN MACRO =================================================================================================


void getXsecExtended(const TString mc_input, int debugMode=0, bool useFEWZ=true, int fineGrid=0)
{
  // check whether it is a calculation
  if (mc_input.Contains("_DebugRun_")) {
    std::cout << "getXsec: _DebugRun_ detected. Terminating the script\n";
    return;
  }

  if (debugMode==1) std::cout << "\n\n\tDEBUG MODE is ON\n\n";
  if (debugMode==-1) std::cout << "\n\n\tPLOT ONLY MODE is ON\n\n";
  std::cout << "DYTools::analysisTag=" << DYTools::analysisTag << "\n";

  if (DYTools::study2D && fineGrid) {
    std::cout << "fineGrid=1 is allowed only for study2D=0\n";
    return;
  }

  // normal calculation

  gBenchmark->Start("getXsec");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================
  
  //Bool_t doSave  = false;    // save plots?
  TString format = "png";   // output file format
  MCInputFileMgr_t inpMgr;
  TString fineGridStr;
  if (fineGrid==1) fineGridStr="fineGrid_";
  else if (fineGrid==2) fineGridStr="summer2011Grid_";
  else if (fineGrid==3) fineGridStr="summer2011specGrid_";
  
  Double_t massLow  = DYTools::massBinLimits[0];
  Double_t massHigh = DYTools::massBinLimits[DYTools::nMassBins];

  // fine grid: 0, 1, 2, ..., (fineGrid_1GeVstop-1), fineGrid_1GeVstop, fineGrid_1GeVstop+fineGridLargerStep, fineGrid_1GeVstop+2*fineGridLargerStep, ..., 1500
  //
  double fineGrid_1GeVstop=200.;
  double fineGrid_LargerStep=20.;

  int locMassBinCount=DYTools::nMassBins;
  double *massRangeEdges=NULL;
  if (!fineGrid) {
    massRangeEdges=new double[locMassBinCount+1];
    for (int i=0; i<=DYTools::nMassBins; ++i) {
      massRangeEdges[i]=DYTools::massBinLimits[i];
    }
  }
  else if (fineGrid==1) {
    // 1GeV grid
    locMassBinCount=int(fineGrid_1GeVstop-massLow+1e-3);
    std::cout << "1GeV grid size=" << locMassBinCount << "\n";
    // fineGrid_LargerStep
    double upperRange= (massHigh-fineGrid_1GeVstop);
    double upperCount= upperRange/fineGrid_LargerStep + 1e-3;
    if (upperCount - trunc(upperCount) > 1e-3) {
      std::cout << "(upperCount -1e-3)=" << (upperCount -1e-3) << "\n";
      std::cout << "should be integer value\n";
      return;
    }
    locMassBinCount += int(upperCount);
    std::cout << "mass grid[" << locMassBinCount << "]\n";

    massRangeEdges=new double[locMassBinCount+1];
    double m=massLow, dm=1;
    int count=0;
    while (m<=massHigh) {
      massRangeEdges[count]=m;
      if (1) {
	if (m<massHigh) std::cout << "count=" << count << ", massRange=" << m << " .. " << (m+dm) << "\n";
	else std::cout << "last edge=" << m << "\n";
      }
      count++;
      m+=dm;
      if (m==fineGrid_1GeVstop) dm=fineGrid_LargerStep;
    }
  }
  else if (fineGrid==2) {
    const int nMassBinTh=518;
    Double_t mxl = 14.0;
    locMassBinCount=nMassBinTh;
    massRangeEdges=new double[nMassBinTh+1];

    for(int iTh=0; iTh<=nMassBinTh; iTh++){
      // mass limits
      if     ( iTh >=   0 && iTh <  11 ) {mxl += 1.0;}
      else if( iTh >=  11 && iTh <  18 ) {mxl += 5.0;}
      else if( iTh >=  18 && iTh < 118 ) {mxl += 1.0;}
      else if( iTh >= 118 && iTh < 340 ) {mxl += 2.0;}
      else if( iTh >= 340 && iTh < nMassBinTh)   {mxl += 5.0; }
      else if( iTh == nMassBinTh)                {mxl = 1500; }
      massRangeEdges[iTh]=mxl;
    }
  }
  else if (fineGrid==3) {
    const int nMassBinTh=518;
    Double_t mxl = 14.0;
    locMassBinCount=nMassBinTh;
    massRangeEdges=new double[nMassBinTh+1];

    int iTh_new=0, tmpCounter=0;
    for(int iTh=0; iTh<=nMassBinTh; iTh++, iTh_new++){
      // mass limits
      if     ( iTh >=   0 && iTh <  11 ) {mxl += 1.0;}
      else if( iTh >=  11 && iTh <  18 ) {mxl += 5.0;}
      else if( iTh >=  18 && iTh < 118 ) {mxl += 1.0;}
      else if( iTh >= 118 && iTh < 340 ) {
	if (iTh>=139) {
	  std::cout << "iTh=" << iTh << ", current mxl=" << mxl << " +2\n";
	  if (iTh==139) tmpCounter=10;
	  tmpCounter--;
	  if (tmpCounter>0) iTh_new--;
	  else if (tmpCounter==0) {
	    tmpCounter=10;
	    std::cout << "iTh=" << iTh << ", iTh_new=" << iTh_new << ", current mxl=" << mxl << " +10\n";
	  }
	}
	mxl += 2.0;
      }
      else if( iTh >= 340 && iTh < nMassBinTh)   {
	std::cout << "iTh=" << iTh << ", current mxl=" << mxl << " +5\n";
	if (iTh<342) iTh_new--;
	else {
	  if ((iTh-342)%4>0) iTh_new--;
	  else std::cout << "iTh=" << iTh << ", iTh_new=" << iTh_new << ", current mxl=" << mxl << " +10\n";
	}
	mxl += 5.0; 
      }
      else if( iTh == nMassBinTh)                {mxl = 1500; }
      massRangeEdges[iTh_new]=mxl;
    }
    massRangeEdges[iTh_new-2]=1500.;
    std::cout << "iTh_new=" << iTh_new << "\n";
    locMassBinCount=iTh_new-2;
  }

  TVectorD massGrid(locMassBinCount+1);
  for (int i=0; i<=locMassBinCount; ++i) massGrid[i]=massRangeEdges[i];
  TH1F hMassIdx("h_massIdx","h_massIdx",locMassBinCount-1,massRangeEdges);
  //delete massRangeEdges;

  if (0) {
    printHisto(&hMassIdx);
    if (0) {
      for (int i=0; i<int(massHigh); ++i) {
	double m=i+0.5;
	std::cout << "i=" << i << ", m=" << m << ", idx=" << (hMassIdx.FindBin(m)-1) << "\n";
      }
    }
    return;
  }

  //std::cout << "massHigh=" << massHigh << "\n";
  //std::cout << "locMassBinCount=" << locMassBinCount << "\n";

  if (!inpMgr.Load(mc_input)) {
    return;
  }
  
  //--------------------------------------------------------------------------------------------------------------
  // Main code 
  //==============================================================================================================
  
  //  
  // Set up histograms
  //

  //vector<TH1D*> hZMassv;
  
  Double_t   nZv = 0;

  TMatrixD nEvents (locMassBinCount,DYTools::nYBinsMax);    // number of weigthed events
  TMatrixD nEventsDET (locMassBinCount,DYTools::nYBinsMax); // number of weighted events in the detector acceptance
  TMatrixD nEventsDETrecoPostIdx (locMassBinCount,DYTools::nYBinsMax);    // number of weigthed events
  TMatrixD w2Events (locMassBinCount,DYTools::nYBinsMax);
  TMatrixD w2EventsDET (locMassBinCount,DYTools::nYBinsMax);
  TMatrixD w2EventsDETrecoPostIdx (locMassBinCount,DYTools::nYBinsMax);
  Double_t nZpeak=0, w2Zpeak=0;
  double nZpeakDET=0, w2ZpeakDET=0;
  double nZpeakDETrecoPostIdx=0, w2ZpeakDETrecoPostIdx=0;

  nEvents = 0;
  w2Events = 0;
  nEventsDET = 0;
  w2EventsDET  = 0;
  nEventsDETrecoPostIdx = 0;
  w2EventsDETrecoPostIdx  = 0;

  //char hname[100];
  //for(UInt_t ifile = 0; ifile<fnamev.size(); ifile++) {
  //  sprintf(hname,"hZMass_%i",ifile); hZMassv.push_back(new TH1F(hname,"",500,0,500)); hZMassv[ifile]->Sumw2();
  //}

  // 
  // Read weights from a file
  //
  const bool useFewzWeights = useFEWZ;
  const bool cutZPT100 = true;
  FEWZ_t fewz(useFewzWeights,cutZPT100);
  if (useFewzWeights && !fewz.isInitialized()) {
    std::cout << "failed to prepare FEWZ correction\n";
    throw 2;
  }

  //
  // Access samples and fill histograms
  //  
  TFile *infile=0;
  TTree *eventTree=0;  
    
  // Data structures to store info from TTrees
  mithep::TGenInfo *gen  = new mithep::TGenInfo();

  // loop over samples  
  double lumi0=0;
  if (debugMode!=-1) {
  for(UInt_t ifile=0; ifile<inpMgr.fileNames().size(); ifile++) {
  
    // Read input file
    cout << "Processing " << inpMgr.fileName(ifile) << "..." << endl;
    infile = new TFile(inpMgr.fileName(ifile));
    assert(infile);

    // Get the TTrees
    eventTree = (TTree*)infile->Get("Events"); assert(eventTree);

    // Find weight for events for this file
    // The first file in the list comes with weight 1,
    // all subsequent ones are normalized to xsection and luminosity
    double lumi  = eventTree->GetEntries()/inpMgr.xsec(ifile);
    if (ifile==0) lumi0=lumi;
    double scale = lumi0/lumi;
    std::cout << "       -> sample weight is " << scale << endl;

    // Set branch address to structures that will store the info  
    eventTree->SetBranchAddress("Gen",&gen);
    TBranch *genBr = eventTree->GetBranch("Gen");
 
    // loop over events    
    nZv += scale * eventTree->GetEntries();

    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
      if (debugMode && (ientry>100000)) break;

      genBr->GetEntry(ientry);
      if (ientry%1000000==0) printProgress("ientry=",ientry,eventTree->GetEntriesFast());

      double massPreFsr = gen->vmass;   // pre-FSR
      double yPreFsr = gen->vy;    // pre-FSR
      double massPostFsr = gen->mass;   // post-FSR
      double yPostFsr = gen->y;    // post-FSR

      if ((massPreFsr < massLow) || (massPreFsr > massHigh)) continue;
      if ((fabs(yPreFsr) < DYTools::yRangeMin) || 
	  (fabs(yPreFsr) > DYTools::yRangeMax)) continue;

      int ibinMassPreFsr=-1, ibinYPreFsr=-1;
      int ibinMassPostFsr=-1, ibinYPostFsr=-1;

      if (!fineGrid) {
	ibinMassPreFsr = DYTools::findMassBin(massPreFsr);
	ibinMassPostFsr= DYTools::findMassBin(massPostFsr);
	ibinYPreFsr = DYTools::findAbsYBin(ibinMassPreFsr, yPreFsr);
	ibinYPostFsr= DYTools::findAbsYBin(ibinMassPostFsr, yPostFsr);
      }
      else {
	ibinMassPreFsr=hMassIdx.FindBin(massPreFsr)-1;
	ibinMassPostFsr=hMassIdx.FindBin(massPostFsr)-1;
	ibinYPreFsr=0;
	ibinYPostFsr=0;
	//printf("massPreFsr=%8.4lf, idx=%3d;  massPostFsr=%8.4lf, idx=%3d\n", massPreFsr,ibinMassPreFsr, massPostFsr,ibinMassPostFsr);
      }

      // We are only interested in the events, reconstructed with 
      // good mass and rapidity 
      if (ibinMassPreFsr==-1 || 
	  ibinMassPreFsr>=locMassBinCount || 
	  ibinYPreFsr==-1) {
	printf(".. skipping mass=%6.4lf, y=%6.4lf. ibinMass=%d, ibinY=%d\n",massPreFsr,yPreFsr,ibinMassPreFsr,ibinYPreFsr);
	continue;
      }

      // Find FEWZ-powheg reweighting factor 
      // that depends on pre-FSR Z/gamma* rapidity, pt, and mass
      double fewz_weight = 1.0;

      if(useFewzWeights) {
	fewz_weight=fewz.getWeight(gen->vmass,gen->vpt,gen->vy);
      }

      //std::cout << "weight=scale*gen->weight*fewz: " << scale << " * " << gen->weight << " * " << fewz_weight << "\n";
      double fullWeight = scale * gen->weight * fewz_weight;
      nEvents(ibinMassPreFsr,ibinYPreFsr) += fullWeight;
      w2Events(ibinMassPreFsr,ibinYPreFsr) += fullWeight*fullWeight;

      // Pre-FSR cross section
      if( DYTools::goodEtPair(gen->vpt_1, gen->vpt_2) &&
	  DYTools::goodEtaPair(gen->veta_1, gen->veta_2) ) {
	nEventsDET(ibinMassPreFsr,ibinYPreFsr) += fullWeight;
	w2EventsDET(ibinMassPreFsr,ibinYPreFsr) += fullWeight*fullWeight;
      }

      // Post-FSR cross section
      if(  (ibinMassPostFsr!=-1) && (ibinYPostFsr!=-1) &&
	   DYTools::goodEtPair(gen->pt_1, gen->pt_2) &&
	   DYTools::goodEtaPair(gen->eta_1, gen->eta_2) ) {
	nEventsDETrecoPostIdx(ibinMassPostFsr,ibinYPostFsr) += fullWeight;
	w2EventsDETrecoPostIdx(ibinMassPostFsr,ibinYPostFsr) += fullWeight*fullWeight;
      }
    }
    delete infile;
    infile=0, eventTree=0;
  }
  delete gen;

  // Determine Z-peak event count
  for (int i=0; i<locMassBinCount; i++) {
    int isZpeak=0;
    // bin idx is (i+1)
    if ((hMassIdx.GetBinLowEdge(i+1)>=60-1e-3) 
	&& (hMassIdx.GetBinLowEdge(i+1+1)<=120+1e-3)) isZpeak=1;
    if (isZpeak) {
      int yiMax=(fineGrid) ? 1:DYTools::nYBins[i];
      for (int yi=0; yi<yiMax; ++yi) {
	nZpeak += nEvents(i,yi);
	w2Zpeak += w2Events(i,yi);
	nZpeakDET += nEventsDET(i,yi);
	w2ZpeakDET += w2EventsDET(i,yi);
	nZpeakDETrecoPostIdx += nEventsDETrecoPostIdx(i,yi);
	w2ZpeakDETrecoPostIdx += w2EventsDETrecoPostIdx(i,yi);
      }
    }
  }
  std::cout << "\n";
  std::cout << "nZpeak=" << nZpeak << ", w2Zpeak=" << w2Zpeak << "\n";
  std::cout << "nZpeakDET=" << nZpeakDET << ", w2ZpeakDET=" << w2ZpeakDET << "\n";
  std::cout << "nZpeakDETrecoPostIdx=" << nZpeakDETrecoPostIdx << ", w2ZpeakDETrecoPostIdx=" << w2ZpeakDETrecoPostIdx << "\n";
  std::cout << "\n";


  if (nZpeak==0) {
    std::cout << "no events in the Z-peak region\n";
    return ;
  }
  } // if (debugMode!=-1)


  // Containers of the normalized event counts
  TMatrixD nEventsNorm (locMassBinCount,DYTools::nYBinsMax);    // number of weigthed events, normalized to Z-peak
  TMatrixD nEventsDETNorm (locMassBinCount,DYTools::nYBinsMax); // number of weighted events in the detector acceptance, normalized to Z-peak
  TMatrixD nEventsDETrecoPostIdxNorm (locMassBinCount,DYTools::nYBinsMax); // number of weighted events in the detector acceptance, normalized to Z-peak
  TMatrixD nEventsNormErr (locMassBinCount,DYTools::nYBinsMax);
  TMatrixD nEventsDETNormErr (locMassBinCount,DYTools::nYBinsMax);
  TMatrixD nEventsDETrecoPostIdxNormErr (locMassBinCount,DYTools::nYBinsMax);

  TMatrixD nEventsErr (locMassBinCount,DYTools::nYBinsMax);
  TMatrixD nEventsDETErr (locMassBinCount,DYTools::nYBinsMax);
  TMatrixD nEventsDETrecoPostIdxErr (locMassBinCount,DYTools::nYBinsMax);

  nEventsNorm=0;
  nEventsDETNorm=0;
  nEventsDETrecoPostIdxNorm=0;
  nEventsNormErr=0;
  nEventsDETNormErr=0;
  nEventsDETrecoPostIdxNormErr=0;

  nEventsErr=0;
  nEventsDETErr=0;

  if (debugMode!=-1) {
  for(int i=0; i<locMassBinCount; i++) {
    int yiMax=(fineGrid) ? 1:DYTools::nYBins[i];
    for (int j=0; j<yiMax; j++) {
      nEventsErr(i,j)=sqrt(w2Events(i,j));
      nEventsDETErr(i,j)=sqrt(w2EventsDET(i,j));
      nEventsDETrecoPostIdxErr(i,j)=sqrt(w2EventsDETrecoPostIdx(i,j));

      nEventsNorm(i,j) = nEvents(i,j)/nZpeak;
      nEventsNormErr(i,j) = 
	getErrorOnRatio(nEvents(i,j),nEventsErr(i,j),
			nZpeak, sqrt(w2Zpeak));

      nEventsDETNorm(i,j) = nEventsDET(i,j)/nZpeakDET;
      nEventsDETNormErr(i,j) =
	getErrorOnRatio(nEventsDET(i,j),nEventsDETErr(i,j),
			nZpeakDET, sqrt(w2ZpeakDET));

      nEventsDETrecoPostIdxNorm(i,j) = nEventsDETrecoPostIdx(i,j)/nZpeakDETrecoPostIdx;
      nEventsDETrecoPostIdxNormErr(i,j) =
	getErrorOnRatio(nEventsDETrecoPostIdx(i,j),nEventsDETrecoPostIdxErr(i,j),
			nZpeakDETrecoPostIdx, sqrt(w2ZpeakDETrecoPostIdx));
    }
  }
  }



  TString outFile= TString("../root_files/xSecThExt_");
  //outFile.Append("2MCfiles_");
  if (!useFewzWeights) outFile.Append("noFEWZ_");
  if (fineGridStr.Length()) outFile.Append(fineGridStr);
  if (debugMode==1) outFile.Append("debug_");
  outFile.Append( DYTools::analysisTag + TString("_tmp.root") );


  TVectorD rapidityGrid(massGrid.GetNoElements()-1);
  rapidityGrid=1;

  if (debugMode!=-1) {
    TFile thFile(outFile,"recreate");
    massGrid.Write("massBinEdges");
    if (fineGrid) rapidityGrid.Write("rapidityBinCount");
    nEvents.Write("nGenEvents");
    nEventsErr.Write("nGenEventsErr");
    nEventsDET.Write("nGenEventsDET");
    nEventsDETErr.Write("nGenEventsDETErr");
    nEventsDETrecoPostIdx.Write("nGenEventsDETrecoPostIdx");
    nEventsDETrecoPostIdxErr.Write("nGenEventsDETRecoPostIdxErr");
    nEventsNorm.Write("nGenEventsNorm");
    nEventsNormErr.Write("nGenEventsNormErr");
    nEventsDETNorm.Write("nGenEventsDETNorm");
    nEventsDETNormErr.Write("nGenEventsDETNormErr");
    nEventsDETrecoPostIdxNorm.Write("nGenEventsDETrecoPostIdxNorm");
    nEventsDETrecoPostIdxNormErr.Write("nGenEventsDETrecoPostIdxNormErr");
    TVectorD zPeakInfo(6);
    zPeakInfo(0)=nZpeak; zPeakInfo(1)=sqrt(w2Zpeak);
    zPeakInfo(2)=nZpeakDET; zPeakInfo(3)=sqrt(w2ZpeakDET);
    zPeakInfo(4)=nZpeakDETrecoPostIdx; zPeakInfo(5)=sqrt(w2ZpeakDETrecoPostIdx);
    zPeakInfo.Write("zPeakCountAndErr");
    thFile.Close();
    std::cout << "file <" << outFile << "> created\n";
  }
  else {
    TFile thFile(outFile);
    nEvents.Read("nGenEvents");
    nEventsErr.Read("nGenEventsErr");
    nEventsDET.Read("nGenEventsDET");
    nEventsDETErr.Read("nGenEventsDETErr");
    nEventsDETrecoPostIdx.Read("nGenEventsDETrecoPostIdx");
    nEventsDETrecoPostIdxErr.Read("nGenEventsDETRecoPostIdxErr");
    nEventsNorm.Read("nGenEventsNorm");
    nEventsNormErr.Read("nGenEventsNormErr");
    nEventsDETNorm.Read("nGenEventsDETNorm");
    nEventsDETNormErr.Read("nGenEventsDETNormErr");
    nEventsDETrecoPostIdxNorm.Read("nGenEventsDETrecoPostIdxNorm");
    nEventsDETrecoPostIdxNormErr.Read("nGenEventsDETrecoPostIdxNormErr");
    TVectorD zPeakInfo(6);
    zPeakInfo.Read("zPeakCountAndErr");
    thFile.Close();
    
    nZpeak=zPeakInfo[0]; w2Zpeak=SQR(zPeakInfo[1]);
    nZpeakDET=zPeakInfo[2]; w2ZpeakDET=SQR(zPeakInfo[3]);
    nZpeakDETrecoPostIdx=zPeakInfo[4]; w2ZpeakDETrecoPostIdx=SQR(zPeakInfo[5]);
    std::cout << "file <" << outFile << "> loaded\n";
  }

  //--------------------------------------------------------------------------------------------------------------
  // Make plots 
  //==============================================================================================================  
  CPlot::sOutDir="plots" + DYTools::analysisTag;

  TString fileNamePlots=outFile;
  fileNamePlots.ReplaceAll("_tmp.root","_plots_tmp.root");
  TFile *filePlots=new TFile(fileNamePlots,"recreate");
  if (!filePlots) {
    std::cout << "failed to create file <" << fileNamePlots << ">\n";
    throw 2;
  }
  // string buffers
  //char ylabel[50];   // y-axis label


  if (DYTools::study2D) {
    TString c2Dname="canvXsectTh_2D";
    if (useFewzWeights) c2Dname.Append("_FEWZ");
    TCanvas *c2D = MakeCanvas(c2Dname,c2Dname,600,600+300*DYTools::study2D);
    std::vector<TH1F*> hXsecV;
    std::vector<CPlot*> cpV;
    hXsecV.reserve(DYTools::nMassBins);
    cpV.reserve(DYTools::nMassBins);
    for (int iM=0; iM<DYTools::nMassBins; ++iM) {
      double massMin=DYTools::massBinLimits[iM];
      double massMax=DYTools::massBinLimits[iM+1];
      CPlot *cp=new CPlot(Form("cpMass%2.0lf_%2.0lf",massMin,massMax),
			  "","|Y|","counts");
      cpV.push_back(cp);
      hXsecV.push_back( plotXsec2D("xsec",iM,nEvents,nEventsErr,cp,kBlack,1) );
    }
    c2D->Divide(2,3);
    for (int ipad=1; ipad<=6; ++ipad) {
      cpV[ipad]->Draw(c2D,false,"png",ipad);
    }
    c2D->Update();
    SaveCanvas(c2D,c2Dname);
  }

  {
    TString canvName=TString("canvXsectTh_") + fineGridStr + TString("1D");
    if (useFewzWeights) canvName.Append("_FEWZ");
    TCanvas *c1D = MakeCanvas(canvName,canvName,600,600);
    CPlot *cp=new CPlot("cplot_1D","","M_{ee} (GeV)","counts");
    TH1F *hXsec= plotXsec1D("xsec1D",nEvents,nEventsErr,cp,kBlue, fineGrid);
    hXsec->SetDirectory(0);
    cp->SetLogx();
    cp->Draw(c1D,false,"png");
    c1D->Update();
    SaveCanvas(c1D,canvName);
  }
	
  {
    TString canvName=TString("canvXsectThNorm_") + fineGridStr + TString("1D");
    if (useFewzWeights) canvName.Append("_FEWZ");
    TCanvas *c1D = MakeCanvas(canvName,canvName,600,600);
    CPlot *cp=new CPlot("cplotNorm_1D","","M_{ee} (GeV)","normalized counts");
    TH1F *hXsec= plotXsec1D("xsec1D",nEventsNorm,nEventsNormErr,cp,kBlue, fineGrid);
    hXsec->SetDirectory(0);
    cp->SetLogx();
    cp->Draw(c1D,false,"png");
    c1D->Update();
    SaveCanvas(c1D,canvName);
  }
	
  /* 

 // Z mass
  sprintf(ylabel,"a.u. / %.1f GeV/c^{2}",hZMassv[0]->GetBinWidth(1));
  CPlot plotZMass1("zmass1","","m(Z) [GeV/c^{2}]",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    plotZMass1.AddHist1D(hZMassv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotZMass1.SetLogy();
  plotZMass1.Draw(c);
  SaveCanvas(c, "zmass1");

  PlotMatrixVariousBinning(accv, "acceptance", "LEGO2",filePlots);
  filePlots->Close();
  if (DYTools::study2D==0)
    Plot1D(accv,accErrv,"acceptance1D","acceptance");
  //delete filePlots;
  
  */
  //--------------------------------------------------------------------------------------------------------------
  // Summary print out
  //==============================================================================================================

  if (!fineGrid) {
    const int printSystErr=0;
    TMatrixD zeroErr=nEventsErr; zeroErr=0;
    printYields("nGenEvents", nEvents,nEventsErr,zeroErr, printSystErr);
    printYields("nGenEventsDET", nEventsDET,nEventsDETErr,zeroErr, printSystErr);
    printYields("nGenEventsDETrecoPostIdx",nEventsDETrecoPostIdx,nEventsDETrecoPostIdxErr, zeroErr,printSystErr);
  }

  gBenchmark->Show("getXsec");
}


//=== FUNCTION DEFINITIONS ======================================================================================

//--------------------------------------------------------------------------------------------------
