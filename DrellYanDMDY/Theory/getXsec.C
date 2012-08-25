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
  if (nTotalErr<nPassErr) {
    std::cout << "getErrorOnRatio: nTotalErr=" << nTotalErr << ", nPassErr=" << nPassErr << "\n";
  }
  double nFailErr = sqrt( fabs(nTotalErr*nTotalErr - nPassErr*nPassErr ) );
  double err = sqrt( (nFail*nFail * nPassErr*nPassErr + 
		      nPass*nPass * nFailErr*nFailErr)
                     / (nTotal*nTotal*nTotal*nTotal) );
  return err;
}

// -----------------------------------------------------------------------

TH1F* plotXsec1D(const char *hname, const TMatrixD &nEvents, const TMatrixD &nEventsErr, CPlot *cp=NULL, int color=kBlack) {
  TH1F *h=new TH1F(hname,hname, DYTools::nMassBins, DYTools::massBinLimits);
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


//=== MAIN MACRO =================================================================================================


void getXsec(const TString mc_input, int debugMode=0, int plotOnly=0)
{
  // check whether it is a calculation
  if (mc_input.Contains("_DebugRun_")) {
    std::cout << "getXsec: _DebugRun_ detected. Terminating the script\n";
    return;
  }

  if (debugMode) std::cout << "\n\n\tDEBUG MODE is ON\n\n";
  std::cout << "DYTools::analysisTag=" << DYTools::analysisTag << "\n";

  // normal calculation

  gBenchmark->Start("getXsec");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================
  
  //Bool_t doSave  = false;    // save plots?
  TString format = "png";   // output file format
  MCInputFileMgr_t inpMgr;
  
  Double_t massLow  = DYTools::massBinLimits[0];
  Double_t massHigh = DYTools::massBinLimits[DYTools::nMassBins];

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

  TMatrixD nEvents (DYTools::nMassBins,DYTools::nYBinsMax);    // number of weigthed events
  TMatrixD nEventsDET (DYTools::nMassBins,DYTools::nYBinsMax); // number of weighted events in the detector acceptance
  TMatrixD w2Events (DYTools::nMassBins,DYTools::nYBinsMax);
  TMatrixD w2EventsDET (DYTools::nMassBins,DYTools::nYBinsMax);
  Double_t nZpeak=0, w2Zpeak=0;
  double nZpeakDET=0, w2ZpeakDET=0;
  
  nEvents = 0;
  w2Events = 0;
  nEventsDET = 0;
  w2EventsDET  = 0;

  //char hname[100];
  //for(UInt_t ifile = 0; ifile<fnamev.size(); ifile++) {
  //  sprintf(hname,"hZMass_%i",ifile); hZMassv.push_back(new TH1F(hname,"",500,0,500)); hZMassv[ifile]->Sumw2();
  //}

  // 
  // Read weights from a file
  //
  const bool useFewzWeights = false;
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
  if (plotOnly==0) {
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
    cout << "       -> sample weight is " << scale << endl;

    // Set branch address to structures that will store the info  
    eventTree->SetBranchAddress("Gen",&gen);
    TBranch *genBr = eventTree->GetBranch("Gen");
 
    // loop over events    
    nZv += scale * eventTree->GetEntries();

    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
      if (debugMode && (ientry>100000)) break;

      genBr->GetEntry(ientry);
      if (ientry%1000000==0) {
	double r=trunc(ientry/double(eventTree->GetEntries())*1000)*0.1;
	std::cout << "ientry=" << ientry << "/" << eventTree->GetEntries() << " (" << r << "%)\n";
      }

      double massPreFsr = gen->vmass;   // pre-FSR
      double yPreFsr = gen->vy;    // pre-FSR

      if ((massPreFsr < massLow) || (massPreFsr > massHigh)) continue;
      if ((fabs(yPreFsr) < DYTools::yRangeMin) || 
	  (fabs(yPreFsr) > DYTools::yRangeMax)) continue;

      int ibinMassPreFsr = DYTools::findMassBin(massPreFsr);
      int ibinYPreFsr = DYTools::findAbsYBin(ibinMassPreFsr, yPreFsr);

      // We are only interested in the events, reconstructed with 
      // good mass and rapidity 
      if (ibinMassPreFsr==-1 || 
	  ibinMassPreFsr>=DYTools::nMassBins || 
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

      // Asymmetric Et cut scheme for DY analysis
      if( ( (gen->pt_1 > DYTools::etMinLead && gen->pt_2 > DYTools::etMinTrail) 
         || (gen->pt_1 > DYTools::etMinTrail && gen->pt_2 > DYTools::etMinLead) )
	  && ((fabs(gen->eta_1)<DYTools::kECAL_GAP_LOW) || (fabs(gen->eta_1)>DYTools::kECAL_GAP_HIGH))
         && ((fabs(gen->eta_2)<DYTools::kECAL_GAP_LOW) || (fabs(gen->eta_2)>DYTools::kECAL_GAP_HIGH))   
	  && (fabs(gen->eta_1)<2.5) && (fabs(gen->eta_2)<2.5)) {
        
	nEventsDET(ibinMassPreFsr,ibinYPreFsr) += fullWeight;
	w2EventsDET(ibinMassPreFsr,ibinYPreFsr) += fullWeight*fullWeight;
      }
    }
    delete infile;
    infile=0, eventTree=0;
  }
  delete gen;

  // Determine Z-peak event count
  for (int i=0; i<DYTools::nMassBins; i++) {
    if ( (DYTools::massBinLimits[i]>=60-1e-3) 
	 && (DYTools::massBinLimits[i+1]<=120+1e-3)) {
      for (int yi=0; yi<DYTools::nYBins[i]; ++yi) {
	nZpeak += nEvents(i,yi);
	w2Zpeak += w2Events(i,yi);
	nZpeakDET += nEventsDET(i,yi);
	w2ZpeakDET += w2EventsDET(i,yi);
      }
    }
  }
  std::cout << "\n";
  std::cout << "nZpeak=" << nZpeak << ", w2Zpeak=" << w2Zpeak << "\n";
  std::cout << "nZpeakDET=" << nZpeakDET << ", w2ZpeakDET=" << w2ZpeakDET << "\n";
  std::cout << "\n";


  if (nZpeak==0) {
    std::cout << "no events in the Z-peak region\n";
    return ;
  }
  } // if (plotOnly==0)


  // Containers of the normalized event counts
  TMatrixD nEventsNorm (DYTools::nMassBins,DYTools::nYBinsMax);    // number of weigthed events, normalized to Z-peak
  TMatrixD nEventsDETNorm (DYTools::nMassBins,DYTools::nYBinsMax); // number of weighted events in the detector acceptance, normalized to Z-peak
  TMatrixD nEventsNormErr (DYTools::nMassBins,DYTools::nYBinsMax);
  TMatrixD nEventsDETNormErr (DYTools::nMassBins,DYTools::nYBinsMax);

  TMatrixD nEventsErr (DYTools::nMassBins,DYTools::nYBinsMax);
  TMatrixD nEventsDETErr (DYTools::nMassBins,DYTools::nYBinsMax);

  nEventsNorm=0;
  nEventsDETNorm=0;
  nEventsNormErr=0;
  nEventsDETNormErr=0;

  nEventsErr=0;
  nEventsDETErr=0;

  if (plotOnly==0) {
  for(int i=0; i<DYTools::nMassBins; i++) {
    for (int j=0; j<DYTools::nYBins[i]; j++) {
      nEventsErr(i,j)=sqrt(w2Events(i,j));
      nEventsDETErr(i,j)=sqrt(w2EventsDET(i,j));

      nEventsNorm(i,j) = nEvents(i,j)/nZpeak;
      nEventsNormErr(i,j) = 
	getErrorOnRatio(nEvents(i,j),nEventsErr(i,j),
			nZpeak, sqrt(w2Zpeak));

      nEventsDETNorm(i,j) = nEventsDET(i,j)/nZpeakDET;
      nEventsDETNormErr(i,j) =
	getErrorOnRatio(nEventsDET(i,j),nEventsDETErr(i,j),
			nZpeakDET, sqrt(w2ZpeakDET));
    }
  }
  }



  TString outFile= TString("../root_files/xSecTh_");
  if (!useFewzWeights) outFile.Append("noFEWZ_");
  outFile.Append( DYTools::analysisTag + TString("_tmp.root") );

  if (plotOnly==0) {
    TFile thFile(outFile,"recreate");
    nEvents.Write("nGenEvents");
    nEventsErr.Write("nGenEventsErr");
    nEventsDET.Write("nGenEventsDET");
    nEventsDETErr.Write("nGenEventsDETErr");
    nEventsNorm.Write("nGenEventsNorm");
    nEventsNormErr.Write("nGenEventsNormErr");
    nEventsDETNorm.Write("nGenEventsDETNorm");
    nEventsDETNormErr.Write("nGenEventsDETNormErr");
    TVectorD zPeakInfo(4);
    zPeakInfo(0)=nZpeak; zPeakInfo(1)=sqrt(w2Zpeak);
    zPeakInfo(2)=nZpeakDET; zPeakInfo(3)=sqrt(w2ZpeakDET);
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
    nEventsNorm.Read("nGenEventsNorm");
    nEventsNormErr.Read("nGenEventsNormErr");
    nEventsDETNorm.Read("nGenEventsDETNorm");
    nEventsDETNormErr.Read("nGenEventsDETNormErr");
    TVectorD zPeakInfo(4);
    zPeakInfo.Read("zPeakCountAndErr");
    thFile.Close();
    
    nZpeak=zPeakInfo[0]; w2Zpeak=SQR(zPeakInfo[1]);
    nZpeakDET=zPeakInfo[2]; w2ZpeakDET=SQR(zPeakInfo[3]);
    std::cout << "file <" << outFile << "> loaded\n";
  }

  //--------------------------------------------------------------------------------------------------------------
  // Make plots 
  //==============================================================================================================  
  CPlot::sOutDir="plots" + DYTools::analysisTag;

  TString fileNamePlots=TString("../root_files/xSecTh_") + DYTools::analysisTag + TString("_plots_tmp.root");
  TFile *filePlots=new TFile(fileNamePlots,"recreate");
  if (!filePlots) {
    std::cout << "failed to create file <" << fileNamePlots << ">\n";
    throw 2;
  }
  // string buffers
  //char ylabel[50];   // y-axis label


  if (DYTools::study2D) {
    TCanvas *c2D = MakeCanvas("canvXsectTh_2D","canvXsectTh_2D",600,600+300*DYTools::study2D);
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
  }

  {
    TCanvas *c1D = MakeCanvas("canvXsectTh_1D","canvXsectTh_1D",600,600);
    CPlot *cp=new CPlot("cplot_1D","","M_{ee} (GeV)","counts");
    TH1F *hXsec= plotXsec1D("xsec1D",nEvents,nEventsErr,cp,kBlue);
    hXsec->SetDirectory(0);
    cp->SetLogx();
    cp->Draw(c1D,false,"png");
    c1D->Update();
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
  /*
  cout << endl;
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  cout << endl; 
  
  cout << labelv[0] << " file: " << fnamev[0] << endl;
  cout << "     Number of generated events: " << nZv << endl;

  if (DYTools::study2D==0)
    {
      printf(" mass bin    preselected      passed     total_A_GEN      BB-BB_A_GEN      EB-BB_A_GEN      EB-EB_A_GEN\n");
      for(int i=0; i<DYTools::nMassBins; i++){
        printf(" %4.0f-%4.0f   %10.1f   %10.1f   %7.4f+-%6.4f  %7.4f+-%6.4f  %7.4f+-%6.4f  %7.4f+-%6.4f \n",
	   DYTools::massBinLimits[i], DYTools::massBinLimits[i+1],
	   nEventsv(i,0), nPassv(i,0),
	   accv(i,0), accErrv(i,0),
	   accBBv(i,0), accErrBBv(i,0),
	   accBEv(i,0), accErrBEv(i,0),
	   accEEv(i,0), accErrEEv(i,0));
      }
    }
  else
    printf("printout format for 2D not chosen");
  cout << endl;

  //sanity check printout
  printSanityCheck(accv, accErrv, "acc");

  */ 
  gBenchmark->Show("getXsec");
}


//=== FUNCTION DEFINITIONS ======================================================================================

//--------------------------------------------------------------------------------------------------
