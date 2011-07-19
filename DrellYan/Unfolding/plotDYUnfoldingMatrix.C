#if !defined(__CINT__) || defined(__MAKECINT__)
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

// define classes and constants to read in ntuple
#include "../Include/EWKAnaDefs.hh"
#include "../Include/TEventInfo.hh"
#include "../Include/TGenInfo.hh"
#include "../Include/TDielectron.hh"   

// Helper functions for Electron ID selection
#include "../Include/EleIDCuts.hh"
#endif

//=== FUNCTION DECLARATIONS ======================================================================================

// void BayesEfficiency(double passed, double total,
//                      double& eff, double& errlow, double& errhigh);
void EfficiencyDivide(double passed, double total,
		  double& eff, double& errlow, double& errhigh);
void SimpleDivide(double passed, double total,
		  double& eff, double& errlow, double& errhigh);

void calculateInvertedMatrixErrors(TMatrixD &T, TMatrixD &TErrPos, TMatrixD &TErrNeg,
				   TMatrixD &TinvErr);

double extraSmearingSigma(double eta);

//=== MAIN MACRO =================================================================================================

void plotDYUnfoldingMatrix(const TString input) 
{
  gBenchmark->Start("plotDYUnfoldingMatrix");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================
  
  Bool_t doSave  = false;    // save plots?
  TString format = "png";   // output file format
  
  vector<TString> fnamev;   // file names   
  vector<TString> labelv;   // legend label
  vector<Int_t>   colorv;   // color in plots
  vector<Int_t>   linev;    // line style
  vector<Double_t> xsecv;
  vector<Double_t> lumiv;
  TString          dirTag;

  ifstream ifs;
  ifs.open(input.Data());
  assert(ifs.is_open());
  string line;
  Int_t state=0;
  while(getline(ifs,line)) {
    if(line[0]=='#') continue;
    if(state == 0){
      dirTag = TString(line);
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
  
  const Double_t kGAP_LOW  = 1.4442;
  const Double_t kGAP_HIGH = 1.566;


  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================

  TRandom random;
  //  
  // Set up histograms
  //
  vector<TH1F*> hZMassv;//, hZMass2v, hZPtv, hZPt2v, hZyv, hZPhiv;  
  
  char hname[100];
  for(UInt_t ifile = 0; ifile<fnamev.size(); ifile++) {
    sprintf(hname,"hZMass_%i",ifile); hZMassv.push_back(new TH1F(hname,"",500,0,500)); hZMassv[ifile]->Sumw2();
  }

  TH1F *hMassDiff   = new TH1F("hMassDiff","", 100, -30, 30);
  TH1F *hMassDiffBB = new TH1F("hMassDiffBB","", 100, -30, 30);
  TH1F *hMassDiffEB = new TH1F("hMassDiffEB","", 100, -30, 30);
  TH1F *hMassDiffEE = new TH1F("hMassDiffEE","", 100, -30, 30);

  TH1F *hMassDiffV[DYTools::nMassBins];
  for(int i=0; i<DYTools::nMassBins; i++){
    sprintf(hname,"hMassDiffV_%d",i);
    hMassDiffV[i] = new TH1F(hname,"",100,-50,50);
  }

  // MC spectra for storage in ROOT file
  // The FSR of RECO means that this is the spectrum of generator-level
  // post-FSR mass for events that were actually reconstructed (i.e. includes
  // efficiency and acceptances losses)
  TVectorD yieldsMcFsrOfRec    (DYTools::nMassBins);
  TVectorD yieldsMcFsrOfRecErr (DYTools::nMassBins);
  TVectorD yieldsMcRec         (DYTools::nMassBins);
  TVectorD yieldsMcRecErr      (DYTools::nMassBins);
  yieldsMcFsrOfRec     = 0;
  yieldsMcFsrOfRecErr  = 0;
  yieldsMcRec          = 0;
  yieldsMcRecErr       = 0;
  // The yields at generator level with mass bins defined by post-FSR di-leptons
  TVectorD yieldsMcFsr    (DYTools::nMassBins);
  yieldsMcFsr          = 0;

  // Vectors for bin to bin corrections
  TVectorD DetCorrFactorNumerator  (DYTools::nMassBins);
  TVectorD DetCorrFactorDenominator(DYTools::nMassBins);
  TVectorD DetCorrFactor           (DYTools::nMassBins);
  TVectorD DetCorrFactorErrPos     (DYTools::nMassBins);
  TVectorD DetCorrFactorErrNeg     (DYTools::nMassBins);
  DetCorrFactorNumerator   = 0;
  DetCorrFactorDenominator = 0;
  DetCorrFactor            = 0;
  DetCorrFactorErrPos      = 0;
  DetCorrFactorErrNeg      = 0;
  
  // Matrices for unfolding
  TMatrixD DetMigration(DYTools::nMassBins, DYTools::nMassBins);
  TMatrixD DetResponse(DYTools::nMassBins, DYTools::nMassBins);
  TMatrixD DetResponseErrPos(DYTools::nMassBins, DYTools::nMassBins);
  TMatrixD DetResponseErrNeg(DYTools::nMassBins, DYTools::nMassBins);
  for(int i=0; i<DYTools::nMassBins; i++){
    for(int j=0; j<DYTools::nMassBins; j++){
      DetMigration(i,j) = 0;
      DetResponse(i,j) = 0;
      DetResponseErrPos(i,j) = 0;
      DetResponseErrNeg(i,j) = 0;
    }
  }

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
    lumiv[ifile] = eventTree->GetEntries()/xsecv[ifile];
    double scale = lumiv[0]/lumiv[ifile];
    cout << "       -> sample weight is " << scale << endl;

    // Set branch address to structures that will store the info  
    eventTree->SetBranchAddress("Info",&info);                TBranch *infoBr       = eventTree->GetBranch("Info");
    eventTree->SetBranchAddress("Gen",&gen);                  TBranch *genBr = eventTree->GetBranch("Gen");
    eventTree->SetBranchAddress("Dielectron",&dielectronArr); TBranch *dielectronBr = eventTree->GetBranch("Dielectron");
  
    // loop over events    
    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
      genBr->GetEntry(ientry);
      infoBr->GetEntry(ientry);

      // Use post-FSR generator level mass in unfolding
      int ibinFsr = DYTools::findMassBin(gen->mass);
      if(ibinFsr != -1 && ibinFsr < yieldsMcFsr.GetNoElements()){
	yieldsMcFsr[ibinFsr] += scale * gen->weight;
      }

      // For EPS2011 for both data and MC (starting from Summer11 production)
      // we use an OR of the twi triggers below. Both are unpresecaled.
      UInt_t eventTriggerBit = kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL 
	| kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL;
      UInt_t leadingTriggerObjectBit = kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele1Obj
	| kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele1Obj;
      UInt_t trailingTriggerObjectBit = kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele2Obj
	| kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele2Obj;
      
      if(!(info->triggerBits & eventTriggerBit)) continue;  // no trigger accept? Skip to next event...                                   

      // loop through dielectrons
      dielectronArr->Clear();
      dielectronBr->GetEntry(ientry);    
      for(Int_t i=0; i<dielectronArr->GetEntriesFast(); i++) {
        const mithep::TDielectron *dielectron = (mithep::TDielectron*)((*dielectronArr)[i]);
	
	// Apply selection
	// Eta cuts
        if((fabs(dielectron->scEta_1)>kGAP_LOW) && (fabs(dielectron->scEta_1)<kGAP_HIGH)) continue;
        if((fabs(dielectron->scEta_2)>kGAP_LOW) && (fabs(dielectron->scEta_2)<kGAP_HIGH)) continue;
	if((fabs(dielectron->scEta_1) > 2.5)       || (fabs(dielectron->scEta_2) > 2.5))       continue;  // outside eta range? Skip to next event...
	
	// Asymmetric SC Et cuts
	if( ! ( ( dielectron->scEt_1 > DYTools::etMinLead  && dielectron->scEt_2 > DYTools::etMinTrail)
		|| ( dielectron->scEt_1 > DYTools::etMinTrail  && dielectron->scEt_2 > DYTools::etMinLead) )) continue;
    	
	// Both electrons must match trigger objects. At least one ordering
	// must match
	if( ! ( 
	       (dielectron->hltMatchBits_1 & leadingTriggerObjectBit && 
		dielectron->hltMatchBits_2 & trailingTriggerObjectBit )
	       ||
	       (dielectron->hltMatchBits_1 & trailingTriggerObjectBit && 
		dielectron->hltMatchBits_2 & leadingTriggerObjectBit ) ) ) continue;
	
	// The Smurf electron ID package is the same as used in HWW analysis
	// and contains cuts like VBTF WP80 for pt>20, VBTF WP70 for pt<10
	// with some customization, plus impact parameter cuts dz and dxy
	if(!passSmurf(dielectron)) continue;  

        /******** We have a Z candidate! HURRAY! ********/

	// Apply extra smearing to MC reconstructed dielectron mass
	// to better resemble the data
	double smear1 = extraSmearingSigma(dielectron->scEta_1);
        double smear2 = extraSmearingSigma(dielectron->scEta_2);
        double smearTotal = sqrt(smear1*smear1 + smear2*smear2);
        double massResmeared = dielectron->mass + random.Gaus(0.0,smearTotal);

	hZMassv[ifile]->Fill(massResmeared,scale * gen->weight);

	// Fill structures for response matrix and bin by bin corrections

	if(ibinFsr != -1 && ibinFsr < yieldsMcFsrOfRec.GetNoElements()){
	    yieldsMcFsrOfRec[ibinFsr] += scale * gen->weight;
	    DetCorrFactorDenominator(ibinFsr) += scale * gen->weight;
	}else if(ibinFsr >= yieldsMcFsrOfRec.GetNoElements())
	  cout << "ERROR: binning problem" << endl;

	int ibinRec = DYTools::findMassBin(massResmeared);
	if(ibinRec != -1 && ibinRec < yieldsMcRec.GetNoElements()){
	  yieldsMcRec[ibinRec] += scale * gen->weight;
	  DetCorrFactorNumerator  (ibinRec) += scale * gen->weight;
	}else if(ibinRec >= yieldsMcRec.GetNoElements())
	  cout << "ERROR: binning problem" << endl;
	
	if( ibinRec != -1 && ibinRec < yieldsMcRec.GetNoElements() 
	    && ibinFsr != -1 && ibinFsr < yieldsMcRec.GetNoElements() ){
	  DetMigration(ibinFsr,ibinRec) += scale * gen->weight;
	}
	
        Bool_t isB1 = (fabs(dielectron->scEta_1)<kGAP_LOW);
        Bool_t isB2 = (fabs(dielectron->scEta_2)<kGAP_LOW);         
	hMassDiff->Fill(massResmeared - gen->mass);
	if( isB1 && isB2 )
	  hMassDiffBB->Fill(massResmeared - gen->mass);
	if( (isB1 && !isB2) || (!isB1 && isB2) )
	  hMassDiffEB->Fill(massResmeared - gen->mass);
	if( !isB1 && !isB2 )
	  hMassDiffEE->Fill(massResmeared - gen->mass);
	
	if(ibinFsr != -1){
	  hMassDiffV[ibinFsr]->Fill(massResmeared - gen->mass);
	}
// 	if(ibinRec == -1)
// 	  printf("Missed bin:  M_fsr=%f   M_reco=%f\n", gen->mass, massResmeared);
	
      } // end loop over dielectrons
    } // end loop over events 
    delete infile;
    infile=0, eventTree=0;
  } // end loop over files
  delete gen;
  
  // Find bin by bin corrections and the errors
  double tCentral, tErrNeg, tErrPos;
  for(int i=0; i<DYTools::nMassBins; i++){
    if ( DetCorrFactorDenominator(i) != 0 ){
      // This method does not take into account correlation between numerator
      // and denominator in calculation of errors. This is a flaw to be corrected
      // in the future.
      SimpleDivide( DetCorrFactorNumerator(i), DetCorrFactorDenominator(i), tCentral, tErrNeg, tErrPos);
      DetCorrFactor(i) = tCentral;
      DetCorrFactorErrPos(i) = tErrPos;
      DetCorrFactorErrNeg(i) = tErrNeg;
    } else {
      DetCorrFactor(i) = 0;
      DetCorrFactorErrPos(i) = 0;
      DetCorrFactorErrNeg(i) = 0;
    }
  }
  
  // Find response matrix, which is simply the normalized migration matrix
  for(int ifsr = 0; ifsr < DetMigration.GetNrows(); ifsr++){
    double nEventsInFsrMassBin = 0;
    for(int ireco = 0; ireco < DetMigration.GetNcols(); ireco++)
      nEventsInFsrMassBin += DetMigration(ifsr,ireco);
    // Now normalize each element and find errors
    for(int ireco = 0; ireco < DetMigration.GetNcols(); ireco++){
      tCentral = 0;
      tErrPos  = 0;
      tErrNeg  = 0;
      // Note: the event counts supposedly are dominated with events with weight "1"
      // coming from the primary MC sample, so the error is assumed Poissonian in
      // the call for efficiency-calculating function below.
      if( nEventsInFsrMassBin != 0 ){
	EfficiencyDivide(DetMigration(ifsr,ireco), nEventsInFsrMassBin, tCentral, tErrNeg, tErrPos);
      // BayesEfficiency does not seem to be working in newer ROOT versions, 
      // so it is replaced by simler method
//         BayesEfficiency( DetMigration(ifsr,ireco), nEventsInFsrMassBin, tCentral, tErrNeg, tErrPos);
      }
      DetResponse      (ifsr,ireco) = tCentral;
      DetResponseErrPos(ifsr,ireco) = tErrPos;
      DetResponseErrNeg(ifsr,ireco) = tErrNeg;
    }
  }

  // Find inverted response matrix
  TMatrixD DetInvertedResponse = DetResponse;
  Double_t det;
  DetInvertedResponse.Invert(&det);
  TMatrixD DetInvertedResponseErr(DetInvertedResponse.GetNrows(), DetInvertedResponse.GetNcols());
  calculateInvertedMatrixErrors(DetResponse, DetResponseErrPos, DetResponseErrNeg, DetInvertedResponseErr);

  // Package bin limits into TVectorD for storing in a file
  TVectorD BinLimitsArray(DYTools::nMassBins+1);
  for(int i=0; i<DYTools::nMassBins+1; i++)
    BinLimitsArray(i) = DYTools::massBinLimits[i];

  // Store constants in the file
  TString outputDir(TString("../root_files/constants/")+dirTag);
  gSystem->mkdir(outputDir,kTRUE);
  TString unfoldingConstFileName(outputDir+TString("/unfolding_constants.root"));
  TFile fConst(unfoldingConstFileName, "recreate" );
  DetResponse             .Write("DetResponse");
  DetInvertedResponse     .Write("DetInvertedResponse");
  DetInvertedResponseErr  .Write("DetInvertedResponseErr");
  BinLimitsArray          .Write("BinLimitsArray");
  DetCorrFactor           .Write("DetCorrFactor");
  DetCorrFactorErrPos     .Write("DetCorrFactorErrPos");
  DetCorrFactorErrNeg     .Write("DetCorrFactorErrNeg");
  fConst.Close();
  
  // Store reference MC arrays in a file
  TString refFileName(outputDir+TString("/yields_MC_unfolding_reference.root"));
  TFile fRef(refFileName, "recreate" );
  BinLimitsArray    .Write("BinLimitsArray");
  yieldsMcFsr       .Write("yieldsMcFsr");
  yieldsMcFsrOfRec  .Write("yieldsMcFsrOfRec");
  yieldsMcRec       .Write("yieldsMcRec");
  fRef.Close();


  //--------------------------------------------------------------------------------------------------------------
  // Make plots 
  //==============================================================================================================  
  TCanvas *c = MakeCanvas("c","c",800,600);

  // string buffers
  char ylabel[50];   // y-axis label

  // Z mass
  sprintf(ylabel,"a.u. / %.1f GeV/c^{2}",hZMassv[0]->GetBinWidth(1));
  CPlot plotZMass1("zmass1","","m(ee) [GeV/c^{2}]",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    plotZMass1.AddHist1D(hZMassv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotZMass1.SetLogy();
  plotZMass1.Draw(c,doSave,format);

  // Create plots of how reco mass looks and how unfolded mass should look like
  TVectorD massBinCentral     (DYTools::nMassBins);
  TVectorD massBinHalfWidth   (DYTools::nMassBins);
  TVectorD yieldMcFsrOfRecErr (DYTools::nMassBins);
  TVectorD yieldMcRecErr      (DYTools::nMassBins);
  for(int i=0; i < DYTools::nMassBins; i++){
    massBinCentral  (i) = (DYTools::massBinLimits[i+1] + DYTools::massBinLimits[i])/2.0;
    massBinHalfWidth(i) = (DYTools::massBinLimits[i+1] - DYTools::massBinLimits[i])/2.0;
    yieldsMcFsrOfRecErr(i) = sqrt(yieldsMcFsrOfRec[i]);
    yieldsMcRecErr     (i) = sqrt(yieldsMcRec[i]);
  }
  TGraphErrors *grFsrOfRec = new TGraphErrors(massBinCentral, yieldsMcFsrOfRec, 
					      massBinHalfWidth, yieldsMcFsrOfRecErr);
  TGraphErrors *grRec      = new TGraphErrors(massBinCentral, yieldsMcRec, 
					      massBinHalfWidth, yieldsMcRecErr);
  TCanvas *d = MakeCanvas("d","d",800,600);
  CPlot plotZMass2("zmass2","","m(ee) [GeV/c^{2}]","events");
  plotZMass2.AddGraph(grFsrOfRec,"no detector effects","PE",kBlack);
  plotZMass2.AddGraph(grRec,"reconstructed","PE",kBlue);
  plotZMass2.Draw(d);

//   double xsize = 600;
//   double ysize = 600;
  double xsize = 600;
  double ysize = 400;

  // Create the plot of the response matrix
  TH2F *hResponse = new TH2F("hResponse","",DYTools::nMassBins, DYTools::massBinLimits,
			     DYTools::nMassBins, DYTools::massBinLimits);
  for(int i=1; i<=DYTools::nMassBins; i++)
    for(int j=1; j<=DYTools::nMassBins; j++)
      hResponse->SetBinContent(i,j,DetResponse(i-1,j-1));
  TCanvas *e = MakeCanvas("e","e",xsize,ysize);
  CPlot plotResponse("response","",
		     "generator post-FSR m(ee) [GeV/c^{2}]",
		     "reconstructed  m(ee) [GeV/c^{2}]");
  plotResponse.AddHist2D(hResponse,"COLZ");
  e->cd();
  plotResponse.SetLogx();
  plotResponse.SetLogy();
  gPad->SetRightMargin(0.15);
  gStyle->SetPalette(1);
  hResponse->GetXaxis()->SetMoreLogLabels();
  hResponse->GetXaxis()->SetNoExponent();
  hResponse->GetYaxis()->SetNoExponent();
  plotResponse.Draw(e);

  // Create the plot of the inverted response matrix
  TH2F *hInvResponse = new TH2F("hInvResponse","",DYTools::nMassBins, DYTools::massBinLimits,
			     DYTools::nMassBins, DYTools::massBinLimits);
  for(int i=1; i<=DYTools::nMassBins; i++)
    for(int j=1; j<=DYTools::nMassBins; j++)
      hInvResponse->SetBinContent(i,j,DetInvertedResponse(i-1,j-1));
  TCanvas *f = MakeCanvas("f","f",xsize,ysize);
  CPlot plotInvResponse("inverted response","",
			"reconstructed  m(ee) [GeV/c^{2}]",
			"generator post-FSR m(ee) [GeV/c^{2}]");
  plotInvResponse.AddHist2D(hInvResponse,"COLZ");
  f->cd();
  plotInvResponse.SetLogx();
  plotInvResponse.SetLogy();
  gPad->SetRightMargin(0.15);
  gStyle->SetPalette(1);
  hInvResponse->GetXaxis()->SetMoreLogLabels();
  hInvResponse->GetXaxis()->SetNoExponent();
  hInvResponse->GetYaxis()->SetNoExponent();
  plotInvResponse.Draw(f);

  // Create a plot of detector resolution without mass binning
  TCanvas *g = MakeCanvas("g","g",600,600);
  CPlot plotMassDiff("massDiff","","reco mass - gen post-FSR mass [GeV/c^{2}]","a.u.");
  hMassDiffBB->Scale(1.0/hMassDiffBB->GetSumOfWeights());
  hMassDiffEB->Scale(1.0/hMassDiffEB->GetSumOfWeights());
  hMassDiffEE->Scale(1.0/hMassDiffEE->GetSumOfWeights());
  plotMassDiff.AddHist1D(hMassDiffBB,"EB-EB","hist",kBlack);
  plotMassDiff.AddHist1D(hMassDiffEB,"EE-EB","hist",kBlue);
  plotMassDiff.AddHist1D(hMassDiffEE,"EE-EE","hist",kRed);
  plotMassDiff.Draw(g);

  // Create a plot of reco - gen post-FSR mass difference for several mass bins
  TCanvas *h = MakeCanvas("h","h",600,600);
  CPlot plotMassDiffV("massDiffV","","reco mass - gen post-FSR mass [GeV/c^{2}]","a.u.");
  hMassDiffV[7]->Scale(1.0/hMassDiffV[7]->GetSumOfWeights());
  hMassDiffV[6]->Scale(1.0/hMassDiffV[6]->GetSumOfWeights());
  hMassDiffV[5]->Scale(1.0/hMassDiffV[5]->GetSumOfWeights());
  hMassDiffV[4]->Scale(1.0/hMassDiffV[4]->GetSumOfWeights());
  plotMassDiffV.AddHist1D(hMassDiffV[4],"50-60 GeV/c^{2}","hist",kBlack);
  plotMassDiffV.AddHist1D(hMassDiffV[5],"60-76 GeV/c^{2}","hist",kBlue);
  plotMassDiffV.AddHist1D(hMassDiffV[6],"76-86 GeV/c^{2}","hist",kRed);
  plotMassDiffV.AddHist1D(hMassDiffV[7],"86-96 GeV/c^{2}","hist",kMagenta);
  plotMassDiffV.SetXRange(-20,50);
  plotMassDiffV.Draw(h);

  // Create graph of bin-to-bin corrections
  TGraphAsymmErrors *grCorrFactor 
    = new TGraphAsymmErrors(massBinCentral, DetCorrFactor,
			    massBinHalfWidth, massBinHalfWidth,
			    DetCorrFactorErrNeg, DetCorrFactorErrPos);
  TCanvas *c11 = MakeCanvas("c11","c11",800,600);
  CPlot plotCorrFactor("corrFactor","","m(ee) [GeV/c^{2}]","correction factor");
  plotCorrFactor.AddGraph(grCorrFactor,"PE",kBlue);
  plotCorrFactor.SetLogx();
  plotCorrFactor.SetYRange(0,2.0);
  plotCorrFactor.Draw(c11);

  //--------------------------------------------------------------------------------------------------------------
  // Summary print out
  //==============================================================================================================
  cout << endl;
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  cout << endl; 
  
  // Printout of all constants, uncomment if needed
//   DetCorrFactor.Print();
//   DetResponse.Print();
  
//   printf("Detector corr factor numerator:\n");
//   DetCorrFactorNumerator.Print();
//   printf("yieldsMcRec:\n");
//   yieldsMcRec.Print();

//   printf("Detector corr factor denominator:\n");
//   DetCorrFactorDenominator.Print();
//   printf("yieldsMcFsrOfRec:\n");
//   yieldsMcFsrOfRec.Print();

  gBenchmark->Show("plotDYUnfoldingMatrix");
}


//=== FUNCTION DEFINITIONS ======================================================================================

// //--------------------------------------------------------------------------------------------------
// void BayesEfficiency(double passed, double total,
//                      double& eff, double& errlow, double& errhigh){

//   // WARNING: this way of calculation of bayesian efficiency errors
//   // is not correct in the case if "passed" and "total" are not plain event
//   // counts, but weighted event counts. Still, hopefully it is not too badly off.

//   // WARNING: in newer ROOT versions this does not quite work for some reason

//   // This function finds Bayesian efficiency using methods from ROOT.
//   // There is a trick since the actualy Efficiency function in ROOT
//   // is not declared public.

//   if(total == 0) {
//     printf("NOT GOOD! can't calculate efficiency with zero denominator\n");
//     eff=0;
//     errlow=0;
//     errhigh=0;
//     return;
//   }

//   // We simply define two 1-bin histograms, put in the total and passed
//   // counts, and use the available method from root to find 
//   // TGraphAsymmErrors for the efficiency, also with 1 bin.
//   TH1F *hTmpTotal = new TH1F("hTmpTotal","",1,0,1);
//   TH1F *hTmpPassed = new TH1F("hTmpPassed","",1,0,1);
  
// //   for(int i=0; i<total; i++)
//   hTmpTotal->Fill(0.5,total);
// //   for(int i=0; i<passed; i++)
//   hTmpPassed->Fill(0.5,passed);

//   TGraphAsymmErrors * hTmpEff = new TGraphAsymmErrors;
//   if( total >= passed){
//     // The functions have changed at certain point, for later ROOT
//     // use plane Divide with option "b"
//     hTmpEff->BayesDivide(hTmpPassed,hTmpTotal);
// //     hTmpEff->Divide(hTmpPassed,hTmpTotal,"b");
  
//     eff = (hTmpEff->GetY())[0];
//     errlow = hTmpEff->GetErrorYlow(0);
//     errhigh = hTmpEff->GetErrorYhigh(0);
//   }else{
//     printf("Bayes efficiency warning: passed > total, using plain calculation\n");
//     SimpleDivide(passed,total,eff,errlow,errhigh);
//   }
  
//   delete hTmpTotal;
//   delete hTmpPassed;
//   delete hTmpEff;
  
//   return;
// }

void EfficiencyDivide(double passed, double total, 
		      double& eff, double& errlow, double& errhigh){
  
  if(total == 0) {
    printf("NOT GOOD! can't divide with zero denominator\n");
    eff=0;
    errlow=0;
    errhigh=0;
    return;
  }
  
  eff = passed/total;
  if(passed<total){
    if(passed != 0)
      errlow = sqrt(eff*(1-eff)/total);
    else
      errlow = 0;
    errhigh = errlow;
  }else{
    printf("Bayes efficiency warning: passed > total, using plain calculation\n");
    SimpleDivide(passed,total,eff,errlow,errhigh);
  }

  return;
}

void SimpleDivide(double passed, double total, 
                     double& eff, double& errlow, double& errhigh){

  if(total == 0) {
    printf("NOT GOOD! can't divide with zero denominator\n");
    eff=0;
    errlow=0;
    errhigh=0;
    return;
  }

  eff = passed/total;
  if(passed != 0)
    errlow = eff*sqrt( 1/passed + 1/total );
  else
    errlow = 0;
  errhigh = errlow;
  return;
}



void calculateInvertedMatrixErrors(TMatrixD &T, TMatrixD &TErrPos, TMatrixD &TErrNeg,
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
  int N = 100000;
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

double extraSmearingSigma(double eta){

  double smear = 0.0;

  // Energy scale corrections from Duncan Ralph derived July 2011
  // based on 204 pb-1 May ReReco and ~650 pb-1 of prompt reco
  // as well as 41X powheg MC  
  const int nEtaBins = 10;
  double smearEtaBinLimits[nEtaBins+1] = 
    {-2.50001, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.50001};
  double smearValues[nEtaBins] = 
    {1.99, 1.36, 0.80, 0.63, 0.55, 0.59, 0.62, 1.27, 1.28, 2.00};
  double smearErrors[nEtaBins] =
    {0.08, 0.14, 0.07, 0.05, 0.09, 0.05, 0.03, 0.09, 0.07, 0.15};

  double feta = fabs(eta);
  for(int i=0; i<nEtaBins; i++){
    if(feta >= smearEtaBinLimits[i] && feta < smearEtaBinLimits[i+1] ){
      smear = smearValues[i];
      break;
    }
  }

  return smear;
  
}


