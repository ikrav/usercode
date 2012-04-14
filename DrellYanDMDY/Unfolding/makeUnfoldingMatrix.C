#if !defined(__CINT__) || defined(__MAKECINT__)

// TODO:
//
// Are trigger bits implemented correctly for full 2011 dataset?
//
// Sort out EfficiencyDivide and other functions!
//
// Review systematics for unfolding. Can it be done more
// efficiently?
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

// define classes and constants to read in ntuple
#include "../Include/EWKAnaDefs.hh"
#include "../Include/TEventInfo.hh"
#include "../Include/TGenInfo.hh"
#include "../Include/TDielectron.hh"   

// Helper functions for Electron ID selection
#include "../Include/EleIDCuts.hh"

#include "../Include/ElectronEnergyScale.hh" //extra smearing


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
//=== MAIN MACRO =================================================================================================

void makeUnfoldingMatrix(const TString input, int systematicsMode = DYTools::NORMAL, int randomSeed = 1, double reweightFsr = 1.0, double massLimit = -1.0, int debugMode=0)
//systematicsMode 0 (NORMAL) - no systematic calc
//1 (RESOLUTION_STUDY) - systematic due to smearing, 2 (FSR_STUDY) - systematics due to FSR, reweighting
//check mass spectra with reweightFsr = 0.95; 1.00; 1.05  
//mass value until which do reweighting
{
  gBenchmark->Start("makeUnfoldingMatrix");

  if (systematicsMode==DYTools::NORMAL)
    std::cout<<"Running script in the NORMAL mode"<<std::endl;
  else if (systematicsMode==DYTools::RESOLUTION_STUDY)
    std::cout<<"Running script in the RESOLUTION_STUDY mode"<<std::endl;
  else if (systematicsMode==DYTools::FSR_STUDY)
    std::cout<<"Running script in the FSR_STUDY mode"<<std::endl;
  else { 
    std::cout<<"requested mode not recognized"<<std::endl;
    assert(0);
  }

  if (debugMode) std::cout << "\n\n\tDEBUG MODE is ON\n\n";
  
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
  TString          escaleTag; // Energy scale calibrations tag

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
  
  // 
  // Set up energy scale corrections
  //
  ElectronEnergyScale escale(escaleTag);
  escale.print();

  if( !escale.isInitialized()) {
    printf("Failed to match escale calibration. Tag: >>%s<<\n", escaleTag.Data());
    assert(0);
  }

  const Double_t kGAP_LOW  = 1.4442;
  const Double_t kGAP_HIGH = 1.566;


  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================

  TRandom random;
  // The random seeds are needed only if we are running this script in systematics mode
  int seed = randomSeed;
  random.SetSeed(seed);
  gRandom->SetSeed(seed);
//   // In the case of systematic studies, generate an array of random offsets
//   TVectorD shift(escale::nEtaBins);
//   shift = 0;
//   if(systematicsMode==DYTools::RESOLUTION_STUDY)
//     for(int i=0; i<escale::nEtaBins; i++)
//       shift[i] = gRandom->Gaus(0,1);

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

  int nUnfoldingBins = DYTools::getNumberOf2DBins();

  TH1F *hMassDiffV[nUnfoldingBins];
  for(int i=0; i<nUnfoldingBins; i++){
    sprintf(hname,"hMassDiffV_%d",i);
    hMassDiffV[i] = new TH1F(hname,"",100,-50,50);
  }

  // MC spectra for storage in ROOT file
  // The FSR of RECO means that this is the spectrum of generator-level
  // post-FSR mass for events that were actually reconstructed (i.e. includes
  // efficiency and acceptances losses)
  TVectorD yieldsMcFsrOfRec    (nUnfoldingBins);
  TVectorD yieldsMcFsrOfRecErr (nUnfoldingBins);
  TVectorD yieldsMcRec         (nUnfoldingBins);
  TVectorD yieldsMcRecErr      (nUnfoldingBins);
  yieldsMcFsrOfRec     = 0;
  yieldsMcFsrOfRecErr  = 0;
  yieldsMcRec          = 0;
  yieldsMcRecErr       = 0;
  // The yields at generator level with mass bins defined by post-FSR di-leptons
  TVectorD yieldsMcFsr    (nUnfoldingBins);
  yieldsMcFsr          = 0;

  // Vectors for bin to bin corrections
  TVectorD DetCorrFactorNumerator  (nUnfoldingBins);
  TVectorD DetCorrFactorDenominator(nUnfoldingBins);
  TVectorD DetCorrFactor           (nUnfoldingBins);
  TVectorD DetCorrFactorErrPos     (nUnfoldingBins);
  TVectorD DetCorrFactorErrNeg     (nUnfoldingBins);
  DetCorrFactorNumerator   = 0;
  DetCorrFactorDenominator = 0;
  DetCorrFactor            = 0;
  DetCorrFactorErrPos      = 0;
  DetCorrFactorErrNeg      = 0;
  
  // Matrices for unfolding
  TMatrixD DetMigration(nUnfoldingBins, nUnfoldingBins);
  TMatrixD DetResponse(nUnfoldingBins, nUnfoldingBins);
  TMatrixD DetResponseErrPos(nUnfoldingBins, nUnfoldingBins);
  TMatrixD DetResponseErrNeg(nUnfoldingBins, nUnfoldingBins);
  for(int i=0; i<nUnfoldingBins; i++){
    for(int j=0; j<nUnfoldingBins; j++){
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
      if (debugMode && (ientry>100000)) break;

      genBr->GetEntry(ientry);
      infoBr->GetEntry(ientry);

      double reweight;
      if (systematicsMode!=DYTools::FSR_STUDY) reweight=1.0;
      else if (((gen->mass)-(gen->vmass))>massLimit) reweight=1.0;
      else reweight=reweightFsr;

      if (ientry<20) std::cout<<"reweight="<<reweight<<std::endl;


      // Use post-FSR generator level mass and rapidity in unfolding
      int iMassGenPostFsr = DYTools::findMassBin2D(gen->mass);
      int iYGenPostFsr = DYTools::findAbsYBin2D(iMassGenPostFsr, gen->y);
      int iIndexFlatGenPostFsr = DYTools::findIndexFlat(iMassGenPostFsr, iYGenPostFsr);
      if(iIndexFlatGenPostFsr != -1 && iIndexFlatGenPostFsr < yieldsMcFsr.GetNoElements()){
         yieldsMcFsr[iIndexFlatGenPostFsr] += reweight * scale * gen->weight;
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

        // We have a Z candidate! HURRAY! 

	// 
	// THIS NEEDS TO BE REVISED, SYSTEMATICS FOR NOW IS COMMENTED OUT
	//
// 	// Apply extra smearing to MC reconstructed dielectron mass
// 	// to better resemble the data
// 	double smear1 = escale::extraSmearingSigma(dielectron->scEta_1);
//         double smear2 = escale::extraSmearingSigma(dielectron->scEta_2);
// 	// In systematics mode, overwrite the smear values with
// 	// shifted ones.
// 	if(systematicsMode==DYTools::RESOLUTION_STUDY){
// 	  smear1 = escale::extraSmearingSigmaShifted(dielectron->scEta_1,shift);
// 	  smear2 = escale::extraSmearingSigmaShifted(dielectron->scEta_2,shift);
// 	}
//         double smearTotal = sqrt(smear1*smear1 + smear2*smear2);
//         double massResmeared = dielectron->mass + random.Gaus(0.0,smearTotal);
	double smearingCorrection = escale.generateMCSmear(dielectron->scEta_1, 
							   dielectron->scEta_2);
	double massResmeared = dielectron->mass + smearingCorrection;

	hZMassv[ifile]->Fill(massResmeared,scale * gen->weight);

	// Fill structures for response matrix and bin by bin corrections

	if(iIndexFlatGenPostFsr != -1 && iIndexFlatGenPostFsr < yieldsMcFsrOfRec.GetNoElements()){
	    yieldsMcFsrOfRec[iIndexFlatGenPostFsr] += reweight * scale * gen->weight;
	    DetCorrFactorDenominator(iIndexFlatGenPostFsr) += reweight * scale * gen->weight;
          }

	else if(iIndexFlatGenPostFsr >= yieldsMcFsrOfRec.GetNoElements())
	  cout << "ERROR: binning problem" << endl;

	// Find bin of the reconstructed mass and rapidity
	int iMassReco = DYTools::findMassBin2D(massResmeared);
	int iYReco = DYTools::findAbsYBin2D(iMassReco, dielectron->y);
	int iIndexFlatReco = DYTools::findIndexFlat(iMassReco, iYReco);
	if(iIndexFlatReco != -1 && iIndexFlatReco < yieldsMcRec.GetNoElements()){
	    yieldsMcRec[iIndexFlatReco] += reweight * scale * gen->weight;
	    DetCorrFactorNumerator(iIndexFlatReco) += reweight * scale * gen->weight;
          }

	else if(iIndexFlatReco >= yieldsMcRec.GetNoElements())
	  cout << "ERROR: binning problem" << endl;
	
	if( iIndexFlatReco != -1 && iIndexFlatReco < yieldsMcRec.GetNoElements() 
	    && iIndexFlatGenPostFsr != -1 && iIndexFlatGenPostFsr < yieldsMcRec.GetNoElements() ){
          DetMigration(iIndexFlatGenPostFsr,iIndexFlatReco) += reweight * scale * gen->weight;
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
	
	if(iIndexFlatGenPostFsr != -1){
	  hMassDiffV[iIndexFlatGenPostFsr]->Fill(massResmeared - gen->mass);
	}
// 	if(iIndexFlatReco == -1)
// 	  printf("Missed bin:  M_fsr=%f   M_reco=%f\n", gen->mass, massResmeared);
	
      } // end loop over dielectrons

    } // end loop over events 
    delete infile;
    infile=0, eventTree=0;
  } // end loop over files
  delete gen;


  //finish reweighting of mass spectra

  
  // Find bin by bin corrections and the errors
  double tCentral, tErrNeg, tErrPos;
  for(int i=0; i<nUnfoldingBins; i++){
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
  //Double_t det;
//   DetInvertedResponse.Invert(&det);
  TMatrixD DetInvertedResponseErr(DetInvertedResponse.GetNrows(), DetInvertedResponse.GetNcols());
//   calculateInvertedMatrixErrors(DetResponse, DetResponseErrPos, DetResponseErrNeg, DetInvertedResponseErr);

  // Package bin limits into TVectorD for storing in a file
//   TVectorD BinLimitsArray(nUnfoldingBins+1);
//   for(int i=0; i<nUnfoldingBins+1; i++)
//     BinLimitsArray(i) = DYTools::massBinLimits[i];

  // Store constants in the file
  TString outputDir(TString("../root_files/constants/")+dirTag);
  if((systematicsMode==DYTools::RESOLUTION_STUDY) || (systematicsMode==DYTools::FSR_STUDY))
    outputDir = TString("../root_files/systematics/")+dirTag;
  gSystem->mkdir(outputDir,kTRUE);
  TString unfoldingConstFileName(outputDir+TString("/unfolding_constants.root"));
  if(systematicsMode==DYTools::RESOLUTION_STUDY){
    unfoldingConstFileName = outputDir+TString("/unfolding_constants_seed_");
    unfoldingConstFileName += seed;
    unfoldingConstFileName += ".root";
  }
  if(systematicsMode==DYTools::FSR_STUDY){
    unfoldingConstFileName = outputDir+TString("/unfolding_constants_reweight_");
    unfoldingConstFileName += int(100*reweightFsr);
    unfoldingConstFileName += ".root";
  }
  TFile fConst(unfoldingConstFileName, "recreate" );
  DetResponse             .Write("DetResponse");
  DetInvertedResponse     .Write("DetInvertedResponse");
  DetInvertedResponseErr  .Write("DetInvertedResponseErr");
//   BinLimitsArray          .Write("BinLimitsArray");
  DetCorrFactor           .Write("DetCorrFactor");
  DetCorrFactorErrPos     .Write("DetCorrFactorErrPos");
  DetCorrFactorErrNeg     .Write("DetCorrFactorErrNeg");
  fConst.Close();
  
  // Store reference MC arrays in a file
  TString refFileName(outputDir+TString("/yields_MC_unfolding_reference.root"));
  TFile fRef(refFileName, "recreate" );
//   BinLimitsArray    .Write("BinLimitsArray");
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
  TVectorD massBinCentral     (nUnfoldingBins);
  TVectorD massBinHalfWidth   (nUnfoldingBins);
  TVectorD yieldMcFsrOfRecErr (nUnfoldingBins);
  TVectorD yieldMcRecErr      (nUnfoldingBins);
  for(int i=0; i < nUnfoldingBins; i++){
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
  TH2F *hResponse = new TH2F("hResponse","",nUnfoldingBins, DYTools::massBinLimits,
			     nUnfoldingBins, DYTools::massBinLimits);
  for(int i=1; i<=nUnfoldingBins; i++)
    for(int j=1; j<=nUnfoldingBins; j++)
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
  TH2F *hInvResponse = new TH2F("hInvResponse","",nUnfoldingBins, DYTools::massBinLimits,
			     nUnfoldingBins, DYTools::massBinLimits);
  for(int i=1; i<=nUnfoldingBins; i++)
    for(int j=1; j<=nUnfoldingBins; j++)
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

  gBenchmark->Show("makeUnfoldingMatrix");
}


//=== FUNCTION DEFINITIONS ======================================================================================
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
