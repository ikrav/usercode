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

#endif


//=== FUNCTION DECLARATIONS ======================================================================================

void computeNormalizedBinContent(double subset, double subsetErr,
				 double total, double totalErr,
				 double& ratio, double& ratioErr);
void calculateInvertedMatrixErrors(TMatrixD &T, TMatrixD &TErrPos, TMatrixD &TErrNeg,
				   TMatrixD &TinvErr);

//=== MAIN MACRO =================================================================================================

void makeUnfoldingMatrix(const TString input, 
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

  if (debugMode) std::cout << "\n\n\tDEBUG MODE is ON\n\n";
  
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

  TriggerConstantSet constantsSet=DetermineTriggerSet(triggerSetString);  
  assert ( constantsSet != TrigSet_UNDEFINED );

  // For MC the trigger does not depend on run number
  const bool isData=kFALSE;
  TriggerSelection requiredTriggers(constantsSet, isData, 0);

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================

  TRandom random;
  // The random seeds are needed only if we are running this script in systematics mode
  int seed = randomSeed;
  random.SetSeed(seed);
  gRandom->SetSeed(seed);
  // In the case of systematic studies, generate an array of random offsets
  TVectorD shift(escale._nEtaBins); // temporary: this vector is outdated by the new features in the escale obj.class
  shift = 0;
  if(systematicsMode==DYTools::RESOLUTION_STUDY) {
    escale.randomizeSmearingWidth(seed);
    for(int i=0; i<escale._nEtaBins; i++)
      shift[i] = gRandom->Gaus(0,1);
  }

  // prepare tools for ESCALE_RESIDUAL
  /*
  TH1F *shapeWeights=NULL;
  if (systematicsMode==DYTools::ESCALE_RESIDUAL) {
    TString shapeFName=TString("../root_files/yields/") + dirTag + TString("/shape_weights.root");
    std::cout << "Obtaining shape_weights.root from <" << shapeFName << ">\n";
    TFile fshape(shapeFName);
    if (!fshape.IsOpen()) {
      std::cout << "failed to open a file <" << shapeFName << ">\n";
      throw 2;
    }
    shapeWeights = (TH1F*)fshape.Get("weights");
    shapeWeights->SetDirectory(0);
    dirTag += TString("_escale_residual");
    std::cout << "changing dirTag to <" << dirTag << ">\n";
  }
  */

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

  int nUnfoldingBins = DYTools::getTotalNumberOfBins();

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

  // MC spectra for storage in ROOT file
  //int nYBinsMax = DYTools::findMaxYBins(); // commented out - rely on DYTools.hh
  TMatrixD yieldsMcPostFsrGen(DYTools::nMassBins,nYBinsMax);
  TMatrixD yieldsMcPostFsrRec(DYTools::nMassBins,nYBinsMax);
  TMatrixD yieldsMcGen(DYTools::nMassBins,nYBinsMax);
  // The errors 2D arrays are not filled at the moment. It needs
  // to be done carefully since events are weighted.
  // For each bin, the error would be sqrt(sum weights^2).
  TMatrixD yieldsMcPostFsrGenErr(DYTools::nMassBins,nYBinsMax);
  TMatrixD yieldsMcPostFsrRecErr(DYTools::nMassBins,nYBinsMax);

  // Matrices for unfolding
  TMatrixD DetMigration(nUnfoldingBins, nUnfoldingBins);
  TMatrixD DetMigrationErr(nUnfoldingBins, nUnfoldingBins);
  TMatrixD DetResponse(nUnfoldingBins, nUnfoldingBins);
  TMatrixD DetResponseErrPos(nUnfoldingBins, nUnfoldingBins);
  TMatrixD DetResponseErrNeg(nUnfoldingBins, nUnfoldingBins);
  for(int i=0; i<nUnfoldingBins; i++){
    for(int j=0; j<nUnfoldingBins; j++){
      DetMigration   (i,j) = 0;
      DetMigrationErr(i,j) = 0;
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
      if (debugMode && (ientry>100000)) break;

      genBr->GetEntry(ientry);
      infoBr->GetEntry(ientry);

      double reweight;
      if (systematicsMode!=DYTools::FSR_STUDY) reweight=1.0;
      else if (((gen->mass)-(gen->vmass))>massLimit) reweight=1.0;
      else reweight=reweightFsr;

      if (ientry<20) std::cout<<"reweight="<<reweight<<std::endl;

      int iMassBinGen = DYTools::findMassBin(gen->mass);
      int iYBinGen = DYTools::findAbsYBin(iMassBinGen, gen->y);
      if ( (iMassBinGen!=-1) && (iYBinGen!=-1) ) {
	yieldsMcGen(iMassBinGen,iYBinGen) += reweight * scale * gen->weight;
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
        if((fabs(dielectron->scEta_1)>kECAL_GAP_LOW) && (fabs(dielectron->scEta_1)<kECAL_GAP_HIGH)) continue;
        if((fabs(dielectron->scEta_2)>kECAL_GAP_LOW) && (fabs(dielectron->scEta_2)<kECAL_GAP_HIGH)) continue;
	if((fabs(dielectron->scEta_1) > 2.5)       || (fabs(dielectron->scEta_2) > 2.5))       continue;  // outside eta range? Skip to next event...
	
	// Asymmetric SC Et cuts
	if( ! ( ( dielectron->scEt_1 > DYTools::etMinLead  && dielectron->scEt_2 > DYTools::etMinTrail)
		|| ( dielectron->scEt_1 > DYTools::etMinTrail  && dielectron->scEt_2 > DYTools::etMinLead) )) continue;
    	
	// Both electrons must match trigger objects. At least one ordering
	// must match
	if( ! requiredTriggers.matchTwoTriggerObjectsAnyOrder( dielectron->hltMatchBits_1,
							       dielectron->hltMatchBits_2,
							       info->runNum) ) continue;
	
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
	int iMassGenPostFsr = DYTools::findMassBin(gen->mass);
	int iYGenPostFsr = DYTools::findAbsYBin(iMassGenPostFsr, gen->y);
	if( iMassGenPostFsr != -1 && iYGenPostFsr != -1)
	  yieldsMcPostFsrGen(iMassGenPostFsr, iYGenPostFsr) += reweight * scale * gen->weight;
      
	// Fill the matrix of the reconstruction level mass and rapidity
	int iMassReco = DYTools::findMassBin(massResmeared);
	int iYReco = DYTools::findAbsYBin(iMassReco, dielectron->y);
	if( iMassReco != -1 && iYReco != -1)
	  yieldsMcPostFsrRec(iMassReco, iYReco) += reweight * scale * gen->weight;
	
        // Unlike the mass vs Y reference yields matrices, to prepare the
	// migration matrix we flatten (mass,Y) into a 1D array, and then
	// store (mass,Y in 1D)_gen vs (mass,Y in 1D)_rec
	int iIndexFlatGen  = DYTools::findIndexFlat(iMassGenPostFsr, iYGenPostFsr);
 	int iIndexFlatReco = DYTools::findIndexFlat(iMassReco, iYReco);
	if( iIndexFlatReco != -1 && iIndexFlatReco < nUnfoldingBins
	    && iIndexFlatGen != -1 && iIndexFlatGen < nUnfoldingBins ){
	  double fullWeight = reweight * scale * gen->weight;
          DetMigration(iIndexFlatGen,iIndexFlatReco) += fullWeight;
	  // Accumulate sum of weights squared, sqrt of the sum is computed later
          DetMigrationErr(iIndexFlatGen,iIndexFlatReco) += fullWeight*fullWeight;
	}
	
        Bool_t isB1 = (fabs(dielectron->scEta_1)<kECAL_GAP_LOW);
        Bool_t isB2 = (fabs(dielectron->scEta_2)<kECAL_GAP_LOW);         
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
  delete gen;

  // Compute the errors on the elements of migration matrix
  // by simply taking the square root over the accumulated sum(w^2)
  for(int i=0; i < DetMigration.GetNrows(); i++)
    for(int j=0; j < DetMigration.GetNcols(); j++)
      if( DetMigrationErr(i,j) >=0 )
	DetMigrationErr(i,j) = sqrt( DetMigrationErr(i,j) );
      else
	printf("makeUnfoldingMatrix::Error: negative weights in DetMigrationErr\n");
  
  // Find response matrix, which is simply the normalized migration matrix
  std::cout << "find response matrix" << std::endl;
  double tCentral, tErr;
  for(int igen = 0; igen < DetMigration.GetNrows(); igen++){
    // First find the normalization for the given generator level slice
    double nEventsInGenBin = 0;
    double nEventsInGenBinErr = 0;
    for(int ireco = 0; ireco < DetMigration.GetNcols(); ireco++){
      nEventsInGenBin += DetMigration(igen,ireco);
      nEventsInGenBinErr += (DetMigrationErr(igen,ireco)*
				 DetMigrationErr(igen,ireco));
    }
    nEventsInGenBinErr = sqrt(nEventsInGenBinErr);

    // Now normalize each element and find errors
    for(int ireco = 0; ireco < DetMigration.GetNcols(); ireco++){
      tCentral = 0;
      tErr     = 0;
      computeNormalizedBinContent(DetMigration(igen,ireco),
				  DetMigrationErr(igen,ireco),
				  nEventsInGenBin,
				  nEventsInGenBinErr,
				  tCentral, tErr);
      DetResponse      (igen,ireco) = tCentral;
      DetResponseErrPos(igen,ireco) = tErr;
      DetResponseErrNeg(igen,ireco) = tErr;
    }
  }

  std::cout << "find inverted response matrix" << std::endl;

  // Find inverted response matrix
  TMatrixD DetInvertedResponse = DetResponse;
  Double_t det;
  DetInvertedResponse.Invert(&det);
  TMatrixD DetInvertedResponseErr(DetInvertedResponse.GetNrows(), DetInvertedResponse.GetNcols());
  calculateInvertedMatrixErrors(DetResponse, DetResponseErrPos, DetResponseErrNeg, DetInvertedResponseErr);


  TVectorD DetResponseArr(nUnfoldingBins);
  TVectorD DetInvertedResponseArr(nUnfoldingBins), DetInvertedResponseErrArr(nUnfoldingBins);
  TVectorD yieldsMcPostFsrGenArr(nUnfoldingBins), yieldsMcPostFsrRecArr(nUnfoldingBins);

  int resFlatten=
    (unfolding::flattenMatrix(DetResponse, DetResponseArr) == 1) &&
    (unfolding::flattenMatrix(DetInvertedResponse, DetInvertedResponseArr) == 1) &&
    (unfolding::flattenMatrix(DetInvertedResponseErr, DetInvertedResponseErrArr) == 1) &&
    (unfolding::flattenMatrix(yieldsMcPostFsrGen, yieldsMcPostFsrGenArr) == 1) &&
    (unfolding::flattenMatrix(yieldsMcPostFsrRec, yieldsMcPostFsrRecArr) == 1);
  if (!resFlatten) {
    std::cout << "Error : failed to flatten the arrays\n";
    assert(0);
  }

  std::cout << "store constants in a file" << std::endl;

  //
  // Store constants and reference arrays in files
  //
  TString outputDir(TString("../root_files/constants/")+dirTag);
  if((systematicsMode==DYTools::RESOLUTION_STUDY) || (systematicsMode==DYTools::FSR_STUDY))
    outputDir = TString("../root_files/systematics/")+dirTag;
  gSystem->mkdir(outputDir,kTRUE);
  //TString unfoldingConstFileName(outputDir+TString("/unfolding_constants.root"));
  TString unfoldingConstFileName=outputDir+
    TString("/unfolding_constants") + analysisTag 
    + TString(".root");
  if(systematicsMode==DYTools::RESOLUTION_STUDY){
    unfoldingConstFileName = outputDir+TString("/unfolding_constants_seed_");
    unfoldingConstFileName += analysisTag;
    unfoldingConstFileName += seed;
    unfoldingConstFileName += ".root";
  }
  if(systematicsMode==DYTools::FSR_STUDY){
    unfoldingConstFileName = outputDir+TString("/unfolding_constants_reweight_");
    unfoldingConstFileName += analysisTag;
    unfoldingConstFileName += int(100*reweightFsr);
    unfoldingConstFileName += ".root";
  }
  TFile fConst(unfoldingConstFileName, "recreate" );
  DetResponse             .Write("DetResponse");
  DetInvertedResponse     .Write("DetInvertedResponse");
  DetInvertedResponseErr  .Write("DetInvertedResponseErr");
  DetResponseArr          .Write("DetResponseFIArray");
  DetInvertedResponseArr  .Write("DetInvertedResponseFIArray");
  DetInvertedResponseErrArr.Write("DetInvertedResponseErrFIArray");
  unfolding::writeBinningArrays(fConst);
  fConst.Close();

  // Store reference MC arrays in a file
  TString refFileName(outputDir+TString("/yields_MC_unfolding_reference_") + analysisTag + TString(".root"));
  TFile fRef(refFileName, "recreate" );
  yieldsMcPostFsrGen.Write("yieldsMcPostFsrGen");
  yieldsMcPostFsrRec.Write("yieldsMcPostFsrRec");
  yieldsMcPostFsrGenArr.Write("yieldsMcPostFsrGenFIArray");
  yieldsMcPostFsrRecArr.Write("yieldsMcPostFsrRecFIArray");
  unfolding::writeBinningArrays(fRef);
  fRef.Close();


  //--------------------------------------------------------------------------------------------------------------
  // Make plots 
  //==============================================================================================================  

  std::cout << "making plots" << std::endl;

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
  TMatrixD resolutionEffect(DYTools::nMassBins,nYBinsMax);
  resolutionEffect = 0;
  for(int i=0; i < resolutionEffect.GetNrows(); i++){
    for(int j=0; j < resolutionEffect.GetNcols(); j++){
      double ngen = yieldsMcPostFsrGen(i,j);
      double nrec = yieldsMcPostFsrRec(i,j);
      if( ngen != 0 )
	resolutionEffect(i,j) = (ngen-nrec)/ngen;
    }
  }
  PlotMatrixVariousBinning(resolutionEffect, "resolution_effect", "LEGO2", NULL);

  // Plot response and inverted response matrices
  TH2F *hResponse = new TH2F("hResponse","",nUnfoldingBins, -0.5, nUnfoldingBins-0.5,
			     nUnfoldingBins, -0.5, nUnfoldingBins-0.5);
  TH2F *hInvResponse = new TH2F("hInvResponse","",nUnfoldingBins, -0.5, nUnfoldingBins-0.5,
				nUnfoldingBins, -0.5, nUnfoldingBins-0.5);
  for(int i=0; i<DetResponse.GetNrows(); i++){
    for(int j=0; j<DetResponse.GetNcols(); j++){
      hResponse->SetBinContent(i,j, DetResponse(i,j));
      hInvResponse->SetBinContent(i,j, DetInvertedResponse(i,j));
    }
  }
  TCanvas *e1 = MakeCanvas("canvResponse","canvResponse",600,600);
  CPlot plotResponse("response","",
		     "flat index gen",
		     "flat index reco");
  plotResponse.AddHist2D(hResponse,"COLZ");
  plotResponse.Draw(e1);
  SaveCanvas(e1,"hResponse");

  TCanvas *e2 = MakeCanvas("canvInvResponse","canvInvResponse",600,600);
  CPlot plotInvResponse("invResponse","",
		     "flat index gen",
		     "flat index reco");
  plotInvResponse.AddHist2D(hInvResponse,"COLZ");
  plotInvResponse.Draw(e2);
  SaveCanvas(e2,"hInvResponse");
 
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

  //--------------------------------------------------------------------------------------------------------------
  // Summary print out
  //==============================================================================================================
  cout << endl;
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  cout << endl; 
  
  // Printout of all constants, uncomment if needed
  //printf("DetCorrFactor:\n"); DetCorrFactor.Print();
  //printf("DetMigration:\n"); DetMigration.Print();
  //printf("DetResponse:\n"); DetResponse.Print();

  //printf("DetInvertedResponse:\n"); DetInvertedResponse.Print();
  //printf("DetInvertedResponseErr:\n"); DetInvertedResponseErr.Print();
  //printf("DetResponseArr:\n"); DetResponseArr.Print();
  //printf("DetInvertedResponseArr:\n"); DetInvertedResponseArr.Print();
  //printf("DetInvertedResonseErrArr:\n"); DetInvertedResponseErrArr.Print();

//   printf("Detector corr factor numerator:\n");
//   DetCorrFactorNumerator.Print();
  //printf("yieldsMcPostFsrRec:\n");
  //yieldsMcPostFsrRec.Print();

  //printf("yieldsMcPostFsrGen:\n");
  //yieldsMcPostFsrGen.Print();

//   printf("Detector corr factor denominator:\n");
//   DetCorrFactorDenominator.Print();
//   printf("yieldsMcPostFsrRecArr:\n");
//   yieldsMcPostFsrRecArr.Print();

  //printf("yieldsMcGen:\n");
  //yieldsMcGen.Print();

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
