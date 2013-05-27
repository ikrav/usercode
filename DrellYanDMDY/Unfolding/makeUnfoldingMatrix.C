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
#include "../Include/InputFileMgr.hh"

//for getting matrix condition number
#include <TDecompLU.h>

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
  TMatrixD yieldsMcPostFsrGen(DYTools::nMassBins,DYTools::nYBinsMax);
  TMatrixD yieldsMcPostFsrRec(DYTools::nMassBins,DYTools::nYBinsMax);
  TMatrixD yieldsMcGen(DYTools::nMassBins,DYTools::nYBinsMax);       // to compare with DrellYan1D
  // The errors 2D arrays are not filled at the moment. It needs
  // to be done carefully since events are weighted.
  // For each bin, the error would be sqrt(sum weights^2).
  TMatrixD yieldsMcPostFsrGenErr(DYTools::nMassBins,DYTools::nYBinsMax);
  TMatrixD yieldsMcPostFsrRecErr(DYTools::nMassBins,DYTools::nYBinsMax);

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
      if (debugMode && (ientry>10)) break;

      genBr->GetEntry(ientry);
      infoBr->GetEntry(ientry);

      double reweight;
      if (systematicsMode!=DYTools::FSR_STUDY) reweight=1.0;
      else if (((gen->mass)-(gen->vmass))>massLimit) reweight=1.0;
      else reweight=reweightFsr;

      if (ientry<20) {
	printf("reweight=%4.2lf, dE_fsr=%+6.4lf\n",reweight,(gen->mass-gen->vmass));
      }

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
	// Et and eta cuts
	if( ! DYTools::goodEtEtaPair( dielectron->scEt_1, dielectron->scEta_1,
				      dielectron->scEt_2, dielectron->scEta_2 ) ) continue;

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
	if( DYTools::energy8TeV == 1){
	  if(!passEGMID2012(dielectron, WP_MEDIUM, info->rhoLowEta)) continue;  
	}else{
	  if(!passEGMID2011(dielectron, WP_MEDIUM, info->rhoLowEta)) continue;  
	}

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
	int iMassGenPostFsr = DYTools::findMassBin(gen->mass);
	int iYGenPostFsr = DYTools::findAbsYBin(iMassGenPostFsr, gen->y);
	if( iMassGenPostFsr != -1 && iYGenPostFsr != -1)
	  yieldsMcPostFsrGen(iMassGenPostFsr, iYGenPostFsr) += reweight * scale * gen->weight;
      
	// Fill the matrix of the reconstruction level mass and rapidity
	int iMassReco = DYTools::findMassBin(massResmeared);
	int iYReco = DYTools::findAbsYBin(iMassReco, dielectron->y);
	double shape_weight = 1.0;
	if( iMassReco != -1 && iYReco != -1) {
	  if (shapeWeights) {
	    shape_weight = (*shapeWeights)[iMassReco][iYReco];
	    //std::cout << "massResmeared=" << massResmeared << ", iMassReco=" << iMassReco << ", shapeWeight=" << shape_weight << "\n";
	  }
	  yieldsMcPostFsrRec(iMassReco, iYReco) += reweight * scale * gen->weight;
	}
	
        // Unlike the mass vs Y reference yields matrices, to prepare the
	// migration matrix we flatten (mass,Y) into a 1D array, and then
	// store (mass,Y in 1D)_gen vs (mass,Y in 1D)_rec
	int iIndexFlatGen  = DYTools::findIndexFlat(iMassGenPostFsr, iYGenPostFsr);
 	int iIndexFlatReco = DYTools::findIndexFlat(iMassReco, iYReco);
	if( iIndexFlatReco != -1 && iIndexFlatReco < nUnfoldingBins
	    && iIndexFlatGen != -1 && iIndexFlatGen < nUnfoldingBins ){
	  double fullWeight = reweight * scale * gen->weight * shape_weight;
	  //std::cout << "adding DetMig(" << iIndexFlatGen << "," << iIndexFlatReco << ") = " << reweight << "*" << scale << "*" << gen->weight << "*" << shape_weight << " = "  << (reweight * scale * gen->weight * shape_weight) << "\n";
          DetMigration(iIndexFlatGen,iIndexFlatReco) += fullWeight;
	  // Accumulate sum of weights squared, sqrt of the sum is computed later
          DetMigrationErr(iIndexFlatGen,iIndexFlatReco) += fullWeight*fullWeight;
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
  delete gen;

  //return;

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

  //Calculation of Unfolding matrix errors using different method
  TMatrixD DetInvertedResponseErr2(nUnfoldingBins,nUnfoldingBins);
  DetInvertedResponseErr2=DetInvertedResponse;
  DetInvertedResponseErr2*=DetInvertedResponse;
  DetInvertedResponseErr2*=DetResponseErrNeg;
  for (int i=0; i<nUnfoldingBins; i++)
    for (int j=0; j<nUnfoldingBins; j++)
      {
         if (DetInvertedResponseErr2(i,j)<0) DetInvertedResponseErr2(i,j)=-DetInvertedResponseErr2(i,j);
      }

  std::cout << "store constants in a file" << std::endl;

  //
  // Store constants and reference arrays in files
  //
  TString outputDir(TString("../root_files/constants/")+dirTag);
  if((systematicsMode==DYTools::RESOLUTION_STUDY) || (systematicsMode==DYTools::FSR_STUDY))
    outputDir = TString("../root_files/systematics/")+dirTag;
  gSystem->mkdir(outputDir,kTRUE);

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

  //TString unfoldingConstFileName(outputDir+TString("/unfolding_constants.root"));
  TString unfoldingConstFileName=outputDir+
    TString("/unfolding_constants") + fnameTag + TString(".root");
  std::cout << "unfoldingConstFileName=<" << unfoldingConstFileName << ">\n";

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
  TString refFileName(outputDir+TString("/yields_MC_unfolding_reference_") + DYTools::analysisTag + TString(".root"));
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
  TMatrixD resolutionEffect(DYTools::nMassBins,DYTools::nYBinsMax);
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

  //draw errors of Unfolding matrix
  TCanvas *cErrors = new TCanvas("cErrors","DetInvertedResponseErr");
  cErrors->Divide(2,2);
  cErrors->cd(1);
  DetInvertedResponseErr.Draw("LEGO2");
  cErrors->cd(2);
  DetInvertedResponseErr2.Draw("LEGO2");


  

  //--------------------------------------------------------------------------------------------------------------
  // Summary print out
  //==============================================================================================================
  cout << endl;
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  cout << endl; 

  //matrix condition number
  TDecompLU lu(DetResponse);
  double condLU=lu.Condition();
  std::cout << " condition number from TDecompLU condLU= " << condLU << std::endl;
  std::cout << " condition number ||DetResponse||*||DetResponseInv||=" << DetResponse.Norm1()*DetInvertedResponse.Norm1() << std::endl;
  std::cout << " chk ROOT bug: -condLU*||DetResponse||=" << (-condLU*DetResponse.Norm1()) << "\n" << std::endl;

  //Print errors of the Unfolding matrix when they exceed 0.1
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
