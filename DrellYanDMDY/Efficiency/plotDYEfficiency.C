#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TClonesArray.h>           // ROOT array class
#include <TCanvas.h>                // class for drawing
#include <TH1F.h>                   // 1D histograms
#include <TCanvas.h>                // class for drawing
#include <TH1F.h>                   // 1D histograms
#include <TGraphErrors.h>           // graphs
#include <TProfile.h>               // profile histograms
#include <TBenchmark.h>             // class to track macro running statistics
#include <TLorentzVector.h>         // class for Lorentz vector computations
#include <TVectorD.h>
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
#include "../Include/latexPrintouts.hh"

// define classes and constants to read in ntuple
#include "../Include/EWKAnaDefs.hh"
#include "../Include/TEventInfo.hh"
#include "../Include/TGenInfo.hh"
#include "../Include/TDielectron.hh"   
#include "../Include/TriggerSelection.hh"
#include "../Include/TVertex.hh"

// Helper functions for Electron ID selection
#include "../Include/EleIDCuts.hh"

#include "../Include/FEWZ.hh"
#include "../Include/EventSelector.hh"
#include "../Include/PUReweight.hh"
#include "../Include/InputFileMgr.hh"



#endif

#define usePUReweight  // Whether apply PU reweighting


//=== FUNCTION DECLARATIONS ======================================================================================

//=== MAIN MACRO =================================================================================================

void plotDYEfficiency(const TString input, 
		      const TString triggerSetString="Full2011DatasetTriggers",
		      int debugMode=0) 
{
  gBenchmark->Start("plotDYEfficiency");

  if (debugMode) std::cout << "\n\n\tDEBUG MODE is ON\n\n";

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================
  
  //Bool_t doSave  = false;    // save plots?
  TString format = "png";   // output file format
  
  vector<TString> fnamev;   // file names   
  vector<TString> labelv;   // legend label
  vector<Int_t>   colorv;   // color in plots
  vector<Int_t>   linev;    // line style
  vector<Double_t> xsecv;
  vector<Double_t> lumiv;
  TString          dirTag;

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
    //escaleTag=mcInp.escaleTag();
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
  
  //for the FSR case
  const bool useFewzWeights = true;
  const bool cutZPT100 = true;
  FEWZ_t fewz(useFewzWeights,cutZPT100);
  if (useFewzWeights && !fewz.isInitialized()) {
    std::cout << "failed to prepare FEWZ correction\n";
    throw 2;
  }

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================
  
  //  
  // Set up histograms
  //
  vector<TH1F*> hZMassv;//, hZMass2v, hZPtv, hZPt2v, hZyv, hZPhiv;  
  
  Double_t   nZv = 0;
  Double_t nZv_puUnweighted=0, nZv_puWeighted=0;

  //int nYBinsMax=findMaxYBins();

  TMatrixD nEventsv (DYTools::nMassBins,DYTools::nYBinsMax);  
  TMatrixD nPassv (DYTools::nMassBins,DYTools::nYBinsMax);  
  TMatrixD effv (DYTools::nMassBins,DYTools::nYBinsMax);  
  TMatrixD effErrv (DYTools::nMassBins,DYTools::nYBinsMax);  

  TMatrixD nEventsBBv (DYTools::nMassBins,DYTools::nYBinsMax);
  TMatrixD nEventsBEv (DYTools::nMassBins,DYTools::nYBinsMax);
  TMatrixD nEventsEEv (DYTools::nMassBins,DYTools::nYBinsMax);
 
  TMatrixD nPassBBv(DYTools::nMassBins,DYTools::nYBinsMax), 
    effBBv(DYTools::nMassBins,DYTools::nYBinsMax), effErrBBv(DYTools::nMassBins,DYTools::nYBinsMax); 
  TMatrixD nPassBEv(DYTools::nMassBins,DYTools::nYBinsMax), 
    effBEv(DYTools::nMassBins,DYTools::nYBinsMax), effErrBEv(DYTools::nMassBins,DYTools::nYBinsMax); 
  TMatrixD nPassEEv(DYTools::nMassBins,DYTools::nYBinsMax), 
    effEEv(DYTools::nMassBins,DYTools::nYBinsMax), effErrEEv(DYTools::nMassBins,DYTools::nYBinsMax);

  TVectorD nEventsZPeakPU(DYTools::nPVBinCount), nPassZPeakPU(DYTools::nPVBinCount);
  TVectorD nEventsZPeakPURaw(100), nPassZPeakPURaw(100);
  TVectorD effZPeakPU(DYTools::nPVBinCount), effErrZPeakPU(DYTools::nPVBinCount);

  nEventsv   = 0;
  nEventsBBv = 0;
  nEventsBEv = 0;
  nEventsEEv = 0;
  nPassv   = 0;
  nPassBBv = 0;
  nPassBEv = 0;
  nPassEEv = 0;

  TMatrixD sumWeightsPassSq (DYTools::nMassBins,DYTools::nYBinsMax);
  TMatrixD sumWeightsTotaSq (DYTools::nMassBins,DYTools::nYBinsMax);
  sumWeightsPassSq = 0;
  sumWeightsTotaSq = 0;


  char hname[100];
  for(UInt_t ifile = 0; ifile<fnamev.size(); ifile++) {
    sprintf(hname,"hZMass_%i",ifile); 
    hZMassv.push_back(new TH1F(hname,"",500,0,1500)); 
    hZMassv[ifile]->Sumw2();
  }

#ifdef usePUReweight
  PUReweight_t puReweight;
  // The commented out code below is not necessary for the Hildreth pile-up reweighting
//   int res=puReweight.setDefaultFile(dirTag,DYTools::analysisTag_USER,0);
//   assert(res);
//   TString outNamePV=puReweight.fileName();
//   TString refPUDistribution="hNGoodPV_data";
//   TString currentPUDistribution="hNGoodPV_zee";
//   res= puReweight.setReference(refPUDistribution) && 
//     puReweight.setActiveSample(currentPUDistribution);
//   if (!res) {
//     std::cout << "failed to locate needed distributions in <" << outNamePV << ">\n";
//     assert(res);
//   }
#endif

  //
  // Access samples and fill histograms
  //  
  TFile *infile=0;
  TTree *eventTree=0;  

  // Data structures to store info from TTrees
  mithep::TEventInfo    *info = new mithep::TEventInfo();
  mithep::TGenInfo *gen  = new mithep::TGenInfo();
  TClonesArray *dielectronArr = new TClonesArray("mithep::TDielectron");
  TClonesArray *pvArr         = new TClonesArray("mithep::TVertex");

  // loop over samples  
  int countMismatch = 0;

  //number of events for which we have binning problem
  int binProblem=0;

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
    AdjustXSectionForSkim(infile,xsecv[ifile],eventTree->GetEntries(),1);
    lumiv[ifile] = eventTree->GetEntries()/xsecv[ifile];
    double extraScale=1; // 4839*1666/27166257.; MC Zee weight in selectEvents
    double scale = extraScale * lumiv[0]/lumiv[ifile];
    cout << "       -> sample weight is " << scale << endl;

    // Set branch address to structures that will store the info  
    eventTree->SetBranchAddress("Info",&info);                
    TBranch *infoBr       = eventTree->GetBranch("Info");
    eventTree->SetBranchAddress("Gen",&gen);                  
    TBranch *genBr = eventTree->GetBranch("Gen");
    eventTree->SetBranchAddress("Dielectron",&dielectronArr); 
    TBranch *dielectronBr = eventTree->GetBranch("Dielectron");
    eventTree->SetBranchAddress("PV", &pvArr);                
    TBranch *pvBr = eventTree->GetBranch("PV");

    // loop over events    
    nZv += scale * eventTree->GetEntries();

    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
      if (debugMode && (ientry>10000)) break;
      //if (ientry>100) break;
      genBr->GetEntry(ientry);
      infoBr->GetEntry(ientry);

      // For EPS2011 for both data and MC (starting from Summer11 production)
      // we use an OR of the twi triggers below. Both are unpresecaled.
      // ULong_t eventTriggerBit = kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL 
      // | kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL;
      // ULong_t leadingTriggerObjectBit = kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele1Obj
      // | kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele1Obj;
      // ULong_t trailingTriggerObjectBit = kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele2Obj
      // | kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele2Obj;

      // old way, to be deleted
      // TriggerConstantSet constantsSet = Full2011DatasetTriggers; // Enum from TriggerSelection.hh

      const bool isData=kFALSE;
      TriggerConstantSet constantsSet = DetermineTriggerSet(triggerSetString);
      assert ( constantsSet != TrigSet_UNDEFINED );

      // Construct the trigger object
      TriggerSelection requiredTriggers(constantsSet, isData, info->runNum);
      assert(requiredTriggers.isDefined());

      // The statements below are commented out: now trigger bit cuts are done
      // in TriggerSelection
//       ULong_t eventTriggerBit = requiredTriggers.getEventTriggerBit(info->runNum);
//       ULong_t leadingTriggerObjectBit = requiredTriggers.getLeadingTriggerObjectBit(info->runNum);
//       ULong_t trailingTriggerObjectBit = requiredTriggers.getTrailingTriggerObjectBit(info->runNum);

      // The A_FSR starting point: gen level quantities
      // The FSR subscript means we work with post-FSR generator level quantities
      double et1 = gen->pt_1;
      double et2 = gen->pt_2;
      double eta1 = gen->eta_1;
      double eta2 = gen->eta_2;

      // Apply acceptance requirement
      if( ! DYTools::goodEtEtaPair( et1, eta1, et2, eta2 ) ) continue;

      // These events are in acceptance, use them for efficiency denominator
      Bool_t isBGen1 = DYTools::isBarrel(eta1);
      Bool_t isBGen2 = DYTools::isBarrel(eta2);
      // determine number of good vertices
      pvBr->GetEntry(ientry);
      int iPUBin=-1;
      int nGoodPV=-1;
      double puWeight=1.0;
	nGoodPV=0;
	const int new_pv_count_code=1;
	if (new_pv_count_code) {
	  nGoodPV=countGoodVertices(pvArr);
	}
	else {
	for (Int_t ipv=0; ipv<pvArr->GetEntriesFast(); ipv++) {
	  const mithep::TVertex *pv = (mithep::TVertex*)((*pvArr)[ipv]);
	  if(pv->nTracksFit                        < 1)  continue;
	  if(pv->ndof                              < 4)  continue;
	  if(fabs(pv->z)                           > 24) continue;
	  if(sqrt((pv->x)*(pv->x)+(pv->y)*(pv->y)) > 2)  continue;
	  nGoodPV++;
	}
	}
      if ((gen->mass>=60) && (gen->mass<=120)) {
	if (nGoodPV>0) {
           if (nGoodPV<=nEventsZPeakPURaw.GetNoElements()) 
	        nEventsZPeakPURaw[nGoodPV] += scale * gen->weight;
	   iPUBin=DYTools::findPUBin(nGoodPV);
	   //std::cout << "iPUBin=" << iPUBin << ", nGoodPV=" << nGoodPV << "\n";
	   if ((iPUBin!=-1) && (iPUBin < nEventsZPeakPU.GetNoElements())) {
	     nEventsZPeakPU[iPUBin] += scale * gen->weight;
	   }
	   //else {
	   //  std::cout << "error in PU bin indexing iPUBin=" << iPUBin << ", nGoodPV=" << nGoodPV << "\n";
	   //}
	}
	//else std::cout << "nGoodPV=" << nGoodPV << "\n";
      }

#ifdef usePUReweight
      puWeight = puReweight.getWeightHildreth(info->nPU);
      nZv_puUnweighted += scale * gen->weight;
      nZv_puWeighted += scale * gen->weight * puWeight;
#endif
    
      // Use post-FSR generator level mass for binning
      int ibinGenM = DYTools::findMassBin(gen->mass);
      int ibinGenY = DYTools::findAbsYBin(ibinGenM,gen->y);
      double totalWeight= scale * gen->weight * puWeight;
      if (useFewzWeights) totalWeight *= fewz.getWeight(gen->vmass,gen->vpt,gen->vy);

      // Accumulate denominator for efficiency calculations
      if(ibinGenM != -1 && ibinGenY != -1 && ibinGenM < DYTools::nMassBins && ibinGenY < DYTools::nYBins[ibinGenM]){
	nEventsv(ibinGenM,ibinGenY) += totalWeight;
        sumWeightsTotaSq(ibinGenM,ibinGenY) += totalWeight*totalWeight;
	// Split events barrel/endcap using matched supercluster or particle eta
	if(isBGen1 && isBGen2)                                  { nEventsBBv(ibinGenM,ibinGenY) += totalWeight; } 
	else if(!isBGen1 && !isBGen2)                           { nEventsEEv(ibinGenM,ibinGenY) += totalWeight; } 
	else if((isBGen1 && !isBGen2) || (!isBGen1 && isBGen2)) { nEventsBEv(ibinGenM,ibinGenY) += totalWeight; }
      }else
        binProblem++;

      // The line below is replaced by the superseeding method, can be cleaned up
      // if(!(info->triggerBits & eventTriggerBit)) continue;  // no trigger accept? Skip to next event...                                   

      if( !(requiredTriggers.matchEventTriggerBit(info->triggerBits, 
						  info->runNum))) 
	continue;
      
      // loop through dielectrons
      dielectronArr->Clear();
      dielectronBr->GetEntry(ientry);

      for(Int_t i=0; i<dielectronArr->GetEntriesFast(); i++) {
        const mithep::TDielectron *dielectron = (mithep::TDielectron*)((*dielectronArr)[i]);
	
	// Apply selection
	// Eta cuts and Et cuts
	if( ! DYTools::goodEtEtaPair( dielectron->scEt_1, dielectron->scEta_1,
				      dielectron->scEt_2, dielectron->scEta_2 ) ) continue;
        Bool_t isB1 = DYTools::isBarrel(dielectron->scEta_1);
        Bool_t isB2 = DYTools::isBarrel(dielectron->scEta_2);
	
	if( !( (isB1 == isBGen1 && isB2 == isBGen2 ) 
	       || (isB1 == isBGen2 && isB2 == isBGen1 ) ) )
	  countMismatch++;

	// Both electrons must match trigger objects. At least one ordering
	// must match
	if( ! requiredTriggers.matchTwoTriggerObjectsAnyOrder( dielectron->hltMatchBits_1,
							       dielectron->hltMatchBits_2,
							       info->runNum) ) continue;

	// The clause below can be deleted, it is superseeded by new methods
//  	if( ! ( 
//  	       (dielectron->hltMatchBits_1 & leadingTriggerObjectBit && 
//  		dielectron->hltMatchBits_2 & trailingTriggerObjectBit )
//  	       ||
//  	       (dielectron->hltMatchBits_1 & trailingTriggerObjectBit && 
//  		dielectron->hltMatchBits_2 & leadingTriggerObjectBit ) ) ) continue;

	// *** Smurf ID is superseeded by new selection ***
// 	// The Smurf electron ID package is the same as used in HWW analysis
// 	// and contains cuts like VBTF WP80 for pt>20, VBTF WP70 for pt<10
// 	// with some customization, plus impact parameter cuts dz and dxy
//  	if(!passSmurf(dielectron)) continue;

	// The selection below is for the EGM working points from spring 2012
	// recommended for both 2011 and 2012 data
	if( DYTools::energy8TeV == 1){
	  if(!passEGMID2012(dielectron, WP_MEDIUM, info->rhoLowEta)) continue;
	}else{
	  if(!passEGMID2011(dielectron, WP_MEDIUM, info->rhoLowEta)) continue;
	}

        // ******** We have a Z candidate! HURRAY! ******** /

	hZMassv[ifile]->Fill(gen->mass,totalWeight);

	// DEBUG
// 	if(ibinGen == 12)
// 	  printf("Gen mass %f  reco mass %f  scEt_1= %f  scEt_2= %f  scEta_1= %f  scEta_2= %f\n",
// 		 gen->mass, dielectron->mass, dielectron->scEt_1, dielectron->scEt_2,
// 		 dielectron->scEta_1, dielectron->scEta_2);
	
	// Accumulate numerator for efficiency calculations
	if ((nGoodPV>0) && (iPUBin!=-1)) { // -1 may also indicate that the mass was not in Z-peak range
	  if ((nGoodPV>=0) && (nGoodPV<=nEventsZPeakPURaw.GetNoElements())) nPassZPeakPURaw[nGoodPV] += totalWeight;
	  if (iPUBin < nPassZPeakPU.GetNoElements()) {
	    nPassZPeakPU[iPUBin] += totalWeight;
	  }
	  //else {
	  //  std::cout << "error in PU bin indexing\n";
	  //}
	}
	if(ibinGenM != -1 && ibinGenY != -1 && ibinGenM < DYTools::nMassBins && ibinGenY < DYTools::nYBins[ibinGenM]){
	  nPassv(ibinGenM,ibinGenY) += totalWeight;
          sumWeightsPassSq(ibinGenM,ibinGenY) += totalWeight*totalWeight;
	  if(isB1 && isB2)                            { nPassBBv(ibinGenM,ibinGenY) += totalWeight; } 
	  else if(!isB1 && !isB2)                     { nPassEEv(ibinGenM,ibinGenY) += totalWeight; } 
	  else if((isB1 && !isB2) || (!isB1 && isB2)) { nPassBEv(ibinGenM,ibinGenY) += totalWeight; }
	}

      } // end loop over dielectrons
    } // end loop over events 
    delete infile;
    infile=0, eventTree=0;
  } // end loop over files
  delete gen;

  cout << "ERROR: binning problem (" << binProblem <<" events in ECAL gap)"<<endl;

  effv      = 0;
  effErrv   = 0;
  effBBv    = 0;
  effErrBBv = 0;
  effBEv    = 0;
  effErrBEv = 0;
  effEEv    = 0;
  effErrEEv = 0;
  for(int i=0; i<DYTools::nMassBins; i++)
    for(int j=0; j<DYTools::nYBins[i]; j++){
      if(nEventsv(i,j) != 0){
        double nPass, nFail, nPassErr, nFailErr;
        nPass=nPassv(i,j);
        nFail=nEventsv(i,j)-nPassv(i,j); 
        nPassErr=sqrt(sumWeightsPassSq(i,j));
        nFailErr=sqrt(sumWeightsTotaSq(i,j)-sumWeightsPassSq(i,j));
        effv(i,j) = nPassv(i,j)/nEventsv(i,j);
        //effErrv(i,j) = sqrt(effv(i,j)*(1-effv(i,j))/nEventsv(i,j));
        effErrv(i,j) = sqrt(( nFail*nFail * nPassErr*nPassErr + nPass*nPass * nFailErr*nFailErr)) / (nEventsv(i,j)*nEventsv(i,j));
      }
    
      if (nEventsBBv(i,j) != 0) {
        effBBv(i,j) = nPassBBv(i,j)/nEventsBBv(i,j);
        effErrBBv(i,j) = sqrt(effBBv(i,j)*(1-effBBv(i,j))/nEventsBBv(i,j));
      }

      if (nEventsBEv(i,j) != 0) {
        effBEv(i,j) = nPassBEv(i,j)/nEventsBEv(i,j);
        effErrBEv(i,j) = sqrt(effBEv(i,j)*(1-effBEv(i,j))/nEventsBEv(i,j));
      }
      
      if (nEventsEEv(i,j) != 0) {
        effEEv(i,j) = nPassEEv(i,j)/nEventsEEv(i,j);
        effErrEEv(i,j) = sqrt(effEEv(i,j)*(1-effEEv(i,j))/nEventsEEv(i,j));
      }
    };

  effZPeakPU=0; effErrZPeakPU=0;
  for (int i=0; i<DYTools::nPVBinCount; ++i) {
    effZPeakPU[i]= nPassZPeakPU[i]/nEventsZPeakPU[i];
    effErrZPeakPU[i]= sqrt(effZPeakPU[i]*(1-effZPeakPU[i])/nEventsZPeakPU[i]);
  }

  printf("Sanity check: gen vs reco barrel-endcap assignment mismatches: %d\n",countMismatch);

  //--------------------------------------------------------------------------------------------------------------
  // Make plots 
  //==============================================================================================================  

  // destination dir
  TString extraTag= DYTools::analysisTag;
  if (!useFewzWeights) extraTag.Append("_noFEWZ");
  CPlot::sOutDir="plots" + extraTag;
  TString outputDir(TString("../root_files/constants/")+dirTag);
  gSystem->mkdir(outputDir,kTRUE);
  TString fnamePlots=outputDir + TString("/event_efficiency_plots") + extraTag + TString(".root");
  TFile *filePlots=new TFile(fnamePlots,"recreate");
  if (!filePlots || !filePlots->IsOpen()) {
    std::cout << "failed to create a file <" << fnamePlots << ">\n";
    throw 2;
  }

  TCanvas *c = MakeCanvas("canvEfficiency","canvEfficiency",800,600);

  // string buffers
  char ylabel[50];   // y-axis label

  // Z mass
  sprintf(ylabel,"a.u. / %.1f GeV/c^{2}",hZMassv[0]->GetBinWidth(1));
  CPlot plotZMass1("zmass1","","m(Z) [GeV/c^{2}]",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    plotZMass1.SetYRange(1, 10000000);
    plotZMass1.AddHist1D(hZMassv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotZMass1.SetLogy();
  plotZMass1.Draw(c);
  SaveCanvas(c, "zmass1");

  PlotMatrixVariousBinning(effv, "efficiency", "LEGO2", filePlots);
  filePlots->Close();
  if (DYTools::study2D==0)
    Plot1D(effv,effErrv,"efficiency1D","efficiency");

  // Store constants in the file
  //TString effConstFileName(outputDir+TString("/event_efficiency_constants.root"));
  TString effConstFileName(outputDir+TString("/event_efficiency_constants") + extraTag + TString(".root"));

  TMatrixD nEventsvErr=nEventsv;
  TMatrixD nEventsBBvErr=nEventsBBv;
  TMatrixD nEventsBEvErr=nEventsBEv;
  TMatrixD nEventsEEvErr=nEventsEEv;
  TMatrixD nPassBBvErr=nPassBBv;
  TMatrixD nPassBEvErr=nPassBEv;
  TMatrixD nPassEEvErr=nPassEEv;

  TMatrixD nPassvErr=nPassv;
  for (int i=0; i<DYTools::nMassBins; ++i) {
    for (int j=0; j<DYTools::nYBins[i]; ++j) {
      nEventsvErr(i,j)=sqrt(sumWeightsTotaSq(i,j));
    }
  }
  for (int i=0; i<DYTools::nMassBins; ++i) {
    for (int j=0; j<DYTools::nYBins[i]; ++j) {
      nPassvErr(i,j)=sqrt(sumWeightsPassSq(i,j));
    }
  }
  nEventsBBvErr=0;
  nEventsBEvErr=0;
  nEventsEEvErr=0;
  nPassBBvErr=0;
  nPassBEvErr=0;
  nPassEEvErr=0;

   TFile fa(effConstFileName,"recreate");
   effv.Write("efficiencyArray");
   effErrv.Write("efficiencyErrArray");

   /*
   nPassv.Write("effEval_nPass");
   nPassvErr.Write("effEval_nPassErr");
   nEventsv.Write("effEval_nTotal");
   nEventsvErr.Write("effEval_nTotalErr");

   nEventsBBv.Write("effEval_nTotalBB");
   nEventsBBvErr.Write("effEval_nTotalBBErr");
   nEventsBEv.Write("effEval_nTotalBE");
   nEventsBEvErr.Write("effEval_nTotalBEErr");
   nEventsEEv.Write("effEval_nTotalEE");
   nEventsEEvErr.Write("effEval_nTotalEEErr");
   nPassBBv.Write("effEval_nPassBB");
   nPassBBvErr.Write("effEval_nPassBBErr");
   nPassBEv.Write("effEval_nPassBE");
   nPassBEvErr.Write("effEval_nPassBEErr");
   nPassEEv.Write("effEval_nPassEE");
   nPassEEvErr.Write("effEval_nPassEEErr");

   effBBv.Write("efficiencyBB");
   effErrBBv.Write("efficiencyBBErr");
   effBEv.Write("efficiencyBE");
   effErrBEv.Write("efficiencyBEErr");
   effEEv.Write("efficiencyEE");
   effErrEEv.Write("efficiencyEEErr");
   */
   //effZPeakPU.Write("efficiencyZPeakPUArray");
   //effErrZPeakPU.Write("efficiencyErrZPeakPUArray");
   //nEventsZPeakPU.Write("nEventsZPeakPUArray");
   //nPassZPeakPU.Write("nPassZPeakPUArray");
   //nEventsZPeakPURaw.Write("nEventsZPeakPURawArray");
   //nPassZPeakPURaw.Write("nPassZPeakPURawArray");
   //unfolding::writeBinningArrays(fa);
   fa.Close();
     
  //--------------------------------------------------------------------------------------------------------------
  // Summary print out
  //==============================================================================================================
  cout << endl;
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  cout << endl; 
  
  cout << labelv[0] << " file: " << fnamev[0] << endl;
  printf("     Number of generated events: %8.1lf\n",nZv);
#ifdef usePUReweight
  printf("     Number of selected events (puUnweighted): %8.1lf\n",nZv_puUnweighted);
  printf("     Number of selected events (puWeighted)  : %8.1lf\n",nZv_puWeighted);
#endif

  const char *yRangeStr=(DYTools::study2D) ? "rapidity range" : "";
  printf(" mass range  %s preselected      passed     total_Eff        BB-BB_Eff        EB-BB_Eff        EB-EB_Eff\n",yRangeStr);
  for(int i=0; i<DYTools::nMassBins; i++){
    double *rapidityBinLimits=DYTools::getYBinLimits(i);
    for (int yi=0; yi<DYTools::nYBins[i]; ++yi) {
      printf(" %4.0f-%4.0f ", DYTools::massBinLimits[i], DYTools::massBinLimits[i+1]);
      if (DYTools::study2D!=0) printf(" %4.2f-%4.2f ", rapidityBinLimits[yi], rapidityBinLimits[yi+1]);
      printf("    %10.0f   %10.0f   %7.4f+-%6.4f  %7.4f+-%6.4f  %7.4f+-%6.4f  %7.4f+-%6.4f \n",
	     nEventsv(i,yi), nPassv(i,yi),
	     effv(i,yi), effErrv(i,yi),
	     effBBv(i,yi), effErrBBv(i,yi),
	     effBEv(i,yi), effErrBEv(i,yi),
	     effEEv(i,yi), effErrEEv(i,yi));
    }
    delete rapidityBinLimits;
  }

  if (1) {
    printf(" mass_range  %s total_preselected  total_passed     BB-BB_preselected  passed    EB-BB_preselected  passed        EB-EB_preselected  passed\n",yRangeStr);
    for(int i=0; i<DYTools::nMassBins; i++){
      double *rapidityBinLimits=DYTools::getYBinLimits(i);
      for (int yi=0; yi<DYTools::nYBins[i]; ++yi) {
	printf(" %4.0f-%4.0f ", DYTools::massBinLimits[i], DYTools::massBinLimits[i+1]);
	if (DYTools::study2D!=0) printf(" %4.2f-%4.2f ", rapidityBinLimits[yi], rapidityBinLimits[yi+1]);
	printf("    %8.1f %8.1f   %8.1f %8.1f  %8.1f %8.1f  %8.1f %8.1f\n",
	       nEventsv(i,yi), nPassv(i,yi),
	       nEventsBBv(i,yi), nPassBBv(i,yi),
	       nEventsBEv(i,yi), nPassBEv(i,yi),
	       nEventsEEv(i,yi), nPassEEv(i,yi));
      }
      delete rapidityBinLimits;
    }
  }

  //printout to the txtFile in latex format
  if (DYTools::study2D)
    {
      latexPrintoutEfficiency2D(effv,effErrv,"Efficiency/plotDYEfficiency.C");
    }
  else if (DYTools::study2D==0)
    {
      latexPrintoutEfficiency1D(effv,effErrv,"Efficiency/plotDYEfficiency.C");
    }

  cout << endl;
  std::cout<<"printout in the Latex format is saved to the text file"<<std::endl;

  printf("\n\nZ-peak efficiency\n");
  printf(" PU bin    preselected      passed     total_Eff\n");
  for(int i=0; i<DYTools::nPVBinCount; i++){
    printf(" %4.0f-%4.0f   %10.0f   %10.0f   %7.4f+-%6.4f\n",
	   DYTools::nPVLimits[i], DYTools::nPVLimits[i+1],
	   nEventsZPeakPU[i], nPassZPeakPU[i],
	   effZPeakPU[i], effErrZPeakPU[i]);
  }

  //sanity check printout
  printSanityCheck(effv, effErrv, "eff");

  cout << endl;
  
  gBenchmark->Show("plotDYEfficiency");
}


//=== FUNCTION DEFINITIONS ======================================================================================

//--------------------------------------------------------------------------------------------------
