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

// define classes and constants to read in ntuple
#include "../Include/EWKAnaDefs.hh"
#include "../Include/TEventInfo.hh"
#include "../Include/TGenInfo.hh"
#include "../Include/TDielectron.hh"   
#include "../Include/TriggerSelection.hh"
#include "../Include/TVertex.hh"

// Helper functions for Electron ID selection
#include "../Include/EleIDCuts.hh"
#endif

//=== FUNCTION DECLARATIONS ======================================================================================

//=== MAIN MACRO =================================================================================================

void plotDYEfficiency(const TString input) 
{
  gBenchmark->Start("plotDYEfficiency");

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
  
  //const Double_t kGAP_LOW  = 1.4442; in DYTools
  //const Double_t kGAP_HIGH = 1.566;  in DYTools


  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================
  
  //  
  // Set up histograms
  //
  vector<TH1F*> hZMassv;//, hZMass2v, hZPtv, hZPt2v, hZyv, hZPhiv;  
  
  Double_t   nZv = 0;

  int nYBinsMax=findMaxYBins();

  TMatrixD nEventsv (DYTools::nMassBins2D,nYBinsMax);  
  TMatrixD nPassv (DYTools::nMassBins2D,nYBinsMax);  
  TMatrixD effv (DYTools::nMassBins2D,nYBinsMax);  
  TMatrixD effErrv (DYTools::nMassBins2D,nYBinsMax);  

  TMatrixD nEventsBBv (DYTools::nMassBins2D,nYBinsMax);
  TMatrixD nEventsBEv (DYTools::nMassBins2D,nYBinsMax);
  TMatrixD nEventsEEv (DYTools::nMassBins2D,nYBinsMax);
 
  TMatrixD nPassBBv(DYTools::nMassBins2D,nYBinsMax), effBBv(DYTools::nMassBins2D,nYBinsMax), effErrBBv(DYTools::nMassBins2D,nYBinsMax); 
  TMatrixD nPassBEv(DYTools::nMassBins2D,nYBinsMax), effBEv(DYTools::nMassBins2D,nYBinsMax), effErrBEv(DYTools::nMassBins2D,nYBinsMax); 
  TMatrixD nPassEEv(DYTools::nMassBins2D,nYBinsMax), effEEv(DYTools::nMassBins2D,nYBinsMax), effErrEEv(DYTools::nMassBins2D,nYBinsMax);

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

  char hname[100];
  for(UInt_t ifile = 0; ifile<fnamev.size(); ifile++) {
    sprintf(hname,"hZMass_%i",ifile); hZMassv.push_back(new TH1F(hname,"",500,0,1500)); hZMassv[ifile]->Sumw2();
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
    lumiv[ifile] = eventTree->GetEntries()/xsecv[ifile];
    double scale = lumiv[0]/lumiv[ifile];
    cout << "       -> sample weight is " << scale << endl;

    // Set branch address to structures that will store the info  
    eventTree->SetBranchAddress("Info",&info);                TBranch *infoBr       = eventTree->GetBranch("Info");
    eventTree->SetBranchAddress("Gen",&gen);                  TBranch *genBr = eventTree->GetBranch("Gen");
    eventTree->SetBranchAddress("Dielectron",&dielectronArr); TBranch *dielectronBr = eventTree->GetBranch("Dielectron");
    eventTree->SetBranchAddress("PV", &pvArr);                TBranch *pvBr = eventTree->GetBranch("PV");

    // loop over events    
    nZv += scale * eventTree->GetEntries();

    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
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

      const bool isData=kFALSE;
      TriggerConstantSet constantsSet = Full2011DatasetTriggers; // Enum from TriggerSelection.hh
      TriggerSelection requiredTriggers(constantsSet, isData, info->runNum);
      ULong_t eventTriggerBit = requiredTriggers.getEventTriggerBit(info->runNum);
      ULong_t leadingTriggerObjectBit = requiredTriggers.getLeadingTriggerObjectBit(info->runNum);
      ULong_t trailingTriggerObjectBit = requiredTriggers.getTrailingTriggerObjectBit(info->runNum);

      // The A_FSR starting point: gen level quantities
      // The FSR subscript means we work with post-FSR generator level quantities
      double et1 = gen->pt_1;
      double et2 = gen->pt_2;
      double eta1 = gen->eta_1;
      double eta2 = gen->eta_2;

      // Apply acceptance requirement
      if( ! ( ( et1 > DYTools::etMinLead  && et2 > DYTools::etMinTrail)
	      || ( et1 > DYTools::etMinTrail  && et2 > DYTools::etMinLead) )) continue;
      if((fabs(eta1) >= 2.5)      || (fabs(eta2) >= 2.5))       continue;
      if((fabs(eta1) >= DYTools::kGAP_LOW) && (fabs(eta1) <= DYTools::kGAP_HIGH)) continue;
      if((fabs(eta2) >= DYTools::kGAP_LOW) && (fabs(eta2) <= DYTools::kGAP_HIGH)) continue;

      // These events are in acceptance, use them for efficiency denominator
      Bool_t isBGen1 = (fabs(eta1)<DYTools::kGAP_LOW);
      Bool_t isBGen2 = (fabs(eta2)<DYTools::kGAP_LOW);
      // determine number of good vertices
      pvBr->GetEntry(ientry);
      int iPUBin=-1;
      int nGoodPV=-1;
      if ((gen->mass>=60) && (gen->mass<=120)) {
	nGoodPV=0;
	for (Int_t ipv=0; ipv<pvArr->GetEntriesFast(); ipv++) {
	  const mithep::TVertex *pv = (mithep::TVertex*)((*pvArr)[ipv]);
	  if(pv->nTracksFit                        < 1)  continue;
	  if(pv->ndof                              < 4)  continue;
	  if(fabs(pv->z)                           > 24) continue;
	  if(sqrt((pv->x)*(pv->x)+(pv->y)*(pv->y)) > 2)  continue;
	  nGoodPV++;
	}
	if (nGoodPV>0) {
           if (nGoodPV<=nEventsZPeakPURaw.GetNoElements()) 
	        nEventsZPeakPURaw[nGoodPV] += scale * gen->weight;
	iPUBin=DYTools::findMassBin(double(nGoodPV),DYTools::nPVBinCount,DYTools::nPVLimits);
	//std::cout << "iPUBin=" << iPUBin << ", nGoodPV=" << nGoodPV << "\n";
	if ((iPUBin!=-1) && (iPUBin < nEventsZPeakPU.GetNoElements())) {
	  nEventsZPeakPU[iPUBin] += scale * gen->weight;
	}
	else {
	  std::cout << "error in PU bin indexing iPUBin=" << iPUBin << ", nGoodPV=" << nGoodPV << "\n";
	}
      } 
	//else std::cout << "nGoodPV=" << nGoodPV << "\n";
      }

      // Use post-FSR generator level mass for binning
      int ibinGenM = DYTools::findMassBin2D(gen->mass);
      int ibinGenY = DYTools::findAbsYBin2D(ibinGenM,gen->y);

      // Accumulate denominator for efficiency calculations
      if(ibinGenM != -1 && ibinGenY != -1 && ibinGenM <= DYTools::nMassBins2D && ibinGenY <= nYBins[ibinGenM]){
	nEventsv(ibinGenM,ibinGenY) += scale * gen->weight;
	// Split events barrel/endcap using matched supercluster or particle eta
	if(isBGen1 && isBGen2)                                  { nEventsBBv(ibinGenM,ibinGenY) += scale * gen->weight; } 
	else if(!isBGen1 && !isBGen2)                           { nEventsEEv(ibinGenM,ibinGenY) += scale * gen->weight; } 
	else if((isBGen1 && !isBGen2) || (!isBGen1 && isBGen2)) { nEventsBEv(ibinGenM,ibinGenY) += scale * gen->weight; }
      }else
        binProblem++;

      if(!(info->triggerBits & eventTriggerBit)) continue;  // no trigger accept? Skip to next event...                                   

      // loop through dielectrons
      dielectronArr->Clear();
      dielectronBr->GetEntry(ientry);

      for(Int_t i=0; i<dielectronArr->GetEntriesFast(); i++) {
        const mithep::TDielectron *dielectron = (mithep::TDielectron*)((*dielectronArr)[i]);
	
	// Apply selection

	// Eta cuts
        if((fabs(dielectron->scEta_1)>DYTools::kGAP_LOW) && (fabs(dielectron->scEta_1)<DYTools::kGAP_HIGH)) continue;
        if((fabs(dielectron->scEta_2)>DYTools::kGAP_LOW) && (fabs(dielectron->scEta_2)<DYTools::kGAP_HIGH)) continue;
	if((fabs(dielectron->scEta_1) > 2.5)       || (fabs(dielectron->scEta_2) > 2.5))       continue;  // outside eta range? Skip to next event...
	
        Bool_t isB1 = (fabs(dielectron->scEta_1)<DYTools::kGAP_LOW);
        Bool_t isB2 = (fabs(dielectron->scEta_2)<DYTools::kGAP_LOW);         

	if( !( (isB1 == isBGen1 && isB2 == isBGen2 ) 
	       || (isB1 == isBGen2 && isB2 == isBGen1 ) ) )
	  countMismatch++;

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

        // ******** We have a Z candidate! HURRAY! ******** /

	hZMassv[ifile]->Fill(gen->mass,scale * gen->weight);

	// DEBUG
// 	if(ibinGen == 12)
// 	  printf("Gen mass %f  reco mass %f  scEt_1= %f  scEt_2= %f  scEta_1= %f  scEta_2= %f\n",
// 		 gen->mass, dielectron->mass, dielectron->scEt_1, dielectron->scEt_2,
// 		 dielectron->scEta_1, dielectron->scEta_2);
	
	// Accumulate numerator for efficiency calculations
	if ((nGoodPV>0) && (iPUBin!=-1)) { // -1 may also indicate that the mass was not in Z-peak range
	  if ((nGoodPV>=0) && (nGoodPV<=nEventsZPeakPURaw.GetNoElements())) nPassZPeakPURaw[nGoodPV] += scale * gen->weight;
	  if (iPUBin < nPassZPeakPU.GetNoElements()) {
	    nPassZPeakPU[iPUBin] += scale * gen->weight;
	  }
	  else {
	    std::cout << "error in PU bin indexing\n";
	  }
	}
	if(ibinGenM != -1 && ibinGenY != -1 && ibinGenM <= DYTools::nMassBins2D && ibinGenY <= nYBins[ibinGenM]){
	  nPassv(ibinGenM,ibinGenY) += scale * gen->weight;
	  if(isB1 && isB2)                            { nPassBBv(ibinGenM,ibinGenY) += scale * gen->weight; } 
	  else if(!isB1 && !isB2)                     { nPassEEv(ibinGenM,ibinGenY) += scale * gen->weight; } 
	  else if((isB1 && !isB2) || (!isB1 && isB2)) { nPassBEv(ibinGenM,ibinGenY) += scale * gen->weight; }
	}

      } // end loop over dielectrons
    } // end loop over events 
    delete infile;
    infile=0, eventTree=0;
  } // end loop over files
  delete gen;
  cout << "ERROR: binning problem. (for " << binProblem <<" events)"<<endl;

  effv      = 0;
  effErrv   = 0;
  effBBv    = 0;
  effErrBBv = 0;
  effBEv    = 0;
  effErrBEv = 0;
  effEEv    = 0;
  effErrEEv = 0;
  for(int i=0; i<DYTools::nMassBins2D; i++)
    for(int j=0; j<DYTools::nYBins[i]; j++){
      if(nEventsv(i,j) != 0){
        effv(i,j) = nPassv(i,j)/nEventsv(i,j);
        effErrv(i,j) = sqrt(effv(i,j)*(1-effv(i,j))/nEventsv(i,j)); 
    
        effBBv(i,j) = nPassBBv(i,j)/nEventsBBv(i,j);
        effErrBBv(i,j) = sqrt(effBBv(i,j)*(1-effBBv(i,j))/nEventsBBv(i,j));
      
        effBEv(i,j) = nPassBEv(i,j)/nEventsBEv(i,j);
        effErrBEv(i,j) = sqrt(effBEv(i,j)*(1-effBEv(i,j))/nEventsBEv(i,j));
      
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
  TCanvas *c = MakeCanvas("c","c",800,600);

  // string buffers
  char ylabel[50];   // y-axis label

  // Z mass
  sprintf(ylabel,"a.u. / %.1f GeV/c^{2}",hZMassv[0]->GetBinWidth(1));
  CPlot plotZMass1("zmass1","","m(Z) [GeV/c^{2}]",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    plotZMass1.AddHist1D(hZMassv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotZMass1.SetLogy();
  plotZMass1.Draw(c);
  SaveCanvas(c, "zmass1");

  PlotMatrixVariousBinning(effv, "efficiency", "LEGO2");

  // Store constants in the file
  TString outputDir(TString("../root_files/constants/")+dirTag);
  gSystem->mkdir(outputDir,kTRUE);
  TString effConstFileName(outputDir+TString("/event_efficiency_constants.root"));

   TFile fa(effConstFileName,"recreate");
   effv.Write("efficiencyArray");
   effErrv.Write("efficiencyErrArray");
   effZPeakPU.Write("efficiencyZPeakPUArray");
   effErrZPeakPU.Write("efficiencyErrZPeakPUArray");
   nEventsZPeakPU.Write("nEventsZPeakPUArray");
   nPassZPeakPU.Write("nPassZPeakPUArray");
   nEventsZPeakPURaw.Write("nEventsZPeakPURawArray");
   nPassZPeakPURaw.Write("nPassZPeakPURawArray");
   fa.Close();
     
  //--------------------------------------------------------------------------------------------------------------
  // Summary print out
  //==============================================================================================================
  //cout << endl;
  //cout << "*" << endl;
  //cout << "* SUMMARY" << endl;
  //cout << "*--------------------------------------------------" << endl;
  //cout << endl; 
  
  //cout << labelv[0] << " file: " << fnamev[0] << endl;
  //printf("     Number of generated events: %8.1lf",nZv);
  //printf(" mass bin    preselected      passed     total_Eff        BB-BB_Eff        EB-BB_Eff        EB-EB_Eff   \n");
  //for(int i=0; i<DYTools::nMassBins; i++){
  //  printf(" %4.0f-%4.0f   %10.0f   %10.0f   %7.4f+-%6.4f  %7.4f+-%6.4f  %7.4f+-%6.4f  %7.4f+-%6.4f \n",
	 //  DYTools::massBinLimits[i], DYTools::massBinLimits[i+1],
	  // nEventsv[i], nPassv[i],
	  // effv[i], effErrv[i],
	  // effBBv[i], effErrBBv[i],
	  // effBEv[i], effErrBEv[i],
	  // effEEv[i], effErrEEv[i]);
  //}

  printf("\n\nZ-peak efficiency\n");
  printf(" PU bin    preselected      passed     total_Eff\n");
  for(int i=0; i<DYTools::nPVBinCount; i++){
    printf(" %4.0f-%4.0f   %10.0f   %10.0f   %7.4f+-%6.4f\n",
	   DYTools::nPVLimits[i], DYTools::nPVLimits[i+1],
	   nEventsZPeakPU[i], nPassZPeakPU[i],
	   effZPeakPU[i], effErrZPeakPU[i]);
  }

  cout << endl;
  
  gBenchmark->Show("plotDYEfficiency");
}


//=== FUNCTION DEFINITIONS ======================================================================================

//--------------------------------------------------------------------------------------------------