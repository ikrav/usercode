#include <TFile.h>
#include <TChain.h>
#include <TString.h>
#include <TBenchmark.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TRandom.h>
#include <TVectorD.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TTimeStamp.h>
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <vector>                   // STL vector class
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings

#include <vector>

#include "../Include/CPlot.hh"
#include "../Include/MitStyleRemix.hh"

#include "../Include/DYTools.hh"
#include "../Include/EleIDCuts.hh"

#include "../Include/EWKAnaDefs.hh"
#include "../Include/TGenInfo.hh"
#include "../Include/TEventInfo.hh"
#include "../Include/TDielectron.hh"
#include "../Include/TElectron.hh"
#include "../Include/TriggerSelection.hh"
#include "../Include/cutFunctions.hh"

using namespace mithep;

//=== FUNCTION DECLARATIONS ======================================================================================

void fillEfficiencyConstants( const TriggerSelection &triggerSet );
void fillOneEfficiency(const TString filename, double *eff, double *effArr, int etaRange);

Bool_t matchedToGeneratorLevel(const TGenInfo *gen, const TDielectron *dielectron);

double findEventScaleFactor(const TElectron *leading, const TElectron *trailing);
double findRecoScaleFactor(const TElectron *ele);
double findIdScaleFactor  (const TElectron *ele);
double findHltScaleFactor (const TElectron *ele);

double findEventScaleFactorSmeared(const TElectron *leading, const TElectron *trailing, int iexp);
double findRecoScaleFactorSmeared(const TElectron *ele, int iexp);
double findIdScaleFactorSmeared  (const TElectron *ele, int iexp);
double findHltScaleFactorSmeared (const TElectron *ele, int iexp);

void drawEfficiencies();
void drawEfficiencyGraphs(TGraphErrors *grData, TGraphErrors *grMc,
			  TString yAxisTitle, TString text, TString plotName);
void drawEfficiencyGraphsAsymmErrors(TGraphAsymmErrors *grData, TGraphErrors *grMc,
			  TString yAxisTitle, TString text, TString plotName);
void drawScaleFactors();
void drawScaleFactorGraphs(TGraphErrors *gr, TString yAxisTitle, TString text, 
			   TString plotName);

void drawEventScaleFactors(TVectorD scaleGsfV, TVectorD scaleGsfErrV,
			   TVectorD scaleIdV , TVectorD scaleIdErrV ,
			   TVectorD scaleHltV, TVectorD scaleHltErrV,
			   TVectorD scaleV   , TVectorD scaleErrV    );
void drawEventScaleFactorGraphs(TGraphErrors *gr, TString yAxisTitle, 
				TString plotName);

double errOnRatio(double a, double da, double b, double db);

//=== Constants ==========================

const bool savePlots = false;

// File names for efficiency measurements from tag and probe
TString          dirTag;

/*
  retire static file definitions
const TString effDataGsfFile = "efficiency_TnP_data_gsf_fit-fit_bins-et5-eta2.root";
const TString effMcGsfFile   = "efficiency_TnP_mc_gsf_count-count_bins-et5-eta2.root";

const TString effDataIdFile  = "efficiency_TnP_data_id_fit-fit_bins-et5-eta2.root";
const TString effMcIdFile    = "efficiency_TnP_mc_id_count-count_bins-et5-eta2.root";

const TString effDataHltFile = "efficiency_TnP_data_hlt_count-count_bins-et5-eta2.root";
const TString effMcHltFile   = "efficiency_TnP_mc_hlt_count-count_bins-et5-eta2.root";
*/

// Define the method used to obtain the efficiencies
const int dataGsfEffMethod = FITnFIT;
const int mcGsfEffMethod   = COUNTnCOUNT;
const int dataIdEffMethod  = FITnFIT;
const int mcIdEffMethod    = COUNTnCOUNT;
const int dataHltEffMethod = COUNTnCOUNT;
const int mcHltEffMethod   = COUNTnCOUNT;

// Declaration of arrays into which efficiencies will be loaded
const int etBinning = DYTools::ETBINS5;
const int etBinCount = DYTools::nEtBins5;
const double *etBinLimits= DYTools::etBinLimits5;
const int etaBinning = DYTools::ETABINS2;
// Reconstruction
double GsfBarrelDataEff   [etBinCount], GsfBarrelDataEffErr[etBinCount];
double GsfEndcapDataEff   [etBinCount], GsfEndcapDataEffErr[etBinCount];
double GsfBarrelMcEff     [etBinCount], GsfBarrelMcEffErr  [etBinCount];
double GsfEndcapMcEff     [etBinCount], GsfEndcapMcEffErr  [etBinCount];
// Identification
double IdBarrelDataEff   [etBinCount], IdBarrelDataEffErr[etBinCount];
double IdEndcapDataEff   [etBinCount], IdEndcapDataEffErr[etBinCount];
double IdBarrelMcEff     [etBinCount], IdBarrelMcEffErr  [etBinCount];
double IdEndcapMcEff     [etBinCount], IdEndcapMcEffErr  [etBinCount];
// HLT trigger
double HltBarrelDataEff   [etBinCount], HltBarrelDataEffErr[etBinCount];
double HltEndcapDataEff   [etBinCount], HltEndcapDataEffErr[etBinCount];
double HltBarrelMcEff     [etBinCount], HltBarrelMcEffErr  [etBinCount];
double HltEndcapMcEff     [etBinCount], HltEndcapMcEffErr  [etBinCount];

// Get the values from TriggerSelection.hh
const double *runLumiV = luminositiesInRunSections;

// Global variables
const int nexp = 100;
double ro_D_B_reco[nexp];
double ro_D_E_reco[nexp];
double ro_M_B_reco[nexp];
double ro_M_E_reco[nexp];

double ro_D_B_id[nexp];
double ro_D_E_id[nexp];
double ro_M_B_id[nexp];
double ro_M_E_id[nexp];

double ro_D_B_hlt[nexp];
double ro_D_E_hlt[nexp];
double ro_M_B_hlt[nexp];
double ro_M_E_hlt[nexp];

//=== MAIN MACRO =================================================================================================

void calcEventEff(const TString input, TString triggerSetString)
{

  gBenchmark->Start("calcEventEff");
  
   // fast check
  TriggerConstantSet triggerSet=DetermineTriggerSet(triggerSetString);  
  assert ( triggerSet != TrigSet_UNDEFINED );
 
  printf("Et bin count is %d\n", etBinCount);
  if (etaBinning!=ETABINS2) {
    std::cout << "this code assumes that eta regions are only (barrel,endcap)\n";
    throw 1;
  }

  CPlot::sOutDir = "plots";

  vector<TString> fnamev;   // file names   
  vector<TString> labelv;   // legend label
  vector<Int_t>   colorv;   // color in plots
  vector<Int_t>   linev;    // line style
  vector<Double_t> xsecv;
  vector<Double_t> lumiv;

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

  // Construct the trigger object
  TriggerSelection triggers(triggerSetString, false, 0); // we work with MC files

  // Read efficiency constants from ROOT files
  // This has to be done AFTER configuration file is parsed
  fillEfficiencyConstants( triggers );

  TH1F *hScale = new TH1F("hScale", "", 150, 0.0, 1.5);
  TH1F *hScaleGsf = new TH1F("hScaleGsf", "", 150, 0.0, 1.5);
  TH1F *hScaleId  = new TH1F("hScaleId" , "", 150, 0.0, 1.5);
  TH1F *hScaleHlt = new TH1F("hScaleHlt", "", 150, 0.0, 1.5);

  TH1F *hZpeakEt = new TH1F("hZpeakEt", "", etBinCount, etBinLimits);
  vector<TH1F*> hLeadingEtV;
  vector<TH1F*> hTrailingEtV;
  vector<TH1F*> hElectronEtV;

  vector<TH1F*> hScaleV;
  vector<TH1F*> hScaleGsfV;
  vector<TH1F*> hScaleIdV;
  vector<TH1F*> hScaleHltV;
  for(int i=0; i<nMassBins; i++){
    TString base = "hScaleV_";
    base += i;
    hScaleV.push_back(new TH1F(base,base,150,0.0,1.5));
    hScaleGsfV.push_back(new TH1F(base+TString("_gsf"),base+TString("_gsf"),150,0.0,1.5));
    hScaleIdV .push_back(new TH1F(base+TString("_id" ),base+TString("_id" ),150,0.0,1.5));
    hScaleHltV.push_back(new TH1F(base+TString("_hlt"),base+TString("_hlt"),150,0.0,1.5));
    base = "hLeadingEt_";
    base += i;
    hLeadingEtV.push_back(new TH1F(base,base,etBinCount, etBinLimits));
    base = "hTrailingEt_";
    base += i;
    hTrailingEtV.push_back(new TH1F(base,base,etBinCount, etBinLimits));
    base = "hElectronEt_";
    base += i;
    hElectronEtV.push_back(new TH1F(base,base,etBinCount, etBinLimits));
  }

  // Create Gaussian-distributed random offsets for each pseudo-experiment
  for(int i=0; i<nexp; i++){
    ro_D_B_reco[i] = gRandom->Gaus(0.0,1.0);
    ro_D_E_reco[i] = gRandom->Gaus(0.0,1.0);
    ro_M_B_reco[i] = gRandom->Gaus(0.0,1.0);
    ro_M_E_reco[i] = gRandom->Gaus(0.0,1.0);
    
    ro_D_B_id[i] = gRandom->Gaus(0.0,1.0);
    ro_D_E_id[i] = gRandom->Gaus(0.0,1.0);
    ro_M_B_id[i] = gRandom->Gaus(0.0,1.0);
    ro_M_E_id[i] = gRandom->Gaus(0.0,1.0);
    
    ro_D_B_hlt[i] = gRandom->Gaus(0.0,1.0);
    ro_D_E_hlt[i] = gRandom->Gaus(0.0,1.0);
    ro_M_B_hlt[i] = gRandom->Gaus(0.0,1.0);
    ro_M_E_hlt[i] = gRandom->Gaus(0.0,1.0);
  }
  // Create container for data for error estimates based on pseudo-experiments
  TH1F *systScale[nMassBins][nexp];
  TH1F *systScaleGsf[nMassBins][nexp];
  TH1F *systScaleId[nMassBins][nexp];
  TH1F *systScaleHlt[nMassBins][nexp];
  for(int i=0; i<nMassBins; i++)
    for(int j=0; j<nexp; j++){
      TString base = "hScaleM_mass";
      base += i;
      base += "_exp";
      base += j;
      systScale[i][j] = new TH1F(base,base,150,0.0,1.5);
      systScaleGsf[i][j] = new TH1F(base+TString("_gsf"),base+TString("_gsf"),150,0.0,1.5);
      systScaleId [i][j] = new TH1F(base+TString("_id" ),base+TString("_id" ),150,0.0,1.5);
      systScaleHlt[i][j] = new TH1F(base+TString("_hlt"),base+TString("_hlt"),150,0.0,1.5);
    }
  
  
  int eventsInNtuple = 0;
  double weightedEventsInNtuple = 0;
  int eventsAfterTrigger = 0;
  int totalCand = 0;
  int totalCandInMassWindow = 0;
  int totalCandInEtaAcceptance = 0;
  int totalCandEtAbove10GeV = 0;
  int totalCandMatchedToGen = 0;
  int totalCandFullSelection = 0;

  // Loop over files
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++){

    //
    // Access samples and fill histograms
    //  
    TFile *infile = 0;
    TTree *eventTree = 0;
        
    // Data structures to store info from TTrees
    mithep::TEventInfo *info = new mithep::TEventInfo();
    mithep::TGenInfo   *gen  = new mithep::TGenInfo();
    TClonesArray *dielectronArr = new TClonesArray("mithep::TDielectron");
    
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
    eventTree->SetBranchAddress("Dielectron",&dielectronArr); TBranch *dielectronBr = eventTree->GetBranch("Dielectron");
    eventTree->SetBranchAddress("Gen",&gen);                  TBranch *genBr = eventTree->GetBranch("Gen");

    // loop over events    
    eventsInNtuple         += eventTree->GetEntries();
    weightedEventsInNtuple += scale * eventTree->GetEntries();
    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
//       for(UInt_t ientry=0; ientry<100000; ientry++) { // This is for faster turn-around in testing
      
      genBr->GetEntry(ientry);
      infoBr->GetEntry(ientry);
      
      /* old trigger defs
      // For EPS2011 for both data and MC (starting from Summer11 production)
      // we use an OR of the twi triggers below. Both are unpresecaled.
      ULong_t eventTriggerBit = kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL 
	| kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL;
      ULong_t leadingTriggerObjectBit = kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele1Obj
	| kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele1Obj;
      ULong_t trailingTriggerObjectBit = kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele2Obj
	| kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele2Obj;
      */
      ULong_t eventTriggerBit = triggers.getEventTriggerBit(info->runNum);
      ULong_t leadingTriggerObjectBit  = triggers.getLeadingTriggerObjectBit(info->runNum);
      ULong_t trailingTriggerObjectBit = triggers.getTrailingTriggerObjectBit(info->runNum);
      
      if(!(info->triggerBits & eventTriggerBit)) continue;  // no trigger accept? Skip to next event...
      eventsAfterTrigger++;
      
      // loop through dielectrons
      dielectronArr->Clear();
      dielectronBr->GetEntry(ientry);    

      for(Int_t i=0; i<dielectronArr->GetEntriesFast(); i++) {
	
	totalCand++;
	const mithep::TDielectron *dielectron = (mithep::TDielectron*)((*dielectronArr)[i]);

	// Consider only events in the mass range of interest
	// Use generator level post-FSR mass.
 	if( gen->mass < massBinLimits[0] || gen->mass > massBinLimits[nMassBins]) continue;
	totalCandInMassWindow++;

	// Exclude ECAL gap region (should already be done for ntuple, but just to make sure...)
	if((fabs(dielectron->scEta_1)>ECAL_GAP_LOW) && (fabs(dielectron->scEta_1)<ECAL_GAP_HIGH)) continue;
	if((fabs(dielectron->scEta_2)>ECAL_GAP_LOW) && (fabs(dielectron->scEta_2)<ECAL_GAP_HIGH)) continue;
	// ECAL acceptance cut on supercluster Et
	if((fabs(dielectron->scEta_1) > 2.5)       || (fabs(dielectron->scEta_2) > 2.5)) continue;  // outside eta range? Skip to next event...
	totalCandInEtaAcceptance++;
	// None of the electrons should be below 10 GeV
	if((dielectron->scEt_1 < 10)               || (dielectron->scEt_2 < 10))	      continue;  // below supercluster ET cut? Skip to next event...
	totalCandEtAbove10GeV++;

	// For MC-only, do generator level matching
	if( ! matchedToGeneratorLevel(gen, dielectron) ) continue;
	totalCandMatchedToGen++;

	// ECAL-driven reconstruction is not required.

	// No cut on opposite charges to avoid systematics related to charge mis-ID	
	//  	if( (dielectron->q_1 == dielectron->q_2 )) continue;

	TElectron *ele1 = extractElectron(dielectron, 1);
	TElectron *ele2 = extractElectron(dielectron, 2);

	// ID cuts
	if( !( passSmurf(ele1) && passSmurf(ele2) ) ) continue;

	// ET and trigger cut on the leading electron
	const TElectron *leading = ele1;
	const TElectron *trailing = ele2;
	if( ele1->scEt < ele2->scEt ){
	  leading = ele2;
	  trailing = ele1;
	}

	// The individual electron trigger match should be done
	// exactly the same way as it is done in the signal selection.
	// (At the moment of this writing it is a bit different, needs to
	// be fixed).
	if( !( leading->scEt  > 20 && (leading ->hltMatchBits & leadingTriggerObjectBit) ) ) continue;
	if( !( trailing->scEt > 10 && (trailing->hltMatchBits & trailingTriggerObjectBit) ) ) continue;
	totalCandFullSelection++;

	double scaleFactor = findEventScaleFactor(leading, trailing);
	double scaleFactorGsf = sqrt(findRecoScaleFactor(leading)*findRecoScaleFactor(trailing));
	double scaleFactorId  = sqrt(findIdScaleFactor(leading)*findIdScaleFactor(trailing));
	double scaleFactorHlt = sqrt(findHltScaleFactor(leading)*findHltScaleFactor(trailing));
	hScale->Fill(scaleFactor, scale);
	hScaleGsf->Fill( scaleFactorGsf, scale);
	hScaleId ->Fill( scaleFactorId, scale);
	hScaleHlt->Fill( scaleFactorHlt, scale);
	// Use generator-level post-FSR mass 
	int ibin = findMassBin(gen->mass);
	hScaleGsfV[ibin]->Fill( scaleFactorGsf, scale);
	hScaleIdV [ibin]->Fill( scaleFactorId, scale);
	hScaleHltV[ibin]->Fill( scaleFactorHlt, scale);
	hScaleV   [ibin]->Fill( scaleFactor, scale);

	hLeadingEtV [ibin]->Fill( leading->scEt, scale);
	hTrailingEtV[ibin]->Fill( trailing->scEt, scale);
	hElectronEtV[ibin]->Fill( leading->scEt, scale);
	hElectronEtV[ibin]->Fill( trailing->scEt, scale);
	if( dielectron->mass > 60 && dielectron->mass < 120){
	  hZpeakEt->Fill(leading->scEt, scale);
	  hZpeakEt->Fill(trailing->scEt, scale);
	}

	// Acumulate pseudo-experiments for error estimate
	for(int iexp = 0; iexp<nexp; iexp++){
 	  scaleFactor = findEventScaleFactorSmeared(leading, trailing, iexp);
	  scaleFactorGsf = sqrt(findRecoScaleFactorSmeared(leading,iexp)*findRecoScaleFactorSmeared(trailing,iexp));
	  scaleFactorId  = sqrt(findIdScaleFactorSmeared(leading,iexp)*findIdScaleFactorSmeared(trailing,iexp));
 	  scaleFactorHlt = sqrt(findHltScaleFactorSmeared(leading,iexp)*findHltScaleFactorSmeared(trailing,iexp));
	  systScale   [ibin][iexp]->Fill(scaleFactor, scale);
	  systScaleGsf[ibin][iexp]->Fill(scaleFactorGsf, scale);
	  systScaleId[ibin][iexp]->Fill(scaleFactorId, scale);
	  systScaleHlt[ibin][iexp]->Fill(scaleFactorHlt, scale);
	}

// 	if(scaleFactor>1.3)
// 	  printf("  leading:   %f    %f      trailing:   %f   %f     mass: %f\n",
// 		 leading->scEt, leading->scEta, trailing->scEt, trailing->scEta, dielectron->mass);


      } // end loop over dielectrons
    } // end loop over events
    
    delete infile;
    infile = 0;
    eventTree = 0;
    delete gen;
    delete info;
    delete dielectronArr;
  } // end loop over files
  
  // Calculate errors on the scale factors
  // The "Mean" are the mean among all pseudo-experiments, very close to the primary scale factor values
  TVectorD scaleMeanV(nMassBins);
  TVectorD scaleMeanErrV(nMassBins);
  TVectorD scaleMeanGsfV(nMassBins);
  TVectorD scaleMeanGsfErrV(nMassBins);
  TVectorD scaleMeanIdV(nMassBins);
  TVectorD scaleMeanIdErrV(nMassBins);
  TVectorD scaleMeanHltV(nMassBins);
  TVectorD scaleMeanHltErrV(nMassBins);
  // Put into these vectors the content of the mean of the primary scale factor distributions
  TVectorD scaleV(nMassBins);
  TVectorD scaleGsfV(nMassBins);
  TVectorD scaleIdV(nMassBins);
  TVectorD scaleHltV(nMassBins);
  for(int ibin = 0; ibin < nMassBins; ibin++){
    scaleMeanV[ibin] = 0;
    scaleMeanErrV[ibin] = 0;
    scaleMeanGsfV[ibin] = 0;
    scaleMeanGsfErrV[ibin] = 0;
    scaleMeanIdV[ibin] = 0;
    scaleMeanIdErrV[ibin] = 0;
    scaleMeanHltV[ibin] = 0;
    scaleMeanHltErrV[ibin] = 0;
    for(int iexp = 0; iexp < nexp; iexp++){
      scaleMeanV[ibin] += systScale[ibin][iexp]->GetMean();
      scaleMeanErrV[ibin] += systScale[ibin][iexp]->GetMean() * systScale[ibin][iexp]->GetMean();

      scaleMeanGsfV[ibin] += systScaleGsf[ibin][iexp]->GetMean();
      scaleMeanGsfErrV[ibin] += systScaleGsf[ibin][iexp]->GetMean() * systScaleGsf[ibin][iexp]->GetMean();

      scaleMeanIdV[ibin] += systScaleId[ibin][iexp]->GetMean();
      scaleMeanIdErrV[ibin] += systScaleId[ibin][iexp]->GetMean() * systScaleId[ibin][iexp]->GetMean();

      scaleMeanHltV[ibin] += systScaleHlt[ibin][iexp]->GetMean();
      scaleMeanHltErrV[ibin] += systScaleHlt[ibin][iexp]->GetMean() * systScaleHlt[ibin][iexp]->GetMean();
    }
    scaleGsfV[ibin] = hScaleGsfV[ibin]->GetMean();
    scaleIdV [ibin] = hScaleIdV [ibin] ->GetMean();
    scaleHltV[ibin] = hScaleHltV[ibin]->GetMean();
    scaleV   [ibin] = hScaleV   [ibin]->GetMean();

    scaleMeanV[ibin] = scaleMeanV[ibin]/nexp;
    scaleMeanErrV[ibin] = sqrt( scaleMeanErrV[ibin] / nexp 
				- scaleMeanV[ibin]*scaleMeanV[ibin] ); 
				
    scaleMeanGsfV[ibin] = scaleMeanGsfV[ibin]/nexp;
    scaleMeanGsfErrV[ibin] = sqrt( scaleMeanGsfErrV[ibin] / nexp 
				- scaleMeanGsfV[ibin]*scaleMeanGsfV[ibin] ); 
				
    scaleMeanIdV[ibin] = scaleMeanIdV[ibin]/nexp;
    scaleMeanIdErrV[ibin] = sqrt( scaleMeanIdErrV[ibin] / nexp 
				- scaleMeanIdV[ibin]*scaleMeanIdV[ibin] ); 
				
    scaleMeanHltV[ibin] = scaleMeanHltV[ibin]/nexp;
    scaleMeanHltErrV[ibin] = sqrt( scaleMeanHltErrV[ibin] / nexp 
				- scaleMeanHltV[ibin]*scaleMeanHltV[ibin] ); 
				
  }

  // Store constants in the file
  TString outputDir(TString("../root_files/constants/")+dirTag);
  gSystem->mkdir(outputDir,kTRUE);
  TString sfConstFileName(outputDir+TString("/scale_factors_") + triggers.triggerConditionsName() + TString(".root"));

  TFile fa(sfConstFileName, "recreate");
  scaleV.Write("scaleFactorArray");
  scaleMeanErrV.Write("scaleFactorErrArray");
  fa.Close();

  printf("Total events in ntuple                                       %15d\n",eventsInNtuple);
  printf("    number of weighted events in ntuple                      %17.1lf\n",weightedEventsInNtuple);
  printf("    events after event level trigger cut                     %15d\n",eventsAfterTrigger);
  printf("\nTotal candidates (no cuts)                                 %15d\n",totalCand);
  printf("        candidates in 15-600 mass window                     %15d\n",totalCandInMassWindow);
  printf("        candidates with eta 0-1.4442, 1.566-2.5              %15d\n",totalCandInEtaAcceptance);
  printf("        candidates, both electrons above 10 GeV              %15d\n",totalCandEtAbove10GeV);
  printf("        candidates matched to GEN level (if MC)              %15d\n",totalCandMatchedToGen);
  printf("        candidates, full selection                           %15d\n",totalCandFullSelection);

  printf("\nTotal average scale factor per event:     %5.3f\n", hScale->GetMean());
  printf(  "     rho_reco per electron                %5.3f\n", hScaleGsf->GetMean());
  printf(  "     rho_id   per electron                %5.3f\n", hScaleId ->GetMean());
  printf(  "     rho_hlt  per electron                %5.3f\n", hScaleHlt->GetMean());

  TCanvas *c1 = new TCanvas("c1","c1",10,10,500,500);
  c1->Divide(2,2);
 
  c1->cd(1);
  hScale->Draw();

  c1->cd(2);
  hScaleGsf->Draw();

  c1->cd(3);
  hScaleId->Draw();

  c1->cd(4);
  hScaleHlt->Draw();

  printf("\nScale factors as a function of mass bin\n");
  printf("    mass          rho_reco          rho_id       rho_hlt        rho_total\n");
  for(int i=0; i<nMassBins; i++){
    printf("   %3.0f - %3.0f     %5.3f +- %5.3f    %5.3f +- %5.3f    %5.3f +- %5.3f     %5.3f +- %5.3f\n",
	   massBinLimits[i], massBinLimits[i+1],
	   hScaleGsfV[i]->GetMean(), scaleMeanGsfErrV[i],
	   hScaleIdV[i]->GetMean(),  scaleMeanIdErrV[i],
	   hScaleHltV[i]->GetMean(), scaleMeanHltErrV[i],
	   hScaleV[i]->GetMean()   , scaleMeanErrV[i]);
  }

  drawEfficiencies();
  drawScaleFactors();
  drawEventScaleFactors(scaleGsfV, scaleMeanGsfErrV,
			scaleIdV , scaleMeanIdErrV ,
			scaleHltV, scaleMeanHltErrV,
			scaleV   , scaleMeanErrV    );

  //
  // Make plots of Et spectra
  //
  // Normalize first
  for(int i=0; i<nMassBins; i++){
    printf("Total events in mass bin %3d     %10.0f\n", i, hLeadingEtV[i]->GetSumOfWeights());
    hLeadingEtV  [i]->Sumw2();
    hTrailingEtV [i]->Sumw2();
    hElectronEtV [i]->Sumw2();
    hLeadingEtV [i]->Scale(1.0/hLeadingEtV[i]->GetSumOfWeights());
    hTrailingEtV[i]->Scale(1.0/hTrailingEtV[i]->GetSumOfWeights());
    hElectronEtV[i]->Scale(1.0/hElectronEtV[i]->GetSumOfWeights());
  }
  printf("Total events around Z peak    %10.0f\n", hZpeakEt->GetSumOfWeights()/2.0);
  hZpeakEt->Sumw2();
  hZpeakEt->Scale(1.0/hZpeakEt->GetSumOfWeights());

  TCanvas *c3 = MakeCanvas("ET canvas 1", "ET canvas 1");
  CPlot etplot1("etplot_bin0", "","E_{T} [GeV]", "N_{ele}, normalized");
  int cbin = 0;
  TString label = "mass bin 0";
  etplot1.SetLogx();
  etplot1.AddHist1D(hElectronEtV[cbin] , label , "PE", kRed);
  etplot1.AddHist1D(hZpeakEt       , "60-120 mass range", "hist,f", kBlack);
  etplot1.SetYRange(0.0,0.6);
  hElectronEtV[cbin]->GetXaxis()->SetMoreLogLabels();
  hElectronEtV[cbin]->GetXaxis()->SetNoExponent();
  hZpeakEt->SetFillStyle(3001);
  hElectronEtV[cbin]->SetMarkerColor(kRed);
  etplot1.Draw(c3,savePlots,"png");

  TCanvas *c4 = MakeCanvas("ET canvas 2", "ET canvas 2");
  CPlot etplot2("etplot_bin3", "","E_{T} [GeV]", "N_{ele}, normalized");
  cbin = 3;
  label = "mass bin 3";
  etplot2.SetLogx();
  etplot2.AddHist1D(hElectronEtV[cbin] , label , "PE", kRed);
  etplot2.AddHist1D(hZpeakEt       , "60-120 mass range", "hist,f", kBlack);
  etplot2.SetYRange(0.0,0.6);
  hElectronEtV[cbin]->GetXaxis()->SetMoreLogLabels();
  hElectronEtV[cbin]->GetXaxis()->SetNoExponent();
  hZpeakEt->SetFillStyle(3001);
  hZpeakEt->SetMarkerColor(kBlack);
  hElectronEtV[cbin]->SetMarkerColor(kRed);
  etplot2.Draw(c4,savePlots,"png");

  TCanvas *c5 = MakeCanvas("ET canvas 3", "ET canvas 3");
  CPlot etplot3("etplot_bin5", "","E_{T} [GeV]", "N_{ele}, normalized");
  cbin = 5;
  label = "mass bin 5";
  etplot3.SetLogx();
  etplot3.AddHist1D(hElectronEtV[cbin] , label , "PE", kRed);
  etplot3.AddHist1D(hZpeakEt       , "60-120 mass range", "hist,f", kBlack);
  etplot3.SetYRange(0.0,0.6);
  hElectronEtV[cbin]->GetXaxis()->SetMoreLogLabels();
  hElectronEtV[cbin]->GetXaxis()->SetNoExponent();
  hZpeakEt->SetFillStyle(3001);
  hZpeakEt->SetMarkerColor(kBlack);
  hElectronEtV[cbin]->SetMarkerColor(kRed);
  etplot3.Draw(c5,savePlots,"png");

  gBenchmark->Show("calcEventEff");
  return;
}

Bool_t matchedToGeneratorLevel(const TGenInfo *gen, const TDielectron *dielectron){

  Bool_t result = kTRUE;
  // In the generator branch of this ntuple, first particle is always
  // negative, and second always positive. In the Dielectron block
  // of the ntuple, the first particle is always the one with larger Pt.
  double dR1=999, dR2=999;
  TLorentzVector v1reco, v2reco, v1gen, v2gen;
  v1reco.SetPtEtaPhiM(dielectron->pt_1, dielectron->eta_1, dielectron->phi_1, 0.000511);
  v2reco.SetPtEtaPhiM(dielectron->pt_2, dielectron->eta_2, dielectron->phi_2, 0.000511);
  v1gen .SetPtEtaPhiM(gen->pt_1, gen->eta_1, gen->phi_1, 0.000511);
  v2gen .SetPtEtaPhiM(gen->pt_2, gen->eta_2, gen->phi_2, 0.000511);
  if( dielectron->q_1 < 0 ){
    dR1 = v1reco.DeltaR(v1gen);
    dR2 = v2reco.DeltaR(v2gen);
  }else{
    dR1 = v1reco.DeltaR(v2gen);
    dR2 = v2reco.DeltaR(v1gen);
  }
  // Require that both are within loose dR of 0.4, otherwise bail out
  if( fabs(dR1) > 0.4 || fabs(dR2) > 0.4 ) result = kFALSE; 
  
  return result;
}

double findEventScaleFactor(const TElectron *leading, const TElectron *trailing){

  double result = findRecoScaleFactor(leading)
    * findRecoScaleFactor(trailing)
    * findIdScaleFactor(leading)
    * findIdScaleFactor(trailing)
    * findHltScaleFactor(leading)
    * findHltScaleFactor(trailing);

  return result;
}

double findRecoScaleFactor(const TElectron *ele){

  double result = 0;

  // Assume 2 eta bins
  int etBin = findEtBin(ele->scEt, etBinning);

  if( etBin == -1){
    // Found bin outside of calibrated range, return 1.0
    result = 1.0;
    return result;
  }

  if( fabs(ele->scEta) < 1.479 ){
    // Barrel
    result = GsfBarrelDataEff[etBin] / GsfBarrelMcEff[etBin];
  }else{
    // Barrel
    result = GsfEndcapDataEff[etBin] / GsfEndcapMcEff[etBin];
  }

  return result;
}

double findIdScaleFactor  (const TElectron *ele){
 
  double result = 0;

  // Assume 2 eta bins
  int etBin = findEtBin(ele->scEt, etBinning);

  if( etBin == -1){
    // Found bin outside of calibrated range, return 1.0
    result = 1.0;
    return result;
  }

  if( fabs(ele->scEta) < 1.479 ){
    // Barrel
    result = IdBarrelDataEff[etBin] / IdBarrelMcEff[etBin];
  }else{
    // Barrel
    result = IdEndcapDataEff[etBin] / IdEndcapMcEff[etBin];
  }

  return result;
}

double findHltScaleFactor (const TElectron *ele){

  double result = 0;

  // Assume 2 eta bins
  int etBin = findEtBin(ele->scEt, etBinning);

  if( etBin == -1){
    // Found bin outside of calibrated range, return 1.0
    result = 1.0;
    return result;
  }

  if( fabs(ele->scEta) < 1.479 ){
    // Barrel
    result = HltBarrelDataEff[etBin] / HltBarrelMcEff[etBin];
  }else{
    // Barrel
    result = HltEndcapDataEff[etBin] / HltEndcapMcEff[etBin];
  }

  return result;
}

// ---------------------- all scale factors smeared ------------------------------

double findEventScaleFactorSmeared(const TElectron *leading, const TElectron *trailing, int iexp){

  double result = findRecoScaleFactorSmeared(leading,iexp)
    * findRecoScaleFactorSmeared(trailing,iexp)
    * findIdScaleFactorSmeared(leading,iexp)
    * findIdScaleFactorSmeared(trailing,iexp)
    * findHltScaleFactorSmeared(leading,iexp)
    * findHltScaleFactorSmeared(trailing,iexp);

  return result;
}

double findRecoScaleFactorSmeared(const TElectron *ele, int iexp){

  double result = 0;

  // Assume 2 eta bins
  int etBin = findEtBin(ele->scEt, etBinning);

  if( etBin == -1){
    // Found bin outside of calibrated range, return 1.0
    result = 1.0;
    return result;
  }

  double effData = 100, effMC = 100;
  if( fabs(ele->scEta) < 1.479 ){
    // Barrel
    effData = GsfBarrelDataEff[etBin] + ro_D_B_reco[iexp] * GsfBarrelDataEffErr[etBin]; 
    effMC   = GsfBarrelMcEff[etBin]   + ro_M_B_reco[iexp] * GsfBarrelMcEffErr  [etBin];
  }else{
    // Barrel
    effData = GsfEndcapDataEff[etBin] + ro_D_E_reco[iexp] * GsfEndcapDataEffErr[etBin];
    effMC   = GsfEndcapMcEff  [etBin] + ro_M_E_reco[iexp] * GsfEndcapMcEffErr  [etBin];
  }

  // In several cases for this efficiency, we have sometthing like 100% +0 -10%,
  // this is not quite correct, but quick to implement, method to fix it.
  if( effData > 100 ) 
    effData = 100;

  result = effData / effMC;
  return result;
}

double findIdScaleFactorSmeared  (const TElectron *ele, int iexp){
 
  double result = 0;

  // Assume 2 eta bins
  int etBin = findEtBin(ele->scEt, etBinning);

  if( etBin == -1){
    // Found bin outside of calibrated range, return 1.0
    result = 1.0;
    return result;
  }

  if( fabs(ele->scEta) < 1.479 ){
    // Barrel
    result = 
      (  IdBarrelDataEff[etBin] + ro_D_B_id[iexp] * IdBarrelDataEffErr[etBin])
      / ( IdBarrelMcEff [etBin] + ro_M_B_id[iexp] * IdBarrelMcEffErr  [etBin]);
  }else{
    // Barrel
    result = 
      (   IdEndcapDataEff[etBin] + ro_D_E_id[iexp] * IdEndcapDataEffErr[etBin])
      / ( IdEndcapMcEff  [etBin] + ro_M_E_id[iexp] * IdEndcapMcEffErr  [etBin]);
  }

  return result;
}

double findHltScaleFactorSmeared (const TElectron *ele, int iexp){

  double result = 0;

  // Assume 2 eta bins
  int etBin = findEtBin(ele->scEt, etBinning);

  if( etBin == -1){
    // Found bin outside of calibrated range, return 1.0
    result = 1.0;
    return result;
  }

  if( fabs(ele->scEta) < 1.479 ){
    // Barrel
    result = 
      (   HltBarrelDataEff[etBin] + ro_D_B_hlt[iexp] * HltBarrelDataEffErr[etBin])
      / ( HltBarrelMcEff  [etBin] + ro_M_B_hlt[iexp] * HltBarrelMcEffErr  [etBin]);
  }else{
    // Barrel
    result = 
      (   HltEndcapDataEff[etBin] + ro_D_E_hlt[iexp] * HltEndcapDataEffErr[etBin])
      / ( HltEndcapMcEff  [etBin] + ro_M_E_hlt[iexp] * HltEndcapMcEffErr  [etBin]);
  }

  return result;
}


 void drawEfficiencyGraphs(TGraphErrors *grData, TGraphErrors *grMc,
			   TString yAxisTitle, TString text, TString plotName){
   
   // Generate "random" canvas name
//    TTimeStamp time;
//    TString c = "c";
//    c += time.GetNanoSec();
   TString c = plotName;
//    printf("Canvas name %s\n", c.Data());

   TCanvas *c2 = MakeCanvas(c,c);
   CPlot plot1(c,"","E_{T} [GeV]", yAxisTitle);
   plot1.SetLogx(); 
   plot1.AddGraph(grData,"data","PE2", kBlue);
   plot1.AddGraph(grMc  ,"MC"  ,"PE", kBlack);
   plot1.Draw(c2);
   plot1.SetYRange(0.2,1.1);
   grData->GetXaxis()->SetTitle("E_{T} [GeV]");
   grData->GetXaxis()->SetMoreLogLabels();
   grData->GetXaxis()->SetNoExponent();
   grData->SetFillStyle(3001);
   grData->SetFillColor(kBlue);
   grMc->SetMarkerStyle(24);
   plot1.TransLegend(0.0, -0.5);
   plot1.AddTextBox(text, 0.6,0.4,0.8,0.5, 0);
   plot1.Draw(c2, savePlots, "png");

   TLine *line = new TLine(10,1.0,500,1.0);
   line->SetLineStyle(kDashed);
   line->Draw("same");

   return;
 }

void drawEfficiencyGraphsAsymmErrors(TGraphAsymmErrors *grData, TGraphErrors *grMc,
				     TString yAxisTitle, TString text, TString plotName){
  
   // Generate "random" canvas name
//    TTimeStamp time;
//    TString c = "c";
//    c += time.GetNanoSec();
   TString c = plotName;
//    printf("Canvas name %s\n", c.Data());

   TCanvas *c2 = MakeCanvas(c,c);
   CPlot plot1(c,"","E_{T} [GeV]", yAxisTitle);
   plot1.SetLogx(); 
   plot1.AddGraph(grData,"data","PE2", kBlue);
   plot1.AddGraph(grMc  ,"MC"  ,"PE", kBlack);
   plot1.Draw(c2);
   plot1.SetYRange(0.2,1.1);
   grData->GetXaxis()->SetTitle("E_{T} [GeV]");
   grData->GetXaxis()->SetMoreLogLabels();
   grData->GetXaxis()->SetNoExponent();
   grData->SetFillStyle(3001);
   grData->SetFillColor(kBlue);
   grMc->SetMarkerStyle(24);
   plot1.TransLegend(0.0, -0.5);
   plot1.AddTextBox(text, 0.6,0.4,0.8,0.5, 0);
   plot1.Draw(c2, savePlots, "png");

   TLine *line = new TLine(10,1.0,500,1.0);
   line->SetLineStyle(kDashed);
   line->Draw("same");

   return;
 }

void drawScaleFactorGraphs(TGraphErrors *gr, TString yAxisTitle, TString text,
			   TString plotName){
   
   // Generate "random" canvas name
//    TTimeStamp time;
//    TString c = "c";
//    c += time.GetNanoSec();
  TString c = plotName;
//   printf("Canvas name %s\n", c.Data());

   TCanvas *c2 = MakeCanvas(c,c);
   CPlot plot1(c,"","E_{T} [GeV]", yAxisTitle);
   plot1.SetLogx(); 
   plot1.AddGraph(gr,"present E_{T}/#eta dep.","PE", kBlack);
   plot1.Draw(c2);
   plot1.SetYRange(0.5,1.5);
   gr->GetXaxis()->SetTitle("E_{T} [GeV]");
   gr->GetXaxis()->SetMoreLogLabels();
   gr->GetXaxis()->SetNoExponent();
   plot1.AddTextBox(text, 0.6,0.4,0.8,0.5, 0);
   plot1.Draw(c2,savePlots, "png");
   
   TLine *line = new TLine(10,1.0,500,1.0);
   line->SetLineStyle(kDashed);
   line->Draw("same");

   return;
 }

void drawEventScaleFactorGraphs(TGraphErrors *gr, TString yAxisTitle, 
				TString plotName){
  
  // Generate "random" canvas name
//   TTimeStamp time;
//   TString c = "c";
//   c += time.GetNanoSec();
  TString c = plotName;
//   printf("Canvas name %s\n", c.Data());
  
  TCanvas *c2 = MakeCanvas(c,c);
  CPlot plot1(c,"","m(e^{+}e^{-}) [GeV]", yAxisTitle);
  plot1.SetLogx(); 
  plot1.AddGraph(gr,"present E_{T}/#eta dep.","PE", kBlack);
  plot1.Draw(c2);
  plot1.SetYRange(0.0,1.5);
  gr->GetXaxis()->SetTitle("m(e^{+}e^{-}) [GeV]");
  gr->GetXaxis()->SetMoreLogLabels();
  gr->GetXaxis()->SetNoExponent();
//   cout << "From main progam CPlot::sOutDir " << CPlot::sOutDir << endl;
  CPlot::sOutDir = "plots";
  plot1.Draw(c2, savePlots, "png");
  
  TLine *line = new TLine(15,1.0,500,1.0);
  line->SetLineStyle(kDashed);
  line->Draw("same");

  return;
 }

void drawEfficiencies(){

  // Make grpahs
  double x[etBinCount];
  double dx[etBinCount];
  for(int i=0; i<etBinCount; i++){
    x[i]  = 0.5*(etBinLimits[i] + etBinLimits[i+1]);
    dx[i] = 0.5*(etBinLimits[i+1] - etBinLimits[i]);
  }

  // Graphs for SC->GSF efficiency
  TGraphErrors *grGsfBarrelDataEff 
//     = new TGraphAsymmErrors(etBinCount, x,  GsfBarrelDataEff, dx, dx, 
// 		       GsfBarrelDataEffErr, GsfBarrelDataEffErrPos);
    = new TGraphErrors(etBinCount, x,  GsfBarrelDataEff, dx, GsfBarrelDataEffErr);
  
  TGraphErrors *grGsfEndcapDataEff 
//     = new TGraphAsymmErrors(etBinCount, x,  GsfEndcapDataEff, dx, dx,
// 		       GsfEndcapDataEffErr, GsfEndcapDataEffErrPos);
    = new TGraphErrors(etBinCount, x,  GsfEndcapDataEff, dx, GsfEndcapDataEffErr);
  
  TGraphErrors *grGsfBarrelMcEff 
    = new TGraphErrors(etBinCount, x,  GsfBarrelMcEff, dx, GsfBarrelMcEffErr);
  
  TGraphErrors *grGsfEndcapMcEff 
    = new TGraphErrors(etBinCount, x,  GsfEndcapMcEff, dx, GsfEndcapMcEffErr);
  
  // Graphs for GSF->IDed efficiency
  TGraphErrors *grIdBarrelDataEff 
    = new TGraphErrors(etBinCount, x,  IdBarrelDataEff, dx, IdBarrelDataEffErr);
  
  TGraphErrors *grIdEndcapDataEff 
    = new TGraphErrors(etBinCount, x,  IdEndcapDataEff, dx, IdEndcapDataEffErr);
  
  TGraphErrors *grIdBarrelMcEff 
    = new TGraphErrors(etBinCount, x,  IdBarrelMcEff, dx, IdBarrelMcEffErr);
  
  TGraphErrors *grIdEndcapMcEff 
    = new TGraphErrors(etBinCount, x,  IdEndcapMcEff, dx, IdEndcapMcEffErr);
  
  // Graphs for IDed->HLT efficiency
  TGraphErrors *grHltBarrelDataEff 
    = new TGraphErrors(etBinCount, x,  HltBarrelDataEff, dx, HltBarrelDataEffErr);
  
  TGraphErrors *grHltEndcapDataEff 
    = new TGraphErrors(etBinCount, x,  HltEndcapDataEff, dx, HltEndcapDataEffErr);
  
  TGraphErrors *grHltBarrelMcEff 
    = new TGraphErrors(etBinCount, x,  HltBarrelMcEff, dx, HltBarrelMcEffErr);
  
  TGraphErrors *grHltEndcapMcEff 
    = new TGraphErrors(etBinCount, x,  HltEndcapMcEff, dx, HltEndcapMcEffErr);
  
  // Draw all graphs
  TString plotName;
  plotName = "plot_eff_gsf_barrel";
  drawEfficiencyGraphs(grGsfBarrelDataEff, grGsfBarrelMcEff,
 		       "efficiency_{RECO}", "Barrel",plotName);

  plotName = "plot_eff_gsf_endcap";
  drawEfficiencyGraphs(grGsfEndcapDataEff, grGsfEndcapMcEff,
 		       "efficiency_{RECO}", "Endcap",plotName);

  plotName = "plot_eff_id_barrel";
  drawEfficiencyGraphs(grIdBarrelDataEff, grIdBarrelMcEff,
 		       "efficiency_{ID}", "Barrel",plotName);

  plotName = "plot_eff_id_endcap";
  drawEfficiencyGraphs(grIdEndcapDataEff, grIdEndcapMcEff,
 		       "efficiency_{ID}", "Endcap",plotName);

  plotName = "plot_eff_hlt_barrel";
  drawEfficiencyGraphs(grHltBarrelDataEff, grHltBarrelMcEff,
 		       "efficiency_{HLT}", "Barrel",plotName);

  plotName = "plot_eff_hlt_endcap";
  drawEfficiencyGraphs(grHltEndcapDataEff, grHltEndcapMcEff,
 		       "efficiency_{HLT}", "Endcap",plotName);

  return;
}

void drawScaleFactors(){

  double x[etBinCount];
  double dx[etBinCount];
  double scaleGsfBarrel   [etBinCount];
  double scaleGsfBarrelErr[etBinCount];
  double scaleGsfEndcap   [etBinCount];
  double scaleGsfEndcapErr[etBinCount];
  double scaleIdBarrel    [etBinCount];
  double scaleIdBarrelErr [etBinCount];
  double scaleIdEndcap    [etBinCount];
  double scaleIdEndcapErr [etBinCount];
  double scaleHltBarrel   [etBinCount];
  double scaleHltBarrelErr[etBinCount];
  double scaleHltEndcap   [etBinCount];
  double scaleHltEndcapErr[etBinCount];
  for(int i=0; i<etBinCount; i++){
    x[i]  = (etBinLimits[i] + etBinLimits[i+1])/2.0;
    dx[i] = (etBinLimits[i+1] - etBinLimits[i])/2.0;
    scaleGsfBarrel   [i] = GsfBarrelDataEff[i] / GsfBarrelMcEff[i];
    scaleGsfBarrelErr[i] = errOnRatio( GsfBarrelDataEff[i], GsfBarrelDataEffErr[i], 
				       GsfBarrelMcEff[i] , GsfBarrelMcEffErr[i]);

    scaleGsfEndcap   [i] = GsfEndcapDataEff[i] / GsfEndcapMcEff[i];
    scaleGsfEndcapErr[i] = errOnRatio( GsfEndcapDataEff[i], GsfEndcapDataEffErr[i], 
				       GsfEndcapMcEff[i] , GsfEndcapMcEffErr[i]);

    scaleIdBarrel    [i] = IdBarrelDataEff[i] / IdBarrelMcEff[i];
    scaleIdBarrelErr [i] = errOnRatio( IdBarrelDataEff[i], IdBarrelDataEffErr[i], 
				       IdBarrelMcEff[i] , IdBarrelMcEffErr[i]);

    scaleIdEndcap    [i] = IdEndcapDataEff[i] / IdEndcapMcEff[i];
    scaleIdEndcapErr [i] = errOnRatio( IdEndcapDataEff[i], IdEndcapDataEffErr[i], 
				       IdEndcapMcEff[i] , IdEndcapMcEffErr[i]);

    scaleHltBarrel   [i] = HltBarrelDataEff[i] / HltBarrelMcEff[i];
    scaleHltBarrelErr[i] = errOnRatio( HltBarrelDataEff[i], HltBarrelDataEffErr[i], 
				       HltBarrelMcEff[i] , HltBarrelMcEffErr[i]);

    scaleHltEndcap   [i] = HltEndcapDataEff[i] / HltEndcapMcEff[i];
    scaleHltEndcapErr[i] = errOnRatio( HltEndcapDataEff[i], HltEndcapDataEffErr[i], 
				       HltEndcapMcEff[i] , HltEndcapMcEffErr[i]);
  }
//   // Zero out HLT scale factor below 20 GeV (i.e. in the very first bin)
//   scaleHltBarrel   [0] = 0;
//   scaleHltBarrelErr[0] = 0;
//   scaleHltEndcap   [0] = 0;
//   scaleHltEndcapErr[0] = 0;

  TGraphErrors *grGsfBarrel 
    = new TGraphErrors(etBinCount, x, scaleGsfBarrel, dx, scaleGsfBarrelErr);

  TGraphErrors *grGsfEndcap 
    = new TGraphErrors(etBinCount, x, scaleGsfEndcap, dx, scaleGsfEndcapErr);

  TGraphErrors *grIdBarrel 
    = new TGraphErrors(etBinCount, x, scaleIdBarrel, dx, scaleIdBarrelErr);

  TGraphErrors *grIdEndcap 
    = new TGraphErrors(etBinCount, x, scaleIdEndcap, dx, scaleIdEndcapErr);

  TGraphErrors *grHltBarrel 
    = new TGraphErrors(etBinCount, x, scaleHltBarrel, dx, scaleHltBarrelErr);

  TGraphErrors *grHltEndcap 
    = new TGraphErrors(etBinCount, x, scaleHltEndcap, dx, scaleHltEndcapErr);

  TString plotName;
  plotName = "plot_scale_gsf_barrel";
  drawScaleFactorGraphs(grGsfBarrel, "scale factor RECO", "Barrel", plotName);
  plotName = "plot_scale_gsf_endcap";
  drawScaleFactorGraphs(grGsfEndcap, "scale factor RECO", "Endcap", plotName);
  plotName = "plot_scale_id_barrel";
  drawScaleFactorGraphs(grIdBarrel , "scale factor ID"  , "Barrel", plotName);
  plotName = "plot_scale_id_endcap";
  drawScaleFactorGraphs(grIdEndcap , "scale factor ID"  , "Endcap", plotName);
  plotName = "plot_scale_hlt_barrel";
  drawScaleFactorGraphs(grHltBarrel, "scale factor HLT" , "Barrel", plotName);
  plotName = "plot_scale_hlt_endcap";
  drawScaleFactorGraphs(grHltEndcap, "scale factor HLT" , "Endcap", plotName);

}

double errOnRatio(double a, double da, double b, double db){

  double result = 0;
  if(a == 0 || b == 0)
    return result;
 
  result = (a/b)*sqrt( (da/a)*(da/a) + (db/b)*(db/b) );
  return result;
}

void drawEventScaleFactors(TVectorD scaleGsfV, TVectorD scaleGsfErrV,
			   TVectorD scaleIdV , TVectorD scaleIdErrV ,
			   TVectorD scaleHltV, TVectorD scaleHltErrV,
			   TVectorD scaleV   , TVectorD scaleErrV    )
{

  // repackage into arrays
  double x[nMassBins];
  double dx[nMassBins];
  double scaleGsfA   [nMassBins];
  double scaleIdA    [nMassBins];
  double scaleHltA   [nMassBins];
  double scaleA      [nMassBins];
  double scaleGsfErrA[nMassBins];
  double scaleIdErrA [nMassBins];
  double scaleHltErrA[nMassBins];
  double scaleErrA   [nMassBins];
  for(int i=0; i<nMassBins; i++){
    x[i] = (massBinLimits[i] + massBinLimits[i+1])/2.0;
    dx[i]= (massBinLimits[i+1] - massBinLimits[i])/2.0;
    scaleGsfA       [i] = scaleGsfV   [i];
    scaleIdA        [i] = scaleIdV    [i];
    scaleHltA       [i] = scaleHltV   [i];
    scaleA          [i] = scaleV      [i];
    scaleGsfErrA    [i] = scaleGsfErrV[i];
    scaleIdErrA     [i] = scaleIdErrV [i];
    scaleHltErrA    [i] = scaleHltErrV[i];
    scaleErrA       [i] = scaleErrV   [i];
  }

  TGraphErrors *grScale = new TGraphErrors(nMassBins, x, scaleA, dx, scaleErrA);
  TGraphErrors *grScaleGsf = new TGraphErrors(nMassBins, x, scaleGsfA, dx, scaleGsfErrA);
  TGraphErrors *grScaleId  = new TGraphErrors(nMassBins, x, scaleIdA , dx, scaleIdErrA );
  TGraphErrors *grScaleHlt = new TGraphErrors(nMassBins, x, scaleHltA, dx, scaleHltErrA);

  TString plotName;
  plotName = "plot_event_scale_gsf";
  drawEventScaleFactorGraphs(grScaleGsf, "RECO scale factor" , plotName);
  plotName = "plot_event_scale_id";
  drawEventScaleFactorGraphs(grScaleId , "ID scale factor"   , plotName);
  plotName = "plot_event_scale_hlt";
  drawEventScaleFactorGraphs(grScaleHlt, "HLT scale factor"  , plotName);
  plotName = "plot_event_scale_full";
  drawEventScaleFactorGraphs(grScale   , "event scale factor", plotName);

}

// This method reads all ROOT files that have efficiencies from
// tag and probe in TMatrixD form and converts the matrices into 
// more simple arrays.
void fillEfficiencyConstants(  const TriggerSelection &triggers ) {
  /*
  TString effDataGsfFile = "efficiency_TnP_data_gsf_fit-fit_bins-et5-eta2.root";
  TString effMcGsfFile   = "efficiency_TnP_mc_gsf_count-count_bins-et5-eta2.root";

  TString effDataIdFile  = "efficiency_TnP_data_id_fit-fit_bins-et5-eta2.root";
  TString effMcIdFile    = "efficiency_TnP_mc_id_count-count_bins-et5-eta2.root";

  TString effDataHltFile = "efficiency_TnP_data_hlt_count-count_bins-et5-eta2.root";
  TString effMcHltFile   = "efficiency_TnP_mc_hlt_count-count_bins-et5-eta2.root";
  */
  TString fnStart="efficiency_TnP_";
  TString fnEnd=".root";
  TString effDataGsfFile = fnStart + getLabel(DATA,GSF,dataGsfEffMethod,etBinning,etaBinning,triggers) + fnEnd;
  TString effMcGsfFile   = fnStart + getLabel(MC  ,GSF,  mcGsfEffMethod,etBinning,etaBinning,triggers) + fnEnd;
  TString effDataIdFile  = fnStart + getLabel(DATA, ID, dataIdEffMethod,etBinning,etaBinning,triggers) + fnEnd;
  TString effMcIdFile    = fnStart + getLabel(MC  , ID,   mcIdEffMethod,etBinning,etaBinning,triggers) + fnEnd;
  TString effDataHltFile = fnStart + getLabel(DATA,HLT,dataHltEffMethod,etBinning,etaBinning,triggers) + fnEnd;
  TString effMcHltFile   = fnStart + getLabel(MC  ,HLT,  mcHltEffMethod,etBinning,etaBinning,triggers) + fnEnd;

  // Continue assuming 2 eta bins
  // Last parameter is 0=barrel, 1=endcap
  fillOneEfficiency(effDataGsfFile, GsfBarrelDataEff, GsfBarrelDataEffErr, 0);
  fillOneEfficiency(effDataGsfFile, GsfEndcapDataEff, GsfEndcapDataEffErr, 1);
  fillOneEfficiency(effMcGsfFile, GsfBarrelMcEff, GsfBarrelMcEffErr, 0);
  fillOneEfficiency(effMcGsfFile, GsfEndcapMcEff, GsfEndcapMcEffErr, 1);

  fillOneEfficiency(effDataIdFile, IdBarrelDataEff, IdBarrelDataEffErr, 0);
  fillOneEfficiency(effDataIdFile, IdEndcapDataEff, IdEndcapDataEffErr, 1);
  fillOneEfficiency(effMcIdFile, IdBarrelMcEff, IdBarrelMcEffErr, 0);
  fillOneEfficiency(effMcIdFile, IdEndcapMcEff, IdEndcapMcEffErr, 1);

  fillOneEfficiency(effMcHltFile, HltBarrelMcEff, HltBarrelMcEffErr, 0);
  fillOneEfficiency(effMcHltFile, HltEndcapMcEff, HltEndcapMcEffErr, 1);

  if ( ! triggers.hltEffMethodIs2011New() ) {
    fillOneEfficiency(effDataHltFile, HltBarrelDataEff, HltBarrelDataEffErr, 0);
    fillOneEfficiency(effDataHltFile, HltEndcapDataEff, HltEndcapDataEffErr, 1);
  }
  else {
    std::cout << " loading files for luminosity reweighting for data HLT-efficiency\n";
    double barrelEff[etBinCount], barrelEffErr[etBinCount];
    double endcapEff[etBinCount], endcapEffErr[etBinCount];
    TriggerSelection locTrig(triggers);

    // clear the array
    for (int i=0; i<etBinCount; ++i) HltBarrelDataEff[i]=0.;
    for (int i=0; i<etBinCount; ++i) HltBarrelDataEffErr[i]=0.;
    // iterate over run chunks
    const double totLumi=runLumiV[0]+runLumiV[1]+runLumiV[2];
    const double invTotLumi=1/totLumi;
    for (int lumiIdx=0; lumiIdx<3; ++lumiIdx) {
      switch(lumiIdx) {
      case 0: locTrig.triggerSet(TrigSet_2011A_SingleEG); break;
      case 1: locTrig.triggerSet(TrigSet_2011A_DoubleEG); break;
      case 2: locTrig.triggerSet(TrigSet_2011B_DoubleEG); break;
      default:
	std::cout << "error in the code\n";
	assert(0);
      }
      effDataHltFile = fnStart + getLabel(DATA,HLT,dataHltEffMethod,etBinning,etaBinning,locTrig) + fnEnd;
      fillOneEfficiency(effDataHltFile, barrelEff, barrelEffErr, 0);
      for (int i=0; i<etBinCount; ++i) HltBarrelDataEff[i]    += invTotLumi * runLumiV[lumiIdx] * barrelEff[i];
      for (int i=0; i<etBinCount; ++i) HltBarrelDataEffErr[i] += invTotLumi * runLumiV[lumiIdx] * barrelEffErr[i];
      fillOneEfficiency(effDataHltFile, endcapEff, endcapEffErr, 1);
      for (int i=0; i<etBinCount; ++i) HltEndcapDataEff[i]    += invTotLumi * runLumiV[lumiIdx] * endcapEff[i];
      for (int i=0; i<etBinCount; ++i) HltEndcapDataEffErr[i] += invTotLumi * runLumiV[lumiIdx] * endcapEffErr[i];
    }
  }

}

void fillOneEfficiency(const TString filename, double *eff, double *effErr, int etaRange){

  TFile f(TString("../root_files/tag_and_probe/")+dirTag+TString("/")+filename);
  if(!f.IsOpen()) assert(0);

  TMatrixD *effMatrix        = (TMatrixD*)f.Get("effArray2D");
  TMatrixD *effMatrixErrLow  = (TMatrixD*)f.Get("effArrayErrLow2D");
  TMatrixD *effMatrixErrHigh = (TMatrixD*)f.Get("effArrayErrHigh2D");

  // Make sure that the objects are present
  if( !(effMatrix && effMatrixErrLow && effMatrixErrHigh) ) assert(0);

  // Make sure that there are only two eta bins and appropriate number of ET bins
  if( effMatrix->GetNcols() != 2 ) {
    printf("The number of eta bins stored in constants files is not 2, crashing\n");
    assert(0);
  }
  if( effMatrix->GetNrows() != getNEtBins(etBinning) ) {
    printf("The number of ET bins stored in constants files is different form expeted, crashing\n");
    printf(" Matrix %d    expect %d\n", effMatrix->GetNrows(), getNEtBins(etBinning));
    assert(0);
  }

  for(int i=0; i<getNEtBins(etBinning); i++){
    eff[i] = (*effMatrix)(i,etaRange);
    // For the errors, take for now the average between asymmetric errors
    // Really, in most cases the errors are identical. 
    // Handling asymmetric errors is rather non-trivial.
    effErr[i] = ( (*effMatrixErrLow)(i,etaRange) + (*effMatrixErrHigh)(i,etaRange))/2.0;
  }
  f.Close();

}
