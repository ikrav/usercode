
/*
   Selection to worry about for e-mu:
    - trigger bits for the whole event
       (depends on whether it is Mu or Electron datset?)
    - for elelectron:
       - which Et cut?
       - trigger matching or not?
    - for muon:
       - which Pt cut?
       - trigger matching or not?

   HOW WE DETERMINE IF THIS IS DATA OR MC?

To do's
1) Find the bin in the fake rate histo associated with the relevant Et (use lower_bound algo)
2) Bin value should be a key in a map where the value in the map is the rate
3) apply fake rate to the event and plot the mass


*/
//================================================================================================
//
// X->e mu selection macro
//
//  * plots distributions associated with selected events
//  * prints list of selected events from data
//  * outputs ROOT files of events passing selection foreach sample, 
//    which can be processed by plotSelect.C
//
//________________________________________________________________________________________________
 
//C-M-\ justify text in emacs

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TBranch.h>                // class to access branches in TTree
#include <TClonesArray.h>           // ROOT array class
#include <TCanvas.h>                // class for drawing
#include <TH1F.h>                   // 1D histograms
#include <TBenchmark.h>             // class to track macro running statistics
#include <TLorentzVector.h>         // 4-vector class
#include <TVector3.h>               // 3D vector class
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <fstream>                  // file stream
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings
#include <assert.h>
#include <boost/scoped_ptr.hpp>     // Use boost sole ownership pointer
#include <boost/shared_ptr.hpp>     // Use boost shared ownership pointer
#include <map>                      // use stl map
#include <set>                      // use stl set
#include <utility>                  // use make_pair for inserting pairs into a map 
#include <algorithm>                //stl algorithms

#include "../../Include/MitStyleRemix.hh"  // style settings for drawing

// define structures to read in ntuple
#include "../../Include/EWKAnaDefs.hh"
#include "../../Include/TEventInfo.hh"
#include "../../Include/TElectron.hh"
#include "../../Include/TPhoton.hh"
#include "../../Include/TVertex.hh"
#include "../../Include/TMuon.hh"

#include "TreeQueue.hh"

// lumi section selection with JSON files
// #include "RunLumiRangeMap.h"
#include "../../Include/JsonParser.hh"

// Helper functions for Electron ID selection
#include "../../Include/EleIDCuts.hh" 

#endif
//=== FUNCTION DECLARATIONS ======================================================================================

const Double_t kGAP_LOW  = 1.4442;
const Double_t kGAP_HIGH = 1.566;
const Double_t REGION_2 = 2.00;
const Double_t REGION_3 = 2.50;

const int nMassBins = 39;
const double massBinLimits[nMassBins+1] = 
  {15,20,25,30,35,40,45,50,55,60,64,72,76,81,86,91,96,101,106,110,115,120,126,133,141,150,160,171,185,200,220,243,273,320,380,440,510,600,1000,1500};
const int eTBins(15);
const double eTBinLimits[eTBins+1] =  {0,5,10,15,20,25,30,35,40,50,60,70,90,150,200,300};
 

// fill ntuple of selected events
/*void fillData(ZemuData *data, 
	      const mithep::TEventInfo *info, 
	      const double mass, const double pt,
	      const mithep::TMuon *muon, 
	      const mithep::TElectron *electron, 
              const UInt_t npv, const Double_t weight);
*/

void correctEnergyScale(mithep::TElectron *electron, mithep::TEventInfo *info);

void TH1FCopy(const TH1F* original, TH1F* copy);

Bool_t passSmurf_elec1(const mithep::TDielectron *dielectron);
Bool_t passSmurf_elec2(const mithep::TDielectron *dielectron);

//===Pull objects and functions from namespace===================================================================

using std::string;
using std::cout;
using std::fstream;
using std::map;
using std::make_pair;
using std::pair;
using std::set;
using std::for_each;
using boost::scoped_ptr;
using boost::shared_ptr;


//============================
//Simple class to hold two values
//Class value does not change if value order is changed
// ie a->b and b->a
//=============================

template <typename T>
class Couple{

public:
  Couple(const T& a, const T& b): first(a), second(b){}
  T alpha() const {return first;}
  T beta() const {return second;}
  bool operator==(const Couple& rhs) const;
  bool operator!=(const Couple& rhs) const;

private:
  T first;
  T second;

};

template <typename T> bool Couple<T>::operator==(const Couple& rhs) const{
  bool pass(false);
  if ((this->first == rhs.first) && (this->second == rhs.second)) {
    pass = true;
    return pass;//stop here then
  }
  if ((this->first == rhs.second) && (this->second == rhs.first)) pass = true;
  return pass;
}

template <typename T> bool Couple<T>::operator!=(const Couple& rhs) const{
  return (!(*this==rhs));
}

struct check_etDiff_functor{
  check_etDiff_functor(const Float_t& inputScET, bool& inputBool):  photonScEt(inputScET), testBool(&inputBool){}
  void operator()(pair<Float_t,UInt_t> scEtPair){
    if (fabs(scEtPair.first - photonScEt) < 0.02) {
      *testBool = true; //assuming it is the same cluster
    }
  }
  Float_t photonScEt;
  bool* testBool;
};

//=== MAIN MACRO =================================================================================================

void fakeRate() 
{  
  gBenchmark->Start("fakeRate");

  
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================
  
  TString  outputDir;         // output directory
  TString  trigName;          // which trigger to use
  TString  format;            // plot format

  //CPlot::sOutDir        = outputDir + TString("/plots");   gSystem->mkdir(CPlot::sOutDir,kTRUE);
  //const TString ntupDir = outputDir + TString("/ntuples"); gSystem->mkdir(ntupDir,kTRUE);

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
  
  //
  // Access samples and fill histograms
  //  
  TFile *infile=0;
  TTree *eventTree=0;  
  
  bool leadingTrigCut(false);

  // Data structures to store info from TTrees
  mithep::TEventInfo *info    = new mithep::TEventInfo();

  TClonesArray *electronArr = new TClonesArray("mithep::TElectron");
  TClonesArray *diElectronArr = new TClonesArray("mithep::TDielectron");
 
  //Set up histos to store events
  const int nBins(30),xmin(0), xmax(300);
  const int rap_Bins = 24;
  const int rap_Bins_highMass = 12;
  float rap_Min = 0;
  float rap_Max = 2.4;

  scoped_ptr<TH1F> cluster_et = scoped_ptr<TH1F>(new TH1F("cluster_et","",nBins,xmin,xmax));
  scoped_ptr<TH1F> cluster_et_reg = scoped_ptr<TH1F>(new TH1F("cluster_et_reg","",eTBins,eTBinLimits));
  scoped_ptr<TH1F> cluster_et_reg2 = scoped_ptr<TH1F>(new TH1F("cluster_et_reg2","",eTBins,eTBinLimits));
  scoped_ptr<TH1F> cluster_et_reg3 = scoped_ptr<TH1F>(new TH1F("cluster_et_reg3","",eTBins,eTBinLimits));
  scoped_ptr<TH1F> electron_et = scoped_ptr<TH1F>(new TH1F("electron_et","",nBins,xmin,xmax));
  scoped_ptr<TH1F> electron_et_reg = scoped_ptr<TH1F>(new TH1F("electron_et_reg","",eTBins,eTBinLimits));
  scoped_ptr<TH1F> electron_et_reg2 = scoped_ptr<TH1F>(new TH1F("electron_et_reg2","",eTBins,eTBinLimits));
  scoped_ptr<TH1F> electron_et_reg3 = scoped_ptr<TH1F>(new TH1F("electron_et_reg3","",eTBins,eTBinLimits));
  scoped_ptr<TH1F> electron_iso_et_reg = scoped_ptr<TH1F>(new TH1F("electron_iso_et_reg","",eTBins,eTBinLimits));
  scoped_ptr<TH1F> electron_iso_et_reg2 = scoped_ptr<TH1F>(new TH1F("electron_iso_et_reg2","",eTBins,eTBinLimits));
  scoped_ptr<TH1F> electron_iso_et_reg3 = scoped_ptr<TH1F>(new TH1F("electron_iso_et_reg3","",eTBins,eTBinLimits));
  scoped_ptr<TH1F> cand_electron_et = scoped_ptr<TH1F>(new TH1F("cand_electron_et","",nBins,xmin,xmax));
  scoped_ptr<TH1F> fRate = scoped_ptr<TH1F>(new TH1F("fRate","",nBins,xmin,xmax));
  shared_ptr<TH1F> fRate_reg(new TH1F("fRate_reg","",eTBins,eTBinLimits));
  shared_ptr<TH1F> fRate_reg2(new TH1F("fRate_reg2","",eTBins,eTBinLimits));
  shared_ptr<TH1F> fRate_reg3(new TH1F("fRate_reg3","",eTBins,eTBinLimits));
  scoped_ptr<TH1F> pfMet = scoped_ptr<TH1F>(new TH1F("pfMet","",50,xmin,500));
  scoped_ptr<TH1F> pfMetx = scoped_ptr<TH1F>(new TH1F("pfMetx","",50,xmin,500));
  scoped_ptr<TH1F> pfMety = scoped_ptr<TH1F>(new TH1F("pfMety","",50,xmin,500));
  scoped_ptr<TH1F> caloMet = scoped_ptr<TH1F>(new TH1F("caloMet","",50,xmin,500));
  scoped_ptr<TH1F> etDiff = scoped_ptr<TH1F>(new TH1F("etDiff","",500,0,0.05));
  scoped_ptr<TH1F> etaDiff = scoped_ptr<TH1F>(new TH1F("etaDiff","",500,0,0.01));
  scoped_ptr<TH1F> testFail = scoped_ptr<TH1F>(new TH1F("testFail","",nMassBins,massBinLimits));

  scoped_ptr<TH1F> testEGsfFMass = scoped_ptr<TH1F>(new TH1F("testEGsfFMass","",nMassBins,massBinLimits));
  scoped_ptr<TH1F> testEGsfF = scoped_ptr<TH1F>(new TH1F("testEGsfF","",nMassBins,massBinLimits));
  scoped_ptr<TH1F> rapMass20to30 = scoped_ptr<TH1F>(new TH1F("rapMass20to30","",rap_Bins, rap_Min, rap_Max));
  scoped_ptr<TH1F> rapMass30to45 = scoped_ptr<TH1F>(new TH1F("rapMass30to45","",rap_Bins, rap_Min, rap_Max));
  scoped_ptr<TH1F> rapMass45to60 = scoped_ptr<TH1F>(new TH1F("rapMass45to60","",rap_Bins, rap_Min, rap_Max));
  scoped_ptr<TH1F> rapMass60to120 = scoped_ptr<TH1F>(new TH1F("rapMass60to120","",rap_Bins, rap_Min, rap_Max));
  scoped_ptr<TH1F> rapMass120to200 = scoped_ptr<TH1F>(new TH1F("rapMass120to200","",rap_Bins, rap_Min, rap_Max));
  scoped_ptr<TH1F> rapMass200to1500 = scoped_ptr<TH1F>(new TH1F("rapMass200to1500","",rap_Bins_highMass, rap_Min, rap_Max));

  if (leadingTrigCut){
    fRate_reg = shared_ptr<TH1F>(new TH1F("fRate_reg","",eTBins,xmin,75));
    fRate_reg2 = shared_ptr<TH1F>(new TH1F("fRate_reg2","",eTBins,xmin,75));
    fRate_reg3 = shared_ptr<TH1F>(new TH1F("fRate_reg3","",eTBins,xmin,75));
  }

  //
  // Set up event dump to file
  //
  ofstream evtfile;
  unsigned int arraysize;
  {
    //create a large char string in a local scope
    //it will be removed from the stack after we know
    //how large it really needs to be
    char buffer[1000]; 
    arraysize = sprintf(buffer, "Output.txt");    
  }

  char evtfname[arraysize];
  sprintf(evtfname, "Output.txt");

  evtfile.open(evtfname);
  assert(evtfile.is_open());
  
  //
  // Set up output ntuple file for the sample
  //
  TString outName = TString("plots.root");
  TFile *outFile = new TFile(outName,"RECREATE");

  //
  // loop through files
  //
  TString ts_evtfname(evtfname);
  
  TreeQueue fakeTrees("FakeDateSource.txt");
  //Get the Tree
  //loop over files
  cout << "start of nextfile\n";
  while( fakeTrees.nextFile() ){

    eventTree = fakeTrees.getTree("Events");
    fakeTrees.status();

    // Set branch address to structures that will store the info  
    eventTree->SetBranchAddress("Info",       &info);          TBranch *infoBr       = eventTree->GetBranch("Info");
    eventTree->SetBranchAddress("Electron" ,  &electronArr  ); TBranch *electronBr   = eventTree->GetBranch("Electron");
  
    UInt_t maxEvents = eventTree->GetEntries();
     
    // loop through events

    for(UInt_t ientry=0; ientry<maxEvents; ++ientry) { 
      if(ientry >= maxEvents) break;
	
      //info->Clear(); should add this back
      infoBr->GetEntry(ientry);
		
      electronArr->Clear(); 
      electronBr->GetEntry(ientry);

      //max number of electrons and photons in event
      UInt_t maxElectronEvts = electronArr->GetEntriesFast();
      if (maxElectronEvts > 1) continue; 
 
      UInt_t photonTrigOr =  kHLT_Photon30_CaloIdVL | kHLT_Photon50_CaloIdVL 
      	| kHLT_Photon75_CaloIdVL | kHLT_Photon90_CaloIdVL | kHLT_Photon125 | kHLT_Photon135; 

      if(!(info->triggerBits & photonTrigOr)) continue;
     
      /* 
      //reject events with 2 or more electrons
      UInt_t nElec(0);
      for(UInt_t ielec=0; ielec < maxElectronEvts; ++ielec) {
	mithep::TElectron *elec = (mithep::TElectron*)((*electronArr)[ielec]);
	if (elec->scEt > 10) ++nElec; 
      }

      if (nElec > 1) continue; */

      double metVal = info->pfMET;
      pfMet->Fill(metVal);
      //cout << "pfMET is: " << info->pfSumET << "\n";
      if (metVal >= 10.0) continue; //reject events with pf MET < 10

      // loop through photons      
      for(UInt_t i=0; i<maxElectronEvts; ++i) {
	mithep::TElectron *elec = (mithep::TElectron*)((*electronArr)[i]);
	if(    elec->scEt < 10      )   continue; //ignore low energy electrons

        //make sure elec did not cause the trigger
	//if (elec->hltMatchBits & triggerEle) continue;
	if (leadingTrigCut){ // This is ignored
	  if (elec->hltMatchBits & photonTrigOr) continue;
	}

	// loop through electrons

	//is ecaldriven
	if ( !(elec->isEcalDriven) ) continue;
	if (elec->nExpHitsInner > 0 ) continue;
              
	if((fabs(elec->scEta)>kGAP_LOW) && (fabs(elec->scEta)<kGAP_HIGH)) continue; // outside eta range? Skip to next event...
        if( fabs(elec->scEta)>REGION_3 ) continue; // outside eta range? Skip to next event...
        
	//isolation cuts
	if ((elec->emIso03/elec->pt) > 0.2) continue;
        if ((elec->hadIso03/elec->pt) > 0.2) continue;
        if ((elec->trkIso03/elec->pt) > 0.2) continue;

        //cuts in regions of eta

	if (fabs(elec->scEta) < kGAP_LOW) {
          if( elec->sigiEtaiEta	> 0.011) continue;
          if (elec->HoverE > 0.10) continue;//H/E cut for photons 
	  if (fabs(elec->deltaPhiIn) > 0.15) continue;	 
	  if (fabs(elec->deltaEtaIn) > 0.01) continue;
	} else {
	  if( elec->sigiEtaiEta	> 0.031) continue;
          if (elec->HoverE > 0.075) continue;//H/E cut for photons 
          if (fabs(elec->deltaPhiIn) > 0.10) continue;	 
	  if (fabs(elec->deltaEtaIn) > 0.01) continue;
	}

        //apply cuts to electron associated with cluster
	//correctEnergyScale(elec, info);	  

	cluster_et->Fill(elec->scEt);//fill with clusters that pass jet cuts
	
	//calculate fake rate in 3 eta regions

	if (fabs(elec->scEta) < kGAP_LOW){
	  cluster_et_reg->Fill(elec->scEt);
	} else if ((fabs(elec->scEta) > kGAP_HIGH) && (fabs(elec->scEta) < REGION_2)){
	  cluster_et_reg2->Fill(elec->scEt);
	} else if ((fabs(elec->scEta) > REGION_2) && (fabs(elec->scEta) < REGION_3)){
	  cluster_et_reg3->Fill(elec->scEt);
	}
	
        //look at rate of just the isolation cut
	
        if (fabs(elec->scEta) < kGAP_LOW){
	  if (elec->trkIso03 < 5.0){
	    electron_iso_et_reg->Fill(elec->scEt);
	  }
	} else if ((fabs(elec->scEta) > kGAP_HIGH) && (fabs(elec->scEta) < REGION_2)){
	  if (elec->trkIso03 < 5.0){
	    electron_iso_et_reg2->Fill(elec->scEt);
	  }
	} else if ((fabs(elec->scEta) > REGION_2) && (fabs(elec->scEta) < REGION_3)){
	  if (elec->trkIso03 < 5.0){
	    electron_iso_et_reg3->Fill(elec->scEt);
	  }
	}
	
	//
	// Apply full electron cuts
	//
	if( ! passSmurf(elec)) continue; ;

	electron_et->Fill(elec->scEt);// are there matches to other clusters?

	//calculate fake rate in 3 eta regions

	if (fabs(elec->scEta) < kGAP_LOW){
	  electron_et_reg->Fill(elec->scEt);
	} else if ((fabs(elec->scEta) > kGAP_HIGH) && (fabs(elec->scEta) < REGION_2)){
	  electron_et_reg2->Fill(elec->scEt);
	} else if ((fabs(elec->scEta) > REGION_2) && (fabs(elec->scEta) < REGION_3)){
	  electron_et_reg3->Fill(elec->scEt);
	}

      } // end loop over photon      
    } // end loop over events
  }//end of loop over files


  //Make a fake rate histo
  //Divide electrons by clusters(Jets)

  TH1FCopy(electron_et.get(), fRate.get());
  fRate->Divide(cluster_et.get());
  TH1FCopy(electron_et_reg.get(), fRate_reg.get());
  fRate_reg->Divide(cluster_et_reg.get());
  TH1FCopy(electron_et_reg2.get(), fRate_reg2.get());
  fRate_reg2->Divide(cluster_et_reg2.get());
  TH1FCopy(electron_et_reg3.get(), fRate_reg3.get());
  fRate_reg3->Divide(cluster_et_reg3.get()); 

  //Put fake rates into a map. Change to use boost unordered map for better performance, actually may not be able to do this as there is no upper or lower bound method for unordered_maps

  map<double,double> rate;
  map<double,double> rate2;
  map<double,double> rate3;

  // Put bin and fake rate values in a map
  
  //  double deltaBin = (xmax - xmin)/nBins;

  for (int i = 1; i < eTBins+1; ++i){
    rate.insert(make_pair(fRate_reg->GetBinLowEdge(i+1),fRate_reg->GetBinContent(i)));//eta region 1
    rate2.insert(make_pair(fRate_reg2->GetBinLowEdge(i+1),fRate_reg2->GetBinContent(i)));// eta region 2
    rate3.insert(make_pair(fRate_reg3->GetBinLowEdge(i+1),fRate_reg3->GetBinContent(i)));// eta region 3 
  }

  //Put 3 eta regions maps into a map
  map<double, map<double,double> > fRate_map;
  fRate_map.insert(make_pair(kGAP_LOW,rate));
  fRate_map.insert(make_pair(REGION_2,rate2));
  fRate_map.insert(make_pair(REGION_3,rate3));
  
  /*  
  Fake rate can be accessed in three steps
  1> Find the correct eta map for a particular value of eta
  2> In this map find the pair with the fake rate for the particular value of et
  3> access the fake rate in that pair

  Test example, find fake rate for eta of 1.2 and et of 20

  map<double, map<double,double> >::iterator rate_iter;
  map<double,double>::iterator rate_iter2;
  double fRateVal;

  rate_iter = fRate_map.upper_bound(1.2); (step 1)
  rate_iter2  = rate_iter->second.upper_bound(20); (step 2)
  fRateVal = rate_iter2->second;    (step 3)       

  Can write the above three lines in one.
  double fRateVal = fRate_map.upper_bound(1.2)->second.upper_bound(20)->second;
  avoids creating temporaries but not easy to decipher, hence the explanation */

  // look at fakes
  
   vector<TString> jsonList;
   
   /*
   jsonList.push_back("/tmp/olaiya/EE/Cert_160404-163869_7TeV_May10ReReco_Collisions11_JSON_v3.txt");
   jsonList.push_back("/tmp/olaiya/EE/Cert_160404-180252_7TeV_PromptReco_Collisions11_JSON.txt");
   jsonList.push_back("/tmp/olaiya/EE/Cert_160404-173692_7TeV_PromptReco_Collisions11_JSON.txt");
   jsonList.push_back("/tmp/olaiya/EE/Cert_170249-172619_7TeV_ReReco5Aug_Collisions11_JSON_v3.txt");
   jsonList.push_back("/tmp/olaiya/EE/Cert_160404-180252_7TeV_PromptReco_Collisions11_JSON.txt");
   */

   jsonList.push_back("/tmp/olaiya/EE/Cert_160404-180252_7TeV_PromptReco_Collisions11_JSON.txt");
   jsonList.push_back("/tmp/olaiya/EE/Cert_160404-180252_7TeV_PromptReco_Collisions11_JSON.txt");
   jsonList.push_back("/tmp/olaiya/EE/Cert_160404-180252_7TeV_PromptReco_Collisions11_JSON.txt");
   jsonList.push_back("/tmp/olaiya/EE/Cert_160404-180252_7TeV_PromptReco_Collisions11_JSON.txt");
   jsonList.push_back("/tmp/olaiya/EE/Cert_160404-180252_7TeV_PromptReco_Collisions11_JSON.txt");

   unsigned int jsonListSize = jsonList.size();

  //start from first file again
  fakeTrees.reset();

  TreeQueue DataTrees("DataFiles.txt");

  unsigned int jsonIndex(0);
  //loop over files
  //while (false){
  while( DataTrees.nextFile() ){
    eventTree = DataTrees.getTree("Events");
    DataTrees.status();

    //really need to write this better!
    JsonParser jsonParser;
    jsonParser.Initialize(jsonList.at(jsonIndex)); 
    cout << "Running over json file: " << jsonList.at(jsonIndex) << "\n";

    if (jsonIndex < jsonListSize ) ++jsonIndex;


    // Set branch address to structures that will store the info  
    eventTree->SetBranchAddress("Info",       &info);          TBranch *infoBr       = eventTree->GetBranch("Info");
    eventTree->SetBranchAddress("Dielectron" ,  &diElectronArr  ); TBranch *diElectronBr   = eventTree->GetBranch("Dielectron");    
    eventTree->SetBranchAddress("Electron" ,  &electronArr  ); TBranch *electronBr   = eventTree->GetBranch("Electron");

    UInt_t maxEvents = eventTree->GetEntries();

    for(UInt_t ientry=0; ientry<maxEvents; ++ientry) {           
      if(ientry >= maxEvents) break;	
    
      infoBr->GetEntry(ientry);		
      diElectronArr->Clear(); 
      diElectronBr->GetEntry(ientry); 
      electronArr->Clear(); 
      electronBr->GetEntry(ientry);

      if(!jsonParser.HasRunLumi(info->runNum, info->lumiSec)) continue;  // not certified run? Skip to next event...

      UInt_t  triggerEle = kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL
            | kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL;

      if(!(info->triggerBits & triggerEle)) continue;	

      vector<int> photonIndex;

      UInt_t maxDiElectrons = diElectronArr->GetEntriesFast();

      map<Float_t,UInt_t> elecEMap;
      map<Float_t,UInt_t> elecCandEMap;
      map<Float_t,UInt_t> elecEFailMap;
      pair<map<Float_t, UInt_t>::iterator, bool> insertPair;
      vector<UInt_t> elecCand;//record electrons for second fake rate method (systematic test)
      vector<UInt_t> photonCand;//record photon cands for seconf fake rate  method (systematic test)

      //We are now looking for an electron pair that passes the cluster cuts but electrons fail 
      //the final electron cuts

      
      for(UInt_t j=0; j<maxDiElectrons; ++j) { //loop over electron pairs
        mithep::TDielectron *dielectron = (mithep::TDielectron*)((*diElectronArr)[j]);

        //et cuts
        Double_t scEt1 = dielectron->scEt_1;
	Double_t scEt2 = dielectron->scEt_2;
	//skip this pairing
        if( ! ( (scEt1>20 && scEt2>10) || (scEt1>10 && scEt2>20) ) ) continue;

	//apply cuts to electron not associated with photon cluster
 
        bool passEle1(true);
        bool passEle2(true);
	bool photonCand1(false);
	bool photonCand2(false);
	bool gsfFtrkiso1(false);
	bool gsfFtrkiso2(false);

        //cout << dielectron->scID_1 << "\t" << dielectron->scID_2 << "\n"; 
        if (dielectron->scEt_1 < 10)  passEle1 = false;	 //reject electron if less than 20GeV
	if( !(dielectron->isEcalDriven_1)   )   passEle1 = false;  // not ECAL seeded electrons? Skip to next event...
        if (dielectron->nExpHitsInner_1 > 0)  passEle1 = false;    
        if((fabs(dielectron->scEta_1)>kGAP_LOW) && (fabs(dielectron->scEta_1)<kGAP_HIGH)) passEle1 = false;
	if( fabs(dielectron->scEta_1) > REGION_3 )   passEle1 = false;  // outside eta range? Skip to next event...
        if ((dielectron->emIso03_1/dielectron->pt_1) > 0.2) passEle1 = false; //iso cut
	if ((dielectron->hadIso03_1/dielectron->pt_1) > 0.2) passEle1 = false; //iso cut
        if ((dielectron->trkIso03_1/dielectron->pt_1) > 0.2) passEle1 = false; //iso cut
        //cuts in regions of eta

	if (fabs(dielectron->scEta_1) < kGAP_LOW) {
          if( dielectron->sigiEtaiEta_1	> 0.011) passEle1 = false;
          if (dielectron->HoverE_1 > 0.10) passEle1 = false;//H/E cut for photons 
	  if (fabs(dielectron->deltaPhiIn_1) > 0.15) passEle1 = false;	 
	  if (fabs(dielectron->deltaEtaIn_1) > 0.01) passEle1 = false;
	} else {
	  if( dielectron->sigiEtaiEta_1	> 0.031) passEle1 = false;
          if (dielectron->HoverE_1 > 0.075) passEle1 = false;//H/E cut for photons
	  if (fabs(dielectron->deltaPhiIn_1) > 0.10) passEle1 = false;	 
	  if (fabs(dielectron->deltaEtaIn_1) > 0.01) passEle1 = false; 
	}
        
        if (dielectron->scEt_2 < 10)  passEle2 = false;	 //reject electron if less than 20GeV
	if( !(dielectron->isEcalDriven_2)   )   passEle2 = false;  // not ECAL seeded electrons? Skip to next event...
        if (dielectron->nExpHitsInner_2 > 0)  passEle2 = false;    
        if((fabs(dielectron->scEta_2)>kGAP_LOW) && (fabs(dielectron->scEta_2)<kGAP_HIGH)) passEle2 = false;
	if( fabs(dielectron->scEta_2) > REGION_3 )   passEle2 = false;  // outside eta range? Skip to next event...
	if ((dielectron->emIso03_2/dielectron->pt_2) > 0.3) passEle2 = false; //iso cut
	if ((dielectron->hadIso03_2/dielectron->pt_2) > 0.3) passEle2 = false; //iso cut

        //cuts in regions of eta

	if (fabs(dielectron->scEta_2) < kGAP_LOW) {
          if( dielectron->sigiEtaiEta_2	> 0.01) passEle2 = false;
          if (dielectron->HoverE_2 > 0.04) passEle2 = false;//H/E cut for photons 
	  if (fabs(dielectron->deltaPhiIn_2) > 0.15) passEle2 = false;	 
	  if (fabs(dielectron->deltaEtaIn_2) > 0.01) passEle2 = false;
	} else {
	  if( dielectron->sigiEtaiEta_2	> 0.03) passEle2 = false;
          if (dielectron->HoverE_2 > 0.10) passEle2 = false;//H/E cut for photons 
	  if (fabs(dielectron->deltaPhiIn_2) > 0.10) passEle2 = false;	 
	  if (fabs(dielectron->deltaEtaIn_2) > 0.01) passEle2 = false;
	}

	if (passEle1){
          insertPair = elecCandEMap.insert(make_pair(dielectron->scEt_1, j));
	  photonCand1 = true;
	}
	if (passEle2){
          insertPair = elecCandEMap.insert(make_pair(dielectron->scEt_2, j));
	  photonCand2 = true;
	}

	//pass gsf electron cut but fail trk isolation
        if (passEle1 && passEle2){ 	 
	  if (dielectron->trkIso03_1 > 5.0) gsfFtrkiso1 = true; //iso cut
	  if (dielectron->trkIso03_2 > 5.0) gsfFtrkiso2 = true; //iso cut	  
	}

        if ( ! (passSmurf_elec1(dielectron))) passEle1 = false; //pass electron id smurf cuts
        if ( ! (passSmurf_elec2(dielectron))) passEle2 = false; //pass electron id smurf cuts
	
	if (photonCand1){
	  if (passEle1){
	    insertPair = elecEMap.insert(make_pair(dielectron->scEt_1, j));
	  } else {
	    elecEFailMap.insert(make_pair(dielectron->scEt_1, j));
	  }
	}
	if (photonCand2){
	  if (passEle2){
	    insertPair = elecEMap.insert(make_pair(dielectron->scEt_2, j));
	  } else {
	    elecEFailMap.insert(make_pair(dielectron->scEt_2, j));
	  }
	}

	//select electron pairs that passed "cluster" cuts for both electrons but then
	//failed full electron cuts for both electrons
	

	//select electron that pass cuts and gsf electron that failed trk iso cut
	if (photonCand1 && photonCand2){//both passed gsfelectron cuts
	  if (passEle1 && gsfFtrkiso2 ){
            TLorentzVector photon1_4v, photon2_4v, ee4v; 
	      photon1_4v.SetPtEtaPhiM(dielectron->scEt_1,dielectron->eta_1,dielectron->phi_1,0.000511);
	      photon2_4v.SetPtEtaPhiM(dielectron->scEt_2, dielectron->eta_2, dielectron->phi_2, 0.000511);
	      ee4v = photon1_4v + photon2_4v;
	      double mass = ee4v.M();
	      double rapidity = ee4v.Rapidity();
         
	      //double fRateVal = fRate_map.upper_bound(dielec->eta_1)->second.upper_bound(dielec->scEt_1)->second;
	      double fRateVal2 = fRate_map.upper_bound(dielectron->eta_2)->second.upper_bound(dielectron->scEt_2)->second;
	      //weight by the fake rate/(1-fake rate) * fake rate/(1-fake rate)
	      double fakeFactor = fRateVal2/(1-fRateVal2);

	      testEGsfF->Fill(mass,fakeFactor);
	      testEGsfFMass->Fill(mass);

	      if ((mass > 20) && (mass <= 30)){
		rapMass20to30->Fill(rapidity,fakeFactor);
	      }
	      if ((mass > 30) && (mass <= 45)){
		rapMass30to45->Fill(rapidity,fakeFactor);
	      }
	      if ((mass > 45) && (mass <= 60)){
		rapMass45to60->Fill(rapidity,fakeFactor);
	      }
	      if ((mass > 60) && (mass <= 120)){
		rapMass60to120->Fill(rapidity,fakeFactor);
	      }
	      if ((mass > 120) && (mass <= 200)){
		rapMass120to200->Fill(rapidity,fakeFactor);
	      }
	      if ((mass > 200) && (mass <= 1500)){
		rapMass200to1500->Fill(rapidity,fakeFactor);
	      }

	  }
	  if (passEle2 && gsfFtrkiso1 ){
	    TLorentzVector photon1_4v, photon2_4v, ee4v;
	      photon1_4v.SetPtEtaPhiM(dielectron->scEt_1,dielectron->eta_1,dielectron->phi_1,0.000511);
	      photon2_4v.SetPtEtaPhiM(dielectron->scEt_2, dielectron->eta_2, dielectron->phi_2, 0.000511);
	      ee4v = photon1_4v + photon2_4v;
	      double mass = ee4v.M();
	      double rapidity = ee4v.Rapidity();

	      //double fRateVal = fRate_map.upper_bound(dielec->eta_1)->second.upper_bound(dielec->scEt_1)->second;
	      double fRateVal = fRate_map.upper_bound(dielectron->eta_1)->second.upper_bound(dielectron->scEt_1)->second;
	      //weight by the fake rate/(1-fake rate) * fake rate/(1-fake rate)

	      double fakeFactor = fRateVal/(1-fRateVal);

	      testEGsfF->Fill(mass,fRateVal/(1-fRateVal));
	      testEGsfFMass->Fill(mass);

              if ((mass > 20) && (mass <= 30)){
		rapMass20to30->Fill(rapidity,fakeFactor);
	      }
	      if ((mass > 30) && (mass <= 45)){
		rapMass30to45->Fill(rapidity,fakeFactor);
	      }
	      if ((mass > 45) && (mass <= 60)){
		rapMass45to60->Fill(rapidity,fakeFactor);
	      }
	      if ((mass > 60) && (mass <= 120)){
		rapMass60to120->Fill(rapidity,fakeFactor);
	      }
	      if ((mass > 120) && (mass <= 200)){
		rapMass120to200->Fill(rapidity,fakeFactor);
	      }
	      if ((mass > 200) && (mass <= 1500)){
		rapMass200to1500->Fill(rapidity,fakeFactor);
	      }
 

	    //elecFailGsfIndex1.push_back(j);
	  }
	}

        if (passEle1 || passEle2) { 
	  elecCand.push_back(j);
        //correctEnergyScale(elec, info); need to apply energy scale corrections
	}


      }//end of loop over electron pairs

    }//end of loop over events
  }//end of loop over files
  
  //Set current output file
  outFile->cd();

  //Write histos to file
  cluster_et->Write();
  electron_et->Write();
  fRate->Write();
  fRate_reg->Write();
  fRate_reg2->Write();
  fRate_reg3->Write();
  pfMet->Write();
  pfMetx->Write();
  pfMety->Write();
  caloMet->Write();
  etDiff->Write();
  etaDiff->Write();
  cand_electron_et->Write();
  testEGsfFMass->Write();
  testEGsfF->Write();
  rapMass20to30->Write();
  rapMass30to45->Write();
  rapMass45to60->Write();
  rapMass60to120->Write();
  rapMass120to200->Write();
  rapMass200to1500->Write();     
 
  //Write out file
  outFile->Write();
  //delete outTree;
  outFile->Close();        
  delete outFile;

  delete infile;
  infile=0, eventTree=0;

  delete info;
  evtfile.close();

  cout << "\n";
  cout << " <> Output saved in " << outputDir << "/" << "\n";
  cout << "\n";
        
  gBenchmark->Show("fakeRate");       
} 


//=== FUNCTION DEFINITIONS ======================================================================================

//--------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------
/*void fillData(ZemuData *data, 
	      const mithep::TEventInfo *info, 
	      const double mass, const double pt,
	      const mithep::TMuon *muon, 
	      const mithep::TElectron *electron, 
              const UInt_t npv, const Double_t weight)
{
  data->runNum         = info->runNum;
  data->evtNum         = info->evtNum;
  data->lumiSec        = info->lumiSec;
  data->nPV            = npv;
  data->mass           = mass;
  data->pt             = pt;

  // Muon
  data->pt_mu           = muon->pt;
  data->eta_mu          = muon->eta;
  data->hltMatchBits_mu = muon->hltMatchBits;
  data->q_mu            = muon->q;

  // Electron
  data->pt_e           = electron->pt;
  data->eta_e          = electron->eta;
  data->scEt_e         = electron->scEt;
  data->scEta_e        = electron->scEta;
  data->hltMatchBits_e = electron->hltMatchBits;
  data->q_e            = electron->q;

  data->weight         = weight;
}
*/
//--------------------------------------------------------------------------------------------------
 void correctEnergyScale(mithep::TElectron *electron, mithep::TEventInfo *info){

   //
   // Energy scale corrections for data
   //

   // The factors below are for WP80 electrons and (almost certain) 
   // have been derived for Sep17 ReReco of 38X.
   
   Double_t scaleEB = 1.0045;
   Double_t scaleEE = 1.0171;
   if(info->runNum > 144114) {
     scaleEB = 1.0071;
     scaleEE = 1.0528;
   }

   Bool_t   isB  = (fabs(electron->scEta) < kGAP_LOW );
   TLorentzVector ele; 
   ele.SetPtEtaPhiM(electron->pt,electron->eta,electron->phi,0.000511);

   electron->scEt = (isB) ? scaleEB*(electron->scEt) : scaleEE*(electron->scEt);

   if(isB) { ele *= scaleEB; } else { ele *= scaleEE; }
   electron->pt  = ele.Pt();
   electron->eta = ele.Eta();
   electron->phi = ele.Phi();
            
   return;
 }

//===================================================================================================
//Function to copy TH1F histos

 void TH1FCopy(const TH1F* orig, TH1F* copy){
    int numXbins = orig->GetNbinsX();
    //check if histos are the same size. They need to be! 
    if (numXbins != copy->GetNbinsX()){
      throw "Histograms have different number of bins";
    }

    for (int i = 1; i <= numXbins; ++i ){
      copy->SetBinContent(i,orig->GetBinContent(i));
      copy->SetBinError(i,orig->GetBinError(i));
    }
  }
//===============================================================================
//Function to apply smurf cuts to a single electron leg of a dielectron object

Bool_t passSmurf_elec1(const mithep::TDielectron *dielectron)
{
  if(fabs(dielectron->d0_1) > 0.02) return kFALSE;
  if(fabs(dielectron->dz_1) > 0.1)  return kFALSE;  
  
  // conversion rejection
  if(dielectron->nExpHitsInner_1 > 0) return kFALSE;
  if(dielectron->isConv_1)            return kFALSE;  
  
  // barrel/endcap dependent requirments      
  if(fabs(dielectron->scEta_1)<1.479) {
    // barrel
    if(dielectron->pfIso04_1 > 0.13*(dielectron->pt_1)) return kFALSE;
    
    if(dielectron->pt_1>20) {
      if(dielectron->sigiEtaiEta_1	> 0.01)  return kFALSE;
      if(fabs(dielectron->deltaPhiIn_1) > 0.06)  return kFALSE;
      if(fabs(dielectron->deltaEtaIn_1) > 0.004) return kFALSE;
      if(dielectron->HoverE_1	        > 0.04)  return kFALSE;
      
    } else {
      if(dielectron->sigiEtaiEta_1	> 0.01)  return kFALSE;
      if(fabs(dielectron->deltaPhiIn_1) > 0.03)  return kFALSE;
      if(fabs(dielectron->deltaEtaIn_1) > 0.004) return kFALSE;
      if(dielectron->HoverE_1	        > 0.025) return kFALSE;    
    }
    
  } else {
    // endcap
    if(dielectron->pfIso04_1 > 0.09*(dielectron->pt_1)) return kFALSE;
    
    if(dielectron->pt_1>20) {
      if(dielectron->sigiEtaiEta_1	> 0.03)  return kFALSE;
      if(fabs(dielectron->deltaPhiIn_1) > 0.03)  return kFALSE;
      if(fabs(dielectron->deltaEtaIn_1) > 0.007) return kFALSE;
      if(dielectron->HoverE_1	        > 0.10)  return kFALSE;
      
    } else {
      if(dielectron->sigiEtaiEta_1	> 0.03)  return kFALSE;
      if(fabs(dielectron->deltaPhiIn_1) > 0.02)  return kFALSE;
      if(fabs(dielectron->deltaEtaIn_1) > 0.005) return kFALSE;
      if(dielectron->HoverE_1	        > 0.10)  return kFALSE;      
    }
  }
  
  if(dielectron->pt_1 < 20 &&
     !((dielectron->fBrem_1>0.15) || (fabs(dielectron->scEta_1)<1 && dielectron->EoverP_1>0.95)))
    return kFALSE;
  
  return kTRUE;
}

//===============================================================================
//Function to apply smurf cuts to a single electron leg of a dielectron object

Bool_t passSmurf_elec2(const mithep::TDielectron *dielectron)
{
  if(fabs(dielectron->d0_2) > 0.02) return kFALSE;
  if(fabs(dielectron->dz_2) > 0.1)  return kFALSE;
  
  // conversion rejection
  if(dielectron->nExpHitsInner_2 > 0) return kFALSE;
  if(dielectron->isConv_2)            return kFALSE;
    
  if(fabs(dielectron->scEta_2)<1.479) {
    // barrel
    if(dielectron->pfIso04_2 > 0.13*(dielectron->pt_2)) return kFALSE;
    
    if(dielectron->pt_2>20) {
      if(dielectron->sigiEtaiEta_2	> 0.01)  return kFALSE;
      if(fabs(dielectron->deltaPhiIn_2) > 0.06)  return kFALSE;
      if(fabs(dielectron->deltaEtaIn_2) > 0.004) return kFALSE;
      if(dielectron->HoverE_2	        > 0.04)  return kFALSE;
      
    } else {
      if(dielectron->sigiEtaiEta_2	> 0.01)  return kFALSE;
      if(fabs(dielectron->deltaPhiIn_2) > 0.03)  return kFALSE;
      if(fabs(dielectron->deltaEtaIn_2) > 0.004) return kFALSE;
      if(dielectron->HoverE_2	        > 0.025) return kFALSE;    
    }
    
  } else {
    // endcap
    if(dielectron->pfIso04_2 > 0.09*(dielectron->pt_2)) return kFALSE;
    
    if(dielectron->pt_2>20) {
      if(dielectron->sigiEtaiEta_2	> 0.03)  return kFALSE;
      if(fabs(dielectron->deltaPhiIn_2) > 0.03)  return kFALSE;
      if(fabs(dielectron->deltaEtaIn_2) > 0.007) return kFALSE;
      if(dielectron->HoverE_2	        > 0.10)  return kFALSE;
      
    } else {
      if(dielectron->sigiEtaiEta_2	> 0.03)  return kFALSE;
      if(fabs(dielectron->deltaPhiIn_2) > 0.02)  return kFALSE;
      if(fabs(dielectron->deltaEtaIn_2) > 0.005) return kFALSE;
      if(dielectron->HoverE_2	        > 0.10)  return kFALSE;      
    }
  }
  
  if(dielectron->pt_2 < 20 && 
     !((dielectron->fBrem_2>0.15) || (fabs(dielectron->scEta_2)<1 && dielectron->EoverP_2>0.95)))
    return kFALSE;
  
  return kTRUE;
}
