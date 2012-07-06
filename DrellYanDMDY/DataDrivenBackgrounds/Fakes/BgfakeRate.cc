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
//  * outputs ROOT files of events passing selection for each sample, 
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
#include <set>
#include <utility>                  // use make_pair for inserting pairs into a map 

#include "../../Include/MitStyleRemix.hh"  // style settings for drawing
//#include "MyTools.hh"        // miscellaneous helper functions
//#include "CSample.hh"        // helper class for organizing input ntuple files

// define structures to read in ntuple
#include "../../Include/EWKAnaDefs.hh"
#include "../../Include/TEventInfo.hh"
#include "../../Include/TElectron.hh"
#include "../../Include/TVertex.hh"

#include "TreeQueue.hh"

// lumi section selection with JSON files
// #include "RunLumiRangeMap.h"
//#include "JsonParser.hh"

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

void correctEnergyScale(mithep::TElectron *electron, mithep::TEventInfo *info);

void TH1FCopy(const TH1F* orig, TH1F* copy);
Bool_t passSmurf_elec1(const mithep::TDielectron *dielectron);
Bool_t passSmurf_elec2(const mithep::TDielectron *dielectron);


//===Pull objects and functions from namespace===================================================================

using std::string;
using std::cout;
using std::fstream;
using std::map;
using std::make_pair;
using boost::scoped_ptr;
using boost::shared_ptr;

//===============================================
//conflist class
//=================================================
//need to put this in its own .cc and .hh file

class ConfList : public TreeQueue{
  public:
  ConfList(string listOfRootFiles);
  bool nextFile(); //method to jump to next file
  void reset();
  double getVal() const {return *current_xsec;}

  private:
  vector<double> xseclist;
  vector<double>::iterator xsec_iter;
  bool storeFilenames();
  vector<double>::iterator current_xsec;
};


ConfList::ConfList(string listOfRootFiles): TreeQueue(listOfRootFiles, true){
  if (storeFilenames()){//anyfilenames stored?
    reset();//Then point to the first file
  } else {
    //should I throw an exception
    throw "empty list of files";
  }
}

bool ConfList::storeFilenames(){  
  string line;
  xseclist.clear();//clear vector
  filelist.clear();//clear vector
  while (rootFileStream.good()){
    getline(rootFileStream,line);
    if (line[0] != '#'){//ignore comment lines
      if (line.size() != 0) {//ignore empty lines
	istringstream concat(line);
	string rootfile;
	double xsec;
	concat >> rootfile >> xsec;
	filelist.push_back(rootfile);
        xseclist.push_back(xsec);
      }
    }
  }
  //need to clear file for reuse
  rootFileStream.clear();

  if (filelist.size() != xseclist.size()){
    throw "size mismatch";
  }

  if (!filelist.empty()) {//same as checking if xseclist is empty
    return true;
  } else {
    return false;
  }
}

bool ConfList::nextFile(){
  bool outcome;
  if (TreeQueue::nextFile()){
    if (xsec_iter != xseclist.end()){
      current_xsec = xsec_iter;
      ++xsec_iter;
      outcome = true;
    } else {
      outcome = false;
    }
  } else {
    outcome = false;
  }

  return outcome;
}

void ConfList::reset(){
  TreeQueue::reset();
  xsec_iter = xseclist.begin();
  current_xsec =  xsec_iter;
}

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

//=== MAIN MACRO =================================================================================================

void BgfakeRate() 
{  
  gBenchmark->Start("BgfakeRate");

  
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================
  
  TString  outputDir;         // output directory
  TString  trigName;          // which trigger to use
  TString  format;            // plot format

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
  
  //
  // Access samples and fill histograms
  //  
  TFile *infile=0;
  TTree *eventTree=0;  
  
  // Data structures to store info from TTrees
  mithep::TEventInfo *info    = new mithep::TEventInfo();
  TClonesArray *diElectronArr = new TClonesArray("mithep::TDielectron");

  //Set up histos to store events
  //  const int nBins(30),xmin(0), xmax(300);
  

  double lumi(4679.7);

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
  TString outName = TString("bgplots.root");
  TFile *outFile = new TFile(outName,"RECREATE");

  //
  // loop through files
  //
  TString ts_evtfname(evtfname);

  TString inName = TString("plots.root");
  infile = new TFile(inName);
  
  const int rap_Bins = 24;
  const int rap_Bins_highMass = 12;
  float rap_Min = 0;
  float rap_Max = 2.4;

  TH1F* fRate_reg = (TH1F*)infile->Get("fRate_reg");
  TH1F* fRate_reg2 = (TH1F*)infile->Get("fRate_reg2");
  TH1F* fRate_reg3 = (TH1F*)infile->Get("fRate_reg3");
  TH1F* testEGsfFMass = (TH1F*)infile->Get("testEGsfFMass");
  TH1F* testEGsfF = (TH1F*)infile->Get("testEGsfF");

  scoped_ptr<TH1F> eEt(new TH1F("eEt","",30,0,300));
  scoped_ptr<TH1F> eeHist(new TH1F("eeHist","",nMassBins,massBinLimits));
  scoped_ptr<TH1F> eeHist2 = scoped_ptr<TH1F>(new TH1F("eeHist2","",nMassBins,massBinLimits));
  scoped_ptr<TH1F> eeHist3 = scoped_ptr<TH1F>(new TH1F("eeHist3","",nMassBins,massBinLimits));
    
  scoped_ptr<TH1F> BkgEGsfFMass = scoped_ptr<TH1F>(new TH1F("BkgEGsfFMass","",nMassBins,massBinLimits));
  scoped_ptr<TH1F> BkgEGsfF = scoped_ptr<TH1F>(new TH1F("BkgEGsfF","",nMassBins,massBinLimits));

  scoped_ptr<TH1F> BkgRapMass20to30 = scoped_ptr<TH1F>(new TH1F("BkgRapMass20to30","",rap_Bins, rap_Min, rap_Max));
  scoped_ptr<TH1F> BkgRapMass30to45 = scoped_ptr<TH1F>(new TH1F("BkgRapMass30to45","",rap_Bins, rap_Min, rap_Max));
  scoped_ptr<TH1F> BkgRapMass45to60 = scoped_ptr<TH1F>(new TH1F("BkgRapMass45to60","",rap_Bins, rap_Min, rap_Max));
  scoped_ptr<TH1F> BkgRapMass60to120 = scoped_ptr<TH1F>(new TH1F("BkgRapMass60to120","",rap_Bins, rap_Min, rap_Max));
  scoped_ptr<TH1F> BkgRapMass120to200 = scoped_ptr<TH1F>(new TH1F("BkgRapMass120to200","",rap_Bins, rap_Min, rap_Max));
  scoped_ptr<TH1F> BkgRapMass200to1500 = scoped_ptr<TH1F>(new TH1F("BkgRapMass200to1500","",rap_Bins_highMass, rap_Min, rap_Max));

  //Put fake rates into a map. Change to use boost unordered map for better performance

  map<double,double> rate;
  map<double,double> rate2;
  map<double,double> rate3;

  // Put bin and fake rate values in a map
  const int eTBins = fRate_reg->GetNbinsX();
  cout << "Number of bins is: " << eTBins << "\n"; 

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
 
  // look at MC events identified as fakes
  ConfList fakeTrees("BgfakeR.txt");  

   //start from first file again
  fakeTrees.reset();
  //loop over files
  while( fakeTrees.nextFile() ){
    eventTree = fakeTrees.getTree("Events");
    fakeTrees.status();

    double xsec = fakeTrees.getVal();    

    // Set branch address to structures that will store the info  
    eventTree->SetBranchAddress("Info",       &info);          TBranch *infoBr       = eventTree->GetBranch("Info");
    eventTree->SetBranchAddress("Dielectron" ,  &diElectronArr  ); TBranch *diElectronBr   = eventTree->GetBranch("Dielectron");    

    double mcScaleFactor(1.0);
    UInt_t maxEvents = eventTree->GetEntries();
    double weight = (xsec * lumi * mcScaleFactor)/maxEvents;

    for(UInt_t ientry=0; ientry<maxEvents; ++ientry) {           
      if(ientry >= maxEvents) break;	
    
      infoBr->GetEntry(ientry);	
      diElectronArr->Clear(); 
      diElectronBr->GetEntry(ientry);	

      UInt_t maxDiElectrons = diElectronArr->GetEntriesFast();

      map<Float_t,UInt_t> elecEMap;
      map<Float_t,UInt_t> elecCandEMap;
      map<Float_t,UInt_t> elecEFailMap;
      pair<map<Float_t, UInt_t>::iterator, bool> insertPair;
      vector<UInt_t> elecCand;//record electrons for second fake rate method (systematic test)
      vector<UInt_t> photonCand;//record photon cands for seconf fake rate  method (systematic test)

       for(UInt_t j=0; j<maxDiElectrons; ++j) { //loop over electron pairs
        mithep::TDielectron *dielectron = (mithep::TDielectron*)((*diElectronArr)[j]);
	
        //et cuts
        Double_t scEt1 = dielectron->scEt_1;
	Double_t scEt2 = dielectron->scEt_2;
	//skip this pairing
        if( ! ( (scEt1>20 && scEt2>10) || (scEt1>10 && scEt2>20) ) ) continue;

        //apply cuts to electron not associated with photon cluster
        eEt->Fill(dielectron->scEt_1,weight);
        eEt->Fill(dielectron->scEt_2,weight); 

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
	if (photonCand1 && photonCand2){
	  if ( (!(passEle1)) && (!(passEle2)) ){
	    //elecPairFailIndex.push_back(j);
	  }
	}

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

	    BkgEGsfF->Fill(mass,(fakeFactor)*weight );
	    BkgEGsfFMass->Fill(mass,weight);

	    if ((mass > 20) && (mass <= 30)){
	      BkgRapMass20to30->Fill(rapidity,fakeFactor*weight);
	    }
	    if ((mass > 30) && (mass <= 45)){
	      BkgRapMass30to45->Fill(rapidity,fakeFactor*weight);
	    }
	    if ((mass > 45) && (mass <= 60)){
	      BkgRapMass45to60->Fill(rapidity,fakeFactor*weight);
	    }
	    if ((mass > 60) && (mass <= 120)){
	      BkgRapMass60to120->Fill(rapidity,fakeFactor*weight);
	    }
	    if ((mass > 120) && (mass <= 200)){
	      BkgRapMass120to200->Fill(rapidity,fakeFactor*weight);
	    }
	    if ((mass > 200) && (mass <= 1500)){
	      BkgRapMass200to1500->Fill(rapidity,fakeFactor*weight);
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

	    BkgEGsfF->Fill(mass,(fakeFactor)*weight);
	    BkgEGsfFMass->Fill(mass,weight);

	    if ((mass > 20) && (mass <= 30)){
	      BkgRapMass20to30->Fill(rapidity,fakeFactor*weight);
	    }
	    if ((mass > 30) && (mass <= 45)){
	      BkgRapMass30to45->Fill(rapidity,fakeFactor*weight);
	    }
	    if ((mass > 45) && (mass <= 60)){
	      BkgRapMass45to60->Fill(rapidity,fakeFactor*weight);
	    }
	    if ((mass > 60) && (mass <= 120)){
	      BkgRapMass60to120->Fill(rapidity,fakeFactor*weight);
	    }
	    if ((mass > 120) && (mass <= 200)){
	      BkgRapMass120to200->Fill(rapidity,fakeFactor*weight);
	    }
	    if ((mass > 200) && (mass <= 1500)){
	      BkgRapMass200to1500->Fill(rapidity,fakeFactor*weight);
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

  // subtract MC ee histos from qcd histo
  testEGsfFMass->Add(BkgEGsfFMass.get(), -1.0);
  testEGsfF->Add(BkgEGsfF.get(), -1.0); 

  //Write histos to file

  testEGsfFMass->Write();
  testEGsfF->Write();
  BkgEGsfFMass->Write();
  BkgEGsfF->Write();

  BkgRapMass20to30->Write();
  BkgRapMass30to45->Write();
  BkgRapMass45to60->Write();
  BkgRapMass60to120->Write();
  BkgRapMass120to200->Write();
  BkgRapMass200to1500->Write();

  //Write out file
  outFile->Write();

  outFile->Close();        
  //delete outFile;

  //delete infile;
  infile=0, eventTree=0;

  //delete info;
  evtfile.close();

  cout << "\n";
  cout << " <> Output saved in " << outputDir << "/" << "\n";
  cout << "\n";
        
  gBenchmark->Show("BgfakeRate");       
} 


//=== FUNCTION DEFINITIONS ======================================================================================

//--------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------

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



