#ifndef tnpSelectEvents_HH
#define tnpSelectEvents_HH

#define tnpSelectEventsIsObject
#define esfSelectEventsIsObject

#if defined(tnpSelectEventsIsObject) || defined(esfSelectEventsIsObject)
#include <TObject.h>
#endif


#include <TTree.h>
#include <TBranch.h>

#include "../Include/PUReweight.hh"

// ------------------------------------------------------

#ifdef tnpSelectEventsIsObject
class tnpSelectEvent_t : public TObject 
#else
struct tnpSelectEvent_t
#endif
{
public:
  typedef enum { _skipWeight, _dontSkipWeight } TCreateBranchesOption_t;
public:

  Double_t mass, yZ; 
  Double_t et,eta; // probe
  UInt_t nGoodPV;
  Double_t eventWeight,weight;

  void assign(Double_t mass_1, Double_t yZee1, Double_t et_1, Double_t eta_1, 
	      UInt_t nGoodPV_1, Double_t event_weight1, Double_t weight1) {
    mass=mass_1; yZ=yZee1; et=et_1; eta=eta_1; 
    nGoodPV=nGoodPV_1;
    eventWeight=event_weight1; weight=weight1;
  }

  void createBranches(TTree *tree, TCreateBranchesOption_t opt) {
    tree->Branch("mass",&this->mass,"mass/D");
    tree->Branch("yZee",&this->yZ,"yZee/D");
    tree->Branch("et",&this->et,"et/D");
    tree->Branch("eta",&this->eta,"eta/D");
    tree->Branch("nGoodPV",&this->nGoodPV,"nGoodPV/i");
    tree->Branch("eventWeight",&this->eventWeight,"eventWeight/D");
    if (opt!=_skipWeight) tree->Branch("weight",&this->weight,"weight/D");
  }

  void setBranchAddress(TTree *tree) {
    tree->SetBranchAddress("mass",&this->mass);
    tree->SetBranchAddress("yZee",&this->yZ);
    tree->SetBranchAddress("et",&this->et);
    tree->SetBranchAddress("eta",&this->eta);
    tree->SetBranchAddress("nGoodPV",&this->nGoodPV);
    tree->SetBranchAddress("eventWeight",&this->eventWeight);
    tree->SetBranchAddress("weight",&this->weight);
  }

  bool insideMassWindow(double mass_low, double mass_high) const {
    return ((mass>=mass_low) && (mass<mass_high));
  }

  bool insideRapidityWindow(double rapidity_min, double rapidity_max) const {
    return ((yZ>=rapidity_min) && (yZ<rapidity_max));
  }

#ifdef tnpSelectEventsIsObject
  ClassDef(tnpSelectEvent_t,1)
#endif
};

// ------------------------------------------------------

#ifdef esfSelectEventsIsObject
class esfSelectEvent_t : public TObject 
#else
struct esfSelectEvent_t
#endif
{
public:

  Double_t genMass,genY,mass,y;
  Double_t et_1,eta_1,et_2,eta_2;
  Double_t weight;
  UInt_t nGoodPV;

  void assign(double _genMass, double _genY,
	      double _mass, double _y1, 
	      double _et1, double _eta1,
	      double _et2, double _eta2, double _weight, UInt_t _nGoodPV) {
    genMass=_genMass; genY=_genY;
    mass=_mass; y=_y1;
    et_1=_et1; eta_1=_eta1;
    et_2=_et2; eta_2=_eta2;
    weight=_weight;
    nGoodPV=_nGoodPV;
  }

  void createBranches(TTree *tree) {
    tree->Branch("genMass",&this->genMass,"genMass/D");
    tree->Branch("genY",&this->genY,"genY/D");
    tree->Branch("mass",&this->mass,"mass/D");
    tree->Branch("y",&this->y,"y/D");
    tree->Branch("et_1",&this->et_1,"et_1/D");
    tree->Branch("eta_1",&this->eta_1,"eta_1/D");
    tree->Branch("et_2",&this->et_2,"et_2/D");
    tree->Branch("eta_2",&this->eta_2,"eta_2/D");
    tree->Branch("weight",&this->weight,"weight/D");
    tree->Branch("nGoodPV",&this->nGoodPV,"nGoodPV/i");
  }

  void setBranchAddress(TTree *tree) {
    tree->SetBranchAddress("genMass",&this->genMass);
    tree->SetBranchAddress("genY",&this->genY);
    tree->SetBranchAddress("mass",&this->mass);
    tree->SetBranchAddress("y",&this->y);
    tree->SetBranchAddress("et_1",&this->et_1);
    tree->SetBranchAddress("eta_1",&this->eta_1);
    tree->SetBranchAddress("et_2",&this->et_2);
    tree->SetBranchAddress("eta_2",&this->eta_2);
    tree->SetBranchAddress("weight",&this->weight);
    tree->SetBranchAddress("nGoodPV",&this->nGoodPV);
  }

  bool insideMassWindow(double mass_low, double mass_high) const {
    return ((mass>=mass_low) && (mass<=mass_high));
  }

#ifdef esfSelectEventsIsObject
  ClassDef(esfSelectEvent_t,1)
#endif
};

// ------------------------------------------------------
// ------------------------------------------------------

// Load default PU distribution (data) from file <puReferenceFName> and use it to create
// weight branch in tag-and-probe selected events file <fname>.
// The selected events file is read to determine PU distribution (which is saved to file
// <savePUFName> in two histograms {savePUHistoNameBase}_pass and {savePUHistoNameBase}_fail

int CreatePUWeightedBranch(const TString &fName, 
			   const TString &puReferenceFName, const TString &puRefDistrName,
			   const TString &savePUFName, const TString &savePUHistoNameBase) {
  std::cout << "entered CreatePUWeightedBranch (" << fName << ")" << std::endl;
  // Open PU reference distribution file
  PUReweight_t puRef;
  int res=puRef.setFile(puReferenceFName) &&
    puRef.setReference(puRefDistrName);
  assert(res && puRef.getHRef());

  PUReweight_t puReweight;
  res=puReweight.setFile(savePUFName,1); // create a new file
  assert(res);

  TFile *selectedEventsFile= new TFile(fName,"UPDATE");
  assert(selectedEventsFile && selectedEventsFile->IsOpen());
  //std::cout << "selectedEventsFile tree entry counts: " << ((TTree*)selectedEventsFile->Get("passTree"))->GetEntries() << ", " << ((TTree*)selectedEventsFile->Get("failTree"))->GetEntries() << "\n";


  // loop over the trees (pass/fail)
  for (int fail=0; res && (fail<2); ++fail) {
    TString pass_fail_str= (fail) ? "_fail" : "_pass";
    TString treeName= (fail) ? "failTree" : "passTree";
    // get the tree
    TTree *tree = (TTree*)selectedEventsFile->Get(treeName); assert(tree);
    // set addresses and create a new branch
    double evWeight,pvWeight;
    UInt_t nGoodPV;
    tree->SetBranchAddress("eventWeight",&evWeight);   
    tree->SetBranchAddress("nGoodPV",&nGoodPV);
    tree->Branch("weight",&pvWeight,"weight/D");
    TBranch *evWeightBr= tree->GetBranch("eventWeight");
    TBranch *pvCountBr= tree->GetBranch("nGoodPV");
    TBranch *weightBr= tree->GetBranch("weight");
    assert(evWeightBr && pvCountBr && weightBr);
    // determine the PU distribution
    TString puDistrName=savePUHistoNameBase + pass_fail_str;
    res= puReweight.setActiveSample(puDistrName);  // create a histo container
    assert(res);
    for (UInt_t i=0; i<tree->GetEntries(); ++i) {
      pvCountBr->GetEntry(i);
      evWeightBr->GetEntry(i);
      puReweight.Fill(nGoodPV,evWeight);
    }
    // set weights
    res=puReweight.setReference(puRef.getHRef()) &&  // set reference distribution
      puReweight.prepareWeights(1);  // prepare weights and save histo
    assert(res);
    // fill the weight branch
    for (UInt_t i=0; i<tree->GetEntries(); ++i) {
      pvCountBr->GetEntry(i);
      pvWeight=evWeight * puReweight.getWeight(nGoodPV);
      weightBr->Fill();
    }
    selectedEventsFile->cd();
    tree->Write();
  }
  selectedEventsFile->Write();
  //std::cout << "selectedEventsFile tree entry counts: " << ((TTree*)selectedEventsFile->Get("passTree"))->GetEntries() << ", " << ((TTree*)selectedEventsFile->Get("failTree"))->GetEntries() << "\n";
  delete selectedEventsFile;
  return 1;
}

// ------------------------------------------------------

#endif
