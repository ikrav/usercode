#ifndef tnpSelectEvents_HH
#define tnpSelectEvents_HH

#define tnpSelectEventsIsObject
#define esfSelectEventsIsObject

#if defined(tnpSelectEventsIsObject) || defined(esfSelectEventsIsObject)
#include <TObject.h>
#endif


#include <TTree.h>
#include <TBranch.h>

// ------------------------------------------------------

#ifdef tnpSelectEventsIsObject
class tnpSelectEvent_t : public TObject 
#else
struct tnpSelectEvent_t
#endif
{
public:

  Double_t mass, yZ; 
  Double_t et,eta; // probe
  UInt_t nGoodPV;
  Double_t weight;

  void assign(Double_t mass_1, Double_t yZee1, Double_t et_1, Double_t eta_1, UInt_t nGoodPV_1, Double_t event_weight) {
    mass=mass_1; yZ=yZee1; et=et_1; eta=eta_1; 
    nGoodPV=nGoodPV_1;
    weight=event_weight;
  }

  void createBranches(TTree *tree) {
    tree->Branch("mass",&this->mass,"mass/D");
    tree->Branch("yZee",&this->yZ,"yZee/D");
    tree->Branch("et",&this->et,"et/D");
    tree->Branch("eta",&this->eta,"eta/D");
    tree->Branch("nGoodPV",&this->nGoodPV,"nGoodPV/i");
    tree->Branch("weight",&this->weight,"weight/D");
  }

  void setBranchAddress(TTree *tree) {
    tree->SetBranchAddress("mass",&this->mass);
    tree->SetBranchAddress("yZee",&this->yZ);
    tree->SetBranchAddress("et",&this->et);
    tree->SetBranchAddress("eta",&this->eta);
    tree->SetBranchAddress("nGoodPV",&this->nGoodPV);
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

  Double_t genMass,mass,et_1,eta_1,et_2,eta_2;
  Double_t weight;
  UInt_t nGoodPV;

  void assign(double _genMass, double _mass, double _et1, double _eta1,
	      double _et2, double _eta2, double _weight, UInt_t _nGoodPV) {
    genMass=_genMass; mass=_mass;
    et_1=_et1; eta_1=_eta1;
    et_2=_et2; eta_2=_eta2;
    weight=_weight;
    nGoodPV=_nGoodPV;
  }

  void createBranches(TTree *tree) {
    tree->Branch("genMass",&this->genMass,"genMass/D");
    tree->Branch("mass",&this->mass,"mass/D");
    tree->Branch("et_1",&this->et_1,"et_1/D");
    tree->Branch("eta_1",&this->eta_1,"eta_1/D");
    tree->Branch("et_2",&this->et_2,"et_2/D");
    tree->Branch("eta_2",&this->eta_2,"eta_2/D");
    tree->Branch("weight",&this->weight,"weight/D");
    tree->Branch("nGoodPV",&this->nGoodPV,"nGoodPV/i");
  }

  void setBranchAddress(TTree *tree) {
    tree->SetBranchAddress("genMass",&this->genMass);
    tree->SetBranchAddress("mass",&this->mass);
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

#endif
