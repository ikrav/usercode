#ifndef storeData_HH
#define storeData_HH

#include <TObject.h>
#include "../Include/ElectronEnergyScale.hh"
#include "../Include/ZeeData.hh"

#define selectedEventDataIsTObject
// -------------------------------------------------------------------------
// -------------------------------------------------------------------------

class SelectedEventData_t 
#ifdef selectedEventDataIsTObject
  : public TObject 
#endif
{
public:
  UInt_t runNum, evtNum, lumiSec;
  double weight;               // event weight
  double mass;                 // dielectron mass
  double scEta_1,scEta_2;      // scEta values of the dielectron candidate
  int nPV,nGoodPV;             // number of vertices, number of good vertices
public:
  SelectedEventData_t() : 
#ifdef selectedEventDataIsTObject
    TObject(), 
#endif
    runNum(0), evtNum(0), lumiSec(0), weight(1.), mass(0.), scEta_1(0.), scEta_2(0.), nPV(1), nGoodPV(1) 
  {}

  SelectedEventData_t(const SelectedEventData_t &a) : 
#ifdef selectedEventDataIsTObject
    TObject(a), 
#endif
    runNum(a.runNum), evtNum(a.evtNum), lumiSec(a.lumiSec), weight(a.weight), mass(a.mass), scEta_1(a.scEta_1), scEta_2(a.scEta_2), nPV(a.nPV), nGoodPV(a.nGoodPV) 
  {}

  SelectedEventData_t(const ZeeData &z) : 
#ifdef selectedEventDataIsTObject
    TObject(), 
#endif
    runNum(0), evtNum(0), lumiSec(0), weight(1.), mass(0.), scEta_1(0.), scEta_2(0.), nPV(1), nGoodPV(1) 
  { this->assign(z); }
  
  void assign(const ZeeData &z) {
    runNum=z.runNum; evtNum=z.evtNum; lumiSec=z.lumiSec;
    weight=z.weight;
    mass=z.mass;
    scEta_1=z.scEta_1; scEta_2=z.scEta_2;
    nPV=z.nPV; nGoodPV=z.nGoodPV;
  }

  int massInsideRange(double mass_min, double mass_max) const { return ((mass>=mass_min) && (mass<=mass_max)) ? 1:0; }

  int massInsideRange(double mass_min, double mass_max, const ElectronEnergyScale &escale, int applySmear) const {
    smearMass(escale,applySmear);
    return massInsideRange(mass_min,mass_max);
  }

  double smearMass(const ElectronEnergyScale &escale, int applySmear) const { 
    double mass_loc=mass;
    if (applySmear) mass_loc += escale.generateMCSmear(scEta_1,scEta_2); 
    return mass_loc;
  }

#ifdef selectedEventDataIsTObject
  ClassDef(SelectedEventData_t,1)
#endif
};

// -------------------------------------------------------------------------
// -------------------------------------------------------------------------

class LumiInfo_t {
public:
  UInt_t runNumMin,runNumMax;
  double lumiWeight;
public:
  LumiInfo_t() : runNumMin(0), runNumMax(0), lumiWeight(0.) {}
  LumiInfo_t(UInt_t run_min, UInt_t run_max, double lumi_weight) : runNumMin(run_min), runNumMax(run_max), lumiWeight(lumi_weight) {}
  LumiInfo_t(const LumiInfo_t &a) : runNumMin(a.runNumMin), runNumMax(a.runNumMax), lumiWeight(a.lumiWeight) {}

  int insideRange(UInt_t run) { return ((run>=runNumMin) && (run<=runNumMax)) ? 1:0; }

  friend std::ostream& operator<<(std::ostream& out, const LumiInfo_t &a) {
    out << " " << a.runNumMin << " -- " << a.runNumMax << "  " << a.lumiWeight;
    return out;
  }

};

// -------------------------------------------------------------------------
// -------------------------------------------------------------------------
/*
class PVDistr_t : public TH1F {
  void Reset() {
    for (int i=1; i<this->GetNbinsX(); ++i) {
      this->SetBinContent(i,0); this->SetBinError(i,0.);
    }
  }
  ClassDef(PVDistr_t,1)
};
*/

// -------------------------------------------------------------------------
// -------------------------------------------------------------------------




#endif
