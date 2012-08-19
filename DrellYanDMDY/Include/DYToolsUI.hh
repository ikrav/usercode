#ifndef DYToolsUI_HH
#define DYToolsUI_HH

#include "../Include/DYTools.hh"

// ------------------------------------------------------------------
//   Names and recognition
// ------------------------------------------------------------------

TString MassBinningName(DYTools::TMassBinning_t set) {
  using namespace DYTools;
  TString name;
  switch(set) {
  case _MassBins_Undefined: name="UNDEFINED"; break;
  case _MassBins_2010: name="MassBins2010"; break;
  case _MassBins_2011: name="MassBins2011"; break;
  case _MassBins_2011_Lumi: name="MassBins2011Lumi"; break;
  case _MassBins_2011_2D: name="MassBins2011_2D"; break;
  default: name="UNKNOWN_MASS_BINNING";
  }
  return name;
}

// ------------------------------------------------------------------

TString CurrentMassBinningName() {
  return MassBinningName(DYTools::massBinningSet);
}

// ------------------------------------------------------------------
// ------------------------------------------------------------------

TString SystematicsStudyName(DYTools::TSystematicsStudy_t study) {
  using namespace DYTools;
  TString name;
  switch(study) {
  case NORMAL: name="NormalRun"; break;
  case RESOLUTION_STUDY: name="ResolutionStudy"; break;
  case FSR_STUDY: name="FSRStudy"; break;
  case ESCALE_RESIDUAL: name="EScale_residual"; break;
  case ESCALE_STUDY: name="EScale_study"; break;
  case ESCALE_STUDY_RND: name="EScale_study_randomized"; break;
  default: name="UNKNOWN_SYSTEMATICS_NAME";
  }
  return name;
}

// ------------------------------------------------------------------

TString TnPMethodName(DYTools::TTnPMethod_t method) { 
  using namespace DYTools;
  TString name;
  switch(method) {
  case COUNTnCOUNT: name="Count&Count"; break;
  case COUNTnFIT: name="Count&Fit"; break;
  case FITnFIT: name="Fit&Fit"; break;
  default: name="Unknown_TnP_method";
  }
  return name;
}

// ------------------------------------------------------------------

TString EfficiencyKindName(DYTools::TEfficiencyKind_t kind) {
  using namespace DYTools;
  TString name;
  switch(kind) {
  case RECO: name="RECO"; break;
  case ID: name="ID"; break;
  case HLT: name="HLT"; break;
  case HLT_leg1: name="HLTleg1"; break;
  case HLT_leg2: name="HLTleg2"; break;
  default: name="Unknown_Efficiency_Kind";
  }
  return name;
}

// ------------------------------------------------------------------

TString EtBinSetName(DYTools::TEtBinSet_t set) {
  using namespace DYTools;
  TString name;
  switch(set) {
  case ETBINS_UNDEFINED: name="EtBinsUndefined"; break;
  case ETBINS1: name="EtBins1"; break;
  case ETBINS5: name="EtBins5"; break;
  case ETBINS6: name="EtBins6"; break;
  case ETBINS7: name="EtBins7"; break;
  case ETBINS8: name="EtBins8"; break;
  case ETBINS9: name="EtBins9"; break;
  default: name="Unknown_Et_BinSet";
    std::cout << "EtBinSetName=<" << name << ">\n";
    assert(0);
  }
  return name;
}

// ------------------------------------------------------------------

TString EtaBinSetName(DYTools::TEtaBinSet_t set) {
  using namespace DYTools;
  TString name;
  switch(set) {
  case ETABINS_UNDEFINED: name="EtaBinsUndefined"; break;
  case ETABINS1: name="EtaBins1"; break;
  case ETABINS2: name="EtaBins2"; break;
  case ETABINS2Negs: name="EtaBins2Negs"; break;
  case ETABINS5: name="EtaBins5"; break;
  case ETABINS5Negs: name="EtaBins5Negs"; break;
  case ETABINS4test: name="EtaBins4test"; break;
  case ETABINS4testNegs: name="EtaBins4testNegs"; break;
  case ETABINS4alt: name="EtaBins4alt"; break;
  case ETABINS4altNegs: name="EtaBins4altNegs"; break;
  case ETABINS5alt: name="EtaBins5alt"; break;
  case ETABINS5altNegs: name="EtaBins5altNegs"; break;
  default: name="Unknown_Eta_BinSet";
    std::cout << "EtaBinSetName=<" << name << ">\n";
    assert(0);
  }
  return name;
}

// ------------------------------------------------------------------

TString DataKindName(DYTools::TDataKind_t kind) {
  using namespace DYTools;
  TString name;
  switch(kind) {
  case DATA: name="Data"; break;
  case MC: name="MC"; break;
  default: name="Unknown_DataKind";
  }
  return name;
}

// ------------------------------------------------------------------
// ------------------------------------------------------------------

DYTools::TSystematicsStudy_t DetermineSystematicsStudy(const TString &str) {
  using namespace DYTools;
  DYTools::TSystematicsStudy_t study=NORMAL;
  if (str.Contains("NORMAL") || str.Contains("NormalRun")) { study=NORMAL; }
  else if (str.Contains("RESOLUTION_STUDY") || str.Contains("ResolutionStudy")) study=RESOLUTION_STUDY;
  else if (str.Contains("FSR_STUDY") || str.Contains("FSRStudy")) study=FSR_STUDY;
  else {
    std::cout << "DetermineSystematicsStudy: failure at <" << str << ">\n";
    assert(0);
  }
  return study;
}

// ------------------------------------------------------------------

DYTools::TTnPMethod_t DetermineTnPMethod(const TString &str) {
  using namespace DYTools;
  DYTools::TTnPMethod_t method=COUNTnCOUNT;
  if (str.Contains("COUNTnCOUNT") || str.Contains("Count&Count")) method=COUNTnCOUNT;
  else if (str.Contains("COUNTnFIT") || str.Contains("Count&Fit")) method=COUNTnFIT;
  else if (str.Contains("FITnFIT") || str.Contains("Fit&Fit")) method=FITnFIT;
  else {
    std::cout << "DetermineTnPMethod failed at <" << str << ">\n";
    assert(0);
  }
  return method;
}

// ------------------------------------------------------------------

DYTools::TEfficiencyKind_t DetermineEfficiencyKind(const TString &str) {
  using namespace DYTools;
  DYTools::TEfficiencyKind_t kind=RECO;
  if (str.Contains("GSF") || str.Contains("Reco") || str.Contains("RECO")) kind=RECO;
  else if (str.Contains("ID")) kind=ID;
  else if (str.Contains("HLTleg1")) kind=HLT_leg1;
  else if (str.Contains("HLTleg2")) kind=HLT_leg2;
  else if (str.Contains("HLT")) kind=HLT;
  else {
    std::cout << "DetermineEfficiencyKind failed at <" << str << ">\n";
    assert(0);
  }
  return kind;
}

// ------------------------------------------------------------------

DYTools::TEtBinSet_t DetermineEtBinSet(const TString& str) {
  using namespace DYTools;
  DYTools::TEtBinSet_t kind=ETBINS1;
  if (str.Contains("ETBINS1") || str.Contains("EtBins1")) kind=ETBINS1;
  else if (str.Contains("ETBINS5") || str.Contains("EtBins5")) kind=ETBINS5;
  else if (str.Contains("ETBINS6") || str.Contains("EtBins6")) kind=ETBINS6;
  else if (str.Contains("ETBINS7") || str.Contains("EtBins7")) kind=ETBINS7;
  else if (str.Contains("ETBINS8") || str.Contains("EtBins8")) kind=ETBINS8;
  else if (str.Contains("ETBINS9") || str.Contains("EtBins9")) kind=ETBINS9;
  else {
    std::cout << "DetermineEtBinSet failed at <" << str  << ">\n";
    assert(0);
  }
  return kind;
}

// ------------------------------------------------------------------

DYTools::TEtaBinSet_t DetermineEtaBinSet(const TString& str) {
  using namespace DYTools;
  DYTools::TEtaBinSet_t kind=ETABINS1;
  if (str.Contains("ETABINS1") || str.Contains("EtaBins1")) kind=ETABINS1;
  else if (str.Contains("ETABINS2Negs") || str.Contains("EtaBins2Negs")) kind=ETABINS2Negs;
  else if (str.Contains("ETABINS2") || str.Contains("EtaBins2")) kind=ETABINS2;
  else if (str.Contains("ETABINS5Negs") || str.Contains("EtaBins5Negs")) kind=ETABINS5Negs;
  else if (str.Contains("ETABINS5altNegs") || str.Contains("EtaBins5altNegs")) kind=ETABINS5altNegs;
  else if (str.Contains("ETABINS5alt") || str.Contains("EtaBins5alt")) kind=ETABINS5alt;
  else if (str.Contains("ETABINS5") || str.Contains("EtaBins5")) kind=ETABINS5;
  else if (str.Contains("ETABINS4altNegs") || str.Contains("EtaBins4altNegs")) kind=ETABINS4altNegs;
  else if (str.Contains("ETABINS4alt") || str.Contains("EtaBins4alt")) kind=ETABINS4alt;
  else if (str.Contains("ETABINS4testNegs") || str.Contains("EtaBins4testNegs")) kind=ETABINS4testNegs;
  else if (str.Contains("ETABINS4test") || str.Contains("EtaBins4test")) kind=ETABINS4test;
  else {
    std::cout << "DetermineEtaBinSet failed at <" << str  << ">\n";
    assert(0);
  }
  return kind;
}

// ------------------------------------------------------------------

DYTools::TDataKind_t DetermineDataKind(const TString &str) {
  using namespace DYTools;
  DYTools::TDataKind_t kind=DATA;
  if ((str.Contains("DATA") || str.Contains("data")) &&
      !(str.Contains("MC") || str.Contains("mc"))) kind=DATA;
  else if ((str.Contains("MC") || str.Contains("mc")) &&
	   !(str.Contains("DATA") || str.Contains("data"))) kind=MC;
  else {
    std::cout << "DetermineDataKind failed at <" << str << ">\n";
    assert(0);
  }
  return kind;
}

// ------------------------------------------------------------------



// ------------------------------------------------------------------
//     Printing
// ------------------------------------------------------------------

std::ostream& operator<<(std::ostream& out, DYTools::TMassBinning_t massBinning) { out << MassBinningName(massBinning); return out; }
std::ostream& operator<<(std::ostream& out, DYTools::TSystematicsStudy_t study) { out << SystematicsStudyName(study); return out; }
std::ostream& operator<<(std::ostream& out, DYTools::TTnPMethod_t method) { out << TnPMethodName(method); return out; }
std::ostream& operator<<(std::ostream& out, DYTools::TEfficiencyKind_t kind) { out << EfficiencyKindName(kind); return out; }
std::ostream& operator<<(std::ostream& out, DYTools::TEtBinSet_t set) { out << EtBinSetName(set); return out; }
std::ostream& operator<<(std::ostream& out, DYTools::TEtaBinSet_t set) { out << EtaBinSetName(set); return out; }
std::ostream& operator<<(std::ostream& out, DYTools::TDataKind_t kind) { out << DataKindName(kind); return out; }


// ------------------------------------------------------------------



#endif
