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
  case _MassBins_test4: name="MassBinsTest4"; break;
  case _MassBins_Zpeak: name="MassBinsZpeak"; break;
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
  case HLT_rndTag: name="HLTrandomTag"; break;
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
  case ETBINS6alt: name="EtBins6alt"; break;
  case ETBINS6altB: name="EtBins6altB"; break;
  case ETBINS7: name="EtBins7"; break;
  case ETBINS7alt: name="EtBins7alt"; break;
  case ETBINS7altB: name="EtBins7altB"; break;
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
  case ETABINS3: name="EtaBins3"; break;
  case ETABINS3Negs: name="EtaBins3Negs"; break;
  case ETABINS5: name="EtaBins5"; break;
  case ETABINS5Negs: name="EtaBins5Negs"; break;
  case ETABINS4test: name="EtaBins4test"; break;
  case ETABINS4testNegs: name="EtaBins4testNegs"; break;
  case ETABINS4alt: name="EtaBins4alt"; break;
  case ETABINS4altNegs: name="EtaBins4altNegs"; break;
  case ETABINS5alt: name="EtaBins5alt"; break;
  case ETABINS5altNegs: name="EtaBins5altNegs"; break;
  case ETABINS8alt: name="EtaBins8alt"; break;
  case ETABINS8altNegs: name="EtaBins8altNegs"; break;
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
  else if (str.Contains("HLTrandomTag")) kind=HLT_rndTag;
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
  else if (str.Contains("ETBINS6altB") || str.Contains("EtBins6altB")) kind=ETBINS6altB;
  else if (str.Contains("ETBINS6alt") || str.Contains("EtBins6alt")) kind=ETBINS6alt;
  else if (str.Contains("ETBINS6") || str.Contains("EtBins6")) kind=ETBINS6;
  else if (str.Contains("ETBINS7altB") || str.Contains("EtBins7altB")) kind=ETBINS7altB;
  else if (str.Contains("ETBINS7alt") || str.Contains("EtBins7alt")) kind=ETBINS7alt;
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
  else if (str.Contains("ETABINS3Negs") || str.Contains("EtaBins3Negs")) kind=ETABINS3Negs;
  else if (str.Contains("ETABINS3") || str.Contains("EtaBins3")) kind=ETABINS3;
  else if (str.Contains("ETABINS5altNegs") || str.Contains("EtaBins5altNegs")) kind=ETABINS5altNegs;
  else if (str.Contains("ETABINS5alt") || str.Contains("EtaBins5alt")) kind=ETABINS5alt;
  else if (str.Contains("ETABINS5") || str.Contains("EtaBins5")) kind=ETABINS5;
  else if (str.Contains("ETABINS4altNegs") || str.Contains("EtaBins4altNegs")) kind=ETABINS4altNegs;
  else if (str.Contains("ETABINS4alt") || str.Contains("EtaBins4alt")) kind=ETABINS4alt;
  else if (str.Contains("ETABINS4testNegs") || str.Contains("EtaBins4testNegs")) kind=ETABINS4testNegs;
  else if (str.Contains("ETABINS4test") || str.Contains("EtaBins4test")) kind=ETABINS4test;
  else if (str.Contains("ETABINS8altNegs") || str.Contains("EtaBins8altNegs")) kind=ETABINS8altNegs;
  else if (str.Contains("ETABINS8alt") || str.Contains("EtaBins8alt")) kind=ETABINS8alt;
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

inline
TString CrossSectionKindName(DYTools::TCrossSectionKind_t kind) {
  using namespace DYTools;
  TString name;
  switch(kind) {
  case _cs_None: name="none"; break;
  case _cs_preFsr: name="preFsr"; break;
  case _cs_preFsrNorm: name="preFsrNorm"; break;
  case _cs_preFsrDet: name="preFsrDet"; break;
  case _cs_preFsrDetNorm: name="preFsrDetNorm"; break;
  case _cs_preFsrDetErr: name="preFsrDetErr"; break;
  case _cs_preFsrDetNormErr: name="preFsrDetNormErr"; break;
  case _cs_preFsrDetSystErr: name="preFsrDetSystErr"; break;
  case _cs_preFsrDetNormSystErr: name="preFsrDetNormSystErr"; break;
  case _cs_postFsr: name="postFsr"; break;
  case _cs_postFsrNorm: name="postFsrNorm"; break;
  case _cs_postFsrDet: name="postFsrDet"; break;
  case _cs_postFsrDetNorm: name="postFsrDetNorm"; break;
  default:
    std::cout << "CrossSectionKindName: cannot determine the name\n";
    assert(0);
  }
  return name;
}

// ------------------------------------------------------------------

inline
DYTools::TCrossSectionKind_t DetermineCrossSectionKind(const TString &str) {
  using namespace DYTools;
  DYTools::TCrossSectionKind_t kind = _cs_None;
  if (str.Contains("preFsrNorm")) kind=_cs_preFsrNorm;
  else if (str.Contains("preFsrDetNormSystErr")) kind=_cs_preFsrDetNormSystErr;
  else if (str.Contains("preFsrDetNormErr")) kind=_cs_preFsrDetNormErr;
  else if (str.Contains("preFsrDetNorm")) kind=_cs_preFsrDetNorm;
  else if (str.Contains("preFsrDetSystErr")) kind=_cs_preFsrDetSystErr;
  else if (str.Contains("preFsrDetErr")) kind=_cs_preFsrDetErr;
  else if (str.Contains("preFsrDet")) kind=_cs_preFsrDet;
  else if (str.Contains("preFsr")) kind=_cs_preFsr;
  else if (str.Contains("postFsrNorm")) kind=_cs_postFsrNorm;
  else if (str.Contains("postFsrDetNorm")) kind=_cs_postFsrDetNorm;
  else if (str.Contains("postFsrDet")) kind=_cs_postFsrDet;
  else if (str.Contains("postFsr")) kind=_cs_postFsr;
  else if (str.Contains("none") || str.Contains("None")) kind=_cs_None;
  else {
    std::cout << "DetermineCrossSectionKind: cannod identify <" << str << ">\n";
    assert(0);
  }
  return kind;
}

// ------------------------------------------------------------------

inline
TString CrossSectionKindLongName(DYTools::TCrossSectionKind_t kind) {
  using namespace DYTools;
  TString name;
  switch(kind) {
  case _cs_None: name="none"; break;
  case _cs_preFsr: name="counts in full space"; break;
  case _cs_preFsrNorm: name="normalized cross section"; break;
  case _cs_preFsrDet: name="counts in DET space"; break;
  case _cs_preFsrDetNorm: name="normalized cross section (DET)"; break;
  case _cs_preFsrDetErr: name="error in DET space"; break;
  case _cs_preFsrDetNormErr: name="error of norm.cross section (DET)"; break;
  case _cs_preFsrDetSystErr: name="syst.error in DET space"; break;
  case _cs_preFsrDetNormSystErr: name="syst.error of norm.cross section (DET)"; break;
  case _cs_postFsr: name="postFsr counts in full space"; break;
  case _cs_postFsrNorm: name="normalized postFsr cross section"; break;
  case _cs_postFsrDet: name="postFsr counts in DET space"; break;
  case _cs_postFsrDetNorm: name="normalized postFsr cross section (DET)"; break;
  default:
    std::cout << "CrossSectionKindLongName: cannot determine the name\n";
    assert(0);
  }
  return name;
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
std::ostream& operator<<(std::ostream& out, DYTools::TCrossSectionKind_t kind) { out << CrossSectionKindName(kind); return out; }


// ------------------------------------------------------------------



#endif
