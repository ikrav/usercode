#ifndef TRIGGERSELECTION_HH
#define TRIGGERSELECTION_HH

#include "../Include/EWKAnaDefs.hh"
#include "../Include/DYTools.hh"
#include <TString.h>
#include <iostream>
#include <sstream>

// -----------------------------------------
//
//  TriggerConstantSet is influenced by the L1 seeding of the trigger
//  Triggers2011::kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL
//  In runs 160329-170759 it was SingleEG (Run2011A), in runs 170826-175770 (Run2011A) 
//  and 175832- (Run2011B) it was DoubleEG
//
// -----------------------------------------

const int excludeJuly2011BadRuns=1;

enum TriggerConstantSet 
  { TrigSet_UNDEFINED =0,
    Full2011DatasetTriggers =10,   //  includes all periods (1+2+4=7)
    TrigSet_2011A_SingleEG  =1, 
    TrigSet_2011A_DoubleEG  =2,
    TrigSet_2011B_DoubleEG  =4,
    Full2012DatasetTriggers =20   //  includes all 2012 periods
  };

// Note: the HLT efficiency calc constants are used for deciding
// whether random tag and probe should be used for the periods
// of DoubleEG L1 seeds. However, there is no need to use the
// random tag and probe according to EGM POG. Thus, no new
// constants are defined for 2012.
enum HLTEfficiencyCalcDef {
  HLTEffCalc_UNDEFINED =0,
  HLTEffCalc_2011Old =1,
  HLTEffCalc_2011New =2,  
  HLTEffCalc_2011HWW =3 // not implemented
};

// run number constants, used in triggerSelection.validRun(runNo)
const UInt_t cFirstEvent2011ASingleEG = 160329;
const UInt_t cLastEvent2011ASingleEG =  170759;
const UInt_t cFirstEvent2011ADoubleEG = 170826;
const UInt_t cLastEvent2011ADoubleEG = 175770;
const UInt_t cFirstEvent2011B = 175832;

// Recorded luminosities in run sections 2011A_SingleEG, 2011A_DoubleEG, 2011B_DoubleEG
const double luminositiesInRunSections[3] = { 1.146, 1.027, 2.445 }; 

// -----------------------------------------------
//           TriggerConstantSet conversions
// -----------------------------------------------

TString TriggerSetName(TriggerConstantSet ts) {
  TString name;
  switch(ts) {
  case TrigSet_UNDEFINED:       name="UNDEFINED"; break;
  case Full2011DatasetTriggers: name="Full2011"; break;
  case TrigSet_2011A_SingleEG:  name="2011A_SingleEG"; break;
  case TrigSet_2011A_DoubleEG:  name="2011A_DoubleEG"; break;
  case TrigSet_2011B_DoubleEG:  name="2011B_DoubleEG"; break;
  case Full2012DatasetTriggers: name="Full2012"; break;
  default: name="<TriggerConstantSet name is unknown>";
  }
  return name;
}

// -----------------------------------------------

TriggerConstantSet DetermineTriggerSet(const TString &str) {
  TriggerConstantSet ts=TrigSet_UNDEFINED;
  if (str.Contains("Full2011")) ts=Full2011DatasetTriggers;
  else if (str.Contains("2011A_SingleEG")) ts=TrigSet_2011A_SingleEG;
  else if (str.Contains("2011A_DoubleEG")) ts=TrigSet_2011A_DoubleEG;
  else if (str.Contains("2011B_DoubleEG")) ts=TrigSet_2011B_DoubleEG;
  else if (str.Contains("Full2012")) ts=Full2012DatasetTriggers;
  return ts;
}

// -----------------------------------------------
// -----------------------------------------------

// No new constants are defined for 2012 data. 
// See comments for "enum HLTEfficiencyCalcDef" above.
TString HLTEfficiencyCalcName(HLTEfficiencyCalcDef ecd) {
  TString name;
  switch(ecd) {
  case HLTEffCalc_UNDEFINED: name="hltEffUndefined"; break;
  case HLTEffCalc_2011Old: name="hltEffOld"; break;
  case HLTEffCalc_2011New: name="hltEffNew"; break;
  case HLTEffCalc_2011HWW: name="hltEffHWW"; break;
  default: name="<HLTEfficiencyCalcName is unknown>";
  }
  return name;
}

// -----------------------------------------------

// No code changes are needed here for 2012 data.
// See comments for "enum HLTEfficiencyCalcDef" above.
HLTEfficiencyCalcDef DetermineHLTEfficiencyCalc(const TString &str) {
  HLTEfficiencyCalcDef ec=HLTEffCalc_UNDEFINED;
  if (str.Contains("HWW")) ec=HLTEffCalc_2011HWW;
  else if (str.Contains("hltEffNew")) ec=HLTEffCalc_2011New;
  else if (DetermineTriggerSet(str)!=TrigSet_UNDEFINED) ec=HLTEffCalc_2011Old;
  return ec;
}

// -----------------------------------------------
// -----------------------------------------------
// -----------------------------------------------

template<class T>
void PrintBits(T n) {
  std::cout << "PrintBits(" << n << ") = ";
  int first=1;
  for (unsigned int i=0; i<8*sizeof(T); ++i) {
    T t=(T(1)<<i);
    if ((t&n) != 0) {
      if (first) first=0; else std::cout << ", ";
      std::cout << i;
    }
  }
  std::cout << "\n";
  return;
}

// -----------------------------------------------

template<class T>
std::string PrintBitsToStdString(T n, const char *msg=NULL, int compact=0) {
  std::stringstream ss(std::stringstream::in | std::stringstream::out);
  if (compact) {
    if (msg) ss << msg;
  }
  else {
    ss << "PrintBits(";
    if (msg) ss << msg << ", ";
    ss << n << ") = ";
  }
  int first=1;
  for (unsigned int i=0; i<8*sizeof(T); ++i) {
    T t=(T(1)<<i);
    if ((t&n) != 0) {
      if (first) first=0; 
      else {
	if (compact) ss << "_";
	else ss << ", ";
      }
      ss << i;
    }
  }
  if (!compact) ss << "\n";
  std::string str;
  getline(ss,str);
  return str;
}

// -----------------------------------------------
// -----------------------------------------------
// -----------------------------------------------

class TriggerSelection{
  
 public:
  TriggerSelection(TriggerConstantSet constantsSet, bool isData, int run, HLTEfficiencyCalcDef hltEffCalc=HLTEffCalc_2011Old):
    _constants(constantsSet),
    _isData(isData),
    _run(run),
    _hltEffCalcAlgo(hltEffCalc)
  {}

  TriggerSelection(const TString& constantsSetString, bool isData, int run):
    _constants(DetermineTriggerSet(constantsSetString)),
    _isData(isData),
    _run(run),
    _hltEffCalcAlgo(DetermineHLTEfficiencyCalc(constantsSetString))
  {}

  TriggerSelection(const TriggerSelection &ts) :
    _constants(ts._constants), _isData(ts._isData), _run(ts._run), _hltEffCalcAlgo(ts._hltEffCalcAlgo)
  {}

  // Access
  TriggerConstantSet triggerSet() const { return _constants; }
  void triggerSet(TriggerConstantSet ts) { _constants=ts; }
  bool actOnData() const { return _isData; }
  void actOnData(bool act_on_data) { _isData = act_on_data; }
  HLTEfficiencyCalcDef hltEffCalcMethod() const { return _hltEffCalcAlgo; }
  void hltEffCalcMethod(HLTEfficiencyCalcDef hltEffCalc) {  _hltEffCalcAlgo = hltEffCalc; }
  bool isDefined() const { return (_constants != TrigSet_UNDEFINED) ? true : false; }
  bool hltEffMethodIsDefined() const { return (_hltEffCalcAlgo != HLTEffCalc_UNDEFINED) ? true : false; }
  bool hltEffMethodIs2011New() const { return (_hltEffCalcAlgo == HLTEffCalc_2011New) ? true : false; }
  //bool hltEffMethodIsHWW() const { return (_hltEffCalcAlgo == HLTEffCalc_2011HWW) ? true : false; }
  TString triggerSetName() const { return TriggerSetName(_constants); }
  TString triggerConditionsName() const { 
    TString name = TriggerSetName(_constants) + TString("_") + HLTEfficiencyCalcName(_hltEffCalcAlgo); 
    return name;
  }
  TString hltEffCalcName() const { return HLTEfficiencyCalcName(_hltEffCalcAlgo); }

  // Filtering
  bool validRun(UInt_t run) const {
    if (!_isData) return true;
    if (excludeJuly2011BadRuns && (run>=171050) && (run<=171578)) return false;
    bool ok=false;
    switch(_constants) {
    case Full2011DatasetTriggers: ok=true; break;
    case TrigSet_2011A_SingleEG: if ((run>=cFirstEvent2011ASingleEG) && (run<=cLastEvent2011ASingleEG)) ok=true; break;
    case TrigSet_2011A_DoubleEG: if ((run>=cFirstEvent2011ADoubleEG) && (run<=cLastEvent2011ADoubleEG)) ok=true; break;
    case TrigSet_2011B_DoubleEG: if (run>=cFirstEvent2011B) ok=true; break;
    case TrigSet_UNDEFINED: 
    case Full2012DatasetTriggers: ok=true; break;
    default:
      ok=false;
    }
    return ok;
  }

  bool validRunRange(UInt_t runNumMin, UInt_t runNumMax) const {
    if (!_isData) return true;
    if (runNumMin>runNumMax) {
      std::cout << "TriggerSelection::validRunRange (runNumMin>runNumMax)\n";
      throw 2;
    }
    UInt_t chkMin=0, chkMax=0;
    switch(_constants) {
    case Full2011DatasetTriggers: chkMin=0; chkMax=UInt_t(1e9); break;
    case TrigSet_2011A_SingleEG: chkMin=cFirstEvent2011ASingleEG; chkMax=cLastEvent2011ASingleEG; break;
    case TrigSet_2011A_DoubleEG: chkMin=cFirstEvent2011ADoubleEG; chkMax=cLastEvent2011ADoubleEG; break;
    case TrigSet_2011B_DoubleEG: chkMin=cFirstEvent2011B; chkMax=UInt_t(1e9); break;
    case Full2012DatasetTriggers: chkMin=0; chkMax=UInt_t(1e9); break;
    case TrigSet_UNDEFINED:       
    default:
      std::cout << "error in TriggerSelection::validRunRange";
      throw 2;
    }
    return ((runNumMax<chkMin) || (runNumMin>chkMax)) ? false : true;
  }

  // Note: Random tag and probe, according to EGM POG, is not
  // needed and usually not used at present.
  //
  bool useRandomTagTnPMethod(UInt_t run) const {
    // For "new" MC use only random tag and probe since Fall11 MC is DoubleEG 
    // for the signal trigger.
    //
    if (!_isData && !(_hltEffCalcAlgo==HLTEffCalc_2011Old) ) {
      std::cout << "trig: randomTag\n";
      return true;
    }
    //
    // If it is data but "old" method is requested do not use random TnP
    // If it is MC and old method is requested, do not use TnP
    if (!_isData || (_hltEffCalcAlgo==HLTEffCalc_2011Old)) {
      std::cout << "trig: not randomTag\n"; 
      return false;
    }
    bool yes=false;
    switch ( _hltEffCalcAlgo ) {
    case HLTEffCalc_2011New:
    case HLTEffCalc_2011HWW:
      if (run>=cFirstEvent2011ADoubleEG) yes=true;
      break;
    default:
      yes=false;
    }
    std::cout << "trig: return randomTag=" << yes << "\n";
    return yes;
  }

  // hand-made catch-all event trigger
  ULong_t getCombinedEventTriggerBit() {
    bool keepIsData=_isData;
    ULong_t bits=0UL;
    _isData=false;
    // Accumulate all MC triggers of interest
    bits |= this->getEventTriggerBit(0);
    bits |= this->getEventTriggerBit_SCtoGSF(0);
    bits |= this->getEventTriggerBit_TagProbe(0,true);
    bits |= this->getEventTriggerBit_TagProbe(0,false);
    _isData=true;
    // Accumulate all data triggers of interest
    if( DYTools::energy8TeV == 1){
      // Trigger bits are same for all runs for 8 TeV data
      bits |= this->getEventTriggerBit(0);
      bits |= this->getEventTriggerBit_SCtoGSF(0);
      bits |= this->getEventTriggerBit_TagProbe(0,true);
      bits |= this->getEventTriggerBit_TagProbe(0,false);
    }else{
      bits |= this->getEventTriggerBit(150000+1);
      bits |= this->getEventTriggerBit(170054+1);
      bits |= this->getEventTriggerBit_SCtoGSF(0);
      bits |= this->getEventTriggerBit_TagProbe(0,true);
      bits |= this->getEventTriggerBit_TagProbe(165088+1,true);
      bits |= this->getEventTriggerBit_TagProbe(0,false);
      bits |= this->getEventTriggerBit_TagProbe(165088+1,false);
    }
    _isData=keepIsData;
    return bits;
  }

  // Trigger bits: main analysis

  ULong_t getEventTriggerBit(UInt_t run) const {
    if (run==0) run=_run;
    ULong_t result = 0;

    // Triggers for 2011 (7 TeV) and 2012 (8 TeV) are different
    if( DYTools::energy8TeV == 1){
      // 8 TeV triggers
      result = Triggers2012::kHLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
    }else{
      // 7 TeV triggers
      //
      // -- old remark:
      // Note: data and MC are the same
      // Note: the trigger
      //     Triggers2011::kHLT_Ele17_CaloIdT_Calo_IsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
      // is packed into the same bit as the trigger
      //     Triggers2011::kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL
      // the difference between the two is only in the name rearrangement
      // -- end of old remark
      if ( !_isData ) {
	result = (Triggers2011::kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL |
		  Triggers2011::kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL);
      }
      else {
	if(validRun(run)) {
	  if( run >= 150000 && run <= 170053)
	    result = Triggers2011::kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL;
	  else if (run >= 170054)
	    result = Triggers2011::kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL;
	}
      }
    } // end if 7 or 8 TeV
    return result;
  };

  ULong_t getLeadingTriggerObjectBit(UInt_t run) const { // no check whether the run is ok!
    if (run==0) run=_run;
    ULong_t result = 0;

    // Triggers for 2011 (7 TeV) and 2012 (8 TeV) are different
    if( DYTools::energy8TeV == 1){
      // 8 TeV triggers
      result = Triggers2012::kHLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele1Obj;
    }else{
      // 7 TeV triggers
      if (!_isData) {
	result = Triggers2011::kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele1Obj |
	  Triggers2011::kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele1Obj;
      }
      else {
	if( run >= 150000 && run <= 170053)
	  result = Triggers2011::kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele1Obj;
	else if(run >= 170054)
	  result = Triggers2011::kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele1Obj;
      }
    } // end if 7 or 8 TeV
    return result;
  };

  ULong_t getTrailingTriggerObjectBit(UInt_t run) const { // no check whether the run is ok!
    if (run==0) run=_run;
    ULong_t result = 0;

    // Triggers for 2011 (7 TeV) and 2012 (8 TeV) are different
    if( DYTools::energy8TeV == 1){
      // 8 TeV triggers
      result = Triggers2012::kHLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele2Obj;
    }else{
      // 7 TeV triggers
      if (!_isData) {
	result = Triggers2011::kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele2Obj |
	  Triggers2011::kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele2Obj;
      }
      else {
	if( run >= 150000 && run <= 170053)
	  result = Triggers2011::kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele2Obj;
	else if(run >= 170054)
	  result = Triggers2011::kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele2Obj;
      }
    } // end if 7 or 8 TeV
    return result;
  };

  // Trigger bits for Tag&Probe analysis

  ULong_t getEventTriggerBit_SCtoGSF(UInt_t run) const {
    if (_isData && !validRun(run)) return 0UL;
    ULong_t bits=0UL;
   
    // Triggers for 2011 (7 TeV) and 2012 (8 TeV) are different
    if( DYTools::energy8TeV == 1){
      // 8 TeV triggers
      bits = Triggers2012::kHLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50;
    }else{
      // 7 TeV triggers
      bits = Triggers2011::kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30;
// 	Triggers2011::kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30 | 
//         Triggers2011::kHLT_Ele32_CaloIdL_CaloIsoVL_SC17  |                // <--- added from eff_Reco.C
//         Triggers2011::kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17;
//
      //Triggers2011::kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30 | 
      //Triggers2011::kHLT_Ele32_CaloIdL_CaloIsoVL_SC17                     // was defined in eff_Reco.C for 2011A(early)

    } // end 7 or 8 TeV
    return bits;
  }

  ULong_t getLeadingTriggerObjectBit_SCtoGSF(int) const { // no check whether the run is ok!
    ULong_t bits=0UL;

    // Triggers for 2011 (7 TeV) and 2012 (8 TeV) are different
    if( DYTools::energy8TeV == 1){
      // 8 TeV triggers
      bits = Triggers2012::kHLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50_EleObj;
    }else{
      // 7 TeV triggers
      
      bits = Triggers2011::kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_EleObj;
// 	Triggers2011::kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_EleObj | 
//         Triggers2011::kHLT_Ele32_CaloIdL_CaloIsoVL_SC17_EleObj   |             //   <--- added from eff_Reco.C
//     	Triggers2011::kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_EleObj;
//
      //Triggers2011::kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_EleObj | 
      //Triggers2011::kHLT_Ele32_CaloIdL_CaloIsoVL_SC17_EleObj                     // was defined in eff_Reco.C for 2011A(early)

    } // end 7 or 8 TeV
    return bits;
  }

  ULong_t getEventTriggerBit_TagProbe(UInt_t run, bool idEffTrigger) const {
    if (_isData && !validRun(run)) return 0UL;
    ULong_t bits=0UL;

    // Triggers for 2011 (7 TeV) and 2012 (8 TeV) are different
    if( DYTools::energy8TeV == 1){
      // 8 TeV triggers
      bits = Triggers2012::kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50
	| Triggers2012::kHLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50;
    }else{
      // 7 TeV triggers
      bits = Triggers2011::kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30;
      //
//       Triggers2011::kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30 | 
//       Triggers2011::kHLT_Ele32_CaloIdL_CaloIsoVL_SC17 |             // <---- added from eff_IdHlt.C
//       Triggers2011::kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17;
      //
      //	Triggers2011::kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_Ele17;      <---------- unknown (Jan 26, 2012)
      //Triggers2011::kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30 | 
      //Triggers2011::kHLT_Ele32_CaloIdL_CaloIsoVL_SC17                     // was defined in eff_IdHlt.C for 2011A(early)
      if (_isData && (((run>=165088) && (run<=170759)) || idEffTrigger)) {
	bits |= Triggers2011::kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30;
      }
    } // end if 7 or 8 TeV
    return bits;
  }

  ULong_t getTagTriggerObjBit(UInt_t run, bool idEffTrigger) const { // no check whether the run is ok!
    ULong_t bits=0UL;

    // Triggers for 2011 (7 TeV) and 2012 (8 TeV) are different
    if( DYTools::energy8TeV == 1){
      // 8 TeV triggers
      bits = Triggers2012::kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50_Ele1Obj
	| Triggers2012::kHLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50_EleObj;
    }else{
      // 7 TeV triggers
      bits = Triggers2011::kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_EleObj;
      //Triggers2011::kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_EleObj;
      //Triggers2011::kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_Ele17_EleObj;  <------------- unknown (Jan 26, 2012)
      //Triggers2011::kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_EleObj | 
      //Triggers2011::kHLT_Ele32_CaloIdL_CaloIsoVL_SC17_EleObj                 // was defined in eff_IdHlt.C for 2011A(early)
      //bits |= Triggers2011::kHLT_Ele32_CaloIdL_CaloIsoVL_SC17_EleObj;   // was defined in ieff_idHlt.C
      if (_isData && (((run>=165088) && (run<=170759)) || idEffTrigger)) {
	bits |= Triggers2011::kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_Ele1Obj;
      }
    } // end if 7 or 8 TeV
    return bits;
  }

  ULong_t getProbeTriggerObjBit_Tight(UInt_t run, bool idEffTrigger) const { // no check whether the run is ok!
    if (0) std::cout << "getProbeTriggerObjBit_Tight(" << run << "," <<  idEffTrigger << ")\n";

    // The probe trigger object bit must match the kind for which
    // we are going to measure the trigger efficiency, i.e. the signal
    // trigger.
    ULong_t bits= 0UL;

    // Triggers for 2011 (7 TeV) and 2012 (8 TeV) are different
    if( DYTools::energy8TeV == 1){
      // 8 TeV triggers
      bits = Triggers2012::kHLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele1Obj;
    }else{
      // 7 TeV triggers

      bits = Triggers2011::kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele1Obj |
	Triggers2011::kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele1Obj;
      
      // THE CODE BELOW IS COMMENTED OUT. The code below appears a mistake.
      // We want to get the bit of the signal trigger in this function, not
      // of the tag and probe trigger.
//     if (_isData && (((run>=165088) && (run<=170759)) || idEffTrigger)) {
//       // Triggers2011::kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30 is DoubleEG in Fall11 MC
//       bits |= Triggers2011::kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_Ele1Obj;
//     }
    } // end if 7 or 8 TeV
    return bits;
  }

  ULong_t getProbeTriggerObjBit_Loose(UInt_t run, bool idEffTrigger) const { // no check whether the run is ok!
    if (0) std::cout << "getProbeTriggerObjBit_Loose(" << run << "," <<  idEffTrigger << ")\n";

    // The probe trigger object bit must match the kind for which
    // we are going to measure the trigger efficiency, i.e. the signal
    // trigger.
    ULong_t bits= 0UL;

    // Triggers for 2011 (7 TeV) and 2012 (8 TeV) are different
    if( DYTools::energy8TeV == 1){
      // 8 TeV triggers
      bits = Triggers2012::kHLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele2Obj;
    }else{
      // 7 TeV triggers
      bits = Triggers2011::kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele2Obj |
	Triggers2011::kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele2Obj;

      // THE CODE BELOW IS COMMENTED OUT. The code below appears a mistake.
      // We want to get the bit of the signal trigger in this function, not
      // of the tag and probe trigger.
//    if (_isData && (((run>=165088) && (run<=170759)) || idEffTrigger)) {
//       // Triggers2011::kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30 is DoubleEG in Fall11 MC
//      bits |= Triggers2011::kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_Ele2Obj;
//     }
    } // end if 7 or 8 TeV
    return bits;
  }

  // 
  // The functions that implement cuts on trigger bits
  //

  Bool_t matchEventTriggerBit(ULong_t bit, UInt_t run) const {
    if( ! (bit & getEventTriggerBit(run) ) )
      return kFALSE;
    return kTRUE;
  }

  Bool_t matchLeadingTriggerObjectBit(ULong_t bit, UInt_t run) const { 
    if( ! (bit & getLeadingTriggerObjectBit(run) ) )
      return kFALSE;
    return kTRUE;    
  }

  Bool_t matchTrailingTriggerObjectBit(ULong_t bit, UInt_t run) const { 
    if( ! (bit & getTrailingTriggerObjectBit(run) ) )
      return kFALSE;
    return kTRUE;    
  }

  Bool_t matchEventTriggerBit_SCtoGSF(ULong_t bit, UInt_t run) const {
    if( ! (bit & getEventTriggerBit_SCtoGSF(run) ) )
      return kFALSE;
    return kTRUE;    
  }

  Bool_t matchLeadingTriggerObjectBit_SCtoGSF(ULong_t bit, UInt_t run) const { 
    if( ! (bit & getLeadingTriggerObjectBit_SCtoGSF(run) ) )
      return kFALSE;
    return kTRUE;    
  }

  Bool_t matchEventTriggerBit_TagProbe(ULong_t bit, UInt_t run, bool idEffTrigger) const {
    if( ! (bit & getEventTriggerBit_TagProbe(run,idEffTrigger) ) )
      return kFALSE;
    return kTRUE;    
  }

  Bool_t matchTagTriggerObjBit(ULong_t bit, UInt_t run, bool idEffTrigger) const { 
    if( ! (bit & getTagTriggerObjBit(run,idEffTrigger) ) )
      return kFALSE;
    return kTRUE;    
  }


  Bool_t matchProbeTriggerObjBit_Tight(ULong_t bit, UInt_t run, bool idEffTrigger) const {
    if( ! (bit & getProbeTriggerObjBit_Tight(run,idEffTrigger) ) )
      return kFALSE;
    return kTRUE;    
  }

  Bool_t matchProbeTriggerObjBit_Loose(ULong_t bit, UInt_t run, bool idEffTrigger) const { 
    if( ! (bit & getProbeTriggerObjBit_Loose(run,idEffTrigger) ) )
      return kFALSE;
    return kTRUE;    
  }

  Bool_t matchTwoTriggerObjectsAnyOrder(ULong_t bit1, ULong_t bit2, UInt_t run) const {
    if( ! ( 
	   (matchLeadingTriggerObjectBit(bit1, run) &&
	    matchTrailingTriggerObjectBit(bit2, run) )
	   ||
	   (matchLeadingTriggerObjectBit(bit2, run) &&
	    matchTrailingTriggerObjectBit(bit1, run) ) ) )
      return kFALSE;
    return kTRUE;
  }
 

 private:
  TriggerConstantSet  _constants;
  bool                _isData;
  int                 _run;             // this is an obsolete data member
  HLTEfficiencyCalcDef   _hltEffCalcAlgo;

};


// -----------------------------------------------

std::ostream& operator<<(std::ostream& out, TriggerConstantSet ts) {
  out << TriggerSetName(ts);
  return out;
}

// -----------------------------------------------

#endif
