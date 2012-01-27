#ifndef TRIGGERSELECTION_HH
#define TRIGGERSELECTION_HH

#include "../Include/EWKAnaDefs.hh"
#include <TString.h>
#include <iostream>

// -----------------------------------------
//
//  TriggerConstantSet is influenced by the L1 seeding of the trigger
//  kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL
//  In runs 170054-170759 it was singleEG (Run2011A), in runs 170826-180252 is was DoubleEG (Run2011A & Run2011B).
//
// -----------------------------------------

enum TriggerConstantSet 
  { TrigSet_UNDEFINED =0,
    Full2011DatasetTriggers =7,   //  Full2011DatasetTriggers includes all periods (1+2+4=7)
    TrigSet_2011A_SingleEG  =1, 
    TrigSet_2011A_DoubleEG  =2,
    TrigSet_2011B_DoubleEG  =4
  };

const UInt_t cFirstEvent2011B = 175770;

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
  return ts;
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
// -----------------------------------------------
// -----------------------------------------------

class TriggerSelection{
  
 public:
  TriggerSelection(TriggerConstantSet constantsSet, bool isData, int run):
    _constants(constantsSet),
    _isData(isData),
    _run(run){};

  bool validRun(UInt_t run) const {
    if (!_isData) return true;
    bool ok=false;
    switch(_constants) {
    case Full2011DatasetTriggers: ok=true; break;
    case TrigSet_2011A_SingleEG: if ((run>=160404)  // lower limit is determined by the 1st 2011 events
				     && (run<=170759)) ok=true; break;
    case TrigSet_2011A_DoubleEG: if ((run>=170826) && (run< cFirstEvent2011B)) ok=true; break;
    case TrigSet_2011B_DoubleEG: if (run>=cFirstEvent2011B) ok=true; break;
    case TrigSet_UNDEFINED: 
    default:
      ok=false;
    }
    return ok;
  }

  bool suitableDataFile(const TString &dataFileName) const {
    bool ok=false;
    switch(_constants) {
    case Full2011DatasetTriggers: ok=true; break;
    case TrigSet_2011A_SingleEG:  
      ok=dataFileName.Contains("r11a-"); 
      break; // all 4 r11a-* files contain such events
    case TrigSet_2011A_DoubleEG:  
      ok=(dataFileName.Contains("r11a-del-a05") || 
	  dataFileName.Contains("r11a-del-o03")) ? true:false;
      break;
    case TrigSet_2011B_DoubleEG:
      ok=dataFileName.Contains("r11b-");
      break;
    default:
      std::cout << "unable to verify the file <" << dataFileName << ">\n";
      assert(0);
    }
    return ok;
  }

  ULong_t getEventTriggerBit(UInt_t run=0) const {
    if (run==0) run=_run;
    ULong_t result = 0;
    // old remark:
      // Note: data and MC are the same
      // Note: the trigger
      //     kHLT_Ele17_CaloIdT_Calo_IsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
      // is packed into the same bit as the trigger
      //     kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL
      // the difference between the two is only in the name rearrangement
    if ( !_isData ) {
      result = (kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL |
		kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL);
    }
    else {
      if(validRun(run)) {
	if( run >= 150000 && run <= 170053)
	  result = kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL;
	else if (run >= 170054)
	  result = kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL;
      }
    }
    return result;
  };

  ULong_t getLeadingTriggerObjectBit(UInt_t run=0) const { // no check whether the run is ok!
    if (run==0) run=_run;
    ULong_t result = 0;
    if (!_isData) {
      result = kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele1Obj |
	kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele1Obj;
    }
    else {
      if( run >= 150000 && run <= 170053)
	result = kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele1Obj;
      else if(run >= 170054)
	result = kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele1Obj;
    }
    return result;
  };

  ULong_t getTrailingTriggerObjectBit(UInt_t run=0) const { // no check whether the run is ok!
    if (run==0) run=_run;
    ULong_t result = 0;
    if (!_isData) {
      result = kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele2Obj |
	kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele2Obj;
    }
    else {
      if( run >= 150000 && run <= 170053)
	result = kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele2Obj;
      else if(run >= 170054)
	result = kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele2Obj;
    }
    return result;
  };


  ULong_t getEventTriggerBit_SCtoGSF(UInt_t run) const {
    if (_isData && !validRun(run)) return 0UL;
    ULong_t bits=
	kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30 | 
	kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17;
      //kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30 | 
      //kHLT_Ele32_CaloIdL_CaloIsoVL_SC17                     // was defined in eff_Reco.C for 2011A(early)
    return bits;
  }

  ULong_t getLeadingTriggerObjectBit_SCtoGSF(int) const { // no check whether the run is ok!
    ULong_t bits=
	kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_EleObj | 
	kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_EleObj;
      //kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_EleObj | 
      //kHLT_Ele32_CaloIdL_CaloIsoVL_SC17_EleObj                     // was defined in eff_Reco.C for 2011A(early)
    return bits;
  }

  ULong_t getEventTriggerBit_TagProbe(UInt_t run) const {
    if (_isData && !validRun(run)) return 0UL;
    ULong_t bits=
      kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30 | 
      kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17;
    //	kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_Ele17;      <---------- unknown
      //kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30 | 
      //kHLT_Ele32_CaloIdL_CaloIsoVL_SC17                     // was defined in eff_IdHlt.C for 2011A(early)
    if (_isData && (run>=165088) && (run<=170759)) {
      bits |= kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30;
    }
    return bits;
  }

  ULong_t getLeadingTriggerObjBit_TagProbe(UInt_t run) const { // no check whether the run is ok!
    ULong_t bits=
      kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_EleObj | 
      kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_EleObj;
      //kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_Ele17_EleObj;  <------------- unknown
      //kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_EleObj | 
      //kHLT_Ele32_CaloIdL_CaloIsoVL_SC17_EleObj                 // was defined in eff_IdHlt.C for 2011A(early)
    if (_isData && (run>=165088) && (run<=170759)) {
      bits |= kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_Ele1Obj;
    }
    return bits;
  }

  ULong_t getTrailingTriggerObjBit_TagProbeTight(UInt_t run) const { // no check whether the run is ok!
    ULong_t bits = 
      kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele1Obj |
      kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele1Obj;
    if (_isData && (run>=165088) && (run<=170759)) {
      bits |= kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_Ele1Obj;
    }
    return bits;
  }

  ULong_t getTrailingTriggerObjBit_TagProbeLoose(UInt_t run) const { // no check whether the run is ok!
    ULong_t bits=
      kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele2Obj |
      kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele2Obj;
    if (_isData && (run>=165088) && (run<=170759)) {
      bits |= kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_Ele2Obj;
    }
    return bits;
  }

 private:
  TriggerConstantSet  _constants;
  bool                _isData;
  int                 _run;

};


// -----------------------------------------------

std::ostream& operator<<(std::ostream& out, TriggerConstantSet ts) {
  out << TriggerSetName(ts);
  return out;
}

// -----------------------------------------------

#endif
