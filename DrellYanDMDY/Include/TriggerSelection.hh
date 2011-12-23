#ifndef TRIGGERSELECTION_HH
#define TRIGGERSELECTION_HH

#include "../Include/EWKAnaDefs.hh"

enum TriggerConstantSet 
  {UNDEFINED,
   Full2011DatasetTriggers
  };

class TriggerSelection{
  
 public:
  TriggerSelection(TriggerConstantSet constantsSet, bool isData, int run):
    _constants(constantsSet),
    _isData(isData),
    _run(run){};

  UInt_t getEventTriggerBit(){
    UInt_t result = -1;
    if(_constants == Full2011DatasetTriggers){
      // Note: data and MC are the same
      // Note: the trigger
      //     kHLT_Ele17_CaloIdT_Calo_IsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
      // is packed into the same bit as the trigger
      //     kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL
      // the difference between the two is only in the name rearrangement
      if( _run >= 150000 && _run <= 170053)
	result = kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL;
      else if(_run >= 170054)
	result = kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL;
    } else {
      printf("TriggerSelection: ERROR unknown constants set requested\n");
    }
    return result;
  };

  UInt_t getLeadingTriggerObjectBit(){
    UInt_t result = -1;
    if(_constants == Full2011DatasetTriggers){
      // See nots in the getEventTriggerBit function
      if( _run >= 150000 && _run <= 170053)
	result = kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele1Obj;
      else if(_run >= 170054)
	result = kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele1Obj;
    } else {
      printf("TriggerSelection: ERROR unknown constants set requested\n");
    }
    return result;
  };

  UInt_t getTrailingTriggerObjectBit(){
    UInt_t result = -1;
    if(_constants == Full2011DatasetTriggers){
      // See nots in the getEventTriggerBit function
      if( _run >= 150000 && _run <= 170053)
	result = kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele2Obj;
      else if(_run >= 170054)
	result = kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele2Obj;
    } else {
      printf("TriggerSelection: ERROR unknown constants set requested\n");
    }
    return result;
  };

 private:
  TriggerConstantSet  _constants;
  bool                _isData;
  int                 _run;

};

#endif
