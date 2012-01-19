#ifndef TRIGGERSELECTION_HH
#define TRIGGERSELECTION_HH

#include "../Include/EWKAnaDefs.hh"
#include <iostream>

enum TriggerConstantSet 
  {UNDEFINED,
   Full2011DatasetTriggers
  };

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

class TriggerSelection{
  
 public:
  TriggerSelection(TriggerConstantSet constantsSet, bool isData, int run):
    _constants(constantsSet),
    _isData(isData),
    _run(run){};

  ULong_t getEventTriggerBit(){
    ULong_t result = 0;
    if(_constants == Full2011DatasetTriggers){
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
      if( _run >= 150000 && _run <= 170053)
	result = kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL;
      else if (_run >= 170054)
	result = kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL;
      }
    } else {
      printf("TriggerSelection: ERROR unknown constants set requested\n");
    }
    return result;
  };

  ULong_t getLeadingTriggerObjectBit(){
    ULong_t result = 0;
    if(_constants == Full2011DatasetTriggers){
      // See nots in the getEventTriggerBit function
      if (!_isData) {
	result = kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele1Obj |
	  kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele1Obj;
      }
      else {
	if( _run >= 150000 && _run <= 170053)
	  result = kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele1Obj;
	else if(_run >= 170054)
	  result = kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele1Obj;
      }
    } else {
      printf("TriggerSelection: ERROR unknown constants set requested\n");
    }
    return result;
  };

  ULong_t getTrailingTriggerObjectBit(){
    ULong_t result = 0;
    if(_constants == Full2011DatasetTriggers){
      // See nots in the getEventTriggerBit function
      if (!_isData) {
	result = kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele2Obj |
	  kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele2Obj;
      }
      else {
	if( _run >= 150000 && _run <= 170053)
	  result = kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele2Obj;
	else if(_run >= 170054)
	result = kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele2Obj;
      }
    }
    else {
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
