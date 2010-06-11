#ifndef SAMPLEBASE_CC
#define SAMPLEBASE_CC

#include "SampleBase.h"
#include "Selection.h"

#include <iostream>
using namespace std;

ClassImp(SampleBase);
  
SampleBase::SampleBase():
  _xsec (0),
  _label ("undefined"),
  _color (1),
  _weight (1.0),
  _currentCandidate (-1),
  _category (NotValidCategory),
  _isInputSourceSetup(false),
  _isNtupleAccessSetup(false)
{}

bool SampleBase::isCandidateReady(int icand){

  if( isNtupleAccessSetup() )
    inputAccessRefresh();
  else
    return false;

  if( _currentCandidate != icand){
    getCandidate(icand);
    if( _currentCandidate != icand )
      return false;
  }
  return true;
}

#endif
