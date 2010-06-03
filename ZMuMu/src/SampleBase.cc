#ifndef SAMPLEBASE_CC
#define SAMPLEBASE_CC

#include "SampleBase.h"

#include <iostream>
using namespace std;

ClassImp(SampleBase);
  
SampleBase::SampleBase(TString filename):
  _filename(filename),
  _infile (0),
  _xsec (0),
  _label ("undefined"),
  _color (1),
  _weight (1.0),
  _currentCandidate (-1),
  _isInitialized (false),
  _isSetup (false)
{}

bool SampleBase::isCandidateReady(int icand){

  if(_infile &&_infile->IsOpen())
    _infile->cd();
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
