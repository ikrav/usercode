#ifndef SAMPLEITALY_CC
#define SAMPLEITALY_CC

#include "SampleItaly.h"

#include <iostream>
using namespace std;

ClassImp(SampleItaly)

SampleItaly::SampleItaly(TString dir, TString base, TString end, int n):
  _dataDir         (dir),
  _fileNameBase    (base),
  _fileNameEnding  (end),
  _nFiles          (n),
  _category        (NotValidCategory),
  _prefix          (""),
  _suffix          ("") 
{
  _chain = new TChain();
  
  if( _nFiles > 0 ){
    
    cout << "SampleItaly: constructing TChain" << endl;
    Long_t *q1=0,*q2=0,*q3=0,*q4=0;
    for(int i = 1; i <= _nFiles; i++) {
      TString filename = _dataDir + _fileNameBase;
      filename += i;
      filename += _fileNameEnding;
      if ( ! gSystem->GetPathInfo(filename,q1,q2,q3,q4) ){
	cout << "  Adding file " << filename << endl;
	_chain->Add(filename);
      }else{
	cout << "  Skipping missing file " << filename << endl;
      }
    }
  }else{
    cout << "Created an empty sample, as requested" << endl;
  }
  return;
}

SampleItaly::~SampleItaly(){
  delete _chain;
}

void SampleItaly::initialize(){

  if( _category == NotValidCategory ) {
    setInitialized(false);
    return;
  }

  if( _category == ZGolden1or2HLT || _category == ZGolden2HLT ||
      _category == ZGolden1HLT    || _category == ZMuMuNonIso ) {
    _prefix = "goodZToMuMuEdmNtupleLoose_zGolden";
  }else if ( _category == ZMuTrk ) {
    _prefix = "goodZToMuMuOneTrackEdmNtupleLoose_zMuTrk";
  }else if ( _category == ZMuTrkMu ) {
    _prefix = "goodZToMuMuOneTrackerMuonEdmNtupleLoose_zMuTrkMu";
  }else if ( _category == ZMuSta ) {
    _prefix = "goodZToMuMuOneStandAloneEdmNtupleLoose_zMuSta";
  }else if ( _category == ZSameCharge ) {
    _prefix = "goodZToMuMuSameChargeEdmNtupleLoose_zSameCharge";
  } else {
    cout << "SampleItaly::initialize: Error: invalod category requested" << endl;
    setInitialized(false);
    return;
  }    
  _suffix = "_ZMuMuSubskim.obj";
  
  setInitialized(true);
  return;
}

#endif
