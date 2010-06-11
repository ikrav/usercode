#ifndef SAMPLEMIT_CC
#define SAMPLEMIT_CC

#include "SampleMIT.h"

#include <iostream>
using namespace std;

ClassImp(SampleMIT);
  
SampleMIT::SampleMIT():
  SampleBase(),
  _filename (""),
  _infile (0),
  _eventInfoTree (0),
  _dimuonTree (0)
{
}

SampleMIT::SampleMIT(TString filename, double xsec, TString label, int color):
  SampleBase(),
  _filename (filename),
  _infile (0),
  _eventInfoTree (0),
  _dimuonTree (0)
{
  _xsec = xsec;
  _label = label;
  _color = color;
  // Open the data file
  setInputSource();
}

SampleMIT::~SampleMIT(){};

bool SampleMIT::setInputSource(){

  bool result = false;
  if( ! _filename.IsNull() )
    _infile = new TFile(_filename);
  
  if(! _filename.IsNull() &&_infile->IsOpen()) {
    setInputSourceSetup(true);
    result = true;
  }else{
    setInputSourceSetup(false);
    if( _filename.IsNull() )
      cout << "SampleMIT: Initialized object with empty file as requested" << endl;
    else
      cout << "Failed to open file " << _filename << endl;
  }

  return result;
}


bool SampleMIT::inputAccessRefresh(){

  if(_infile &&_infile->IsOpen()){
    _infile->cd();
  } else
    return false;

  return true;
}

int SampleMIT::getNEvents(){

  if( ! inputAccessRefresh() ) {
    cout << "SampleMIT::getNEvents: error: data file could not be accessed" << endl;
    return -1;
  }

  if( ! isNtupleAccessSetup() ){
    cout << "SampleMIT::getNEvents: error: ntuple is not set up" << endl;
    return -1;
  }
  return _eventInfoTree->GetEntries();
}

int SampleMIT::getNCandidates(){

  if( ! inputAccessRefresh() )
    return -1;

//   if(_infile &&_infile->IsOpen()){
//     _infile->cd();
//   }else
//     return -1;

  if( ! isNtupleAccessSetup() ){
    cout << "SampleMIT::getNCandidates: error: ntuple is not set up" << endl;
    return -1;
  }
  return _dimuonTree->GetEntries();
}

bool SampleMIT::setNtupleAccess(){

  bool result = false;
  if( ! isInputSourceSetup() ) {
    cout << "SampleMIT::setNtupleAccess: error: can not set up, data file not initialized" << endl;
    return result;
  }

//   if( ! isInitialized() ){
//     return result;
//   }

  _eventInfoTree = (TTree*)_infile->FindObjectAny("EventInfo");
  _dimuonTree    = (TTree*)_infile->FindObjectAny("Dimuon");
  
  if( _eventInfoTree != 0 && _dimuonTree != 0 ){
    setNtupleAccessSetup(true);
  }else{
    setNtupleAccessSetup(false);
    cout << "SampleMIT::setup: error: can not find ntuples in the file" << endl;
    return result;
  }

  // Set branch address to structures that will store the info  
  _eventInfoTree  ->SetBranchAddress("EventInfo",&_eventInfo);
  _dimuonTree     ->SetBranchAddress("Dimuon",   &_dimuon);

  result = true;
  return result;
}

void SampleMIT::getCandidate(int icand){

  if( ! inputAccessRefresh() )
    return;

//   if(_infile &&_infile->IsOpen())
//     _infile->cd();
//   else
//     return;

  if( ! isNtupleAccessSetup() )
    return;

  if(_currentCandidate == icand)
    return;

  _dimuonTree->GetEntry(icand);
  if( _dimuon.eventId > _eventInfoTree->GetEntries() ){
    cout << "SampleMIT::getCandidate: error: mismatch EvengInfo vs Dimuons" << endl;
    _currentCandidate = -1;
    return;
  }

  _eventInfoTree->GetEntry(_dimuon.eventId-1);
  if( _dimuon.eventId != _eventInfo.eventId ){
    cout << "SampleMIT::getCandidate: error: sanity check failed" << endl;
    _currentCandidate = -1;
    return;
  }

  _currentCandidate = icand;
  return;
}

//
// Access variables of an ntuple entry  
//
int SampleMIT::triggerBits(int icand){

  if( !isCandidateReady(icand) )
    return -1;

  return _eventInfo.triggerBits;
}

double SampleMIT::mass(int icand){

  if( !isCandidateReady(icand) )
    return -1;

  return _dimuon.mass;
}
  
double SampleMIT::mu1pt(int icand){

  if( !isCandidateReady(icand) )
    return -1;

  return _dimuon.pt_1;
}

double SampleMIT::mu2pt(int icand){

  if( !isCandidateReady(icand) )
    return -1;

  return _dimuon.pt_2;
}

double SampleMIT::mu1eta(int icand){

  if( !isCandidateReady(icand) )
    return -1;

  return _dimuon.eta_1;
}

double SampleMIT::mu2eta(int icand){

  if( !isCandidateReady(icand) )
    return -1;

  return _dimuon.eta_2;
}

int SampleMIT::mu1hltMatchBits(int icand){

  if( !isCandidateReady(icand) )
    return -1;

  return _dimuon.hltMatchBits_1;
}

int SampleMIT::mu2hltMatchBits(int icand){

  if( !isCandidateReady(icand) )
    return -1;

  return _dimuon.hltMatchBits_2;
}

int SampleMIT::mu1typeBits(int icand){

  if( !isCandidateReady(icand) )
    return -1;

  return _dimuon.typeBits_1;
}

int SampleMIT::mu2typeBits(int icand){

  if( !isCandidateReady(icand) )
    return -1;

  return _dimuon.typeBits_2;
}

double SampleMIT::mu1trkIso(int icand){

  if( !isCandidateReady(icand) )
    return -1;

  return _dimuon.trkIso03_1;
}

double SampleMIT::mu2trkIso(int icand){

  if( !isCandidateReady(icand) )
    return -1;

  return _dimuon.trkIso03_2;
}

int SampleMIT::mu1charge(int icand){

  if( !isCandidateReady(icand) )
    return -1;

  return _dimuon.q_1;
}

int SampleMIT::mu2charge(int icand){

  if( !isCandidateReady(icand) )
    return -1;

  return _dimuon.q_2;
}

int SampleMIT::nJets(int icand){

  if( !isCandidateReady(icand) )
    return -1;

  return _eventInfo.nJets;
}

#endif
