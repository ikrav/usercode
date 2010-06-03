#ifndef SAMPLEMIT_CC
#define SAMPLEMIT_CC

#include "SampleMIT.h"

#include <iostream>
using namespace std;

ClassImp(SampleMIT);
  
SampleMIT::SampleMIT(TString filename):
  SampleBase(filename),
  _eventInfoTree (0),
  _dimuonTree (0)
{
  initializeFile(filename);
}

SampleMIT::SampleMIT(TString filename, double xsec, TString label, int color):
  SampleBase(filename),
  _eventInfoTree (0),
  _dimuonTree (0)
{
  initializeFile(filename);
  _xsec = xsec;
  _label = label;
  _color = color;
}

void SampleMIT::initializeFile(TString filename){
  if(!filename.IsNull())
    _infile = new TFile(filename);
  
  if(!filename.IsNull() &&_infile->IsOpen())
    setInitialized(true);
  else{
    setInitialized(false);
    if(filename.IsNull())
      cout << "SampleMIT: Initialized object with empty file as requested" << endl;
    else
      cout << "Failed to open file " << filename << endl;
  }

  return;
}

int SampleMIT::getNEvents(){

  if(_infile &&_infile->IsOpen()){
    _infile->cd();
  }else
    return -1;

  if( ! isSetup() ){
    cout << "SampleMIT::getNEvents: error: ntuple is not set up" << endl;
    return -1;
  }
  return _eventInfoTree->GetEntries();
}

int SampleMIT::getNCandidates(){

  if(_infile &&_infile->IsOpen()){
    _infile->cd();
  }else
    return -1;

  if( ! isSetup() ){
    cout << "SampleMIT::getNCandidates: error: ntuple is not set up" << endl;
    return -1;
  }
  return _dimuonTree->GetEntries();
}

int SampleMIT::setup(){

  int status = 0;
  if(_infile && _infile->IsOpen())
    _infile->cd();
  else 
    return status;

  if( ! isInitialized() ){
    cout << "SampleMIT::setup: error: can not set up, not initialized" << endl;
    return status;
  }

  _eventInfoTree = (TTree*)_infile->FindObjectAny("EventInfo");
  _dimuonTree    = (TTree*)_infile->FindObjectAny("Dimuon");
  
  if( _eventInfoTree != 0 && _dimuonTree != 0 ){
    setSetup(true);
  }else{
    setSetup(false);
    cout << "SampleMIT::setup: error: can not find ntuples in the file" << endl;
    return status;
  }

  // Set branch address to structures that will store the info  
  _eventInfoTree  ->SetBranchAddress("EventInfo",&_eventInfo);
  _dimuonTree     ->SetBranchAddress("Dimuon",   &_dimuon);

  status=1;
  return status;
}

void SampleMIT::getCandidate(int icand){

  if(_infile &&_infile->IsOpen())
    _infile->cd();
  else
    return;

  if( ! isSetup() )
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
