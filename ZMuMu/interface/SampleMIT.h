#ifndef SAMPLEMIT_H
#define SAMPLEMIT_H

#include "SampleBase.h"
// This header file describes data structures specific to
// MIT ntuples
#include "ZAnaStructDefs.hh"

#include "TROOT.h"
#include "TTree.h"

class SampleMIT : public SampleBase {

 public:
  
  SampleMIT(TString filename);
  SampleMIT(TString filename, double xsec, TString label, int color);

  void initializeFile(TString filename);

  // Methods specific to an ntuple implementation
  int getNEvents();       // return number of MC events
  int getNCandidates();   // return number of candidates
  int setup();       // set up ntuple access
  void getCandidate(int icand); // run GetEntry on relevant trees
  // Access variables of an ntuple entry
  int triggerBits(int icand);
  double mass(int icand);
  double mu1pt(int icand);
  double mu2pt(int icand);
  double mu1eta(int icand);
  double mu2eta(int icand);
  int mu1hltMatchBits(int icand);
  int mu2hltMatchBits(int icand);
  int mu1typeBits(int icand);
  int mu2typeBits(int icand);
  double mu1trkIso(int icand);
  double mu2trkIso(int icand);
  int mu1charge(int icand);  
  int mu2charge(int icand);  

  int nJets(int icand);
  
  ClassDef(SampleMIT,1)


 protected:

  TTree *_eventInfoTree;
  TTree *_dimuonTree;

  // The structures below are containers for SetBranchAddress
  TEventInfo _eventInfo;
  TDimuon    _dimuon;

};

#endif
