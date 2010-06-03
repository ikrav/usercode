#ifndef SAMPLEBASE_H
#define SAMPLEBASE_H

#include "TROOT.h"
#include "TFile.h"

class SampleBase {

 public:
  SampleBase(TString filename);

  inline void setCrossSection(double xsec){_xsec = xsec;}
  inline void setLabel       (TString label){_label = label;} 
  inline void setColor       (int color){_color = color;} 
  inline void setWeight      (double weight){_weight = weight;} 
  inline void setInitialized (bool a){_isInitialized = a;}
  inline void setSetup       (bool a){_isSetup = a;}
  
  inline double  getCrossSection(){return _xsec;} 
  inline TString getLabel       (){return _label;} 
  inline int     getColor       (){return _color;} 
  inline double  getWeight      (){return _weight;} 
  inline bool    isInitialized  (){return _isInitialized;}
  inline bool    isSetup        (){return _isSetup;}
  
  // Methods specific to an ntuple implementation
  virtual int getNEvents() = 0; // return number of MC events
  virtual int getNCandidates() = 0; // return number of candidates
  virtual int setup() = 0;       // set up ntuple access
  virtual void getCandidate(int icand) = 0; // run GetEntry on relevant trees  

  bool isCandidateReady(int icand);

  // Access variables of an ntuple entry  
  virtual int triggerBits(int icand) = 0;  
  virtual double mass(int icand) = 0;  
  virtual double mu1pt(int icand) = 0;  
  virtual double mu2pt(int icand) = 0;  
  virtual double mu1eta(int icand) = 0;  
  virtual double mu2eta(int icand) = 0;  
  virtual int mu1hltMatchBits(int icand) = 0;  
  virtual int mu2hltMatchBits(int icand) = 0;  
  virtual int mu1typeBits(int icand) = 0;  
  virtual int mu2typeBits(int icand) = 0;  
  virtual double mu1trkIso(int icand) = 0;  
  virtual double mu2trkIso(int icand) = 0;  
  virtual int mu1charge(int icand) = 0;  
  virtual int mu2charge(int icand) = 0;  

   ClassDef(SampleBase,1)

 protected:

  TString  _filename;    // filename where ntuple is stored
  TFile   *_infile;      // file with the ntuples
  double   _xsec;        // cross section
  TString  _label;       // label describing the sample, TLatex syntax ok
  int      _color;       // color of the histograms in case we draw them
  double   _weight;      // weight of this sample when several samples are added
  int      _currentCandidate; // index of the candidate presently loaded

  bool     _isInitialized; // Ntuples file(s) found and opened
  bool     _isSetup;       // ntuple TTrees found, non-null
};

#endif
