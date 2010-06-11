#ifndef SAMPLEBASE_H
#define SAMPLEBASE_H

#include "TROOT.h"
#include "TFile.h"

enum CategoriesList {
  NotValidCategory = -1,
  ZGolden1or2HLT   = 0,
  ZGolden2HLT      = 1,
  ZGolden1HLT      = 2,
  ZMuTrk           = 3,
  ZMuTrkMu         = 4,
  ZMuSta           = 5,
  ZMuMuNonIso      = 6,
  ZSameCharge      = 7
};

class SampleBase {

 public:
  SampleBase();
  ~SampleBase(){}; 
/*   SampleBase(TString filename); */
  

  inline void setCrossSection(double xsec){_xsec = xsec;}
  inline void setLabel       (TString label){_label = label;} 
  inline void setColor       (int color){_color = color;} 
  inline void setWeight      (double weight){_weight = weight;} 
/*   inline void setInitialized (bool a){_isInitialized = a;} */
/*   inline void setSetup       (bool a){_isSetup = a;} */
  inline void setCategory    (int cat){_category = cat;};
  inline void setInputSourceSetup(bool a){_isInputSourceSetup = a;};
  inline void setNtupleAccessSetup(bool a){_isNtupleAccessSetup = a;};
  
  inline double  getCrossSection(){return _xsec;} 
  inline TString getLabel       (){return _label;} 
  inline int     getColor       (){return _color;} 
  inline double  getWeight      (){return _weight;} 
/*   inline bool    isInitialized  (){return _isInitialized;} */
/*   inline bool    isSetup        (){return _isSetup;} */
  inline bool    isInputSourceSetup  (){return _isInputSourceSetup;}
  inline bool    isNtupleAccessSetup (){return _isNtupleAccessSetup;}
  inline int     getCategory    (){return _category;};
  
  // Set up input. Both of the functions have to be defined
  // in the derived classes. The idea is that the derived
  // class will have a function setInputSourceXXX(relevant pars)
  // that sets class variables, and calls setInputSource() to
  // actually perform the operation. Simular for ntuple access.
  //
  //  This function defines data file or files where
  // the ntuples are stored, and checks that those exist.
  virtual bool setInputSource() = 0;
  //  This function sets up ntuple ttree or tchain, and
  // defines tree branches.
  virtual bool setNtupleAccess() = 0;
  // In case of multiple files or multiple tchains when several
  // sumples are in use in the same root session, one
  // sometimes needs to "cd" into the file of this sample, or 
  // otherwise pull it to the surface.
  virtual bool inputAccessRefresh() = 0;

  // Methods specific to an ntuple implementation
  virtual int getNEvents() = 0; // return number of MC events
  virtual int getNCandidates() = 0; // return number of candidates
/*   virtual int setup() = 0;       // set up ntuple access */
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

  double   _xsec;        // cross section
  TString  _label;       // label describing the sample, TLatex syntax ok
  int      _color;       // color of the histograms in case we draw them
  double   _weight;      // weight of this sample when several samples are added
  int      _currentCandidate; // index of the candidate presently loaded

  int      _category;

  bool     _isInputSourceSetup; // Data file(s) found and opened
  bool     _isNtupleAccessSetup; // Ntuple access set up

};

#endif
