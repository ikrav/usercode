#ifndef SELECTION_H
#define SELECTION_H

#include "TROOT.h"

enum TriggerBit
  {
    tNoTrigger = 0,
    tHLT_Mu9   = 1,
    tHLT_Mu11  = 2,
    tHLT_Mu15  = 4,
    tHLT_Jet30 = 8
  };

enum MuType 
  { 
    mNoMuon     = 0,
    mGlobal     = 1, 
    mTracker    = 2, 
    mStandAlone = 4
  };

enum MuMuType
  {
    MuMu2HLT              = 1,
    MuMu1HLT              = 2,
    MuTk                  = 3,
    MuSa                  = 4,
    MuMu1HLT1NonIsolated  = 5
  };

class Selection {
  
 public:
  Selection();

  inline bool isInitialized(){return _isInitialized;};
  inline void setInitialized(bool in){_isInitialized = in;}

  void initialize(int trigger = tHLT_Mu9,
		  double minPt = 20,
		  double maxEta1 = 3.0,
		  double maxEta2 = 3.0,
		  double minMass = 50,
		  double maxMass = 200,
		  double trackIso = 3);

  bool passBaselineSelection(int trigger,
			     double mass,
			     double pt1,
			     double pt2,
			     double eta1,
			     double eta2);
			     

  int findCategory(int mu1type,
		   int mu2type,
		   double mu1Iso,
		   double mu2Iso,
		   int mu1HLTmatch,
		   int mu2HLTmatch,
		   int mu1q,
		   int mu2q);
  
  ClassDef(Selection,1)

 private:

  bool _isInitialized;
  int    _trigger;
  double _minPt;
  double _maxEta1;
  double _maxEta2;
  double _minMass;
  double _maxMass;
  double _trackIso;

};

#endif
