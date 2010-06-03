#ifndef SELECTION_CC
#define SELECTION_CC

#include "Selection.h"
#include <math.h>
#include <iostream>

using namespace std;

ClassImp(Selection);

Selection::Selection()
  :_isInitialized(false),
   _trigger(tNoTrigger),
   _minPt(0),
   _maxEta1(5),
   _maxEta2(5),
   _minMass(0.0),
   _maxMass(200.0),
   _trackIso(3.0)
{}

void Selection::initialize(int trigger,
			   double minPt,double maxEta1, double maxEta2,
			   double minMass, double maxMass, double trackIso)
{
  _trigger  = trigger;
  _minPt    = minPt;   
  _maxEta1  = maxEta1;  
  _maxEta2  = maxEta2;  
  _minMass  = minMass; 
  _maxMass  = maxMass; 
  _trackIso = trackIso;
  setInitialized(true);

  return;
}

bool Selection::passBaselineSelection(int trigger, double mass, double pt1, double pt2,
			   double eta1, double eta2)
{
  bool result = false;

  if( !isInitialized() )
    return result;

  if( trigger & _trigger  &&
      mass > _minMass      && mass < _maxMass      &&
      pt1 > _minPt         && pt2 > _minPt         &&
      ( (fabs(eta1) < _maxEta1 && fabs(eta2) < _maxEta2) ||
	(fabs(eta1) < _maxEta2 && fabs(eta2) < _maxEta1) )
      )
    result = true;
  
  
  return result;
}

int Selection::findCategory(int mu1type, int mu2type,
			    double mu1iso, double mu2iso,
			    int mu1HLTmatch, int mu2HLTmatch,
			    int mu1q, int mu2q)
{

  int category = -1;
  // The categories are mutually exclusive
  if(
     (mu1type & mGlobal)      && (mu2type & mGlobal)      &&
     mu1iso < _trackIso       && mu2iso < _trackIso       &&
     (mu1HLTmatch & _trigger) && (mu2HLTmatch & _trigger) &&
     mu1q*mu2q < 0 ) 
    {
      // Two isolated global muons, both HLT-matched, opposite sign
      category = 1;
    }else if (
	      (mu1type & mGlobal)        && (mu2type & mGlobal)        &&
	      mu1iso < _trackIso         && mu2iso < _trackIso         &&
	      ( (mu1HLTmatch & _trigger) || (mu2HLTmatch & _trigger) ) &&
	      mu1q*mu2q < 0 ) 
    {
      // Two isolated global muons, only one HLT matched, opposite sign
      // Note: already the 2-hlt already failed above
      category = 2;
    }else if (
	 ( ( (mu1type & mGlobal) && (mu1HLTmatch & _trigger) 
	     && (mu2type == mTracker || mu2type == mNoMuon) ) ||
	   ( (mu2type & mGlobal) && (mu2HLTmatch & _trigger) 
	     && (mu1type == mTracker || mu1type == mNoMuon) ) ) &&
	 mu1iso < _trackIso && mu2iso < _trackIso &&
	 mu1q*mu2q < 0 )
    {
      // One Global and one Tracker muon. Both isolated. Global muon HLT matched.
      // Opposite sign.
      // Note: using == or & for checking muon type is intentional. Global muons
      //  are normally both Tracker and StandAlone.
      // Note: muon type zero means it is a general track (as opposed to tracker muon).
      category = 3;
    }else if (
	      ( ( (mu1type & mGlobal) && (mu1HLTmatch & _trigger)
	     && (mu2type == mStandAlone) ) ||
		( (mu2type & mGlobal) && (mu2HLTmatch & _trigger)
	     && (mu1type == mStandAlone) ) ) &&
	 mu1iso < _trackIso && mu2iso < _trackIso &&
	 mu1q*mu2q < 0 )
    {
      // One Global and one StandAlone muon. Both isolated. Global muon HLT matched.
      // Opposite sign.
      // Note: using == or & for checking muon type is intentional. Global muons
      //  are normally both Tracker and StandAlone.
      category = 4;
    }else if (
	      (mu1type & mGlobal)        && (mu2type & mGlobal)        &&
	      (mu1iso > _trackIso || mu2iso > _trackIso )              &&
	      ( (mu1HLTmatch & _trigger) || (mu2HLTmatch & _trigger) ) &&
	      mu1q*mu2q < 0 ) 
    {
      // Two global muons. At least one HLT matched. At least one non-isolated.
      category = 5;
    }else if (
	      (mu1type & mGlobal)        && (mu2type & mGlobal)        &&
	      mu1iso < _trackIso         && mu2iso < _trackIso         &&
	      ( (mu1HLTmatch & _trigger) || (mu2HLTmatch & _trigger) ) &&
	      mu1q*mu2q > 0 ) 
    {
      // Two isolated global muons, at least one HLT matched, SAME sign
      category = 6;
    }

  return category;
}

#endif
