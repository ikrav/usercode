#ifndef EventSelector_HH
#define EventSelector_HH

#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <iostream>
#include <ostream>
#include <string>

#include "DYTools.hh"
#include "ElectronEnergyScale.hh"
#include "TDielectron.hh"
#include "TVertex.hh"
#include "EleIDCuts.hh"
#include "eventCounter.h"

// -------------------------------------------------

class DielectronSelector_t {
public:
  typedef enum { _selectNone, _selectDefault } TSelectionType_t;
  typedef enum { _escaleNone, _escaleUncorrected, _escaleData, _escaleDataRnd, _escaleMC } TEScaleCorrection_t;
protected:
  TSelectionType_t fSelection;
  ElectronEnergyScale *fEScale;
  UInt_t fTotalCandidates, fCandidatesGoodEta, fCandidatesGoodEt, fCandidatesHLTMatched, fCandidatesIDPassed, fCandidatesMassAboveMinLimit;
public:
  DielectronSelector_t(TSelectionType_t set_selection, ElectronEnergyScale *escale=NULL);

  bool testDielectron_default(mithep::TDielectron *dielectron, TEScaleCorrection_t applyEScale, 
			      ULong_t leadingTriggerObjectBit, ULong_t trailingTriggerObjectBit, 
			      double rho,
			      eventCounter_t *ec=NULL); // the last argument is rho for PU correction of the isolaiton

  // dielectron may be modified by escale corrections
  bool testDielectron(mithep::TDielectron *dielectron, TEScaleCorrection_t applyEScale, 
		      ULong_t leadingTriggerObjectBit, ULong_t trailingTriggerObjectBit,
		      double rho, 
		      eventCounter_t *ec=NULL) {// the last argument is rho for PU correction of the isolaiton
    bool ok=kFALSE;
    switch(fSelection) {
    case _selectNone: break;
    case _selectDefault: ok=testDielectron_default(dielectron,applyEScale,leadingTriggerObjectBit,trailingTriggerObjectBit,rho,ec); break;
    default:
      ok=kFALSE;
    }
    return ok;
  }

  static std::string selectionName(DielectronSelector_t::TSelectionType_t selection) {
    std::string name;
    switch(selection) {
    case DielectronSelector_t::_selectNone: name="None"; break;
    case DielectronSelector_t::_selectDefault: name="default"; break;
    default: name="UNKNOWN";
    }
    return name;
  }

  bool operator()(mithep::TDielectron *dielectron, TEScaleCorrection_t applyEScale, 
		  ULong_t leadingTriggerObjectBit, ULong_t trailingTriggerObjectBit,
		  double rho, eventCounter_t *ec=NULL) {
    return testDielectron(dielectron,applyEScale,
			  leadingTriggerObjectBit,trailingTriggerObjectBit,
			  rho,ec);
  }

  std::ostream& printCounts(std::ostream&);
};

// -------------------------------------------------
// -------------------------------------------------

// if this is a spec.skim file, rescale xsec
int AdjustXSectionForSkim(TFile *infile, Double_t &xsec, UInt_t numEntries, int verbatim=1) {
  if(xsec>0) { 
    // if this is a spec.skim file, rescale xsec
    TTree *descrTree=(TTree*)infile->Get("Description");
    if (descrTree) {
      UInt_t origNumEntries=0;
      descrTree->SetBranchAddress("origNumEntries",&origNumEntries);
      descrTree->GetEntry(0);
      if (origNumEntries>0) {
	Double_t factor=numEntries/double(origNumEntries);
	if (verbatim) std::cout << " -> rescaling xsec by " << factor << " due to skimming\n";
	xsec*=factor;
      }
      delete descrTree;
    }
    else {
      if (verbatim) std::cout << "descrTree not found\n";
    }
  }
  return 1;
}

// -------------------------------------------------
// -------------------------------------------------

//pvArr: TClonesArray("mithep::TVertex");
UInt_t countGoodVertices(const TClonesArray *pvArr) {
  UInt_t nGoodPV=0;
  for(Int_t ipv=0; ipv<pvArr->GetEntriesFast(); ipv++) {
    const mithep::TVertex *pv = (mithep::TVertex*)((*pvArr)[ipv]);                  
    if(pv->nTracksFit                        < 1)  continue;
    if(pv->ndof                              < 4)  continue;
    if(fabs(pv->z)                           > 24) continue;
    if(sqrt((pv->x)*(pv->x)+(pv->y)*(pv->y)) > 2)  continue;
    nGoodPV++;
  }
  return nGoodPV;
}

// -------------------------------------------------
// -------------------------------------------------


#endif
