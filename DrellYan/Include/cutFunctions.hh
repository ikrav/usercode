#if !defined(__CINT__) || defined(__MAKECINT__)

#include <TLorentzVector.h>         // class for Lorentz vector computations

// define classes and constants to read in ntuple
#include "TGenInfo.hh"
#include "TDielectron.hh"
#include "TElectron.hh"
#include "TPhoton.hh"
#include "DYTools.hh"
#include "EleIDCuts.hh"

#endif

Bool_t dielectronMatchedToGeneratorLevel(const mithep::TGenInfo *gen, const mithep::TDielectron *dielectron);

Bool_t electronMatchedToGeneratorLevel(const mithep::TGenInfo *gen, const mithep::TElectron *electron);

Bool_t scMatchedToGeneratorLevel(const mithep::TGenInfo *gen, const mithep::TPhoton *sc);

bool passID(const mithep::TElectron *electron);

bool isTag(const mithep::TElectron *electron, UInt_t trigger);

TString getLabel(int sample, int effType, int method,  int etBinning, int etaBinning);
