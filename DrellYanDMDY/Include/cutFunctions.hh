#if !defined(__CINT__) || defined(__MAKECINT__)

#include <TLorentzVector.h>         // class for Lorentz vector computations

// define classes and constants to read in ntuple
#include "TGenInfo.hh"
#include "TDielectron.hh"
#include "TElectron.hh"
#include "TPhoton.hh"
#include "DYTools.hh"
#include "DYToolsUI.hh"
#include "EleIDCuts.hh"
#include "TriggerSelection.hh"

#endif

Bool_t dielectronMatchedToGeneratorLevel(const mithep::TGenInfo *gen, const mithep::TDielectron *dielectron);

Bool_t electronMatchedToGeneratorLevel(const mithep::TGenInfo *gen, const mithep::TElectron *electron);

Bool_t scMatchedToGeneratorLevel(const mithep::TGenInfo *gen, const mithep::TPhoton *sc);

bool passID(const mithep::TElectron *electron, double rho);

bool isTag(const mithep::TElectron *electron, ULong_t trigger, double rho);

TString getLabel(int sample, DYTools::TEfficiencyKind_t effType, int method, DYTools::TEtBinSet_t etBinning, DYTools::TEtaBinSet_t etaBinning, const TriggerSelection &trigSet);

