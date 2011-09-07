#ifndef ZHtoEEbbCuts_HH
#define ZHtoEEbbCuts_HH

namespace ZHtoEEbbCuts {

  // Electron cut values
  const double cutEleETMin  = 20.0;
  const double cutEleEtaMax = 2.5;
  const double cutEleIsoMax = 0.1;
  
  // Z cut values
  const double cutZMassMin  = 75.0;
  const double cutZMassMax  = 105.0;
  const double cutZPtMin    = 100.0;
  
  // Jet cut values
  const double cutJetPTMin  = 20.0;
  const double cutJetEtaMax = 2.5;
  const double cutJetTrackCountMin = 2;
  const double cutEMFractionMin = 0.01;
  const double cutHadFractionMin = 0.01;
  const double cutEMFractionMax = 0.99;
  const double cutHadFractionMax = 0.99;
  const double csv1cut = 0.898; // CSV Tight
  const double csv2cut = 0.50;
  
  // Higgs (di-b-jet) cut values)
  // Mass cuts for 115 GeV Higgs
  const double cutHMassMin =  95.0;
  const double cutHMassMax = 125.0;
  const double cutHPtMin   = 100.0;
  
  // Other
  const double cutDRJetLeptonMin = 0.3;
  const double cutDPhiZHMin = 2.90;
  const int cutAdditionalCentralJetsMax = 1;
}

#endif
