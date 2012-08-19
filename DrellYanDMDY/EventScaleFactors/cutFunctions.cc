#if !defined(__CINT__) || defined(__MAKECINT__)
#include "../Include/cutFunctions.hh"
#endif

Bool_t dielectronMatchedToGeneratorLevel(const mithep::TGenInfo *gen, const mithep::TDielectron *dielectron){

  Bool_t result = kTRUE;

  // Matching is done as prescribed in AN-12-116
  //   - dR < 0.2
  //   - match to status 1 (post-FSR)
  //   - no charge matching

  // In the generator branch of this ntuple, first particle is always
  // negative, and second always positive. In the Dielectron block
  // of the ntuple, the first particle is always the one with larger Pt.
  double dR1_1=999, dR2_1=999;
  double dR1_2=999, dR2_2=999;
  TLorentzVector v1reco, v2reco, v1gen, v2gen;
  v1reco.SetPtEtaPhiM(dielectron->pt_1, dielectron->eta_1, dielectron->phi_1, 0.000511);
  v2reco.SetPtEtaPhiM(dielectron->pt_2, dielectron->eta_2, dielectron->phi_2, 0.000511);
  v1gen .SetPtEtaPhiM(gen->pt_1, gen->eta_1, gen->phi_1, 0.000511);
  v2gen .SetPtEtaPhiM(gen->pt_2, gen->eta_2, gen->phi_2, 0.000511);
  // Try both assignments, at least one assignment should match
  dR1_1 = v1reco.DeltaR(v1gen);
  dR2_1 = v2reco.DeltaR(v2gen);
  dR1_2 = v1reco.DeltaR(v2gen);
  dR2_2 = v2reco.DeltaR(v1gen);
  // Require that both are within the required cone
  bool matchAssignment1 = (fabs(dR1_1) < 0.2 && fabs(dR2_1) < 0.2 );
  bool matchAssignment2 = (fabs(dR1_2) < 0.2 && fabs(dR2_2) < 0.2 );
  if( ! (matchAssignment1 || matchAssignment2) ) result = kFALSE; 
  
  return result;
}

Bool_t electronMatchedToGeneratorLevel(const mithep::TGenInfo *gen, const mithep::TElectron *electron){
  
  Bool_t result = kTRUE;
  
  // Matching is done as prescribed in AN-12-116
  //   - dR < 0.2
  //   - match to status 1 (post-FSR)
  //   - no charge matching

  // In the generator branch of this ntuple, first particle is always
  // negative, and second always positive (but this is not used at present
  // as there is no charge matching).
  double dR1=999;
  double dR2=999;
  TLorentzVector vreco, v1gen, v2gen;
  vreco.SetPtEtaPhiM(electron->pt, electron->eta, electron->phi, 0.000511);
  v1gen .SetPtEtaPhiM(gen->pt_1, gen->eta_1, gen->phi_1, 0.000511);
  v2gen .SetPtEtaPhiM(gen->pt_2, gen->eta_2, gen->phi_2, 0.000511);
  
  dR1 = vreco.DeltaR(v1gen);
  dR2 = vreco.DeltaR(v2gen);

  if( !( fabs(dR1) < 0.2 || fabs(dR2) < 0.2 ) ) result = kFALSE; 
  
  return result;
}

Bool_t scMatchedToGeneratorLevel(const mithep::TGenInfo *gen, const mithep::TPhoton *sc){

  Bool_t result = kTRUE;
  // We do not know which of the gen electrons possibly
  // produced this supercluster, so we check both.
  double dR1=999, dR2=999;
  TLorentzVector vreco, v1gen, v2gen;
  vreco.SetPtEtaPhiM(sc->pt, sc->eta, sc->phi, 0.000511);
  v1gen .SetPtEtaPhiM(gen->pt_1, gen->eta_1, gen->phi_1, 0.000511);
  v2gen .SetPtEtaPhiM(gen->pt_2, gen->eta_2, gen->phi_2, 0.000511);
  dR1 = vreco.DeltaR(v1gen);
  dR2 = vreco.DeltaR(v2gen);
  // Require that at least one is  within loose dR of 0.4, otherwise bail out
  if( fabs(dR1) > 0.4 && fabs(dR2) > 0.4 ) result = kFALSE; 
  
  return result;
}

bool passID(const mithep::TElectron *electron, double rho){

  bool result = passEGM2011(electron, WP_TIGHT, rho);
  return result;
}

bool passIDTag(const mithep::TElectron *electron, double rho){

  bool result = passEGM2011(electron, WP_TIGHT, rho);
  return result;
}

bool isTag(const mithep::TElectron *electron, ULong_t trigger, double rho){

  bool elePassID  = passIDTag(electron, rho);
  bool elePassHLT =  (electron ->hltMatchBits & trigger);

  bool result = ( elePassID && elePassHLT && (electron->pt > 25) );

  return result;
}

// -------------------------------------------------------------------

TString getLabel(int sample, int effType, int method,  int etBinning, int etaBinning, const TriggerSelection &trigSet){
  using namespace DYTools;

  TString label = analysisTag;
  if (analysisTag.Length()>0) label.Append("_");

  assert ( trigSet.isDefined() );
  if (sample != -1111) {
    label+= trigSet.triggerSetName();
  }

  if(sample == DATA)
    label += "_data";
  else if((sample == MC) || (sample == -1111))
    label += "_mc";
  else
    assert(0);

  if( effType == RECO )
    label += "_reco";
  else if( effType == ID )
    label += "_id";
  else if(effType == HLT ) {
    if ((sample==DATA) && trigSet.hltEffMethodIs2011New()) label += "_hlt2011new";
    else label += "_hlt";
  }
  else
    assert(0);

  if (sample != -1111) {
    if(method == COUNTnCOUNT)
      label += "_count-count";
    else if( method == COUNTnFIT ) 
      label += "_count-fit";
    else if( method == FITnFIT ) 
      label += "_fit-fit";
    else
      assert(0);
  }

  label += "_bins-et";
  label += getNEtBins(etBinning);
  label += "-eta";
  label += getNEtaBins(etaBinning);

  return label;
}

// -------------------------------------------------------------------
