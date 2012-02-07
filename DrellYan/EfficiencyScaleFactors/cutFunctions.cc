#if !defined(__CINT__) || defined(__MAKECINT__)
#include "../Include/cutFunctions.hh"
#endif

Bool_t dielectronMatchedToGeneratorLevel(const mithep::TGenInfo *gen, const mithep::TDielectron *dielectron){

  Bool_t result = kTRUE;
  // In the generator branch of this ntuple, first particle is always
  // negative, and second always positive. In the Dielectron block
  // of the ntuple, the first particle is always the one with larger Pt.
  double dR1=999, dR2=999;
  TLorentzVector v1reco, v2reco, v1gen, v2gen;
  v1reco.SetPtEtaPhiM(dielectron->pt_1, dielectron->eta_1, dielectron->phi_1, 0.000511);
  v2reco.SetPtEtaPhiM(dielectron->pt_2, dielectron->eta_2, dielectron->phi_2, 0.000511);
  v1gen .SetPtEtaPhiM(gen->pt_1, gen->eta_1, gen->phi_1, 0.000511);
  v2gen .SetPtEtaPhiM(gen->pt_2, gen->eta_2, gen->phi_2, 0.000511);
  if( dielectron->q_1 < 0 ){
    dR1 = v1reco.DeltaR(v1gen);
    dR2 = v2reco.DeltaR(v2gen);
  }else{
    dR1 = v1reco.DeltaR(v2gen);
    dR2 = v2reco.DeltaR(v1gen);
  }
  // Require that both are within loose dR of 0.4, otherwise bail out
  if( fabs(dR1) > 0.4 || fabs(dR2) > 0.4 ) result = kFALSE; 
  
  return result;
}

Bool_t electronMatchedToGeneratorLevel(const mithep::TGenInfo *gen, const mithep::TElectron *electron){
  
  Bool_t result = kTRUE;
  
  // In the generator branch of this ntuple, first particle is always
  // negative, and second always positive. 
  double dR=999;
  TLorentzVector vreco, v1gen, v2gen;
  vreco.SetPtEtaPhiM(electron->pt, electron->eta, electron->phi, 0.000511);
  v1gen .SetPtEtaPhiM(gen->pt_1, gen->eta_1, gen->phi_1, 0.000511);
  v2gen .SetPtEtaPhiM(gen->pt_2, gen->eta_2, gen->phi_2, 0.000511);
  
  if( electron->q < 0 ){
    dR = vreco.DeltaR(v1gen);
  }else{
    dR = vreco.DeltaR(v2gen);
  }
  // Require that both are within loose dR of 0.4, otherwise bail out
  if( fabs(dR) > 0.4 ) result = kFALSE; 
  
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

bool passID(const mithep::TElectron *electron){

  bool result = passSmurf(electron);
  return result;
}

bool isTag(const mithep::TElectron *electron, ULong_t trigger){

  bool elePassID  = passID(electron);
  bool elePassHLT =  (electron ->hltMatchBits & trigger);

  bool result = ( elePassID && elePassHLT && (electron->scEt > 20) );

  return result;
}

// -------------------------------------------------------------------

TString getLabel(int sample, int effType, int method,  int etBinning, int etaBinning, const TriggerSelection &trigSet){

  TString label = "";

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

  if( effType == GSF )
    label += "_gsf";
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
