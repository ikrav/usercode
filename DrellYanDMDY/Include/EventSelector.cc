#include "EventSelector.hh"
#include <TLorentzVector.h>

// ---------------------------------------------------------------

DielectronSelector_t::DielectronSelector_t(TSelectionType_t set_selection, ElectronEnergyScale *escale) :
  fSelection(set_selection),
  fEScale(escale),
  fTotalCandidates(0), fCandidatesGoodEta(0), fCandidatesGoodEt(0),
  fCandidatesHLTMatched(0), fCandidatesIDPassed(0), 
  fCandidatesMassAboveMinLimit(0)
{
}

// ---------------------------------------------------------------

bool DielectronSelector_t::testDielectron_default(mithep::TDielectron *dielectron, TEScaleCorrection_t applyEScale, ULong_t leadingTriggerObjectBit, ULong_t trailingTriggerObjectBit) {

  fTotalCandidates++;

  bool ok=kFALSE;
  // assume 1 iteration loop
  for (int dummy_iter=0; dummy_iter<1; ++dummy_iter) {
    // Exclude ECAL gap region and cut out of acceptance electrons
    if((fabs(dielectron->scEta_1)>kECAL_GAP_LOW) && (fabs(dielectron->scEta_1)<kECAL_GAP_HIGH)) continue;
    if((fabs(dielectron->scEta_2)>kECAL_GAP_LOW) && (fabs(dielectron->scEta_2)<kECAL_GAP_HIGH)) continue;
    if((fabs(dielectron->scEta_1) > 2.5)       || (fabs(dielectron->scEta_2) > 2.5))       continue;  // outside eta range? Skip to next event...
  
    fCandidatesGoodEta++;

    //
    // Energy scale corrections for data
    // NOTE: the electrons and dielectron 4-vectors are updated, the supercluster quantities are not
    //
    Double_t scEt1 = dielectron->scEt_1;
    Double_t scEt2 = dielectron->scEt_2;
    // Electron energy scale correction
    if(applyEScale==_escaleData) {
      if (!fEScale) {
	std::cout << "Error in testDielectron_default: escaleData is requested, but the pointer is null\n";
	throw 2;
      }

      double corr1 = fEScale->getEnergyScaleCorrection(dielectron->scEta_1);
      double corr2 = fEScale->getEnergyScaleCorrection(dielectron->scEta_2);
      scEt1 = dielectron->scEt_1 * corr1;
      scEt2 = dielectron->scEt_2 * corr2;

      TLorentzVector ele1; 
      ele1.SetPtEtaPhiM(dielectron->pt_1,dielectron->eta_1,dielectron->phi_1,0.000511);
      ele1 *= corr1;
      dielectron->pt_1  = ele1.Pt();
      dielectron->eta_1 = ele1.Eta();
      dielectron->phi_1 = ele1.Phi();
      
      TLorentzVector ele2; 
      ele2.SetPtEtaPhiM(dielectron->pt_2,dielectron->eta_2,dielectron->phi_2,0.000511);
      ele2 *= corr2;
      dielectron->pt_2  = ele2.Pt();
      dielectron->eta_2 = ele2.Eta();
      dielectron->phi_2 = ele2.Phi();
      
      TLorentzVector vDiEle = ele1+ele2;            
      dielectron->mass = vDiEle.M();
      dielectron->pt   = vDiEle.Pt();
      dielectron->y    = vDiEle.Rapidity();
      dielectron->phi  = vDiEle.Phi(); 
    }
       	  
    // requirements on BOTH electrons
    // For DY ET cuts are asymmetric:
    if( ! ( (scEt1>20 && scEt2>10) || (scEt1>10 && scEt2>20) ) ) continue;
    
    fCandidatesGoodEt++;

    // Both electrons must match trigger objects. At least one ordering
    // must match
    if( ! ( 
	   (dielectron->hltMatchBits_1 & leadingTriggerObjectBit && 
	    dielectron->hltMatchBits_2 & trailingTriggerObjectBit )
	   ||
	   (dielectron->hltMatchBits_1 & trailingTriggerObjectBit && 
	    dielectron->hltMatchBits_2 & leadingTriggerObjectBit ) ) ) continue;
    
    fCandidatesHLTMatched++;

    // Other cuts to both electrons
    
    // The Smurf electron ID package is the same as used in HWW analysis
    // and contains cuts like VBTF WP80 for pt>20, VBTF WP70 for pt<10
    // with some customization, plus impact parameter cuts dz and dxy
    if(!passSmurf(dielectron)) continue;

    fCandidatesIDPassed++;

    // loose mass window 
    double minMass= (applyEScale==_escaleUncorrected) ? 5 : 10;
      if( dielectron->mass < minMass ) continue;

    fCandidatesMassAboveMinLimit++;

    ok=kTRUE; // selection PASSED

  } while(0); // fictitious loop
  return ok;
}


// ---------------------------------------------------------------

std::ostream& DielectronSelector_t::printCounts(std::ostream &out) {
  char buf[20];
  out << "DielectronSelector (selection=" << DielectronSelector_t::selectionName(fSelection) << ")\n";
  const char *format="%10u";
  sprintf(buf,format,fTotalCandidates);
  out << "Total number of dielectron candidates  " << buf << "\n";
  sprintf(buf,format,fCandidatesGoodEta);
  out << "Total candidates with good eta         " << buf << "\n";
  sprintf(buf,format,fCandidatesGoodEt);
  out << "Total candidates with good Et          " << buf << "\n";
  sprintf(buf,format,fCandidatesHLTMatched);
  out << "Total candidates HLT matched           " << buf << "\n";
  sprintf(buf,format,fCandidatesIDPassed);
  out << "Total candidates ID passed             " << buf << "\n";
  sprintf(buf,format,fCandidatesMassAboveMinLimit);
  out << "Total candidates with mass above limit " << buf << "\n";
  return out;
}

// ---------------------------------------------------------------

// ---------------------------------------------------------------
