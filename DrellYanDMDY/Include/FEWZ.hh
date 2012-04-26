#ifndef FEWZ_HH
#define FEWZ_HH

#include "DYTools.hh"
#include <TH2D.h>

//  Definition of a FEWZ_t class which knows the internal structure
//  of FEWZ correction file and provides the correction factor for 
// the provided values of mass and rapidity y

const int _nMassBinsFEWZ = 40;
const double _massBinLimitsFEWZ[_nMassBinsFEWZ+1] = 
  {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 
   81, 86, 91, 96, 101, 106, 110, 115, 120, 126, 133, 141, 
   150, 160, 171, 185, 200, 220, 243, 273, 320, 380, 440, 
   510, 600, 1000, 1500}; // 40 bins

// --------------------------------------

class FEWZ_t {
protected:
  TH2D *weights[_nMassBinsFEWZ], *weightErrors[_nMassBinsFEWZ];
  bool fInitialized, fCutZPT100;
public:
  FEWZ_t(bool loadWeights, bool do_cutZPT100);
  bool isInitialized() const { return fInitialized; }
  bool cutZPT100() const { return fCutZPT100; }

  Double_t getValue(TH2D **values, Double_t gen_ee_mass, Double_t gen_pt, Double_t gen_rapidity) const {
    if (!fInitialized) {
      std::cout << "Error: FEWZ::getValue is requested on uninitialized object\n";
      return Double_t(0);
    }
    Double_t value=1.0;
    int ibinPreFsr = DYTools::_findMassBin(gen_ee_mass,_nMassBinsFEWZ,_massBinLimitsFEWZ);
    // If mass is larger than the highest bin boundary
    // (last bin), use the last bin.
    if(ibinPreFsr == -1 && gen_ee_mass >= _massBinLimitsFEWZ[_nMassBinsFEWZ] )
      ibinPreFsr = _nMassBinsFEWZ-1;
    // Find FEWZ-powheg reweighting factor 
    // that depends on pre-FSR Z/gamma* rapidity, pt, and mass
    if(ibinPreFsr != -1 && ibinPreFsr < _nMassBinsFEWZ){
      int ptBin = values[ibinPreFsr]->GetXaxis()->FindBin( gen_pt );
      int yBin = values[ibinPreFsr]->GetYaxis()->FindBin( gen_rapidity );
      // In case if pt or y are outside of the weight maps,
      // set them to the closest bin.
      if(ptBin == values[ibinPreFsr]->GetNbinsX() + 1)
	ptBin = values[ibinPreFsr]->GetNbinsX();
      if(ptBin == 0)
	ptBin = 1;
      if(yBin == values[ibinPreFsr]->GetNbinsY() + 1)
	yBin = values[ibinPreFsr]->GetNbinsY();
      if(yBin == 0)
	yBin = 1;
      // Apply PT cut if needed
      if( fCutZPT100 ) 
	if( ptBin == values[ibinPreFsr]->GetNbinsX() )
	  ptBin = values[ibinPreFsr]->GetNbinsX() - 1;

      value = values[ibinPreFsr]->GetBinContent( ptBin, yBin);
    }
    else {
      // Error printout is commented out: the maps now go down
      // to 15 GeV. Events with generator level mass below 15 GeV 
      // may contribute to 
      // reconstructed events when reconstructed mass is above 15 GeV
// 	      cout << "Error: vmass outside of FEWZ weight maps, vmass=" 
// 		   << gen->vmass << endl;
    }
    return value;
  }

  Double_t getWeight(Double_t gen_ee_mass, Double_t gen_pt, Double_t gen_rapidity) {
    return getValue(weights,gen_ee_mass,gen_pt,gen_rapidity);
  }

  Double_t getWeightError(Double_t gen_ee_mass, Double_t gen_pt, Double_t gen_rapidity) {
    return getValue(weightErrors,gen_ee_mass,gen_pt,gen_rapidity);
  }

};

// --------------------------------------


#endif
