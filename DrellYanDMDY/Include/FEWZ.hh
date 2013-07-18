#ifndef FEWZ_HH
#define FEWZ_HH

#include "DYTools.hh"
#include <TH2D.h>

//  Definition of a FEWZ_t class which knows the internal structure
//  of FEWZ correction file and provides the correction factor for 
// the provided values of mass and rapidity y

// 1) weights used for 7 TeV analysis are provided as histograms.
// 2) weights used for 8 TeV analysis are provided as a macro from Alexey. 
// This class interfaces the macro. The definitions need adaptation in case
// the macro is changed

#ifdef _8TeV_analysis_
// 8 TeV FEWZ map is based on file from Alexey
const int _nMassBinsFEWZ_8TeV = 14;
const double _massBinLimitsFEWZ_8TeV[_nMassBinsFEWZ_8TeV + 1] = 
  { 15, 20, 30, 45, 60, 72, 106, 120, 133, 150, 160, 171, 200, 400, 1500 };
const int _nMassBinsFEWZ = _nMassBinsFEWZ_8TeV;
const double *_massBinLimitsFEWZ = _massBinLimitsFEWZ_8TeV;
#else
// 7 TeV FEWZ map is based on file with histograms
const int _nMassBinsFEWZ = 40;
const double _massBinLimitsFEWZ[_nMassBinsFEWZ+1] = 
  {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 
   81, 86, 91, 96, 101, 106, 110, 115, 120, 126, 133, 141, 
   150, 160, 171, 185, 200, 220, 243, 273, 320, 380, 440, 
   510, 600, 1000, 1500}; // 40 bins
#endif


#ifdef _8TeV_analysis_
const int nptbins_FEWZ8TeV = 21;
const int nrapbins_FEWZ8TeV = 8;
const double pt_bin_FEWZ8TeV[nptbins_FEWZ8TeV] = {0.0,20.0,30.,35.,40.,45.,50.,55.,60.,65.,70.,75.,80.,85.,90.,100.,120.,150.,300.,600.,1000.};
const double rap_bin_FEWZ8TeV[nrapbins_FEWZ8TeV] = {0.,0.2,0.4,0.7,1.1,1.9,2.4,1000.0};
#endif


// --------------------------------------

#ifdef _8TeV_analysis_
//int Find_Index_FEWZ8TeV( double muon_input, bool isPt ); // A.Svyatkovskiy function
double weight_FEWZ8TeV( double muon_pt, double muon_rap, double mass); // A.Svyatkovskiy function
#endif


// --------------------------------------

class FEWZ_t {
protected:
  TH2D *weights[_nMassBinsFEWZ], *weightErrors[_nMassBinsFEWZ];
  bool fInitialized, fCutZPT100;
public:
  FEWZ_t(bool loadWeights=true, bool do_cutZPT100=true);
  bool isInitialized() const { return fInitialized; }
  bool cutZPT100() const { return fCutZPT100; }

  // ------
  // this method is used for 7 TeV analysis
  // ------
#ifndef _8TeV_analysis_
  Double_t getValue_7TeV(TH2D **values, Double_t gen_ee_mass, Double_t gen_pt, Double_t gen_rapidity) const {
    if (!fInitialized) {
      std::cout << "Error: FEWZ::getValue_7TeV is requested on uninitialized object\n";
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
#endif

  // ------
  // mass bin index
  // ------

  Int_t getIndex(Double_t gen_ee_mass) const {
    Int_t ibinPreFsr = DYTools::_findMassBin(gen_ee_mass,_nMassBinsFEWZ,_massBinLimitsFEWZ);
    return ibinPreFsr;
  }

  // ------
  // helping methods
  // ------

  double getMassCenter(Int_t mass_bin) const {
    double mass=0.5*(_massBinLimitsFEWZ[mass_bin] + _massBinLimitsFEWZ[mass_bin+1]);
    return mass;
  }

  double getLowMassEdge(Int_t mass_bin) const {
    double mass=_massBinLimitsFEWZ[mass_bin];
    return mass;
  }

  double getHighMassEdge(Int_t mass_bin) const {
    double mass=_massBinLimitsFEWZ[mass_bin+1];
    return mass;
  }

#ifdef _8TeV_analysis_
  int prepareHisto_8TeV(Int_t mass_bin, int abs_rapidity=0);
#endif

  // ------
  // access to histograms
  // in 7TeV analysis, the pointer is returned
  // in 8TeV analysis, the histogram is also created!
  // ------
  
  const TH2D *getHistoPtr(Int_t mass_bin, int abs_rapidity=0) {
#ifdef _8TeV_analysis_
    if (!weights[mass_bin]) {
      // the histogram has to be created
      prepareHisto_8TeV(mass_bin,abs_rapidity);
    }
#endif
    return weights[mass_bin];
  }


  // ------
  // main method
  // ------

  Double_t getWeight(Double_t gen_ee_mass, Double_t gen_pt, Double_t gen_rapidity) {
#ifdef _8TeV_analysis_
    return weight_FEWZ8TeV(gen_pt,fabs(gen_rapidity),gen_ee_mass);
#else
    return getValue(weights,gen_ee_mass,gen_pt,gen_rapidity);
#endif
  }

  // ------
  // error is available only for 7 TeV analysis
  // ------


  Double_t getWeightError(Double_t gen_ee_mass, Double_t gen_pt, Double_t gen_rapidity) {
#ifdef _8TeV_analysis_
    if (0) std::cout << "getWeightError(mass=" << gen_ee_mass << ", pt=" << gen_pt << ", gen_rapidity=" << gen_rapidity << ")\n"; // get rid of a compiler warning
    return 0.;
#else
    return getValue(weightErrors,gen_ee_mass,gen_pt,gen_rapidity);
#endif
  }

};

// --------------------------------------


#endif
