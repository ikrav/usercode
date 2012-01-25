#ifndef ElectronEnergyScaleAdv_HH
#define ElectronEnergyScaleAdv_HH


#include "../../Include/ElectronEnergyScale.hh"
#include "../../Include/EtaEtaMass.hh"

class ElectronEnergyScaleAdv_t : public ElectronEnergyScale {

public:

  // Constructor
  ElectronEnergyScaleAdv_t(CalibrationSet calibrationSet) : ElectronEnergyScale(calibrationSet) {}
  ElectronEnergyScaleAdv_t(const TString &escaleTagName) : ElectronEnergyScale(escaleTagName) {}

  int getEtaEtaIdx(double eta1, double eta2) const;

  int LoadEEMFile(const TString &eemFileName, vector<vector<double>*> &eemData) const;

};


#endif
