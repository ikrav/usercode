#ifndef ElectronEnergyScale_HH
#define ElectronEnergyScale_HH

#include <TString.h>
#include <TF1.h>
#include <TRandom.h>

class ElectronEnergyScale {

public:

  enum CalibrationSet {
    UNDEFINED,
    Date20110901_EPS11_default
  };

  // Constructor
  ElectronEnergyScale(CalibrationSet calibrationSet);

  // Initialization
  bool initializeAllConstants();
  bool initializeExtraSmearingFunction();

  // Access
  double getEnergyScaleCorrection(double eta);
  double generateMCSmear(double eta1, double eta2);

  // Several functions for systematic studies. These
  // randomize constants for energy scale corrections
  // within their errors.
  void   randomizeEnergyScaleCorrections(int seed);
  double getEnergyScaleCorrectionRandomized(double eta);

  void   randomizeSmearingWidth(int seed);
  double generateMCSmearRandomized(double eta1, double eta2);

  void print();

  // Internal functions, not for general use
  double getEnergyScaleCorrectionAny(double eta, bool randomize);
  double generateMCSmearAny(double eta1, double eta2, bool randomize);

private:

  CalibrationSet     _calibrationSet;
  bool               _isInitialized;

  // The data structure assumes the following:
  //  - the binning is in eta only
  //  - data constants contain a single scale factor 
  //  - MC constants contain several variable describing the shape
  //         of the smearing function
  //  - each constant, data or MC, comes with an error

  // Data constants
  int                    _nEtaBins;
  double *               _etaBinLimits;
  double *               _dataConst;
  double *               _dataConstErr;
  // randomized for systematics studies
  double *               _dataConstRandomized;
  bool                   _energyScaleCorrectionRandomizationDone;
  
  // MC constants. The number of constants and the names depend
  // on particular calibration set. For now, maximum possible is four.
  int                    _nMCConstants;
  TString                _mcConst1Name;
  TString                _mcConst2Name;
  TString                _mcConst3Name;
  TString                _mcConst4Name;
  double *               _mcConst1;
  double *               _mcConst2;
  double *               _mcConst3;
  double *               _mcConst4;
  double *               _mcConst1Err;
  double *               _mcConst2Err;
  double *               _mcConst3Err;
  double *               _mcConst4Err;
  bool                   _smearingWidthRandomizationDone;

  // Functions to be used for extra smearing
  static const int nMaxFunctions = 50;
  TF1 *smearingFunctionGrid[nMaxFunctions][nMaxFunctions];
  TF1 *smearingFunctionGridRandomized[nMaxFunctions][nMaxFunctions];

};

#endif
