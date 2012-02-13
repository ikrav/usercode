#ifndef ElectronEnergyScale_HH
#define ElectronEnergyScale_HH

#include <TString.h>
#include <TF1.h>
#include <TRandom.h>
#include <TH1F.h>

class ElectronEnergyScale {

public:

  enum CalibrationSet {
    UNDEFINED=0,
    UNCORRECTED,
    Date20110901_EPS11_default,
    Date20120101_default,
    CalSet_File_Gauss
  };

  // Constructor
  ElectronEnergyScale(CalibrationSet calibrationSet);
  ElectronEnergyScale(const TString &escaleTagName);
  // Destructor
  //~ElectronEnergyScale() { this->clear(); }
  // Clean-up
  void clear();

  // Initialization
  void init(CalibrationSet calSet);
  void init(const TString &stringWithEScaleTagName);
  bool initializeAllConstants();
  bool initializeExtraSmearingFunction();
  bool isInitialized() const { return _isInitialized; }
  bool isScaleRandomized() const { return _energyScaleCorrectionRandomizationDone; }
  bool isSmearRandomized() const { return _smearingWidthRandomizationDone; }

  bool loadInputFile(const TString &fileName, int debug=0);
  bool assignConstants(const std::vector<string> &lines, int debug=0);
  static bool AssignConstants(const std::vector<string> &lines, int count, double *eta, double *scale, double *scaleErr, double *smear, double *smearErr, int debug=0);

  // Access
  CalibrationSet getCalibrationSet() const { return _calibrationSet; }
  bool setCalibrationSet(CalibrationSet calSet); // in specific cases allow to change the calibration set
  int numberOfEtaBins() const { return _nEtaBins; }
  int numberOfEtaEtaBins() const { return _nEtaBins*(_nEtaBins+1); }

  // Eta bin index in the same way as TH1F bin index
  int getEtaBinIdx(double eta) const {
    int idx=-1;
    for (int i=0; (idx<0) && (i<_nEtaBins); ++i) {
      if ((eta>_etaBinLimits[i]) && (eta<=_etaBinLimits[i+1])) idx=i+1;
    }
    return idx;
  }

  double getEtaBinCenter(int i) const { return 0.5*(_etaBinLimits[i-1]+_etaBinLimits[i]); }

  double getEnergyScaleCorrection(double eta) const;
  // old-style smear (event shift)
  double generateMCSmear(double eta1, double eta2) const;
  // updated smear (distribution) : one event
  bool addSmearedWeight(TH1F *hMassDestination, int eta1Bin, int eta2Bin, double mass, double weight) const {
    assert(this->isInitialized());
    return addSmearedWeightAny(hMassDestination,eta1Bin,eta2Bin,mass,weight,kFALSE);
  }
  // updated smear (distribution) : smear collection
  void smearDistribution(TH1F *destination, int eta1Bin, int eta2Bin, const TH1F *source) const {
    assert(this->isInitialized());
    smearDistributionAny(destination,eta1Bin,eta2Bin,source,kFALSE);
  }

  // Several functions for systematic studies. These
  // randomize constants for energy scale corrections
  // within their errors.
  void   randomizeEnergyScaleCorrections(int seed);
  double getEnergyScaleCorrectionRandomized(double eta) const;

  void   randomizeSmearingWidth(int seed);

  // old-style smear (event shift)
  double generateMCSmearRandomized(double eta1, double eta2) const;
  // updated smear (distribution) : one event
  bool addSmearedWeightRandomized(TH1F *hMassDestination, int eta1Bin, int eta2Bin, double mass, double weight) const {
    assert(this->isInitialized()); assert(this->isSmearRandomized());
    return addSmearedWeightAny(hMassDestination,eta1Bin,eta2Bin,mass,weight,kTRUE);
  }
  // updated smear (distribution) : smear collection
  void smearDistributionRandomized(TH1F *destination, int eta1Bin, int eta2Bin, const TH1F *source) const {
    assert(this->isInitialized()); assert(this->isSmearRandomized());
    smearDistributionAny(destination,eta1Bin,eta2Bin,source,kTRUE);
  }

  void print() const;

  // Useful functions
  static CalibrationSet DetermineCalibrationSet(const TString &escaleTagName, TString *inputFileName=NULL);  // TString -> CalibrationSet
  static TString CalibrationSetName(CalibrationSet escaleTag, const TString *fileName);   // CalibrationSet -> TString
  static TString CalibrationSetFunctionName(CalibrationSet escaleTag); // name of the calibrating function

  TString calibrationSetName() const { return ElectronEnergyScale::CalibrationSetName(this->_calibrationSet, &this->_inpFileName); }
  TString calibrationSetFunctionName() const { return ElectronEnergyScale::CalibrationSetFunctionName(this->_calibrationSet); }

  TH1F* createScaleHisto(const TString &namebase) const;
  TH1F* createSmearHisto(const TString &namebase, int parameterNo) const;


protected:
  // Internal functions, not for general use
  double getEnergyScaleCorrectionAny(double eta, bool randomize) const;

  // old-style smear (event shift)
  double generateMCSmearAny(double eta1, double eta2, bool randomize) const;
  // updated smear (distribution) : one event
  bool addSmearedWeightAny(TH1F *hMassDestination, int eta1Bin, int eta2Bin, double mass, double weight, bool randomize) const;
  // updated smear (distribution) : smear collection
  void smearDistributionAny(TH1F *destination, int eta1Bin, int eta2Bin, const TH1F *source, bool randomize) const;

  TH1F* createParamHisto(const TString &namebase, const TString &nameTag, const double *params, const double *paramErrs) const;


private:

  CalibrationSet     _calibrationSet;
  TString            _inpFileName;
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

protected:
  // Functions to be used for extra smearing
  static const int nMaxFunctions = 50;
  TF1 *smearingFunctionGrid[nMaxFunctions][nMaxFunctions];
  TF1 *smearingFunctionGridRandomized[nMaxFunctions][nMaxFunctions];

};


#endif
