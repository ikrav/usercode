//
// This file contains implementation of methods related to applying
// electron energy scale corrections.
//

#include "../Include/ElectronEnergyScale.hh"

ElectronEnergyScale::ElectronEnergyScale(CalibrationSet calibrationSet):
  _calibrationSet(calibrationSet),
  _isInitialized(false),
  _energyScaleCorrectionRandomizationDone(false),
  _smearingWidthRandomizationDone(false)
{

  if( calibrationSet == UNDEFINED )
    return;

  if( !initializeAllConstants())
    return;

  if( !initializeExtraSmearingFunction())
    return;

  _isInitialized = true;
  return;
}

bool ElectronEnergyScale::initializeAllConstants(){
  
  bool success = true;
  if( _calibrationSet == Date20110901_EPS11_default ){
    //
    // Constants from energy scale calibrations
    // done for Summer 11 result. Note that the 
    // constants are symmetric about eta = 0.
    //
    const int nEtaBins = 12;
    const double etaBinLimits[nEtaBins+1] = 
      {-2.50001, -2.0, -1.5, -1.2, -0.8, -0.4, 0.0, 0.4, 0.8, 1.2, 1.5, 2.0, 2.50001};
    
    const double corrValues[nEtaBins] = 
      {1.04642, 1.00187, 1.01556, 1.00500, 1.00093, 1.00149, 1.00149, 1.00093, 1.00500, 1.01556, 1.00187, 1.04642};
    const double corrErrors[nEtaBins] = 
      {4.28928e-04,3.39718e-04,4.89342e-04,5.80480e-05,1.21192e-05,1.27489e-04,1.27489e-04,1.21192e-05,5.80480e-05,4.89342e-04,3.39718e-04,4.28928e-04};
    
    const double smearValues[nEtaBins] = 
      {2.05888e+00,1.46747e+00,1.14861e+00,7.63770e-01,5.04140e-01,5.27258e-01,5.27258e-01,5.04140e-01,7.63770e-01,1.14861e+00,1.46747e+00,2.05888e+00};
    const double smearErrors[nEtaBins] =
      {2.85889e-02,3.85260e-02,4.26451e-02,3.22979e-02,3.76972e-02,3.32377e-02,3.32377e-02,3.76972e-02,3.22979e-02,4.26451e-02,3.85260e-02,2.85889e-02};
    //
    // Prepare those constants that are used in this parameterization for
    // assigning the values.
    // 
    double dummyArray[nEtaBins];
    double dummyArray2[nEtaBins+1];
    _nEtaBins = nEtaBins;
    _etaBinLimits = (double *)malloc(sizeof(dummyArray2));
    _dataConst    = (double *)malloc(sizeof(dummyArray));
    _dataConstErr = (double *)malloc(sizeof(dummyArray));
    _dataConstRandomized = (double *)malloc(sizeof(dummyArray));
    _nMCConstants = 1;
    _mcConst1Name = "smear";
    _mcConst2Name = "<none>";
    _mcConst3Name = "<none>";
    _mcConst4Name = "<none>";
    _mcConst1 = (double *)malloc(sizeof(dummyArray));
    _mcConst2 = 0;
    _mcConst3 = 0;
    _mcConst4 = 0;
    _mcConst1Err = (double *)malloc(sizeof(dummyArray));
    _mcConst2Err = 0;
    _mcConst3Err = 0;
    _mcConst4Err = 0;
    //
    for(int i=0; i<nEtaBins; i++){
      _etaBinLimits[i] = etaBinLimits[i];
      _dataConst   [i] = corrValues[i];
      _dataConstErr[i] = corrErrors[i];
      _mcConst1    [i] = smearValues[i];
      _mcConst1Err [i] = smearErrors[i];
    }
    _etaBinLimits[nEtaBins] = etaBinLimits[nEtaBins];
  }
  
  return success;
}

// The extra smearing function is to provide smearing for
// the mass of an event based on the individual parameters
// of two electrons. Thus the function is actually an 2D
// array of functions, with each pair of (i,j) eta bins for
// a given dielectron candidate corresponding to its unique
// smearing function.
bool ElectronEnergyScale::initializeExtraSmearingFunction(){

  bool success = true;
  // A sanity check. The function that initializes constants
  // should have been run by now.
  if( _nEtaBins <= 0 && _nEtaBins > nMaxFunctions)
    return false;

  for( int i=0; i<_nEtaBins; i++){
    for( int j=0; j<_nEtaBins; j++){
      TString fname = TString::Format("smearing_function_%03d_%03d", i, j);
      if( _calibrationSet == Date20110901_EPS11_default ){
	if(_mcConst1 == 0) continue;
	smearingFunctionGrid[i][j] = new TF1(fname, "gaus(0)", -10, 10);
	smearingFunctionGrid[i][j]->SetNpx(500);
	double si = _mcConst1[i];
	double sj = _mcConst1[j];
	smearingFunctionGrid[i][j]->SetParameters(1.0,0.0,sqrt(si*si+sj*sj));
      } else {
	success = false;
      }
    } // end inner loop ove eta bins
  } // end outer loop over eta bins

  return success;
}

void   ElectronEnergyScale::randomizeEnergyScaleCorrections(int seed){

  if( !_isInitialized ){
    printf("ElectronEnergyScale ERROR: the object is not properly initialized\n");
    return;
  }

  TRandom rand;
  rand.SetSeed(seed);
  _energyScaleCorrectionRandomizationDone = true;
  for(int i=0; i<_nEtaBins; i++){
    _dataConstRandomized[i] = _dataConst[i] + rand.Gaus(0.0,_dataConstErr[i]);
  }

  return;
}

double ElectronEnergyScale::getEnergyScaleCorrection(double eta){

  double result = 1.0;
  bool randomize = false;
  result = getEnergyScaleCorrectionAny(eta,randomize);

  return result;
}

double ElectronEnergyScale::getEnergyScaleCorrectionRandomized(double eta){

  double result = 1.0;
  if( !_energyScaleCorrectionRandomizationDone ){
    printf("ElectronEnergyScale ERROR: can not get randomized escale, randomization is not done\n");
    return result;
  }

  bool randomize = true;
  result = getEnergyScaleCorrectionAny(eta,randomize);

  return result;
}

double ElectronEnergyScale::getEnergyScaleCorrectionAny(double eta, bool randomize){

  double result = 1.0;
  if( !_isInitialized ){
    printf("ElectronEnergyScale ERROR: the object is not properly initialized\n");
    return result;
  }

  for(int i=0; i<_nEtaBins; i++){
    if(eta >= _etaBinLimits[i] && eta < _etaBinLimits[i+1] ){
      if( !randomize )
	result = _dataConst[i];
      else
	result = _dataConstRandomized[i];
      break;
    }
  }

  return result;
}

void ElectronEnergyScale::randomizeSmearingWidth(int seed){

  if( !_isInitialized ){
    printf("ElectronEnergyScale ERROR: the object is not properly initialized\n");
    return;
  }

  TRandom rand;
  rand.SetSeed(seed);
  _smearingWidthRandomizationDone = true;
  if( _calibrationSet == Date20110901_EPS11_default ){

    for( int i=0; i<_nEtaBins; i++){
      for( int j=0; j<_nEtaBins; j++){
	TString fname = TString::Format("smearing_function_randomized_%03d_%03d", i, j);
	smearingFunctionGridRandomized[i][j] = new TF1(fname, "gaus(0)", -10, 10);
	smearingFunctionGridRandomized[i][j]->SetNpx(500);
	double si = _mcConst1[i] + rand.Gaus(0.0,_mcConst1Err[i]);
	double sj = _mcConst1[j] + rand.Gaus(0.0,_mcConst1Err[j]);
	smearingFunctionGridRandomized[i][j]->SetParameters(1.0,0.0,sqrt(si*si+sj*sj));
      } // end inner loop ove eta bins
    } // end outer loop over eta bins

  } else {
    // This place should be never reached. This is just a sanity check.
    printf("ElectronEnergyScale ERROR: failed to created randomized smearing function\n");
  }

  return;
}

double ElectronEnergyScale::generateMCSmear(double eta1, double eta2){
  
  bool randomize = false;
  double result = generateMCSmearAny(eta1, eta2, randomize);

  return result;
}

double ElectronEnergyScale::generateMCSmearRandomized(double eta1, double eta2){

  double result = 0.0;

  if( !_smearingWidthRandomizationDone ){
    printf("ElectronEnergyScale ERROR: can not get randomized smear, randomization is not done\n");
    return result;
  }
  
  bool randomize = true;
  result = generateMCSmearAny(eta1, eta2, randomize);

  return result;
}

double ElectronEnergyScale::generateMCSmearAny(double eta1, double eta2, bool randomize){
  
  double result = 0;
  if( !_isInitialized ){
    printf("ElectronEnergyScale ERROR: the object is not properly initialized\n");
    return result;
  }
  
  int count = 0;
  int ibin = 0;
  int jbin = 0;
  for(int i=0; i<_nEtaBins; i++){
    if(eta1 >= _etaBinLimits[i] && eta1 < _etaBinLimits[i+1] ){
      ibin = i;
      count++;
    }
    if(eta2 >= _etaBinLimits[i] && eta2 < _etaBinLimits[i+1] ){
      jbin = i;
      count++;
    }
  }
  if(count != 2) printf("ElectronEnergyScale: Smear function ERROR\n");
 
  if( !randomize)
    result = smearingFunctionGrid[ibin][jbin]->GetRandom();
  else
    result = smearingFunctionGridRandomized[ibin][jbin]->GetRandom();

  return result;
}

void ElectronEnergyScale::print(){

  printf("\nEnergy scale corrections used:\n");
  printf("   Calibration set (%d): ",_calibrationSet);
  if( _calibrationSet == Date20110901_EPS11_default ){
    //
    printf("Date20110901_EPS11_default\n");
    printf("   Smearing function: Gaussian\n");
  } else {
    printf("   UNKNOWN. No further details available!\n");
    return;
  }
  printf("   Constants:\n");
  printf("     eta-bin      Escale-const      MC-const-1          MC-const-2          MC-const-3          MC-const-4\n");
  printf("                              %16s    %16s    %16s   %16s\n",
	 _mcConst1Name.Data(), _mcConst2Name.Data(), _mcConst3Name.Data(), _mcConst4Name.Data());
  for(int i=0; i<_nEtaBins; i++){
    printf("   %5.2f- %5.2f  %6.4f+-%5.4f",
	   _etaBinLimits[i], _etaBinLimits[i+1],
	   _dataConst[i], _dataConstErr[i]);
    if( _mcConst1 != 0 )
      printf("   %3.1e+- %3.1e", _mcConst1[i], _mcConst1Err[i]);
    if( _mcConst2 != 0 )
      printf("   %3.1e+- %3.1e", _mcConst2[i], _mcConst2Err[i]);
    if( _mcConst3 != 0 )
      printf("   %3.1e+- %3.1e", _mcConst3[i], _mcConst3Err[i]);
    if( _mcConst4 != 0 )
      printf("   %3.1e+- %3.1e", _mcConst4[i], _mcConst4Err[i]);
    printf("\n");
  }
  printf("\n");
  
  return;
}
