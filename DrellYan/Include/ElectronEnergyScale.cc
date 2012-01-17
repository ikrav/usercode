//
// This file contains implementation of methods related to applying
// electron energy scale corrections.
//

#include "../Include/ElectronEnergyScale.hh"

//------------------------------------------------------

ElectronEnergyScale::ElectronEnergyScale(CalibrationSet calibrationSet):
  _calibrationSet(calibrationSet),
  _isInitialized(false),
  _energyScaleCorrectionRandomizationDone(false),
  _smearingWidthRandomizationDone(false)
{

  if(_calibrationSet == UNDEFINED )
    return;

  if( !initializeAllConstants())
    return;

  if( !initializeExtraSmearingFunction())
    return;

  _isInitialized = true;
  return;
}

//------------------------------------------------------

ElectronEnergyScale::ElectronEnergyScale(const TString &escaleTagName):
  _calibrationSet(UNDEFINED),
  _isInitialized(false),
  _energyScaleCorrectionRandomizationDone(false),
  _smearingWidthRandomizationDone(false)
{

  _calibrationSet=ElectronEnergyScale::DetermineCalibrationSet(escaleTagName);

  if(_calibrationSet == UNDEFINED )
    return;

  std::cout << "_calibrationSet=" << ElectronEnergyScale::CalibrationSetName(_calibrationSet) << "\n";

  if( !initializeAllConstants()) {
    std::cout << "failed to initialize\n";
    return;
  }

  if( !initializeExtraSmearingFunction()) {
    std::cout << "failed to prepare extra smearing function\n";
    return;
  }

  _isInitialized = true;
  return;
}

//------------------------------------------------------

bool AssignConstants(const std::vector<string> &lines, int count, double *scale, double *scaleErr, double *smear, double *smearErr, int debug=0) {
  int etaDivCount=0;
  for (unsigned int i=0; !etaDivCount && (i<lines.size()); ++i) {
    if (lines[i].find("EtaDivisionCount")!=std::string::npos) {
      etaDivCount=atoi(lines[i].c_str()+lines[i].find('=')+1);
    }
  }
  if ((etaDivCount!=count) && (etaDivCount*2!=count)) {
    std::cout << "AssignConstants: got lines with etaDivCount=" << etaDivCount << ", while allocation states count=" << count << "\n";
    assert(0);
  }
  double *d,*derr;
  for (unsigned int i=0; i<lines.size(); ++i) {
    d=NULL; derr=NULL;
    if (lines[i].find("scale_")!=std::string::npos) {
      d=scale; derr=scaleErr;
    }
    else if (lines[i].find("smear_")!=std::string::npos) {
      d=smear; derr=smearErr;
    }
    if (d) {
      const char *s=lines[i].c_str();
      int idx=atoi(s + lines[i].find('_') + 1);
      double val=atof(s + lines[i].find(' '));
      double valErr=atof(s + lines[i].find(' ',lines[i].find('.')));
      if (etaDivCount==count) {
	d[idx]=val; derr[idx]=valErr;
      }
      else if (etaDivCount*2==count) {
	d[etaDivCount-idx-1]=val;
	derr[etaDivCount-idx-1]=valErr;
	d[etaDivCount+idx]=val;
	derr[etaDivCount+idx]=valErr;
      }
      else assert(0);
    }
  }
  if (debug) {
    std::cout << "got \n";
    for (unsigned int i=0; i<lines.size(); ++i) std::cout << "--> " << lines[i] << "\n";
    std::cout << "derived \n";
    for (int i=0; i<count; ++i) {
      std::cout << "scale_" << i << "  " << scale[i] << " " << scaleErr[i] << "\n";
    }
    for (int i=0; i<count; ++i) {
      std::cout << "smear_" << i << "  " << smear[i] << " " << smearErr[i] << "\n";
    }
  }
  return true;
}

//------------------------------------------------------

bool ElectronEnergyScale::initializeAllConstants(){
  
  bool success = true;
  int nEtaBins1=0;
  //
  // Determine the number of bins
  //
  switch(_calibrationSet) {
  case UNCORRECTED: nEtaBins1=1; break;
  case Date20110901_EPS11_default: nEtaBins1=12; break;
  case Date20120101_Gauss_6bins: nEtaBins1=12; break;
  case Date20120101_Gauss_6binNegs: nEtaBins1=12; break;
  default:
    std::cout << "ElectronEnergyScale::initializeAllConstants: is not ready for the _calibrationSet= " << ElectronEnergyScale::CalibrationSetName(_calibrationSet) << " (" << int(_calibrationSet) << ")\n";
    return false;
  }

  //
  // Allocate memory
  //
  _nEtaBins = nEtaBins1;
  _etaBinLimits = (double *)malloc(sizeof(double)*(nEtaBins1+1));
  _dataConst    = (double *)malloc(sizeof(double)*nEtaBins1);
  _dataConstErr = (double *)malloc(sizeof(double)*nEtaBins1);
  _dataConstRandomized = (double *)malloc(sizeof(double)*nEtaBins1);
  switch(_calibrationSet) {
  default:
    _nMCConstants = 1;
    _mcConst1Name = "smear";
    _mcConst2Name = "<none>";
    _mcConst3Name = "<none>";
    _mcConst4Name = "<none>";
    _mcConst1 = (double *)malloc(sizeof(double)*nEtaBins1);
    _mcConst2 = 0;
    _mcConst3 = 0;
    _mcConst4 = 0;
    _mcConst1Err = (double *)malloc(sizeof(double)*nEtaBins1);
    _mcConst2Err = 0;
    _mcConst3Err = 0;
    _mcConst4Err = 0;
  }
  
  // 
  // Assign values
  // 
  switch( _calibrationSet ) { 
  case UNCORRECTED: {    
    if (nEtaBins1!=1) assert(0);
    assert(_etaBinLimits); 
    assert(_dataConst); assert(_dataConstErr);
    assert(_mcConst1); assert(_mcConst1Err);
    int i=0;
    _etaBinLimits[i] = -2.5; _etaBinLimits[i+1]= 2.5;
    _dataConst   [i] = 1.0;
    _dataConstErr[i] = 0.0;
    _mcConst1    [i] = 0.;
    _mcConst1Err [i] = 0.;
  }
    break;
  case Date20110901_EPS11_default: {
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

    if (nEtaBins1!=nEtaBins) assert(0);
    assert(_etaBinLimits); 
    assert(_dataConst); assert(_dataConstErr);
    assert(_mcConst1); assert(_mcConst1Err);
    for(int i=0; i<nEtaBins; i++){
      _etaBinLimits[i] = etaBinLimits[i];
      _dataConst   [i] = corrValues[i];
      _dataConstErr[i] = corrErrors[i];
      _mcConst1    [i] = smearValues[i];
      _mcConst1Err [i] = smearErrors[i];
    }
    _etaBinLimits[nEtaBins] = etaBinLimits[nEtaBins];
  }
    break;

  case Date20120101_Gauss_6bins: {
    //
    // Data of full 2011. Eta bins like for
    // for Summer 11 result. Note that the 
    // constants are symmetric about eta = 0.
    //
    const int nEtaBins = 12;
    const double etaBinLimits[nEtaBins+1] = 
      {-2.50001, -2.0, -1.5, -1.2, -0.8, -0.4, 0.0, 0.4, 0.8, 1.2, 1.5, 2.0, 2.50001};
    std::vector<string> lines;
    lines.push_back("! g_esfWorkCase=6 (6 bins)");
    lines.push_back("! g_esfFitModel=1 (fit model Gauss)");
    lines.push_back("scaling sqrt");
    lines.push_back("! bins  0.00 0.40 0.80 1.20 1.50 2.00 2.50");
    lines.push_back("MCOverData=53382.626340");
    lines.push_back("EtaDivisionCount=6");
    lines.push_back("ScalingFactorsCount=6");
    lines.push_back("SmearingFactorsCount=6");
    lines.push_back("scale_0      1.00427 -9.33125e-05");
    lines.push_back("scale_1      1.00487  -9.2255e-05");
    lines.push_back("scale_2      1.01121  -0.00012073");
    lines.push_back("scale_3      1.02144 -0.000199873");
    lines.push_back("scale_4     0.989138 -0.000184884");
    lines.push_back("scale_5      1.01263 -0.000218029");
    lines.push_back("smear_0     0.485588   -0.0169862");
    lines.push_back("smear_1     0.435977   -0.0204939");
    lines.push_back("smear_2     0.739175   -0.0163163");
    lines.push_back("smear_3      1.17939   -0.0205025");
    lines.push_back("smear_4      1.50625   -0.0194033");
    lines.push_back("smear_5      1.99576   -0.0155175");

    if (nEtaBins1!=nEtaBins) assert(0);

    assert(_etaBinLimits); 
    assert(_dataConst); assert(_dataConstErr);
    assert(_mcConst1); assert(_mcConst1Err);
    for(int i=0; i<nEtaBins+1; i++) _etaBinLimits[i] = etaBinLimits[i];
    if (!AssignConstants(lines, nEtaBins,_dataConst,_dataConstErr,_mcConst1,_mcConst1Err)) assert(0);
  }
    break;


  case Date20120101_Gauss_6binNegs: {
    //
    // Data of full 2011. Eta bins like for
    // for Summer 11 result. The
    // constants are NOT symmetric about eta = 0.
    //
    const int nEtaBins = 12;
    const double etaBinLimits[nEtaBins+1] = 
      {-2.50001, -2.0, -1.5, -1.2, -0.8, -0.4, 0.0, 0.4, 0.8, 1.2, 1.5, 2.0, 2.50001};
    std::vector<string> lines;
    lines.push_back("! g_esfWorkCase=12 (6 bins on each eta side)");
    lines.push_back("! g_esfFitModel=1 (fit model Gauss)");
    lines.push_back("scaling sqrt");
    lines.push_back("! bins  -2.50 -2.00 -1.50 -1.20 -0.80 -0.40 0.00 0.40 0.80 1.20 1.50 2.00 2.50");
    lines.push_back("MCOverData=15565.640556");
    lines.push_back("EtaDivisionCount=12");
    lines.push_back("ScalingFactorsCount=12");
    lines.push_back("SmearingFactorsCount=12");
    lines.push_back("scale_0      1.01223 -0.000283215");
    lines.push_back("scale_1     0.989088 -0.000210782");
    lines.push_back("scale_2      1.02209 -0.000250671");
    lines.push_back("scale_3       1.0114 -0.000155979");
    lines.push_back("scale_4      1.00601 -0.000155017");
    lines.push_back("scale_5       1.0057 -0.000109885");
    lines.push_back("scale_6      1.00293 -0.000137715");
    lines.push_back("scale_7      1.00371 -2.22722e-05");
    lines.push_back("scale_8      1.01091  -8.8298e-05");
    lines.push_back("scale_9      1.02074 -0.000232632");
    lines.push_back("scale_10     0.989139 -0.000200525");
    lines.push_back("scale_11      1.01299 -0.000299865");
    lines.push_back("smear_0      1.94993   -0.0222418");
    lines.push_back("smear_1      1.55407   -0.0263524");
    lines.push_back("smear_2      1.19115   -0.0280501");
    lines.push_back("smear_3      0.76452   -0.0219793");
    lines.push_back("smear_4     0.429524   -0.0282913");
    lines.push_back("smear_5      0.52105   -0.0239832");
    lines.push_back("smear_6     0.454759   -0.0256199");
    lines.push_back("smear_7     0.470073   -0.0273218");
    lines.push_back("smear_8     0.736159   -0.0226141");
    lines.push_back("smear_9      1.19079   -0.0283186");
    lines.push_back("smear_10      1.48027   -0.0271734");
    lines.push_back("smear_11       2.0404   -0.0214981");

    if (nEtaBins1!=nEtaBins) assert(0);

    assert(_etaBinLimits); 
    assert(_dataConst); assert(_dataConstErr);
    assert(_mcConst1); assert(_mcConst1Err);
    for(int i=0; i<nEtaBins+1; i++) _etaBinLimits[i] = etaBinLimits[i];
    if (!AssignConstants(lines, nEtaBins,_dataConst,_dataConstErr,_mcConst1,_mcConst1Err)) assert(0);
  }
    break;

  default:
    std::cout << "ElectronEnergyScale::initializeAllConstants: is not ready for the _calibrationSet= " << ElectronEnergyScale::CalibrationSetName(_calibrationSet) << " (" << int(_calibrationSet) << ") [3]\n";
    return false;
  }
  
  return success;
}

//------------------------------------------------------

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
      switch(_calibrationSet) {
      case UNCORRECTED: break;
      case Date20110901_EPS11_default:
      case Date20120101_Gauss_6bins:
      case Date20120101_Gauss_6binNegs: {
	if(_mcConst1 == 0) continue;
	smearingFunctionGrid[i][j] = new TF1(fname, "gaus(0)", -10, 10);
	smearingFunctionGrid[i][j]->SetNpx(500);
	double si = _mcConst1[i];
	double sj = _mcConst1[j];
	smearingFunctionGrid[i][j]->SetParameters(1.0,0.0,sqrt(si*si+sj*sj));
      }
	break;
      default:
	success = false;
      }
    } // end inner loop over eta bins
  } // end outer loop over eta bins

  return success;
}

//------------------------------------------------------

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

//------------------------------------------------------

double ElectronEnergyScale::getEnergyScaleCorrection(double eta) const {

  double result = 1.0;
  bool randomize = false;
  result = getEnergyScaleCorrectionAny(eta,randomize);

  return result;
}

//------------------------------------------------------

double ElectronEnergyScale::getEnergyScaleCorrectionRandomized(double eta) const {

  double result = 1.0;
  if( !_energyScaleCorrectionRandomizationDone ){
    printf("ElectronEnergyScale ERROR: can not get randomized escale, randomization is not done\n");
    return result;
  }

  bool randomize = true;
  result = getEnergyScaleCorrectionAny(eta,randomize);

  return result;
}

//------------------------------------------------------

double ElectronEnergyScale::getEnergyScaleCorrectionAny(double eta, bool randomize) const {

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

//------------------------------------------------------

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
	if (i>j) { smearingFunctionGridRandomized[i][j]->SetParameters(smearingFunctionGridRandomized[j][i]->GetParameters()); }
      } // end inner loop over eta bins
    } // end outer loop over eta bins

  } else {
    // This place should be never reached. This is just a sanity check.
    printf("ElectronEnergyScale ERROR: failed to created randomized smearing function\n");
  }

  return;
}

//------------------------------------------------------

double ElectronEnergyScale::generateMCSmear(double eta1, double eta2) const {
  
  bool randomize = false;
  double result = generateMCSmearAny(eta1, eta2, randomize);

  return result;
}

//------------------------------------------------------

double ElectronEnergyScale::generateMCSmearRandomized(double eta1, double eta2) const {

  double result = 0.0;

  if( !_smearingWidthRandomizationDone ){
    printf("ElectronEnergyScale ERROR: can not get randomized smear, randomization is not done\n");
    return result;
  }
  
  bool randomize = true;
  result = generateMCSmearAny(eta1, eta2, randomize);

  return result;
}

//------------------------------------------------------

double ElectronEnergyScale::generateMCSmearAny(double eta1, double eta2, bool randomize) const {
  
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

//------------------------------------------------------

void ElectronEnergyScale::print() const {

  printf("\nEnergy scale corrections used:\n");
  printf("   Calibration set (%d): %s\n", _calibrationSet, ElectronEnergyScale::CalibrationSetName(this->_calibrationSet).Data());
  printf("   Smearing function: %s\n",ElectronEnergyScale::CalibrationSetFunctionName(this->_calibrationSet).Data());
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

//------------------------------------------------------
//------------------------------------------------------

ElectronEnergyScale::CalibrationSet ElectronEnergyScale::DetermineCalibrationSet(const TString &escaleTagName) {
  ElectronEnergyScale::CalibrationSet calibrationSet  = ElectronEnergyScale::UNDEFINED;
  if ( escaleTagName == TString("UNCORRECTED")) {
    calibrationSet = ElectronEnergyScale::UNCORRECTED;
  }
  else if ( escaleTagName == TString("Date20110901_EPS11_default")) {
    calibrationSet = ElectronEnergyScale::Date20110901_EPS11_default;
  }
  else if ( escaleTagName == TString("Date20120101_Gauss_6bins")) {
    calibrationSet = ElectronEnergyScale::Date20120101_Gauss_6bins;
  }
  else if ( escaleTagName == TString("Date20120101_Gauss_6binNegs")) {
    calibrationSet = ElectronEnergyScale::Date20120101_Gauss_6binNegs;
  }
  else{
    printf("Failed to match escale calibration. Tag: >>%s<<\n", escaleTagName.Data());
    assert(0); 
  }
  //std::cout << "DetermineCalibrationSet(" << escaleTagName << ") returns " << ElectronEnergyScale::CalibrationSetName(calibrationSet) << " (" << int(calibrationSet) << ")\n";
  return calibrationSet;
}

//------------------------------------------------------

TString ElectronEnergyScale::CalibrationSetName(ElectronEnergyScale::CalibrationSet escaleTag) {
  TString name="UNDEFINED";
  switch(escaleTag) {
  case ElectronEnergyScale::UNDEFINED: break;
  case ElectronEnergyScale::UNCORRECTED: name="UNCORRECTED"; break;
  case ElectronEnergyScale::Date20110901_EPS11_default: name="Date20110910_EPS11_default"; break;
  case ElectronEnergyScale::Date20120101_Gauss_6bins: name="Date20120101_Gauss_6bins"; break;
  case ElectronEnergyScale::Date20120101_Gauss_6binNegs: name="Date20120101_Gauss_6binNegs"; break;
  default:
    name="CalibrationSetName_undetermined";
  }
  return name;
}

//------------------------------------------------------

TString ElectronEnergyScale::CalibrationSetFunctionName(ElectronEnergyScale::CalibrationSet escaleTag) {
  TString name="Gauss";
  switch (escaleTag) {
  case ElectronEnergyScale::UNDEFINED: name="undefined"; break;
  case ElectronEnergyScale::UNCORRECTED: name="uncorrected"; break;
  case ElectronEnergyScale::Date20110901_EPS11_default: break; // Gauss
  case ElectronEnergyScale::Date20120101_Gauss_6bins: break; // Gauss
  case ElectronEnergyScale::Date20120101_Gauss_6binNegs: break; // Gauss
  default:
    name="CalibrationSetFunctionName_undetermined";
  }
  return name;
}

//------------------------------------------------------
