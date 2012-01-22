//
// This file contains implementation of methods related to applying
// electron energy scale corrections.
//

#include "../Include/ElectronEnergyScale.hh"
#include <fstream>
#include <sstream>
#include "MyTools.hh"

//------------------------------------------------------

ElectronEnergyScale::ElectronEnergyScale(CalibrationSet calibrationSet):
  _calibrationSet(calibrationSet),
  _inpFileName(),
  _isInitialized(false),
  _energyScaleCorrectionRandomizationDone(false),
  _smearingWidthRandomizationDone(false)
{
  this->init(calibrationSet);
}

//------------------------------------------------------

ElectronEnergyScale::ElectronEnergyScale(const TString &escaleTagName):
  _calibrationSet(UNDEFINED),
  _inpFileName(),
  _isInitialized(false),
  _energyScaleCorrectionRandomizationDone(false),
  _smearingWidthRandomizationDone(false)
{
  this->init(escaleTagName);
}

//------------------------------------------------------

void ElectronEnergyScale::clear() {
  if (_isInitialized) {
    if (_etaBinLimits) delete[] _etaBinLimits;
    if (_dataConst) delete[] _dataConst;
    if (_dataConstErr) delete[] _dataConstErr;
    if (_dataConstRandomized) delete[] _dataConstRandomized;
    if (_mcConst1) delete[] _mcConst1;
    if (_mcConst2) delete[] _mcConst2;
    if (_mcConst3) delete[] _mcConst3;
    if (_mcConst4) delete[] _mcConst4;
    if (_mcConst1Err) delete[] _mcConst1Err;
    if (_mcConst2Err) delete[] _mcConst2Err;
    if (_mcConst3Err) delete[] _mcConst3Err;
    if (_mcConst4Err) delete[] _mcConst4Err;
  }
}

//------------------------------------------------------

void ElectronEnergyScale::init(const TString &stringWithEScaleTagName) {
  this->clear();
  this->init(ElectronEnergyScale::DetermineCalibrationSet(stringWithEScaleTagName,&_inpFileName));
  return;
}

//------------------------------------------------------

void ElectronEnergyScale::init(CalibrationSet calSet) {
  _calibrationSet=calSet;
  _isInitialized = false;
   
  if( _calibrationSet == UNDEFINED )
    return;

  std::cout << "_calibrationSet=" << ElectronEnergyScale::CalibrationSetName(_calibrationSet,&_inpFileName) << "\n";

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

bool ElectronEnergyScale::loadInputFile(const TString &fileName, int debug) {
  ifstream fin;
  fin.open(fileName);
  if (!fin) {
    std::cout << "failed to open a file <" << fileName << ">\n";
    return kFALSE;
  }
  vector<string> lines;
  string s;
  while (getline(fin,s)) {
    lines.push_back(s);
  }
  fin.close();
  if (debug) { std::cout << "from <" << fileName << "> loaded "; PrintVec("loaded ",lines,1); }
  bool res=assignConstants(lines,debug);
  if (!res) std::cout << "error from loadInputFile(" << fileName << ")\n";
  return res;
}

//------------------------------------------------------

bool ElectronEnergyScale::assignConstants(const std::vector<string> &lines, int debug) {
  clear();
  int linesHaveNegativeEtas=-1;
  int etaDivCount=0;
  for (unsigned int i=0; i<lines.size(); i++) {
    if ((lines[i].find("g_esfWorkCase")!=std::string::npos) && 
	(lines[i].find("g_esfWorkCaseShortName")==std::string::npos)) {
      linesHaveNegativeEtas= (lines[i].find("bins on each eta side")!=std::string::npos) ? 1:0;
    }
    else if (lines[i].find("EtaDivisionCount")!=std::string::npos) {
      etaDivCount=atoi(lines[i].c_str() + lines[i].find('=') + 1);
    }
  }
  // check that lines have the expected info
  if (debug) std::cout << "linesHaveNegativeEtas=" << linesHaveNegativeEtas << ", etaDivCount=" << etaDivCount << "\n"; 
  if ((linesHaveNegativeEtas==-1) || (etaDivCount==0)) {
    std::cout << "assignConstants: the provided lines do not have a keyword 'g_esfWorkCase' or 'EtaDivisionCount'\n";
    return kFALSE;
  }
  // allocate memory
  if (linesHaveNegativeEtas==0) etaDivCount*=2;
  if (debug) std::cout << "etaDivCount=" << etaDivCount << "\n";
  _nEtaBins = etaDivCount;
  _etaBinLimits = new double[_nEtaBins+1];
  _dataConst = new double[_nEtaBins];
  _dataConstErr = new double[_nEtaBins];
  _nMCConstants=1;
  _mcConst1 = new double[_nEtaBins];
  _mcConst2 = 0;
  _mcConst3 = 0;
  _mcConst4 = 0;
  _mcConst1Err = new double[_nEtaBins];
  _mcConst2Err = 0;
  _mcConst3Err = 0;
  _mcConst4Err = 0;
  _mcConst1Name = "smear";

  assert(_etaBinLimits); 
  assert(_dataConst); assert(_dataConstErr);
  assert(_mcConst1); assert(_mcConst1Err);
  int res=ElectronEnergyScale::assignConstants(lines,_nEtaBins,_etaBinLimits,_dataConst,_dataConstErr,_mcConst1,_mcConst1Err,debug);
  if (!res) std::cout << "failed in this->assignConstants(lines)\n";
  return res;
}

//------------------------------------------------------

bool ElectronEnergyScale::assignConstants(const std::vector<string> &lines, int count, double *eta, double *scale, double *scaleErr, double *smear, double *smearErr, int debug) {
  int etaDivCount=0;
  for (unsigned int i=0; !etaDivCount && (i<lines.size()); ++i) {
    if (lines[i].find("EtaDivisionCount")!=std::string::npos) {
      etaDivCount=atoi(lines[i].c_str()+lines[i].find('=')+1);
    }
  }
  // check etaDivCount count
  if ((etaDivCount!=count) && (etaDivCount*2!=count)) {
    std::cout << "assignConstants: got lines with etaDivCount=" << etaDivCount << ", while allocation states count=" << count << "\n";
    assert(0);
  }

  // Assign etaBins
  for (unsigned int i=0; i<lines.size(); ++i) {
    if (lines[i].find("! bins")!=std::string::npos) {
      std::stringstream s(lines[i].c_str()+7);
      int idx=(etaDivCount==count) ? 0 : etaDivCount;
      for ( ; idx<=count; ++idx) {
	s >> eta[idx];
      }
      if (etaDivCount*2==count) {
	idx=2*etaDivCount;
	for (int ii=0; ii<etaDivCount; ++ii) {
	  eta[ii] = - eta[idx-ii];
	}
      }
      if (fabs(eta[0]+2.5)<1e-6) eta[0]-=1e-5;
      if (fabs(eta[count]-2.5)<1e-6) eta[count]+=1e-5;
      break;
    }
  }

  // Assign scale/smear constants
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
  case Date20120101_default: nEtaBins1=12; break;
  case CalSet_File_Gauss: nEtaBins1=-1; break;
  default:
    std::cout << "ElectronEnergyScale::initializeAllConstants: is not ready for the _calibrationSet= " << ElectronEnergyScale::CalibrationSetName(_calibrationSet,&_inpFileName) << " (" << int(_calibrationSet) << ")\n";
    return false;
  }

  if (nEtaBins1>0) {
    //
    // Allocate memory
    //
    _nEtaBins = nEtaBins1;
    _etaBinLimits = new double[nEtaBins1+1];
    _dataConst    = new double[nEtaBins1];
    _dataConstErr = new double[nEtaBins1];
    _dataConstRandomized = new double[nEtaBins1];
    switch(_calibrationSet) {
    default:
      _nMCConstants = 1;
      _mcConst1Name = "smear";
      _mcConst2Name.Clear();
      _mcConst3Name.Clear();
      _mcConst4Name.Clear();
      _mcConst1 = new double[nEtaBins1];
      _mcConst2 = 0;
      _mcConst3 = 0;
      _mcConst4 = 0;
      _mcConst1Err = new double[nEtaBins1];
      _mcConst2Err = 0;
      _mcConst3Err = 0;
      _mcConst4Err = 0;
    }
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
    _etaBinLimits[i] = -2.50001; _etaBinLimits[i+1]= 2.50001;
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

  case Date20120101_default: {
    //
    // Data of full 2011. Eta bins like for
    // for Summer 11 result. The
    // constants are NOT symmetric about eta = 0.
    //
    const int nEtaBins = 12;
    //const double etaBinLimits[nEtaBins+1] = 
    //  {-2.50001, -2.0, -1.5, -1.2, -0.8, -0.4, 0.0, 0.4, 0.8, 1.2, 1.5, 2.0, 2.50001};
    std::vector<string> lines;
    lines.push_back("! Date 2012 Jan 19\n");
    lines.push_back("! g_esfWorkCase=15 (6 bins on each eta side)");
    lines.push_back("! g_esfWorkCaseShortName= 6binNegs");
    lines.push_back("! g_esfFitModel=1 (fit model Gauss)");
    lines.push_back("scaling sqrt");
    lines.push_back("! bins  -2.50 -2.00 -1.50 -1.20 -0.80 -0.40 0.00 0.40 0.80 1.20 1.50 2.00 2.50");
    lines.push_back("MCOverData=15503.114069");
    lines.push_back("EtaDivisionCount=12");
    lines.push_back("ScalingFactorsCount=12");
    lines.push_back("SmearingFactorsCount=12");
    lines.push_back("scale_0      1.01214 -0.000303694");
    lines.push_back("scale_1     0.989212 -0.000268224");
    lines.push_back("scale_2      1.02211 -0.000280481");
    lines.push_back("scale_3      1.01145 -0.000168768");
    lines.push_back("scale_4      1.00602  -0.00014565");
    lines.push_back("scale_5      1.00573 -0.000144208");
    lines.push_back("scale_6      1.00293   -0.0001399");
    lines.push_back("scale_7      1.00375 -0.000146988");
    lines.push_back("scale_8      1.01098 -0.000170304");
    lines.push_back("scale_9      1.02083 -0.000283934");
    lines.push_back("scale_10     0.989319 -0.000270072");
    lines.push_back("scale_11      1.01284 -0.000304156");
    lines.push_back("smear_0      1.94897   -0.0226282");
    lines.push_back("smear_1      1.56441   -0.0266465");
    lines.push_back("smear_2       1.1956   -0.0284747");
    lines.push_back("smear_3     0.768533   -0.0221135");
    lines.push_back("smear_4     0.432683    -0.028314");
    lines.push_back("smear_5     0.525586   -0.0239693");
    lines.push_back("smear_6     0.455353   -0.0257488");
    lines.push_back("smear_7     0.471485   -0.0277377");
    lines.push_back("smear_8     0.742029   -0.0229296");
    lines.push_back("smear_9      1.19892   -0.0287477");
    lines.push_back("smear_10      1.49028   -0.0275794");
    lines.push_back("smear_11      2.03342   -0.0219498");

    if (nEtaBins1!=nEtaBins) assert(0);

    assert(_etaBinLimits); 
    assert(_dataConst); assert(_dataConstErr);
    assert(_mcConst1); assert(_mcConst1Err);
    //for(int i=0; i<nEtaBins+1; i++) _etaBinLimits[i] = etaBinLimits[i];
    if (!assignConstants(lines, nEtaBins,_etaBinLimits,_dataConst,_dataConstErr,_mcConst1,_mcConst1Err)) assert(0);
  }
    break;

  case CalSet_File_Gauss: {
    if (_inpFileName.Length()==0) {
      std::cout << "ElectronEnergyScale::initializeAllConstants. Calibration set CalSet_File_Gauss requires input file to be set\n";
      assert(0);
    }
    int debug=0;
    success = loadInputFile(_inpFileName,debug);
    break;
  }

  default:
    std::cout << "ElectronEnergyScale::initializeAllConstants: is not ready for the _calibrationSet= " << ElectronEnergyScale::CalibrationSetName(_calibrationSet,&_inpFileName) << " (" << int(_calibrationSet) << ") [3]\n";
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
      case Date20120101_default:
      case CalSet_File_Gauss: {
	if(_mcConst1 == 0) continue;
	smearingFunctionGrid[i][j] = new TF1(fname, "gaus(0)", -10, 10);
	smearingFunctionGrid[i][j]->SetNpx(500);
	double si = _mcConst1[i];
	double sj = _mcConst1[j];
	double sij=sqrt(si*si+sj*sj);
	smearingFunctionGrid[i][j]->SetParameters(1.0/(sij*sqrt(8*atan(1))),0.0,sij);
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
  if (_calibrationSet==UNCORRECTED) return;

  for(int i=0; i<_nEtaBins; i++){
    _dataConstRandomized[i] = _dataConst[i] + rand.Gaus(0.0,_dataConstErr[i]);
  }

  return;
}

//------------------------------------------------------

double ElectronEnergyScale::getEnergyScaleCorrection(double eta) const {

  double result = 1.0;
  bool randomize = false;
  if (_calibrationSet != UNCORRECTED) {
    result = getEnergyScaleCorrectionAny(eta,randomize);
  }

  return result;
}

//------------------------------------------------------

double ElectronEnergyScale::getEnergyScaleCorrectionRandomized(double eta) const {

  double result = 1.0;
  if (_calibrationSet == UNCORRECTED) return result;

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
  if (_calibrationSet == UNCORRECTED) return result;

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
  switch( _calibrationSet ) {
  case UNCORRECTED: break;

  case Date20110901_EPS11_default:
  case Date20120101_default:
  case CalSet_File_Gauss: {

    for( int i=0; i<_nEtaBins; i++){
      for( int j=0; j<_nEtaBins; j++){
	TString fname = TString::Format("smearing_function_randomized_%03d_%03d", i, j);
	smearingFunctionGridRandomized[i][j] = new TF1(fname, "gaus(0)", -10, 10);
	smearingFunctionGridRandomized[i][j]->SetNpx(500);
	double si = _mcConst1[i] + rand.Gaus(0.0,_mcConst1Err[i]);
	double sj = _mcConst1[j] + rand.Gaus(0.0,_mcConst1Err[j]);
	double sij= sqrt(si*si+sj*sj);
	smearingFunctionGridRandomized[i][j]->SetParameters(1.0/(sij*sqrt(8*atan(1))),0.0,sij);
	if (i>j) { smearingFunctionGridRandomized[i][j]->SetParameters(smearingFunctionGridRandomized[j][i]->GetParameters()); }
      } // end inner loop over eta bins
    } // end outer loop over eta bins
  }
    break;
  default:
    // This place should be never reached. This is just a sanity check.
    printf("ElectronEnergyScale ERROR: failed to created randomized smearing function\n");
  }

  return;
}

//------------------------------------------------------

double ElectronEnergyScale::generateMCSmear(double eta1, double eta2) const {
  
  bool randomize = false;
  if (_calibrationSet == UNCORRECTED) return 0.0;

  double result = generateMCSmearAny(eta1, eta2, randomize);

  return result;
}

//------------------------------------------------------

double ElectronEnergyScale::generateMCSmearRandomized(double eta1, double eta2) const {

  double result = 0.0;
  if (_calibrationSet == UNCORRECTED) return 0.0;

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
  if (_calibrationSet == UNCORRECTED) return result;

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
  printf("   Calibration set (%d): %s\n", _calibrationSet, ElectronEnergyScale::CalibrationSetName(this->_calibrationSet,&this->_inpFileName).Data());
  printf("   Smearing function: %s\n",ElectronEnergyScale::CalibrationSetFunctionName(this->_calibrationSet).Data());
  printf("   Constants:\n");
  printf("     eta-bin      Escale-const      MC-const-1");
  if ( _mcConst2 ) printf("          MC-const-2");
  if ( _mcConst3 ) printf("          MC-const-3");
  if ( _mcConst4 ) printf("          MC-const-4");
  printf("\n");
  printf("                              %16s    %16s    %16s   %16s\n",
	 _mcConst1Name.Data(), _mcConst2Name.Data(), _mcConst3Name.Data(), _mcConst4Name.Data());
  for(int i=0; i<_nEtaBins; i++){
    printf("   %5.2f- %5.2f  %6.4f+-%5.4f",
    //printf("   %8.5f- %8.5f  %6.4f+-%5.4f",
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

ElectronEnergyScale::CalibrationSet ElectronEnergyScale::DetermineCalibrationSet(const TString &escaleTagName_orig, TString *inputFileName) {
  ElectronEnergyScale::CalibrationSet calibrationSet  = ElectronEnergyScale::UNDEFINED;
  TString fileName;
  TString escaleTagName;
  {
    //std::cout << "entered DetermineCalibrationSet(\"" << escaleTagName_orig << "\", " << ((inputFileName) ? "ptr":"NULL") << ")\n";
    Ssiz_t pos = escaleTagName_orig.Index("#");
    if (pos>0) escaleTagName=escaleTagName_orig(0,pos-1); 
    else escaleTagName=escaleTagName_orig;
    //std::cout << "escaleTagName_orig={" << escaleTagName_orig << "}, escaleTagName={" << escaleTagName << "}\n";
  }
  if ( escaleTagName.Contains("UNCORRECTED")) {
    calibrationSet = ElectronEnergyScale::UNCORRECTED;
  }
  else if ( escaleTagName.Contains("Date20110901_EPS11_default")) {
    calibrationSet = ElectronEnergyScale::Date20110901_EPS11_default;
  }
  else if ( escaleTagName.Contains("Date20120101_default")) {
    calibrationSet = ElectronEnergyScale::Date20120101_default;
  }
  else if ( escaleTagName.Contains("Date20120101_Gauss_6bins")) {
    calibrationSet = ElectronEnergyScale::CalSet_File_Gauss;
    fileName="../root_files/constants/testESF_6bins_Gauss_20120119.inp";
  }
  else if ( escaleTagName.Contains("Date20120101_Gauss_6binNegs")) {
    calibrationSet = ElectronEnergyScale::CalSet_File_Gauss;
    fileName="../root_files/constants/testESF_6binNegs_Gauss_20120119.inp";
  }
  else if (escaleTagName.Contains("File")) {
    if (escaleTagName.Contains("File_Gauss") ||
	escaleTagName.Contains("FileGauss")) {      
      Ssiz_t pos=escaleTagName.Index("File");
      while ((pos<escaleTagName.Length()) && (escaleTagName[pos]!=' ')) pos++;
      while ((pos<escaleTagName.Length()) && (escaleTagName[pos]==' ')) pos++;
      if ((pos==escaleTagName.Length()) || (escaleTagName[pos]=='\n')) {
	std::cout << "since escaleTagName contains 'FileGauss', the keyword has to be followed by a file name\n";
      }
      else {
	calibrationSet = ElectronEnergyScale::CalSet_File_Gauss;
	fileName=escaleTagName(pos,escaleTagName.Length()-pos);
      }
    }
  }
  else{
    printf("Failed to match escale calibration. Tag: >>%s<<\n", escaleTagName.Data());
    assert(0); 
  }
    
  if (fileName.Length()) {
    if (inputFileName) *inputFileName=fileName;
    else {
      std::cout << "warning: ElectronEnergyScale::DetermineCalibrationSet located a keyword meaning the file name will be provided, yet no pointer to the container was supplied\n";
    }
  }
  //std::cout << "DetermineCalibrationSet(" << escaleTagName << ") returns " << ElectronEnergyScale::CalibrationSetName(calibrationSet,&fileName) << " (" << int(calibrationSet) << ")\n";
  return calibrationSet;
}

//------------------------------------------------------

TString ElectronEnergyScale::CalibrationSetName(ElectronEnergyScale::CalibrationSet escaleTag, const TString *fileName) {
  TString name="UNDEFINED";
  switch(escaleTag) {
  case ElectronEnergyScale::UNDEFINED: break;
  case ElectronEnergyScale::UNCORRECTED: name="UNCORRECTED"; break;
  case ElectronEnergyScale::Date20110901_EPS11_default: name="Date20110910_EPS11_default"; break;
  case ElectronEnergyScale::Date20120101_default: name="Date20120101_default"; break;
  case ElectronEnergyScale::CalSet_File_Gauss: 
    name="FileGauss(";
    if (fileName) name+=(*fileName);
    name+=")";
    break;
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
  case ElectronEnergyScale::Date20120101_default: break; // Gauss
  case ElectronEnergyScale::CalSet_File_Gauss: break; // Gauss
  default:
    name="CalibrationSetFunctionName_undetermined";
  }
  return name;
}

//------------------------------------------------------
