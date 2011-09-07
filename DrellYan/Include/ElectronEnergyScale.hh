#ifndef ElectronEnergyScale
#define ElectronEnergyScale

namespace escale{

  // Energy scale corrections and MC extra smearing from Andrius drived July 2011
  // on 204 pb-1 May 10 ReReco and 800 pb-1 PromptReco, matched
  // against Summer11 powheg MC
  // Note that Andrius has 6 eta bins in absolute  eta, so value "i"
  // below is equal to value "N-i"
  const int nEtaBins = 12;
  const double escaleEtaBinLimits[nEtaBins+1] = 
    {-2.50001, -2.0, -1.5, -1.2, -0.8, -0.4, 0.0, 0.4, 0.8, 1.2, 1.5, 2.0, 2.50001};
  const double corrValues[nEtaBins] = 
    {1.04642, 1.00187, 1.01556, 1.00500, 1.00093, 1.00149, 1.00149, 1.00093, 1.00500, 1.01556, 1.00187, 1.04642};
  const double corrErrors[nEtaBins] = 
    {4.28928e-04,3.39718e-04,4.89342e-04,5.80480e-05,1.21192e-05,1.27489e-04,1.27489e-04,1.21192e-05,5.80480e-05,4.89342e-04,3.39718e-04,4.28928e-04};
  const double smearValues[nEtaBins] = 
    {2.05888e+00,1.46747e+00,1.14861e+00,7.63770e-01,5.04140e-01,5.27258e-01,5.27258e-01,5.04140e-01,7.63770e-01,1.14861e+00,1.46747e+00,2.05888e+00};
  const double smearErrors[nEtaBins] =
    {2.85889e-02,3.85260e-02,4.26451e-02,3.22979e-02,3.76972e-02,3.32377e-02,3.32377e-02,3.76972e-02,3.22979e-02,4.26451e-02,3.85260e-02,2.85889e-02};

  double findEnergyScaleCorrection(double eta){
    
    double corr = 1.0;
    
    for(int i=0; i<nEtaBins; i++){
      if(eta >= escaleEtaBinLimits[i] && eta < escaleEtaBinLimits[i+1] ){
	corr = corrValues[i];
	break;
      }
    }
    
    return corr;
    
  }
  
  double extraSmearingSigma(double eta){
    
    double smear = 0.0;
    
    for(int i=0; i<nEtaBins; i++){
      if(eta >= escaleEtaBinLimits[i] && eta < escaleEtaBinLimits[i+1] ){
	smear = smearValues[i];
	break;
      }
    }
    
    return smear;
    
  }
  
  // This function is meant to be used for systematic
  // studies where we want to smear the extra smearing constant
  // within its errors.
  double extraSmearingSigmaShifted(double eta, TVectorD &shift){
    
    double smear = 0.0;
    
    for(int i=0; i<nEtaBins; i++){
      if(eta >= escaleEtaBinLimits[i] && eta < escaleEtaBinLimits[i+1] ){
	smear = smearValues[i] + shift[i] * smearErrors[i];
	break;
      }
    }
    
    return smear;
    
  }
  
} // end of namespace
#endif
