double findEnergyScaleCorrection(double eta){

  double corr = 1.0;

  // Energy scale corrections from Duncan Ralph derived July 2011
  // based on 204 pb-1 May ReReco and ~650 pb-1 of prompt reco
  // as well as 41X powheg MC  
//   const int nEtaBins = 10;
//   double corrEtaBinLimits[nEtaBins+1] = 
//     {-2.50001, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.50001};
//   double corrValues[nEtaBins] = 
//     {0.9890, 1.0250, 0.9932, 1.0029, 0.9986, 1.0027, 1.0036, 0.9887, 1.0209, 0.9812};
//   double corrErrors[nEtaBins] = 
//     {0.0006, 0.0007, 0.0003, 0.0003, 0.0007, 0.0000, 0.0002, 0.0012, 0.0002, 0.0028};

  // Energy scale corrections from Andrius drived July 2011
  // on 204 pb-1 May 10 ReReco and 800 pb-1 PromptReco, matched
  // against Summer11 powheg MC
  // Note that Andrius has 6 eta bins in absolute  eta, so value "i"
  // below is equal to value "N-i"
  const int nEtaBins = 12;
  double corrEtaBinLimits[nEtaBins+1] = 
    {-2.50001, -2.0, -1.5, -1.2, -0.8, -0.4, 0.0, 0.4, 0.8, 1.2, 1.5, 2.0, 2.50001};
  double corrValues[nEtaBins] = 
    {1.04642, 1.00187, 1.01556, 1.00500, 1.00093, 1.00149, 1.00149, 1.00093, 1.00500, 1.01556, 1.00187, 1.04642};
  double corrErrors[nEtaBins] = 
    {4.28928e-04,3.39718e-04,4.89342e-04,5.80480e-05,1.21192e-05,1.27489e-04,1.27489e-04,1.21192e-05,5.80480e-05,4.89342e-04,3.39718e-04,4.28928e-04};

  for(int i=0; i<nEtaBins; i++){
    if(eta >= corrEtaBinLimits[i] && eta < corrEtaBinLimits[i+1] ){
      corr = corrValues[i];
      break;
    }
  }
  
  return corr;

}

double extraSmearingSigma(double eta){

  double smear = 0.0;

  // Energy scale corrections from Duncan Ralph derived July 2011
  // based on 204 pb-1 May ReReco and ~650 pb-1 of prompt reco
  // as well as 41X powheg MC  
//   const int nEtaBins = 10;
//   double smearEtaBinLimits[nEtaBins+1] = 
//     {-2.50001, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.50001};
//   double smearValues[nEtaBins] = 
//     {1.99, 1.36, 0.80, 0.63, 0.55, 0.59, 0.62, 1.27, 1.28, 2.00};
//   double smearErrors[nEtaBins] =
//     {0.08, 0.14, 0.07, 0.05, 0.09, 0.05, 0.03, 0.09, 0.07, 0.15};

  // Energy scale corrections+extra smear from Andrius drived July 2011
  // on 204 pb-1 May 10 ReReco and 800 pb-1 PromptReco, matched
  // against Summer11 powheg MC
  // Note that Andrius has 6 eta bins in absolute  eta, so value "i"
  // below is equal to value "N-i"
  const int nEtaBins = 12;
  double smearEtaBinLimits[nEtaBins+1] = 
    {-2.50001, -2.0, -1.5, -1.2, -0.8, -0.4, 0.0, 0.4, 0.8, 1.2, 1.5, 2.0, 2.50001};
  double smearValues[nEtaBins] = 
    {2.05888e+00,1.46747e+00,1.14861e+00,7.63770e-01,5.04140e-01,5.27258e-01,5.27258e-01,5.04140e-01,7.63770e-01,1.14861e+00,1.46747e+00,2.05888e+00};
  double smearErrors[nEtaBins] =
    {2.85889e-02,3.85260e-02,4.26451e-02,3.22979e-02,3.76972e-02,3.32377e-02,3.32377e-02,3.76972e-02,3.22979e-02,4.26451e-02,3.85260e-02,2.85889e-02};

  for(int i=0; i<nEtaBins; i++){
    if(eta >= smearEtaBinLimits[i] && eta < smearEtaBinLimits[i+1] ){
      smear = smearValues[i];
      break;
    }
  }

  return smear;
  
}

