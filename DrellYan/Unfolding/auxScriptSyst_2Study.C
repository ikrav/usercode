// The script will be superceded by evaluateUnfoldingSyst.sh

{

  //gROOT->ProcessLine(".L plotDYUnfoldingMatrix.C+");

  //plotDYUnfoldingMatrix("../config_files/summer11mc.input",DYTools::FSR_STUDY,1,1.05,-1);
  //plotDYUnfoldingMatrix("../config_files/summer11mc.input",DYTools::FSR_STUDY,1,0.95,-1);

  //for (int i=1; i<=20; i++)
  //{
  // plotDYUnfoldingMatrix("../config_files/summer11mc.input",DYTools::RESOLUTION_STUDY,1000+i);
  //}

  gROOT->ProcessLine(".L calcUnfoldingSystematics.C+");
  calcUnfoldingSystematics("../config_files/xsecCalc.conf");


  //void plotDYUnfoldingMatrix(const TString input, int systematicsMode = 0, int randomSeed = 1, int reweightInt = 100, double massLimit = -1)
  //systematicsMode 0 - no systematic calc, no reweighting
  //1 - systematic mode, 2 - (reweighting of mass diff < -1 GeV) mode
  //check mass spectra with reweight = 95%; 100%; 105%  
  //mass value until which do reweighting
}
