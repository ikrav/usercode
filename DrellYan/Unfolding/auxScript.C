{
  gROOT->ProcessLine(".L plotDYUnfoldingMatrixNew.C+");
  //plotDYUnfoldingMatrix("../config_files/summer11mc.input",2,1,100,-1);
  plotDYUnfoldingMatrix("../config_files/summer11mc.input",2,1,1.05,-1);
  plotDYUnfoldingMatrix("../config_files/summer11mc.input",2,1,0.95,-1);

  gROOT->ProcessLine(".L calcUnfoldingSystematicsNew.C+");
  calcUnfoldingSystematics("../config_files/xsecCalc.conf");


  //void plotDYUnfoldingMatrix(const TString input, int systematicsMode = 0, int randomSeed = 1, int reweightInt = 100, double massLimit = -1)
  //systematicsMode 0 - no systematic calc, no reweighting
  //1 - systematic mode, 2 - (reweighting of mass diff < -1 GeV) mode
  //check mass spectra with reweight = 95%; 100%; 105%  
  //mass value until which do reweighting
}
