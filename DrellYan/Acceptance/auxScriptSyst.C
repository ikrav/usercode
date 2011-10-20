{
  gROOT->ProcessLine(".L plotDYAcceptance.C+");


  plotDYAcceptance("../config_files/summer11mc.input",DYTools::FSR_STUDY,1.05);
  plotDYAcceptance("../config_files/summer11mc.input",DYTools::FSR_STUDY,0.95);


  gROOT->ProcessLine(".L calcAcceptanceSystematics.C+");
  calcAcceptanceSystematics("../config_files/xsecCalc.conf");

  

}
