{
  gROOT->ProcessLine(".L plotDYAcceptance.C+");


  plotDYAcceptance("../config_files/summer11mc.input");
  plotDYAcceptance("../config_files/summer11mc.input",2,1.05);
  plotDYAcceptance("../config_files/summer11mc.input",2,0.95);


  gROOT->ProcessLine(".L calcAcceptanceSystematics.C+");
  calcAcceptanceSystematics("../config_files/xsecCalc.conf");

  

}
