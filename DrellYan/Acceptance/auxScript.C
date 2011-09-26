{
  gROOT->ProcessLine(".L plotDYAcceptanceNew.C+");
  plotDYAcceptance("../config_files/summer11mc.input");
  plotDYAcceptance("../config_files/summer11mc.input",2,1.05);
  plotDYAcceptance("../config_files/summer11mc.input",2,0.95);
}
