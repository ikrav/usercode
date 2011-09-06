{
  gROOT->ProcessLine(".L plotDYEfficiency.C+");
  plotDYEfficiency("../config_files/summer11mc.input");
}
