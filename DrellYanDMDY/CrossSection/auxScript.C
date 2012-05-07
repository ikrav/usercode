{
  gROOT->ProcessLine(".L calcCrossSection.C+");
  calcCrossSection("../config_files/xsecCalc.conf");
}
