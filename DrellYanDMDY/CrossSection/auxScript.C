{
  gROOT->ProcessLine(".L calcCrossSectionFsr.C+");
  calcCrossSectionFsr("../config_files/xsecCalc.conf");
}
