{
  //gROOT->ProcessLine(".L prepareYields.C+");
  //prepareYields("../config_files/dataT3.conf");
  gROOT->ProcessLine(".L subtractBackground.C+");
  subtractBackground("../config_files/dataT3.conf");
}
