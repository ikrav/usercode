{
  gROOT->ProcessLine(".L subtractBackground.C+");
  subtractBackground("../config_files/data.conf");
}
