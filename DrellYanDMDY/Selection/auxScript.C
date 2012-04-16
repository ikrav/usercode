{
  gROOT->ProcessLine(".L selectEvents.C+");
  selectEvents("../config_files/dataT3.conf"); 
  //selectEvents("../config_files/dataReduced.conf"); 
}
