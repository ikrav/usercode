{
  gROOT->ProcessLine(".L selectEvents.C+");
  selectEvents("../config_files/dataT3.conf","Full2011DatasetTriggers",DYTools::NORMAL); 
  //selectEvents("../config_files/dataReduced.conf"); 
}
