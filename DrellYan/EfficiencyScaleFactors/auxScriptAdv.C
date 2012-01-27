{
  const int calcMCReco=1;
  const int calcMCID  =0;
  const int calcMCHLT =0;
  const int calcExpReco=1;
  const int calcExpID  =0;
  const int calcExpHLT =0;
  const int calcEffScaleFactors=0;
  TString triggerSet="Full2011";


  if (calcMCReco) {
  gROOT->ProcessLine(".L eff_Reco.C+");
  eff_Reco("../config_files/sfFall11_mc_RECO.conf",triggerSet);
  
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
  std::cout<<"DONE: eff_Reco(\"../config_files/sfFall11_mc_RECO.conf\"," << triggerSet << ")"<<std::endl;
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
  }
  
  if (calcExpReco) {
  gROOT->ProcessLine(".L eff_Reco.C+");
  eff_Reco("../config_files/sfFall11_data_RECO.conf",triggerSet);
  
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
  std::cout<<"DONE: eff_Reco(\"../config_files/sfFall11_data_RECO.conf\"," << triggerSet << ")"<<std::endl;
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;  
  }


  if (calcMCID) {
  gROOT->ProcessLine(".L eff_IdHlt.C+");
  eff_IdHlt("../config_files/sfFall11_mc_ID.conf",triggerSet);
  
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
  std::cout<<"DONE: eff_IdHlt(\"../config_files/sfFall11_mc_ID.conf\"," << triggerSet << ")"<<std::endl;
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;  
  }
   
  if (calcExpID) {
  gROOT->ProcessLine(".L eff_IdHlt.C+");
  eff_IdHlt("../config_files/sfFall11_data_ID.conf",triggerSet);
  
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
  std::cout<<"DONE: eff_IdHlt(\"../config_files/sfFall11_data_ID.conf\"," << triggerSet << ")"<<std::endl;
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;  
  }
   
  if (calcMCHLT) {
  gROOT->ProcessLine(".L eff_IdHlt.C+");
  eff_IdHlt("../config_files/sfFall11_mc_HLT.conf",triggerSet);
  
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
  std::cout<<"DONE: eff_IdHlt(\"../config_files/sfFall11_mc_HLT.conf\"," << triggerSet << ")"<<std::endl;
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;  
  }
  
  if (calcExpHLT) {
  gROOT->ProcessLine(".L eff_IdHlt.C+");
  eff_IdHlt("../config_files/sfFall11_data_HLT.conf",triggerSet);
  
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
  std::cout<<"DONE: eff_IdHlt(\"../config_files/sfFall11_data_HLT.conf\"," << triggerSet << ")"<<std::endl;
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;  
  }
  
  if (calcEffScaleFactors) {
  gROOT->ProcessLine(".L calcEventEff.C+");
  calcEventEff("../config_files/summer11mc.input",triggerSet);
  
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
  std::cout<<"DONE: calcEventEff(\"../config_files/summer11mc.input\"," << triggerSet << ")"<<std::endl;
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl; 
  }


}
