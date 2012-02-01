{  // auxScriptAdv allows to turn on/off different calculations

  TString triggerSet="Full2011_hltEffNew";
  TString mcMainInputFile="../config_files/fall11mc.input"; // used in CalcEventEff.C
  TString mcFileStart="../config_files/sf";   // beginning of the file used in eff_*.C 
  TString inpFile;
  //mcFileStart="../config_files/sfFall11";

  gROOT->ProcessLine(".L eff_Reco.C+");
  inpFile = mcFileStart + TString("_mc_RECO.conf");
  eff_Reco(inpFile,triggerSet);
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
  std::cout<<"DONE: eff_Reco(\""<< inpFile << "\",\"" << triggerSet <<"\")"<<std::endl;
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;


  gROOT->ProcessLine(".L eff_Reco.C+");
  inpFile = mcFileStart + TString("_data_RECO.conf");
  eff_Reco(inpFile,triggerSet);
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
  std::cout<<"DONE: eff_Reco(\""<< inpFile << "\",\"" << triggerSet <<"\")"<<std::endl;
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;  
}
