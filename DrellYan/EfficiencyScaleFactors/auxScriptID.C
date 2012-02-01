{  // auxScriptAdv allows to turn on/off different calculations

  TString triggerSet="Full2011_hltEffNew";
  TString mcMainInputFile="../config_files/fall11mc.input"; // used in CalcEventEff.C
  TString mcFileStart="../config_files/sf";   // beginning of the file used in eff_*.C 
  TString inpFile;
  //mcFileStart="../config_files/sfFall11";

  gROOT->ProcessLine(".L eff_IdHlt.C+");
  inpFile = mcFileStart + TString("_mc_ID.conf");
  eff_IdHlt(inpFile,triggerSet);
  
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
  std::cout<<"DONE: eff_IdHlt(\"" << inpFile << "\",\"" << triggerSet <<"\")"<<std::endl;
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;  

   
  gROOT->ProcessLine(".L eff_IdHlt.C+");
  inpFile = mcFileStart + TString("_data_ID.conf");
  eff_IdHlt(inpFile,triggerSet);
  
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
  std::cout<<"DONE: eff_IdHlt(\"" << inpFile << "\",\"" << triggerSet <<"\")"<<std::endl;
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;  

}
