{
  const int calcMCReco=0;
  const int calcMCID  =1;
  const int calcMCHLT =0;
  const int calcExpReco=0;
  const int calcExpID  =0;
  const int calcExpHLT =0;
  const int calcEffScaleFactors=0;
  TString triggerSet="Full2011";

  TString mcMainInputFile="../config_files/fall11mc.input"; // used in CalcEventEff.C
  TString mcFileStart="../config_files/sf";   // beginning of the file used in eff_*.C 
  TString inpFile;
  mcFileStart="../config_files/sfFall11";


  if (calcMCReco) {
  gROOT->ProcessLine(".L eff_Reco.C+");
  inpFile = mcFileStart + TString("_mc_RECO.conf");
  eff_Reco(inpFile,triggerSet);
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
  std::cout<<"DONE: eff_Reco(\""<< inpFile << "\",\"" << triggerSet <<"\")"<<std::endl;
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
  }
  
  if (calcExpReco) {
  gROOT->ProcessLine(".L eff_Reco.C+");
  inpFile = mcFileStart + TString("_data_RECO.conf");
  eff_Reco(inpFile,triggerSet);
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
  std::cout<<"DONE: eff_Reco(\""<< inpFile << "\",\"" << triggerSet <<"\")"<<std::endl;
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;  
  }


  if (calcMCID) {
  gROOT->ProcessLine(".L eff_IdHlt.C+");
  inpFile = mcFileStart + TString("_mc_ID.conf");
  eff_IdHlt(inpFile,triggerSet);
  
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
  std::cout<<"DONE: eff_IdHlt(\"" << inpFile << "\",\"" << triggerSet <<"\")"<<std::endl;
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;  
  }
   
  if (calcExpID) {
  gROOT->ProcessLine(".L eff_IdHlt.C+");
  inpFile = mcFileStart + TString("_data_ID.conf");
  eff_IdHlt(inpFile,triggerSet);
  
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
  std::cout<<"DONE: eff_IdHlt(\"" << inpFile << "\",\"" << triggerSet <<"\")"<<std::endl;
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;  
  }
   
  if (calcMCHLT) {
  gROOT->ProcessLine(".L eff_IdHlt.C+");
  inpFile = mcFileStart + TString("_mc_HLT.conf");
  eff_IdHlt(inpFile,triggerSet);
  
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
  std::cout<<"DONE: eff_IdHlt(\"" << inpFile << "\",\"" << triggerSet <<"\")"<<std::endl;
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;  
  }
  
  if (calcExpHLT) {
  gROOT->ProcessLine(".L eff_IdHlt.C+");
  inpFile = mcFileStart + TString("_data_HLT.conf");
  eff_IdHlt(inpFile,triggerSet);
  
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
  std::cout<<"DONE: eff_IdHlt(\"" << inpFile << "\",\"" << triggerSet <<"\")"<<std::endl;
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;  
  }
  
  if (calcEffScaleFactors) {
  gROOT->ProcessLine(".L calcEventEff.C+");
  calcEventEff(mcMainInputFile,triggerSet);

  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
  std::cout<<"DONE: calcEventEff(\"" << mcMainInputFile << "\",\"" << triggerSet << "\")"<<std::endl;
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl; 
  }


}
