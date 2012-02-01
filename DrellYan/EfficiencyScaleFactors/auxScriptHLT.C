{  // auxScriptAdv allows to turn on/off different calculations

  TString triggerSet="Full2011_hltEffNew";
  TString mcMainInputFile="../config_files/fall11mc.input"; // used in CalcEventEff.C
  TString mcFileStart="../config_files/sf";   // beginning of the file used in eff_*.C 
  TString inpFile;
  //mcFileStart="../config_files/sfFall11";

  if ( ! triggerSet.Contains("hltEffNew") || ! triggerSet.Contains("Full2011") ) {
  gROOT->ProcessLine(".L eff_IdHlt.C+");
  inpFile = mcFileStart + TString("_mc_HLT.conf");
  eff_IdHlt(inpFile,triggerSet);
  
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
  std::cout<<"DONE: eff_IdHlt(\"" << inpFile << "\",\"" << triggerSet <<"\")"<<std::endl;
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
  std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;  
  }
  
  gROOT->ProcessLine(".L eff_IdHlt.C+");

  if ( ! triggerSet.Contains("hltEffNew") || ! triggerSet.Contains("Full2011") ) {
    inpFile = mcFileStart + TString("_data_HLT.conf");
    eff_IdHlt(inpFile,triggerSet);
  
    std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
    std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
    std::cout<<"DONE: eff_IdHlt(\"" << inpFile << "\",\"" << triggerSet <<"\")"<<std::endl;
    std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
    std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;  
  }
  else {
    TString mcInpFile = mcFileStart + TString("_mc_HLT.conf");
    inpFile = mcFileStart + TString("_data_HLT.conf");
    TString *trigSetArr= new TString[3];
    trigSetArr[0]="2011A_SingleEG_hltEffNew"; trigSetArr[1]="2011A_DoubleEG_hltEffNew"; trigSetArr[2]="2011B_DoubleEG_hltEffNew";
    for (int k=0; k<3; ++k) {
      eff_IdHlt(mcInpFile,trigSetArr[k]); // create templates
      std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
      std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
      std::cout<<"DONE: eff_IdHlt(\"" << mcInpFile << "\",\"" << trigSetArr[k] <<"\")"<<std::endl;
      std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
      std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;  

      eff_IdHlt(inpFile,trigSetArr[k]);
      std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
      std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
      std::cout<<"DONE: eff_IdHlt(\"" << inpFile << "\",\"" << trigSetArr[k] <<"\")"<<std::endl;
      std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;
      std::cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"<<std::endl;  
    }
  }

}
