{  
  
  gROOT->ProcessLine(".x ../Include/rootlogon.C");

  gROOT->ProcessLine(".L cutFunctions.cc+");
  gROOT->ProcessLine(".L fitFunctionsCore.cc+");
  gROOT->ProcessLine(".L fitFunctions.cc+");

}

