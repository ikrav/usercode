{  
 
  gROOT->ProcessLine(".L ../../Include/TDielectron.hh+");
  gROOT->ProcessLine(".L ../../Include/TElectron.hh+");
  gROOT->ProcessLine(".L ../../Include/DYTools.hh+");
  gROOT->ProcessLine(".L ../../Include/TMuon.hh+");
  gROOT->ProcessLine(".L ../../Include/TEventInfo.hh+");
  gROOT->ProcessLine(".L ../../Include/TGenInfo.hh+");
  gROOT->ProcessLine(".L ../../Include/TPhoton.hh+");
  gROOT->ProcessLine(".L ../../Include/TJet.hh+");
  gROOT->ProcessLine(".L ../../Include/TVertex.hh+");
  gROOT->ProcessLine(".L ../../Include/EleIDCuts.hh+");
  gROOT->ProcessLine(".L ../../Include/ZeeData.hh+");

  gROOT->ProcessLine(".L ../../Include/TriggerSelection.hh+");

  gROOT->ProcessLine(".L ../../Include/JsonParser.cc+");
  gROOT->ProcessLine(".L ../../Include/ElectronEnergyScale.cc+");
  gROOT->ProcessLine(".L ../../Include/EtaEtaMass.hh+");
  gROOT->ProcessLine(".L ../../Include/FEWZ.cc+");
  gROOT->ProcessLine(".L ../../Include/EventSelector.cc+");
  gROOT->ProcessLine(".L ../../Include/InputFileMgr.cc+");
  gROOT->ProcessLine(".L ../../Include/PUReweight.cc+");

  //gROOT->ProcessLine(".L ../../Unfolding/UnfoldingTools.C+");
  //gROOT->ProcessLine(".L ../../Include/plotFunctions.cc+");

   //gROOT->ProcessLine(".x ../../Include/rootlogon.C");
  gROOT->Macro("TreeQueue.cc+");
  gROOT->Macro("ConfList.cc+");
}
