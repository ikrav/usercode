{  
  if (gSystem->Getenv("CMSSW_VERSION")) {
    TString rfitpath("/afs/cern.ch/cms/sw/slc5_ia32_gcc434/lcg/roofit/5.26.00-cms5/include");
    TString path = gSystem->GetIncludePath();
    path += "-I. -I$ROOTSYS/src -I";
    path += rfitpath;
    gSystem->SetIncludePath(path.Data());
    
    TString str = gSystem->GetMakeSharedLib();
    if (str.Contains("-m32")==0 && str.Contains("-m64")==0) {
      str.ReplaceAll("g++", "g++ -m32");
      gSystem->SetMakeSharedLib(str);
    }      
    
  }
  
// Load "MIT Style" plotting
  gROOT->Macro("../Include/CPlot.cc+");
  gROOT->Macro("../Include/MitStyleRemix.cc+");

  // Load structures for ntuple analysis
  // The ones below have class definitions and need to be compiled
  gROOT->ProcessLine(".L ../Include/TDielectron.hh+");
  gROOT->ProcessLine(".L ../Include/TElectron.hh+");
  gROOT->ProcessLine(".L ../Include/TMuon.hh+");
  gROOT->ProcessLine(".L ../Include/TEventInfo.hh+");
  gROOT->ProcessLine(".L ../Include/TGenInfo.hh+");
  gROOT->ProcessLine(".L ../Include/TPhoton.hh+");
  gROOT->ProcessLine(".L ../Include/TJet.hh+");
  gROOT->ProcessLine(".L ../Include/TVertex.hh+");
  gROOT->ProcessLine(".L ../Include/EleIDCuts.hh+");
  gROOT->ProcessLine(".L ../Include/JsonParser.cc+");
             
}
