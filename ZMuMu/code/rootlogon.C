{
  gROOT->ProcessLine(".x CPlot.cc+");
  gROOT->ProcessLine(".x MitStyleRemix.cc+");
       
  gROOT->ProcessLine(".L Selection.cc+");
  gROOT->ProcessLine(".L SampleBase.cc+");
  gROOT->ProcessLine(".L SampleMIT.cc+");

}
