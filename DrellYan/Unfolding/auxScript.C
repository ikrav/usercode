{
  gROOT->ProcessLine(".L plotDYUnfoldingMatrix.C+");
  plotDYUnfoldingMatrix("../config_files/summer11mc.input");
  for (int i=1; i<=20; i++)
  {
    plotDYUnfoldingMatrix("../config_files/summer11mc.input",true,1000+i);
  }
  gROOT->ProcessLine(".L calcUnfoldingSystematics.C+");
  calcUnfoldingSystematics("../config_files/xsecCalc.conf")
}
