{  
  gSystem->Exec("cd ../; ln -s ../Include . ; ln -s ../Unfolding . ; cd Demos");
  gROOT->ProcessLine(".x ../Interface/rootlogon.C");

  //gROOT->ProcessLine(".L eesSmearEventDemo.C+");

}
