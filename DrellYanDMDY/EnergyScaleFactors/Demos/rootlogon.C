{  
  gSystem->Exec("cd ../; ln -s ../Include . ; ln -s ../Unfolding . ; cd Demos");
  gROOT->ProcessLine(".x ../Include/rootlogon.C");

}
