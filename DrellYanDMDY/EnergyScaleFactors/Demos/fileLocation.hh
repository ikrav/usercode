#ifndef fileLocation_HH
#define fileLocation_HH

void setEEMFileLocation(TString &eemFileExp, TString &eemFileMC) {
  TString path="../../root_files/selected_events/DY_m10+pr+a05+o03+pr_4680pb/ntuples/";
  path="/media/spektras/DYee2011/2ndProd-ntuple/";
  eemFileExp= path; 
  eemFileExp.Append("data_bpVer1_EtaEtaM.root");
  eemFileMC=path; 
  eemFileMC.Append("zee_bpVer1_EtaEtaM.root");
}

#endif
