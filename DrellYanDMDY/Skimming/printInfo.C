#include <TROOT.h>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include "../Include/InputFileMgr.hh"
#include <iostream>

void printInfo(TString fname) {
  TFile infile(fname);
  if (!infile.IsOpen()) {
    std::cout << "failed to open a file <" << fname << ">\n";
    return;
  }
  std::cout << "Keys : ";
  infile.ls();

  TTree *descrTree=(TTree*)infile.Get("Description");
  if (descrTree) {
    TDescriptiveInfo_t *info=new TDescriptiveInfo_t();
    descrTree->SetBranchAddress("description",&info);
    TBranch *br=descrTree->GetBranch("description");
    br->GetEntry(0);
    info->print();
    delete br;
    delete info;
  }
  delete descrTree;
  infile.Close();
}
