//
// This file contains implementation of methods related to applying
// electron energy scale corrections.
//

#include "ElectronEnergyScaleAdv.hh"
#include <fstream>
#include <sstream>
#include "../../Include/MyTools.hh"

#include <TCanvas.h>
#include <TFile.h>
#include <TBranch.h>
#include <TTree.h>


//------------------------------------------------------

int ElectronEnergyScaleAdv_t::getEtaEtaIdx(double eta1, double eta2) const {
  int idx1=getEtaBinIdx(eta1)-1;
  int idx2=getEtaBinIdx(eta2)-1;
  int idx=-1;
  if ((idx1>=0) && (idx2>=0)) {
    if (idx1>idx2) { idx=idx1; idx1=idx2; idx2=idx; }
    idx=(2*_nEtaBins-1-idx1)*idx1/2 + idx2;
  }
  return idx;
}

//------------------------------------------------------

int ElectronEnergyScaleAdv_t::LoadEEMFile(const TString &eemFileName, vector<vector<double>*> &eemData) const {
  EtaEtaMassData_t *eem = new EtaEtaMassData_t();
  int res=1;
  int etaEtaCount = this->numberOfEtaEtaBins();
  eemData.clear(); eemData.reserve(etaEtaCount+1);

  // Read file twice. First time we will determine the number of masses
  // in each (eta1,eta2) bin and allocate memory. The second time we will
  // store the mass value
  const int optimize=1; // optimize memory
  for (int loop=0; res && (loop<1+optimize); ++loop) {
    //std::cout << "loop=" << loop << std::endl;
    std::vector<int> counts(etaEtaCount);
    const char *fname=eemFileName.Data();
    TFile *fin = new TFile(fname);
    assert(fin);
    TTree *tree = (TTree*)fin->Get("Data"); assert(tree);
    tree->SetBranchAddress("Data",&eem);
    TBranch *branch = tree->GetBranch("Data"); assert(branch);
    UInt_t entryMax = tree->GetEntries();
    for (UInt_t ientry=0; ientry < entryMax; ++ientry) {
      branch->GetEntry(ientry);
      int idx=this->getEtaEtaIdx(eem->eta1(),eem->eta2());
      //std::cout << loop << " got " << (*eem) << ", idx=" << idx << "\n";
      if ((idx>=0) && (idx < etaEtaCount)) {
	if (optimize && (loop==0)) {
	  // loops 0 -- count the number of events
	  counts[idx]++;
	}
	else {
	  // loops 1 -- store the memory
	  std::vector<double>* store=eemData[idx];
	  store->push_back(eem->mass());
	}
      }
    }
    //delete branch;
    //delete tree;
    delete fin;
    if (optimize && (loop==0)) {
      // allocate memory
      unsigned int k=0;
      for (int i=0; i<this->numberOfEtaBins(); ++i) {
	for (int j=i; j<this->numberOfEtaBins(); ++j, ++k) {
	  std::vector<double>* tmp=new std::vector<double>();
	  assert(tmp);
	  //std::cout << "counts[k]=" << counts[k] << "\n";
	  tmp->reserve(counts[k]);
	  eemData.push_back(tmp);
	}
      }
    }
  }
  return 1;
}

//------------------------------------------------------
