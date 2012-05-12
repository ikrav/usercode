#ifndef PUReweight_HH
#define PUReweight_HH

#include <TString.h>
#include <TFile.h>
#include <TH1F.h>
#include <iostream>

class PUReweight_t {
protected:
  TString FName; // file name
  TFile *FFile; // pointer to a file
  TH1F *hRef;   // reference histogram
  TH1F *hActive; // active histogram
  TH1F *hWeight; // histogram of weights
public:
  PUReweight_t() : FName(), FFile(NULL), hRef(NULL), hActive(NULL), hWeight(NULL) { }
  ~PUReweight_t() { this->clear(); }

  void clear() {
    if (hRef) delete hRef;
    if (hActive) delete hActive;
    if (hWeight) delete hWeight;
    if (FFile) delete FFile;
  }

  // access
  const TH1F* getHRef() const { return hRef; }
  const TH1F* getHActive() const { return hActive; }
  const TH1F* getHWeigth() const { return hWeight; }


  double getWeight(int nGoodPV) const {
    double w=0;
    if (!hWeight) {
      std::cout << "PUReweight::getWeight: call setActiveSample first\n";
    }
    else {
      int idx=hWeight->FindBin( nGoodPV );
      w=hWeight->GetBinContent( idx );
    }
    return w;
  }


  int setFile(const TString &fname);
  int setReference(const TString &setName);
  int setActiveSample(const TString &setName);

  int printWeights(std::ostream&) const;
  int printActiveDistr_and_Weights(std::ostream& out) const;
};


#endif
