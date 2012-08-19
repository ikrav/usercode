#ifndef PUReweight_HH
#define PUReweight_HH

#include <TROOT.h>
#include <TString.h>
#include <TFile.h>
#include <TH1F.h>
#include <iostream>
//#include "../Include/DYTools.hh"
//#include "../Include/MyTools.hh"


class PUReweight_t {
public:
  typedef enum { maxPVs=45 } TConst_t;
protected:
  TString FName; // file name
  TFile *FFile; // pointer to a file
  TH1F *hRef;   // reference histogram
  TH1F *hActive; // active histogram
  TH1F *hWeight; // histogram of weights
  int FCreate; // whether a file is being created (1) or updated (2), otherwise - reading (0)
public:
  PUReweight_t() : FName(), FFile(NULL), hRef(NULL), hActive(NULL), hWeight(NULL), FCreate(0) { }
  ~PUReweight_t() { this->clear(); }

  void clear() {
    if (hRef) { delete hRef; hRef=0; }
    if ((FCreate!=0) && FFile && hActive) { FFile->cd(); hActive->Write(); }
    if (hActive) { delete hActive; hActive=0; }
    if (hWeight) { delete hWeight; hWeight=0; }
    if (FFile) { delete FFile; FFile=0; }
    FCreate=0;
  }

  // access
  const TString& fileName() const { return FName; }
  const TH1F* getHRef() const { return hRef; }
  const TH1F* getHActive() const { return hActive; }
  const TH1F* getHWeigth() const { return hWeight; }
  int getCreate() const { return FCreate; }


  double getWeight(int nGoodPV) const {
    double w=0;
    if (!hWeight) {
      std::cout << "PUReweight::getWeight: call setActiveSample first\n";
    }
    else {
      if (nGoodPV> PUReweight_t::maxPVs+1) nGoodPV=PUReweight_t::maxPVs+1;
      int idx=hWeight->FindBin( nGoodPV );
      w=hWeight->GetBinContent( idx );
    }
    return w;
  }

  int setDefaultFile(const TString &dirTag, const TString &analysisTag, 
		     int create=0) {
    TString fname=TString("../root_files/selected_events/") + dirTag + 
      TString("/npv") + analysisTag + TString(".root");
    int res=setFile(fname,create);
    if (!res) std::cout << "Error in PUReweight::setDefaultFile\n";
    return res;
  }

  int setFile(const TString &fname, int create=0);
  int setReference(const TString &setName);
  int setReference(const TH1F *new_hRef) { 
    hRef=(TH1F*)new_hRef->Clone(new_hRef->GetName()+TString("_clone")); 
    hRef->SetDirectory(0);
    return 1; 
  }
  int setActiveSample(const TString &setName); // set active sample and calculate weights
  int prepareWeights(int save_weights); // needs to be called if setReference(TH1F) was called
  
  int Fill(UInt_t nGoodPV, double weight) {
    if ((FCreate==0) || !hActive) {
      std::cout << "cannot fill the histogram:\n";
      if (FCreate==0) std::cout << " - file is not opened for creation\n";
      if (!hActive) std::cout << " - active histogram is not set\n";
      return 0;
    }
    if (nGoodPV> PUReweight_t::maxPVs+1) nGoodPV=PUReweight_t::maxPVs+1;
    hActive->Fill(nGoodPV,weight);
    return 1;
  }

  int printActiveDistr_and_Weights(std::ostream& out) const;

  int printWeights(std::ostream &out) const { 
    int res=this->printHisto(out,hWeight,"hWeight");
    if (!res) out << "in printWeights\n";
    return res;
  }

  int printActive(std::ostream &out) const { 
    int res=this->printHisto(out,hActive,"hActive");
    if (!res) out << "in printActive\n";
    return res;
  }

protected:
  TH1F *newHisto(const TString &name) const {
    // histogram: PUs from 0 to maxPVs with an overflow bin (maxPVs+1) 
    TH1F *h= new TH1F(name,name,PUReweight_t::maxPVs+2,-0.5,Double_t(PUReweight_t::maxPVs+1.5));
    h->Sumw2();
    h->GetXaxis()->SetTitle("nGoodPVs"); h->GetYaxis()->SetTitle("weight (a.u.)");
    return h;
  }

  int printHisto(std::ostream& out, const TH1F* histo, const TString &name) const;
};


#endif
