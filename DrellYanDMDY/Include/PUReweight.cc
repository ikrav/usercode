#include "../Include/PUReweight.hh"

// --------------------------------------------------------------

int PUReweight_t::setFile(const TString &fname) {
  if (fname == FName) {
    std::cout << "PUReweight::setFile warning: method repeatedly called for the same file name <" << fname << ">\n";
    return 0;
  }
  this->clear();
  FFile = new TFile(fname);
  if (!FFile || !FFile->IsOpen()) {
    std::cout << "failed to open a file <" << fname << ">\n";
    return 0;
  }
  FName=fname;
  return 1;
}


// --------------------------------------------------------------

int PUReweight_t::setReference(const TString &name) {
  if (!FFile) {
    std::cout << "PUReweight::setReference(" << name << "): file is not open\n";
    return 0;
  }
  if (hRef) delete hRef;
  hRef = (TH1F*) FFile->Get(name);
  if (!hRef) {
    std::cout << "PUReweight::setReference(" << name << "): failed to set reference\n";
    return 0;
  }
  hRef->SetDirectory(0);
  return 1;
}

// --------------------------------------------------------------

int PUReweight_t::setActiveSample(const TString &name) {
  if (!FFile) {
    std::cout << "PUReweight::setActiveSample(" << name << "): file is not open\n";
    return 0;
  }
  if (!hRef) {
    std::cout << "PUReweight::setActiveSample(" << name << "): hRef is not set (call setReference first)\n";
    return 0;
  }
  if (hActive) delete hActive;
  if (hWeight) delete hWeight;

  hActive = (TH1F*) FFile->Get(name);
  if (!hActive) {
    std::cout << "PUReweight::setActiveSample(" << name << "): failed to set active sample\n";
    return 0;
  }
  hActive->SetDirectory(0);
  hWeight= (TH1F*)hRef->Clone(name+TString("_puWeights"));
  hWeight->SetDirectory(0);
  hWeight->Scale( hActive->GetSumOfWeights() / hRef->GetSumOfWeights() );
  hWeight->Divide(hActive);
  return 1;
}

// --------------------------------------------------------------

int PUReweight_t::printActiveDistr_and_Weights(std::ostream& out) const {
  if (!hActive || !hWeight) {
    out << "PUReweight::printActiveDistr_and_Weights: hActive or hWeight is not set (call setActiveSample first)\n";
    return 0;
  }
  char buf[100];
  out << "values of " << hActive->GetName() << " and " << hWeight->GetName() << "\n";
  for(int i=1; i<=hWeight->GetNbinsX(); i++) {
    sprintf(buf," %5.2f    %6.4f  %6.4f     %6.4f  %6.4f\n",hWeight->GetBinCenter(i),hActive->GetBinContent(i),hActive->GetBinError(i),hWeight->GetBinContent(i),hWeight->GetBinError(i));
    out << buf;
  }
  return 1;
}

// --------------------------------------------------------------

int PUReweight_t::printWeights(std::ostream& out) const {
  if (!hWeight) {
    out << "PUReweight::printWeights: hWeight is not set (call setActiveSample first)\n";
    return 0;
  }
  char buf[100];
  out << "values of " << hWeight->GetName() << "\n";
  for(int i=1; i<=hWeight->GetNbinsX(); i++) {
    sprintf(buf," %5.2f    %f    %f\n",hWeight->GetBinCenter(i),hWeight->GetBinContent(i),hWeight->GetBinError(i));
    out << buf;
  }
  return 1;
}

// --------------------------------------------------------------
// --------------------------------------------------------------
