#include "../Include/PUReweight.hh"
#include "assert.h"

// --------------------------------------------------------------
PUReweight_t::PUReweight_t(): FName(), FFile(NULL), hRef(NULL), 
			      hActive(NULL), hWeight(NULL), FCreate(0) 
{ 

  initializeHildrethWeights();

}

// --------------------------------------------------------------

void PUReweight_t::initializeHildrethWeights(){

  // The pile-up reweighting setup according to the Hildreth's method
  // is done below.
  //   For general information on this, see
  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupMCReweightingUtilities
  //
  //   The "target" (data) PU distribution comes
  // from luminosity and inelastic pp cross section convolution,
  // and is obtained according to the prescription from 
  //     https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJSONFileforData
  // It is independent on analysis selection, but depends on the JSON file.
  // This distribution comes from running an official tool (takes ~1 hour
  // on full 2011 JSON), and it returns a ROOT file with a TH1D with a certain
  // binning. One needs to repack it then to the standard DY binning and TH1F
  // separately before using it in this class.
  // 
  //   The "source" (MC) PU distribution is the generator-level number 
  // of simulated pile-up. It is independent on the physical process
  // in the MC sample, or selection applied in the analysis. It depends
  // on the MC production, i.e. on the pile-up scenario (such as "S6" or S7"),
  // and on the vintage (Fall11, or Summer12).
  //     To prepare it with Drell-Yan ntuples, one needs to simply take
  // any MC sample and plot the "Info.nPU" variable which, for MC, contains
  // the number of simulated pile-up at generator level. The binning *has* to be
  // the same as that of the target distribution. Note: the comment
  // in the TEventInfo.hh header says it is "reconstructed vertices", but
  // it is a misinformation, according to the author Kevin Sung, it is really
  // the gen-level quantity.

  TString ftargetName = "../root_files/pileup/dataPileupHildreth_full2011_20121110_repacked.root";
  TFile f1(ftargetName);
  if( ! f1.IsOpen()){
    printf("Failed to find the target for Hildreth's PU reweighting\n");
    printf("  failed to open the file %s\n", ftargetName.Data());
    printf("  See PUReweight.cc for more detail\n");
    assert(0);
  }
  TH1F *target = (TH1F*)f1.Get("pileup_lumibased_data");
  if( target == 0 ){
    printf("Failed to find the histogram for pileup in the file %s\n", 
	   ftargetName.Data());
    printf("  See PUReweight.cc for more detail\n");
    assert(0);
  }

  TString fsourceName = "../root_files/pileup/mcPileupHildreth_full2011_20121110_repacked.root";
  TFile f2(fsourceName);
  if( ! f2.IsOpen()){
    printf("Failed to find the source for Hildreth's PU reweighting\n");
    printf("  failed to open the file %s\n", fsourceName.Data());
    printf("  See PUReweight.cc for more detail\n");
    assert(0);
  }
  TH1F *source = (TH1F*)f2.Get("pileup_simulevel_mc");
  if( source == 0 ){
    printf("Failed to find the histogram for pileup in the file %s\n", 
	   fsourceName.Data());
    printf("  See PUReweight.cc for more detail\n");
    assert(0);
  }

  // Make sure histograms have the same binning
  bool nBinsMatch = (source->GetNbinsX() == target->GetNbinsX());
  bool lowBoundaryMatch = 
    ( source->GetXaxis()->GetBinLowEdge(1) == target->GetXaxis()->GetBinLowEdge(1));
  bool upBoundaryMatch = 
    ( source->GetXaxis()->GetBinUpEdge(source->GetNbinsX())
      == target->GetXaxis()->GetBinUpEdge(target->GetNbinsX()));
  if( !( nBinsMatch && lowBoundaryMatch && upBoundaryMatch ) ){
    printf("Failed to find weights in PUReweight: the source and the target\n");
    printf(" for the Hildreth's method reweighting have different binning\n");
    assert(0);
  }
  
  // Normalize the source and the target
  target->Scale(1/target->GetSumOfWeights());
  source->Scale(1/source->GetSumOfWeights());

  // Find the weights distribution
  hWeightHildreth = (TH1F*)target->Clone("hWeightHildreth");
  hWeightHildreth->SetDirectory(0);
  hWeightHildreth->Divide(source);
  
  f1.Close();
  f2.Close();
  
//   printf("Pileup weights for the Hildreth method are constructed.\n");
//   for(int i=1; i<= hWeightHildreth->GetNbinsX(); i++){
//     printf("PU=%2d  weight=%f\n", i, hWeightHildreth->GetBinContent(i));
//   }
  return;
}

// --------------------------------------------------------------

int printHisto_local(std::ostream& out, const TH1F* histo) {
  if (!histo) {
    out << "printHisto: histo is null\n";
    return 0;
  }
  char buf[100];
  out << "values of " << histo->GetName() << "\n";
  for(int i=1; i<=histo->GetNbinsX(); i++) {
    double x=histo->GetBinLowEdge(i);
    double w=histo->GetBinWidth(i);
    sprintf(buf," %5.2f-%5.2f    %f    %f\n",
	    x,x+w,histo->GetBinContent(i),histo->GetBinError(i));
    out << buf;
  }
  return 1;
}

 
// --------------------------------------------------------------
// --------------------------------------------------------------

int PUReweight_t::setFile(const TString &fname, int create) {
  if ((fname == FName) && (FCreate==create)) {
    std::cout << "PUReweight::setFile warning: method repeatedly called for the same file name <" << fname << ">, create=" << create << "\n";
    return 0;
  }
  this->clear(); // save active histo if needed and close file
  TString opt;
  switch(create) {
  case 0: opt="READ"; break;
  case 1: opt="RECREATE"; break;
  case 2: opt="UPDATE"; break;
  default:
    std::cout << "PUReweight::setFile: create=" << create << " option is not recognized\n";
    return 0;
  }
  FCreate=create;
  FFile = new TFile(fname,opt.Data());
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
  if (FCreate) {
    std::cout << "File is opened for writing (create=" << FCreate << "), it's not possible to setReference\n";
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

  if (FCreate==0) {
    // reading mode

    if (!hRef) {
      std::cout << "PUReweight::setActiveSample(" << name << "): hRef is not set (call setReference first)\n";
      return 0;
    }

    if (hActive) { delete hActive; hActive=0; }
    if (hWeight) { delete hWeight; hWeight=0; }
    
    hActive = (TH1F*) FFile->Get(name);
    if (!hActive) {
      std::cout << "PUReweight::setActiveSample(" << name << "): failed to set active sample\n";
      return 0;
    }
    if (!prepareWeights(0)) {
      std::cout << "in method setActiveSample\n";
      return 0;
    }
  }
  else {
    // writing mode
    if (hActive) {
      FFile->cd(); hActive->Write();
      delete hActive;
    }
    hActive=this->newHisto(name);
    hActive->SetDirectory(FFile);
  }
  return 1;
}

// --------------------------------------------------------------

int PUReweight_t::prepareWeights(int save_weight) {
  if (!hActive) {
    std::cout << "PUReweight::prepareWeights: hActive is not set (maybe call setActiveSample)\n";
    return 0;
  }
  if (!hRef) {
    std::cout << "PUReweight::prepareWeights: hRef is not set (call setReference first)\n";
    return 0;
  }
  
  if (hWeight) delete hWeight;
  
  // check that the PU division is the same
  int ok=(hActive->GetNbinsX() == hRef->GetNbinsX()) ? 1:0;
  for (int ibin=1; ok && (ibin<=hActive->GetNbinsX()); ++ibin) {
    if ( (hActive->GetBinLowEdge(ibin) != hRef->GetBinLowEdge(ibin)) ||
	 (hActive->GetBinWidth(ibin) != hRef->GetBinWidth(ibin)) ) {
      ok=0;
    }
  }
  if (!ok) {
    std::cout << "prepareWeights: hRef and hActive have different binnings\n";
    std::cout << "hRef: "; printHisto_local(std::cout, hRef);
    std::cout << "hActive: "; printHisto_local(std::cout, hActive);
    return 0;
  }
  

  hWeight= (TH1F*)hRef->Clone( hActive->GetName() + TString("_puWeights") );
  assert(hWeight);
  hWeight->SetDirectory(0);
  hWeight->Scale( hActive->GetSumOfWeights() / hRef->GetSumOfWeights() );
  hWeight->Divide(hActive);
  if ((hWeight->GetBinLowEdge(1)==-0.5) && (hWeight->GetBinWidth(1)==1.)) {
    hWeight->SetBinContent(1,0.); hWeight->SetBinError(1,0.);
  }
  if (save_weight) {
    if (FCreate!=0) { FFile->cd(); hWeight->Write(); }
    else std::cout << "PUReweight::prepareWeights: save is requested but the file is in reading mode\n";
  }
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

int PUReweight_t::printHisto(std::ostream& out, const TH1F* histo, const TString &name) const {
  if (!histo) {
    out << "PUReweight::printHisto: " << name << " is not set\n";
    return 0;
  }
  char buf[100];
  out << "values of " << histo->GetName() << "\n";
  for(int i=1; i<=histo->GetNbinsX(); i++) {
    sprintf(buf," %5.2f    %f    %f\n",
	    histo->GetBinCenter(i),histo->GetBinContent(i),histo->GetBinError(i));
    out << buf;
  }
  return 1;
}

// --------------------------------------------------------------
// --------------------------------------------------------------
