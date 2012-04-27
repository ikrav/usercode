#include "FEWZ.hh"
#include <TFile.h>
#include <TString.h>

// ----------------------------------------------------------

FEWZ_t::FEWZ_t(bool loadWeights, bool do_cutZPT100) : fInitialized(kFALSE), fCutZPT100(do_cutZPT100) {
  if (loadWeights) {
    TFile fweights("../root_files/fewz/weights_stepwise_prec10-5_fine12.root");
    if(do_cutZPT100)
      cout << "FEWZ_t NOTE: in MC, for Z/gamma* PT>100 the FEWZ weights for 80<PT<100 GeV are used!" << endl;
    if( !fweights.IsOpen() ) assert(0);
    bool ok=kTRUE;
    for(int i=0; ok && (i<_nMassBinsFEWZ); i++){
      TString hnames = TString::Format("weight_%02d",i+1);
      weights[i] = (TH2D*)fweights.Get(hnames);
      hnames = TString::Format("h_weighterror_%02d",i+1);
      weightErrors[i] = (TH2D*)fweights.Get(hnames);
      if (!weights[i] || !weightErrors[i]) ok=kFALSE;
      else {
	weights[i]->SetDirectory(0);
	weightErrors[i]->SetDirectory(0);
      }
    }
    fInitialized=ok;
  }
}

// ----------------------------------------------------------
