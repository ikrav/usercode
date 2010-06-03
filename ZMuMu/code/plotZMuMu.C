//================================================================================================
//
// Plots for reconstructed Z candidates
//
// Z candidate categories:
//   0. standard selection: 2 isolated global muons at least one matched to an HLT primitive
//   1. mu-mu-2HLT  : 2 isolated global muons both matched to HLT primitives
//   2. mu-mu-1HLT  : 2 isolated global muons only one compatible with an HLT primitive
//   3. mu-tk       : 1 isolated global muon + 1 isolated tracker muon
//   4. mu-sa       : 1 isolated global muon + 1 isolated standalone muon
//   5. mu-mu-noIso : 2 global muons, at least one compatible with HLT, at least one FAILS isolation
//
//________________________________________________________________________________________________

#include "SampleMIT.h"
#include "Selection.h"

#include <TBenchmark.h>
#include <TH1F.h>

#include "CPlot.hh"              // helper class for plots
#include "MitStyleRemix.hh"      // style settings for drawing

#include <iostream>
#include <iomanip>
using namespace std;

void plotZMuMu(Int_t category=0){

  gBenchmark->Start("plotZMuMu");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================

  double lumi = 10.0; // luminosity normalization [pb^-1]

  //---------------------------------------------------------
  // Define all samples
  //---------------------------------------------------------
  vector<SampleMIT*> samples;
  
  // Data. Leave blank if we do not have any.
  SampleMIT *sampleData = new SampleMIT("");
  samples.push_back(sampleData);
  
  TString mcLocation = "/server/03a/schorlem/ntuples3/";
  cout << "Create sample objects" << endl;
  // Signal Monte Carlo
  SampleMIT *sampleMCZmm = new SampleMIT(mcLocation+TString("s09-zmm-7-mc3_ntuple.root"),
			1606.6,                 // cross-section
			"Z #rightarrow #mu#mu", // label
			kBlue);
  samples.push_back(sampleMCZmm);


  // Assorted backgrounds
  SampleMIT *sampleMCwm = new SampleMIT(mcLocation+TString("s09-wm-7-mc3_ntuple.root"),
		     9679.9*0.742,                 // cross-section
		     "W #rightarrow #mu#nu", // label
		     kOrange+7);
  samples.push_back(sampleMCwm);

  SampleMIT *sampleMCttbar = new SampleMIT(mcLocation+TString("s09-ttbar-7-mc3_ntuple.root"),
		     165.0,                 // cross-section
		     "t#bar{t}", // label
		     kGreen+2);
  samples.push_back(sampleMCttbar);

  SampleMIT *sampleMCZtt = new SampleMIT(mcLocation+TString("s09-ztt-7-mc3_ntuple.root"),
		     1606.6,                 // cross-section
		     "Z #rightarrow #tau#tau", // label
		     kMagenta+2);
  samples.push_back(sampleMCZtt);

  SampleMIT *sampleMCincmu15 = new SampleMIT(mcLocation+TString("s09-incmu15-7-mc3_ntuple.root"),
		     0.2969*0.00037*1e+09,                 // cross-section
		     "QCD", // label
		     kRed);
  samples.push_back(sampleMCincmu15);

  SampleMIT *sampleMCww = new SampleMIT(mcLocation+TString("s09-ww-7-mc3_ntuple.root"),
		     42.9,                 // cross-section
		     "WW", // label
		     kYellow+1);
  samples.push_back(sampleMCww);

  // Ntuple initialization and some constants
  cout << "Initialize ntuples" << endl;
  for(int isample = 0; isample < samples.size(); isample++){
    if (isample == 0 && !samples[isample]->isInitialized() )continue;
    samples[isample]->setup();
    if(samples[isample]->getNEvents() != 0)
      samples[isample]->setWeight(lumi*samples[isample]->getCrossSection()/
 				 (double)samples[isample]->getNEvents());
    else
      samples[isample]->setWeight(0.0);
  }

  //---------------------------------------------------------
  // Set up selection
  //---------------------------------------------------------

  int    trigger  = tHLT_Mu9;
  double minMass  = 60;
  double maxMass  = 120;
  double minPt    = 20;
  double maxEta1  = 2.4;
  double maxEta2  = 2.1;
  double trackIso = 3;
  Selection selection;
  selection.initialize(trigger,minPt,maxEta1,maxEta2,minMass,maxMass,trackIso);

  //---------------------------------------------------------
  // Loop over samples and candidates
  //---------------------------------------------------------

  //
  // Set up histograms<D-c>
  //  
  vector <TH1F*> hZMassv;
  vector <TH1F*> hNJetsv;
  TH1F *hbkg = new TH1F("hbkg", "all bkgs", 180, 20, 200);
 
  for(UInt_t isample=0; isample<samples.size(); isample++) {
    char hname[100];    
    sprintf(hname,"hZMass_%i",isample); hZMassv.push_back(new TH1F(hname,"",150,50,200)); hZMassv[isample]->Sumw2();
    sprintf(hname,"hNJets_%i",isample); hNJetsv.push_back(new TH1F(hname,"",20,-0.5,19.5)); hNJetsv[isample]->Sumw2();
  }


  // flag for if there is a data sample to process
  Bool_t hasData = samples[0]->isInitialized();
  
  // loop over samples
  for(UInt_t isample=0; isample<samples.size(); isample++) {  
    SampleMIT *s = samples[isample];
    cout << endl << "Process sample " << s->getLabel() << endl;
    cout << "total events "     << s->getNEvents() << endl;
    cout << "total candidates " << s->getNCandidates() << endl;
    for( int ientry=0; ientry<s->getNCandidates(); ientry++){
      if(ientry%10000 == 0) {
	printf("."); fflush(stdout);
      }
      //
      // Apply canidate selection
      //
      // General cuts
      if( ! selection.passBaselineSelection( s->triggerBits(ientry),
					     s->mass(ientry),
					     s->mu1pt(ientry),
					     s->mu2pt(ientry),
					     s->mu1eta(ientry),
					     s->mu2eta(ientry)) )
	  continue;
      // Apply category cuts
      int foundCategory = selection.findCategory( s->mu1typeBits(ientry),
						  s->mu2typeBits(ientry),
						  s->mu1trkIso(ientry),
						  s->mu2trkIso(ientry),
						  s->mu1hltMatchBits(ientry),
						  s->mu2hltMatchBits(ientry),
						  s->mu1charge(ientry),
						  s->mu2charge(ientry) );
      if( category == 0 && !(foundCategory == 1 || foundCategory == 2) ) 
	continue;
      else if ( category != 0 && ! (category == foundCategory) )
	continue;	   

      // Tests and additional cuts if needed:
      double pt1 = s->mu1pt(ientry);
      double pt2 = s->mu2pt(ientry);
      int nJets = s->nJets(ientry);
      int mu1iso = s->mu1trkIso(ientry);
      int mu2iso = s->mu2trkIso(ientry);
      int mu1type = s->mu1typeBits(ientry);
      int mu2type = s->mu2typeBits(ientry);
      //if( ! (pt1 > 30 && pt2 > 30) ) continue;
      //if( (mu1iso > 3 && mu2iso > 3) ) continue;
      //if(nJets>0) continue;
      //if( !(mu1type == mNoMuon || mu2type == mNoMuon) )continue;

//       // Require opposite charge (if not using categories, otherwise
//       // the cut is applied in the category requirements)
//       if( s->mu1charge(ientry)*s->mu2charge(ientry) >0 )
// 	continue;

      // Fill histograms
      double mass = s->mass(ientry);
      hZMassv[isample]->Fill(s->mass(ientry),s->getWeight());
      hNJetsv[isample]->Fill(s->nJets(ientry),s->getWeight());
      if(isample>1) // excluding 0=data and 1=signalMC
	hbkg->Fill(s->mass(ientry),s->getWeight());

      // If we do not have data, then in the "data" histogram add all MC contributions
      if( !hasData)
	hZMassv[0]->Fill(s->mass(ientry),s->getWeight());
    
    } // end loop over entries

  } // end loop over samples

  //---------------------------------------------------------
  // Display mass plots
  //---------------------------------------------------------
  TCanvas *c = MakeCanvas("c","c",800,600);
  
  // string buffers
  char pname[50];    // plot name
  char ylabel[50];   // y-axis label
  char text[100];    // text box string
  
  // plot title
  TString title = "";

  // set how many decimal places to print depending on "lumi" variable
  Int_t precision = 2;
  if(lumi >= 1) precision = 1;
  if(lumi >= 10) precision = 0;
  
  // dimuon mass
  sprintf(pname,"zmass%i",category);
  sprintf(ylabel,"Events / %.1f GeV/c^{2}",hZMassv[0]->GetBinWidth(1));
  CPlot plotZMass(pname,title,"m(#mu#mu) [GeV/c^{2}]",ylabel);
  plotZMass.AddHist1D(hZMassv[0],
		      samples[0]->getLabel(),"E",
		      samples[0]->getColor());
  for(UInt_t ibg=2; ibg<hZMassv.size(); ibg++)
  {
    plotZMass.AddToStack(hZMassv[ibg],
			 samples[ibg]->getLabel(),
			 samples[ibg]->getColor());
  }
  plotZMass.AddToStack(hZMassv[1],
		       samples[1]->getLabel(),
		       samples[1]->getColor() );
  sprintf(text,"#int#font[12]{L}dt = %.*f pb^{-1}",precision,lumi);
  plotZMass.AddTextBox(text,0.21,0.85,0.41,0.8,0);
  plotZMass.SetYRange(0.01,1000);
  plotZMass.SetLogy();
  plotZMass.Draw(c,false);    
 

    //--------------------------------------------------------------------------------------------------------------
  // Summary print out
  //==============================================================================================================
  cout << endl;
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  cout << "  Luminosity = " << setprecision(2) << lumi << " / pb" << endl;
  cout << endl;
 
  if(hasData) 
    cout << "       Data files: " << samples[0]->getLabel() << endl;
  
  if(category==0)      { cout << " *** Category: mu-mu-1.or.2HLT" << endl; }
  else if(category==1) { cout << " *** Category: mu-mu-2HLT" << endl; }
  else if(category==2) { cout << " *** Category: mu-mu-1HLT" << endl; }
  else if(category==3) { cout << " *** Category: mu-tk" << endl; }
  else if(category==4) { cout << " *** Category: mu-sa" << endl; }
  else if(category==5) { cout << " *** Category: mu-mu-noIso" << endl; }
  else if(category==6) { cout << " *** Category: mu-mu-same-sign" << endl; }
  
  cout << setw(20) << "Total events =";
  cout << setw(8) << setprecision(2) << fixed << hZMassv[0]->Integral() << endl;

  if(samples.size()>1) {  
    cout << setw(20) << "MC Signal =";
    cout << setw(8) << setprecision(2) << fixed << hZMassv[1]->Integral();
    cout << " \u00B1 "; 
    cout << setw(4) << setprecision(2) << sqrt(hZMassv[1]->GetEffectiveEntries())*samples[1]->getWeight();
    cout << ",  xsec = " << setw(9) << setprecision(2) << fixed << samples[1]->getCrossSection() << " pb";
    cout << ",  files: " << samples[1]->getLabel() << endl; 
  
    for(UInt_t ibg=2; ibg<hZMassv.size(); ibg++) {    
      cout << setw(17) << "MC Bkgd " << ibg-1 << " =";
      cout << setw(8) << setprecision(2) << fixed << hZMassv[ibg]->Integral();
      cout << " \u00B1 ";
      cout << setw(4) << setprecision(2) << sqrt(hZMassv[ibg]->GetEffectiveEntries())*samples[ibg]->getWeight();
      cout << ",  xsec = " << setw(9) << setprecision(2) << fixed << samples[ibg]->getCrossSection() << " pb"; 
      cout << ",  files: " << samples[ibg]->getLabel() << endl;     
    }
    cout << endl;  
  }


  gBenchmark->Show("plotZMuMu");
  return;
}
