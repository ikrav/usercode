#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TCanvas.h>                // class for drawing
#include <TH1F.h>                   // 1D histograms
#include <TGraphErrors.h>           // graphs
#include <TProfile.h>               // profile histograms
#include <TBenchmark.h>             // class to track macro running statistics
#include <TLorentzVector.h>         // class for Lorentz vector computations
#include <TVectorD.h>
#include <TGraphErrors.h>
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <vector>                   // STL vector class
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings

#include "../Include/CPlot.hh"          // helper class for plots
#include "../Include/MitStyleRemix.hh"  // style settings for drawing
#include "../Include/MyTools.hh"        // miscellaneous helper functions
#include "../Include/DYTools.hh"

// define classes and constants to read in ntuple
#include "../Include/EWKAnaDefs.hh"
#include "../Include/TGenInfo.hh"
#endif

//=== FUNCTION DECLARATIONS ======================================================================================

//=== MAIN MACRO =================================================================================================

void plotDYAcceptance(const TString input, int systematicsMode = DYTools::NORMAL, double reweightFsr = 1.0, double massLimit=-1)
//systematicsMode 0 (NORMAL) - no systematic calc
//2 (FSR_STUDY) - systematics due to FSR, reweighting
//check mass spectra with reweight = 95%; 100%; 105%  
//mass value until which do reweighting 
{
  // check whether it is a calculation
  if (input.Contains("_DebugRun_")) {
    std::cout << "plotDYAcceptance: _DebugRun_ detected. Terminating the script\n";
    return;
  }

  // normal calculation

  gBenchmark->Start("plotDYAcceptance");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================
  
  if (systematicsMode==DYTools::NORMAL)
    std::cout<<"Running script in the NORMAL mode"<<std::endl;
  else if (systematicsMode==DYTools::FSR_STUDY)
    std::cout<<"Running script in the FSR_STUDY mode"<<std::endl;
  else { 
    std::cout<<"requested mode not recognized"<<std::endl;
    assert(0);
  }

  Bool_t doSave  = false;    // save plots?
  TString format = "png";   // output file format
  
  vector<TString> fnamev;   // file names   
  vector<TString> labelv;   // legend label
  vector<Int_t>   colorv;   // color in plots
  vector<Int_t>   linev;    // line style
  vector<Double_t> xsecv;
  vector<Double_t> lumiv;
  TString          dirTag;

  Double_t massLow  = DYTools::massBinLimits[0];
  Double_t massHigh = DYTools::massBinLimits[nMassBins];
  
  ifstream ifs;
  ifs.open(input.Data());
  assert(ifs.is_open());
  string line;
  Int_t state=0;
  while(getline(ifs,line)) {
    if(line[0]=='#') continue;
    if(state == 0){
      dirTag = TString(line);
      state++;
      continue;
    }else{
      string fname;
      Int_t color, linesty;
      stringstream ss(line);
      Double_t xsec;
      ss >> fname >> xsec >> color >> linesty;
      string label = line.substr(line.find('@')+1);
      fnamev.push_back(fname);
      labelv.push_back(label);
      colorv.push_back(color);
      linev.push_back(linesty);
      xsecv.push_back(xsec);
      lumiv.push_back(0);
    }
  }
  ifs.close();
  
  const Double_t kGAP_LOW  = 1.4442;
  const Double_t kGAP_HIGH = 1.566;


  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================
  
  //  
  // Set up histograms
  //
  vector<TH1F*> hZMassv;//, hZMass2v, hZPtv, hZPt2v, hZyv, hZPhiv;  
  
  Double_t   nZv = 0;
  TVectorD nEventsv (DYTools::nMassBins);  
  TVectorD nPassv   (DYTools::nMassBins);
  TVectorD accv     (DYTools::nMassBins);
  TVectorD accErrv  (DYTools::nMassBins);

  TVectorD nPassBBv(DYTools::nMassBins), accBBv(DYTools::nMassBins), accErrBBv(DYTools::nMassBins); 
  TVectorD nPassBEv(DYTools::nMassBins), accBEv(DYTools::nMassBins), accErrBEv(DYTools::nMassBins); 
  TVectorD nPassEEv(DYTools::nMassBins), accEEv(DYTools::nMassBins), accErrEEv(DYTools::nMassBins);

  // Vectors for calculation of errors with weighted sums
  TVectorD w2Eventsv (DYTools::nMassBins);  
  TVectorD w2Passv   (DYTools::nMassBins);


  nEventsv = 0;
  nPassv   = 0;
  nPassBBv = 0;
  nPassBEv = 0;
  nPassEEv = 0;

  w2Eventsv = 0;
  w2Passv   = 0;


  char hname[100];
  for(UInt_t ifile = 0; ifile<fnamev.size(); ifile++) {
    sprintf(hname,"hZMass_%i",ifile); hZMassv.push_back(new TH1F(hname,"",500,0,500)); hZMassv[ifile]->Sumw2();
  }

  // 
  // Read weights from a file
  //
  const bool useFewzWeights = true;
  const bool cutZPT100 = true;
  if(cutZPT100)
    cout << "NOTE: in MC, for Z/gamma* PT>100 the FEWZ weights for 80<PT<100 GeV are used!" << endl;
  TH2D *weights[DYTools::nMassBins];
  TH2D *weightErrors[DYTools::nMassBins];
  TFile fweights("../root_files/fewz/weights_stepwise_prec10-5_fine12.root");
  if( !fweights.IsOpen() ) assert(0);
  for(int i=0; i<DYTools::nMassBins; i++){
    TString hname_loc = TString::Format("weight_%02d",i+1);
    weights[i] = (TH2D*)fweights.Get(hname_loc);
    hname_loc = TString::Format("h_weighterror_%02d",i+1);
    weightErrors[i] = (TH2D*)fweights.Get(hname_loc);
  }

  //
  // Access samples and fill histograms
  //  
  TFile *infile=0;
  TTree *eventTree=0;  
    
  // Data structures to store info from TTrees
  mithep::TGenInfo *gen  = new mithep::TGenInfo();
  
  // loop over samples  
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
  
    // Read input file
    cout << "Processing " << fnamev[ifile] << "..." << endl;
    infile = new TFile(fnamev[ifile]); 
    assert(infile);

    // Get the TTrees
    eventTree = (TTree*)infile->Get("Events"); assert(eventTree);

    // Find weight for events for this file
    // The first file in the list comes with weight 1,
    // all subsequent ones are normalized to xsection and luminosity
    double xsec=xsecv[ifile];
    AdjustXSectionForSkim(infile,xsec,eventTree->GetEntries(),1);
    lumiv[ifile] = eventTree->GetEntries()/xsec;
    double scale = lumiv[0]/lumiv[ifile];
    //     if(ifile != 0) scale *= 0.87;
    cout << "       -> sample weight is " << scale << endl;

    // Set branch address to structures that will store the info  
    eventTree->SetBranchAddress("Gen",&gen);
    TBranch *genBr = eventTree->GetBranch("Gen");
  
    // loop over events    
    nZv += scale * eventTree->GetEntries();
    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
      //if (ientry>100000) break;

      genBr->GetEntry(ientry);

      // Which mass is used?
      double massPreFsr = gen->vmass;   // pre-FSR
      double mass = gen->mass;    // post-FSR
      if((mass < massLow) || (mass > massHigh)) continue;

      double reweight;
      if (systematicsMode!=DYTools::FSR_STUDY) reweight=1.0;
      else if ((mass-massPreFsr)>massLimit) reweight=1.0;
      else reweight=reweightFsr;

      
      int ibin = DYTools::findMassBin(mass);
      int ibinPreFsr = DYTools::findMassBin(massPreFsr);
      // If mass is larger than the highest bin boundary
      // (last bin), use the last bin.
      if(ibinPreFsr == -1 && mass >= massBinLimits[nMassBins] )
	ibinPreFsr = nMassBins-1;
      // Find FEWZ-powheg reweighting factor 
      // that depends on pre-FSR Z/gamma* rapidity, pt, and mass
      double fewz_weight = 1.0;
      if(useFewzWeights){
	if(ibinPreFsr != -1 && ibinPreFsr < DYTools::nMassBins){
	  int ptBin = weights[ibinPreFsr]->GetXaxis()->FindBin( gen->vpt );
	  int yBin = weights[ibinPreFsr]->GetYaxis()->FindBin( gen->vy );
	  // In case if pt or y are outside of the weight maps,
	  // set them to the closest bin.
 	  if(ptBin == weights[ibinPreFsr]->GetNbinsX() + 1)
 	    ptBin = weights[ibinPreFsr]->GetNbinsX();
 	  if(ptBin == 0)
 	    ptBin = 1;
 	  if(yBin == weights[ibinPreFsr]->GetNbinsY() + 1)
 	    yBin = weights[ibinPreFsr]->GetNbinsY();
 	  if(yBin == 0)
 	    yBin = 1;
	  // Apply PT cut if needed
	  if( cutZPT100 ) 
	    if( ptBin == weights[ibinPreFsr]->GetNbinsX() )
	      ptBin = weights[ibinPreFsr]->GetNbinsX() - 1;
	  fewz_weight = weights[ibinPreFsr]->GetBinContent( ptBin, yBin);
						      
	}else
	  cout << "Error: binning problem" << endl;
      }
//       printf("mass= %f   pt= %f    Y= %f     weight= %f\n",gen->mass, gen->vpt, gen->vy, fewz_weight);

      if(ibin != -1 && ibin < nEventsv.GetNoElements()){
	double fullWeight = reweight * scale * gen->weight * fewz_weight;
	nEventsv[ibin] += fullWeight;
	w2Eventsv[ibin] += fullWeight*fullWeight;
      }else if(ibin >= nEventsv.GetNoElements())
	cout << "ERROR: binning problem" << endl;

      Bool_t isB1 = (fabs(gen->eta_1)<kGAP_LOW);
      Bool_t isB2 = (fabs(gen->eta_2)<kGAP_LOW);
      Bool_t isE1 = (fabs(gen->eta_1)>kGAP_HIGH);
      Bool_t isE2 = (fabs(gen->eta_2)>kGAP_HIGH);
      // Asymmetric Et cut scheme for DY analysis
      if( ( (gen->pt_1 > DYTools::etMinLead && gen->pt_2 > DYTools::etMinTrail) 
	    || (gen->pt_1 > DYTools::etMinTrail && gen->pt_2 > DYTools::etMinLead) )
         && ((fabs(gen->eta_1)<kGAP_LOW) || (fabs(gen->eta_1)>kGAP_HIGH))
         && ((fabs(gen->eta_2)<kGAP_LOW) || (fabs(gen->eta_2)>kGAP_HIGH))   
	 && (fabs(gen->eta_1)<2.5) && (fabs(gen->eta_2)<2.5)) {
        
	if(ibin != -1 && ibin < nPassv.GetNoElements()){
	  double fullWeight = reweight * scale * gen->weight * fewz_weight;
	  nPassv[ibin] += fullWeight;
	  w2Passv[ibin] += fullWeight*fullWeight;
	  if(isB1 && isB2)                          { nPassBBv[ibin] += reweight * scale * gen->weight * fewz_weight; } 
	  else if(isE1 && isE2)                     { nPassEEv[ibin] += reweight * scale * gen->weight * fewz_weight; } 
	  else if((isB1 && isE2) || (isE1 && isB2)) { nPassBEv[ibin] += reweight * scale * gen->weight * fewz_weight; }
	}
      }
      hZMassv[ifile]->Fill(mass,reweight * scale * gen->weight * fewz_weight);
    }   
    delete infile;
    infile=0, eventTree=0;
  }
  delete gen;
  
  accv      = 0;
  accErrv   = 0;
  accBBv    = 0;
  accErrBBv = 0;
  accBEv    = 0;
  accErrBEv = 0;
  accEEv    = 0;
  accErrEEv = 0;
  for(int i=0; i<DYTools::nMassBins; i++){
    if(nEventsv[i] != 0){
      accv[i] = nPassv[i]/nEventsv[i];
      // The commented out piece does not have correct error calculation
      // for weighted events.
      //
      // accErrv[i] = sqrt(accv[i]*(1-accv[i])/nEventsv[i]); 
      //
      // The correct error calculation when Npass and Nfail
      // are accumulated as sum of weights:
      //    Npass = sum(weight) over all passed events
      //    NpassErr = sqrt( sum( weight^2 ) ) over all passed events
      //    Similarly, for Nfail.
      //    After that, we propagate errors in formula Npass/(Npass+Nfail)
      //    where the errors are computed as above.
      //
      double nTotal = nEventsv[i];
      double nPass = nPassv[i];
      double nPassErr = sqrt( w2Passv[i] );
      double nFail = nTotal - nPass;
      double nFailErr = sqrt( w2Eventsv[i] - w2Passv[i] );
      accErrv[i] = sqrt( ( nFail*nFail * nPassErr*nPassErr + nPass*nPass * nFailErr*nFailErr)
			 / (nTotal*nTotal*nTotal*nTotal));

      accBBv[i] = nPassBBv[i]/nEventsv[i];
      accErrBBv[i] = sqrt(accBBv[i]*(1-accBBv[i])/nEventsv[i]);
      
      accBEv[i] = nPassBEv[i]/nEventsv[i];
      accErrBEv[i] = sqrt(accBEv[i]*(1-accBEv[i])/nEventsv[i]);
      
      accEEv[i] = nPassEEv[i]/nEventsv[i];
      accErrEEv[i] = sqrt(accEEv[i]*(1-accEEv[i])/nEventsv[i]);
    }
  }
  
  //--------------------------------------------------------------------------------------------------------------
  // Make plots 
  //==============================================================================================================  
  TCanvas *c = MakeCanvas("c","c",800,600);

  // string buffers
  char ylabel[50];   // y-axis label

  // Z mass
  sprintf(ylabel,"a.u. / %.1f GeV/c^{2}",hZMassv[0]->GetBinWidth(1));
  CPlot plotZMass1("zmass1","","m(Z) [GeV/c^{2}]",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    plotZMass1.AddHist1D(hZMassv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotZMass1.SetLogy();
  plotZMass1.Draw(c,doSave,format);

  // Prepare the overall acceptance plot
  double x [DYTools::nMassBins];
  double dx[DYTools::nMassBins];
  double y [DYTools::nMassBins];
  double dy[DYTools::nMassBins];
  for(int i=0; i<DYTools::nMassBins; i++){
    x[i]  = (DYTools::massBinLimits[i  ] + DYTools::massBinLimits[i+1])/2.0;
    dx[i] = (DYTools::massBinLimits[i+1] - DYTools::massBinLimits[i  ])/2.0;
    y[i] = accv[i];
    dy[i] = accErrv[i];
  }
  TGraphErrors *acceptanceGraph = 
    new TGraphErrors(DYTools::nMassBins,x,y,dx,dy);

  TCanvas *c1 = MakeCanvas("c1","c1",1000,600);
  CPlot plotAcceptance("Acceptance","","M_{ee} [GeV/c^{2}]","acceptance");
  plotAcceptance.SetLogx();
  plotAcceptance.AddGraph((TGraph*)acceptanceGraph,"PE",600,kFullDotMedium,1); 
  plotAcceptance.SetYRange(0,1);
  plotAcceptance.SetXRange(10,1500.0);
  acceptanceGraph->GetYaxis()->SetTitleOffset(1.0);
  acceptanceGraph->GetXaxis()->SetMoreLogLabels();
  acceptanceGraph->GetXaxis()->SetNoExponent();
  acceptanceGraph->SetMarkerStyle(20);
  acceptanceGraph->SetMarkerSize(1);
  plotAcceptance.Draw(c1,doSave,format);
          
  // Store constants in the file
  TString accConstFileName(TString("../root_files/"));
  if (systematicsMode==DYTools::NORMAL){
    accConstFileName+=TString("constants/")+dirTag;
    gSystem->mkdir(accConstFileName,kTRUE);
    accConstFileName+=TString("/acceptance_constants.root");
  }
  else if (systematicsMode==DYTools::FSR_STUDY){
    accConstFileName+=TString("systematics/")+dirTag;
    gSystem->mkdir(accConstFileName,kTRUE);
    accConstFileName+=TString("/acceptance_constants_reweight_");
    accConstFileName+=int(reweightFsr*100);
    accConstFileName+=TString(".root");
  }
  TFile fa(accConstFileName,"recreate");
  acceptanceGraph->Write("A_FSR");
  accv.Write("acceptanceArray");
  accErrv.Write("acceptanceErrArray");
  fa.Close();

  //--------------------------------------------------------------------------------------------------------------
  // Summary print out
  //==============================================================================================================
  cout << endl;
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  cout << endl; 
  
  cout << labelv[0] << " file: " << fnamev[0] << endl;
  cout << "     Number of generated events: " << nZv << endl;
  printf(" mass bin    preselected      passed     total_A_GEN      BB-BB_A_GEN      EB-BB_A_GEN      EB-EB_A_GEN\n");
  for(int i=0; i<DYTools::nMassBins; i++){
    printf(" %4.0f-%4.0f   %10.0f   %10.0f   %7.4f+-%6.4f  %7.4f+-%6.4f  %7.4f+-%6.4f  %7.4f+-%6.4f \n",
	   DYTools::massBinLimits[i], DYTools::massBinLimits[i+1],
	   nEventsv[i], nPassv[i],
	   accv[i], accErrv[i],
	   accBBv[i], accErrBBv[i],
	   accBEv[i], accErrBEv[i],
	   accEEv[i], accErrEEv[i]);
  }

  cout << endl;
  
  gBenchmark->Show("plotDYAcceptance");
}


//=== FUNCTION DEFINITIONS ======================================================================================

//--------------------------------------------------------------------------------------------------
