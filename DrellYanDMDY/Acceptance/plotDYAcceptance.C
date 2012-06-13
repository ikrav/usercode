#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TStyle.h>
#include <TCanvas.h>                // class for drawing
#include <TH1F.h>                   // 1D histograms
#include <TGraphErrors.h>           // graphs
#include <TGraph2DErrors.h>           // graphs
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
#include "../Include/plotFunctions.hh"

// define classes and constants to read in ntuple
#include "../Include/EWKAnaDefs.hh"
#include "../Include/TGenInfo.hh"

#include "../Include/EventSelector.hh"
#include "../Include/FEWZ.hh"
#include "../Include/UnfoldingTools.hh"
#endif

//=== FUNCTION DECLARATIONS ======================================================================================

//=== MAIN MACRO =================================================================================================

void plotDYAcceptance(const TString input, int systematicsMode = DYTools::NORMAL, double reweightFsr = 1.0, double massLimit=-1, int debugMode=0)
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

  if (debugMode) std::cout << "\n\n\tDEBUG MODE is ON\n\n";

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

  //Bool_t doSave  = false;    // save plots?
  TString format = "png";   // output file format
  
  vector<TString> fnamev;   // file names   
  vector<TString> labelv;   // legend label
  vector<Int_t>   colorv;   // color in plots
  vector<Int_t>   linev;    // line style
  vector<Double_t> xsecv;
  vector<Double_t> lumiv;
  TString          dirTag;

  Double_t massLow  = DYTools::massBinLimits[0];
  Double_t massHigh = DYTools::massBinLimits[DYTools::nMassBins];
  
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
  
  //const Double_t kGAP_LOW  = 1.4442; - shifted to DYTools.hh
  //const Double_t kGAP_HIGH = 1.566;  - shifted to DYTools.hh


  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================
  
  //  
  // Set up histograms
  //

  //int nYBinsMax=DYTools::findMaxYBins();


  vector<TH1F*> hZMassv;//, hZMass2v, hZPtv, hZPt2v, hZyv, hZPhiv;  
  
  Double_t   nZv = 0;

  TMatrixD nEventsv (DYTools::nMassBins,nYBinsMax);  
  TMatrixD nPassv   (DYTools::nMassBins,nYBinsMax);
  TMatrixD accv     (DYTools::nMassBins,nYBinsMax);
  TMatrixD accErrv  (DYTools::nMassBins,nYBinsMax);

  TMatrixD nPassBBv(DYTools::nMassBins,nYBinsMax), accBBv(DYTools::nMassBins,nYBinsMax), accErrBBv(DYTools::nMassBins,nYBinsMax); 
  TMatrixD nPassBEv(DYTools::nMassBins,nYBinsMax), accBEv(DYTools::nMassBins,nYBinsMax), accErrBEv(DYTools::nMassBins,nYBinsMax); 
  TMatrixD nPassEEv(DYTools::nMassBins,nYBinsMax), accEEv(DYTools::nMassBins,nYBinsMax), accErrEEv(DYTools::nMassBins,nYBinsMax);

  // Vectors for calculation of errors with weighted sums
  TMatrixD w2Eventsv (DYTools::nMassBins,nYBinsMax);  
  TMatrixD w2Passv   (DYTools::nMassBins,nYBinsMax);


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
  FEWZ_t fewz(useFewzWeights,cutZPT100);
  if (useFewzWeights && !fewz.isInitialized()) {
    std::cout << "failed to prepare FEWZ correction\n";
    throw 2;
  }
  const int new_fewz_code=1;
  TH2D *weights[DYTools::nMassBins]; // temporary
  TH2D *weightErrors[DYTools::nMassBins]; // temporary
  if (!new_fewz_code) {
  if(cutZPT100)
    cout << "NOTE: in MC, for Z/gamma* PT>100 the FEWZ weights for 80<PT<100 GeV are used!" << endl;
  TFile fweights("../root_files/fewz/weights_stepwise_prec10-5_fine12.root");
  if( !fweights.IsOpen() ) assert(0);
  for(int i=0; i<DYTools::nMassBins; i++){
    TString hname_loc = TString::Format("weight_%02d",i+1);
    weights[i] = (TH2D*)fweights.Get(hname_loc);
    hname_loc = TString::Format("h_weighterror_%02d",i+1);
    weightErrors[i] = (TH2D*)fweights.Get(hname_loc);
    weights[i]->SetDirectory(0); weightErrors[i]->SetDirectory(0);
  }
  }
  //
  // Access samples and fill histograms
  //  
  TFile *infile=0;
  TTree *eventTree=0;  
    
  // Data structures to store info from TTrees
  mithep::TGenInfo *gen  = new mithep::TGenInfo();

  int noFewz=0;//counter of events for which fewz weight was not found
 
  // loop over samples  
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
  
    // Read input file
    //std::cout<<"fnamev["<<ifile<<"]="<<fnamev[ifile]<<std::endl;
    cout << "Processing " << fnamev[ifile] << "..." << endl;
    infile = new TFile(fnamev[ifile]); 
    assert(infile);

    // Get the TTrees
    eventTree = (TTree*)infile->Get("Events"); assert(eventTree);

    // Find weight for events for this file
    // The first file in the list comes with weight 1,
    // all subsequent ones are normalized to xsection and luminosity
    AdjustXSectionForSkim(infile,xsecv[ifile],eventTree->GetEntries(),1);
    lumiv[ifile] = eventTree->GetEntries()/xsecv[ifile];
    double scale = lumiv[0]/lumiv[ifile];
    //     if(ifile != 0) scale *= 0.87;
    cout << "       -> sample weight is " << scale << endl;

    // Set branch address to structures that will store the info  
    eventTree->SetBranchAddress("Gen",&gen);
    TBranch *genBr = eventTree->GetBranch("Gen");
 
    // loop over events    
    nZv += scale * eventTree->GetEntries();

    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
      if (debugMode && (ientry>10000)) break;

      genBr->GetEntry(ientry);

      // Which mass is used?
      double massPreFsr = gen->vmass;   // pre-FSR
      double mass = gen->mass;    // post-FSR
      //double yPreFsr = gen->vy;   // pre-FSR
      double y = gen->y;    // post-FSR

      if ((mass < massLow) || (mass > massHigh)) continue;
      if ((fabs(y) < DYTools::yRangeMin) || (fabs(y) > DYTools::yRangeMax)) continue;

      double reweight;
      if (systematicsMode!=DYTools::FSR_STUDY) reweight=1.0;
      else if ((mass-massPreFsr)>massLimit) reweight=1.0;
      else reweight=reweightFsr;

      //ibinM1DPreFsr is used to determine fewz_weights
      //when new_fewz_code=0
      int ibinM1DPreFsr = DYTools::_findMassBin2011(massPreFsr);
      // If mass is larger than the highest bin boundary
      // (last bin), use the last bin.
      if(ibinM1DPreFsr == -1 && 
	 massPreFsr >= DYTools::_massBinLimits2011[DYTools::_nMassBins2011] ) {
	ibinM1DPreFsr = DYTools::_nMassBins2011-1;
      }

      int ibinMass = DYTools::findMassBin(mass);
      //int ibinMassPreFsr = DYTools::findMassBin(massPreFsr);
      int ibinY = DYTools::findAbsYBin(ibinMass, y);
      //int ibinYPreFsr = DYTools::findAbsYBin(ibinMassPreFsr, yPreFsr);

      // Find FEWZ-powheg reweighting factor 
      // that depends on pre-FSR Z/gamma* rapidity, pt, and mass
      double fewz_weight = 1.0;

      if(useFewzWeights) {
	if (new_fewz_code) fewz_weight=fewz.getWeight(gen->vmass,gen->vpt,gen->vy);
	else
        {

	  if(ibinM1DPreFsr != -1 && ibinM1DPreFsr < DYTools::_nMassBins2011)
            {
	      int ptBin = weights[ibinM1DPreFsr]->GetXaxis()->FindBin( gen->vpt );
	      int yBin = weights[ibinM1DPreFsr]->GetYaxis()->FindBin( gen->vy );
	      // In case if pt or y are outside of the weight maps,
	      // set them to the closest bin.
 	      if(ptBin == weights[ibinM1DPreFsr]->GetNbinsX() + 1)
 	        ptBin = weights[ibinM1DPreFsr]->GetNbinsX();
 	      if(ptBin == 0)
 	        ptBin = 1;
 	      if(yBin == weights[ibinM1DPreFsr]->GetNbinsY() + 1)
 	        yBin = weights[ibinM1DPreFsr]->GetNbinsY();
 	      if(yBin == 0)
 	        yBin = 1;
	      // Apply PT cut if needed
	      if( cutZPT100 ) 
	        if( ptBin == weights[ibinM1DPreFsr]->GetNbinsX() )
	          ptBin = weights[ibinM1DPreFsr]->GetNbinsX() - 1;
	      fewz_weight = weights[ibinM1DPreFsr]->GetBinContent( ptBin, yBin);

            }
          else
	    //cout << "Error: binning problem" << endl;
            noFewz++;

	}
//       printf("mass= %f   pt= %f    Y= %f     weight= %f\n",gen->mass, gen->vpt, gen->vy, fewz_weight);
      }

      if(ibinMass != -1 && ibinMass < nEventsv.GetNrows()){
	double fullWeight = reweight * scale * gen->weight * fewz_weight;
	nEventsv(ibinMass,ibinY) += fullWeight;
	w2Eventsv(ibinMass,ibinY) += fullWeight*fullWeight;
      }else if(ibinMass >= nEventsv.GetNrows())
	cout << "ERROR: binning problem" << endl;

      Bool_t isB1 = (fabs(gen->eta_1)<kECAL_GAP_LOW);
      Bool_t isB2 = (fabs(gen->eta_2)<kECAL_GAP_LOW);
      Bool_t isE1 = (fabs(gen->eta_1)>kECAL_GAP_HIGH);
      Bool_t isE2 = (fabs(gen->eta_2)>kECAL_GAP_HIGH);
      // Asymmetric Et cut scheme for DY analysis
      if( ( (gen->pt_1 > DYTools::etMinLead && gen->pt_2 > DYTools::etMinTrail) 
	    || (gen->pt_1 > DYTools::etMinTrail && gen->pt_2 > DYTools::etMinLead) )
         && ((fabs(gen->eta_1)<kECAL_GAP_LOW) || (fabs(gen->eta_1)>kECAL_GAP_HIGH))
         && ((fabs(gen->eta_2)<kECAL_GAP_LOW) || (fabs(gen->eta_2)>kECAL_GAP_HIGH))   
	 && (fabs(gen->eta_1)<2.5) && (fabs(gen->eta_2)<2.5)) {
        
	if(ibinMass != -1 && ibinMass < nPassv.GetNrows()){
	  double fullWeight = reweight * scale * gen->weight * fewz_weight;
	  nPassv(ibinMass,ibinY) += fullWeight;
	  w2Passv(ibinMass,ibinY) += fullWeight*fullWeight;
	  if(isB1 && isB2)                          { nPassBBv(ibinMass,ibinY) += reweight * scale * gen->weight * fewz_weight; } 
	  else if(isE1 && isE2)                     { nPassEEv(ibinMass,ibinY) += reweight * scale * gen->weight * fewz_weight; } 
	  else if((isB1 && isE2) || (isE1 && isB2)) { nPassBEv(ibinMass,ibinY) += reweight * scale * gen->weight * fewz_weight; }
	}
      }
      hZMassv[ifile]->Fill(mass,reweight * scale * gen->weight * fewz_weight);

    }   
    delete infile;
    infile=0, eventTree=0;

  }
  delete gen;

  std::cout<<"for "<<noFewz<<" events fewz-weight was not found"<<std::endl;  

  accv      = 0;
  accErrv   = 0;
  accBBv    = 0;
  accErrBBv = 0;
  accBEv    = 0;
  accErrBEv = 0;
  accEEv    = 0;
  accErrEEv = 0;

  for(int i=0; i<DYTools::nMassBins; i++)
    for (int j=0; j<nYBins[i]; j++)
      {
        if(nEventsv(i,j) != 0)
          {
             accv(i,j) = nPassv(i,j)/nEventsv(i,j);
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
             double nTotal = nEventsv(i,j);
             double nPass = nPassv(i,j);
             double nPassErr = sqrt( w2Passv(i,j) );
             double nFail = nTotal - nPass;
             double nFailErr = sqrt( w2Eventsv(i,j) - w2Passv(i,j) );
             accErrv(i,j) = sqrt( ( nFail*nFail * nPassErr*nPassErr + nPass*nPass * nFailErr*nFailErr)
			 / (nTotal*nTotal*nTotal*nTotal));

             accBBv(i,j) = nPassBBv(i,j)/nEventsv(i,j);
             accErrBBv(i,j) = sqrt(accBBv(i,j)*(1-accBBv(i,j))/nEventsv(i,j));
      
             accBEv(i,j) = nPassBEv(i,j)/nEventsv(i,j);
             accErrBEv(i,j) = sqrt(accBEv(i,j)*(1-accBEv(i,j))/nEventsv(i,j));
      
             accEEv(i,j) = nPassEEv(i,j)/nEventsv(i,j);
             accErrEEv(i,j) = sqrt(accEEv(i,j)*(1-accEEv(i,j))/nEventsv(i,j));
          }
      }

  
  //--------------------------------------------------------------------------------------------------------------
  // Make plots 
  //==============================================================================================================  

  CPlot::sOutDir="plots" + analysisTag;
  if (systematicsMode==DYTools::FSR_STUDY) {
    CPlot::sOutDir += "_reweight";
    CPlot::sOutDir += int(reweightFsr*100);
  }

  TString outDir= TString("../root_files/");
  if (systematicsMode==DYTools::NORMAL)  outDir+=TString("constants/");
  else if (systematicsMode==DYTools::FSR_STUDY)  outDir+=TString("systematics/");
  outDir+=dirTag;
  gSystem->mkdir(outDir,kTRUE);

  TString fileNamePlots=outDir + "/acceptance_plots" + analysisTag + ".root";
  TFile *filePlots=new TFile(fileNamePlots,"recreate");
  if (!filePlots) {
    std::cout << "failed to create file <" << fileNamePlots << ">\n";
    throw 2;
  }
  TCanvas *c = MakeCanvas("canvAcceptance","canvAcceptance",800,600);

  // string buffers
  char ylabel[50];   // y-axis label

  // Z mass
  sprintf(ylabel,"a.u. / %.1f GeV/c^{2}",hZMassv[0]->GetBinWidth(1));
  CPlot plotZMass1("zmass1","","m(Z) [GeV/c^{2}]",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    plotZMass1.AddHist1D(hZMassv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotZMass1.SetLogy();
  plotZMass1.Draw(c);
  SaveCanvas(c, "zmass1");

  PlotMatrixVariousBinning(accv, "acceptance", "LEGO2",filePlots);
  filePlots->Close();
  //delete filePlots;
  

  //          ------ can be adjusted to use outDir
  // Store constants in the file
  TString accConstFileName(TString("../root_files/"));
  if (systematicsMode==DYTools::NORMAL){
    accConstFileName+=TString("constants/")+dirTag;
    gSystem->mkdir(accConstFileName,kTRUE);
    //accConstFileName+=TString("/acceptance_constants.root");
    accConstFileName+=TString("/acceptance_constants" + analysisTag + ".root");
  }
  else if (systematicsMode==DYTools::FSR_STUDY){
    accConstFileName+=TString("systematics/")+dirTag;
    gSystem->mkdir(accConstFileName,kTRUE);
    accConstFileName+=TString("/acceptance_constants_reweight_") + analysisTag;
    accConstFileName+=TString("_");
    accConstFileName+=int(reweightFsr*100);
    accConstFileName+=TString(".root");
  }
  TFile fa(accConstFileName,"recreate");
  accv.Write("acceptanceMatrix");
  accErrv.Write("acceptanceErrMatrix");
  unfolding::writeBinningArrays(fa);
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

  if (DYTools::study2D==0)
    {
      printf(" mass bin    preselected      passed     total_A_GEN      BB-BB_A_GEN      EB-BB_A_GEN      EB-EB_A_GEN\n");
      for(int i=0; i<DYTools::nMassBins; i++){
        printf(" %4.0f-%4.0f   %10.1f   %10.1f   %7.4f+-%6.4f  %7.4f+-%6.4f  %7.4f+-%6.4f  %7.4f+-%6.4f \n",
	   DYTools::massBinLimits[i], DYTools::massBinLimits[i+1],
	   nEventsv(i,0), nPassv(i,0),
	   accv(i,0), accErrv(i,0),
	   accBBv(i,0), accErrBBv(i,0),
	   accBEv(i,0), accErrBEv(i,0),
	   accEEv(i,0), accErrEEv(i,0));
      }
    }
  else
    printf("printout format for 2D not chosen");
  cout << endl;
 
  gBenchmark->Show("plotDYAcceptance");
}


//=== FUNCTION DEFINITIONS ======================================================================================

//--------------------------------------------------------------------------------------------------
