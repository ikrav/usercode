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
#include "../Include/plotFunctions.hh"

// define classes and constants to read in ntuple
#include "../Include/EWKAnaDefs.hh"
#include "../Include/TGenInfo.hh"
#endif

//=== FUNCTION DECLARATIONS ======================================================================================

//=== MAIN MACRO =================================================================================================

void plotDYFSRCorrections(const TString input) 
{
  gBenchmark->Start("plotDYFSRCorrections");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================
  
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
  
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================
  
  //  
  // Set up histograms
  //
  vector<TH1F*> hZMassv;//, hZMass2v, hZPtv, hZPt2v, hZyv, hZPhiv;  
  TH1F *hMassPreFsr = new TH1F("hMassPreFsr","",500,0,1500);
  TH1F *hMassPostFsr = new TH1F("hMassPostFsr","",500,0,1500);
  
  UInt_t   nZv = 0;
  TMatrixD nEventsv (DYTools::nMassBins,findMaxYBins());  
  TMatrixD nPassv   (DYTools::nMassBins,findMaxYBins());
  TMatrixD nCorrelv  (DYTools::nMassBins,findMaxYBins());
  TMatrixD corrv     (DYTools::nMassBins,findMaxYBins());
  TMatrixD corrErrv  (DYTools::nMassBins,findMaxYBins());

  nEventsv = 0;
  nPassv = 0;
  nCorrelv =0;

  char hname[100];
  for(UInt_t ifile = 0; ifile<fnamev.size(); ifile++) {
    sprintf(hname,"hZMass_%i",ifile); hZMassv.push_back(new TH1F(hname,"",500,0,1500)); hZMassv[ifile]->Sumw2();
  }

  // 
  // Read weights from a file
  //
  const bool useFewzWeights = true;
  const bool cutZPT100 = true;
  if(cutZPT100)
    cout << "NOTE: in MC, for Z/gamma* PT>100 the FEWZ weights for 80<PT<100 GeV are used!" << endl;

  TH2D *weights[DYTools::nMassBins2011];
  TH2D *weightErrors[DYTools::nMassBins2011];

  TFile fweights("../root_files/fewz/weights_stepwise_prec10-5_fine12.root");


  if( !fweights.IsOpen() ) assert(0);

  for(int i=0; i<DYTools::nMassBins2011; i++){
    TString hname1 = TString::Format("weight_%02d",i+1);
    weights[i] = (TH2D*)fweights.Get(hname1);
    hname1 = TString::Format("h_weighterror_%02d",i+1);
    weightErrors[i] = (TH2D*)fweights.Get(hname1);
  }


  //
  // Access samples and fill histograms
  //  
  TFile *infile=0;
  TTree *eventTree=0;  
    
  // Data structures to store info from TTrees
  mithep::TGenInfo *gen  = new mithep::TGenInfo();
  
  // loop over samples  

  int noReweight=0;
  int binProblemPreFsr=0;
  int binProblemPostFsr=0;
  int binProblemCorrel=0;


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
      //if (ientry>1000) break;
      genBr->GetEntry(ientry);

      double mass = gen->vmass;    // pre-FSR
      double massPostFsr = gen->mass;    // post-FSR
      if((mass < massLow) || (mass > massHigh)) continue;
      double y = gen->vy;    // pre-FSR
      double yPostFsr = gen->y;    // post-FSR
      if((y < yRangeMin) || (y > yRangeMax)) continue;

      int ibinM1D = DYTools::findMassBin2011(mass);
      // If mass is larger than the highest bin boundary
      // (last bin), use the last bin.
      if(ibinM1D == -1 && mass >= massBinLimits2011[nMassBins2011] )
	ibinM1D = nMassBins2011-1;
      int ibinM1DPostFsr = DYTools::findMassBin2011(massPostFsr); 

      int ibinM = DYTools::findMassBin(mass);
      // If mass is larger than the highest bin boundary
      // (last bin), use the last bin.
      if(ibinM == -1 && mass >= massBinLimits[nMassBins] )
	ibinM = nMassBins-1;
      int ibinMPostFsr = DYTools::findMassBin(massPostFsr);

      int ibinY = DYTools::findAbsYBin2D(ibinM, y);
      int ibinYPostFsr = DYTools::findAbsYBin2D(ibinMPostFsr, yPostFsr);

      // Find FEWZ-powheg reweighting factor 
      // that depends on pre-FSR Z/gamma* rapidity and pt
      double fewz_weight = 1.0;
      if(useFewzWeights){
	if(ibinM1D != -1 && ibinM1D < DYTools::nMassBins2011){
	  // Check that the virtual Z has Pt and Y within the
	  // weight map range. If it is in the underflow or overflow,
	  // move the index to point to the appropriate edge bin
	  int binX = weights[ibinM1D]->GetXaxis()->FindBin( gen->vpt );
	  int binY = weights[ibinM1D]->GetYaxis()->FindBin( gen->vy );
	  if(binX == 0) binX += 1;
	  if(binX == weights[ibinM1D]->GetNbinsX() + 1) binX -= 1;
	  if(binY == 0) binY += 1;
	  if(binY == weights[ibinM1D]->GetNbinsY() + 1) binY -= 1;
	  // Apply PT cut if needed
	  if( cutZPT100 ) 
	    if( binX == weights[ibinM1D]->GetNbinsX() )
	      binX = weights[ibinM1D]->GetNbinsX() - 1;
	  fewz_weight = weights[ibinM1D]->GetBinContent( binX, binY);
	}else
          noReweight++;
      }
    //       printf("mass= %f   pt= %f    Y= %f     weight= %f\n",gen->mass, gen->vpt, gen->vy, fewz_weight);

      if(ibinM != -1 && ibinM<nMassBins && ibinY!=-1 && ibinY<nYBins[ibinM])
	nEventsv(ibinM,ibinY) += scale * gen->weight * fewz_weight;
      else if(ibinM >= nMassBins)
        binProblemPreFsr++;


      if(ibinMPostFsr != -1 && ibinMPostFsr < nMassBins &&  ibinYPostFsr!=-1 && ibinYPostFsr<nYBins[ibinMPostFsr])
	nPassv(ibinMPostFsr,ibinYPostFsr) += scale * gen->weight * fewz_weight;
      else if(ibinMPostFsr >= nMassBins || ibinYPostFsr>=nYBins[ibinMPostFsr]) {
	// Do nothing: post-fsr mass could easily be below the lowest edge of mass range
// 	cout << "ERROR: binning problem post-FSR, bin=" << ibinPostFsr << "  mass=" << massPostFsr << endl;
        binProblemPostFsr++;
      }

      if (ibinM==ibinMPostFsr && ibinY==ibinYPostFsr)
      {
      if(ibinM != -1 && ibinM < nMassBins && ibinY!=-1 && ibinY<nYBins[ibinM])
          nCorrelv(ibinM,ibinY) += scale * gen->weight * fewz_weight;
      else if(ibinM >= nMassBins || ibinY>=nYBins[ibinM])
	binProblemCorrel++;
      }

      hZMassv[ifile]->Fill(mass,scale * gen->weight * fewz_weight);
      hMassPreFsr->Fill(mass, scale*gen->weight * fewz_weight);
      hMassPostFsr->Fill(massPostFsr, scale*gen->weight * fewz_weight);

    }   
    delete infile;
    infile=0, eventTree=0;
  }
  delete gen;
  cout << "Error: binning problem FEWZ bins. For " << noReweight << "  events"<< endl;
  cout << "ERROR: binning problem pre-FSR. For " << binProblemPreFsr << "  events" << endl;
  cout << "ERROR: binning problem post-FSR. For " << binProblemPostFsr << "  events" << endl;
  cout << "ERROR: binning problem Correlation. For " << binProblemCorrel << "  events" << endl;

  corrv      = 0;
  corrErrv   = 0;
  for(int i=0; i<DYTools::nMassBins; i++)
    for (int j=0; j<DYTools::nYBins[i]; j++){
      if(nEventsv(i,j) != 0){
        corrv(i,j) = nPassv(i,j)/nEventsv(i,j);
        corrErrv(i,j) = corrv(i,j) * sqrt( 1.0/nPassv(i,j) + 1.0/nEventsv(i,j)- 2*nCorrelv(i,j)/(nPassv(i,j)*nEventsv(i,j)) );
//        corrErrv[i] = corrv[i] * sqrt( 1.0/nPassv[i] + 1.0/nEventsv[i] );
//       corrErrv[i] = sqrt(corrv[i]*(1-corrv[i])/nEventsv[i]);
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
  plotZMass1.Draw(c);
  SaveCanvas(c, "zmass1");

  // Pre FSR vs post-FSR plots
  TCanvas *c2 = MakeCanvas("c2","c2",800,600);
  CPlot plotOverlay("overlay","","m(ee) [GeV/c^{2}]","Events");
  plotOverlay.AddHist1D(hMassPreFsr,"pre-FSR","hist",kBlue); 
  plotOverlay.AddHist1D(hMassPostFsr,"post-FSR","hist",kRed); 
  plotOverlay.SetLogy();
  plotOverlay.Draw(c2);
  SaveCanvas(c2, "overlay");

  PlotMatrixVariousBinning(corrv, "N_PosrFsr/N_PreFsr", "LEGO2");


          
  // Store constants in the file
  TString outputDir(TString("../root_files/constants/")+dirTag);
  gSystem->mkdir(outputDir,kTRUE);
  TString fsrConstFileName(outputDir+TString("/fsr_constants.root"));

  TFile fa(fsrConstFileName,"recreate");
  //correctionGraph->Write("CORR_FSR");
  corrv.Write("fsrCorrectionArray");
  corrErrv.Write("fsrCorrectionErrArray");
  fa.Close();

  //--------------------------------------------------------------------------------------------------------------
  // Summary print out
  //==============================================================================================================
  //cout << endl;
  //cout << "*" << endl;
  //cout << "* SUMMARY" << endl;
  //cout << "*--------------------------------------------------" << endl;
  //cout << endl; 
  
  //cout << labelv[0] << " file: " << fnamev[0] << endl;
  //cout << "     Number of generated events: " << nZv << endl;
  //printf(" mass bin    preselected      passed     FSR correction\n");
  //for(int i=0; i<DYTools::nMassBins; i++){
   // printf(" %4.0f-%4.0f   %10.0f   %10.0f   %7.4f+-%6.4f \n",
	//   DYTools::massBinLimits[i], DYTools::massBinLimits[i+1],
	 //  nEventsv[i], nPassv[i], 
          // corrv[i], corrErrv[i]);
  //}

  //cout << endl;

  gBenchmark->Show("plotDYFSRCorrections");
}
