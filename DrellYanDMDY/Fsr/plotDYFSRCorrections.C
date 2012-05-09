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

#include "../Include/EventSelector.hh"
#include "../Include/FEWZ.hh"
#include "../Include/UnfoldingTools.hh"

#endif

//=== FUNCTION DECLARATIONS ======================================================================================

//=== MAIN MACRO =================================================================================================

void plotDYFSRCorrections(const TString input, bool sansAcc=0, int debugMode=0) 
{
  gBenchmark->Start("plotDYFSRCorrections");

  if (debugMode) std::cout << "\n\n\tDEBUG MODE is ON\n\n";

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================
   
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
  
  UInt_t   nZ = 0;
  Double_t nZweighted=0;
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
    sprintf(hname,"hZMass_%i",ifile); 
    hZMassv.push_back(new TH1F(hname,"",500,0,1500)); 
    hZMassv[ifile]->Sumw2();
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
  TH2D *weights[DYTools::_nMassBins2011]; // temporary
  TH2D *weightErrors[DYTools::_nMassBins2011]; // temporary
  if (!new_fewz_code) {
  if(cutZPT100)
    cout << "NOTE: in MC, for Z/gamma* PT>100 the FEWZ weights for 80<PT<100 GeV are used!" << endl;
  TFile fweights("../root_files/fewz/weights_stepwise_prec10-5_fine12.root");
  if( !fweights.IsOpen() ) assert(0);
  for(int i=0; i<DYTools::_nMassBins2011; i++){
    TString hname1 = TString::Format("weight_%02d",i+1);
    weights[i] = (TH2D*)fweights.Get(hname1);
    hname1 = TString::Format("h_weighterror_%02d",i+1);
    weightErrors[i] = (TH2D*)fweights.Get(hname1);
  }
  }

  //
  // Access samples and fill histograms
  //  
  TFile *infile=0;
  TTree *eventTree=0;  
    
  // Data structures to store info from TTrees
  mithep::TGenInfo *gen  = new mithep::TGenInfo();
  int binProblemFEWZ=0;
  int binProblemPreFsr=0;
  int binProblemPostFsr=0;
  int binProblemCorrel=0;

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
    AdjustXSectionForSkim(infile,xsecv[ifile],eventTree->GetEntries(),1);
    lumiv[ifile] = eventTree->GetEntries()/xsecv[ifile];
    double scale = lumiv[0]/lumiv[ifile];
//     if(ifile != 0) scale *= 0.87;
    cout << "       -> sample weight is " << scale << endl;

    // Set branch address to structures that will store the info  
    eventTree->SetBranchAddress("Gen",&gen);
    TBranch *genBr = eventTree->GetBranch("Gen");
  
    // loop over events    
    nZ += eventTree->GetEntries();
    nZweighted += scale * eventTree->GetEntries();
    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
      if (debugMode && (ientry>10000)) break; // debug option
      genBr->GetEntry(ientry);

      // if sansAcc mode, Discard events that are not in kinematics acceptance
      if (sansAcc)
        {
           if( ! (DYTools::isBarrel(gen->eta_1) || 
		  DYTools::isEndcap(gen->eta_1) ) ) continue;
           if( ! (DYTools::isBarrel(gen->eta_2) || 
		  DYTools::isEndcap(gen->eta_2) ) ) continue;
           // both above 10 GeV and at least one above 20 GeV
           if( ! (gen->pt_1 > 10 && gen->pt_2 > 10) ) continue;
           if( ! (gen->pt_1 > 20 || gen->pt_2 > 20) ) continue;
        }

      double mass = gen->vmass;    // pre-FSR
      double massPostFsr = gen->mass;    // post-FSR
      if((mass < massLow) || (mass > massHigh)) continue;
      double y = gen->vy;    // pre-FSR
      double yPostFsr = gen->y;    // post-FSR
      if((y < yRangeMin) || (y > yRangeMax)) continue;

      int ibinM1D = DYTools::_findMassBin2011(mass);      // temporary
      // If mass is larger than the highest bin boundary
      // (last bin), use the last bin.
      if(ibinM1D == -1 && mass >= _massBinLimits2011[_nMassBins2011] )
	ibinM1D = _nMassBins2011-1;
      //int ibinM1DPostFsr = DYTools::_findMassBin2011(massPostFsr);  // temporary

      int ibinM = DYTools::findMassBin(mass);
      // If mass is larger than the highest bin boundary
      // (last bin), use the last bin.
      if(ibinM == -1 && mass >= massBinLimits[nMassBins] )
	ibinM = nMassBins-1;
      int ibinMPostFsr = DYTools::findMassBin(massPostFsr);

      int ibinY = DYTools::findAbsYBin(ibinM,y);
      int ibinYPostFsr = DYTools::findAbsYBin(ibinMPostFsr,yPostFsr);

      // Find FEWZ-powheg reweighting factor 
      // that depends on pre-FSR Z/gamma* rapidity and pt
      double fewz_weight = 1.0;
      if(useFewzWeights){
	if (new_fewz_code) {
	  fewz_weight=fewz.getWeight(gen->vmass,gen->vpt,gen->vy);
	}
	else {
	if(ibinM1D != -1 && ibinM1D < DYTools::_nMassBins2011){
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
          binProblemFEWZ++;
      }
     //       printf("mass= %f   pt= %f    Y= %f     weight= %f\n",gen->mass, gen->vpt, gen->vy, fewz_weight);
      }

      if(ibinM != -1 && ibinM<nMassBins && ibinY!=-1 && ibinY<nYBins[ibinM])
	nEventsv(ibinM,ibinY) += scale * gen->weight * fewz_weight;
      else if(ibinM >= nMassBins || ibinY>=nYBins[ibinM])
        binProblemPreFsr++;

      if(ibinMPostFsr!=-1 && ibinMPostFsr<nMassBins &&  ibinYPostFsr!=-1 
	 && ibinYPostFsr<nYBins[ibinMPostFsr])
	nPassv(ibinMPostFsr,ibinYPostFsr) += scale * gen->weight * fewz_weight;
      else if(ibinMPostFsr >= nMassBins || 
	      ibinYPostFsr>=nYBins[ibinMPostFsr]) {
	// Do nothing: post-fsr mass could easily be below the lowest edge of mass range
        binProblemPostFsr++;
      }

      if (ibinM==ibinMPostFsr && ibinY==ibinYPostFsr)
        {
          if(ibinM != -1 && ibinM<nMassBins && ibinY!=-1 && 
	     ibinY<nYBins[ibinM])
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
  cout << "Error: binning problem FEWZ bins, " << binProblemFEWZ 
       <<"  events"<< endl;
  cout << "ERROR: binning problem, " << binProblemPreFsr<<" events"<< endl;
  cout << "ERROR: binning problem Correlation, " << binProblemCorrel 
       <<"  events" << endl;

  corrv      = 0;
  corrErrv   = 0;
  for(int i=0; i<DYTools::nMassBins; i++)
    for (int j=0; j<DYTools::nYBins[i]; j++)
      {
        if(nEventsv(i,j) != 0)
          {
            corrv(i,j) = nPassv(i,j)/nEventsv(i,j);
            corrErrv(i,j) = corrv(i,j) * 
	      sqrt( 1.0/nPassv(i,j) + 1.0/nEventsv(i,j)- 
		    2*nCorrelv(i,j)/(nPassv(i,j)*nEventsv(i,j)) );
            //corrErrv[i] = corrv[i] * sqrt( 1.0/nPassv[i] + 1.0/nEventsv[i] );
            //corrErrv[i] = sqrt(corrv[i]*(1-corrv[i])/nEventsv[i]);
          }
       }

  //--------------------------------------------------------------------------------------------------------------
  // Make plots 
  //==============================================================================================================  

  // prepare tag
  TString addStr=TString("_") + analysisTag;
  if (sansAcc) addStr.Append("_sans_acc");
  //else addStr=""; 

  // Prepare directories
  TString outputDir(TString("../root_files/constants/")+dirTag);
  gSystem->mkdir(outputDir,kTRUE);
  TString fsr_constants="/fsr_constants";
  fsr_constants+=addStr;   // addStr contains analysisTag
  fsr_constants+=".root";

  TString fileNamePlots=outputDir + fsr_constants;
  fileNamePlots.Replace(fileNamePlots.Index(".root"),
			sizeof(".root"),"_plots.root");
  TFile *filePlots=new TFile(fileNamePlots,"recreate");
  if (!filePlots || !filePlots->IsOpen()) {
    std::cout << "failed to create file <" << fileNamePlots << ">\n";
    throw 2;
  }

  TString canvName=TString("canvFsr_zmass1") + addStr;
  TCanvas *c = MakeCanvas(canvName,canvName,800,600);

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
  TString zmass1="zmass1";
  zmass1+=addStr;
  SaveCanvas(c, zmass1);
  c->Write();

  // Pre FSR vs post-FSR plots
  canvName = TString("canv_preFsr_vs_postFSR") + addStr;
  TCanvas *c2 = MakeCanvas(canvName,canvName,800,600);
  CPlot plotOverlay("overlay","","m(ee) [GeV/c^{2}]","Events");
  plotOverlay.AddHist1D(hMassPreFsr,"pre-FSR","hist",kBlue); 
  plotOverlay.AddHist1D(hMassPostFsr,"post-FSR","hist",kRed); 
  plotOverlay.SetLogy();
  TString overlay="overlay";
  overlay+=addStr;
  plotOverlay.Draw(c2);
  SaveCanvas(c2, overlay);
  c2->Write();

  plotOverlay.SetLogx();
  canvName.Append("_logM");
  overlay.Append("_logM");
  plotOverlay.SetXRange(1,DYTools::massBinLimits[DYTools::nMassBins]);
  TCanvas *c3 = MakeCanvas(canvName,canvName,800,600);
  plotOverlay.Draw(c3); // prepare legend
  plotOverlay.TransLegend(-1.,0.);        // incorrect placement?!
  plotOverlay.Draw(c3);
  SaveCanvas(c3, overlay);
  c3->Write();
  

  TString NoverN="N_PosrFsr_over_N_PreFsr";
  NoverN+=addStr;
  std::cout<<"I'M going to save plots"<<std::endl;
  PlotMatrixVariousBinning(corrv, NoverN, "LEGO2",filePlots);

  if (filePlots) filePlots->Close();


  // Store constants in the file
  TString fsrConstFileName(outputDir+fsr_constants);
  TFile fa(fsrConstFileName,"recreate");
  //correctionGraph->Write("CORR_FSR");
  corrv.Write("fsrCorrectionMatrix");
  corrErrv.Write("fsrCorrectionErrMatrix");
  unfolding::writeBinningArrays(fa);
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
  //cout << "     Number of generated events:    " << nZ << endl;
  //char buf[30];
  //sprintf(buf,"%3.1lf",nZweighted);
  //cout << "     Number of weighted gen.events: " << buf << endl;
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
