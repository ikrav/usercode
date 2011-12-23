//================================================================================================
//
// Prepare binned histograms with signal and background events for further analysis.
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TCanvas.h>                // class for drawing
#include <TH1F.h>                   // 1D histograms
#include <TBenchmark.h>             // class to track macro running statistics
#include <TVector3.h>               // 3D vector class
#include <TArrayD.h>
#include <TVectorD.h>
#include <TLorentzVector.h>         // 4-vector class
#include <TRandom.h>
#include <TDatime.h>                // time stamp
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings

#include "../Include/CPlot.hh"          // helper class for plots
#include "../Include/MitStyleRemix.hh"  // style settings for drawing
#include "../Include/MyTools.hh"        // miscellaneous helper functions
#include "../Include/CSample.hh"        // helper class for organizing input ntuple files

// define structures to read in ntuple
#include "../Include/ZeeData.hh"

#include "../Include/ElectronEnergyScale.hh"        // energy scale correction

#include "../Include/DYTools.hh"

#endif

// Forward declarations
void drawRapidityInMassSlices(vector<TMatrixD*> &yields, 
			      vector<TMatrixD*> &yieldsSumw2,
			      vector<TString>   &sampleTags,
			      int nhist, TCanvas *canvas);

//=== MAIN MACRO =================================================================================================

void prepareYields(const TString conf  = "data_plot.conf") 
{  
  gBenchmark->Start("prepareYields");

  
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================

  TString  outputDir;         // output directory
  Double_t lumi;              // luminosity (pb^-1)
  Bool_t   doWeight;          // weight events?
  TString  escaleTag;         // Energy scale calibrations tag
  TString  format;            // plot format

  vector<TString>  snamev;    // sample name (for output file)
  vector<CSample*> samplev;   // data/MC samples
    
  //
  // parse .conf file
  //
  ifstream ifs;
  ifs.open(conf.Data());
  assert(ifs.is_open());
  string line;
  Int_t state=0;
  while(getline(ifs,line)) {
    if(line[0]=='#') continue;
    if(line[0]=='%') { 
      state++; 
      continue; 
    }
    if(line[0]=='$') {
      samplev.push_back(new CSample());
      stringstream ss(line);
      string chr;
      string sname;
      Int_t color;
      ss >> chr >> sname >> color;
      string label = line.substr(line.find('@')+1);
      snamev.push_back(sname);
      samplev.back()->label = label;
      samplev.back()->color = color;
      continue;
    }
    
    if(state==0) {  // general settings
      stringstream ss1(line); ss1 >> lumi;
      getline(ifs,line);
      stringstream ss2(line); ss2 >> doWeight;
      getline(ifs,line);
      outputDir = TString(line);
      getline(ifs,line);
      stringstream ss3(line); ss3 >> escaleTag;
      getline(ifs,line);
      format = TString(line);
      
    } else if(state==1) {  // define data sample
      string fname;
      Double_t xsec;
      string json;
      stringstream ss(line);
      ss >> fname >> xsec >> json;
      samplev.back()->fnamev.push_back(fname);
      samplev.back()->xsecv.push_back(xsec);
      samplev.back()->jsonv.push_back(json);
    
    } else if(state==2) {  // define MC samples
      string fname;
      Double_t xsec;
      stringstream ss(line);
      ss >> fname >> xsec;
      samplev.back()->fnamev.push_back(fname);
      samplev.back()->xsecv.push_back(xsec);
    }
  }
  ifs.close();
  
  ElectronEnergyScale::CalibrationSet calibrationSet 
    = ElectronEnergyScale::UNDEFINED;
  if( escaleTag == TString("Date20110901_EPS11_default")){
    calibrationSet = ElectronEnergyScale::Date20110901_EPS11_default;
  }else{
    printf("Failed to match escale calibration. Tag: >>%s<<\n", escaleTag.Data());
    assert(0);
  }

  // sOutDir is a static data member in the CPlot class.
  // There is a strange crash of the whole ROOT session well after
  // this script is executed when one attempts to exit ROOT, with 
  // a dump of memory map. This happens only on UNL Tier3, but
  // there is no sign of a problem on any other computer.
  //   The consequence of this variable is not set is that the plots
  // will be created in the local directory rather than the
  // one configured through sOutDir.
//   CPlot::sOutDir = outputDir + TString("/plots");

  Bool_t hasData = (samplev[0]->fnamev.size()>0);
    
  //
  // Canvas dimensions
  //
  Int_t canw=800, canh=600;
  
  
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //============================================================================================================== 

  vector<TMatrixD*> yields;
  vector<TMatrixD*> yieldsSumw2;
  int maxYBins = DYTools::findMaxYBins();
  vector<TString> sampleTags;

  vector<TH1F*>    hMassv;
  vector<TH1F*>    hMassBinsv;
  vector<Double_t> nSelv;
  vector<Double_t> nSelVarv;  

  char hname[100];
  for(UInt_t isam=0; isam<samplev.size(); isam++) {
    TMatrixD *yieldsMatrix = new TMatrixD(DYTools::nMassBins2D, maxYBins);
    (*yieldsMatrix) = 0;
    yields.push_back(yieldsMatrix);
    //
    TMatrixD *yieldsSumw2Matrix = new TMatrixD(DYTools::nMassBins2D, maxYBins);
    (*yieldsSumw2Matrix) = 0;
    yieldsSumw2.push_back(yieldsSumw2Matrix);
    //
    sampleTags.push_back(snamev[isam]);

    sprintf(hname,"hMass_%i",isam);  
    hMassv.push_back(new TH1F(hname,"",1490,10,1500));  
    hMassv[isam]->Sumw2();
    sprintf(hname,"hMassBins_%i",isam);  
    hMassBinsv.push_back(new TH1F(hname,"",
				  DYTools::nMassBins2D, 
				  DYTools::massBinLimits2D));
    hMassBinsv[isam]->Sumw2();
			 
    nSelv.push_back(0);
    nSelVarv.push_back(0);  
  }
  
  // 
  // Set up energy scale corrections
  //
  ElectronEnergyScale escale(calibrationSet);
  escale.print();

  ZeeData data;
  TRandom random;

  // Open file with number of PV distributions for pile-up reweighting
  const TString fnamePV = outputDir+TString("/npv.root");
  TFile *pvfile = new TFile(fnamePV);
  assert(pvfile);
  TH1F *hPVData = 0;
  if(hasData){ 
    hPVData = (TH1F*)pvfile->Get("hNGoodPV_data"); assert(hPVData); 
  }

  //
  //  Diboson backgrounds need to be saved separately, but plotted
  //  together. Detect whether there are separate ww/wz/zz contribution,
  //  and whether merging is needed later. Of course, this relies on
  // the fact that the file data.conf has names ww, wz, zz for those
  // contributions.
  int nDibosonSamples = 0;
  bool mergeDibosons = false;
  for(UInt_t isam=0; isam<samplev.size(); isam++) {
    if( snamev[isam] == "ww" || snamev[isam] == "wz" || snamev[isam] == "zz")
      nDibosonSamples++;
  }
  if(nDibosonSamples==3)
    mergeDibosons = true;

  //
  // Access samples and fill histograms
  //  
  TFile *infile=0;
  TTree *eventTree=0; 
  
  for(UInt_t isam=0; isam<samplev.size(); isam++) {
    if((isam==0) && !hasData) continue;
    
    const TString fname = outputDir + TString("/ntuples/") + snamev[isam] + TString("_select.root");
    cout << "Processing " << fname << "..." << endl;   
    infile = new TFile(fname);
    assert(infile); 

    // Prepare weights for pile-up reweighting for MC
    TH1F *hPVThis = (TH1F*) pvfile->Get(TString("hNGoodPV_")+snamev[isam]); assert(hPVThis);
    // Normalize or not? Not clear
    hPVThis->Scale( hPVData->GetSumOfWeights()/hPVThis->GetSumOfWeights());
    TH1F *puWeights = (TH1F*)hPVData->Clone("puWeights");
    puWeights->Divide(hPVThis);
//     for(int i=1; i<=puWeights->GetNbinsX(); i++)
//       printf(" %f    %f    %f\n",puWeights->GetBinCenter(i),puWeights->GetBinContent(i),puWeights->GetBinError(i));

    // Get the TTree and set branch address
    eventTree = (TTree*)infile->Get("Events"); assert(eventTree); 
    eventTree->SetBranchAddress("Events",&data);


    TMatrixD *thisSampleYields = yields.at(isam);
    TMatrixD *thisSampleYieldsSumw2 = yieldsSumw2.at(isam);

    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
      eventTree->GetEntry(ientry);
      Double_t weight = data.weight;
      
      // Any extra weight factors:
      double weightPU = puWeights->GetBinContent( puWeights->FindBin( data.nPV ));
      weight *= weightPU;
      
      // If This is MC, add extra smearing to the mass
      // We apply extra smearing to all MC samples: it is may be
      // not quite right for fake electron backgrounds, but these
      // are not dominant, and in any case we do not have corrections
      // for fake electrons.
      if(isam!=0) {            	    
	double smearingCorrection = escale.generateMCSmear(data.scEta_1, data.scEta_2);
	data.mass = data.mass + smearingCorrection;
      }

      // Find the 2D bin for this event:
      int massBin = findMassBin2D(data.mass);
      int yBin    = findAbsYBin2D(massBin, data.y);
      (*thisSampleYields)(massBin,yBin) += weight;
      (*thisSampleYieldsSumw2)(massBin,yBin) += weight*weight;

      hMassv[isam]->Fill(data.mass,weight);
      hMassBinsv[isam]->Fill(data.mass,weight);

      nSelv[isam] += weight;
      nSelVarv[isam] += weight*weight;
      
    }
    delete infile;
    infile=0, eventTree=0;
  }

  //--------------------------------------------------------------------------------------------------------------
  // Make plots
  //==============================================================================================================

  // Buffers for labels and comments
  char ylabel[100];   // y-axis label
  char lumitext[50];
  if(lumi>0) {
    if(lumi<1) { sprintf(lumitext,"#int#font[12]{L}dt = %.0f nb^{-1}",1000.*lumi); }
    if(lumi<1000) { sprintf(lumitext,"#int#font[12]{L}dt = %.0f nb^{-1}",1000.*lumi); }
    else       { sprintf(lumitext,"#int#font[12]{L}dt = %.3g fb^{-1}",lumi/1000.0); }
  }

  // Merge diboson histograms if needed
  TH1F *hMassBinsDibosons = (TH1F*)hMassBinsv[1]->Clone("hMassBinsDibosons");
  TH1F *hMassDibosons = (TH1F*)hMassv[1]->Clone("hMassDibosons");
  hMassDibosons->Reset();
  Int_t colorDibosons = 1;
  TString labelDibosons = "WW/WZ/ZZ";
  for(UInt_t isam=0; isam<samplev.size(); isam++) {
    if( snamev[isam] == "ww" || snamev[isam] == "wz" || snamev[isam] == "zz"){
      hMassDibosons->Add(hMassv[isam]);
      hMassBinsDibosons->Add(hMassBinsv[isam]);
      // Use color of the last diboson entry
      colorDibosons = samplev[isam]->color;
    }
  }

  // Additional normalizatoin for MC:
  //
  // Ideally, we would normalize all MC samples to data luminosity
  // In practice, however, it is not easy because of two reasons:
  //  - the luminosity is known with an error (systematic shift of 6% is
  //       for example suspected in mid-2011)
  //  - data/MC scale factors for efficiency to select events may
  //       move normalization off by another 5%
  // Therefore, we normalize total MC to the data Z peak. This gives us
  // the scale factor that is applied to all samples. 
  // In the following calculation it is assumed that the first histogram
  // is data and the last is signal MC.  
  TH1F *totalMCMass = (TH1F*)hMassv[0]->Clone("totalMCMass");
  totalMCMass->Reset();
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    totalMCMass->Add(hMassv[isam]);
  }
  double massNormMin = 60.0;
  double massNormMax = 120.0;
  double dataOverMc = hMassv[0]->Integral(hMassv[0]->FindBin(massNormMin+0.001),
					  hMassv[0]->FindBin(massNormMax-0.001)) /
    totalMCMass->Integral(totalMCMass->FindBin(massNormMin+0.001),
			  totalMCMass->FindBin(massNormMax-0.001));
  printf("data to MC extra correction from Z peak normalization: %f\n",dataOverMc);
  
  // Rescale all MC samples. This is not totally proper for fake lepton
  // backgrounds, but ok for backgrounds with true leptons, and those are dominant
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    hMassv[isam]->Scale(dataOverMc);
    hMassBinsv[isam]->Scale(dataOverMc);
    printf("  MC %s IS RESCALED for plotting\n", snamev[isam].Data());
  }
  
  //
  // Draw mass spectrum without rapidity binning
  //
  // First, draw the histograms with fine binning.
  TCanvas *c1 = MakeCanvas("c1","c1",canw,canh);
  sprintf(ylabel,"Events / %.1f GeV/c^{2}",hMassv[0]->GetBinWidth(1));
  CPlot plotMass("mass","","m(e^{+}e^{-}) [GeV/c^{2}]",ylabel);
  plotMass.SetLogx();
  if(hasData) { plotMass.AddHist1D(hMassv[0],samplev[0]->label,"E"); }
  // Do not draw separately dibosons, but draw the merged histogram if needed
  if(mergeDibosons)
    plotMass.AddToStack(hMassDibosons, labelDibosons, colorDibosons);
  for(UInt_t isam=1; isam<samplev.size(); isam++){
    if( !(mergeDibosons && (snamev[isam]=="ww" || snamev[isam]=="wz" || snamev[isam]=="zz")))
      plotMass.AddToStack(hMassv[isam],samplev[isam]->label,samplev[isam]->color);
  }
  plotMass.SetLegend(0.75,0.55,0.98,0.9);
  if(lumi>0) plotMass.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
  if(hasData){
    hMassv[0]->GetXaxis()->SetMoreLogLabels();
    hMassv[0]->GetXaxis()->SetNoExponent();
  }
  plotMass.SetLogy();
  plotMass.SetLogx();
  plotMass.Draw(c1);
  plotMass.SetYRange(1.0,300000);  
  plotMass.Draw(c1,kFALSE,format);

  // Second, draw with the mass binning used in the analysis
  TCanvas *c2 = MakeCanvas("c2","c2",canw,canh);
  CPlot plotMassBins("mass","","m(e^{+}e^{-}) [GeV/c^{2}]","Events");
  plotMassBins.SetLogx();
  if(hasData) { plotMassBins.AddHist1D(hMassBinsv[0],samplev[0]->label,"E"); }
  // Do not draw separately dibosons, but draw the merged histogram if needed
  if(mergeDibosons)
    plotMassBins.AddToStack(hMassDibosons, labelDibosons, colorDibosons);
  for(UInt_t isam=1; isam<samplev.size(); isam++){
    if( !(mergeDibosons && (snamev[isam]=="ww" || snamev[isam]=="wz" || snamev[isam]=="zz")))
      plotMassBins.AddToStack(hMassBinsv[isam],samplev[isam]->label,samplev[isam]->color);
  }
  plotMassBins.SetLegend(0.75,0.55,0.98,0.9);
  if(lumi>0) plotMassBins.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
  if(hasData){
    hMassBinsv[0]->GetXaxis()->SetMoreLogLabels();
    hMassBinsv[0]->GetXaxis()->SetNoExponent();
  }
  plotMassBins.SetLogy();
  plotMassBins.SetLogx();
  plotMassBins.Draw(c2);
  plotMassBins.SetYRange(1.0,2000000);  
  plotMassBins.Draw(c2,kFALSE,format);
 
  //
  // Draw flattened distribution
  //
  // Create the histograms from the yields arrays
  int flatIndexMax = DYTools::getNumberOf2DBins();
  TH1F *hFlattened[samplev.size()];
  for(UInt_t isam=0; isam < samplev.size(); isam++){
    sprintf(hname, "hFlattanedMass_%i", isam);
    hFlattened[isam] = new TH1F(hname,"",flatIndexMax, 0.0, 1.0*flatIndexMax);
    TMatrixD *thisSampleYields = yields.at(isam);
    TMatrixD *thisSampleYieldsSumw2 = yieldsSumw2.at(isam);
    for(int im = 0; im < DYTools::nMassBins2D; im++){
      for(int iy = 0; iy < DYTools::nYBins[im]; iy++){
	int iflat = findIndexFlat(im, iy);
	hFlattened[isam]->SetBinContent(iflat, (*thisSampleYields)(im,iy) );
	hFlattened[isam]->SetBinError(iflat, sqrt((*thisSampleYieldsSumw2)(im,iy)) );
      }
    }
  }
  // Merge dibosons
  TH1F *hFlattenedDibosons = (TH1F*)hFlattened[1]->Clone("hFlattenedDibosons");
  hFlattenedDibosons->Reset();
  for(UInt_t isam=0; isam<samplev.size(); isam++) {
    if( snamev[isam] == "ww" || snamev[isam] == "wz" || snamev[isam] == "zz"){
      hFlattenedDibosons->Add(hFlattened[isam]);
    }
  }

  // Draw the flattened figure.
  TCanvas *c3 = MakeCanvas("c3","c3",canw,canh);
  CPlot plotFlattened("mass","","m(e^{+}e^{-}) [GeV/c^{2}]","Events");
  if(hasData) { plotFlattened.AddHist1D(hFlattened[0],samplev[0]->label,"E"); }
  // Do not draw separately dibosons, but draw the merged histogram if needed
  if(mergeDibosons)
    plotFlattened.AddToStack(hFlattenedDibosons, labelDibosons, colorDibosons);
  for(UInt_t isam=1; isam<samplev.size(); isam++){
    if( !(mergeDibosons && (snamev[isam]=="ww" || snamev[isam]=="wz" || snamev[isam]=="zz")))
      plotFlattened.AddToStack(hFlattened[isam],samplev[isam]->label,samplev[isam]->color);
  }
  plotFlattened.SetLegend(0.75,0.55,0.98,0.9);
  if(lumi>0) plotFlattened.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
  if(hasData){
    hFlattened[0]->GetXaxis()->SetMoreLogLabels();
    hFlattened[0]->GetXaxis()->SetNoExponent();
  }
  plotFlattened.SetLogy();
  plotFlattened.Draw(c3);
  plotFlattened.SetYRange(1.0,200000);  
  plotFlattened.Draw(c3,kFALSE,format);

  //
  // Draw rapidity in mass slices
  //

  TCanvas *c4 = MakeCanvas("c4","c4",canw,canh);
  drawRapidityInMassSlices(yields, yieldsSumw2, sampleTags, DYTools::nMassBins2D, c4);

  //--------------------------------------------------------------------------------------------------------------
  // Save the results
  //==============================================================================================================

  // Pack info into writable objects
  TVectorD massBinLimits(nMassBins2D+1);
  TVectorD rapidityBinning(nMassBins2D+1);
  for(int i=0; i <= nMassBins2D; i++){
    massBinLimits(i) = DYTools::massBinLimits2D[i];
    rapidityBinning(i) = DYTools::nYBins[i];
  }
  // This dummy object is only needed to convey the number
  // of samples considered. The method below is rather convoluted,
  // but I do not know a better one. Ideally, we would just store
  // a list of strings, each string containing the sample name.
  TVectorD dummySampleCount(sampleTags.size());
  dummySampleCount = 0;

  TString outputDirYields(outputDir.Data());
  outputDirYields.ReplaceAll("selected_events","yields");
  gSystem->mkdir(outputDirYields,kTRUE);
  TString fNameOutYields(outputDirYields+TString("/yields"));
  fNameOutYields += ".root";
  TFile fYields( fNameOutYields, "recreate" );
  massBinLimits      .Write("massBinLimits");
  rapidityBinning    .Write("rapidityBinning");
  dummySampleCount   .Write("dummySampleCount");
  char objName[100];
  for(UInt_t isam = 0; isam < yields.size(); isam++){
    sprintf(objName,"sample_name_%i",isam);
    TObjString *sampleNameStorable = new TObjString( sampleTags.at(isam) );
    sampleNameStorable->Write(objName);
    sprintf(objName,"yields_%s",sampleTags.at(isam).Data());
    yields[isam]       ->Write(objName);
    sprintf(objName,"yieldsSumw2_%s",sampleTags.at(isam).Data());
    yieldsSumw2[isam]  ->Write(objName);
  }
  fYields.Close();
  
  //--------------------------------------------------------------------------------------------------------------
  // Summary print out
  //==============================================================================================================
  cout << endl;
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  cout << endl;
  
  if(hasData) {
    cout << "   Data: " << setprecision(1) << fixed << nSelv[0] << " Drell-Yan events!" << endl;
    for(UInt_t ifile=0; ifile<samplev[0]->fnamev.size(); ifile++)
      cout << "     " << samplev[0]->fnamev[ifile] << endl;
  }
  cout << endl;
  if(samplev.size()>1) {
    cout << "   MC:" << endl;
    for(UInt_t isam=1; isam<samplev.size(); isam++) {
      cout << "      " << snamev[isam] << setw(8) << setprecision(3) << fixed << nSelv[isam] << " +/- ";
      cout << setw(4) << setprecision(3) << fixed << sqrt(nSelVarv[isam]) << endl;
      for(UInt_t ifile=0; ifile<samplev[isam]->fnamev.size(); ifile++) {
        cout << "         " << samplev[isam]->fnamev[ifile] << endl;
      }
    }
    cout << endl;
  }
  cout << endl;

  //
  // Summary printout in mass bins, integrated over rapidity
  //
  // Add yields over rapidity bins
  double totalData            [DYTools::nMassBins2D];
  double totalSignalMC        [DYTools::nMassBins2D];
  double totalSignalMCError   [DYTools::nMassBins2D];
  double totalBg     [DYTools::nMassBins2D];
  double totalBgError[DYTools::nMassBins2D];
  for( int im=0; im<DYTools::nMassBins2D; im++){
    totalData         [im] = 0;
    totalSignalMC     [im] = 0;
    totalSignalMCError[im] = 0;
    totalBg           [im] = 0;
    totalBgError      [im] = 0;
    for(int iy = 0; iy < DYTools::nYBins[im]; iy++){
      for( UInt_t isam = 0; isam < yields.size(); isam++){
	if( sampleTags.at(isam) == TString("data") ){
	  totalData[im] += (*yields.at(isam))(im,iy);
	}else if ( sampleTags.at(isam) == TString("zee") ){
	  totalSignalMC[im] += (*yields.at(isam))(im,iy);
	  totalSignalMCError[im] += (*yieldsSumw2.at(isam))(im,iy);
	}else{
	  // what remains are background samples
	  totalBg[im] += (*yields.at(isam))(im,iy);
	  totalBgError[im] += (*yieldsSumw2.at(isam))(im,iy);
	}
      } // end loop over samples
    } // end loop over rapidity bins
    totalBgError[im] = sqrt( totalBgError[im] );
    totalSignalMCError[im] = sqrt( totalSignalMCError[im] );
  } // end loop over mass bins

  printf("Printout of the data, MC signal and MC backgrounds integrated over rapidity\n");
  printf("     mass bin        data      MC-signal     MC-backgr\n");
  for(int im = 0; im < DYTools::nMassBins2D; im++){
    printf("%5.0f-%5.0f GeV: ", DYTools::massBinLimits2D[im],
	   DYTools::massBinLimits2D[im+1]);
    printf(" %7.0f", totalData[im]);
    printf(" %9.1f+-%6.1f", totalSignalMC[im], totalSignalMCError[im]);
    printf(" %7.2f+-%5.2f", totalBg[im], totalBgError[im]);
    printf("\n");
  }

  gBenchmark->Show("prepareYields");      
}

void drawRapidityInMassSlices(vector<TMatrixD*> &yields, 
			      vector<TMatrixD*> &yieldsSumw2,
			      vector<TString>   &sampleTags,
			      int nhist, TCanvas *canvas){

  if(nhist != 7){
    printf("ERROR: this function can draw only 6 histograms on one canvas\n");
    printf("       after underflow is dropped\n");
    printf("       Drawing is aborted\n");
    return;
  }

  // Create rapidity histograms for data and total MC,
  // drop underflow bin.
  TH1F *histData[6];
  TH1F *histTotalMC[6];
  char hname[100];
  for(int im=0; im<6; im++){ // Loop over mass slices
    // Start from 1, ignoring 0=underflow
    int iMassBin = im + 1;
    sprintf(hname,"rapidity_data_for_massbin_%i",iMassBin);
    histData[im] = new TH1F(hname,"",nYBins[iMassBin], DYTools::yRangeMin, DYTools::yRangeMax);
    sprintf(hname,"rapidity_mc_for_massbin_%i",iMassBin);
    histTotalMC[im] = new TH1F(hname,"",nYBins[iMassBin], DYTools::yRangeMin, DYTools::yRangeMax);
    for(UInt_t isam = 0; isam < yields.size(); isam++){ // Loop over samples
      TMatrixD *thisSampleYields      = yields.at(isam);
      TMatrixD *thisSampleYieldsSumw2 = yieldsSumw2.at(isam);
      for(int iy = 0; iy < DYTools::nYBins[iMassBin]; iy++){ // Loop over rapidity bins
	// In the code below, we are being careful with indices:
	//   the content bins in the histogram start from 1, while the rapidity
	//   bin index in the yields arrays start from zero.
	if(sampleTags.at(isam) == TString("data") ){
	  histData[im]->SetBinContent(iy+1, 
				      histData[im]->GetBinContent(iy+1)
				      + (*thisSampleYields)(iMassBin,iy) );
	  double oldError = histData[im]->GetBinError(iy+1);
	  histData[im]->SetBinError(iy+1,
				   sqrt( oldError*oldError 
					 + (*thisSampleYieldsSumw2)(iMassBin,iy)) );
	} else {
	  histTotalMC[im]->SetBinContent(iy+1, 
					histTotalMC[im]->GetBinContent(iy+1)
					+ (*thisSampleYields)(iMassBin,iy) );
	  double oldError = histTotalMC[im]->GetBinError(iy+1);
	  histTotalMC[im]->SetBinError(iy+1,
				      sqrt( oldError*oldError 
					    + (*thisSampleYieldsSumw2)(iMassBin,iy)) );
	}
      } // End loop over rapidity bins
    } // End loop over samples
  } // End loop over mass slices

  canvas->cd();
  
  double topCanvasMargin = 0.1;
  double bottomCanvasMargin = 0.1;
  const int nPadRows = 3;
  const int nPadCols = 2;
  double plotsize = (1.0 - topCanvasMargin - bottomCanvasMargin)/nPadRows;

  TPad *pads[6];
  for(int irow = 0; irow < nPadRows; irow++){
    for(int icol = 0; icol < nPadCols; icol++){
      int padindex = icol*nPadRows + irow;
      TString padname = "pad_";
      padname += padindex;
      double xlow = icol*0.5;
      double xhigh = (icol+1)*0.5;
      double ylow = bottomCanvasMargin + plotsize * irow;
      if( irow == 0 ) 
	ylow = 0;
      double yhigh = bottomCanvasMargin + plotsize * (irow + 1);
      if( irow == nPadRows - 1)
	yhigh = 1.0;
      pads[padindex] = new TPad(padname,padname,xlow,ylow,xhigh,yhigh);
      if(irow != 0)
	pads[padindex]->SetBottomMargin(0);
      else
	pads[padindex]->SetBottomMargin(bottomCanvasMargin/(plotsize+bottomCanvasMargin));
      if(irow != nPadRows - 1)
	pads[padindex]->SetTopMargin(0);
      else
	pads[padindex]->SetTopMargin(topCanvasMargin/(plotsize+topCanvasMargin));
    }
  }

  printf("Draw slices\n");
  for(int irow = 0; irow < nPadRows; irow++){
    for(int icol = 0; icol < nPadCols; icol++){
      int padindex = icol*nPadRows + irow;
      canvas->cd();
      pads[padindex]->Draw();
      pads[padindex]->cd();
      //
      double lsize = 0.1; //font size
      if(irow != 0 && irow != nPadRows - 1)
 	lsize = lsize * (bottomCanvasMargin+plotsize)/plotsize;
      histData[padindex]->GetYaxis()->SetLabelSize(lsize);
      histData[padindex]->GetYaxis()->SetNdivisions(505);
      histTotalMC[padindex]->SetLineWidth(2);
      // 
      if(irow == 0){
	histData[padindex]->GetXaxis()->SetLabelSize(lsize);
	histData[padindex]->GetXaxis()->SetTitle("rapidity");
	histData[padindex]->GetXaxis()->SetTitleSize(0.15);
	histData[padindex]->GetXaxis()->SetTitleOffset(0.8);
      }      
      histData[padindex]->SetMarkerSize(0.8);
      histData[padindex]->Draw("pe");
      histTotalMC[padindex]->Draw("hist,same");
    }
  }

}
