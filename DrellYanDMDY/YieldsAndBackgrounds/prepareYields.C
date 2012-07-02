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
#include <THStack.h>
#include <TBenchmark.h>             // class to track macro running statistics
#include <TVector3.h>               // 3D vector class
#include <TArrayD.h>
#include <TVectorD.h>
#include <TMatrixD.h>
#include <TLorentzVector.h>         // 4-vector class
#include <TRandom.h>
#include <TDatime.h>                // time stamp

#include <TLatex.h> 

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
#include "../Include/DYToolsUI.hh"
#include "../Include/plotFunctions.hh"
#include "../Include/PUReweight.hh"
#include "../Include/UnfoldingTools.hh"

#endif

// Forward declarations
/*
void DrawMassPeak(vector<TH1F*> hMassv, vector<CSample*> samplev, vector<TString> snamev, TH1F* hMassDibosons, bool hasData, 
                   bool mergeDibosons, TString labelDibosons, Int_t colorDibosons, Double_t lumi, char* lumitext,  bool actualBinning);

void DrawFlattened(vector<TMatrixD*> yields, vector<TMatrixD*> yieldsSumw2, vector<CSample*> samplev, vector<TString> snamev, bool hasData, 
                   bool mergeDibosons, TString labelDibosons, Int_t colorDibosons, Double_t lumi, char* lumitext);

void Draw6Canvases(vector<TMatrixD*> yields, vector<TMatrixD*> yieldsSumw2,
                    vector<CSample*> samplev, vector<TString> snamev, 
                    bool hasData, double dataOverMc, double* dataOverMcEachBin, bool normEachBin=1, bool singleCanvas=0);
void SomeHistAttributes (TH1F* hist, TString samplename);
*/
//void SaveCanvas(TCanvas* canv, TString canvName);

//=== MAIN MACRO =================================================================================================

void prepareYields(const TString conf  = "data_plot.conf",
		   DYTools::TSystematicsStudy_t runMode=DYTools::NORMAL,
		   const TString plotsDirExtraTag="")
{  
  gBenchmark->Start("prepareYields");

  std::cout << "\n\nRun mode: " << SystematicsStudyName(runMode) << "\n";
  switch(runMode) {
  case DYTools::NORMAL:
  case DYTools::ESCALE_STUDY:
  case DYTools::ESCALE_STUDY_RND:
    break;
  default:
    std::cout << "prepareYields is not ready for runMode=" << SystematicsStudyName(runMode) << "\n";
    throw 2;
  }
 
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
  
  // 
  // Set up energy scale corrections
  //
  ElectronEnergyScale escale(escaleTag);
  assert(escale.isInitialized());
  escale.print();

  // sOutDir is a static data member in the CPlot class.
  // There is a strange crash of the whole ROOT session well after
  // this script is executed when one attempts to exit ROOT, with 
  // a dump of memory map. This happens only on UNL Tier3, but
  // there is no sign of a problem on any other computer.
  //   The consequence of this variable is not set is that the plots
  // will be created in the local directory rather than the
  // one configured through sOutDir.
//   CPlot::sOutDir = outputDir + TString("/plots");

  if ((runMode==DYTools::ESCALE_STUDY) || (runMode==DYTools::ESCALE_STUDY_RND)) {
    CPlot::sOutDir = "plots_escale/";
  }
  else CPlot::sOutDir = "plots";
  CPlot::sOutDir += plotsDirExtraTag;

  Bool_t hasData = (samplev[0]->fnamev.size()>0);
    
  //
  // Canvas dimensions
  //
  //Int_t canw=800, canh=600;
  
  
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
    TMatrixD *yieldsMatrix = new TMatrixD(DYTools::nMassBins, maxYBins);
    (*yieldsMatrix) = 0;
    yields.push_back(yieldsMatrix);
    //
    TMatrixD *yieldsSumw2Matrix = new TMatrixD(DYTools::nMassBins, maxYBins);
    (*yieldsSumw2Matrix) = 0;
    yieldsSumw2.push_back(yieldsSumw2Matrix);
    //
    sampleTags.push_back(snamev[isam]);

    sprintf(hname,"hMass_%i",isam);  
    hMassv.push_back(new TH1F(hname,"",1490,10,1500));  
    hMassv[isam]->Sumw2();
    sprintf(hname,"hMassBins_%i",isam);  
    hMassBinsv.push_back(new TH1F(hname,"",
				  DYTools::nMassBins,
				  DYTools::massBinLimits));
    hMassBinsv[isam]->Sumw2();
			 
    nSelv.push_back(0);
    nSelVarv.push_back(0);  
  }
  

#ifdef ZeeData_is_TObject
  ZeeData_t *data = new ZeeData_t();
#else
  ZeeData *data = new ZeeData();
#endif
  TRandom random;

  int puReweight_new_code=1;
  // Open file with number of PV distributions for pile-up reweighting
  // TString("/npv.root");
  PUReweight_t puWeight;
  const TString fnamePV = outputDir+TString("/npv") + analysisTag_USER + TString(".root");
  std::cout << "fnamePV=<" << fnamePV << ">\n";
  TFile *pvfile = NULL;
  TH1F *hPVData = 0;
  if (puReweight_new_code) {
    assert(puWeight.setFile(fnamePV));
    assert(puWeight.setReference("hNGoodPV_data"));
  }
  else {
    pvfile=new TFile(fnamePV);
    assert(pvfile);
    if(hasData){ 
      hPVData = (TH1F*)pvfile->Get("hNGoodPV_data"); assert(hPVData); 
    }
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

    TString fname = outputDir + TString("/ntuples/") + snamev[isam] + analysisTag_USER + TString("_select.root");
    if ((isam==0) && 
	((runMode==DYTools::ESCALE_STUDY) || (runMode==DYTools::ESCALE_STUDY_RND))) {
      // selectEvents corrects only data energies
      TString fnameTag=TString("_select_") + escale.calibrationSetShortName();
      fname.Replace(fname.Index("_select."),sizeof("_select.")-2,fnameTag);
    }
    cout << "Processing " << fname << "..." << endl;   
    infile = new TFile(fname);
    assert(infile); 

    // Prepare weights for pile-up reweighting for MC
    TH1F *puWeights=NULL;
    if (puReweight_new_code) {
      assert(puWeight.setActiveSample(TString("hNGoodPV_")+snamev[isam]));
      //puWeight.printActiveDistr_and_Weights(std::cout);
      //puWeight.printWeights(std::cout);
    }
    else {
      TH1F *hPVThis = (TH1F*) pvfile->Get(TString("hNGoodPV_")+snamev[isam]); assert(hPVThis);
      hPVThis->SetDirectory(0);
      // Normalize or not? Not clear
      hPVThis->Scale( hPVData->GetSumOfWeights()/hPVThis->GetSumOfWeights());
      puWeights = (TH1F*)hPVData->Clone("puWeights");
      puWeights->SetDirectory(0);
      puWeights->Divide(hPVThis);
      for(int i=1; i<=puWeights->GetNbinsX(); i++)
	printf(" %f    %f    %f\n",puWeights->GetBinCenter(i),puWeights->GetBinContent(i),puWeights->GetBinError(i));
    }

    // Get the TTree and set branch address
    eventTree = (TTree*)infile->Get("Events"); assert(eventTree); 
    eventTree->SetBranchAddress("Events",&data);


    TMatrixD *thisSampleYields = yields.at(isam);
    TMatrixD *thisSampleYieldsSumw2 = yieldsSumw2.at(isam);

    std::cout << "here are " << eventTree->GetEntries() << " entries in " << snamev[isam] << " sample\n";
    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
      eventTree->GetEntry(ientry);
      Double_t weight = data->weight;
      
      // Any extra weight factors:
      double weightPU = (puReweight_new_code) ?
	puWeight.getWeight( data->nPV ) :
	puWeights->GetBinContent( puWeights->FindBin( data->nPV ));
      weight *= weightPU;
      
      // If This is MC, add extra smearing to the mass
      // We apply extra smearing to all MC samples: it is may be
      // not quite right for fake electron backgrounds, but these
      // are not dominant, and in any case we do not have corrections
      // for fake electrons.
      if(isam!=0) {            	    
	double smearingCorrection = escale.generateMCSmear(data->scEta_1, data->scEta_2);
	data->mass = data->mass + smearingCorrection;
      }

      // Find the 2D bin for this event:
      int massBin = findMassBin(data->mass);
      int yBin    = findAbsYBin(massBin, data->y);

      if ((massBin==-1) || (yBin==-1)) // out of range
	continue;

      (*thisSampleYields)(massBin,yBin) += weight;
      (*thisSampleYieldsSumw2)(massBin,yBin) += weight*weight;

      hMassv[isam]->Fill(data->mass,weight);
      hMassBinsv[isam]->Fill(data->mass,weight);

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
  //char ylabel[100];   // y-axis label
  char lumitext[50];
  if(lumi>0) {
    if(lumi<1) { sprintf(lumitext,"#int#font[12]{L}dt = %.0f nb^{-1}",1000.*lumi); }
    if(lumi<1000) { sprintf(lumitext,"#int#font[12]{L}dt = %.0f nb^{-1}",1000.*lumi); }
    else       { sprintf(lumitext,"#int#font[12]{L}dt = %.3g fb^{-1}",lumi/1000.0); }
  }

  // Merge diboson histograms if needed
  TH1F *hMassBinsDibosons = (TH1F*)hMassBinsv[1]->Clone("hMassBinsDibosons");
  TH1F *hMassDibosons = (TH1F*)hMassv[1]->Clone("hMassDibosons");
  hMassBinsDibosons->Reset();
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

  // Additional normalization for MC:
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

  double dataOverMcEachBin[nMassBins+1];
  for (int i=0; i<nMassBins; i++)
    {
      dataOverMcEachBin[i] = hMassv[0]->Integral(hMassv[0]->FindBin(massBinLimits[i]+0.001),hMassv[0]->FindBin(massBinLimits[i+1]-0.001)) /
      totalMCMass->Integral(totalMCMass->FindBin(massBinLimits[i]+0.001),totalMCMass->FindBin(massBinLimits[i+1]-0.001));
      printf("data to MC %i bin norm: %f\n",i,dataOverMcEachBin[i]);
    }

  //std::cout << "\n\nsetting dataOverMc=1\n"; dataOverMc=1; // for 1_to_1 comparison with 1D
  
  // Rescale all MC samples. This is not totally proper for fake lepton
  // backgrounds, but ok for backgrounds with true leptons, and those are dominant
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    hMassv[isam]->Scale(dataOverMc);
    hMassBinsv[isam]->Scale(dataOverMc);
    printf("  MC %s IS RESCALED for plotting\n", snamev[isam].Data());
  }

  //
  // Prepare outputDir and the plot file
  //

  TString outputDirYields(outputDir.Data());
  outputDirYields.ReplaceAll("selected_events","yields");
  gSystem->mkdir(outputDirYields,kTRUE);
  TString fNameOutYieldPlots(outputDirYields+TString("/yield_plots") + analysisTag);
  fNameOutYieldPlots += ".root";
  TFile *fYieldPlots = new TFile( fNameOutYieldPlots, "recreate" );
  if (!fYieldPlots) {
    std::cout << "Failed to create a file <" << fNameOutYieldPlots << ">\n";
    throw 2;
  }

  //
  // Draw mass spectrum without rapidity binning
  //


  // First, draw the mass histograms with fine mass binning
  DrawMassPeak(hMassv, samplev, snamev, hMassDibosons, hasData, mergeDibosons, labelDibosons, colorDibosons, lumi, lumitext, 0, fYieldPlots);

  // Second, draw the mass histograms with the mass binning used in the analysis
  DrawMassPeak(hMassBinsv, samplev, snamev, hMassBinsDibosons, hasData, mergeDibosons, labelDibosons, colorDibosons, lumi, lumitext, 1, fYieldPlots);

  // Draw the flattened figure (Y histograms for different mass regions)
  DrawFlattened(yields, yieldsSumw2, samplev, snamev, hasData, mergeDibosons, labelDibosons, colorDibosons, lumi, lumitext, fYieldPlots);

  // Draw rapidity in mass slices 
  
  Draw6Canvases(yields, yieldsSumw2, samplev, snamev, hasData, dataOverMc, dataOverMcEachBin, 1, 0, fYieldPlots);
  Draw6Canvases(yields, yieldsSumw2, samplev, snamev, hasData, dataOverMc, dataOverMcEachBin, 1, 1, fYieldPlots);

  fYieldPlots->Close();

  //--------------------------------------------------------------------------------------------------------------
  // Save the results
  //==============================================================================================================

  // Pack info into writable objects
  TVectorD massBinLimits(nMassBins+1);
  TVectorD rapidityBinning(nMassBins+1);
  for(int i=0; i <= nMassBins; i++){
    massBinLimits(i) = DYTools::massBinLimits[i];
    rapidityBinning(i) = DYTools::nYBins[i];
  }

  // This dummy object is only needed to convey the number
  // of samples considered. The method below is rather convoluted,
  // but I do not know a better one. Ideally, we would just store
  // a list of strings, each string containing the sample name.
  TVectorD dummySampleCount(sampleTags.size());
  dummySampleCount = 0;

  gSystem->mkdir(outputDirYields,kTRUE);
  TString fNameOutYields(outputDirYields+TString("/yields") + analysisTag);
  fNameOutYields += ".root";
  TFile fYields( fNameOutYields, "recreate" );
  massBinLimits      .Write("massBinLimits");
  rapidityBinning    .Write("rapidityBinning");
  dummySampleCount   .Write("dummySampleCount");
  unfolding::writeBinningArrays(fYields);
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
  
  /*
  // Save mass histograms into a separate file for a direct comparison 
  // with DrellYan(1D) package
  TString fNameOutHists(outputDirYields+"/massHist");
  fNameOutHists += analysisTag + TString(".root");
  TFile fMassHists(fNameOutHists,"recreate");
  for(UInt_t isam=0; isam<samplev.size(); isam++) {
    hMassBinsv[isam]->Write(snamev[isam]);
  }
  std::cout << "file <" << fNameOutHists << "> created\n";
  fMassHists.Close();
  */


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
  double totalData            [DYTools::nMassBins];
  double totalSignalMC        [DYTools::nMassBins];
  double totalSignalMCError   [DYTools::nMassBins];
  double totalBg     [DYTools::nMassBins];
  double totalBgError[DYTools::nMassBins];

  for( int im=0; im<DYTools::nMassBins; im++){
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
  for(int im = 0; im < DYTools::nMassBins; im++){
    printf("%5.0f-%5.0f GeV: ", DYTools::massBinLimits[im],
	   DYTools::massBinLimits[im+1]);
    printf(" %7.0f", totalData[im]);
    printf(" %9.1f+-%6.1f", totalSignalMC[im], totalSignalMCError[im]);
    printf(" %7.2f+-%5.2f", totalBg[im], totalBgError[im]);
    printf("\n");
  }
  printf("Note: these MC numbers are not rescaled!\n");

  if (1) {
  // A different view of background table
  printf("\n\nPrintout of the backgrounds for all mass bins, view II\n");

  printf("            ");
  for(UInt_t isam=0; isam<samplev.size(); isam++) {
    printf(" %14s ",snamev[isam].Data());
  }
  printf("           total          fraction\n");
  for(int ibin=0; ibin<nMassBins; ibin++){
    printf("%5.1f-%5.1f GeV: ",
           hMassBinsv[0]->GetXaxis()->GetBinLowEdge(ibin+1),
           hMassBinsv[0]->GetXaxis()->GetBinUpEdge(ibin+1));
    // Data:
    printf(" %7.0f+-%5.0f ",hMassBinsv[0]->GetBinContent(ibin+1),hMassBinsv[0]->GetBinError(ibin+1) );
    // Individual MC samples
    double total=0., totalError=0.;
    for(UInt_t isam=1; isam<samplev.size(); isam++) {
      double thisContent = hMassBinsv[isam]->GetBinContent(ibin+1);
      double thisError = hMassBinsv[isam]->GetBinError(ibin+1);
      printf(" %7.2f+-%5.2f ",thisContent, thisError);
      if ( (isam!=0) && (snamev[isam]!=TString("zee"))) {
	total+= thisContent;
	totalError+=thisError*thisError;
      }
    }
    totalError = sqrt(totalError);
    // Total
    printf("  %8.2f+-%6.2f",total, totalError);
    printf("  %5.1f\n",100*total/hMassBinsv[0]->GetBinContent(ibin+1));
  }
  }



  gBenchmark->Show("prepareYields");      
}

