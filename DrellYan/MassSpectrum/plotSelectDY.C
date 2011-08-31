//================================================================================================
//
// Plots for reconstructed Z candidates
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

//=== MAIN MACRO =================================================================================================

void plotSelectDY(const TString conf  = "data_plot.conf") 
{  
  gBenchmark->Start("plotSelectDY");

  
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================

  TString  outputDir;         // output directory
  Double_t lumi;              // luminosity (pb^-1)
  Bool_t   doWeight;          // weight events?
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
  
  CPlot::sOutDir = outputDir + TString("/plots");

  Bool_t hasData = (samplev[0]->fnamev.size()>0);
    
  //
  // Canvas dimensions
  //
  Int_t canw=800, canh=600;
  
  
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //============================================================================================================== 
  
  vector<TH1F*>    hMassv;
  vector<TH1F*>    hMassBinsv;
  vector<TH1F*>    hPtv;
  vector<TH1F*>    hyv;
  vector<TH1F*>    hPhiv;   
  vector<TH1F*>    hElePtv;
  vector<TH1F*>    hEleEtav;
  vector<TH1F*>    hElePhiv;
  vector<TH1F*>    hCaloMetv;
  vector<TH1F*>    hCaloSumEtv;
  vector<TH1F*>    hTCMetv;
  vector<TH1F*>    hTCSumEtv;
  vector<TH1F*>    hPFMetv;
  vector<TH1F*>    hPFSumEtv;
  vector<TH1F*>    hNPVv;
  vector<Double_t> nSelv;
  vector<Double_t> nSelVarv;  

  char hname[100];
  for(UInt_t isam=0; isam<samplev.size(); isam++) {
    sprintf(hname,"hMass_%i",isam);  hMassv.push_back(new TH1F(hname,"",790,10,800));  hMassv[isam]->Sumw2();
    sprintf(hname,"hMassBins_%i",isam);  hMassBinsv.push_back(new TH1F(hname,"",nMassBins,massBinLimits));  hMassBinsv[isam]->Sumw2();
    sprintf(hname,"hPt_%i",isam);    hPtv.push_back(new TH1F(hname,"",20,0,100));     hPtv[isam]->Sumw2();
    sprintf(hname,"hy_%i",isam);     hyv.push_back(new TH1F(hname,"",20,-3,3));       hyv[isam]->Sumw2();
    sprintf(hname,"hPhi_%i",isam);   hPhiv.push_back(new TH1F(hname,"",20,-3.2,3.2)); hPhiv[isam]->Sumw2();
    
    sprintf(hname,"hElePt_%i",isam);  hElePtv.push_back(new TH1F(hname,"",30,0,150));     hElePtv[isam]->Sumw2();
    sprintf(hname,"hEleEta_%i",isam); hEleEtav.push_back(new TH1F(hname,"",20,-3,3));     hEleEtav[isam]->Sumw2();
    sprintf(hname,"hElePhi_%i",isam); hElePhiv.push_back(new TH1F(hname,"",20,-3.2,3.2)); hElePhiv[isam]->Sumw2(); 

    sprintf(hname,"hCaloMet_%i",isam);   hCaloMetv.push_back(new TH1F(hname,"",20,0,60));    hCaloMetv[isam]->Sumw2();
    sprintf(hname,"hCaloSumEt_%i",isam); hCaloSumEtv.push_back(new TH1F(hname,"",25,0,500)); hCaloSumEtv[isam]->Sumw2();
    sprintf(hname,"hTCMet_%i",isam);     hTCMetv.push_back(new TH1F(hname,"",20,0,60));      hTCMetv[isam]->Sumw2();
    sprintf(hname,"hTCSumEt_%i",isam);   hTCSumEtv.push_back(new TH1F(hname,"",25,0,500));   hTCSumEtv[isam]->Sumw2();
    sprintf(hname,"hPFMet_%i",isam);     hPFMetv.push_back(new TH1F(hname,"",20,0,60));      hPFMetv[isam]->Sumw2();
    sprintf(hname,"hPFSumEt_%i",isam);   hPFSumEtv.push_back(new TH1F(hname,"",25,0,500));   hPFSumEtv[isam]->Sumw2();  
    
    sprintf(hname,"hNPV_%i",isam); hNPVv.push_back(new TH1F(hname,"",10,-0.5,9.5)); hNPVv[isam]->Sumw2();
      
    nSelv.push_back(0);
    nSelVarv.push_back(0);  
  }
  
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
    for(int i=1; i<=puWeights->GetNbinsX(); i++)
      printf(" %f    %f    %f\n",puWeights->GetBinCenter(i),puWeights->GetBinContent(i),puWeights->GetBinError(i));

    // Get the TTree and set branch address
    eventTree = (TTree*)infile->Get("Events"); assert(eventTree); 
    eventTree->SetBranchAddress("Events",&data);

    // DEBUG
    ofstream out;
    if(isam==0){
      out.open("data_event_dump.txt");
      assert(out.is_open());
    }

    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
      eventTree->GetEntry(ientry);
      Double_t weight = data.weight;

      // Any extra weight factors:
      double weightPU = puWeights->GetBinContent( puWeights->FindBin( data.nPV ));
      weight *= weightPU;
      
      // Apply energy scale

      Double_t scEt1 = data.scEt_1;
      Double_t scEt2 = data.scEt_2;
      // Electron energy scale correction
      if(isam==0) {            	    

	// WARNING: Energy scale should have been applied by this point
	// during the original selection. The code below is left in case
	// if related studies need to be done. But the correction really
	// needs to be applied in the original selection script so that 
	// ET cut on electrons is applied properly there.
	//   	double corr1 = findEnergyScaleCorrection(data.scEta_1);
	//   	double corr2 = findEnergyScaleCorrection(data.scEta_2);
 	double corr1 = 1;
 	double corr2 = 1;
	scEt1 = data.scEt_1 * corr1;
	scEt2 = data.scEt_2 * corr2;
	
	TLorentzVector ele1; 
	ele1.SetPtEtaPhiM(data.pt_1,data.eta_1,data.phi_1,0.000511);
	ele1 *= corr1;
	data.pt_1  = ele1.Pt();
	data.eta_1 = ele1.Eta();
	data.phi_1 = ele1.Phi();
        
	TLorentzVector ele2; 
	ele2.SetPtEtaPhiM(data.pt_2,data.eta_2,data.phi_2,0.000511);
	ele2 *= corr2;
	data.pt_2  = ele2.Pt();
	data.eta_2 = ele2.Eta();
	data.phi_2 = ele2.Phi();
	
	TLorentzVector vDiEle = ele1+ele2;            
	data.mass = vDiEle.M();
	data.pt   = vDiEle.Pt();
	data.y    = vDiEle.Rapidity();
	data.phi  = vDiEle.Phi(); 
      }

      // If This is MC, add extra smearing to the mass
      if(isam!=0) {            	    
	double smear1 = extraSmearingSigma(data.scEta_1);
	double smear2 = extraSmearingSigma(data.scEta_2);
	double smearTotal = sqrt(smear1*smear1 + smear2*smear2);
	data.mass = data.mass + random.Gaus(0.0,smearTotal);
      }

      hMassv[isam]->Fill(data.mass,weight);
      hMassBinsv[isam]->Fill(data.mass,weight);
      // DEBUG in case event list needs to be printed
//       if(isam==0){
// 	char line[500];
// 	sprintf(line,"%4i %4i %16i %f %f %f %f %f\n",
// 		data.runNum, data.lumiSec, data.evtNum,
// 		data.mass, data.pt_1, data.pt_2, data.scEt_1, data.scEt_2);
// 	out << line;
//       }

      hPtv[isam]  ->Fill(data.pt,  weight);
      hyv[isam]   ->Fill(data.y,   weight);
      hPhiv[isam] ->Fill(data.phi, weight);
      
      hElePtv[isam] ->Fill(data.pt_1, weight); hElePtv[isam] ->Fill(data.pt_2, weight);
      hEleEtav[isam]->Fill(data.eta_1,weight); hEleEtav[isam]->Fill(data.eta_2,weight);
      hElePhiv[isam]->Fill(data.phi_1,weight); hElePhiv[isam]->Fill(data.phi_2,weight);

      hNPVv[isam]->Fill(data.nPV,weight);
	  
      nSelv[isam] += weight;
      nSelVarv[isam] += weight*weight;
      
      TVector3 met;
      
      if(data.caloMEx!=0 || data.caloMEy!=0) {       
        met.SetXYZ(data.caloMEx, data.caloMEy, 0);
        hCaloMetv[isam]  ->Fill(met.Perp(),         weight);
        hCaloSumEtv[isam]->Fill(data.caloSumET,weight);
      }
      
      if(data.tcMEx!=0 || data.tcMEy!=0) {       
        met.SetXYZ(data.tcMEx, data.tcMEy, 0);
        hTCMetv[isam]  ->Fill(met.Perp(),       weight);
        hTCSumEtv[isam]->Fill(data.tcSumET,weight);
      }
      
      if(data.pfMEx!=0 || data.pfMEy!=0) {       
        met.SetXYZ(data.pfMEx, data.pfMEy, 0);
        hPFMetv[isam]  ->Fill(met.Perp(),       weight);
        hPFSumEtv[isam]->Fill(data.pfSumET,weight);
      }
    }
    delete infile;
    infile=0, eventTree=0;
    if(isam==0)
      out.close();
  }

  //--------------------------------------------------------------------------------------------------------------
  // Make plots
  //==============================================================================================================
  TCanvas *c = MakeCanvas("c","c",canw,canh);
  
  // string buffers
  char ylabel[100];   // y-axis label
  char lumitext[50];
  if(lumi>0) {
    if(lumi<1) { sprintf(lumitext,"#int#font[12]{L}dt = %.0f nb^{-1}",1000.*lumi); }
    else       { sprintf(lumitext,"#int#font[12]{L}dt = %.3g pb^{-1}",lumi); }
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


//   char text[100];  
//   sprintf(text,"%i Events",(Int_t)(hMassv[0]->Integral()));
//   sprintf(ylabel,"a.u. / %.1f GeV/c^{2}",hMassv[0]->GetBinWidth(1));
//   CPlot plotMass("mass0","3_6_X","m(e^{+}e^{-}) [GeV/c^{2}]",ylabel);
//   if(hasData) { 
//     // hMassv[0]->Scale(1.0/(hMassv[0]->Integral()));
//     plotMass.AddHist1D(hMassv[0],samplev[0]->label,"E"); 
//   }
//   for(UInt_t isam=1; isam<samplev.size(); isam++) {
//     // hMassv[isam]->Scale(1.0/(hMassv[isam]->Integral()));
//     plotMass.AddToStack(hMassv[isam],samplev[isam]->label,samplev[isam]->color);
//   }
//   //plotMass.SetYRange(0,0.2);
//   plotMass.TransLegend(0.1,0);
//   plotMass.AddTextBox(text,0.21,0.85,0.41,0.8,0);
//   sprintf(text,"p_{KS} = %.4f",hMassv[0]->KolmogorovTest(hMassv[1]));
//   plotMass.AddTextBox(text,0.21,0.78,0.41,0.73,0);
//   plotMass.Draw(c,kTRUE,format);


  // scaling of MC to have the same area as data in 60-120
//   for(UInt_t isam=1; isam<samplev.size(); isam++){
//     if(snamev[isam] == "zee" ){
//       double scale = 1/1.12;
//       cout << "Scaling signal MC by " << scale << endl;
//       hMassv[isam]->Scale(scale);
//     }else{
//       cout << "sample " << snamev[isam] << " is not scaled" << endl;
//     }
//   }

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
    printf("  MC %s IS RESCALED\n", snamev[isam].Data());
  }

  // dielectron mass
  double massMinPlot = 20;
  double massMaxPlot = 600;
  
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
  plotMass.SetYRange(1.0,35000);  
  plotMass.SetXRange(massMinPlot,massMaxPlot);  
  plotMass.Draw(c,kTRUE,format);
  
  TCanvas *c0 = MakeCanvas("c0","c0",canw,canh);
  plotMass.SetName("masslog");
  plotMass.SetLogy();
  plotMass.SetLogx();
  plotMass.Draw(c0);
  plotMass.SetYRange(1.0,60e3);  
  plotMass.SetXRange(massMinPlot,massMaxPlot);  
  plotMass.Draw(c0,kTRUE,format);

  TCanvas *c2 = MakeCanvas("c2","c2",canw,canh);
  // dielectron mass
//   sprintf(ylabel,"Events / %.1f GeV/c^{2}",hMassBinsv[0]->GetBinWidth(1));
  sprintf(ylabel,"Events");
  CPlot plotMassBins("massBins","","m(e^{+}e^{-}) [GeV/c^{2}]",ylabel);
  plotMassBins.SetLogx();
  if(hasData) { plotMassBins.AddHist1D(hMassBinsv[0],samplev[0]->label,"E"); }
  // Do not draw separately dibosons, but draw the merged histogram if needed
  if(mergeDibosons)
    plotMassBins.AddToStack(hMassBinsDibosons,labelDibosons,colorDibosons);
  for(UInt_t isam=1; isam<samplev.size(); isam++){
    if( !(mergeDibosons && (snamev[isam]=="ww" || snamev[isam]=="wz" || snamev[isam]=="zz")))
      plotMassBins.AddToStack(hMassBinsv[isam],samplev[isam]->label,samplev[isam]->color);
  }
  if(hasData){
    hMassBinsv[0]->GetXaxis()->SetMoreLogLabels();
    hMassBinsv[0]->GetXaxis()->SetNoExponent();
  }
  plotMassBins.SetXRange(massMinPlot,massMaxPlot);  
  plotMassBins.SetYRange(1e-1,300000);  
  plotMassBins.SetLegend(0.75,0.55,0.98,0.9);
  if(lumi>0) plotMassBins.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
  plotMassBins.Draw(c2,kTRUE,format);
  
  TCanvas *c20 = MakeCanvas("c20","c20",canw,canh);
  plotMassBins.SetName("masslogBins");
  plotMassBins.SetLogy();
  plotMassBins.SetLogx();
  plotMassBins.SetXRange(massMinPlot,massMaxPlot);  
  plotMassBins.SetYRange(10,30e4);  
  plotMassBins.Draw(c20,kTRUE,format);




  /*    
  // dielectron pT
  sprintf(ylabel,"Events / %.1f GeV/c",hPtv[0]->GetBinWidth(1));
  CPlot plotPt("pt","","p_{T}(e^{+}e^{-}) [GeV/c]",ylabel);
  if(hasData) { plotPt.AddHist1D(hPtv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++)
    plotPt.AddToStack(hPtv[isam],samplev[isam]->label,samplev[isam]->color);
  if(lumi>0) plotPt.AddTextBox(lumitext,0.41,0.85,0.61,0.8,0);
  plotPt.SetLegend(0.75,0.55,0.98,0.9);
  plotPt.Draw(c,kTRUE,format);
  
  plotPt.SetName("ptlog");  
  plotPt.SetLogy();
  plotPt.SetYRange(1e-2,1e3);  
  plotPt.Draw(c,kTRUE,format);
    
  // dielectron eta
  sprintf(ylabel,"Events / %.2f",hyv[0]->GetBinWidth(1));
  CPlot ploty("y","","y(e^{+}e^{-})",ylabel);
  if(hasData) { ploty.AddHist1D(hyv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++)
    ploty.AddToStack(hyv[isam],samplev[isam]->label,samplev[isam]->color);
  if(lumi>0) ploty.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
  ploty.SetLegend(0.75,0.55,0.98,0.9);
  ploty.SetYRange(0,1.5*(ploty.GetStack()->GetMaximum()));
  ploty.Draw(c,kTRUE,format);
    
  // dielectron phi
  sprintf(ylabel,"Events / %.2f",hPhiv[0]->GetBinWidth(1));
  CPlot plotPhi("phi","","#phi(e^{+}e^{-})",ylabel);
  if(hasData) { plotPhi.AddHist1D(hPhiv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++)
    plotPhi.AddToStack(hPhiv[isam],samplev[isam]->label,samplev[isam]->color);
  if(lumi>0) plotPhi.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
  plotPhi.SetYRange(0,1.5*(plotPhi.GetStack()->GetMaximum()));
  plotPhi.SetLegend(0.75,0.55,0.98,0.9);
  plotPhi.Draw(c,kTRUE,format);
    
  // electron pT
  sprintf(ylabel,"Events / %.1f GeV/c",hElePtv[0]->GetBinWidth(1));
  CPlot plotElePt("elept","","p_{T}(e) [GeV/c]",ylabel);
  if(hasData) { plotElePt.AddHist1D(hElePtv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++)
    plotElePt.AddToStack(hElePtv[isam],samplev[isam]->label,samplev[isam]->color);
  if(lumi>0) plotElePt.AddTextBox(lumitext,0.45,0.85,0.65,0.8,0);
  plotElePt.SetLegend(0.75,0.55,0.98,0.9);
  plotElePt.Draw(c,kTRUE,format);
  
  plotElePt.SetName("eleptlog");  
  plotElePt.SetLogy();
  plotElePt.SetYRange(1e-2,1e3);  
  plotElePt.Draw(c,kTRUE,format);
    
  // electron eta
  sprintf(ylabel,"Events / %.2f",hEleEtav[0]->GetBinWidth(1));
  CPlot plotEleEta("eleeta","","#eta(e) [GeV/c]",ylabel);
  if(hasData) { plotEleEta.AddHist1D(hEleEtav[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++)
    plotEleEta.AddToStack(hEleEtav[isam],samplev[isam]->label,samplev[isam]->color);
  if(lumi>0) plotEleEta.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
  plotEleEta.SetLegend(0.75,0.55,0.98,0.9);
  plotEleEta.SetYRange(0,1.5*(plotEleEta.GetStack()->GetMaximum()));
  plotEleEta.Draw(c,kTRUE,format);
    
  // electron phi
  sprintf(ylabel,"Events / %.2f",hElePhiv[0]->GetBinWidth(1));
  CPlot plotElePhi("elephi","","#phi(e) [GeV/c]",ylabel);
  if(hasData) { plotElePhi.AddHist1D(hElePhiv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++)
    plotElePhi.AddToStack(hElePhiv[isam],samplev[isam]->label,samplev[isam]->color);
  if(lumi>0) plotElePhi.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
  plotElePhi.SetYRange(0,1.5*(plotElePhi.GetStack()->GetMaximum()));
  plotElePhi.Draw(c,kTRUE,format);

  //
  // MET and SumET
  // 
  sprintf(ylabel,"Events / %.1f GeV",hCaloMetv[0]->GetBinWidth(1));
  CPlot plotCaloMet("calomet","","Calo #slash{E}_{T} [GeV]",ylabel);
  if(hasData) { plotCaloMet.AddHist1D(hCaloMetv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++)
    plotCaloMet.AddToStack(hCaloMetv[isam],samplev[isam]->label,samplev[isam]->color);
  if(lumi>0) plotCaloMet.AddTextBox(lumitext,0.35,0.85,0.55,0.8,0);
  plotCaloMet.TransLegend(0.05,0);
  plotCaloMet.Draw(c,kTRUE,format);
  
  sprintf(ylabel,"Events / %.1f GeV/c",hCaloSumEtv[0]->GetBinWidth(1));
  CPlot plotCaloSumEt("calosumet","","Calo #Sigma_{}E_{T} [GeV]",ylabel);
  if(hasData) { plotCaloSumEt.AddHist1D(hCaloSumEtv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++)
    plotCaloSumEt.AddToStack(hCaloSumEtv[isam],samplev[isam]->label,samplev[isam]->color);
  if(lumi>0) plotCaloSumEt.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
  plotCaloSumEt.SetLegend(0.7,0.55,0.9,0.9);
  plotCaloSumEt.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events / %.1f GeV",hTCMetv[0]->GetBinWidth(1));
  CPlot plotTCMet("tcmet","","TC #slash{E}_{T} [GeV]",ylabel);
  if(hasData) { plotTCMet.AddHist1D(hTCMetv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++)
    plotTCMet.AddToStack(hTCMetv[isam],samplev[isam]->label,samplev[isam]->color);
  if(lumi>0) plotTCMet.AddTextBox(lumitext,0.35,0.85,0.55,0.8,0);
  plotTCMet.TransLegend(0.05,0);
  plotTCMet.Draw(c,kTRUE,format);
  
  sprintf(ylabel,"Events / %.1f GeV/c",hTCSumEtv[0]->GetBinWidth(1));
  CPlot plotTCSumEt("tcsumet","","TC #Sigma_{}E_{T} [GeV]",ylabel);
  if(hasData) { plotTCSumEt.AddHist1D(hTCSumEtv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++)
    plotTCSumEt.AddToStack(hTCSumEtv[isam],samplev[isam]->label,samplev[isam]->color);
  if(lumi>0) plotTCSumEt.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
  plotTCSumEt.SetLegend(0.7,0.55,0.9,0.9);
  plotTCSumEt.Draw(c,kTRUE,format);
  
  sprintf(ylabel,"Events / %.1f GeV",hPFMetv[0]->GetBinWidth(1));
  CPlot plotPFMet("pfmet","","PF #slash{E}_{T} [GeV]",ylabel);
  if(hasData) { plotPFMet.AddHist1D(hPFMetv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++)
    plotPFMet.AddToStack(hPFMetv[isam],samplev[isam]->label,samplev[isam]->color);
  if(lumi>0) plotPFMet.AddTextBox(lumitext,0.35,0.85,0.55,0.8,0);
  plotPFMet.TransLegend(0.05,0);
  plotPFMet.Draw(c,kTRUE,format);
  
  sprintf(ylabel,"Events / %.1f GeV/c",hPFSumEtv[0]->GetBinWidth(1));
  CPlot plotPFSumEt("pfsumet","","PF #Sigma_{}E_{T} [GeV]",ylabel);
  if(hasData) { plotPFSumEt.AddHist1D(hPFSumEtv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++)
    plotPFSumEt.AddToStack(hPFSumEtv[isam],samplev[isam]->label,samplev[isam]->color);
  if(lumi>0) plotPFSumEt.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
  plotPFSumEt.SetLegend(0.7,0.55,0.9,0.9);
  plotPFSumEt.Draw(c,kTRUE,format);

  //
  // Primary Vertex Multiplicity
  //
  CPlot plotNPV("npv","","N_{PV}","Events");
  if(hasData) { plotNPV.AddHist1D(hNPVv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++)
    plotNPV.AddToStack(hNPVv[isam],samplev[isam]->label,samplev[isam]->color);
  if(lumi>0) plotNPV.AddTextBox(lumitext,0.41,0.85,0.61,0.8,0);
  plotNPV.TransLegend(0.08,0);
  plotNPV.SetLogy();
  plotNPV.Draw(c,kTRUE,format);
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
    cout << "   Data: " << setprecision(1) << fixed << nSelv[0] << " Z events!" << endl;
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

  // Print backgrounds table for all mass bins
  double totalBg[nMassBins];
  double totalBgError[nMassBins];
  printf("Printout of the backgrounds for all mass bins\n");
  for(int ibin=0; ibin<nMassBins; ibin++){
    printf("%5.1f-%5.1f GeV: ",
	   hMassBinsv[0]->GetXaxis()->GetBinLowEdge(ibin+1),
	   hMassBinsv[0]->GetXaxis()->GetBinUpEdge(ibin+1));

    double total = 0, totalError = 0;
    for(UInt_t isam=0; isam<samplev.size(); isam++) {
      double thisContent = hMassBinsv[isam]->GetBinContent(ibin+1);
      double thisError = hMassBinsv[isam]->GetBinError(ibin+1);
      printf(" %s %7.3f+-%6.3f ",snamev[isam].Data(), thisContent, thisError);
      if( isam != 0 && snamev[isam] != TString("zee") ){
	total += thisContent;
	totalError += thisError*thisError;
      }
    }
    totalError = sqrt(totalError);
    printf("  total = %f+-%f\n",total, totalError);
    totalBg[ibin] = total;
    totalBgError[ibin] = totalError;
  } 
  // A different view of background table
  printf("Printout of the backgrounds for all mass bins, view II\n");

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
    for(UInt_t isam=1; isam<samplev.size(); isam++) {
      double thisContent = hMassBinsv[isam]->GetBinContent(ibin+1);
      double thisError = hMassBinsv[isam]->GetBinError(ibin+1);
      printf(" %7.2f+-%5.2f ",thisContent, thisError);
    }
    // Total
    printf("  %8.2f+-%6.2f",totalBg[ibin], totalBgError[ibin]);
    printf("  %5.1f\n",100*totalBg[ibin]/hMassBinsv[0]->GetBinContent(ibin+1));
  } 

//   double rho = 0.9177; // from 10/264 v3 with WP95
//   printf("Using DUMMY APRIME, do not trust final x-sec values\n");
//   double aprime[nMassBins] = {
//     0.020408,
//     0.005599,
//     0.050018,
//     0.173977,
//     0.257167,
//     0.313543,
//     0.364771,
//     0.386062,
//     0.402452,
//     0.408596,
//     0.417387,
//     0.431701,
//     0.554400}; 

  printf("Printout of the signal and total backgorund, and x-sec for mass bins\n");
  // Also save into array to store
  TVectorD YieldsSignal      (nMassBins);
  TVectorD YieldsSignalErr   (nMassBins);
  TVectorD BinLimitsForYields(nMassBins+1);
  for(int ibin=0; ibin<=nMassBins; ibin++){
    BinLimitsForYields[ibin] = massBinLimits[ibin];
  }
  printf("      mass              selected     background    backgr. fraction, %%\n");
  for(int ibin=0; ibin<nMassBins; ibin++){
    printf("%5.1f-%5.1f GeV: ",
	   hMassBinsv[0]->GetXaxis()->GetBinLowEdge(ibin+1),
	   hMassBinsv[0]->GetXaxis()->GetBinUpEdge(ibin+1));

    double selected      = hMassBinsv[0]->GetBinContent(ibin+1);
    double selectedError = hMassBinsv[0]->GetBinError(ibin+1);

    double signal = selected - totalBg[ibin];
    double signalError = sqrt(selectedError*selectedError +
			      totalBgError[ibin]*totalBgError[ibin]);

    printf(" %7.0f+-%5.0f   %7.1f+-%4.1f  %4.1f\n",
	   selected, selectedError,
	   totalBg[ibin], totalBgError[ibin],
	   totalBg[ibin]/selected*100);

    YieldsSignal   [ibin] = signal;
    YieldsSignalErr[ibin] = signalError;
  } 
  
  // Save yields information into a root file
  TDatime stamp;
  TString outputDirYields(outputDir.Data());
  outputDirYields.ReplaceAll("selected_events","yields");
  gSystem->mkdir(outputDirYields,kTRUE);
  TString fNameOutYields(outputDirYields+TString("/yields"));
//   fNameOutYields += TString::Format("%d_%d",stamp.GetDate(),stamp.GetTime());
  fNameOutYields += ".root";
  TFile fYields( fNameOutYields, "recreate" );
  YieldsSignal      .Write("YieldsSignal");
  YieldsSignalErr   .Write("YieldsSignalErr");
  BinLimitsForYields.Write("BinLimitsForYields");
  fYields.Close();

  // Save mass histograms into a separate file
  TString fNameOutHists(outputDirYields+"/massHist");
//   fNameOutHists += TString::Format("%d_%d",stamp.GetDate(),stamp.GetTime());
  fNameOutHists += ".root";
  TFile fMassHists(fNameOutHists,"recreate");
  for(UInt_t isam=0; isam<samplev.size(); isam++) {
    hMassBinsv[isam]->Write(snamev[isam]);
  }
  fMassHists.Close();

  
  gBenchmark->Show("plotSelectDY");      
}

