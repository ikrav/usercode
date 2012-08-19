#if !defined(__CINT__) || defined(__MAKECINT__)
#include "../Include/fitFunctions.hh"
#include "../Include/MyTools.hh"
#include <TEntryList.h>
#endif

/*
void measurePassAndFail(double &signal, double &signalErr, double &efficiency, double &efficiencyErr,TTree *passTree, TTree *failTree,TCanvas *passCanvas, TCanvas *failCanvas,const char* setBinsType){

  // Define data sets
  RooRealVar mass("mass","mass",60, 120);
  mass.setBins(120);
  RooRealVar pt ("et" ,"et" ,10.0, 1000);
  RooRealVar eta("eta","eta",-10, 10);
  RooDataSet *dataPass = new RooDataSet("dataPass","dataPass",passTree,RooArgSet(mass,pt,eta));
  RooDataSet *dataFail = new RooDataSet("dataFail","dataFail",RooArgList(mass,pt,eta),Import(*failTree));
  RooCategory probeType("probeType","probeType");
  probeType.defineType("pass");
  probeType.defineType("fail");
  RooAbsData *data;
  // If needed do binned fit
  bool unbinnedFit = false;
  if(unbinnedFit){
    data = new RooDataSet("data","data",mass,Index(probeType),
			  Import("pass",*dataPass), Import("fail",*dataFail));
  }else{
    RooDataHist *dataPassBinned = dataPass->binnedClone("dataPassBinned","dataPassBinned");
    RooDataHist *dataFailBinned = dataFail->binnedClone("dataFailBinned","dataFailBinned");
    data = new RooDataHist("data","data",mass,Index(probeType),
			   Import("pass",*dataPassBinned), Import("fail",*dataFailBinned));    
  }
  // Define the PDFs
  //
  // Common pieces
  //
  // True signal model
  RooRealVar zMass ("zMass" ,"zMass" ,91.188);
  RooRealVar zWidth("zWidth","zWidth",2.495);
  RooBreitWigner bwPdf("bwPdf","bwPdf",mass,zMass,zWidth);
  // Total signal events and efficiency
  RooRealVar nsignal("nsignal","nsignal",1000,0.0,1.0e7);
  RooRealVar eff    ("eff"    ,"eff"    ,0.7,0.0,1.0);

  //  PDF for the PASS sample
  // 
  // Background
  RooRealVar lambdaBgPass("lambdaBgPass","lambdaBgPass",-0.1, -0.5, 0.0);
  RooExponential bgPassPdf("bgPassPdf","bgPassPdf",mass,lambdaBgPass);
  // Signal
  //     - resolution function
  RooRealVar cbMeanPass("cbMeanPass","cbMeanPass"   ,0.0, -10.0,10.0);
  RooRealVar cbWidthPass("cbWidthPass","cbWidthPass",1.0,  0.1, 5.0);
  RooRealVar cbAlphaPass("cbAlphaPass","cbAlphaPass",5.0,  0.0,20.0);
  RooRealVar cbNPass("cbNPass","cbNPass"            ,1.0,  0.0,10.0);
  RooCBShape cbPassPdf("cbPassPdf","cbPassPdf",mass,cbMeanPass,cbWidthPass,cbAlphaPass,cbNPass);
  //     - realistic model
  mass.setBins(10000,setBinsType);
  RooFFTConvPdf signalPassPdf("signalPassPdf","signalPassPdf",mass, bwPdf, cbPassPdf);
  // Combine signal and background
  RooFormulaVar nsigPass("nsigPass","nsigPass","@0*@1",RooArgList(nsignal,eff));
  RooRealVar nbgPass ("nbgPass" ,"nbgPass" ,1,0.0,1.0e5);
  RooAddPdf passPdf("passPdf","passPdf",RooArgList(signalPassPdf,bgPassPdf), RooArgList(nsigPass,nbgPass));
  
  //
  //  PDF for the FAIL sample
  // 
  // Background
  RooRealVar lambdaBgFail("lambdaBgFail","lambdaBgFail",-0.1, -0.5, 0.0);
  RooExponential bgFailPdf("bgFailPdf","bgFailPdf",mass,lambdaBgFail);
  // Signal
  //     - resolution function
  RooRealVar cbMeanFail("cbMeanFail","cbMeanFail"   ,0.0, -10.0,10.0);
  RooRealVar cbWidthFail("cbWidthFail","cbWidthFail",1.0,   0.1, 5.0);
  RooRealVar cbAlphaFail("cbAlphaFail","cbAlphaFail",5.0,   0.0,20.0);
  RooRealVar cbNFail("cbNFail","cbNFail"            ,1.0,   0.0,10.0);
  RooCBShape cbFailPdf("cbFailPdf","cbFailPdf",mass,cbMeanFail,cbWidthFail,cbAlphaFail,cbNFail);
  //     - realistic model
  RooFFTConvPdf signalFailPdf("signalFailPdf","signalFailPdf",mass, bwPdf, cbFailPdf);
  // Combine signal and background
  RooFormulaVar nsigFail("nsigFail","nsigFail","@0*(1.0-@1)",RooArgList(nsignal,eff));
  RooRealVar nbgFail ("nbgFail" ,"nbgFail" ,1,0.0,1.0e5);
  RooAddPdf failPdf("failPdf","failPdf",RooArgList(signalFailPdf,bgFailPdf), RooArgList(nsigFail,nbgFail));
  
  // Combine pass and fail
  RooSimultaneous fullPdf("fullPdf","fullPdf",probeType);
  fullPdf.addPdf(passPdf,"pass");
  fullPdf.addPdf(failPdf,"fail");

  // Do the fit

  // Start with a reasonable point and do rough approximation first
  double total = dataPass->numEntries()+dataFail->numEntries();
  nsignal.setVal( 0.99*total);
  eff.setVal(0.90);
  nbgPass.setVal(0.01*total);
  nbgFail.setVal(0.01*total);
  cbAlphaPass.setVal(1.0);
  cbAlphaFail.setVal(0.5);
  cbNPass    .setVal(5.0);
  cbNFail    .setVal(5.0);
  cbAlphaPass.setConstant(kTRUE);
  cbAlphaFail.setConstant(kTRUE);
  cbNPass    .setConstant(kTRUE);
  cbNFail    .setConstant(kTRUE);
  RooFitResult *result = fullPdf.fitTo(*data,Extended(kTRUE),Save(),RooFit::NumCPU(2,true));

  // Release shape parameters and refine the fit
  cbAlphaPass.setConstant(kFALSE);
  cbAlphaFail.setConstant(kFALSE);
  cbNPass    .setConstant(kFALSE);
  cbNFail    .setConstant(kFALSE);
  result = fullPdf.fitTo(*data,Extended(kTRUE),Save(),RooFit::NumCPU(2,true));

  cout << "Fit status 1st iteration " << result->status() << endl;
//   if(!result->status()){
//     result = fullPdf.fitTo(*data,Extended(kTRUE),Save());
//     cout << "Fit status 2d iteration " << result->status() << endl;
//   }

  // Plot
  passCanvas->cd();
  passCanvas->SetWindowPosition(0,0);
  passCanvas->Draw();
  RooPlot *framePass = mass.frame();
  dataPass->plotOn(framePass);
  passPdf.plotOn(framePass);
  passPdf.plotOn(framePass,Components("bgPassPdf"),LineStyle(kDashed));
  framePass->Draw();

  failCanvas->cd();
  failCanvas->SetWindowPosition(0+ failCanvas->GetWindowWidth(),0);
  failCanvas->Draw();
  RooPlot *frameFail = mass.frame();
  dataFail->plotOn(frameFail);
  failPdf.plotOn(frameFail);
  failPdf.plotOn(frameFail,Components("bgFailPdf"),LineStyle(kDashed));
  frameFail->Draw();

  signal        = nsignal.getVal();
  signalErr     = nsignal.getError();
  efficiency    = eff.getVal();
  efficiencyErr = eff.getError();

  return;
}
*/

void measureEfficiency(TTree *passTree, TTree *failTree, 
		       int method, int etBinning, int etaBinning, 
		       TCanvas *canvas, ofstream &effOutput, ofstream &fitLog,
		       bool useTemplates, TFile *templatesFile, 
		       TFile *resultsRootFile, TFile *plotsRootFile,
		       int NsetBins, bool isRECO, const char* setBinsType, 
		       TString dirTag, const TString &picFileExtraTag, 
		       int puBin) {
  // puBin is important for the fit

  // For COUNTnCOUNT method we should write to root file results
  // from measureEfficiencyCountAndCount routine, otherwise
  // from measureEfficiencyWithFit routine.
  bool saveCountingToRootFile = true;
  if( method == DYTools::COUNTnFIT || method == DYTools::FITnFIT )
    saveCountingToRootFile = false;
  
  // Always report counting method results
  measureEfficiencyCountAndCount(passTree, failTree, etBinning, etaBinning, 
				 canvas, effOutput, 
				 saveCountingToRootFile, resultsRootFile, 
				 plotsRootFile);
  
  if( method == DYTools::COUNTnFIT || method == DYTools::FITnFIT )
    measureEfficiencyWithFit(passTree, failTree, 
			     method, etBinning, etaBinning, 
			     canvas, effOutput, fitLog,
			     useTemplates, templatesFile, 
			     resultsRootFile, plotsRootFile,
			     NsetBins, isRECO, setBinsType, dirTag, 
			     picFileExtraTag, puBin);
  

  return;
}

// -------------------------------------------------------------------------

void measureEfficiencyPU(TTree *passTreeFull, TTree *failTreeFull, 
			 int method, int etBinning, int etaBinning, 
			 TCanvas *canvas,ofstream &effOutput, ofstream &fitLog,
			 bool useTemplates, TFile *templatesFile, 
			 const TString &resultRootFileBase,
			 int NsetBins, bool isRECO, const char* setBinsType, 
			 TString dirTag, const TString &picFileExtraTag,
			 int puDependence
			 ) {

  if (!puDependence) {
    fitLog << "\nmeasureEfficiencyPU: no PU-dependence requested\n";
    TString resRootFName=resultRootFileBase + TString(".root");
    TFile *resultsRootFile=new TFile(resRootFName,"recreate");
    TString resPlotsFName=resultRootFileBase + TString("-plots.root");
    TFile *resultPlotsFile=new TFile(resPlotsFName,"recreate");
    measureEfficiency(passTreeFull,failTreeFull,method,etBinning,etaBinning,
		      canvas,effOutput,fitLog,useTemplates,
		      templatesFile,resultsRootFile,resultPlotsFile,
		      NsetBins,isRECO,setBinsType,
		      dirTag,picFileExtraTag);
  }
  else {
    // prepare pu-dependent trees
    vector<TTree*> passTreeV, failTreeV;
    passTreeV.reserve(DYTools::nPVBinCount);
    failTreeV.reserve(DYTools::nPVBinCount);
    Double_t mass,et,eta;
    UInt_t nPV;

    // 1. create trees
    for (int pass=0; pass<2; ++pass) {
      for (int pu_i=0; pu_i<DYTools::nPVBinCount; ++pu_i) {
	char buf[50];
	UInt_t pvMin=UInt_t(DYTools::nPVLimits[pu_i  ]+0.6);
	UInt_t pvMax=UInt_t(DYTools::nPVLimits[pu_i+1]-0.4);
	sprintf(buf,"%sTreePU_%u_%u",(pass)? "pass":"fail",pvMin,pvMax);
	//std::cout << "buf=" << buf << "\n";
	TTree *tree=new TTree(buf,buf);
	tree->SetDirectory(0);
	tree->Branch("mass",&mass,"mass/D");
	tree->Branch("et",&et,"et/D");
	tree->Branch("eta",&eta,"eta/D");
	if (pass) passTreeV.push_back(tree); else failTreeV.push_back(tree);
      }
    }

    // fill trees
    UInt_t excludedPass=0, excludedFail=0;
    for (int pass=0; pass<2; ++pass) {
      vector<TTree*> *treeV = (pass) ? &passTreeV : &failTreeV;
      TTree *treeFull=(pass) ? passTreeFull : failTreeFull;
      treeFull->SetBranchAddress("mass",&mass);
      treeFull->SetBranchAddress("et",&et);
      treeFull->SetBranchAddress("eta",&eta);
      treeFull->SetBranchAddress("nGoodPV",&nPV);
      for (UInt_t i=0; i<treeFull->GetEntries(); ++i) {
	treeFull->GetEntry(i);
	unsigned int pu_idx=(unsigned int)(DYTools::findPUBin(nPV));
	if (pu_idx<treeV->size()) {
	  (*treeV)[pu_idx]->Fill();
	}
	else {
	  if (pass) excludedPass++; else excludedFail++;
	}
      }
    }
    if (1) {
      effOutput << "\nEvent distribution by PU\n";
      effOutput << "  PU range    passCount    failCount\n";
      char buf[100];
      UInt_t nPass=excludedPass;
      UInt_t nFail=excludedFail;
      for (int pu_i=0; pu_i<DYTools::nPVBinCount; ++pu_i) {
	nPass+=passTreeV[pu_i]->GetEntries();
	nFail+=failTreeV[pu_i]->GetEntries();
	sprintf(buf,"  %4.1lf..%4.1lf    %8lld   %8lld", 
		DYTools::nPVLimits[pu_i],DYTools::nPVLimits[pu_i+1],
		passTreeV[pu_i]->GetEntries(), failTreeV[pu_i]->GetEntries());
	effOutput << buf << "\n";
      }
      sprintf(buf,"  excluded      %8u   %8u", excludedPass,excludedFail);
      effOutput << buf << "\n"
		<<"-------------------------------------\n";
      sprintf(buf,"  total         %8u   %8u", nPass,nFail);
      effOutput << buf << "\n";
    }
  
    // measure efficiency
    for (int pu_i=0; pu_i<DYTools::nPVBinCount; ++pu_i) {
      std::cout << " pu_i=" << pu_i << "\n";
      char buf[50];
      UInt_t pvMin=UInt_t(DYTools::nPVLimits[pu_i  ]+0.6);
      UInt_t pvMax=UInt_t(DYTools::nPVLimits[pu_i+1]-0.4);
      sprintf(buf,"_%u_%u.root",pvMin,pvMax);
      TString resRootFName=resultRootFileBase + TString(buf);
      TFile *resultsRootFile=new TFile(resRootFName,"recreate");
      sprintf(buf,"-plots_%u_%u.root",pvMin,pvMax);
      TString resPlotFName=resultRootFileBase + TString(buf);
      TFile *resultPlotFile=new TFile(resPlotFName,"recreate");
      fitLog << "\nmeasureEfficiencyPU: PU range " 
	     << pvMin << " - " << pvMax << "\n";
      effOutput << "\nmeasureEfficiencyPU: PU range " 
		<< pvMin << " - " << pvMax << "\n";
      std::cout << "call measure efficiency" << std::endl;
      measureEfficiency(passTreeV[pu_i],failTreeV[pu_i],
			method,etBinning,etaBinning,
			canvas,effOutput,fitLog,useTemplates,templatesFile,
			resultsRootFile,resultPlotFile,
			NsetBins,isRECO,setBinsType,
			dirTag,picFileExtraTag, pu_i+1);
      std::cout << "done measure efficiency" << std::endl;
    }

    ClearVec(passTreeV); ClearVec(failTreeV);
  }

  return;
}

void measureEfficiencyCountAndCount(TTree *passTree, TTree *failTree, 
			    int etBinning, int etaBinning, 
			    TCanvas *canvas, ofstream &effOutput,
			    bool saveResultsToRootFile, TFile *resultsRootFile,
			    TFile *resultPlotsFile){
  
  int nEt                = DYTools::getNEtBins(etBinning);
  const double *limitsEt = DYTools::getEtBinLimits(etBinning);
  
  int nEta                = DYTools::getNEtaBins(etaBinning);
  const double *limitsEta = DYTools::getEtaBinLimits(etaBinning);
  printf("eta bins %d\n", nEta);

  TMatrixD effArray2D(nEt, nEta);
  TMatrixD effArrayErrLow2D(nEt, nEta);
  TMatrixD effArrayErrHigh2D(nEt, nEta);
  TMatrixD effArray2DWeighted(nEt, nEta);
  TMatrixD effArrayErrLow2DWeighted(nEt, nEta);
  TMatrixD effArrayErrHigh2DWeighted(nEt, nEta);
 
  std::vector<std::string> lines;
  lines.reserve(50);

  effOutput << endl;
  effOutput << "Efficiency, counting method:\n";  
  effOutput << "     SC ET         SC eta           efficiency             pass         fail\n";

  lines.push_back("\nEfficiency, counting method (weighted):\n");
  lines.push_back("     SC ET         SC eta           efficiency             pass         fail\n");
  for(int j=0; j<nEta; j++){
    for(int i=0; i<nEt; i++){
      double effCount, effErrLowCount, effErrHighCount;
      TString etCut = TString::Format(" ( et >=%6.1f && et <%6.1f ) ",
				      limitsEt[i], limitsEt[i+1]);
      TString etaCut = 
	TString::Format(" ( abs(eta) >= %5.3f && abs(eta) < %5.3f ) ",
				       limitsEta[j], limitsEta[j+1]);
      TString cut = etCut + TString(" && ") + etaCut;
      //cout << cut << endl;
      double probesPass = passTree->GetEntries(cut);
      double probesFail = failTree->GetEntries(cut);

      canvas->cd(1 + 2*(i + j*nEt) + 0);
      passTree->Draw("weight >> hwPass",cut);
      canvas->cd(1 + 2*(i + j*nEt) + 1);
      failTree->Draw("weight >> hwFail",cut);
      TH1 *hwPass=(TH1*)gDirectory->Get("hwPass");
      TH1 *hwFail=(TH1*)gDirectory->Get("hwFail");
      double probesPassWeighted = hwPass->GetMean() * hwPass->GetEntries();
      double probesFailWeighted = hwFail->GetMean() * hwFail->GetEntries();
      double effCountWeighted, effErrLowCountWeighted, effErrHighCountWeighted;
      //std::cout << " probesPass=" << probesPass << ", probesFail=" << probesFail << "\n";
      //std::cout << " hwPass->GetEntries=" << hwPass->GetEntries() << ", hwFail->GetEntries=" << hwFail->GetEntries() << "\n";
      //std::cout << " probesPassWeighted=" << probesPassWeighted << ", probesFailWeighted=" << probesFailWeighted << "\n";

      DYTools::calcEfficiency( probesPass, probesPass+probesFail, 
			       DYTools::EFF_CLOPPER_PEARSON,
			       effCount, effErrLowCount, effErrHighCount);
      DYTools::calcEfficiency( probesPassWeighted, 
			       probesPassWeighted+probesFailWeighted,
			       DYTools::EFF_CLOPPER_PEARSON,
			       effCountWeighted, 
			       effErrLowCountWeighted,effErrHighCountWeighted);
      char strOut[200];
      sprintf(strOut, "hweightPass_Et_%1.0f-%1.0f__Eta_%5.3f-%5.3f",
	      limitsEt[i],limitsEt[i+1],
	      limitsEta[j],limitsEta[j+1]);
      hwPass->SetName(strOut);
      if (resultPlotsFile) { resultPlotsFile->cd(); hwPass->Write(); }
      sprintf(strOut, "hweightFail_Et_%1.0f-%1.0f__Eta_%5.3f-%5.3f",
	      limitsEt[i],limitsEt[i+1],
	      limitsEta[j],limitsEta[j+1]);
      hwFail->SetName(strOut);
      if (resultPlotsFile) { resultPlotsFile->cd(); hwFail->Write(); }

      sprintf(strOut, "   %3.0f - %3.0f   %5.3f - %5.3f   %5.1f +%5.1f -%5.1f    %10.0f  %10.0f\n",
	     limitsEt[i], limitsEt[i+1],
	     limitsEta[j], limitsEta[j+1],
	     effCount*100, effErrHighCount*100, effErrLowCount*100,
	     probesPass, probesFail);
      effOutput << strOut;

      sprintf(strOut, "   %3.0f - %3.0f   %5.3f - %5.3f   %5.1f +%5.1f -%5.1f    %10.0f  %10.0f\n",
	      limitsEt[i], limitsEt[i+1],
	      limitsEta[j], limitsEta[j+1],
	      effCountWeighted*100, 
	      effErrHighCountWeighted*100, effErrLowCountWeighted*100,
	      probesPassWeighted, probesFailWeighted);
      lines.push_back(strOut);

      canvas->cd(1 + 2*(i + j*nEt) + 0);
      passTree->Draw("mass",cut);
      canvas->cd(1 + 2*(i + j*nEt) + 1);
      failTree->Draw("mass",cut);
      canvas->Update();
      if (resultPlotsFile) {
	resultPlotsFile->cd();
	TString canvName=TString::Format("canvMassDistr_Et_%2.0f-%2.0f_abs_eta_%5.3f-%5.3f",
					 limitsEt[i],limitsEt[i+1],
					 limitsEta[j],limitsEta[j+1]);
	canvName.ReplaceAll(".","_");
	TCanvas ctemp(canvName,canvName,900,400);
	ctemp.Divide(2,1);
	ctemp.cd(1); passTree->Draw("mass",cut);
	ctemp.cd(2); failTree->Draw("mass",cut);
	ctemp.Write();
      }
      effArray2D(i,j) = effCount;
      effArrayErrLow2D(i,j) = effErrLowCount;
      effArrayErrHigh2D(i,j) = effErrHighCount;
      effArray2DWeighted(i,j) = effCountWeighted;
      effArrayErrLow2DWeighted(i,j) = effErrLowCountWeighted;
      effArrayErrHigh2DWeighted(i,j) = effErrHighCountWeighted;
    }
  }
  effOutput << endl;
  for (unsigned int i=0; i<lines.size(); ++i) {
    effOutput << lines[i];
  }
  effOutput << endl;

  if(saveResultsToRootFile){
    if(resultsRootFile && resultsRootFile->IsOpen()){
      resultsRootFile->cd();
      effArray2D.Write("effArray2D");
      effArrayErrLow2D.Write("effArrayErrLow2D");
      effArrayErrHigh2D.Write("effArrayErrHigh2D");    
      effArray2DWeighted.Write("effArray2DWeighted");
      effArrayErrLow2DWeighted.Write("effArrayErrLow2DWeighted");
      effArrayErrHigh2DWeighted.Write("effArrayErrHigh2DWeighted");    
      resultsRootFile->Close();
    }else assert(0);
  }

  if (resultPlotsFile && resultPlotsFile->IsOpen()) {
    resultPlotsFile->cd();
    canvas->Write();
  }

  return;
}

// --------------------------------------------------

void measureEfficiencyWithFit(TTree *passTree, TTree *failTree, 
		      int method, int etBinning, int etaBinning, 
		      TCanvas *canvas, ofstream &effOutput, ofstream &fitLog,
		      bool useTemplates, TFile *templatesFile, 
		      TFile *resultsRootFile, TFile *resultPlotsFile,
		      int NsetBins, bool isRECO, const char* setBinsType, 
		      TString dirTag, const TString &picFileExtraTag,
		      int puBin){
  
  int nEt                = DYTools::getNEtBins(etBinning);
  const double *limitsEt = DYTools::getEtBinLimits(etBinning);

  int nEta                = DYTools::getNEtaBins(etaBinning);
  const double *limitsEta = DYTools::getEtaBinLimits(etaBinning);
  printf("eta bins %d\n", nEta);
  
  TMatrixD effArray2D(nEt, nEta);
  TMatrixD effArrayErrLow2D(nEt, nEta);
  TMatrixD effArrayErrHigh2D(nEt, nEta);

  TString methodStr;
  switch(method) {
  case DYTools::COUNTnCOUNT: methodStr="Count+Count"; break;
  case DYTools::COUNTnFIT: methodStr="Count+Fit"; break;
  case DYTools::FITnFIT: methodStr="Fit+Fit"; break;
  default:
    std::cout << "unknown method of efficiency calculation\n";
    assert(0);
  }
 
  effOutput << endl;
  effOutput << "Efficiency, " << methodStr << " method:\n";  
  effOutput << "     SC ET         SC eta           efficiency             pass         fail\n";
  for(int j=0; j<nEta; j++){
    for(int i=0; i<nEt; i++){
  //for(int j=0; j<1; j++){
    //for(int i=0; i<1; i++){
       TString etCut = TString::Format(" ( et >=%6.1f && et <%6.1f ) ",
				      limitsEt[i], limitsEt[i+1]);
      TString etaCut = 
	TString::Format(" ( abs(eta) >= %5.3f && abs(eta) < %5.3f ) ",
				       limitsEta[j], limitsEta[j+1]);
      TString cut = etCut + TString(" && ") + etaCut;
      //cout << cut.Data() << endl;
      double probesPass = passTree->GetEntries(cut);
      double probesFail = failTree->GetEntries(cut);
      TPad *passPad = (TPad*)canvas->GetPad(1 + 2*(i + j*nEt) + 0);
      TPad *failPad = (TPad*)canvas->GetPad(1 + 2*(i + j*nEt) + 1);
      double efficiency, efficiencyErrHi, efficiencyErrLo;
      printf("\n ==\n");
      char strOut[200];
      sprintf(strOut,
	   " ==   Start fitting Et: %3.0f - %3.0f  and eta:  %5.3f - %5.3f \n",
	     limitsEt[i], limitsEt[i+1],
	     limitsEta[j], limitsEta[j+1]);
      printf("%s",strOut);
      printf(" ==\n\n");
      fitLog << endl << strOut << endl;
      
      if(!useTemplates){
	fitMass(passTree, failTree, cut, method, 
		efficiency,efficiencyErrHi, efficiencyErrLo, passPad, failPad, 
		resultPlotsFile,
		fitLog, NsetBins, isRECO, setBinsType, dirTag);
      }
      else{
	printf("\nMASS TEMPLATES ARE USED IN THE FIT\n\n");
	// In case templates are used, find the right templates
	TH1F *templatePass = 
	  getPassTemplate(i,j,etaBinning, templatesFile, puBin);
	TH1F *templateFail = 
	  getFailTemplate(i,j,etaBinning, templatesFile, puBin);

	if (0) {
	  TCanvas *cx=new TCanvas("cx","cx",600,600);
	  templatePass->SetLineColor(kGreen);
	  templatePass->DrawCopy("hist");
	  templateFail->SetLineColor(kRed);
	  templateFail->DrawCopy("hist same");
	  cx->Update();
	  double xx; std::cin>>xx;
	}
      
	fitMassWithTemplates(passTree, failTree, cut, method, 
			     efficiency, efficiencyErrHi, efficiencyErrLo,
			     passPad, failPad, resultPlotsFile,
			     fitLog, templatePass, templateFail, 
			     isRECO, setBinsType, dirTag, picFileExtraTag);
	//std::cout << "after fitMassWithTemplates: enter a value\n"; double x; std::cin >> x;
	char buf[50];
	sprintf(buf,"tmp_eta%d_et%d.png",j,i);
	canvas->SaveAs(buf);
 
      }
            

      sprintf(strOut, "   %3.0f - %3.0f   %5.3f - %5.3f   %5.1f +%5.1f -%5.1f    %10.0f  %10.0f\n",
	      limitsEt[i], limitsEt[i+1],
	      limitsEta[j], limitsEta[j+1],
	      efficiency*100, efficiencyErrHi*100, efficiencyErrLo*100,
	      probesPass, probesFail);
      effOutput << strOut;
      effArray2D(i,j) = efficiency;
      effArrayErrLow2D(i,j) = efficiencyErrLo;
      effArrayErrHigh2D(i,j) = efficiencyErrHi;
      
    }
  }

  effOutput << endl;

  if(resultsRootFile && resultsRootFile->IsOpen()){
    resultsRootFile->cd();
    effArray2D.Write("effArray2D");
    effArrayErrLow2D.Write("effArrayErrLow2D");
    effArrayErrHigh2D.Write("effArrayErrHigh2D");    
    resultsRootFile->Close();
  }else assert(0);

  if (resultPlotsFile && resultPlotsFile->IsOpen()) {
    resultPlotsFile->cd();
    canvas->Write();
  }

  //std::cout << "before leaving measureEfficiencyWithFit: enter a value\n"; double x; std::cin >> x;
   
  return;
}

// --------------------------------------------------

int getTemplateBin(int etBin, int etaBin, int etaBinning){

  int templateBin = -1;

  if( etBin != -1 && etaBin != -1)
    templateBin = etBin * DYTools::getNEtaBins(etaBinning) + etaBin;

  return templateBin;

}

// --------------------------------------------------

TString getTemplateName(int etBin, int etaBin, const char *pass_fail_str, 
			int puBin) {
  char buf[50];
  sprintf(buf,"hMassTemplate_Et%d_eta%d",etBin,etaBin);
  int len=0;
  if (puBin>0) {
    len=strlen(buf);
    int puMin=int(DYTools::nPVLimits[puBin-1]+0.6);
    int puMax=int(DYTools::nPVLimits[puBin  ]-0.4);
    sprintf(buf+len,"_%d_%d",puMin,puMax);
  }
  len=strlen(buf);
  sprintf(buf+len,"_%s",pass_fail_str);
  return TString(buf);
}

// --------------------------------------------------

TH1F * getPassTemplate(int etBin, int etaBin, int etaBinning, TFile *file, 
		       int puBin){

  TH1F *hist = 0;
  if(file == 0)
    return hist;

  int templateBin = getTemplateBin(etBin, etaBin, etaBinning);
  if( templateBin == -1 )
    return hist;

  TString name = getTemplateName(etBin,etaBin,"pass",puBin);

  hist = (TH1F*)file->Get(name);
  return hist;
}

// --------------------------------------------------

TH1F * getFailTemplate(int etBin, int etaBin, int etaBinning, TFile *file, 
		       int puBin){

  TH1F *hist = 0;
  if(file == 0)
    return hist;

  int templateBin = getTemplateBin(etBin, etaBin, etaBinning);
  if( templateBin == -1 )
    return hist;

  TString name = getTemplateName(etBin,etaBin,"fail",puBin);

  hist = (TH1F*)file->Get(name);
  return hist;
}
