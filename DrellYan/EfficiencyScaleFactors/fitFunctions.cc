#if !defined(__CINT__) || defined(__MAKECINT__)
#include "../Include/fitFunctions.hh"
#endif

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

void measureEfficiency(TTree *passTree, TTree *failTree, 
		       int method, int etBinning, int etaBinning, 
		       TCanvas *canvas, ofstream &effOutput, ofstream &fitLog,
		       bool useTemplates, TFile *templatesFile, TFile *resultsRootFile,
		       int NsetBins, bool isRECO, const char* setBinsType, TString dirTag){

  // For COUNTnCOUNT method we should write to root file results
  // from measureEfficiencyCountAndCount routine, otherwise
  // from measureEfficiencyWithFit routine.
  bool saveCountingToRootFile = true;
  if( method == COUNTnFIT || method == FITnFIT )
    saveCountingToRootFile = false;

  // Always report counting method results
  measureEfficiencyCountAndCount(passTree, failTree, etBinning, etaBinning, 
				 canvas, effOutput, saveCountingToRootFile, resultsRootFile);
  
  if( method == COUNTnFIT || method == FITnFIT )
    measureEfficiencyWithFit(passTree, failTree, 
			     method, etBinning, etaBinning, 
			     canvas, effOutput, fitLog,
			     useTemplates, templatesFile, resultsRootFile,
			     NsetBins, isRECO, setBinsType, dirTag);
  
  return;
}

void measureEfficiencyCountAndCount(TTree *passTree, TTree *failTree, 
				    int etBinning, int etaBinning, 
				    TCanvas *canvas, ofstream &effOutput,
				    bool saveResultsToRootFile, TFile *resultsRootFile){

  int nEt                = getNEtBins(etBinning);
  const double *limitsEt = getEtBinLimits(etBinning);

  int nEta                = getNEtaBins(etaBinning);
  const double *limitsEta = getEtaBinLimits(etaBinning);
  printf("eta bins %d\n", nEta);

  TMatrixD effArray2D(nEt, nEta);
  TMatrixD effArrayErrLow2D(nEt, nEta);
  TMatrixD effArrayErrHigh2D(nEt, nEta);
 
  effOutput << endl;
  effOutput << "Efficiency, counting method:\n";  
  effOutput << "     SC ET         SC eta           efficiency             pass         fail\n";
  for(int j=0; j<nEta; j++){
    for(int i=0; i<nEt; i++){
      double effCount, effErrLowCount, effErrHighCount;
      TString etCut = TString::Format(" ( et >=%6.1f && et <%6.1f ) ",
				      limitsEt[i], limitsEt[i+1]);
      TString etaCut = TString::Format(" ( abs(eta) >= %5.3f && abs(eta) < %5.3f ) ",
				       limitsEta[j], limitsEta[j+1]);
      TString cut = etCut + TString(" && ") + etaCut;
      //cout << cut.Data() << endl;
      double probesPass = passTree->GetEntries(cut);
      double probesFail = failTree->GetEntries(cut);
      DYTools::calcEfficiency( probesPass, probesPass+probesFail, DYTools::EFF_CLOPPER_PEARSON,
			       effCount, effErrLowCount, effErrHighCount);
      char strOut[200];
      sprintf(strOut, "   %3.0f - %3.0f   %5.3f - %5.3f   %5.1f +%5.1f -%5.1f    %10.0f  %10.0f\n",
	     limitsEt[i], limitsEt[i+1],
	     limitsEta[j], limitsEta[j+1],
	     effCount*100, effErrLowCount*100, effErrHighCount*100,
	     probesPass, probesFail);
      effOutput << strOut;
      canvas->cd(1 + 2*(i + j*nEt) + 0);
      passTree->Draw("mass",cut);
      canvas->cd(1 + 2*(i + j*nEt) + 1);
      failTree->Draw("mass",cut);
      canvas->Update();
      effArray2D(i,j) = effCount;
      effArrayErrLow2D(i,j) = effErrLowCount;
      effArrayErrHigh2D(i,j) = effErrHighCount;
    }
  }
  effOutput << endl;

  if(saveResultsToRootFile){
    if(resultsRootFile && resultsRootFile->IsOpen()){
      resultsRootFile->cd();
      effArray2D.Write("effArray2D");
      effArrayErrLow2D.Write("effArrayErrLow2D");
      effArrayErrHigh2D.Write("effArrayErrHigh2D");    
      resultsRootFile->Close();
    }else assert(0);
  }

  return;
}

void measureEfficiencyWithFit(TTree *passTree, TTree *failTree, 
			      int method, int etBinning, int etaBinning, 
			      TCanvas *canvas, ofstream &effOutput, ofstream &fitLog,
			      bool useTemplates, TFile *templatesFile, TFile *resultsRootFile, 
			      int NsetBins, bool isRECO, const char* setBinsType, TString dirTag){
  
  int nEt                = getNEtBins(etBinning);
  const double *limitsEt = getEtBinLimits(etBinning);

  int nEta                = getNEtaBins(etaBinning);
  const double *limitsEta = getEtaBinLimits(etaBinning);
  printf("eta bins %d\n", nEta);
  
  TMatrixD effArray2D(nEt, nEta);
  TMatrixD effArrayErrLow2D(nEt, nEta);
  TMatrixD effArrayErrHigh2D(nEt, nEta);
 
  effOutput << endl;
  effOutput << "Efficiency, Count+Fit method:\n";  
  effOutput << "     SC ET         SC eta           efficiency             pass         fail\n";
  for(int j=0; j<nEta; j++){
    for(int i=0; i<nEt; i++){
  //for(int j=0; j<1; j++){
    //for(int i=0; i<1; i++){
       TString etCut = TString::Format(" ( et >=%6.1f && et <%6.1f ) ",
				      limitsEt[i], limitsEt[i+1]);
      TString etaCut = TString::Format(" ( abs(eta) >= %5.3f && abs(eta) < %5.3f ) ",
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
      sprintf(strOut," ==   Start fitting Et: %3.0f - %3.0f  and eta:  %5.3f - %5.3f \n",
	     limitsEt[i], limitsEt[i+1],
	     limitsEta[j], limitsEta[j+1]);
      printf("%s",strOut);
      printf(" ==\n\n");
      fitLog << endl << strOut << endl;
      // In case templates are used, find the right templates
      TH1F *templatePass = getPassTemplate(i,j,etaBinning, templatesFile);
      TH1F *templateFail = getFailTemplate(i,j,etaBinning, templatesFile);
      
      
      if(!useTemplates){
	fitMass(passTree, failTree, cut, method, efficiency, efficiencyErrHi, efficiencyErrLo, passPad, failPad, fitLog, NsetBins, isRECO, setBinsType, dirTag);
      }else{
	printf("\nMASS TEMPLATES ARE USED IN THE FIT\n\n");
	fitMassWithTemplates(passTree, failTree, cut, method, 
			     efficiency, efficiencyErrHi, efficiencyErrLo,
			     passPad, failPad, fitLog, templatePass, templateFail, 
			     isRECO, setBinsType, dirTag);
      }
            

      sprintf(strOut, "   %3.0f - %3.0f   %5.3f - %5.3f   %5.1f +%5.1f -%5.1f        %10.0f  %10.0f\n",
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
  
  return;
}

int getTemplateBin(int etBin, int etaBin, int etaBinning){

  int templateBin = -1;

  if( etBin != -1 && etaBin != -1)
    templateBin = etBin * getNEtaBins(etaBinning) + etaBin;

  return templateBin;

}

TH1F * getPassTemplate(int etBin, int etaBin, int etaBinning, TFile *file){

  TH1F *hist = 0;
  if(file == 0)
    return hist;

  int templateBin = getTemplateBin(etBin, etaBin, etaBinning);
  if( templateBin == -1 )
    return hist;

  TString name = "hMassTemplate_Et";
  name += etBin;
  name += "_eta";
  name += etaBin;
  name += "_pass";

  hist = (TH1F*)file->Get(name);
  return hist;
}

TH1F * getFailTemplate(int etBin, int etaBin, int etaBinning, TFile *file){

  TH1F *hist = 0;
  if(file == 0)
    return hist;

  int templateBin = getTemplateBin(etBin, etaBin, etaBinning);
  if( templateBin == -1 )
    return hist;

  TString name = "hMassTemplate_Et";
  name += etBin;
  name += "_eta";
  name += etaBin;
  name += "_fail";

  hist = (TH1F*)file->Get(name);
  return hist;
}
