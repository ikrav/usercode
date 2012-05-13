#if !defined(__CINT__) || defined(__MAKECINT__)
#include "../Include/fitFunctionsCore.hh"
#include <RooEffProd.h>
#endif

int performWeightedFit=1;
int bkgPassContainsErrorFunction=0;
int bkgFailContainsErrorFunction=0;


void printCorrelations(ostream& os, RooFitResult *res)
{
  ios_base::fmtflags flags = os.flags();
  const RooArgList *parlist = res->correlation("eff");
  
  os << "  Correlation Matrix" << endl;
  os << " --------------------" << endl;
  for(Int_t i=0; i<parlist->getSize(); i++) {
    for(Int_t j=0; j<parlist->getSize(); j++) 
      os << "  " << setw(7) << setprecision(4) << fixed << res->correlationMatrix()(i,j);    
    os << endl;
  }
  os.flags(flags);
}

struct RealLimit 
{
  double av;
  double lo;
  double hi;
};

void fitMass(TTree *passTree, TTree *failTree, TString cut, int mode, double &efficiency, double &efficiencyErrHi, double &efficiencyErrLo, TPad *passPad, TPad *failPad, ofstream &fitLog, int NsetBins, bool isRECO, const char* setBinsType, TString dirTag){
  
  // meaningless check, saving from compiler complaints
  if (dirTag.Length() && 0) std::cout << "fitMass : dirTag=" << dirTag << "\n";

  RealLimit lims[19];

  for (int j=0; j<19; j++)
    {
      lims[j].av=0.0;
      lims[j].lo=0.0;
      lims[j].hi=0.0;
    }

  //  [0].name="mass";
  //  [1].name="et";
  //  [2].name="eta";
  //  [3].name="zMass";
  //  [4].name="zWidth";
  //  [5].name="nsignal";
  //  [6].name="eff";
  //  [7].name="lambdaBgPass";
  //  [8].name="cbMeanPass";
  //  [9].name="cbWidthPass";
  //  [10].name="cbAlphaPass";
  //  [11].name="cbNPass";
  //  [12].name="nbgPass";
  //  [13].name="lambdaBgFail";
  //  [14].name="cbMeanFail";
  //  [15].name="cbWidthFail";
  //  [16].name="cbAlphaFail";
  //  [17].name="cbNFail";
  //  [18].name="nbgFail";
  
  if (!isRECO)//for id and hlt
    {
      lims[0].lo=60;    lims[0].hi=120;
      lims[1].lo=10.0;  lims[1].hi=1000.0;
      lims[2].lo=-10;   lims[2].hi=10;

      lims[3].av=91.188;  
      lims[4].av=2.495;

      lims[5].av=1000;  lims[5].lo=0.0;   lims[5].hi=1.0e7;
      lims[6].av=0.7;   lims[6].lo=0.0;   lims[6].hi=1.0;
      lims[7].av=-0.1;  lims[7].lo=-0.5;  lims[7].hi=0.5;
      lims[8].av=0.0;   lims[8].lo=-10.0; lims[8].hi=10.0;
      lims[9].av=1.0;   lims[9].lo=0.1;   lims[9].hi=5.0;
      lims[10].av=5.0;  lims[10].lo=0.0;  lims[10].hi=20.0;
      lims[11].av=1.0;  lims[11].lo=0.0;  lims[11].hi=10.0;
      lims[12].av=1.0;  lims[12].lo=0.0;  lims[12].hi=1.0e5;
      lims[13].av=-0.1; lims[13].lo=-0.5; lims[13].hi=0.5;
      lims[14].av=0.0;  lims[14].lo=-10.0;lims[14].hi=10.0;
      lims[15].av=1.0;  lims[15].lo=0.1;  lims[15].hi=5.0;
      lims[16].av=5.0;  lims[16].lo=0.0;  lims[16].hi=20.0;
      lims[17].av=1.0;  lims[17].lo=0.0;  lims[17].hi=10.0;
      lims[18].av=1.0;  lims[18].lo=0.0;  lims[18].hi=1.0e5;
      
    }
  else if (isRECO)//for reco
    {
      lims[0].lo=60;    lims[0].hi=120;
      lims[1].lo=10.0;  lims[1].hi=1000.0;
      lims[2].lo=-10;   lims[2].hi=10;

      lims[3].av=91.188;  
      lims[4].av=2.495;

      lims[5].av=1000;  lims[5].lo=0.0;   lims[5].hi=1.0e7;
      lims[6].av=0.7;   lims[6].lo=0.0;   lims[6].hi=1.0;
      lims[7].av=-0.1;  lims[7].lo=-0.5;  lims[7].hi=0.5;
      lims[8].av=0.0;   lims[8].lo=-5.0;  lims[8].hi=5.0;
      lims[9].av=1.0;   lims[9].lo=0.1;   lims[9].hi=6.0;
      lims[10].av=5.0;  lims[10].lo=0.0;  lims[10].hi=20.0;
      lims[11].av=1.0;  lims[11].lo=0.0;  lims[11].hi=10.0;
      lims[12].av=1.0;  lims[12].lo=0.0;  lims[12].hi=1.0e5;
      lims[13].av=-0.1; lims[13].lo=-0.5; lims[13].hi=0.5;
      lims[14].av=0.0;  lims[14].lo=-5.0; lims[14].hi=5.0;
      lims[15].av=1.0;  lims[15].lo=0.1;  lims[15].hi=6.0;
      lims[16].av=5.0;  lims[16].lo=0.0;  lims[16].hi=20.0;
      lims[17].av=1.0;  lims[17].lo=0.0;  lims[17].hi=10.0;
      lims[18].av=1.0;  lims[18].lo=0.0;  lims[18].hi=1.0e5;      
    }

  
  // Define data sets
  RooRealVar mass("mass","mass",lims[0].lo, lims[0].hi);
  mass.setBins(NsetBins);
  RooRealVar et ("et" ,"et" ,lims[1].lo, lims[1].hi);
  RooRealVar eta("eta","eta",lims[2].lo, lims[2].hi);
  RooRealVar weight("weight","weight",0.,100.);
  RooFormulaVar rooCut("rooCut","rooCut",cut,RooArgSet(et,eta));
  RooArgSet dsetArgs(mass,et,eta);
  int weightedFit=0;
  if (performWeightedFit && (passTree->GetBranch("weight")!=NULL)) {
    weightedFit=1;
    dsetArgs.add(weight);
  }
  RooDataSet  *dataUnbinnedPass = 
    new RooDataSet("dataUnbinnedPass","dataUnbinnedPass",
		   passTree,dsetArgs, rooCut);
  RooDataSet *dataUnbinnedPassUnweighted=NULL;
  RooDataSet  *dataUnbinnedFail = 
    new RooDataSet("dataUnbinnedFail","dataUnbinnedFail",
						 failTree,dsetArgs, rooCut);
  RooDataSet *dataUnbinnedFailUnweighted=NULL;
  std::string dline(70,'-'); dline+='\n';

  if (weightedFit) {
    std::cout << "\n\n\tweighted Fit\n\n";
    dataUnbinnedPassUnweighted=dataUnbinnedPass;
    RooDataSet *dSet=dataUnbinnedPassUnweighted;
    dataUnbinnedPass= 
      new RooDataSet("dataUnbinnedPassWeighted","dataUnbinnedPassWeighted",
		     dSet, *dSet->get(), 0, "weight");
    dataUnbinnedFailUnweighted=dataUnbinnedFail;
    dSet=dataUnbinnedFailUnweighted;
    dataUnbinnedFail=
      new RooDataSet("dataUnbinnedFailWeighted","dataUnbinnedFailWeighted",
		     dSet, *dSet->get(), 0, "weight");
    std::cout << dline; dataUnbinnedPass->Print("V"); std::cout << dline;
    std::cout << dline; dataUnbinnedFail->Print("V"); std::cout << dline;
    std::cout << std::endl;    
  }
  else std::cout << "\n\n\tunweighted Fit\n\n";

  RooDataHist *dataBinnedPass   = dataUnbinnedPass->binnedClone("dataBinnedPass","dataBinnedPass");
  RooDataHist *dataBinnedFail   = dataUnbinnedFail->binnedClone("dataBinnedFail","dataBinnedFail");
  RooCategory probeType("probeType","probeType");
  probeType.defineType("pass");
  probeType.defineType("fail");
  RooAbsData *data;

  // If needed do binned fit
  bool unbinnedFit = true;
  if(unbinnedFit){
    RooArgSet combiDSetArgs(mass);
    if (weightedFit) dsetArgs.add(weight);
    data = new RooDataSet("data","data",combiDSetArgs,Index(probeType),
			  Import("pass",*dataUnbinnedPass), 
			  Import("fail",*dataUnbinnedFail));
    std::cout << dline; data->Print("V"); std::cout << dline;
    cout << endl << "Setting up UNBINNED fit" << endl << endl;
  }else{
    data = new RooDataHist("data","data",mass,Index(probeType),
			   Import("pass",*dataBinnedPass), 
			   Import("fail",*dataBinnedFail));    
    cout << endl << "Setting up BINNED fit" << endl << endl;
  }
  // Define the PDFs
  //
  // Common pieces
  //
  // True signal model
  RooRealVar zMass ("zMass" ,"zMass" ,lims[3].av);
  RooRealVar zWidth("zWidth","zWidth",lims[4].av);
  RooBreitWigner bwPdf("bwPdf","bwPdf",mass,zMass,zWidth);
  // Total signal events and efficiency
  RooRealVar nsignal("nsignal","nsignal",lims[5].av,lims[5].lo,lims[5].hi);
  RooRealVar eff("eff","eff",lims[6].av,lims[6].lo,lims[6].hi);

  //  PDF for the PASS sample
  // 
  // Background
  RooRealVar lambdaBgPass("lambdaBgPass","lambdaBgPass",lims[7].av, lims[7].lo, lims[7].hi);
  RooExponential bgPassPdf("bgPassPdf","bgPassPdf",mass,lambdaBgPass);
  RooRealVar bgPassErfHalfPos("bgPassErfHalfPos","bgPassErfHalfPos",0, lims[0].hi);
  RooFormulaVar bgPassErf("bgPassErf","0.5*(TMath::Erf((@0-@1)/0.5)+1.)",RooArgList(mass,bgPassErfHalfPos));
  RooEffProd bgPassTotPdf("bgPassTotPdf","pass background model: exp * erf",bgPassPdf,bgPassErf);

  // Signal
  //     - resolution function
  RooRealVar cbMeanPass("cbMeanPass","cbMeanPass"   ,lims[8].av, lims[8].lo, lims[8].hi);
  RooRealVar cbWidthPass("cbWidthPass","cbWidthPass",lims[9].av, lims[9].lo, lims[9].hi);
  RooRealVar cbAlphaPass("cbAlphaPass","cbAlphaPass",lims[10].av, lims[10].lo, lims[10].hi);
  RooRealVar cbNPass("cbNPass","cbNPass"            ,lims[11].av, lims[11].lo, lims[11].hi);
  RooCBShape cbPassPdf("cbPassPdf","cbPassPdf",mass,cbMeanPass,cbWidthPass,cbAlphaPass,cbNPass);
  //     - realistic model
  mass.setBins(10000,setBinsType);
  RooFFTConvPdf signalPassPdf("signalPassPdf","signalPassPdf",mass, bwPdf, cbPassPdf);
  // Combine signal and background
  RooFormulaVar nsigPass("nsigPass","nsigPass","@0*@1",RooArgList(nsignal,eff));
  RooRealVar nbgPass ("nbgPass","nbgPass",lims[12].av,lims[12].lo,lims[12].hi);
  RooAbsPdf *passPdf=NULL;
  if( mode == COUNTnFIT ){
    RooGenericPdf * simpleSignal = new RooGenericPdf("simpleSignal","simpleSignal",
						     "1.0",RooArgList());
    RooExtendPdf * simpleSignalExtended = new RooExtendPdf("passPdf", "passPdf",
							   *simpleSignal, nsigPass);
    passPdf = simpleSignalExtended;
  }else if( mode == FITnFIT ){
    if (bkgPassContainsErrorFunction) {
      cout << "pass background model: exponential times error function\n";
      fitLog << "pass background model: exponential times error function\n";
      passPdf = new RooAddPdf("passPdf","passPdf",RooArgList(signalPassPdf,bgPassTotPdf),RooArgList(nsigPass,nbgPass));
    }
    else {
      passPdf = new RooAddPdf("passPdf","passPdf",RooArgList(signalPassPdf,bgPassPdf), RooArgList(nsigPass,nbgPass));
    }
  }else{
    printf("ERROR: inappropriate mode requested\n");
    return;
  }
  //
  //  PDF for the FAIL sample
  // 
  // Background
  RooRealVar lambdaBgFail("lambdaBgFail","lambdaBgFail",lims[13].av, lims[13].lo, lims[13].hi);
  RooExponential bgFailPdf("bgFailPdf","bgFailPdf",mass,lambdaBgFail);
  RooRealVar bgFailErfHalfPos("bgFailErfHalfPos","bgFailErfHalfPos",0, lims[0].hi);
  RooFormulaVar bgFailErf("bgFailErf","0.5*(TMath::Erf((@0-@1)/0.5)+1.)",RooArgList(mass,bgFailErfHalfPos));
  RooEffProd bgFailTotPdf("bgFailTotPdf","fail background model: exp * erf",bgFailPdf,bgFailErf);

  // Signal
  //     - resolution function
  RooRealVar cbMeanFail("cbMeanFail","cbMeanFail",   lims[14].av,lims[14].lo,lims[14].hi);
  RooRealVar cbWidthFail("cbWidthFail","cbWidthFail",lims[15].av, lims[15].lo,lims[15].hi);
  RooRealVar cbAlphaFail("cbAlphaFail","cbAlphaFail",lims[16].av,lims[16].lo,lims[16].hi);
  RooRealVar cbNFail("cbNFail","cbNFail"            ,lims[17].av,    lims[17].lo,    lims[17].hi);
  RooCBShape cbFailPdf("cbFailPdf","cbFailPdf",mass,cbMeanFail,cbWidthFail,cbAlphaFail,cbNFail);
  //     - realistic model
  RooFFTConvPdf signalFailPdf("signalFailPdf","signalFailPdf",mass, bwPdf, cbFailPdf);
  // Combine signal and background
  RooFormulaVar nsigFail("nsigFail","nsigFail","@0*(1.0-@1)",RooArgList(nsignal,eff));
  RooRealVar nbgFail ("nbgFail" ,"nbgFail" ,lims[18].av,lims[18].lo,lims[18].hi);
  RooAddPdf *failPdf=NULL;
  if (bkgFailContainsErrorFunction) {
    cout << "fail background model: exponential times error function\n";
    fitLog << "fail background model: exponential times error function\n";
    failPdf=new RooAddPdf("failPdf","failPdf",RooArgList(signalFailPdf,bgFailTotPdf), RooArgList(nsigFail,nbgFail));
  }
  else {
    failPdf=new RooAddPdf("failPdf","failPdf",RooArgList(signalFailPdf,bgFailPdf), RooArgList(nsigFail,nbgFail));
  }
  
  // Combine pass and fail
  RooSimultaneous fullPdf("fullPdf","fullPdf",probeType);
  fullPdf.addPdf(*passPdf,"pass");
  fullPdf.addPdf(*failPdf,"fail");

  
  // Do the fit
  // Start with a reasonable point and do rough approximation first
  double total = dataUnbinnedPass->numEntries() + dataUnbinnedFail->numEntries();
  nsignal.setVal( 0.99*total);
  eff.setVal(0.90);
  if( mode == FITnFIT ){
    nbgPass.setVal(0.01*total);
    cbAlphaPass.setVal(1.0);
    cbNPass    .setVal(5.0);
    cbAlphaPass.setConstant(kTRUE);
    cbNPass    .setConstant(kTRUE);
  }
  nbgFail.setVal(0.01*total);
  cbAlphaFail.setVal(0.5);
  cbNFail    .setVal(5.0);
  cbAlphaFail.setConstant(kTRUE);
  cbNFail    .setConstant(kTRUE);
  RooFitResult *result = fullPdf.fitTo(*data,
				       Extended(kTRUE),
				       Save(),
				       NumCPU(2,true));

  // Release shape parameters and refine the fit
  if( mode == FITnFIT ){
    cbAlphaPass.setConstant(kFALSE);
    cbNPass    .setConstant(kFALSE);
  }
  cbAlphaFail.setConstant(kFALSE);
  cbNFail    .setConstant(kFALSE);
  result = fullPdf.fitTo(*data,
			 Extended(kTRUE),
			 Minos(RooArgSet(eff)),
			 Save(),
			 NumCPU(2,true));
  // If minos fails, refit without minos
  if((fabs(eff.getErrorLo())<5e-5) || (eff.getErrorHi()<5e-5)){
    cout << "MINOS FAILS" << endl;
    result = fullPdf.fitTo(*data,
			   Extended(kTRUE),
			   Save(),
			   NumCPU(2,true));
  }

  efficiency       = eff.getVal();
  efficiencyErrHi  = eff.getErrorHi();
  efficiencyErrLo  = fabs(eff.getErrorLo());

  // Draw fit results
  passPad->cd();
  passPad->Clear();
  RooPlot *framePass = mass.frame();
  dataUnbinnedPass->plotOn(framePass);
  if(mode == FITnFIT){
    passPdf->plotOn(framePass);
    passPdf->plotOn(framePass,Components("bgPassPdf"),LineStyle(kDashed));
  }
  framePass->Draw();
  passPad->Update();


  failPad->cd();
  failPad->Clear();
  RooPlot *frameFail = mass.frame();
  dataUnbinnedFail->plotOn(frameFail);
  failPdf->plotOn(frameFail);
  failPdf->plotOn(frameFail,Components("bgFailPdf"),LineStyle(kDashed));
  frameFail->Draw();
  failPad->Update();


  // Print fit outcome into fit log
  result->printStream(fitLog,RooPrintable::kValue,RooPrintable::kVerbose);
  fitLog << endl;
  printCorrelations(fitLog, result);

  //std::cout << "enter a value\n"; double x; std::cin >> x;
  
  return;
}


void fitMassWithTemplates(TTree *passTree, TTree *failTree, TString cut, int mode, double &efficiency, double &efficiencyErrHi, double &efficiencyErrLo, TPad *passPad, TPad *failPad, ofstream &fitLog, TH1F *templatePass, TH1F *templateFail, bool isRECO, const char* setBinsType, TString dirTag, const TString &picFileExtraTag){

    
  RealLimit lims[12];
  
  for (int j=0; j<12; j++)
    {
      lims[j].av=0.0;
      lims[j].lo=0.0;
      lims[j].hi=0.0;
    }

  //  [0].name="mass";
  //  [1].name="et";
  //  [2].name="eta";
  //  [3].name="nsignal";
  //  [4].name="eff";
  //  [5].name="lambdaBgPass";
  //  [6].name="resMeanPass";
  //  [7].name="resSigma";
  //  [8].name="nbgPass";
  //  [9].name="lambdaBgFail";
  //  [10].name="resMeanFail";
  //  [11].name="nbgFail";
  
  if (!isRECO)
    {
      lims[0].lo=60;    lims[0].hi=120;
      lims[1].lo=10.0;  lims[1].hi=1000.0;
      lims[2].lo=-10;   lims[2].hi=10;

      lims[3].av=1000;  lims[3].lo=0.0;   lims[3].hi=1.0e7;
      lims[4].av=0.7;   lims[4].lo=0.0;   lims[4].hi=1.0;
      lims[5].av=-0.1;  lims[5].lo=-0.5;  lims[5].hi=0.5;
      lims[6].av=0.0;   lims[6].lo=-5.0;  lims[6].hi=5.0;
      lims[7].av=1.0;   lims[7].lo=0.1;   lims[7].hi=6.0;
      lims[8].av=1.0;   lims[8].lo=0.0;   lims[8].hi=1.0e5;
      lims[9].av=-0.1;  lims[9].lo=-0.5;  lims[9].hi=0.5;
      lims[10].av=0.0;  lims[10].lo=-5.0; lims[10].hi=5.0;
      lims[11].av=1.0;  lims[11].lo=0.0;  lims[11].hi=1.0e5;
      
    }
  else if (isRECO)
    {
      lims[0].lo=60;    lims[0].hi=120;
      lims[1].lo=10.0;  lims[1].hi=1000.0;
      lims[2].lo=-10;   lims[2].hi=10;

      lims[3].av=1000;  lims[3].lo=0.0;   lims[3].hi=1.0e7;
      lims[4].av=0.7;   lims[4].lo=0.0;   lims[4].hi=1.0;
      lims[5].av=-0.1;  lims[5].lo=-0.5;  lims[5].hi=0.5;
      lims[6].av=0.0;   lims[6].lo=-5.0;  lims[6].hi=5.0;
      lims[7].av=1.0;   lims[7].lo=0.1;   lims[7].hi=6.0;
      lims[8].av=1.0;   lims[8].lo=0.0;   lims[8].hi=1.0e5;
      lims[9].av=-0.1;  lims[9].lo=-0.5;  lims[9].hi=0.5;
      lims[10].av=0.0;  lims[10].lo=-3.0; lims[10].hi=3.0;
      lims[11].av=1.0;  lims[11].lo=0.0;  lims[11].hi=1.0e5;
      
    }
  
  /*  
    std::cout<<"beg-llllllllllllllllllllllimssssssssssssssssssssssss"<<std::endl;  
    
    for (int j=0; j<12; j++)
    {
    std::cout<<lims[j].name<<": "<<lims[j].av<<", "<<lims[j].lo<<", "<<lims[j].hi<<std::endl;
    }
    
    std::cout<<"end-llllllllllllllllllllllimssssssssssssssssssssssss"<<std::endl;  
  */
  
  
  // Define data sets
  RooRealVar mass("mass","mass",lims[0].lo, lims[0].hi);
  mass.setBins(30);
  RooRealVar et ("et" ,"et" ,lims[1].lo, lims[1].hi);
  RooRealVar eta("eta","eta",lims[2].lo, lims[2].hi);
  RooRealVar weight("weight","weight",0.,100.);
  RooFormulaVar rooCut("rooCut","rooCut",cut,RooArgSet(et,eta));
  RooArgSet dsetArgs(mass,et,eta);
  int weightedFit=0;
  if (performWeightedFit && (passTree->GetBranch("weight")!=NULL)) {
    weightedFit=1;
    dsetArgs.add(weight);
  }
  RooDataSet  *dataUnbinnedPass = 
    new RooDataSet("dataUnbinnedPass","dataUnbinnedPass",
		   passTree,dsetArgs, rooCut);
  RooDataSet *dataUnbinnedPassUnweighted=NULL;
  RooDataSet  *dataUnbinnedFail = 
    new RooDataSet("dataUnbinnedFail","dataUnbinnedFail",
		   failTree,dsetArgs, rooCut);
  RooDataSet *dataUnbinnedFailUnweighted=NULL;
  std::string dline(70,'-'); dline+='\n';

  if (weightedFit) {
    std::cout << "\n\n\tweighted Fit\n\n";
    dataUnbinnedPassUnweighted=dataUnbinnedPass;
    RooDataSet *dSet=dataUnbinnedPassUnweighted;
    dataUnbinnedPass= 
      new RooDataSet("dataUnbinnedPassWeighted","dataUnbinnedPassWeighted",
		     dSet, *dSet->get(), 0, "weight");
    dataUnbinnedFailUnweighted=dataUnbinnedFail;
    dSet=dataUnbinnedFailUnweighted;
    dataUnbinnedFail=
      new RooDataSet("dataUnbinnedFailWeighted","dataUnbinnedFailWeighted",
		     dSet, *dSet->get(), 0, "weight");
    std::cout << dline; dataUnbinnedPass->Print("V"); std::cout << dline;
    std::cout << dline; dataUnbinnedFail->Print("V"); std::cout << dline;
    std::cout << std::endl;
  }
  else std::cout << "\n\n\tunweighted Fit (" << dataUnbinnedPass->numEntries() << "p," << dataUnbinnedFail->numEntries() << "f)\n\n";

  RooDataHist *dataBinnedPass   = 
    dataUnbinnedPass->binnedClone("dataBinnedPass","dataBinnedPass");
  RooDataHist *dataBinnedFail   = 
    dataUnbinnedFail->binnedClone("dataBinnedFail","dataBinnedFail");
  RooCategory probeType("probeType","probeType");
  probeType.defineType("pass");
  probeType.defineType("fail");
  RooAbsData *data;

  // If needed do binned fit
  bool unbinnedFit = true;
  if(unbinnedFit){
    RooArgSet combiDSetArgs(mass);
    if (weightedFit) dsetArgs.add(weight);
    data = new RooDataSet("data","data",combiDSetArgs,Index(probeType),
			  Import("pass",*dataUnbinnedPass), 
			  Import("fail",*dataUnbinnedFail));
    std::cout << dline; data->Print("V"); std::cout << dline;
    cout << endl << "Setting up UNBINNED fit" << endl << endl;
  }else{
    data = new RooDataHist("data","data",mass,Index(probeType),
			   Import("pass",*dataBinnedPass), 
			   Import("fail",*dataBinnedFail));    
    cout << endl << "Setting up BINNED fit" << endl << endl;
  }
  // Define the PDFs
  //
  // Common pieces
  //
  // Total signal events and efficiency
  RooRealVar nsignal("nsignal","nsignal",lims[3].av,lims[3].lo,lims[3].hi);
  RooRealVar eff    ("eff"    ,"eff"    ,lims[4].av,lims[4].lo,lims[4].hi);

  //  PDF for the PASS sample
  // 
  // Background
  RooRealVar lambdaBgPass("lambdaBgPass","lambdaBgPass",lims[5].av, lims[5].lo, lims[5].hi);
  RooExponential bgPassPdf("bgPassPdf","bgPassPdf",mass,lambdaBgPass);
  RooRealVar bgPassErfHalfPos("bgPassErfHalfPos","bgPassErfHalfPos",0, lims[0].hi);
  RooFormulaVar bgPassErf("bgPassErf","0.5*(TMath::Erf((@0-@1)/0.5)+1.)",RooArgList(mass,bgPassErfHalfPos));
  RooEffProd bgPassTotPdf("bgPassTotPdf","pass background model: exp * erf",bgPassPdf,bgPassErf);
  // Signal
  //     - resolution function
  RooRealVar resMeanPass("resMeanPass","cbMeanPass"   ,lims[6].av, lims[6].lo, lims[6].hi);
  RooRealVar resSigma   ("resSigma"   ,"resSigma  "   ,lims[7].av, lims[7].lo, lims[7].hi);
  RooGaussian resPassPdf("resPassPdf","resPassPdf", mass, resMeanPass, resSigma);
  //      - mc template
  RooDataHist rooTemplatePass("rooTemplatePass","rooTemplatePass",RooArgList(mass),templatePass);
  RooHistPdf templatePassPdf("templatePassPdf","templatePassPdf",RooArgSet(mass),rooTemplatePass);
  //     - realistic model
  mass.setBins(10000,setBinsType);
//   RooFFTConvPdf signalPassPdf("signalPassPdf","signalPassPdf",mass, templatePassPdf, cbPassPdf);
  RooFFTConvPdf signalPassPdf("signalPassPdf","signalPassPdf",mass, templatePassPdf, resPassPdf);
  // Combine signal and background
  RooFormulaVar nsigPass("nsigPass","nsigPass","@0*@1",RooArgList(nsignal,eff));
  RooRealVar nbgPass ("nbgPass","nbgPass",lims[8].av, lims[8].lo, lims[8].hi);
  RooAbsPdf *passPdf;
  if( mode == COUNTnFIT ){
    RooGenericPdf * simpleSignal = new RooGenericPdf("simpleSignal","simpleSignal",
						     "1.0",RooArgList());
    RooExtendPdf * simpleSignalExtended = new RooExtendPdf("passPdf", "passPdf",
							   *simpleSignal, nsigPass);
    passPdf = simpleSignalExtended;
  }else if( mode == FITnFIT ){
    if (bkgPassContainsErrorFunction) {
      cout << "pass background model: exponential times error function\n";
      fitLog << "pass background model: exponential times error function\n";
      passPdf = new RooAddPdf("passPdf","passPdf",RooArgList(signalPassPdf,bgPassTotPdf),RooArgList(nsigPass,nbgPass));
    }
    else {
      passPdf = new RooAddPdf("passPdf","passPdf",RooArgList(signalPassPdf,bgPassPdf), RooArgList(nsigPass,nbgPass));
    }
  }else{
    printf("ERROR: inappropriate mode requested\n");
    return;
  }
  //
  //  PDF for the FAIL sample
  // 
  // Background
  RooRealVar lambdaBgFail("lambdaBgFail","lambdaBgFail",lims[9].av, lims[9].lo, lims[9].hi);
  RooExponential bgFailPdf("bgFailPdf","bgFailPdf",mass,lambdaBgFail);
  RooRealVar bgFailErfHalfPos("bgFailErfHalfPos","bgFailErfHalfPos",0, lims[0].hi);
  RooFormulaVar bgFailErf("bgFailErf","0.5*(TMath::Erf((@0-@1)/0.5)+1.)",RooArgList(mass,bgFailErfHalfPos));
  RooEffProd bgFailTotPdf("bgFailTotPdf","fail background model: exp * erf",bgFailPdf,bgFailErf);
  // Signal
  //     - resolution function
  // The limits for the "fail" come from stating at the fit results without
  // splitting into Et bins. In some cases, the peak is not there at all. (for reco)
  RooRealVar resMeanFail("resMeanFail","cbMeanFail"   ,lims[10].av, lims[10].lo, lims[10].hi);
  RooGaussian resFailPdf("resFailPdf","resFailPdf", mass, resMeanFail, resSigma);
  //      - mc template
  RooDataHist rooTemplateFail("rooTemplateFail","rooTemplateFail",RooArgList(mass),templateFail);
  RooHistPdf templateFailPdf("templateFailPdf","templateFailPdf",RooArgSet(mass),rooTemplateFail);
  //     - realistic model
  //   RooFFTConvPdf signalFailPdf("signalFailPdf","signalFailPdf",mass, templateFailPdf, cbFailPdf);
  RooFFTConvPdf signalFailPdf("signalFailPdf","signalFailPdf",mass, templateFailPdf, resFailPdf);
  // Combine signal and background
  RooFormulaVar nsigFail("nsigFail","nsigFail","@0*(1.0-@1)",RooArgList(nsignal,eff));
  RooRealVar nbgFail ("nbgFail" ,"nbgFail" ,lims[11].av,lims[11].lo,lims[11].hi);
  RooAddPdf *failPdf=NULL;
  if (bkgFailContainsErrorFunction) {
    cout << "fail background model: exponential times error function\n";
    fitLog << "fail background model: exponential times error function\n";
    failPdf=new RooAddPdf("failPdf","failPdf",RooArgList(signalFailPdf,bgFailTotPdf), RooArgList(nsigFail,nbgFail));
  }
  else {
    failPdf=new RooAddPdf("failPdf","failPdf",RooArgList(signalFailPdf,bgFailPdf), RooArgList(nsigFail,nbgFail));
  }
  
  // Combine pass and fail
  RooSimultaneous fullPdf("fullPdf","fullPdf",probeType);
  fullPdf.addPdf(*passPdf,"pass");
  fullPdf.addPdf(*failPdf,"fail");

  
  // Do the fit
  double total = dataUnbinnedPass->numEntries() + dataUnbinnedFail->numEntries();
  nsignal.setVal( 0.99*total);
  eff.setVal(0.90);
  if( mode == FITnFIT ){
    nbgPass.setVal(0.01*total);
  }
  nbgFail.setVal(0.01*total);

  RooFitResult *result;  
  std::cout << "fitting \n\n";
  result = fullPdf.fitTo(*data, Extended(kTRUE), Minos(RooArgSet(eff)), Save(), NumCPU(2,true));
  
  
  
  // If minos fails, refit without minos
  if((fabs(eff.getErrorLo())<5e-5) || (eff.getErrorHi()<5e-5)){
    cout << "\n\n\tMINOS FAILS\n\n" << endl;
    result = fullPdf.fitTo(*data, Extended(kTRUE), Save(), NumCPU(2,true));
  }
  
  efficiency     = eff.getVal();
  efficiencyErrHi  = eff.getErrorHi();
  efficiencyErrLo  = fabs(eff.getErrorLo());

  // Draw fit results
  passPad->cd();
  passPad->Clear();
  RooPlot *framePass = mass.frame();
  dataUnbinnedPass->plotOn(framePass);
  if(mode == FITnFIT){
    passPdf->plotOn(framePass);
    passPdf->plotOn(framePass,Components("bgPassPdf"),LineStyle(kDashed));
  }
  framePass->Draw();
  passPad->Update();

  TString cutF=cut.ReplaceAll(" ","_");
  cutF=cutF.ReplaceAll(".","_");
  cutF=cutF.ReplaceAll("&&","_");
  cutF=cutF.ReplaceAll(">","_");
  cutF=cutF.ReplaceAll("<","_");
  cutF=cutF.ReplaceAll("=","_");
  cutF=cutF.ReplaceAll("(","_");
  cutF=cutF.ReplaceAll(")","_");
  cutF=cutF.ReplaceAll("____","_");
  cutF=cutF.ReplaceAll("___","_");
  cutF=cutF.ReplaceAll("__","_");

  TString pngFileBase=TString("../root_files/tag_and_probe/") + dirTag + TString("/fit-") + picFileExtraTag;
  if (isRECO) pngFileBase+="-reco-"; else pngFileBase+="-id-";
  pngFileBase += cutF;
  
  TCanvas cPass("cPass","cPass");
  framePass->Draw();
  cPass.Update();

  /*
  if (isRECO) cPass.Print(("../root_files/tag_and_probe/"+(std::string)dirTag+"/fit-reco-"+(std::string)cutF+"-pass.png").c_str());
  else cPass.Print(("../root_files/tag_and_probe/"+(std::string)dirTag+"/fit-id-"+(std::string)cutF+"-pass.png").c_str());
  */
  cPass.Print((pngFileBase + TString("-pass.png")).Data());

  failPad->cd();
  failPad->Clear();
  RooPlot *frameFail = mass.frame();
  dataUnbinnedFail->plotOn(frameFail);
  failPdf->plotOn(frameFail);
  failPdf->plotOn(frameFail,Components("bgFailPdf"),LineStyle(kDashed));
  frameFail->Draw();
  failPad->Update();

  TCanvas cFail("cFail","cFail");
  frameFail->Draw();
  cFail.Update();
  /*
  if (isRECO) cFail.Print(("../root_files/tag_and_probe/"+(std::string)dirTag+"/fit-reco-"+(std::string)cutF+"-fail.png").c_str());
  else cFail.Print(("../root_files/tag_and_probe/"+(std::string)dirTag+"/fit-id-"+(std::string)cutF+"-fail.png").c_str());
  */
  cFail.Print((pngFileBase + TString("-fail.png")).Data());

  // Print fit outcome into fit log
  result->printStream(fitLog,RooPrintable::kValue,RooPrintable::kVerbose);
  fitLog << endl;
  printCorrelations(fitLog, result);
  
  
  return;
}

