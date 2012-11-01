#if !defined(__CINT__) || defined(__MAKECINT__)
#include "../Include/fitFunctionsCore.hh"
#include <RooEffProd.h>
#endif

int performWeightedFit=1;
int debugWeightedFit=0; // extra datasets are created if this is 1
int bkgPassContainsErrorFunction=0;
int bkgFailContainsErrorFunction=0;
bool unbinnedFitDefault=false;
bool unbinnedFit=unbinnedFitDefault;


// helping variables for fitMassWithTemplates
const int cLoadVariableLimitsFromFile=0;
const TString variableLimitsFName="user_fit_params.dat";
const TString variableLimitsFNameTempl="user_fit_params_wTemplates.dat";

struct RealLimit 
{
  std::string name;
  double av;
  double lo;
  double hi;

  void clear() { av=0.; lo=0.; hi=0.; }
};


inline int PosOk(size_t pos) { return (pos==std::string::npos) ? 0 : 1; }

// ---------------------------------------------------------

int loadVariableLimitsFromFile(int isReco, const TString &cut, 
			       int limitCount, RealLimit *limits, 
			       int fitWithTemplates) {
  unbinnedFit=unbinnedFitDefault;
  TString fname=
    (fitWithTemplates) ? variableLimitsFNameTempl : variableLimitsFName;
  std::ifstream fin(fname);
  if (!fin.is_open()) {
    std::cout << "loadVariableLimitsFromFile is on, but the file <" << 
      fname << "> could not be found\n";
    assert(0);
  }

  for (int i=0; i<limitCount; ++i) {
    limits[i].clear();
  }

  std::string s;
  int debugLoad=1;
  int stage=0;
  int found=0;
  while (!fin.eof() && getline(fin,s)) {
    if (debugLoad) std::cout << "got line <" << s << ">\n";
    if (s[0]=='#') {
      if (stage==2+limitCount) stage=0;
      if (debugLoad) std::cout << "... skipping '#'\n";
      continue;
    }
    if (s[0]=='%') {
      if (debugLoad) std::cout << "... file end symbol '%' detected\n";
      break;
    }
    if ((stage==0) && PosOk(s.find("isReco="))) {
      if (debugLoad) std::cout << "isReco detected\n";
      int isReco_in_file=atoi(s.c_str() + 7);
      if (debugLoad) std::cout << "isReco_in_file=" << isReco_in_file << "\n";
      if (isReco==isReco_in_file) {
	if (debugLoad) std::cout << "... next stage\n";
	stage++;
      }
      else {
	if (debugLoad) std::cout << "... not suitable\n";
      }
    }
    else if ((stage==1) && PosOk(s.find("cut="))) { 
      if (debugLoad) std::cout << "'cut' detected. Seeking {" << cut << "}\n";
      size_t pos1=s.find('{');
      size_t pos2=s.find('}');
      if (!PosOk(pos1) || !PosOk(pos2)) {
	std::cout << "could not extract cut from a line <" << s << ">: missing beginning '{' and/or ending '}'\n";
	assert(0);
      }
      TString cut_in_file=s.substr(pos1+1,pos2-pos1-1);
      if (debugLoad) std::cout << "cut in file={" << cut_in_file << "}\n";
      if (cut.Contains(cut_in_file)) {
	if (debugLoad) std::cout << "cut ok. Next stage\n";
	stage++;
      }
      else {
	if (debugLoad) std::cout << "... cut not suitable.Stage=0\n";
	stage=0;
      }
    }
    else if ((stage>=2) && (stage<2+limitCount)) {
      int idx=stage-2;
      std::stringstream ss(s);
      if (idx<3) {
	ss >> limits[idx].lo >> limits[idx].hi;
	if (debugLoad) std::cout << "limits[" << idx <<"].lo=" << limits[idx].lo << ", .hi=" << limits[idx].hi << "\n";
      }
      else {
	ss >> limits[idx].av >> limits[idx].lo >> limits[idx].hi;
	if (debugLoad) std::cout << "limits[" << idx <<"].av=" << limits[idx].av << ", .lo=" << limits[idx].lo << ", .hi=" << limits[idx].hi << "\n";
      }
      stage++;
    }
    else if (stage==2+limitCount) {
      found++;
      if (PosOk(s.find("unbinnedFit="))) {
	if (debugLoad) std::cout << "unbinnedFit key detected\n";
	size_t posMax=s.find("unbinnedFit=") + TString("unbinnedFit=").Length() + 3;
	size_t posYes=s.find("true");
	unbinnedFit=(posYes<posMax) ? true : false;
	if (debugLoad) std::cout << "setting unbinnedFit to " << ((unbinnedFit) ? "true" : "false") << "\n";
      }
      stage=0;
    }
    else {
      if (debugLoad) std::cout << "... skipping\n";
    }
  }
  fin.close();
  
  return found;
}

// ----------------------------------------------------

void printCorrelations(ostream& os, RooFitResult *res)
{
  if (!res) {
    os << "printCorrelations: fit result ptr is null" << endl;
    return;
  }
  ios_base::fmtflags flags = os.flags();
  const RooArgList *parlist = res->correlation("eff");
  
  os << "  Correlation Matrix" << endl;
  os << " --------------------" << endl;
  for(Int_t i=0; i<parlist->getSize(); i++) {
    for(Int_t j=0; j<parlist->getSize(); j++) 
      os << "  " << setw(7) << setprecision(4) << fixed 
	 << res->correlationMatrix()(i,j);    
    os << endl;
  }
  os.flags(flags);
}


void fitMass(TTree *passTree, TTree *failTree, TString cut, int mode, 
	     double &efficiency, double &efficiencyErrHi, 
	     double &efficiencyErrLo, 
	     TPad *passPad, TPad *failPad, TFile *plotsRootFile, 
	     ofstream &fitLog, int NsetBins, bool isRECO, 
	     const char* setBinsType, TString dirTag){

  // meaningless check, saving from compiler complaints
  if (dirTag.Length() && 0) 
    std::cout << "fitMass : dirTag=" << dirTag << "\n";

  RealLimit lims[19];

  for (int j=0; j<19; j++)
    {
      lims[j].av=0.0;
      lims[j].lo=0.0;
      lims[j].hi=0.0;
    }

    lims[0].name="mass";
    lims[1].name="et";
    lims[2].name="eta";
    lims[3].name="zMass";
    lims[4].name="zWidth";
    lims[5].name="nsignal";
    lims[6].name="eff";
    lims[7].name="lambdaBgPass";
    lims[8].name="cbMeanPass";
    lims[9].name="cbWidthPass";
    lims[10].name="cbAlphaPass";
    lims[11].name="cbNPass";
    lims[12].name="nbgPass";
    lims[13].name="lambdaBgFail";
    lims[14].name="cbMeanFail";
    lims[15].name="cbWidthFail";
    lims[16].name="cbAlphaFail";
    lims[17].name="cbNFail";
    lims[18].name="nbgFail";
  
  int foundOnFile=0;
  if (cLoadVariableLimitsFromFile) {
    foundOnFile=loadVariableLimitsFromFile(isRECO,cut,19,lims,0);
    std::cout << "loadVariableLimitsFromFile found " << foundOnFile << " records for cut={" << cut << "}\n";

    std::cout<<"\nbeg-llllllllllllllllllllllimssssssssssssssssssssssss"<<std::endl;  
    
    for (int j=0; j<12; j++)
      {
	std::cout<<lims[j].name<<": "<<lims[j].av<<", "<<lims[j].lo<<", "<<lims[j].hi<<std::endl;
      }
    
    std::cout<<"end-llllllllllllllllllllllimssssssssssssssssssssssss"<<std::endl;   
  }

  if (!foundOnFile) {
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
      //lims[10].av=5.0;  lims[10].lo=0.0;  lims[10].hi=20.0;
      lims[10].av=5.0;  lims[10].lo=-20.0;  lims[10].hi=20.0;
      lims[11].av=1.0;  lims[11].lo=0.0;  lims[11].hi=10.0;
      lims[12].av=1.0;  lims[12].lo=0.0;  lims[12].hi=1.0e5;
      lims[13].av=-0.1; lims[13].lo=-0.5; lims[13].hi=0.5;
      lims[14].av=0.0;  lims[14].lo=-5.0; lims[14].hi=5.0;
      lims[15].av=1.0;  lims[15].lo=0.1;  lims[15].hi=6.0;
      //lims[16].av=5.0;  lims[16].lo=0.0;  lims[16].hi=20.0;
      lims[16].av=5.0;  lims[16].lo=-20.0;  lims[16].hi=20.0;
      lims[17].av=1.0;  lims[17].lo=0.0;  lims[17].hi=10.0;
      lims[18].av=1.0;  lims[18].lo=0.0;  lims[18].hi=1.0e5;      
    }
  }
  
  // Define data sets
  RooRealVar mass("mass","mass",lims[0].lo, lims[0].hi);
  mass.setBins(NsetBins);
  RooRealVar et ("et" ,"et" ,lims[1].lo, lims[1].hi);
  RooRealVar eta("eta","eta",lims[2].lo, lims[2].hi);
  RooRealVar weight("weight","weight",0.,100.);
  RooFormulaVar rooCut("rooCut","rooCut",cut,RooArgSet(et,eta));
  RooArgSet dsetArgs(mass,et,eta);
  if (performWeightedFit && (passTree->GetBranch("weight")==NULL)) {
    std::cout << "fitMass: perforkWeightedFit=1, "
	      << "but the tree has no branch \"weight\"\n";
  }
  RooDataSet  *dataUnbinnedPass = NULL;
  RooDataSet  *dataUnbinnedFail = NULL;
  RooDataSet *dataUnbinnedPassUnweighted= (!debugWeightedFit) ? NULL :
    new RooDataSet("dataUnbinnedPass","dataUnbinnedPass",
		   passTree,dsetArgs, rooCut);
  RooDataSet *dataUnbinnedFailUnweighted= (!debugWeightedFit) ? NULL :
    new RooDataSet("dataUnbinnedFail","dataUnbinnedFail",
		   failTree,dsetArgs, rooCut);
  std::string dline(70,'-'); dline+='\n';

  if (!performWeightedFit) {
    dataUnbinnedPass = 
      new RooDataSet("dataUnbinnedPass","dataUnbinnedPass",
		     passTree,dsetArgs, rooCut);
    dataUnbinnedFail = 
      new RooDataSet("dataUnbinnedFail","dataUnbinnedFail",
		     failTree,dsetArgs, rooCut);
    std::cout << "\n\n\tunweighted Fit" << std::endl;
  }
  else {
    std::cout << "\n\n\tweighted Fit\n\n";
    fitLog << "\n\n\tweighted Fit\n\n";
    dsetArgs.add(weight);
    dataUnbinnedPass = 
      new RooDataSet("dataUnbinnedPass","dataUnbinnedPass",
		     passTree,dsetArgs, rooCut, "weight");
    dataUnbinnedFail = 
      new RooDataSet("dataUnbinnedFail","dataUnbinnedFail",
		     failTree,dsetArgs, rooCut, "weight");

    std::cout << dline; dataUnbinnedPass->Print("V"); std::cout << dline;
    std::cout << dline; dataUnbinnedFail->Print("V"); std::cout << dline;
    std::cout << std::endl;
  }

  RooDataHist *dataBinnedPass   = (unbinnedFit) ? NULL :
    dataUnbinnedPass->binnedClone("dataBinnedPass","dataBinnedPass");
  RooDataHist *dataBinnedFail   = (unbinnedFit) ? NULL :
    dataUnbinnedFail->binnedClone("dataBinnedFail","dataBinnedFail");
  RooDataHist *dataBinnedPassUnweighted   = (!debugWeightedFit) ? NULL :
    dataUnbinnedPassUnweighted->binnedClone("dataBinnedPassNoW","dataBinnedPassNoW");
  RooDataHist *dataBinnedFailUnweighted   = (!debugWeightedFit) ? NULL :
    dataUnbinnedFailUnweighted->binnedClone("dataBinnedFailNoW","dataBinnedFailNoW");
  RooCategory probeType("probeType","probeType");
  probeType.defineType("pass");
  probeType.defineType("fail");
  RooAbsData *data;

  // If needed do binned fit
  //bool unbinnedFit = false;
  if(unbinnedFit){
    RooArgSet combiDSetArgs(mass);
    if (performWeightedFit) {
      combiDSetArgs.add(weight);
      data = new RooDataSet("data","data",combiDSetArgs,Index(probeType),
			    Import("pass",*dataUnbinnedPass), 
			    Import("fail",*dataUnbinnedFail),
			    RooFit::WeightVar("weight"));
    }
    else {
      data = new RooDataSet("data","data",combiDSetArgs,Index(probeType),
			    Import("pass",*dataUnbinnedPass), 
			    Import("fail",*dataUnbinnedFail));
    }
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
  if( mode == DYTools::COUNTnFIT ){
    RooGenericPdf * simpleSignal = new RooGenericPdf("simpleSignal","simpleSignal",
						     "1.0",RooArgList());
    RooExtendPdf * simpleSignalExtended = new RooExtendPdf("passPdf", "passPdf",
							   *simpleSignal, nsigPass);
    passPdf = simpleSignalExtended;
  }else if( mode == DYTools::FITnFIT ){
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
  if( mode == DYTools::FITnFIT ){
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

  /*
 [#0] WARNING:InputArguments -- RooAbsPdf::fitTo(fullPdf) WARNING: a likelihood 
 fit is request of what appears to be weighted data.
 While the estimated values of the parameters will always be calculated taking 
 the weights into account,
 there are multiple ways to estimate the errors on these parameter values. You 
 are advised to make an explicit choice on the error calculation:
 - Either provide SumW2Error(kTRUE), to calculate a sum-of-weights corrected 
 HESSE error matrix (error will be proportional to the number of events)
 - Or provide SumW2Error(kFALSE), to return errors from original HESSE error 
 matrix (which will be proportional to the sum of the weights)
 If you want the errors to reflect the information contained in the provided 
 dataset, choose kTRUE.
 If you want the errors to reflect the precision you would be able to obtain 
 with an unweighted dataset with 'sum-of-weights' events, choose kFALSE.
  */

  std::cout << "fitting \n\n";
  Bool_t sumw2error= kFALSE;
  RooFitResult *result = 
    fullPdf.fitTo(*data,
		  Extended(kTRUE),
		  Save(),
		  NumCPU(2,true)
		  , RooFit::SumW2Error(sumw2error)
		  );

  // Release shape parameters and refine the fit
  if( mode == DYTools::FITnFIT ){
    cbAlphaPass.setConstant(kFALSE);
    cbNPass    .setConstant(kFALSE);
  }
  cbAlphaFail.setConstant(kFALSE);
  cbNFail    .setConstant(kFALSE);
  result = fullPdf.fitTo(*data,
			 Extended(kTRUE),
			 Minos(RooArgSet(eff)),
			 Save(),
			 NumCPU(2,true)
			 , RooFit::SumW2Error(sumw2error)
			 );
  // If minos fails, refit without minos
  if((fabs(eff.getErrorLo())<5e-5) || (eff.getErrorHi()<5e-5) || 
     (fabs(eff.getAsymErrorLo())<5e-5) || (eff.getAsymErrorHi()<5e-5)) {
    cout << "MINOS FAILS" << endl;
    fitLog << "MINOS FAILS" << endl;
    result = fullPdf.fitTo(*data,
			   Extended(kTRUE),
			   Save(),
			   NumCPU(2,true)
			   , RooFit::SumW2Error(sumw2error)
			   );
    if((fabs(eff.getErrorLo())<5e-5) || (eff.getErrorHi()<5e-5) || 
       (fabs(eff.getAsymErrorLo())<5e-5) || (eff.getAsymErrorHi()<5e-5)) {
       cout << "SECOND FIT FAILURE" << endl;
       fitLog << "SECOND FIT FAILURE" << endl;
     }
  }

  efficiency       = eff.getVal();
  efficiencyErrHi  = eff.getErrorHi();
  efficiencyErrLo  = fabs(eff.getErrorLo());

  // Draw fit results
  passPad->cd();
  passPad->Clear();
  RooPlot *framePass = mass.frame();

  if (unbinnedFit) dataUnbinnedPass->plotOn(framePass);
  else dataBinnedPass->plotOn(framePass);

  if(mode == DYTools::FITnFIT){
    passPdf->plotOn(framePass);
    passPdf->plotOn(framePass,Components("bgPassPdf"),LineStyle(kDashed));
  }
  framePass->Draw();
  passPad->Update();
  if (plotsRootFile) { plotsRootFile->cd(); passPad->Write(); }


  failPad->cd();
  failPad->Clear();
  RooPlot *frameFail = mass.frame();

  if (unbinnedFit) dataUnbinnedFail->plotOn(frameFail);
  else dataBinnedFail->plotOn(frameFail);
  
  failPdf->plotOn(frameFail);
  failPdf->plotOn(frameFail,Components("bgFailPdf"),LineStyle(kDashed));
  frameFail->Draw();
  failPad->Update();
  if (plotsRootFile) { plotsRootFile->cd(); failPad->Write(); }


  // Print fit outcome into fit log
  result->printStream(fitLog,RooPrintable::kValue,RooPrintable::kVerbose);
  fitLog << endl;
  printCorrelations(fitLog, result);


 // Clean-up

  if (dataUnbinnedPass) delete dataUnbinnedPass;
  if (dataUnbinnedPassUnweighted) delete dataUnbinnedPassUnweighted;
  if (dataUnbinnedFail) delete dataUnbinnedFail;
  if (dataUnbinnedFailUnweighted) delete dataUnbinnedFailUnweighted;
  if (dataBinnedPass) delete dataBinnedPass;
  if (dataBinnedFail) delete dataBinnedFail;
  if (dataBinnedPassUnweighted) delete dataBinnedPassUnweighted;
  if (dataBinnedFailUnweighted) delete dataBinnedFailUnweighted;
  if (data) delete data;
  //if (simpleSignal) delete simpleSignal;
  //if (simpleSignalExtended) delete simpleSignalExtended;
  if (passPdf) delete passPdf;
  if (failPdf) delete failPdf;
  //if (framePass) delete framePass;
  //if (frameFail) delete frameFail;
  if (result) delete result;
  
  return;
}


void fitMassWithTemplates(TTree *passTree, TTree *failTree, TString cut, 
			  int mode, double &efficiency, 
			  double &efficiencyErrHi, double &efficiencyErrLo,
			  TPad *passPad, TPad *failPad, TFile *plotsRootFile,
			  ofstream &fitLog, 
			  TH1F *templatePass, TH1F *templateFail, bool isRECO,
			  const char* setBinsType, TString dirTag, 
			  const TString &picFileExtraTag){

    
  RealLimit lims[12];
  
  for (int j=0; j<12; j++)
    {
      lims[j].av=0.0;
      lims[j].lo=0.0;
      lims[j].hi=0.0;
    }

    lims[0].name="mass";
    lims[1].name="et";
    lims[2].name="eta";
    lims[3].name="nsignal";
    lims[4].name="eff";
    lims[5].name="lambdaBgPass";
    lims[6].name="resMeanPass";
    lims[7].name="resSigma";
    lims[8].name="nbgPass";
    lims[9].name="lambdaBgFail";
    lims[10].name="resMeanFail";
    lims[11].name="nbgFail";
  
  int foundOnFile=0;
  if (cLoadVariableLimitsFromFile) {
    foundOnFile=loadVariableLimitsFromFile(isRECO,cut,12,lims,1);
    std::cout << "loadVariableLimitsFromFile found " << foundOnFile << " records for cut={" << cut << "}\n";

    std::cout<<"\nbeg-llllllllllllllllllllllimssssssssssssssssssssssss"<<std::endl;  
    
    for (int j=0; j<12; j++)
      {
	std::cout<<lims[j].name<<": "<<lims[j].av<<", "<<lims[j].lo<<", "<<lims[j].hi<<std::endl;
      }
    
    std::cout<<"end-llllllllllllllllllllllimssssssssssssssssssssssss"<<std::endl;   
  }


  if (!foundOnFile) {
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
      lims[7].av=1.5;   lims[7].lo=1.0;   lims[7].hi=6.0; // Low limit at 1 GeV to ensure smoothness of binned templates
      lims[8].av=1.0;   lims[8].lo=0.0;   lims[8].hi=1.0e5;
      lims[9].av=-0.1;  lims[9].lo=-0.5;  lims[9].hi=0.5;
      //lims[10].av=0.0;  lims[10].lo=-3.0; lims[10].hi=3.0; 
      lims[10].av=0.0;  lims[10].lo=-5.0; lims[10].hi=5.0; // 2012.07.10
      lims[11].av=1.0;  lims[11].lo=0.0;  lims[11].hi=1.0e5;
      
    }
  }
  

  // Define data sets
  RooRealVar mass("mass","mass",lims[0].lo, lims[0].hi);
  mass.setBins(30);
  RooRealVar et ("et" ,"et" ,lims[1].lo, lims[1].hi);
  RooRealVar eta("eta","eta",lims[2].lo, lims[2].hi);
  RooRealVar weight("weight","weight",0.,100.);
  RooFormulaVar rooCut("rooCut","rooCut",cut,RooArgSet(et,eta));
  RooArgSet dsetArgs(mass,et,eta);
  if (performWeightedFit && (passTree->GetBranch("weight")==NULL)) {
    std::cout << "fitMassWithTemplates: perforkWeightedFit=1, "
	      << "but the tree has no branch \"weight\"\n";
  }
  RooDataSet  *dataUnbinnedPass = NULL;
  RooDataSet  *dataUnbinnedFail = NULL;
  std::string dline(70,'-'); dline+='\n';

  RooDataSet *dataUnbinnedPassNoWeight = (!debugWeightedFit) ? NULL :
    new RooDataSet("dataUnbinnedPassNoWeight","dataUnbinnedPassNoWeight",
		   passTree,dsetArgs, rooCut);
  RooDataSet *dataUnbinnedFailNoWeight = (!debugWeightedFit) ? NULL :
    new RooDataSet("dataUnbinnedFailNoWeight","dataUnbinnedFailNoWeight",
		   failTree,dsetArgs, rooCut);
 
  if (!performWeightedFit) {
    dataUnbinnedPass = 
      new RooDataSet("dataUnbinnedPass","dataUnbinnedPass",
		     passTree,dsetArgs, rooCut);
    dataUnbinnedFail = 
      new RooDataSet("dataUnbinnedFail","dataUnbinnedFail",
		     failTree,dsetArgs, rooCut);
    std::cout << "\n\n\tunweighted Fit (" << dataUnbinnedPass->numEntries() 
	      << "p," << dataUnbinnedFail->numEntries() << "f)\n\n";
  }
  else {
    std::cout << "\n\n\tweighted Fit\n\n";
    fitLog << "\n\n\tweighted Fit\n\n";
    dsetArgs.add(weight);
    dataUnbinnedPass = 
      new RooDataSet("dataUnbinnedPass","dataUnbinnedPass",
		     passTree,dsetArgs, rooCut, "weight");
    dataUnbinnedFail = 
      new RooDataSet("dataUnbinnedFail","dataUnbinnedFail",
		     failTree,dsetArgs, rooCut, "weight");

    std::cout << dline; dataUnbinnedPass->Print("V"); std::cout << dline;
    std::cout << dline; dataUnbinnedFail->Print("V"); std::cout << dline;
    std::cout << std::endl;
  }
 
  RooDataHist *dataBinnedPass   = (unbinnedFit) ? NULL :
    dataUnbinnedPass->binnedClone("dataBinnedPass","dataBinnedPass");
  RooDataHist *dataBinnedFail   = (unbinnedFit) ? NULL :
    dataUnbinnedFail->binnedClone("dataBinnedFail","dataBinnedFail");
  RooDataHist *dataBinnedPassNoWeight   = (!debugWeightedFit) ? NULL :
    dataUnbinnedPassNoWeight->binnedClone("dataBinnedPassNoWeight","dataBinnedPassNoWeight");
  RooDataHist *dataBinnedFailNoWeight   = (!debugWeightedFit) ? NULL :
    dataUnbinnedFailNoWeight->binnedClone("dataBinnedFailNoWeight","dataBinnedFailNoWeight");
  RooCategory probeType("probeType","probeType");
  probeType.defineType("pass");
  probeType.defineType("fail");
  RooAbsData *data;
  
  if (0 && debugWeightedFit) {
    // test pass
    TCanvas *cTest=new TCanvas("cTestPass","cTestPass",800,800);
    cTest->Divide(2,2);
    RooPlot *frameTest1 = mass.frame(RooFit::Title("unB W pass (r) vs unB noW pass (v)"));
    dataUnbinnedPass->plotOn(frameTest1,RooFit::LineColor(kRed),RooFit::MarkerColor(kRed));
    dataUnbinnedPassNoWeight->plotOn(frameTest1,RooFit::LineColor(kViolet),RooFit::MarkerColor(kViolet),RooFit::MarkerStyle(24));
    cTest->cd(1);
    frameTest1->Draw();
    RooPlot *frameTest2 = mass.frame(RooFit::Title("unB W pass (r) vs B W pass (g)"));
    dataUnbinnedPass->plotOn(frameTest2,RooFit::LineColor(kRed),RooFit::MarkerColor(kRed));
    dataBinnedPass->plotOn(frameTest2,RooFit::LineColor(kGreen+1),RooFit::MarkerColor(kGreen+1),RooFit::MarkerStyle(24));
    cTest->cd(2);
    frameTest2->Draw();
    RooPlot *frameTest3 = mass.frame(RooFit::Title("B noW pass (b) vs unB noW pass (v)"));
    dataBinnedPassNoWeight->plotOn(frameTest3,RooFit::LineColor(kBlue+1),RooFit::MarkerColor(kBlue+1));
    dataUnbinnedPassNoWeight->plotOn(frameTest3,RooFit::LineColor(kViolet),RooFit::MarkerColor(kViolet),RooFit::MarkerStyle(24));
    cTest->cd(3);
    frameTest3->Draw();
    RooPlot *frameTest4 = mass.frame(RooFit::Title("B noW pass (b) vs B W pass (g)"));
    dataBinnedPassNoWeight->plotOn(frameTest4,RooFit::LineColor(kBlue+1),RooFit::MarkerColor(kBlue+1));
    dataBinnedPass->plotOn(frameTest4,RooFit::LineColor(kGreen+1),RooFit::MarkerColor(kGreen+1),RooFit::MarkerStyle(24));
    cTest->cd(4);
    frameTest4->Draw();
    cTest->Update();
    return;
  } // debugWeightedFit

  if (0 && debugWeightedFit) {
    // test fail
    TCanvas *cTest=new TCanvas("cTestFail","cTestFail",800,800);
    cTest->Divide(2,2);
    RooPlot *frameTest1 = mass.frame(RooFit::Title("unB W fail (r) vs unB noW fail (v)"));
    dataUnbinnedFail->plotOn(frameTest1,RooFit::LineColor(kRed),RooFit::MarkerColor(kRed));
    dataUnbinnedFailNoWeight->plotOn(frameTest1,RooFit::LineColor(kViolet),RooFit::MarkerColor(kViolet),RooFit::MarkerStyle(24));
    cTest->cd(1);
    frameTest1->Draw();
    RooPlot *frameTest2 = mass.frame(RooFit::Title("unB W fail (r) vs B W fail (g)"));
    dataUnbinnedFail->plotOn(frameTest2,RooFit::LineColor(kRed),RooFit::MarkerColor(kRed));
    dataBinnedFail->plotOn(frameTest2,RooFit::LineColor(kGreen+1),RooFit::MarkerColor(kGreen+1),RooFit::MarkerStyle(24));
    cTest->cd(2);
    frameTest2->Draw();
    RooPlot *frameTest3 = mass.frame(RooFit::Title("B noW fail (b) vs unB noW fail (v)"));
    dataBinnedFailNoWeight->plotOn(frameTest3,RooFit::LineColor(kBlue+1),RooFit::MarkerColor(kBlue+1));
    dataUnbinnedFailNoWeight->plotOn(frameTest3,RooFit::LineColor(kViolet),RooFit::MarkerColor(kViolet),RooFit::MarkerStyle(24));
    cTest->cd(3);
    frameTest3->Draw();
    RooPlot *frameTest4 = mass.frame(RooFit::Title("B noW fail (b) vs B W fail (g)"));
    dataBinnedFailNoWeight->plotOn(frameTest4,RooFit::LineColor(kBlue+1),RooFit::MarkerColor(kBlue+1));
    dataBinnedFail->plotOn(frameTest4,RooFit::LineColor(kGreen+1),RooFit::MarkerColor(kGreen+1),RooFit::MarkerStyle(24));
    cTest->cd(4);
    frameTest4->Draw();
    cTest->Update();
    return;
  } // debugWeightedFit


  // If needed do binned fit
  //bool unbinnedFit = true;
  if(unbinnedFit){
    RooArgSet combiDSetArgs(mass);
    if (performWeightedFit) {
      combiDSetArgs.add(weight);
      data = new RooDataSet("data","data",combiDSetArgs,Index(probeType),
			    Import("pass",*dataUnbinnedPass), 
			    Import("fail",*dataUnbinnedFail),
			    RooFit::WeightVar("weight"));
    }
    else {
      data = new RooDataSet("data","data",combiDSetArgs,Index(probeType),
			    Import("pass",*dataUnbinnedPass), 
			    Import("fail",*dataUnbinnedFail));
    }
    std::cout << dline; data->Print("V"); std::cout << dline;
    cout << endl << "Setting up UNBINNED fit" << endl << endl;
  }else{
    data = new RooDataHist("data","data",mass,Index(probeType),
			   Import("pass",*dataBinnedPass), 
			   Import("fail",*dataBinnedFail));    
    cout << endl << "Setting up BINNED fit" << endl << endl;
  }

  if (0 && debugWeightedFit) {
    // test combined dset
    RooDataSet *dataNoWeight = 
      new RooDataSet("dataNoWeight","dataNoWeight",
		     RooArgSet(mass),Index(probeType),
		     Import("pass",*dataUnbinnedPassNoWeight), 
		     Import("fail",*dataUnbinnedFailNoWeight));
    RooDataHist *dataBinnedNoWeight =
      new RooDataHist("dataBinnedNoWeight","dataBinnedNoWeight",
		     RooArgSet(mass),Index(probeType),
		     Import("pass",*dataBinnedPassNoWeight), 
		     Import("fail",*dataBinnedFailNoWeight));

    TCanvas *cTest=new TCanvas("cTestCombi","cTestCombi",800,400);
    cTest->Divide(2,1);
    TString title=((unbinnedFit) ? "UNBINNED ":"BINNED ");
    RooPlot *frameTest1 = mass.frame(RooFit::Title(title + TString("W pass (r) vs noW pass (v)")));
    data->plotOn(frameTest1,RooFit::Cut("probeType==probeType::pass"),RooFit::LineColor(kRed),RooFit::MarkerColor(kRed));
    if (unbinnedFit)
      dataNoWeight->plotOn(frameTest1,RooFit::Cut("probeType==probeType::pass"),RooFit::LineColor(kViolet),RooFit::MarkerColor(kViolet),RooFit::MarkerStyle(24));
    else
      dataBinnedNoWeight->plotOn(frameTest1,RooFit::Cut("probeType==probeType::pass"),RooFit::LineColor(kViolet),RooFit::MarkerColor(kViolet),RooFit::MarkerStyle(24));
    cTest->cd(1);
    frameTest1->Draw();

    RooPlot *frameTest2 = mass.frame(RooFit::Title(title + TString("W fail (r) vs noW fail (v)")));
    data->plotOn(frameTest2,RooFit::Cut("probeType==probeType::fail"),RooFit::LineColor(kRed),RooFit::MarkerColor(kRed));
    if (unbinnedFit) 
      dataNoWeight->plotOn(frameTest2,RooFit::Cut("probeType==probeType::fail"),RooFit::LineColor(kViolet),RooFit::MarkerColor(kViolet),RooFit::MarkerStyle(24));

    else 
      dataBinnedNoWeight->plotOn(frameTest2,RooFit::Cut("probeType==probeType::fail"),RooFit::LineColor(kViolet),RooFit::MarkerColor(kViolet),RooFit::MarkerStyle(24));
    cTest->cd(2);
    frameTest2->Draw();
    cTest->Update();
    if (dataNoWeight) delete dataNoWeight;
    if (dataBinnedNoWeight) delete dataBinnedNoWeight;
    return;
  } // debugWeightedFit

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
  if( mode == DYTools::COUNTnFIT ){
    RooGenericPdf * simpleSignal = new RooGenericPdf("simpleSignal","simpleSignal",
						     "1.0",RooArgList());
    RooExtendPdf * simpleSignalExtended = new RooExtendPdf("passPdf", "passPdf",
							   *simpleSignal, nsigPass);
    passPdf = simpleSignalExtended;
  }else if( mode == DYTools::FITnFIT ){
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

  // The old resMeanFail is commented out below. The new one
  // is constructed as resMeanPass+resMeanFailShift so that one can
  // constrain resMeanFail = resMeanPass by fixing the "shift" to zero.
  // This is useful when the failing probes distribution is not sensitive to
  // Gaussian smearing parameters.
  RooRealVar resMeanFailShift("resMeanFailShift","resMeanFailShift", 0.0, 2*lims[10].lo, 2*lims[10].hi);
  RooFormulaVar resMeanFail("resMeanFail", "@0+@1", RooArgList(resMeanPass, resMeanFailShift));
  //   RooRealVar resMeanFail("resMeanFail","cbMeanFail"   ,lims[10].av, lims[10].lo, lims[10].hi);

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

  // Is this a low momentum case of the RECO efficiency? 
  if(isRECO){
    printf("\n Fitting RECO efficiency. Check if special handing for low Et is needed.\n");

    // Establish whether the present Et/Eta bin contains
    // any probes below 30 GeV
    TString lowEtSubsetCut = cut;
    lowEtSubsetCut += Form(" && (%s<30.0)", et.GetName());
    RooDataSet *dataPassLowEt = new RooDataSet("dataPassLowEt", "dataPassLowEt", passTree, 
					       dsetArgs, lowEtSubsetCut, "weight");
    RooDataSet *dataFailLowEt = new RooDataSet("dataFailLowEt", "dataFailLowEt", failTree, 
					       dsetArgs, lowEtSubsetCut, "weight");
    double lowEtEntries = dataPassLowEt->sumEntries() + dataFailLowEt->sumEntries();

    double allEntries = data->sumEntries();
    printf("Found all entries           : %f\n", allEntries);
    printf("Found entries with Et<30 GeV: %f\n", lowEtEntries);

    // Modify the parameter constraints for low Et
    if( lowEtEntries > 1 ) {
      // For low Et data for the SC->GSF efficiency, constrain the mean
      // of the Gaussian smearing for the failed probes to be equal to the
      // mean for the passing probes.
      printf("\n The mean of the smearing for the failing probes is set equal\n");
      printf(" to the mean of the smearing for the passing probes.\n\n"); 
      resMeanFailShift.setVal(0.0);
      resMeanFailShift.setConstant(kTRUE);
    }
  }
  
  
  // Combine pass and fail
  RooSimultaneous fullPdf("fullPdf","fullPdf",probeType);
  fullPdf.addPdf(*passPdf,"pass");
  fullPdf.addPdf(*failPdf,"fail");

  
  // Do the fit
  double total = dataUnbinnedPass->numEntries() + dataUnbinnedFail->numEntries();
  nsignal.setVal( 0.99*total);
  eff.setVal(0.90);
  if( mode == DYTools::FITnFIT ){
    nbgPass.setVal(0.01*total);
  }
  nbgFail.setVal(0.01*total);

  // Print the initial state of the PDF
  printf("\nThe full PDF parameter printout:\n");
  fullPdf.getVariables()->Print("v");
  /*
 [#0] WARNING:InputArguments -- RooAbsPdf::fitTo(fullPdf) WARNING: a likelihood 
 fit is request of what appears to be weighted data.
 While the estimated values of the parameters will always be calculated taking 
 the weights into account,
 there are multiple ways to estimate the errors on these parameter values. You 
 are advised to make an explicit choice on the error calculation:
 - Either provide SumW2Error(kTRUE), to calculate a sum-of-weights corrected 
 HESSE error matrix (error will be proportional to the number of events)
 - Or provide SumW2Error(kFALSE), to return errors from original HESSE error 
 matrix (which will be proportional to the sum of the weights)
 If you want the errors to reflect the information contained in the provided 
 dataset, choose kTRUE.
 If you want the errors to reflect the precision you would be able to obtain 
 with an unweighted dataset with 'sum-of-weights' events, choose kFALSE.
  */

  std::cout << "fitting \n\n";
  Bool_t sumw2error= kFALSE;
  RooFitResult *result=
    fullPdf.fitTo(*data, Extended(kTRUE), 
		  Minos(RooArgSet(eff)), Save(), NumCPU(2,true)
		  , RooFit::SumW2Error(sumw2error)
		  );
  
  
  
  // If minos fails, check if we can remove the background component in 
  // the passing sample
  if((fabs(eff.getErrorLo())<5e-5) || (eff.getErrorHi()<5e-5) || 
     (fabs(eff.getAsymErrorLo())<5e-5) || (eff.getAsymErrorHi()<5e-5)) {
    if (nbgPass.getVal()<1e-5) {
      cout << "\n\n\tMINOS FAILS (removing nbgPass) eff=" 
	   << eff.getVal() << "\n\n" << endl;
      fitLog << "\n\n\tMINOS FAILS (removing nbgPass) eff=" 
	     << eff.getVal() << "\n\n" << endl;
      nbgPass.setVal(0.);
      nbgPass.setConstant(kTRUE);
      lambdaBgPass.setConstant(kTRUE);

      result=
	fullPdf.fitTo(*data, Extended(kTRUE), 
		      Minos(RooArgSet(eff)), Save(), NumCPU(2,true)
		      , RooFit::SumW2Error(sumw2error)
		      );
    }
  }

 // If minos fails, refit without minos
  if((fabs(eff.getErrorLo())<5e-5) || (eff.getErrorHi()<5e-5) || 
     (fabs(eff.getAsymErrorLo())<5e-5) || (eff.getAsymErrorHi()<5e-5)) {
    cout << "\n\n\tMINOS FAILS\n\n" << endl;
    fitLog << "\n\n\tMINOS FAILS\n\n" << endl;
    result = fullPdf.fitTo(*data, Extended(kTRUE), Save(), NumCPU(2,true)
			   , RooFit::SumW2Error(sumw2error)
			   );
    if((fabs(eff.getErrorLo())<5e-5) || (eff.getErrorHi()<5e-5)) {
       cout << "\n\n\tSECOND FIT FAILURE\n\n" << endl;
       fitLog << "\n\n\tSECOND FIT FAILURE\n\n" << endl;
     }
  }
  // After the fit, print all info about parameters, constant and floating,
  // initial and final.
  result->Print("v");
  
  efficiency     = eff.getVal();
  efficiencyErrHi  = eff.getErrorHi();
  efficiencyErrLo  = fabs(eff.getErrorLo());

  // Draw fit results
  passPad->cd();
  passPad->Clear();
  RooPlot *framePass = mass.frame();

  if (unbinnedFit) dataUnbinnedPass->plotOn(framePass);
  else dataBinnedPass->plotOn(framePass);

  if(mode == DYTools::FITnFIT){
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

  TString pngFilePath=TString("../root_files/tag_and_probe/") + dirTag + TString("/");
  TString pngFileBase="fit-";
  if (performWeightedFit) pngFileBase.Append("weighted-");
  if (unbinnedFit) pngFileBase.Append("unbinned-"); else pngFileBase.Append("binned-");
  pngFileBase.Append(picFileExtraTag);
  if (isRECO) pngFileBase+="-reco-"; else pngFileBase+="-id-";
  pngFileBase += cutF;
  
  TString canvasName="canv-" + pngFileBase + TString("-pass");
  TCanvas cPass(canvasName,canvasName);
  framePass->Draw();
  cPass.Update();
  if (plotsRootFile) { plotsRootFile->cd(); cPass.Write(); }

  /*
  if (isRECO) cPass.Print(("../root_files/tag_and_probe/"+(std::string)dirTag+"/fit-reco-"+(std::string)cutF+"-pass.png").c_str());
  else cPass.Print(("../root_files/tag_and_probe/"+(std::string)dirTag+"/fit-id-"+(std::string)cutF+"-pass.png").c_str());
  */
  cPass.Print((pngFilePath + pngFileBase + TString("-pass.png")).Data());

  failPad->cd();
  failPad->Clear();
  RooPlot *frameFail = mass.frame();
  dataUnbinnedFail->plotOn(frameFail);
  failPdf->plotOn(frameFail);
  failPdf->plotOn(frameFail,Components("bgFailPdf"),LineStyle(kDashed));
  frameFail->Draw();
  failPad->Update();

  canvasName="canv-" + pngFileBase + TString("-fail");
  TCanvas cFail(canvasName,canvasName);
  frameFail->Draw();
  cFail.Update();
  if (plotsRootFile) { plotsRootFile->cd(); cFail.Write(); }

  cFail.Print((pngFilePath + pngFileBase + TString("-fail.png")).Data());

  // Print fit outcome into fit log
  result->printStream(fitLog,RooPrintable::kValue,RooPrintable::kVerbose);
  fitLog << endl;
  printCorrelations(fitLog, result);
  
 
  // Clean-up

  if (dataUnbinnedPass) delete dataUnbinnedPass;
  if (dataUnbinnedFail) delete dataUnbinnedFail;
  if (dataBinnedPass) delete dataBinnedPass;
  if (dataBinnedFail) delete dataBinnedFail;
  if (data) delete data;
  //if (simpleSignal) delete simpleSignal;
  //if (simpleSignalExtended) delete simpleSignalExtended;
  if (passPdf) delete passPdf;
  if (failPdf) delete failPdf;
  //if (framePass) delete framePass;
  //if (frameFail) delete frameFail;
  if (result) delete result;

  if (dataUnbinnedPassNoWeight) delete dataUnbinnedPassNoWeight;
  if (dataUnbinnedFailNoWeight) delete dataUnbinnedFailNoWeight;
  if (dataBinnedPassNoWeight) delete dataBinnedPassNoWeight;
  if (dataBinnedFailNoWeight) delete dataBinnedFailNoWeight;

  return;
}

