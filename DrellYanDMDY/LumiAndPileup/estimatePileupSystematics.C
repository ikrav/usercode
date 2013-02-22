{

  // This scriptestimates the overall pile-up systematics
  // by comparing the fully corrected pre-FSR r-shape obtained
  // with the standard PU reweighting and the reweightings
  // where inelastic cross section is faried +-5%.
  //
  // Additionally, the event scale factors (TnP output)
  // and the signal MC full event efficiency are compared
  // for the three reweightings.
  //
  // Note: it appears that a significant fraction of the effect
  // on the r-shape comes from the pile-up effects in detector
  // resolution unfolding (it is not printed here by itself,
  // but the r-shape variation overall, of coruse, includes it).
  //
  
  TString dirDefault = "/home/hep/ikrav/releases/another_UserCode_v2/UserCode/ikravchenko/DrellYanDMDY/";
  TString dirPlus5   = "/home/hep/ikrav/releases/another_UserCode_v2_pileup_syst/UserCode/ikravchenko/DrellYanDMDY/";
  TString dirMinus5   = "/home/hep/ikrav/releases/another_UserCode_v2_pileup_syst_minus/UserCode/ikravchenko/DrellYanDMDY/";

  TFile *feffdef = new TFile(dirDefault 
			     + TString("root_files/constants/DY_m10+pr+a05+o03+pr_4839pb/")
			     + TString("event_efficiency_constants1D.root") );
  TMatrixD *effdef = (TMatrixD*)feffdef->Get("efficiencyArray");
  
  TFile *feffplus5percent = new TFile(dirPlus5
				      + TString("root_files/constants/DY_m10+pr+a05+o03+pr_4839pb/")
				      + TString("event_efficiency_constants1D.root"));
  TMatrixD *effplus5percent = (TMatrixD*)feffplus5percent->Get("efficiencyArray");
  
  TFile *feffminus5percent = new TFile(dirMinus5
				      + TString("root_files/constants/DY_m10+pr+a05+o03+pr_4839pb/")
				      + TString("event_efficiency_constants1D.root"));
  TMatrixD *effminus5percent = (TMatrixD*)feffminus5percent->Get("efficiencyArray");
  
  printf("Efficiency:\n");
  printf("    def         with+5%%      with-5%%   diff with+5%%, %% rel   diff with+5%%, %% rel \n");
  double diffEff[40];
  for(int i=0; i<40; i++){
    double diff1 = fabs( ((*effdef)(i,0) - (*effplus5percent)(i,0)) )/ (*effdef)(i,0);
    double diff2 = fabs( ((*effdef)(i,0) - (*effminus5percent)(i,0)) )/ (*effdef)(i,0);
    printf("%2d  %6.4f       %6.4f       %6.4f          %6.2f          %6.2f\n", i, (*effdef)(i,0), 
	   (*effplus5percent)(i,0), (*effminus5percent)(i,0), 
	   diff1*100, diff2*100);
    diffEff[i] = diff1;
    if( diff2 > diff1 ) diffEff[i] = diff2;
  }
  
  // Scale factors

  TFile *fscdef = new TFile(dirDefault
			    + TString("root_files/constants/DY_m10+pr+a05+o03+pr_4839pb/")
			    + TString("scale_factors_1D_Full2011_hltEffOld_PU.root") );
  TMatrixD *scdef = (TMatrixD*)fscdef->Get("scaleFactor");

  TFile *fscplus5percent = new TFile(dirPlus5
				     + TString("root_files/constants/DY_m10+pr+a05+o03+pr_4839pb/")
				     + TString("scale_factors_1D_Full2011_hltEffOld_PU.root") );
  TMatrixD *scplus5percent = (TMatrixD*)fscplus5percent->Get("scaleFactor");

  TFile *fscminus5percent = new TFile(dirMinus5
				     + TString("root_files/constants/DY_m10+pr+a05+o03+pr_4839pb/")
				     + TString("scale_factors_1D_Full2011_hltEffOld_PU.root") );
  TMatrixD *scminus5percent = (TMatrixD*)fscminus5percent->Get("scaleFactor");

  double diffSF[40];
  printf("\nScale Factor\n"); 
  printf("bin      def        with+5%%      with-5%%   diff with+5%%, %% rel   diff with-5%%, %% rel\n");
  for(int i=0; i<40; i++){
    double diff1 = fabs( ((*scdef)(i,0) - (*scplus5percent)(i,0)) )/ (*scdef)(i,0);
    double diff2 = fabs( ((*scdef)(i,0) - (*scminus5percent)(i,0)) )/ (*scdef)(i,0);
    printf("%2d     %6.4f       %6.4f       %6.4f         %6.2f            %6.2f\n", i, (*scdef)(i,0), 
	   (*scplus5percent)(i,0), (*scminus5percent)(i,0),  
	   diff1*100, diff2*100);
    diffSF[i] = diff1;
    if( diff2>diff1 )
      diffSF[i] = diff2;
  }

  // Cross sections

  TFile *fxsdef = new TFile(dirDefault
			    + TString("root_files/DY_m10+pr+a05+o03+pr_4839pb_fsrUnfGood/")
			    + TString("xSec_results_1D.root") );
  TMatrixD *xsdef = (TMatrixD*)fxsdef->Get("normXSecByBin");
  TVectorD *limits = (TVectorD*)fxsdef->Get("massBinLimits");
  
  TFile *fxsplus5percent = new TFile(dirPlus5
				     + TString("root_files/DY_m10+pr+a05+o03+pr_4839pb_fsrUnfGood/")
				     + TString("xSec_results_1D.root") );
  TMatrixD *xsplus5percent = (TMatrixD*)fxsplus5percent->Get("normXSecByBin");
  
  TFile *fxsminus5percent = new TFile(dirMinus5
				     + TString("root_files/DY_m10+pr+a05+o03+pr_4839pb_fsrUnfGood/")
				     + TString("xSec_results_1D.root") );
  TMatrixD *xsminus5percent = (TMatrixD*)fxsminus5percent->Get("normXSecByBin");
  
  double diffXS[40];
  printf("\nCross section\n");
  printf("bin      def        with+5%%      with-5%%   diff with+5%%, %% rel   diff with-5%%, %% rel\n");
  for(int i=0; i<40; i++){
    double diff1 = fabs( ((*xsdef)(i,0) - (*xsplus5percent)(i,0)) )/ (*xsdef)(i,0);
    double diff2 = fabs( ((*xsdef)(i,0) - (*xsminus5percent)(i,0)) )/ (*xsdef)(i,0);
    printf("%2d     %6.4f       %6.4f       %6.4f         %6.2f            %6.2f\n", i, (*xsdef)(i,0), 
	   (*xsplus5percent)(i,0), (*xsminus5percent)(i,0),  
	   diff1*100, diff2*100);
    diffXS[i] = diff1;
    if(diff2 > diff1 )
      diffXS[i] = diff2;
  }

  // Summary table
  printf("\nSummary table\n");
  printf(" mass    dEff,%%      dSF, %%       dXSEC, %%\n");
  for(int i=0; i<40; i++){
    printf("%4.0f-%4.0f  %7.2f     %7.2f      %7.2f\n", 
	   (*limits)(i), (*limits)(i+1),
	   diffEff[i]*100, diffSF[i]*100, diffXS[i]*100);

  }  

}
