// Given ASCII input create root file with the cross section
//

#include "../Include/DYTools.hh"
#include "../Include/InputFileMgr.hh"
#include <TTree.h>
#include <TBranch.h>
#include <TVectorD.h>
//#include 
#include <sstream>
#include <fstream>

//------------------------------------------------------------------------------------------------------------------------

int debug=0;

typedef enum { _dmitrii=1, _matthew=2 } TFileNaming_t;
TFileNaming_t fileNaming=
  _matthew;
//_dmitrii;

//------------------------------------------------------------------------------------------------------------------------


void createRootFile_code(TString topDir, TString outputRootFile="", int iMassCheck=-1) {

  if (DYTools::study2D==0) {
    std::cout << "DYTools::study2D should be 1\n";
    return;
  }

  if (fileNaming==_matthew) {
    std::cout << "\n Assuming: pdfNegError is not on file\n";
  }


  //
  // Don't write TObject part of the objects
  //
  TDescriptiveInfo_t::Class()->IgnoreTObjectStreamer();

  if (!outputRootFile.Length()) {
    if (topDir.Contains("/")) {
      std::cout << "topDir cannot contain '/' if outputRootFile is empty\n";
      return;
    }
    outputRootFile=TString("xSec_") + topDir + TString(".root");
  }


  TDescriptiveInfo_t *description= new TDescriptiveInfo_t();
  description->_info.reserve(1000);

  std::string s;

  TMatrixD xSect(DYTools::nMassBins,DYTools::nYBinsMax);
  TMatrixD xSectErr(DYTools::nMassBins,DYTools::nYBinsMax);
  TMatrixD xSectSystErrPos(DYTools::nMassBins,DYTools::nYBinsMax);
  TMatrixD xSectSystErrNeg(DYTools::nMassBins,DYTools::nYBinsMax);
  xSect=0; xSectErr=0; 
  xSectSystErrPos=0; xSectSystErrNeg=0;

  for (int iMass=1; iMass<7; ++iMass) {
    if ((iMassCheck>=0) && (iMass!=iMassCheck)) continue;

    int iFileIdx=iMass-1;
    TString inputTextFileBase;
    switch(fileNaming) {
    case _dmitrii: 
      inputTextFileBase=
	Form("/dy%d/NNLO.output_dy%d_NNLO.txt",iFileIdx,iFileIdx);
      break;
    case _matthew:
      inputTextFileBase=
	Form("/namgrid_e_NNLO_%1.0fto%1.0f_abkm/namgrid_e_NNLO_%1.0fto%1.0f_abkm.output",
	     DYTools::massBinLimits[iMass], DYTools::massBinLimits[iMass+1],
	     DYTools::massBinLimits[iMass], DYTools::massBinLimits[iMass+1]);
      break;
    default:
      std::cout << " code is not ready to work with this fileNaming\n";
      assert(0);
    }
	
    TString inputTextFile = topDir + inputTextFileBase;
    std::cout << "reading text file <" << inputTextFile << ">" << std::endl;

    std::ifstream fin(inputTextFile);
    if (!fin.is_open()) {
      std::cout << "Failed to open the file <" << inputTextFile << ">\n";
      return;
    }
    
    int stage=0;
    int downCount=1;

    description->_info.push_back(Form("File: %s",inputTextFile.Data()));
    double massMin=-1, massMax=-1;
    double y,cs,csErr,csSErr1,csSErr2;
    double yPrev=0,csPrev=0;

    while (!fin.eof() && getline(fin,s)) {
      description->_info.push_back(s);
      if (debug) std::cout << "stage=" << stage << ", s=<" << s << ">\n";
      //const char *c=s.c_str();
      if ((stage==0) && PosOk(s,"Lepton-pair invariant mass")) {
	size_t tmpPos=s.find('=');
	if (!PosOk(tmpPos)) {
	  std::cout << "\nError:\nline=<" << s << ">\n";
	  std::cout << "a line with 'Lepton-pair invariant mass' does not contain '='\n";
	  return;
	}
	//std::cout << "chk: <" << s.substr(tmpPos+1,s.size()) << ">\n";
	std::stringstream ss(s.substr(tmpPos+1,s.size()));
	if (PosOk(s,"minimum")) {
	  ss >> massMin;
	}
	else {
	  ss >> massMax;
	  std::cout << "massMin=" << massMin << ", massMax=" << massMax << "\n";
	  stage++;
	  
	  int idxMin=DYTools::findMassBin(massMin);
	  int idxMax=DYTools::findMassBin(massMax);
	  if ((idxMax==-1) && (massMax=1500)) idxMax=DYTools::nMassBins;
	  if ((idxMin<0) || (idxMax<0) || (idxMax!=idxMin+1) 
	      || (iMass!=idxMin)) {
	    std::cout << "\nError:\nline=<" << s << ">\n";
	    printf("  expecting iMass=%d\n",iMass);
	    printf("  massMin=%4.1lf, massMinIdx=%d,  massMax=%4.1lf, massMaxIdx=%d\n",massMin,idxMin,massMax,idxMax);
	    return;
	  }
	}
      }
      else if ((stage==1) && PosOk(s,"Z/W rapidity")) {
	stage++;
      }
      else if (stage==2) {
	if (downCount>0) { downCount--; continue; }
	if (s.size()==0) { stage++; continue; }
	if (iMass<0) {
	  std::cout << "table stage was reached, but iMass<0\n";
	  return;
	}
	if (debug) std::cout << "read <" << s << ">\n";
	std::stringstream ss(s);
	ss >> y;
	if (y<0) continue;
	else {
	  ss >> cs >> csErr >> csSErr1 >> csSErr2;
	  if (fileNaming==_matthew) { csSErr2=csSErr1; }
	  int iY=DYTools::findYBin(iMass,y);
	  if (iY==-1) {
	    std::cout << "iY=" << iY << " from y=" << y << " in massBin=" << iMass << "\n";
	    return;
	  }
	  if (debug) printf(" keeping for (%d,%d): %6.4lf %8.6lf %8.6lf %8.6lf\n",iMass,iY,cs,csErr,csSErr1,csSErr2);
	  //std::cout << "iMass=" << iMass << ", nMassBins=" << DYTools::nMassBins << "\n";
	  if (iMass!=DYTools::nMassBins-1) {
	    xSect[iMass][iY] += cs;
	    xSectErr[iMass][iY] += csErr*csErr;
	    xSectSystErrPos[iMass][iY] += csSErr1*csSErr1;
	    xSectSystErrNeg[iMass][iY] += csSErr2*csSErr2;
	  }
	  else {
	    int iYMin=DYTools::findYBin(iMass,y-0.0495);
	    int iYMax=DYTools::findYBin(iMass,y+0.0495);
	    if ((iYMin==iY) && (iYMax==iY)) {
	      xSect[iMass][iY] += cs;
	      xSectErr[iMass][iY] += csErr*csErr;
	      xSectSystErrPos[iMass][iY] += csSErr1*csSErr1;
	      xSectSystErrNeg[iMass][iY] += csSErr2*csSErr2;
	      yPrev=y; csPrev=cs;
	    }
	    else {
	      std::cout << "cs correction for y=" << y << ", iYMin=" << iYMin << ", iYMax=" << iYMax << ", iY=" << iY << "\n";
	      double y_1=y-0.15;
	      double y0=y-0.05;
	      double y1=y+0.05;
	      double wsigma_1to0=csPrev;
	      double wsigma0to1=cs;
	      double wsigma0to1Err=csErr;
	      double ystar=y;
	      double wSigmaStar1=0, wSigmaStar1Err=0;
	      double wSigmaStar2=0, wSigmaStar2Err=0;
	      subdivideBinWeightByLinearApprox(y_1, wsigma_1to0,
					       y0, wsigma0to1, wsigma0to1Err,
					       y1,
					       0., -1e9,
					       ystar,
					       wSigmaStar1,wSigmaStar1Err,
					       wSigmaStar2,wSigmaStar2Err);
	      xSect[iMass][iYMin] += wSigmaStar1;
	      xSectErr[iMass][iYMin] += wSigmaStar1Err*wSigmaStar1Err;
	      xSect[iMass][iY] += wSigmaStar2;
	      xSectErr[iMass][iY] += wSigmaStar2Err*wSigmaStar2Err;

	      // divide 1st syst error
	      subdivideBinWeightByLinearApprox(y_1, wsigma_1to0,
					       y0, wsigma0to1, csSErr1,
					       y1,
					       0., -1e9,
					       ystar,
					       wSigmaStar1,wSigmaStar1Err,
					       wSigmaStar2,wSigmaStar2Err);

	      xSectSystErrPos[iMass][iYMin] += wSigmaStar1Err*wSigmaStar1Err;
	      xSectSystErrPos[iMass][iY] += wSigmaStar2Err*wSigmaStar2Err;

	      // divide 2nd syst error
	      subdivideBinWeightByLinearApprox(y_1, wsigma_1to0,
					       y0, wsigma0to1, csSErr2,
					       y1,
					       0., -1e9,
					       ystar,
					       wSigmaStar1,wSigmaStar1Err,
					       wSigmaStar2,wSigmaStar2Err);

	      xSectSystErrNeg[iMass][iYMin] += wSigmaStar1Err*wSigmaStar1Err;
	      xSectSystErrNeg[iMass][iY] += wSigmaStar2Err*wSigmaStar2Err;
	    }
	  }
	}
      }
      else if (stage==3) {
	break;
      }
    }
    if (iMass!=5) description->_info.push_back("\n\n");
    fin.close();
  }

  int nUnfoldingBins = DYTools::getTotalNumberOfBins();
  TVectorD xSectFI(nUnfoldingBins);
  TVectorD xSectErrFI(nUnfoldingBins);
  xSectFI=0; xSectErrFI=0;

  int idx=0;
  double zPeakCS=0, zPeakCSErr=0;
  for (int iM=0; iM<DYTools::nMassBins; ++iM) {
    double mass=DYTools::massBinLimits[iM]+0.1;
    int isZpeak= ((mass>=60) && (mass<=120)) ? 1:0;
    std::cout << "mass=" << mass << ", isZpeak=" << isZpeak << "\n";
    for (int iY=0;  iY<DYTools::nYBins[iM]; ++iY, ++idx) {
      if (isZpeak) {
	zPeakCS+= xSect[iM][iY];
	zPeakCSErr+=xSectErr[iM][iY];
      }
      xSectFI[idx]=xSect[iM][iY];
      xSectErr[iM][iY] = sqrt(xSectErr[iM][iY]);
      xSectErrFI[idx]=xSectErr[iM][iY];
      xSectSystErrPos[iM][iY] = sqrt(xSectSystErrPos[iM][iY]);
      xSectSystErrNeg[iM][iY] = sqrt(xSectSystErrNeg[iM][iY]);
    }
  }
  zPeakCSErr= sqrt(zPeakCSErr);
  std::cout << "zPeakCS=" << zPeakCS << "\n";

  TMatrixD xSectNorm=xSect;
  TMatrixD xSectNormErr=xSectErr;
  TMatrixD xSectNormSystErrPos=xSectSystErrPos;
  TMatrixD xSectNormSystErrNeg=xSectSystErrNeg;
  xSectNorm *= (1/zPeakCS);
  xSectNormErr *= (1/zPeakCS);
  xSectNormSystErrPos *= (1/zPeakCS);
  xSectNormSystErrNeg *= (1/zPeakCS);

  // 
  //  Make plots
  // 
  TFile fout(outputRootFile,"recreate");
  xSect.Write("xSectThDET");
  xSectErr.Write("xSectThDETErr");
  xSectSystErrPos.Write("xSectThDETSystErrPos");
  xSectSystErrNeg.Write("xSectThDETSystErrNeg");
  xSectFI.Write("xSectThDETFIArray");
  xSectErrFI.Write("xSectThDETErrFIArray");

  xSectNorm.Write("xSectThDETNorm");
  xSectNormErr.Write("xSectThDETNormErr");
  xSectNormSystErrPos.Write("xSectThDETNormSystErrPos");
  xSectNormSystErrNeg.Write("xSectThDETNormSystErrNeg");

  TTree *descriptionTree = new TTree("Description","description");
  descriptionTree->Branch("description","TDescriptiveInfo_t",&description);
  descriptionTree->Fill();
  descriptionTree->Write();
  fout.Close();

  // 
  //  Make plots
  // 
  TString plotsFName=outputRootFile;
  plotsFName.ReplaceAll(".root","-plots.root");
  TFile fPlots(plotsFName,"recreate");
  for (int iM=1; iM<6; ++iM) {
    int perMassBinWidth=0;
    TString hName=Form("xsec_Mass_%1.0lf_%1.0lf",DYTools::massBinLimits[iM],DYTools::massBinLimits[iM+1]);
    TH1F *h=extractRapidityDependence(hName,"",
				      xSect, xSectErr,
				      iM,perMassBinWidth);
    assert(h);
    h->Write();
    TString hNameNorm=Form("xsecNorm_Mass_%1.0lf_%1.0lf",DYTools::massBinLimits[iM],DYTools::massBinLimits[iM+1]);
    TH1F *hNorm=extractRapidityDependence(hNameNorm,"",
					  xSectNorm, xSectNormErr,
					  iM,perMassBinWidth);
    hNorm->Write();
    assert(hNorm);
  }

  TMatrixD zero=xSectErr;
  zero=0;
  for (int iM=1; iM<7; ++iM) {
    int perMassBinWidth=0;
    TString hName=Form("xsecErr_Mass_%1.0lf_%1.0lf",DYTools::massBinLimits[iM],DYTools::massBinLimits[iM+1]);
    TH1F *h=extractRapidityDependence(hName,"",
				      xSectErr, zero,
				      iM,perMassBinWidth);
    assert(h);
    h->Write();
    TString hNameNorm=Form("xsecNormErr_Mass_%1.0lf_%1.0lf",DYTools::massBinLimits[iM],DYTools::massBinLimits[iM+1]);
    TH1F *hNorm=extractRapidityDependence(hNameNorm,"",
					  xSectNormErr, zero,
					  iM,perMassBinWidth);
    hNorm->Write();
    assert(hNorm);
  }
  fPlots.Close();

  std::cout << "files <" << outputRootFile << "> and <" << plotsFName << "> created\n";
}


// 
// Sample ASCII file NNLO.output_dy0_NNLO.txt:
/*

 CMS collision Energy =    7000.00000000000     
 Collider type = pp
 ======================================
 Factorization scale                =    25.0000000000000     
 Renormalization scale              =    25.0000000000000     
 ==================================
 NNLO Cross-Section
 Strong coupling                    =   0.146798992895555     
 PDF set = CTEQ10W                                                     
 ======================================
 Z pole focus OFF
 Lower M limit for non-Zpole     =    20.0000000000000     
 Upper M limit for non-Zpole     =    30.0000000000000     
 ======================================
 Alpha QED                          =   7.812500000000000E-003
 ======================================
 Z mass (GeV)                       =    91.1876000000000     
 Z width (GeV)                      =    2.49520000000000     
 Z->ll partial width                =   8.398799999999999E-002
 sin^2(theta)                       =   0.222550000000000     
 up quark charge                    =   0.666666700000000     
 down quark charge                  =  -0.333333300000000     
 lepton chage                       =   -1.00000000000000     
 up quark vector coupling           =   0.409100000000000     
 down quark vector coupling         =  -0.704500000000000     
 lepton vector coupling             =  -0.113600000000000     
 up quark axial coupling            =   -1.00000000000000     
 down quark axial coupling          =    1.00000000000000     
 lepton axial coupling              =    1.00000000000000     
 ======================================
 Lepton-pair invariant mass minimum =    20.0000000000000     
 Lepton-pair invariant mass maximum =    30.0000000000000     
 Transverse mass minimum            =    0.00000000000000     
 Transverse mass maximum            =    1000.00000000000     
 Z pT minimum                       =    0.00000000000000     
 Z pT maximum                       =    1000.00000000000     
 Z rapidity minimum                 =   -20.0000000000000     
 Z rapidity maximum                 =    20.0000000000000     
 Lepton pT minimum                  =    0.00000000000000     
 Lepton pT maximum                  =    1000.00000000000     
 Anti-lepton pT minimum             =    0.00000000000000     
 Anti-lepton pT maximum             =    1000.00000000000     
 Lepton eta absolute value?         =   T
 Lepton pseudorapidity minimum      =    0.00000000000000     
 Lepton pseudorapidity maximum      =    100.000000000000     
 Anti-lepton eta absolute value?    =   T
 Anti-lepton pseudorapidity minimum =    0.00000000000000     
 Anti-lepton pseudorapidity maximum =    100.000000000000     
 ======================================
 ======================================
 Jet merging algorithm              =   ktal
 Jet algorithm cone size (deltaR)   =   0.400000000000000     
 parton-parton Rsep (cone algo only)=    1.30000000000000     
 Minimum pT for Observable Jets     =    20.0000000000000     
 Maximum eta for Observable Jets    =    4.50000000000000     
 Minimum Number of Jets             =             0
 Maximum Number of Jets             =             2
 Leading jet pT minimum             =    0.00000000000000     
 ======================================
 Lep-Anti-lep deltaR minimum        =    0.00000000000000     
 Lep-Anti-lep deltaPhi minimum      =    0.00000000000000     
 Lep-Anti-lep deltaPhi maximum      =    4.00000000000000     
 Lep-Jet deltaR minimum             =    0.00000000000000     
 ======================================
 Rapidity of Z cutoff for CS angles =    0.00000000000000     
 
 ===========   VEGAS PARAMETERS   ========
 
 Maximum number of evaluations      =   1000000000
 Actual number of evaluations       =   1001250000
 Requested relative precision (%)   =   1.000000000000000E-003
 Requested absolute precision       =   1.000000000000000E-019
 Total integration time (sec)       =    78400.0000000000     
 Random number seed                 =          115
 
 ===========   RESULT    =================
 
 Sigma (pb)                  =    5.78539
 Integration Error (pb)      =    0.00131464
 PDF Error (pb)              +   0.197178
                             -   0.23181
 
 ===========   HISTOGRAMS    =================
      bin           weight   numerical error       + pdf error      - pdf error

  ----   Z/W pT         ----

     5.00        0.0188143        0.00053315       0.000879686        0.00108957
    15.00         0.519899       0.000730248         0.0234495         0.0284606
    25.00           1.3106        0.00072443         0.0529867         0.0628402
    35.00          1.58019       0.000735487         0.0584937          0.067418
    45.00          0.94733       0.000552241         0.0315707         0.0361062
    55.00         0.524858       0.000435369         0.0156818         0.0179606
    65.00         0.305929       0.000364372        0.00825348        0.00950309
    75.00         0.187546       0.000316267        0.00458128        0.00530729
    85.00         0.119931       0.000280736         0.0026806        0.00313583
    95.00         0.079163       0.000252857        0.00165098        0.00195654
   105.00        0.0539829        0.00023132        0.00104698        0.00125136
   115.00        0.0372821       0.000213273       0.000677827       0.000814824
   125.00        0.0267786       0.000197921       0.000485841       0.000580567
   135.00        0.0190717       0.000183625       0.000336358       0.000402097
   145.00        0.0144675       0.000171316       0.000251737       0.000303894
   155.00        0.0105534       0.000159032       0.000168265       0.000212699
   165.00       0.00790566       0.000148803       0.000147774        0.00017311
   175.00        0.0062648       0.000140362       0.000118872       0.000140128
   185.00       0.00502181       0.000132318       0.000102428       0.000119815
   195.00       0.00396197       0.000124652       8.20816e-05        9.5067e-05

  ----   Z/W rapidity   ----

    -2.35        0.0455816       0.000131884        0.00128592        0.00184215
    -2.25        0.0718099       0.000158095        0.00199934        0.00288739
    -2.15        0.0900085       0.000171997        0.00247109        0.00358097
    -2.05         0.102835       0.000179118        0.00280258        0.00405157
    -1.95         0.110416       0.000189566        0.00303481        0.00435052
    -1.85         0.110587       0.000189425        0.00306046        0.00433126
    -1.75         0.107094       0.000183323          0.003016        0.00417267
    -1.65         0.102556        0.00018038        0.00293366        0.00393218
    -1.55         0.106177       0.000176734        0.00316624        0.00415146
    -1.45         0.106019       0.000176981        0.00325717        0.00416627
    -1.35         0.107274       0.000179214        0.00332173        0.00416848
    -1.25         0.115059       0.000179943          0.003787        0.00458873
    -1.15         0.122407       0.000183696        0.00419322        0.00492338
    -1.05         0.129459       0.000186666        0.00461158         0.0052593
    -0.95         0.135694       0.000188844        0.00506227        0.00561527
    -0.85         0.140129       0.000191094        0.00542354        0.00584177
    -0.75         0.143039       0.000191867        0.00572018        0.00601971
    -0.65         0.145696       0.000192232        0.00600145        0.00619412
    -0.55         0.146928       0.000193048        0.00626909        0.00633833
    -0.45         0.148578       0.000195029        0.00643387        0.00642733
    -0.35         0.148971       0.000193669        0.00658809        0.00649957
    -0.25         0.148943       0.000195355         0.0066843        0.00653344
    -0.15         0.149201       0.000194116        0.00674969         0.0065564
    -0.05         0.149286       0.000173725         0.0068248        0.00660367
     0.05         0.149269       0.000173474        0.00682741        0.00660213
     0.15         0.149196       0.000193866        0.00674678        0.00655645
     0.25         0.148926       0.000195303        0.00668778        0.00653833
     0.35         0.148762       0.000193184        0.00658314        0.00649369
     0.45         0.148393       0.000193752        0.00644139         0.0064273
     0.55         0.147004       0.000191617        0.00626219        0.00633716
     0.65         0.145391       0.000190993        0.00598994        0.00618496
     0.75         0.142986       0.000190758        0.00571786        0.00602201
     0.85         0.139948       0.000189879        0.00542193        0.00584526
     0.95         0.135414       0.000187644        0.00504668        0.00559313
     1.05         0.129243       0.000185384         0.0046124        0.00526133
     1.15         0.122037       0.000181902        0.00417499        0.00490981
     1.25         0.114739       0.000178478        0.00378109        0.00458088
     1.35         0.107195       0.000177854        0.00331831        0.00417347
     1.45         0.105926       0.000175955        0.00325375        0.00416473
     1.55         0.106113       0.000175595         0.0031559        0.00414975
     1.65         0.102243       0.000178483        0.00293518        0.00393551
     1.75         0.106827       0.000181719        0.00300699        0.00416725
     1.85         0.110309       0.000187889        0.00305797        0.00432766
     1.95         0.110059       0.000188463        0.00302222        0.00433481
     2.05         0.102458       0.000177954        0.00279928        0.00404223
     2.15        0.0898419       0.000170922        0.00246543        0.00357472
     2.25        0.0715781       0.000157306        0.00199424        0.00288692
     2.35        0.0453411       0.000130932        0.00128031        0.00183702

  ----   Q_ll Invaria   ----

    20.50         0.582849       0.000472115         0.0187087         0.0218656
    21.50          0.56351       0.000440501         0.0182408         0.0213366
    22.50          0.54914       0.000414894         0.0178796         0.0209174
    23.50          0.53922       0.000393764         0.0178001         0.0208087
    24.50          0.53408       0.000377155         0.0179443         0.0210114
    25.50         0.536555       0.000365601         0.0181046         0.0211964
    26.50         0.547174       0.000361213         0.0187789         0.0220249
    27.50         0.570944       0.000376661         0.0199552         0.0234696
    28.50         0.630773       0.000449368         0.0225685         0.0266626
    29.50         0.730416       0.000498232         0.0274055         0.0325643


*/
