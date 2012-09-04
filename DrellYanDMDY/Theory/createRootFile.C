// Given ASCII input create root file with the cross section
//

#include "../Include/DYTools.hh"
#include "../Include/InputFileMgr.hh"
#include "../Include/MyTools.hh"
#include <TTree.h>
#include <TBranch.h>
#include <TVectorD.h>
//#include 
#include <sstream>
#include <fstream>

void createRootFile(TString inputTextFile="", TString outputRootFile="") {

  //
  // Don't write TObject part of the objects
  //
  TDescriptiveInfo_t::Class()->IgnoreTObjectStreamer();


  TDescriptiveInfo_t *description= new TDescriptiveInfo_t();
  description->_info.reserve(1000);

  std::ifstream fin(inputTextFile);
  if (!fin.is_open()) {
    std::cout << "Failed to open the file <" << inputTextFile << ">\n";
    return;
  }

  std::string s;
  int downCount=-1;

  const int downCountForLabels=1; // number of lines to skip until labels
  const int downCountForTable=1;  // number of lines to skip until table
  const char *chkHeaderString="numerical error";


  TMatrixD xSect(DYTools::nMassBins,DYTools::nYBinsMax);
  TMatrixD xSectErr(DYTools::nMassBins,DYTools::nYBinsMax);
  TMatrixD xSectSystErrPos(DYTools::nMassBins,DYTools::nYBinsMax);
  TMatrixD xSectSystErrNeg(DYTools::nMassBins,DYTools::nYBinsMax);
  xSect=0; xSectErr=0; 
  xSectSystErrPos=0; xSectSystErrNeg=0;

  int iMass=-1;
  int stage=0;

  while (!fin.eof() && getline(fin,s)) {
    description->_info.push_back(s);
    //const char *c=s.c_str();
    if ((stage==0) && PosOk(s,"GeV")) {
      std::stringstream ss(s);
      double massMin=-1, massMax=-1;
      ss >> massMin >> massMax;
      if (massMax<0) massMax=-massMax;
      std::cout << "massMin=" << massMin << ", massMax=" << massMax << "\n";
      int idxMin=DYTools::findMassBin(massMin);
      int idxMax=DYTools::findMassBin(massMax);
      if ((idxMax==-1) && (massMax=1500)) idxMax=DYTools::nMassBins;
      if ((idxMin<0) || (idxMax<0) || (idxMax!=idxMin+1)) {
	std::cout << "\nError:\nline=<" << s << ">\n";
	printf("  massMin=%4.1lf, massMinIdx=%d,  massMax=%4.1lf, massMaxIdx=%d\n",massMin,idxMin,massMax,idxMax);
	return;
      }
      iMass=idxMin;
      downCount=downCountForLabels;
      stage++;
    }
    else if (stage==1) {
      if (downCount>0) { downCount--; continue; }
      if (PosOk(s,chkHeaderString)) {
	stage++;
	downCount=downCountForTable;
      }
      else {
	std::cout << "check failed: could not locate <" 
		  << chkHeaderString << "> in the expected line\n";
	return;
      }
    }
    else if (stage==2) {
      if (downCount>0) { downCount--; continue; }
      if (s.size()==0) { stage=0; continue; }
      if (iMass<0) {
	std::cout << "table stage was reached, but iMass<0\n";
	return;
      }
      std::stringstream ss(s);
      double y,cs,csErr,csSErr1,csSErr2;
      ss >> y;
      if (y<0) continue;
      else {
	ss >> cs >> csErr >> csSErr1 >> csSErr2;
	int iY=DYTools::findYBin(iMass,y);
	if (iY==-1) {
	  std::cout << "iY=" << iY << " from y=" << y << " in massBin=" << iMass << "\n";
	  return;
	}
	printf(" keeping for (%d,%d): %6.4lf %8.6lf %8.6lf %8.6lf\n",iMass,iY,cs,csErr,csSErr1,csSErr2);
	std::cout << "iMass=" << iMass << ", nMassBins=" << DYTools::nMassBins << "\n";
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
	  }
	  else {
	    std::cout << "cs correction for y=" << y << ", iYMin=" << iYMin << ", iYMax=" << iYMax << ", iY=" << iY << "\n";
	    xSect[iMass][iYMin] += 0.5*cs;
	    xSectErr[iMass][iYMin] += 0.5*csErr*csErr;
	    xSectSystErrPos[iMass][iYMin] += 0.5*csSErr1*csSErr1;
	    xSectSystErrNeg[iMass][iYMin] += 0.5*csSErr2*csSErr2;
	    xSect[iMass][iY] += 0.5*cs;
	    xSectErr[iMass][iY] += 0.5*csErr*csErr;
	    xSectSystErrPos[iMass][iY] += 0.5*csSErr1*csSErr1;
	    xSectSystErrNeg[iMass][iY] += 0.5*csSErr2*csSErr2;
	  }
	}
      }
    }
  }
  fin.close();

  int nUnfoldingBins = DYTools::getTotalNumberOfBins();
  TVectorD xSectFI(nUnfoldingBins);
  TVectorD xSectErrFI(nUnfoldingBins);
  xSectFI=0; xSectErrFI=0;

  int idx=0;
  double zPeakCS=0, zPeakCSErr=0;
  for (int iM=0; iM<DYTools::nMassBins; ++iM) {
    double mass=DYTools::massBinLimits[iM]+0.1;
    int isZpeak= ((mass>=60) && (mass<=120)) ? 1:0;
    std::cout << "mass=" << mass << "\n";
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

  TFile fout(outputRootFile,"recreate");
  xSect.Write("xSectThDET");
  xSectErr.Write("xSectThDETErr");
  xSectSystErrPos.Write("xSectThDETSystErrPos");
  xSectSystErrNeg.Write("xSectThDETSystErrNeg");
  xSectFI.Write("xSectThDETFIArray");
  xSectErrFI.Write("xSectThDETErrFIArray");

  xSect *= (1/zPeakCS);
  xSectErr=0;
  xSect.Write("xSectThDETNorm");
  xSectErr.Write("xSectThDETNormErr");

  TTree *descriptionTree = new TTree("Description","description");
  descriptionTree->Branch("description","TDescriptiveInfo_t",&description);
  descriptionTree->Fill();
  descriptionTree->Write();
  fout.Close();
  std::cout << "file <" << outputRootFile << "> created\n";
}


// 
// Sample ASCII file start:
/*

Dimitri Bourilkov  UF  24-Aug-2012

CTEQ10W - di-electrons in CMS acceptance


20-30 GeV

 bin           weight    numerical error       + pdf error       - pdf error

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


30-45 GeV

 bin           weight    numerical error       + pdf error       - pdf error

-2.35          0.14559        0.00167554        0.00561575        0.00754393
-2.25         0.236062        0.00186169        0.00912496         0.0125639
-2.15         0.299846         0.0019937         0.0114233         0.0156402
-2.05          0.34686        0.00207365         0.0131558         0.0180529
...


*/
