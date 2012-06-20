#include "TVectorD.h"
#include "TFile.h"
#include "TString.h"
#include "TSystem.h"
#include <iostream>
#include <fstream>
#include <sstream>

#include "../Include/DYTools.hh"
#include "../Include/UnfoldingTools.hh"
#include "calcUnfoldingSystematics.C"
#include "../Include/ElectronEnergyScale.hh"

// returns 1 - ok, 0 - binning failure, -1 - file failure
int applyUnfoldingLocal(TVectorD &vin, TVectorD &vout, TString matrixFileName, int printLoadedData=1);

// save texTable
int printTexTable(const TString &texFileName, const std::vector<TString>& headers, const std::vector<TVectorD*> &data, const std::vector<double> &factors);


// -----------------------------------------------------------

void calcEscaleSystematics(TString lumiTag="DY_m10+pr+a05+o03+pr_4680pb", 
int saveTexTable=0){

  int nUnfoldingBins = DYTools::getTotalNumberOfBins();
  TVectorD observedYields(nUnfoldingBins);
  TVectorD observedYieldsErr(nUnfoldingBins);
  TVectorD dummyArr(nUnfoldingBins);
  TVectorD unfoldedYields(nUnfoldingBins);

  // 
  // Indicators (and also counters) whether the data was present
  //
  int countEScaleSyst=0;
  int countResidualShape=0;
  int countEtaSyst=0;
  int countFittingShape=0;
  std::vector<TString> usedFiles;

  ElectronEnergyScale escale("Date20120101_default");

  //
  // Calculate error associated with statistical 
  // uncertainty on the energy scale factors
  //
  TVectorD unfoldedYieldsMean(nUnfoldingBins);
  TVectorD unfoldedYieldsRMS(nUnfoldingBins);
  TVectorD unfoldedYieldsSquaredMean(nUnfoldingBins);

  unfoldedYieldsMean = 0;
  unfoldedYieldsRMS = 0;
  unfoldedYieldsSquaredMean = 0;

  TString matrixFileName = TString("../root_files/constants/") + lumiTag + 
    TString("/unfolding_constants") + analysisTag + TString(".root");
  const int nFiles1 = 20;  // expected number of files
  if (1)
  for(int ifile=0; ifile<nFiles1; ifile++){
    int seed = 1000+ifile;
    escale.randomizeEnergyScaleCorrections(seed);
    TString fname = TString("../root_files/yields/") + lumiTag + 
      TString("_escale_randomized/yields_bg-subtracted") + analysisTag + 
      TString("__") + escale.calibrationSetShortName();
    fname += ".root";
    TFile file(fname);
    if (!file.IsOpen()) {
      std::cout << "failed to open a file <" << fname << ">\n";
    }
    else {
      file.Close();

      // register
      usedFiles.push_back(fname);
      // work with data
      int res=	(readData(fname, observedYields,observedYieldsErr,dummyArr) == 1)
	&& (applyUnfoldingLocal(observedYields,unfoldedYields,matrixFileName) == 1);
      if (res==1) {
	countEScaleSyst++;
	//std::cout << "unfoldedYields for seed=" << seed << ": "; unfoldedYields.Print();
	// Accumulate mean and RMS
	for(int idx = 0; idx < nUnfoldingBins; idx++){
	  unfoldedYieldsMean[idx] += unfoldedYields[idx];
	  unfoldedYieldsSquaredMean[idx] += unfoldedYields[idx]*unfoldedYields[idx];
	}
      }
    }
  }

  // Final calculation of the mean and RMS for Smearing
  TVectorD escaleRandomizedSystRelative(nUnfoldingBins);
  if (countEScaleSyst) {
    for(int idx = 0; idx < nUnfoldingBins; idx++){
      unfoldedYieldsMean[idx] = unfoldedYieldsMean[idx]/double(countEScaleSyst);
      unfoldedYieldsSquaredMean[idx] = unfoldedYieldsSquaredMean[idx]/double(countEScaleSyst);
      unfoldedYieldsRMS[idx] = sqrt(unfoldedYieldsSquaredMean[idx] - 
				  unfoldedYieldsMean[idx]*unfoldedYieldsMean[idx]);
      escaleRandomizedSystRelative[idx] = unfoldedYieldsRMS[idx]/unfoldedYieldsMean[idx];
    }
    //std::cout << "unfoldedYieldsMean=" << ": "; unfoldedYieldsMean.Print();
    //std::cout << "escaleRandomizedSystRelative: "; escaleRandomizedSystRelative.Print();
  }
  else {
    std::cout << "escaleRandomizomized files were not found\n";
  }

  //
  // Calculate error related to the difference between
  // data and MC mass distribution shape after all corrections
  // had been applied to data and MC. "Residual differences"
  //
  TVectorD unfoldedYieldsVariation(nUnfoldingBins);
  TVectorD escaleShapeSystRelative(nUnfoldingBins);
  TString fname = TString("../root_files/yields/") + lumiTag + 
    TString("/yields_bg-subtracted") + analysisTag + TString(".root");
  int res=1;
  {
  TFile file(fname);
  if (file.IsOpen()) {
    file.Close();
    if (readData(fname,observedYields,observedYieldsErr,dummyArr) != 1) {
      std::cout << "failed to get data from file <" << fname << ">\n";
      throw 2;
    }
  //
    matrixFileName = TString("../root_files/constants/") + lumiTag + 
      TString("/unfolding_constants") + analysisTag + TString(".root");
    res=applyUnfoldingLocal(observedYields, unfoldedYields, matrixFileName);
    if (res==1) usedFiles.push_back(matrixFileName);

  //
    if (res==1) {
      matrixFileName = TString("../root_files/constants/") + lumiTag + 
	TString("_escale_residual/unfolding_constants") + analysisTag +
	TString("_escaleResidual.root");
      res=applyUnfoldingLocal(observedYields, unfoldedYieldsVariation, matrixFileName);
      if (res==1) {
	countResidualShape++;
	usedFiles.push_back(matrixFileName);
      }
    }
  }
  if (res==1) {
    for(int idx=0; idx<nUnfoldingBins; idx++){
      escaleShapeSystRelative[idx] = 
	fabs(unfoldedYields[idx] - unfoldedYieldsVariation[idx])
	/(unfoldedYields[idx] + unfoldedYieldsVariation[idx]);
    }
  }
  }

  //
  // Calculate error related to extra smearing function shape
  //
  std::vector<TString> shapeNames;
  shapeNames.push_back("_6binNegs_BreitWigner");
  shapeNames.push_back("_6binNegs_Voigtian");
  std::vector<TVectorD*> unfoldedYieldsShape;
  if (1)
  for (unsigned int i=0; i<shapeNames.size(); ++i) {
    TVectorD observedYieldsShape(nUnfoldingBins);
    TVectorD observedYieldsShapeErr(nUnfoldingBins);
    TString shapeFName=TString("../root_files/yields/") + lumiTag + 
      TString("_escale_shape/yields_bg-subtracted") + analysisTag +
      TString("_") + shapeNames[i] + TString(".root");
    TFile fShape(shapeFName);
    if (fShape.IsOpen()) {
      fShape.Close();
      usedFiles.push_back(shapeFName);
      if ( readData(shapeFName,observedYieldsShape,observedYieldsShapeErr,dummyArr) !=1) {
	std::cout << "failed to load data from <" << shapeFName << ">\n";
	throw 2;
      }
      TString matrixFileNameShape = TString("../root_files/constants/") + lumiTag + 
	TString("_escale_shape/unfolding_constants") + analysisTag +
	TString("_") + shapeNames[i] + TString(".root");
      TVectorD *shapeYields=new TVectorD(nUnfoldingBins);
      {
	TFile ftmp(TString("tmp_") + shapeNames[i] + TString(".root"),"recreate");
	observedYieldsShape.Write("YieldsSignal");
	TVectorD tmp(nUnfoldingBins);
	tmp=0;
	tmp.Write("YieldsSignalErr");
	tmp.Write("YieldsSignalSystErr");
	ftmp.Close();
      }
      res=applyUnfoldingLocal(observedYieldsShape, *shapeYields, matrixFileNameShape);
      if (res) {
	countFittingShape++;
	usedFiles.push_back(matrixFileNameShape);
	unfoldedYieldsShape.push_back(shapeYields);
      }
      else delete shapeYields;
    }
  }
  
  //
  TVectorD escaleFitShapeSystRelative(nUnfoldingBins);
  for(int idx=0; idx<nUnfoldingBins; idx++){
    double min = unfoldedYields[idx];
    double max = unfoldedYields[idx];
    for (unsigned int iShape=0; iShape<unfoldedYieldsShape.size(); ++iShape) {
      min = TMath::Min( min, (*unfoldedYieldsShape[iShape])[idx] );
      max = TMath::Max( max, (*unfoldedYieldsShape[iShape])[idx] );
    }
    escaleFitShapeSystRelative[idx] = fabs(max-min)/(max+min);
  }

  //
  // Calculate error related to eta binning
  //
  std::vector<TString> etaBinNames;
  etaBinNames.push_back("_2binNegs_Gauss");
  etaBinNames.push_back("_3EB3EENegs_Gauss");
  etaBinNames.push_back("_4binNegs_Gauss");
  etaBinNames.push_back("_4EB3EENegs_Gauss");
  etaBinNames.push_back("_5binNegs_Gauss");
  etaBinNames.push_back("_6binNegs_Gauss_20120119"); // default
  etaBinNames.push_back("_6bins_Gauss_20120119");
  std::vector<TVectorD*> unfoldedYieldsEta;
  unfoldedYieldsEta.reserve(etaBinNames.size());
  if (1)
  for (unsigned int i=0; i<etaBinNames.size(); ++i) {
    TString fEtaFileName=TString("../root_files/yields/") + lumiTag + 
      TString("_escale_eta/yields_bg-subtracted") + analysisTag +
      TString("_") + etaBinNames[i] + TString(".root");
    TFile fEta(fEtaFileName);
    if (fEta.IsOpen()) {
      fEta.Close();
      usedFiles.push_back(fEtaFileName);
      TVectorD observedYieldsEta(nUnfoldingBins);
      TVectorD observedYieldsEtaErr(nUnfoldingBins);
      if (readData(fEtaFileName,observedYieldsEta,observedYieldsEtaErr,dummyArr)!=1) {
	std::cout << "failed to get info from <" << fEtaFileName << ">\n";
	throw 2;
      }
      matrixFileName = TString("../root_files/constants/") + lumiTag + 
	TString("_escale_eta/unfolding_constants") + analysisTag +
	TString("_") + etaBinNames[i] + TString(".root");
      TVectorD *unfYields = new TVectorD(nUnfoldingBins);
      res=applyUnfoldingLocal(observedYieldsEta, *unfYields, matrixFileName);
      if (res==1) {
	countEtaSyst++;
	usedFiles.push_back(matrixFileName);
	unfoldedYieldsEta.push_back(unfYields);
      }
      else {
	delete unfYields;
      }
    }
  }

  //
  TVectorD escaleEtaBinSystRelative(nUnfoldingBins);
  for(int idx=0; idx<nUnfoldingBins; idx++){
    double min = unfoldedYields[idx];
    double max = unfoldedYields[idx];
    for (unsigned int iShape=0; iShape<unfoldedYieldsEta.size(); ++iShape) {
      min = TMath::Min( min, (*unfoldedYieldsEta[iShape])[idx] );
      max = TMath::Max( max, (*unfoldedYieldsEta[iShape])[idx] );
    }
    escaleEtaBinSystRelative[idx] = fabs(max-min)/(max+min);
  }

  // 
  // Put all errors together
  //
  TVectorD escaleSystPercent(nUnfoldingBins);
  for(int idx = 0; idx < nUnfoldingBins; idx++){
    escaleSystPercent[idx] = 100.0 *
      sqrt(
	   escaleRandomizedSystRelative[idx]*escaleRandomizedSystRelative[idx]
	   + escaleShapeSystRelative[idx]*escaleShapeSystRelative[idx]
	   + escaleFitShapeSystRelative[idx]*escaleFitShapeSystRelative[idx]
	   + escaleEtaBinSystRelative[idx]*escaleEtaBinSystRelative[idx]
	   );
  }

  printf("\n\n");
  printf("mass_range y_range stat_part shape-fit residual binning total,%% \n");
   for(int i=0, idx=0; i<nMassBins; i++){
    double *rapidityBinLimits=DYTools::getYBinLimits(i);
    for (int yi=0; yi<nYBins[i]; ++yi, ++idx) {
      printf("%4.0f-%4.0f  %4.2lf-%4.2lf  %6.2f  %6.2f  %6.2f  %6.2f  %6.2f\n",
	     massBinLimits[i],massBinLimits[i+1],
	     rapidityBinLimits[yi],rapidityBinLimits[yi+1],
	     escaleRandomizedSystRelative[idx]*100.0,
	     escaleFitShapeSystRelative[idx]*100.0,
	     escaleShapeSystRelative[idx]*100.0,
	     escaleEtaBinSystRelative[idx]*100.0,
	     escaleSystPercent[idx]
	     );
    }
    delete rapidityBinLimits;
  }

  std::cout << "\n\n";
  std::cout << std::string(80,'-') << "\n";
  std::cout << " Summary of successful loads:\n";
  std::cout << "  escale systematics   file count=  " << countEScaleSyst << "\n";
  std::cout << "  residual shape syst. file count=  " << countResidualShape << "\n";
  std::cout << "  eta systematics      file count=  " << countEtaSyst << "\n";
  std::cout << "  fitting shape        file count=  " << countFittingShape << "\n";
  std::cout << "\n";
  if (1) {
    std::cout << " Loaded files:\n";
    for (unsigned int i=0; i<usedFiles.size(); ++i) {
      std::cout << " " << usedFiles[i] << "\n";
    }
    std::cout << std::endl;
  }

  if (saveTexTable) {
    std::vector<TString> headers;
    std::vector<TVectorD*> data;
    std::vector<double> factors;
    headers.push_back("mass range");
    headers.push_back("rapidity range");
    headers.push_back("statistical, \\%"); data.push_back(&escaleRandomizedSystRelative); factors.push_back(100.);
    headers.push_back("shape, \\%"); data.push_back(&escaleFitShapeSystRelative); factors.push_back(100.);
    headers.push_back("residual, \\%"); data.push_back(&escaleShapeSystRelative); factors.push_back(100.);
    headers.push_back("$\\eta$ binning, \\%"); data.push_back(&escaleEtaBinSystRelative); factors.push_back(100.);
    headers.push_back("total, \\%"); data.push_back(&escaleSystPercent); factors.push_back(1.);
    
    TString texFName=TString("../root_files/systematics/") + lumiTag + 
      TString("/escale_systematics") + analysisTag + TString("_tmp.tex");
    printTexTable(texFName,headers,data,factors);
  }

  TString finalFName=TString("../root_files/systematics/") + lumiTag + TString("/escale_systematics_tmp.root");
  TFile fout(finalFName,"recreate");
  escaleSystPercent.Write("escaleSystPercent");
  fout.Close();

  return;
}


//-----------------------------------------------------------------
// Unfold
//-----------------------------------------------------------------
int  applyUnfoldingLocal(TVectorD &vin, TVectorD &vout, TString matrixFileName, int printLoadedData)
{

  // Read unfolding constants
  //std::cout << "unfold: Load constants from <" << matrixFileName << ">" << std::endl;

  // Construct file names

  int res=unfolding::unfold(vin, vout, matrixFileName);
  if (res!=1) {
    std::cout << " ... in function applyUnfoldingLocal[calcEscaleSystematics]\n";
    return 0;
  }

  // Print the result. Mainly for debugging purposes
  if (printLoadedData) {
    printf("\nUNFOLD: Results for the data from, yields:\n");
    printf("                   yields observed        after unfolding            \n");
    for(int i=0, idx=0; i<DYTools::nMassBins; i++){
      double *rapidityBinLimits=DYTools::getYBinLimits(i);
      for (int yi=0; yi<nYBins[i]; ++yi, ++idx) {
	printf("%4.0f-%4.0f %4.2lf-%4.2lf   %8.1f       %8.1f\n",
	       DYTools::massBinLimits[i], DYTools::massBinLimits[i+1],
	       rapidityBinLimits[yi],rapidityBinLimits[yi+1],
	       vin[idx], vout[idx]);
      }
      delete rapidityBinLimits;
    }
    printf("\n");
  }

  return 1;
}

//-----------------------------------------------------------------
// save texTable
//-----------------------------------------------------------------

int printTexTable(const TString &texFileName, const std::vector<TString>& headers, const std::vector<TVectorD*> &data, const std::vector<double> &factors) {
  if ((headers.size()!=data.size()+2) || (headers.size()<2)) {
    std::cout << "printTexTable. vector size mismatch: " << headers.size() << " headers and " << data.size() << " data\n";
    return 0;
  }

  std::string s;
  FILE *fout=fopen(texFileName.Data(),"w");
  fprintf(fout,"\n\n\n");
  fprintf(fout,"\\begin{table}[tbhH]\n");
  fprintf(fout,"\\caption{\\label{tbl:escaleSyst} Electron energy scale systematics}\n");
  fprintf(fout,"\\begin{center}\n\\begin{tabular}{|c|");
  for (unsigned int i=0; i<data.size(); ++i) fprintf(fout,"d|");
  fprintf(fout,"}\n");
  fprintf(fout,"\\hline\n");
  fprintf(fout," %s ",headers[0].Data());
  for (unsigned int i=1; i<headers.size(); ++i) fprintf(fout," & %s ",headers[i].Data());
  fprintf(fout,"\\\\\n\\hline\n");
  for(int i=0, idx=0; i<nMassBins; i++) {
    double *rapidityBinLimits=DYTools::getYBinLimits(i);
    for (int yi=0; yi<nYBins[i]; ++yi, ++idx) {
      fprintf(fout," %4.0lf $-$ %4.0lf ",massBinLimits[i],massBinLimits[i+1]);
      fprintf(fout," & %4.2lf $-$ %4.2lf ",rapidityBinLimits[yi],rapidityBinLimits[yi+1]);
      for (unsigned int j=0; j<headers.size()-2; ++j) {
	fprintf(fout," & %5.2lf ",(*data[j])[idx]*factors[j]);
      }
      fprintf(fout,"\\\\\n");
    }
  }
  fprintf(fout,"\\hline\n");
  fprintf(fout,"\\end{tabular}\n\\end{center}\n\\end{table}\n\n");
  fclose(fout);
  std::cout << "texFileName=<" << texFileName << "> saved\n";
  return 1;
}

