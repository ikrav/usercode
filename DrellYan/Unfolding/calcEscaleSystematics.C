#include "TVectorD.h"
#include "TFile.h"
#include "TString.h"
#include "TSystem.h"
#include <iostream>
#include <fstream>
#include <sstream>

#include "../Include/DYTools.hh"
#include "../Include/UnfoldingTools.hh"

// returns 1 - ok, 0 - binning failure, -1 - file failure
int applyUnfoldingShort(TVectorD &vin, TVectorD &vout, TString matrixFileName, int printLoadedData=0);

// save texTable
int printTexTable(const TString &texFileName, int nMassBins, const double *massBins, const std::vector<TString>& headers, const std::vector<TVectorD*> &data, const std::vector<double> &factors);


// -----------------------------------------------------------

void calcEscaleSystematics(TString lumiTag="DY_m10+pr+a05+o03+pr_4680pb", 
int saveTexTable=1){

  TVectorD observedYields(nMassBins);
  TVectorD unfoldedYields(nMassBins);

  // 
  // Indicators (and also counters) whether the data was present
  //
  int flagEScaleSyst=0;
  int flagResidualShape=0;
  int flagEtaSyst=0;
  int flagFittingShape=0;
  std::vector<TString> usedFiles;

  //
  // Calculate error associated with statistical 
  // uncertainty on the energy scale factors
  //
  TVectorD unfoldedYieldsMean(nMassBins);
  TVectorD unfoldedYieldsRMS(nMassBins);
  TVectorD unfoldedYieldsSquaredMean(nMassBins);

  unfoldedYieldsMean = 0;
  unfoldedYieldsRMS = 0;
  unfoldedYieldsSquaredMean = 0;

  TString matrixFileName = TString("../root_files/constants/") + lumiTag + TString("/unfolding_constants.root");
  const int nFiles1 = 20;  // expected number of files
  for(int ifile=0; ifile<nFiles1; ifile++){
    int seed = 1001+ifile;
    TString fname = TString("../root_files/yields/") + lumiTag + TString("_escale_randomized/yields_bg-subtracted_seed");
    fname += seed;
    fname += ".root";
    TFile file(fname);
    if (file.IsOpen()) {
      // register
      usedFiles.push_back(fname);
      // work with data
      observedYields = *(TVectorD*)file.Get("YieldsSignal");
      int res=applyUnfoldingShort(observedYields, unfoldedYields,matrixFileName);
      if (res==1) {
	flagEScaleSyst++;
	//std::cout << "unfoldedYields for seed=" << seed << ": "; unfoldedYields.Print();
	// Accumulate mean and RMS
	for(int iMassBin = 0; iMassBin < nMassBins; iMassBin++){
	  unfoldedYieldsMean[iMassBin] += unfoldedYields[iMassBin];
	  unfoldedYieldsSquaredMean[iMassBin] += unfoldedYields[iMassBin]*unfoldedYields[iMassBin];
	}
      }
      file.Close();
    }
  }

  // Final calculation of the mean and RMS for Smearing
  TVectorD escaleRandomizedSystRelative(nMassBins);
  if (flagEScaleSyst) {
    for(int i = 0; i < nMassBins; i++){
      unfoldedYieldsMean[i] = unfoldedYieldsMean[i]/double(flagEScaleSyst);
      unfoldedYieldsSquaredMean[i] = unfoldedYieldsSquaredMean[i]/double(flagEScaleSyst);
      unfoldedYieldsRMS[i] = sqrt(unfoldedYieldsSquaredMean[i] - 
				  unfoldedYieldsMean[i]*unfoldedYieldsMean[i]);
      escaleRandomizedSystRelative[i] = unfoldedYieldsRMS[i]/unfoldedYieldsMean[i];
    }
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
  TVectorD unfoldedYieldsVariation(nMassBins);
  TVectorD escaleShapeSystRelative(nMassBins);
  TString fname = TString("../root_files/yields/") + lumiTag + TString("/yields_bg-subtracted.root");
  TFile file(fname);
  int res=1;
  if (file.IsOpen()) {
    observedYields = *(TVectorD*)file.Get("YieldsSignal");
  //
    matrixFileName = TString("../root_files/constants/") + lumiTag + TString("/unfolding_constants.root");
    res=applyUnfoldingShort(observedYields, unfoldedYields, matrixFileName);
    if (res==1) usedFiles.push_back(matrixFileName);

  //
    if (res==1) {
      matrixFileName = TString("../root_files/constants/") + lumiTag + TString("_escale_residual/unfolding_constants.root");
      res=applyUnfoldingShort(observedYields, unfoldedYieldsVariation, matrixFileName);
      if (res==1) {
	flagResidualShape++;
	usedFiles.push_back(matrixFileName);
      }
    }
  }
  if (res==1) {
    for(int i=0; i<nMassBins; i++){
      escaleShapeSystRelative[i] = fabs(unfoldedYields[i] - unfoldedYieldsVariation[i])
	/(unfoldedYields[i] + unfoldedYieldsVariation[i]);
    }
  }

  //
  // Calculate error related to extra smearing function shape
  //
  std::vector<TString> shapeNames;
  shapeNames.push_back("_6binNegs_BreitWigner");
  shapeNames.push_back("_6binNegs_Voigtian");
  std::vector<TVectorD*> unfoldedYieldsShape;
  for (unsigned int i=0; i<shapeNames.size(); ++i) {
    TVectorD observedYieldsShape(nMassBins);
    TString shapeFName=TString("../root_files/yields/") + lumiTag + TString("_escale_shape/yields_bg-subtracted") + shapeNames[i] + TString(".root");
    TFile fShape(shapeFName);
    if (fShape.IsOpen()) {
      usedFiles.push_back(shapeFName);
      observedYieldsShape = *(TVectorD*)file.Get("YieldsSignal");
      TString matrixFileNameShape = TString("../root_files/constants/") + lumiTag + TString("_escale_shape/unfolding_constants") + shapeNames[i] + TString(".root");
      TVectorD *shapeYields=new TVectorD(nMassBins);
      res=applyUnfoldingShort(observedYieldsShape, *shapeYields, matrixFileNameShape);
      if (res) {
	flagFittingShape++;
	usedFiles.push_back(matrixFileNameShape);
	unfoldedYieldsShape.push_back(shapeYields);
      }
      else delete shapeYields;
    }
  }
  
  //
  TVectorD escaleFitShapeSystRelative(nMassBins);
  for(int i=0; i<nMassBins; i++){
    double min = unfoldedYields[i];
    double max = unfoldedYields[i];
    for (unsigned int iShape=0; iShape<unfoldedYieldsShape.size(); ++iShape) {
      min = TMath::Min( min, (*unfoldedYieldsShape[iShape])[i] );
      max = TMath::Max( max, (*unfoldedYieldsShape[iShape])[i] );
    }
    escaleFitShapeSystRelative[i] = fabs(max-min)/(max+min);
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
  for (unsigned int i=0; i<etaBinNames.size(); ++i) {
    TString fEtaFileName=TString("../root_files/yields/") + lumiTag + TString("_escale_eta/yields_bg-subtracted") + etaBinNames[i] + TString(".root");
    TFile fEta(fEtaFileName);
    if (fEta.IsOpen()) {
      usedFiles.push_back(fEtaFileName);
      TVectorD observedYieldsEta(nMassBins);
      observedYieldsEta = *(TVectorD*)file.Get("YieldsSignal");
      matrixFileName = TString("../root_files/constants/") + lumiTag + TString("_escale_eta/unfolding_constants") + etaBinNames[i] + TString(".root");
      TVectorD *unfYields = new TVectorD(nMassBins);
      res=applyUnfoldingShort(observedYieldsEta, *unfYields, matrixFileName);
      if (res==1) {
	flagEtaSyst++;
	usedFiles.push_back(matrixFileName);
	unfoldedYieldsEta.push_back(unfYields);
      }
      else {
	delete unfYields;
      }
    }
  }

  //
  TVectorD escaleEtaBinSystRelative(nMassBins);
  for(int i=0; i<nMassBins; i++){
    double min = unfoldedYields[i];
    double max = unfoldedYields[i];
    for (unsigned int iShape=0; iShape<unfoldedYieldsEta.size(); ++iShape) {
      min = TMath::Min( min, (*unfoldedYieldsEta[iShape])[i] );
      max = TMath::Max( max, (*unfoldedYieldsEta[iShape])[i] );
    }
    escaleEtaBinSystRelative[i] = fabs(max-min)/(max+min);
  }

  // 
  // Put all errors together
  //
  TVectorD escaleSystPercent(nMassBins);
  for(int i = 0; i < nMassBins; i++){
    escaleSystPercent[i] = 100.0 *
      sqrt(
	   escaleRandomizedSystRelative[i]*escaleRandomizedSystRelative[i]
	   + escaleShapeSystRelative[i]*escaleShapeSystRelative[i]
	   + escaleFitShapeSystRelative[i]*escaleFitShapeSystRelative[i]
	   + escaleEtaBinSystRelative[i]*escaleEtaBinSystRelative[i]
	   );
  }

  printf("\n\n");
  printf("mass bin     stat_part    shape-fit     residual    binning    total,%% \n");
  for(int i=0; i<nMassBins; i++){
    printf("%4.0f-%4.0f   %6.2f      %6.1f        %6.1f         %6.2f            %6.2f\n",
	   massBinLimits[i],massBinLimits[i+1],
	   escaleRandomizedSystRelative[i]*100.0,
	   escaleFitShapeSystRelative[i]*100.0,
	   escaleShapeSystRelative[i]*100.0,
	   escaleEtaBinSystRelative[i]*100.0,
	   escaleSystPercent[i]
	   );
  }

  std::cout << "\n\n";
  std::cout << std::string(80,'-') << "\n";
  std::cout << " Summary of successful loads:\n";
  std::cout << "  escale systematics   file count=  " << flagEScaleSyst << "\n";
  std::cout << "  residual shape syst. file count=  " << flagResidualShape << "\n";
  std::cout << "  eta systematics      file count=  " << flagEtaSyst << "\n";
  std::cout << "  fitting shape        file count=  " << flagFittingShape << "\n";
  std::cout << "\n";
  if (0) {
    std::cout << " Loaded files:\n";
    for (unsigned int i=0; i<usedFiles.size(); ++i) {
      std::cout << " " << usedFiles[i] << "\n";
    }
  }

  if (saveTexTable) {
    std::vector<TString> headers;
    std::vector<TVectorD*> data;
    std::vector<double> factors;
    headers.push_back("mass bin"); 
    headers.push_back("statistical, \\%"); data.push_back(&escaleRandomizedSystRelative); factors.push_back(100.);
    headers.push_back("shape, \\%"); data.push_back(&escaleFitShapeSystRelative); factors.push_back(100.);
    headers.push_back("residual, \\%"); data.push_back(&escaleShapeSystRelative); factors.push_back(100.);
    headers.push_back("$\\eta$ binning, \\%"); data.push_back(&escaleEtaBinSystRelative); factors.push_back(100.);
    headers.push_back("total, \\%"); data.push_back(&escaleSystPercent); factors.push_back(1.);
    
    TString texFName=TString("../root_files/systematics/") + lumiTag + TString("/escale_systematics_tmp.tex");
    printTexTable(texFName,nMassBins,massBinLimits,headers,data,factors);
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
int  applyUnfoldingShort(TVectorD &vin, TVectorD &vout, TString matrixFileName, int printLoadedData)
{

  // Read unfolding constants
  std::cout << "unfold: Load constants from <" << matrixFileName << ">" << std::endl;


  // Construct file names

  int res=unfolding::unfold(vin, vout, matrixFileName);
  if (res!=1) {
    std::cout << " ... in function applyUnfoldingShort\n";
    return 0;
  }

  // Print the result. Mainly for debugging purposes
  if (printLoadedData) {
    printf("\nUNFOLD: Results for the data, yields:\n");
    printf("                   yields observed        after unfolding            \n");
    for(int i=0; i<DYTools::nMassBins; i++){
      printf("%4.0f-%4.0f   %8.1f       %8.1f\n",
	     DYTools::massBinLimits[i], DYTools::massBinLimits[i+1],
	     vin[i], vout[i]);
    }
    printf("\n");
  }

  return 1;
}


//-----------------------------------------------------------------
// save texTable
//-----------------------------------------------------------------

int printTexTable(const TString &texFileName, int nMassBins, const double *massBinLimits, const std::vector<TString>& headers, const std::vector<TVectorD*> &data, const std::vector<double> &factors) {
  if ((headers.size()!=data.size()+1) || (headers.size()<2)) {
    std::cout << "printTexTable. vector size mismatch: " << headers.size() << " headers and " << data.size() << " data\n";
    return 0;
  }
  if (!massBinLimits) {
    std::cout << "printTexTable: massBinLimits is null\n";
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
  for(int i=0; i<nMassBins; i++) {
    fprintf(fout," %4.0lf $-$ %4.0lf ",massBinLimits[i],massBinLimits[i+1]);
    for (unsigned int j=0; j<headers.size()-1; ++j) {
      fprintf(fout," & %5.2lf ",(*data[j])[i]*factors[j]);
    }
    fprintf(fout,"\\\\\n");
  }
  fprintf(fout,"\\hline\n");
  fprintf(fout,"\\end{tabular}\n\\end{center}\n\\end{table}\n\n");
  fclose(fout);
  std::cout << "texFileName=<" << texFileName << "> saved\n";
  return 1;
}

