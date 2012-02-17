#include "TVectorD.h"
#include "TFile.h"
#include "TString.h"
#include "TSystem.h"
#include <iostream>
#include <fstream>
#include <sstream>

#include "../Include/DYTools.hh"
#include "../Include/UnfoldingTools.hh"

void  applyUnfoldingShort(TVectorD &vin, TVectorD &vout, TString matrixFileName);

void calcEscaleSystematics(){

  TVectorD observedYields(nMassBins);
  TVectorD unfoldedYields(nMassBins);

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

  TString matrixFileName = 
    "../root_files/constants/DY_m10+pr_1142pb/unfolding_constants.root";
  const int nFiles1 = 20;
  for(int ifile=0; ifile<nFiles1; ifile++){
    int seed = 1001+ifile;
    TString fname = "../root_files/yields/DY_m10+pr_1142pb_escale_randomized/yields_bg-subtracted_seed";
    fname += seed;
    fname += ".root";
    TFile file(fname);
    observedYields = *(TVectorD*)file.Get("YieldsSignal");
    applyUnfoldingShort(observedYields, unfoldedYields,matrixFileName);
    unfoldedYields.Print();
    // Accumulate mean and RMS
    for(int iMassBin = 0; iMassBin < nMassBins; iMassBin++){
      unfoldedYieldsMean[iMassBin] += unfoldedYields[iMassBin];
      unfoldedYieldsSquaredMean[iMassBin] += unfoldedYields[iMassBin]*unfoldedYields[iMassBin];
    }
    file.Close();
  }

  // Final calculation of the mean and RMS for Smearing
  TVectorD escaleRandomizedSystRelative(nMassBins);
  for(int i = 0; i < nMassBins; i++){
    unfoldedYieldsMean[i] = unfoldedYieldsMean[i]/(1.0*nFiles1);
    unfoldedYieldsSquaredMean[i] = unfoldedYieldsSquaredMean[i]/(1.0*nFiles1);
    unfoldedYieldsRMS[i] = sqrt(unfoldedYieldsSquaredMean[i] - 
				unfoldedYieldsMean[i]*unfoldedYieldsMean[i]);
    escaleRandomizedSystRelative[i] = unfoldedYieldsRMS[i]/unfoldedYieldsMean[i];
  }

  escaleRandomizedSystRelative.Print();

  //
  // Calculate error related to the difference between
  // data and MC mass distribution shape after all corrections
  // had been applied to data and MC
  //
  TVectorD unfoldedYieldsVariation(nMassBins);
  TString fname = "../root_files/yields/DY_m10+pr_1142pb/yields_bg-subtracted.root";
  TFile file(fname);
  observedYields = *(TVectorD*)file.Get("YieldsSignal");
  //
  matrixFileName = "../root_files/constants/DY_m10+pr_1142pb/unfolding_constants.root";
  applyUnfoldingShort(observedYields, unfoldedYields, matrixFileName);
  //
  matrixFileName = "../root_files/constants/DY_m10+pr_1142pb_shape_syst/unfolding_constants.root";
  applyUnfoldingShort(observedYields, unfoldedYieldsVariation, matrixFileName);
  TVectorD escaleShapeSystRelative(nMassBins);
  for(int i=0; i<nMassBins; i++){
    escaleShapeSystRelative[i] = fabs(unfoldedYields[i] - unfoldedYieldsVariation[i])
      /(unfoldedYields[i] + unfoldedYieldsVariation[i]);
  }

  //
  // Calculate error related to extra smearing function shape
  //
  TVectorD observedYieldsBW(nMassBins);
  TFile fBW("../root_files/yields/DY_m10+pr_1142pb_6binNegs_BreitWigner/yields_bg-subtracted.root");
  observedYieldsBW = *(TVectorD*)file.Get("YieldsSignal");
  TString matrixFileNameBW = "../root_files/constants/DY_m10+pr_1142pb_6binNegs_BreitWigner/unfolding_constants.root";
  TVectorD unfoldedYieldsBW(nMassBins);
  applyUnfoldingShort(observedYieldsBW, unfoldedYieldsBW, matrixFileNameBW);
  //
  TVectorD observedYieldsVT(nMassBins);
  TFile fVT("../root_files/yields/DY_m10+pr_1142pb_6binNegs_Voigtian/yields_bg-subtracted.root");
  observedYieldsVT = *(TVectorD*)file.Get("YieldsSignal");
  TString matrixFileNameVT = "../root_files/constants/DY_m10+pr_1142pb_6binNegs_Voigtian/unfolding_constants.root";
  TVectorD unfoldedYieldsVT(nMassBins);
  applyUnfoldingShort(observedYieldsVT, unfoldedYieldsVT, matrixFileNameVT);
  //
  TVectorD escaleFitShapeSystRelative(nMassBins);
  for(int i=0; i<nMassBins; i++){
    double max = TMath::Max(TMath::Max(unfoldedYields[i],unfoldedYieldsBW[i]),
			    unfoldedYieldsVT[i]);
    double min = TMath::Min(TMath::Min(unfoldedYields[i],unfoldedYieldsBW[i]),
			    unfoldedYieldsVT[i]);
    escaleFitShapeSystRelative[i] = fabs(max-min)/(max+min);
  }

  //
  // Calculate error related to eta binning
  //
  TVectorD observedYields2bin(nMassBins);
  TFile f2bin("../root_files/yields/DY_m10+pr_1142pb_6binNegs_2binNegs_Gauss/yields_bg-subtracted.root");
  observedYields2bin = *(TVectorD*)file.Get("YieldsSignal");
  TString matrixFileName2bin = "../root_files/constants/DY_m10+pr_1142pb_2binNegs_Gauss/unfolding_constants.root";
  TVectorD unfoldedYields2bin(nMassBins);
  applyUnfoldingShort(observedYields2bin, unfoldedYields2bin, matrixFileName2bin);
  //
  TVectorD observedYields4bin(nMassBins);
  TFile f4bin("../root_files/yields/DY_m10+pr_1142pb_4binNegs_Gauss/yields_bg-subtracted.root");
  observedYields4bin = *(TVectorD*)file.Get("YieldsSignal");
  TString matrixFileName4bin = "../root_files/constants/DY_m10+pr_1142pb_4binNegs_Gauss/unfolding_constants.root";
  TVectorD unfoldedYields4bin(nMassBins);
  applyUnfoldingShort(observedYields4bin, unfoldedYields4bin, matrixFileName4bin);
  //
  TVectorD escaleBinSystRelative(nMassBins);
  for(int i=0; i<nMassBins; i++){
    double max = TMath::Max(TMath::Max(unfoldedYields[i],unfoldedYields2bin[i]),
			    unfoldedYields4bin[i]);
    double min = TMath::Min(TMath::Min(unfoldedYields[i],unfoldedYields2bin[i]),
			    unfoldedYields4bin[i]);
    escaleBinSystRelative[i] = fabs(max-min)/(max+min);
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
	   + escaleBinSystRelative[i]*escaleBinSystRelative[i]
	   );
  }

  printf("mass bin     stat_part    shape-fit   shape-reweight    binning    total,%% \n");
  for(int i=0; i<nMassBins; i++){
    printf("%4.0f-%4.0f   %6.2f      %6.1f        %6.1f         %6.2f            %6.2f\n",
	   massBinLimits[i],massBinLimits[i+1],
	   escaleRandomizedSystRelative[i]*100.0,
	   escaleFitShapeSystRelative[i]*100.0,
	   escaleShapeSystRelative[i]*100.0,
	   escaleBinSystRelative[i]*100.0,
	   escaleSystPercent[i]
	   );
  }

  TFile fout("../root_files/systematics/DY_m10+pr_1142pb/escale_systematics_tmp.root","recreate");
  escaleSystPercent.Write("escaleSystPercent");
  fout.Close();

  return;
}

//-----------------------------------------------------------------
// Unfold
//-----------------------------------------------------------------
void  applyUnfoldingShort(TVectorD &vin, TVectorD &vout, TString matrixFileName)
{

  // Read unfolding constants
  printf("unfold: Load constants\n"); fflush(stdout);


  // Construct file names

  unfolding::unfold(vin, vout, matrixFileName);

  // Print the result. Mainly for debugging purposes
  printf("\nUNFOLD: Results for the data, yields:\n");
  printf("                   yields observed        after unfolding            \n");
  for(int i=0; i<DYTools::nMassBins; i++){
    printf("%4.0f-%4.0f   %8.1f       %8.1f\n",
	   DYTools::massBinLimits[i], DYTools::massBinLimits[i+1],
	   vin[i], vout[i]);
  }
  printf("\n");

  return;
}

