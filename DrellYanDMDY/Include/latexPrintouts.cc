#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TMatrixD.h> 
#include <TString.h> 

#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings

#include "../Include/DYTools.hh"        // helper class for organizing input ntuple files
#include "../Include/latexPrintouts.hh"
     
#endif

// -----------------------------------------------------------------------------

void latexPrintoutAcceptance2D(TMatrixD accv, TMatrixD accErrv, TString producedBy)
{
   latexPrintoutOneValue2D(accv, accErrv, producedBy, "acceptance  \\%", "acceptance2D" , "Numerical values of the post-FSR acceptance for %4.0f-%4.0f GeV mass slice of \\DYee candidates");
   //TString baseOfReferenceName="acceptance2D";
   //TString tableName = "Numerical values of the post-FSR acceptance for %4.0f-%4.0f GeV mass slice of \\DYee candidates";
   //TString valueName = "acceptance"
}


void latexPrintoutAcceptance1D(TMatrixD accv, TMatrixD accErrv, TString producedBy)
{
   latexPrintoutOneValue1D(accv, accErrv, producedBy, "acceptance  \\%", "acceptance1D" , "Numerical values of the post-FSR acceptance for 1D measurement of \\DYee candidates" );
   //TString baseOfReferenceName="acceptance1D";
   //TString tableName = "Numerical values of the post-FSR acceptance for 1D measurement of \\DYee candidates";
   //TString valueName = "acceptance"

}

void latexPrintoutEfficiency2D(TMatrixD effv, TMatrixD effErrv, TString producedBy)
{
   latexPrintoutOneValue2D(effv, effErrv, producedBy, "efficiency  \\%", "efficiency2D" , "Reconstruction and selection efficiency $\\epsilon^{mc}$ for %4.0f-%4.0f GeV mass slice of \\DYee candidates");
   //TString baseOfReferenceName="efficiency2D";
   //TString tableName = "Reconstruction and selection efficiency $\\epsilon^{mc}$ for %4.0f-%4.0f GeV mass slice of \\DYee candidates";
   //TString valueName = "efficiency"
}


void latexPrintoutEfficiency1D(TMatrixD effv, TMatrixD effErrv, TString producedBy)
{
   latexPrintoutOneValue1D(effv, effErrv, producedBy, "efficiency \\%", "efficiency1D" , "Reconstruction and selection efficiency $\\epsilon^{mc}$ of \\DYee candidates" );
   //TString baseOfReferenceName="event-sf1D";
   //TString tableName = "Reconstruction and selection efficiency $\\epsilon^{mc}$ of \\DYee candidates";
   //TString valueName = "$\rho_{data/mc}$"
}

void latexPrintoutScaleFactors2D(TMatrixD scalev, TMatrixD scaleErrv, TString producedBy)
{
   latexPrintoutOneValue2D(scalev, scaleErrv, producedBy, "$\rho_{data/mc}$", "event-sf2D" , "Scale factors for correcting MC event efficiency for %4.0f-%4.0f GeV mass slice of \\DYee candidates");
   //TString baseOfReferenceName="event-sf2D";
   //TString tableName = "Scale factors for correcting MC event efficiency for %4.0f-%4.0f GeV mass slice of \\DYee candidates";
   //TString valueName = "$\rho_{data/mc}$"
}

void latexPrintoutScaleFactors1D(TMatrixD scalev, TMatrixD scaleErrv, TString producedBy)
{
   latexPrintoutOneValue1D(scalev, scaleErrv, producedBy, "$\rho_{data/mc}$", "event-sf1D" , "Scale factors for correcting MC event efficiency" );
   //TString baseOfReferenceName="event-sf1D";
   //TString tableName = "Scale factors for correcting MC event efficiency";
   //TString valueName = "$\rho_{data/mc}$"
}

void latexPrintoutOneValue2D(TMatrixD value, TMatrixD valueErr, TString producedBy, TString valueName, TString baseOfReferenceName, TString tableName)
{
   FILE* txtFile;
   TString valueNameForSaving=valueName;
   valueNameForSaving.ReplaceAll(" ","");
   valueNameForSaving.ReplaceAll(".","");
   valueNameForSaving.ReplaceAll("/","");
   valueNameForSaving.ReplaceAll("\\","");
   valueNameForSaving.ReplaceAll("\"","");
   valueNameForSaving.ReplaceAll(",","");
   valueNameForSaving.ReplaceAll(":","");
   valueNameForSaving.ReplaceAll(";","");
   valueNameForSaving.ReplaceAll("%","");
   TString txtFileName="tables2D"+valueNameForSaving+".txt";
   txtFile = fopen(txtFileName,"w");
   TString strComment="% Tables produced by "+ producedBy + " in 2D measurement\n";
   fprintf(txtFile,strComment);  
   fprintf(txtFile,"\n\n\n");
   for (int mslice=0; mslice<DYTools::nMassBins; mslice++) 
     {   
       fprintf(txtFile,"\\begin{table}\n");
       TString strHeader="\\caption{\\label{tab:" + baseOfReferenceName+"-%d} " + tableName + "}\n";
       fprintf(txtFile,strHeader,mslice, DYTools::massBinLimits[mslice], DYTools::massBinLimits[mslice+1] ); 
       fprintf(txtFile,"\\begin{center}\\small{\n");
       fprintf(txtFile,"\\begin{tabular}[2]{|c|c||c|c|}\n");
       fprintf(txtFile,"\\hline\n");
       TString strQuantities = " $|Y|$  &   " + valueName + "  & $|Y|$  &   " + valueName + " \\\\ \n";
       fprintf(txtFile,strQuantities);
       fprintf(txtFile,"\\hline\n");
       int halfBins=DYTools::nYBins[mslice]/2;
       for (int j=0; j<halfBins; j++)
         {
           fprintf(txtFile,"%1.1f-%1.1f &", j*(DYTools::yRangeMax-DYTools::yRangeMin)/DYTools::nYBins[mslice], (j+1)*(DYTools::yRangeMax-DYTools::yRangeMin)/DYTools::nYBins[mslice]);
           fprintf(txtFile,"$ %7.4f \\pm %6.4f $ & ", value(mslice,j), valueErr(mslice,j));
           if (halfBins+j<DYTools::nYBins[mslice])
             {
               fprintf(txtFile,"%1.1f-%1.1f &", (halfBins+j)*(DYTools::yRangeMax-DYTools::yRangeMin)/DYTools::nYBins[mslice], (halfBins+j+1)*(DYTools::yRangeMax-DYTools::yRangeMin)/DYTools::nYBins[mslice]); 
               fprintf(txtFile,"$ %7.4f \\pm %6.4f $ \\\\ \n ", value(mslice,halfBins+j), valueErr(mslice,halfBins+j));
             }
           else fprintf(txtFile,"& \\\\ \n ");
        }
        fprintf(txtFile,"\\hline\n");
        fprintf(txtFile,"\\end{tabular}}\n");
        fprintf(txtFile,"\\end{center}\n");
        fprintf(txtFile,"\\end{table}\n");
        fprintf(txtFile,"\n\n");
     }
   fclose(txtFile);
}


void latexPrintoutOneValue1D(TMatrixD value, TMatrixD valueErr, TString producedBy, TString valueName, TString baseOfReferenceName, TString tableName)
{
   FILE* txtFile;
   TString valueNameForSaving=valueName;
   valueNameForSaving.ReplaceAll(" ","");
   valueNameForSaving.ReplaceAll(".","");
   valueNameForSaving.ReplaceAll("/","");
   valueNameForSaving.ReplaceAll("\\","");
   valueNameForSaving.ReplaceAll("\"","");
   valueNameForSaving.ReplaceAll(",","");
   valueNameForSaving.ReplaceAll(":","");
   valueNameForSaving.ReplaceAll(";","");
   valueNameForSaving.ReplaceAll("%","");
   TString txtFileName="tables1D"+valueNameForSaving+".txt";
   txtFile = fopen(txtFileName,"w");
   TString strComment="% Tables produced by "+ producedBy + " in 1D measurement\n";
   fprintf(txtFile,strComment); 
   fprintf(txtFile,"\n"); 
 
   fprintf(txtFile,"\\begin{table}\n");
   TString strHeader="\\caption{\\label{tab:" + baseOfReferenceName+"} " + tableName + "}\n"; 
   fprintf(txtFile,strHeader); 
   fprintf(txtFile,"\\begin{center}\\small{\n");
   fprintf(txtFile,"\\begin{tabular}[2]{|c|c||c|c|}\n");
   fprintf(txtFile,"\\hline\n");
   TString strQuantities = " $|Y|$  &   " + valueName + " & $|Y|$  &   " + valueName + " \\\\ \n";
   fprintf(txtFile,strQuantities);
   fprintf(txtFile,"\\hline\n");
   int halfBins=DYTools::nMassBins/2;
   for (int mslice=0; mslice<halfBins; mslice++)
     {
       fprintf(txtFile,"%4.0f-%4.0f &", DYTools::massBinLimits[mslice], DYTools::massBinLimits[mslice+1] );
       fprintf(txtFile,"$ %7.4f \\pm %6.4f $ & ", value(mslice,0), valueErr(mslice,0));
       if (halfBins+mslice<DYTools::nMassBins)
         {
           fprintf(txtFile,"%4.0f-%4.0f &", DYTools::massBinLimits[halfBins+mslice], DYTools::massBinLimits[halfBins+mslice+1]); 
           fprintf(txtFile,"$ %7.4f \\pm %6.4f $ \\\\ \n ", value(mslice+halfBins,0), valueErr(mslice+halfBins,0));
         }
       else fprintf(txtFile,"& \\\\ \n ");
     }
   fprintf(txtFile,"\\hline\n");
   fprintf(txtFile,"\\end{tabular}}\n");
   fprintf(txtFile,"\\end{center}\n");
   fprintf(txtFile,"\\end{table}\n");
   fprintf(txtFile,"\n\n");
   fclose(txtFile);
}
