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

void latexPrintoutAcceptance2D(TMatrixD accv, TMatrixD accErrv)
{
   FILE* txtFile;
   TString txtFileName;
   if (DYTools::study2D) txtFileName="tables2D.txt";
   else txtFileName="tables1D.txt";
   txtFile = fopen(txtFileName,"w");
   fprintf(txtFile,"%Table of Acceptance in 2D measurement\n\n\n");  
   for (int mslice=0; mslice<DYTools::nMassBins; mslice++) 
     {   
       fprintf(txtFile,"\\begin{table}\n");
       fprintf(txtFile,"\\caption{\\label{tab:acceptance2D-%d} Numerical values of the post-FSR acceptance for %4.0f-%4.0f GeV mass slice of \\DYee candidates}\n",mslice, DYTools::massBinLimits[mslice], DYTools::massBinLimits[mslice+1] ); 
       fprintf(txtFile,"\\begin{center}\\small{\n");
       fprintf(txtFile,"\\begin{tabular}[2]{|c|c||c|c|}\n");
       fprintf(txtFile,"\\hline\n");
       fprintf(txtFile," $|Y|$  &   acceptance \\%  & $|Y|$  &   acceptance \\% \\\\ \n");
       fprintf(txtFile,"\\hline\n");
       int halfBins=DYTools::nYBins[mslice]/2;
       for (int j=0; j<halfBins; j++)
         {
           fprintf(txtFile,"%1.1f-%1.1f &", j*(DYTools::yRangeMax-DYTools::yRangeMin)/DYTools::nYBins[mslice], (j+1)*(DYTools::yRangeMax-DYTools::yRangeMin)/DYTools::nYBins[mslice]);
           fprintf(txtFile,"$ %7.4f \\pm %6.4f $ & ", accv(mslice,j), accErrv(mslice,j));
           if (halfBins+j<DYTools::nYBins[mslice])
             {
               fprintf(txtFile,"%1.1f-%1.1f &", (halfBins+j)*(DYTools::yRangeMax-DYTools::yRangeMin)/DYTools::nYBins[mslice], (halfBins+j+1)*(DYTools::yRangeMax-DYTools::yRangeMin)/DYTools::nYBins[mslice]); 
               fprintf(txtFile,"$ %7.4f \\pm %6.4f $ \\\\ \n ", accv(mslice,halfBins+j), accErrv(mslice,halfBins+j));
             }
           else fprintf(txtFile,"& \\\\ \n ");
        }
        fprintf(txtFile,"\\hline\n");
        fprintf(txtFile,"\\end{tabular}}\n");
        fprintf(txtFile,"\\end{center}\n");
        fprintf(txtFile,"\\end{table}\n");
        fprintf(txtFile,"\n\n");
     }
}

