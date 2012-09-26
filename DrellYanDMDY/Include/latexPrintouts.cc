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
}


void latexPrintoutAcceptance1D(TMatrixD accv, TMatrixD accErrv, TString producedBy)
{
   latexPrintoutOneValue1D(accv, accErrv, producedBy, "acceptance  \\%", "acceptance1D" , "Numerical values of the post-FSR acceptance for 1D measurement of \\DYee candidates" );
}

void latexPrintoutEfficiency2D(TMatrixD effv, TMatrixD effErrv, TString producedBy)
{
   latexPrintoutOneValue2D(effv, effErrv, producedBy, "efficiency  \\%", "efficiency2D" , "Reconstruction and selection efficiency $\\epsilon^{mc}$ for %4.0f-%4.0f GeV mass slice of \\DYee candidates");
}


void latexPrintoutEfficiency1D(TMatrixD effv, TMatrixD effErrv, TString producedBy)
{
   latexPrintoutOneValue1D(effv, effErrv, producedBy, "efficiency \\%", "efficiency1D" , "Reconstruction and selection efficiency $\\epsilon^{mc}$ of \\DYee candidates" );
}

void latexPrintoutScaleFactors2D(TMatrixD scalev, TMatrixD scaleErrv, TString producedBy)
{
   latexPrintoutOneValue2D(scalev, scaleErrv, producedBy, "$\rho_{data/mc}$", "event-sf2D" , "Scale factors for correcting MC event efficiency for %4.0f-%4.0f GeV mass slice of \\DYee candidates");
}

void latexPrintoutScaleFactors1D(TMatrixD scalev, TMatrixD scaleErrv, TString producedBy)
{
   latexPrintoutOneValue1D(scalev, scaleErrv, producedBy, "$\rho_{data/mc}$", "event-sf1D" , "Scale factors for correcting MC event efficiency" );
}

void latexPrintoutFsr2D(TMatrixD corrv, TMatrixD corrErrv, TString producedBy)
{
   latexPrintoutOneValue2D(corrv, corrErrv, producedBy, "FSR \\%", "fsr-binbybin-2D" , "Numerical values of the Fsr corrections in full phase space for %4.0f-%4.0f GeV mass slice of \\DYee candidates");
}

void latexPrintoutFsr1D(TMatrixD corrv, TMatrixD corrErrv, TString producedBy)
{
   latexPrintoutOneValue1D(corrv, corrErrv, producedBy, "FSR \\%", "fsr-binbybin-1D" , "Numerical values of the Fsr corrections in full phase space of \\DYee candidates");
}

void latexPrintoutFsrInAcceptance2D(TMatrixD corrv, TMatrixD corrErrv, TString producedBy)
{
   latexPrintoutOneValue2D(corrv, corrErrv, producedBy, "FSR in acceptance \\%", "fsrInAcc-binbybin-2D" , "Numerical values of the Fsr corrections in detector phase space  (i.e within acceptance) for %4.0f-%4.0f GeV mass slice of \\DYee candidates");
}

void latexPrintoutFsrInAcceptance1D(TMatrixD corrv, TMatrixD corrErrv, TString producedBy)
{
   latexPrintoutOneValue1D(corrv, corrErrv, producedBy, "FSR in acceptance \\%", "fsrInAcc-binbybin-2D" , "Numerical values of the Fsr corrections in full phase space of \\DYee candidates");
}

void latexPrintoutOneValue2D(TMatrixD value, TMatrixD valueErr, TString producedBy, TString valueName, TString baseOfReferenceName, TString tableName)
{

   int nValues=1;
   TMatrixD* values[1];   values[0]=&value;
   int valuesType[1];   valuesType[0]=1;
   TMatrixD* valuesErr1[1];   valuesErr1[0]=&valueErr;
   TMatrixD* valuesErr2[1];   valuesErr2[0]=0;
   TString valuesName[1];   valuesName[0]=valueName;
   latexPrintoutTwoColumns2D(nValues, valuesType, values, valuesErr1, valuesErr2, producedBy, valuesName, baseOfReferenceName, tableName);

}


void latexPrintoutOneValue1D(TMatrixD value, TMatrixD valueErr, TString producedBy, TString valueName, TString baseOfReferenceName, TString tableName)
{
 
   int nValues=1;
   TMatrixD* values[1];   values[0]=&value;
   int valuesType[1];   valuesType[0]=1;
   TMatrixD* valuesErr1[1];   valuesErr1[0]=&valueErr;
   TMatrixD* valuesErr2[1];   valuesErr2[0]=0;
   TString valuesName[1];   valuesName[0]=valueName;
   latexPrintoutTwoColumns1D(nValues, valuesType, values, valuesErr1, valuesErr2, producedBy, valuesName, baseOfReferenceName, tableName);
}

void latexPrintoutTwoColumns2D(const int nValues, int* valuesType, TMatrixD** values, TMatrixD** valuesErr1, TMatrixD** valuesErr2, TString producedBy, TString* valuesName, TString baseOfReferenceName, TString tableName)
{
   FILE* txtFile;
   TString valueNameForSaving=valuesName[0];
   for (int i=1; i<nValues; i++)
     valueNameForSaving+=valuesName[i];
   valueNameForSaving.ReplaceAll(" ","");
   valueNameForSaving.ReplaceAll(".","");
   valueNameForSaving.ReplaceAll("/","");
   valueNameForSaving.ReplaceAll("\"","");
   valueNameForSaving.ReplaceAll("\\","");
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

       TString strBeginTabular="\\begin{tabular}[2]{|c|";
       for (int i=0; i<nValues; i++)
         strBeginTabular+="c|";
       strBeginTabular+="|c|";
       for (int i=0; i<nValues; i++)
         strBeginTabular+="c|";
       strBeginTabular+="}\n";   
       fprintf(txtFile,strBeginTabular);

       fprintf(txtFile,"\\hline\n");

       TString strQuantities = " $|Y|$ ";
       for (int i=0; i<nValues; i++)
         {
           strQuantities+=" & ";
           strQuantities+=valuesName[i];       
         }   
       strQuantities += " $|Y|$  ";
       for (int i=0; i<nValues; i++)
         {
           strQuantities+=" & ";
           strQuantities+=valuesName[i];       
         } 
       strQuantities +=  " \\\\ \n";      
       fprintf(txtFile,strQuantities);
       fprintf(txtFile,"\\hline\n");

       int halfBins=DYTools::nYBins[mslice]/2;
       for (int j=0; j<halfBins; j++)
         {
           fprintf(txtFile,"%1.1f-%1.1f", j*(DYTools::yRangeMax-DYTools::yRangeMin)/DYTools::nYBins[mslice], (j+1)*(DYTools::yRangeMax-DYTools::yRangeMin)/DYTools::nYBins[mslice]);
           for (int i=0; i<nValues; i++)
             {
                TMatrixD& temp0= *values[i];
                TMatrixD& temp1= *valuesErr1[i];
                TMatrixD& temp2= *valuesErr2[i];
                if (valuesType[i]==0) 
                  fprintf(txtFile," & $ %7.4f$ ", temp0(mslice,j));
                else if (valuesType[i]==1) 
                  fprintf(txtFile," & $ %7.4f \\pm %6.4f $ ", temp0(mslice,j), temp1(mslice,j));
                else if (valuesType[i]==2) 
                  fprintf(txtFile," & $ %7.4f \\pm %6.4f \\pm %6.4f $", temp0(mslice,j), temp1(mslice,j), temp2(mslice,j));
             }

           fprintf(txtFile,"&");
           if (halfBins+j<DYTools::nYBins[mslice])
             {
               fprintf(txtFile,"%1.1f-%1.1f", (halfBins+j)*(DYTools::yRangeMax-DYTools::yRangeMin)/DYTools::nYBins[mslice], (halfBins+j+1)*(DYTools::yRangeMax-DYTools::yRangeMin)/DYTools::nYBins[mslice]); 

               for (int i=0; i<nValues; i++)
                 {
                   TMatrixD& temp0= *values[i];
                   TMatrixD& temp1= *valuesErr1[i];
                   TMatrixD& temp2= *valuesErr2[i];
                   if (valuesType[i]==0) 
                     fprintf(txtFile," & $ %7.4f$ ", temp0(mslice,j));
                   else if (valuesType[i]==1) 
                     fprintf(txtFile," & $ %7.4f \\pm %6.4f $ ", temp0(mslice,j), temp1(mslice,j));
                   else if (valuesType[i]==2) 
                     fprintf(txtFile," & $ %7.4f \\pm %6.4f \\pm %6.4f $", temp0(mslice,j), temp1(mslice,j), temp2(mslice,j));
                 }

               fprintf(txtFile,"& \\\\ \n ");
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

void latexPrintoutTwoColumns1D(const int nValues, int* valuesType, TMatrixD** values, TMatrixD** valuesErr1, TMatrixD** valuesErr2, TString producedBy, TString* valuesName, TString baseOfReferenceName, TString tableName)
{
   //valueTypes: 0 - value; 1 - value+-error; 2 - value+=error1+-error2
   FILE* txtFile;
   TString valueNameForSaving=valuesName[0];
   for (int i=1; i<nValues; i++)
     valueNameForSaving+=valuesName[i];
   valueNameForSaving.ReplaceAll(" ","");
   valueNameForSaving.ReplaceAll(".","");
   valueNameForSaving.ReplaceAll("/","");
   valueNameForSaving.ReplaceAll("\"","");
   valueNameForSaving.ReplaceAll("\\","");
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

   TString strBeginTabular="\\begin{tabular}[2]{|c|";
   for (int i=0; i<nValues; i++)
     strBeginTabular+="c|";
   strBeginTabular+="|c|";
   for (int i=0; i<nValues; i++)
     strBeginTabular+="c|";
   strBeginTabular+="}\n";   
   fprintf(txtFile,strBeginTabular);

   fprintf(txtFile,"\\hline\n");
   TString strQuantities = " mass, GeV  ";
   for (int i=0; i<nValues; i++)
     {
       strQuantities+=" & ";
       strQuantities+=valuesName[i];       
     }   
   strQuantities += " & mass, GeV  ";
   for (int i=0; i<nValues; i++)
     {
       strQuantities+=" & ";
       strQuantities+=valuesName[i];       
     } 
   strQuantities +=  " \\\\ \n";      
   fprintf(txtFile,strQuantities);

   fprintf(txtFile,"\\hline\n");
   int halfBins=DYTools::nMassBins/2;
   for (int mslice=0; mslice<halfBins; mslice++)
     {
       fprintf(txtFile,"%4.0f-%4.0f", DYTools::massBinLimits[mslice], DYTools::massBinLimits[mslice+1] );
       for (int i=0; i<nValues; i++)
         {
            TMatrixD& temp0= *values[i];
            TMatrixD& temp1= *valuesErr1[i];
            TMatrixD& temp2= *valuesErr2[i];
            if (valuesType[i]==0) 
              fprintf(txtFile," & $ %7.4f$ ", temp0(mslice,0));
            else if (valuesType[i]==1) 
              fprintf(txtFile," & $ %7.4f \\pm %6.4f $ ", temp0(mslice,0), temp1(mslice,0));
            else if (valuesType[i]==2) 
              fprintf(txtFile," & $ %7.4f \\pm %6.4f \\pm %6.4f $ ", temp0(mslice,0), temp1(mslice,0), temp2(mslice,0));
         }

       fprintf(txtFile,"&");

       if (halfBins+mslice<DYTools::nMassBins)
         {
           fprintf(txtFile,"%4.0f-%4.0f", DYTools::massBinLimits[halfBins+mslice], DYTools::massBinLimits[halfBins+mslice+1]); 
           for (int i=0; i<nValues; i++)
             {
               TMatrixD& temp0= *values[i];
               TMatrixD& temp1= *valuesErr1[i];
               TMatrixD& temp2= *valuesErr2[i];
               if (valuesType[i]==0) 
                 fprintf(txtFile," & $ %7.4f$ ", temp0(mslice+halfBins,0));
               else if (valuesType[i]==1) 
                 fprintf(txtFile," & $ %7.4f \\pm %6.4f $ ", temp0(mslice+halfBins,0), temp1(mslice+halfBins,0));
               else if (valuesType[i]==2) 
                 fprintf(txtFile," & $ %7.4f \\pm %6.4f \\pm %6.4f $ ", temp0(mslice+halfBins,0), temp1(mslice+halfBins,0), temp2(mslice+halfBins,0));
             }
           fprintf(txtFile,"& \\\\ \n ");

           //fprintf(txtFile,"$ %7.4f \\pm %6.4f $ \\\\ \n ", value(mslice+halfBins,0), valueErr(mslice+halfBins,0));
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
