#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>
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
#include "../Include/MyTools.hh"
     
#endif

// -----------------------------------------------------------------------------

void latexPrintoutAcceptance2D(TMatrixD accv, TMatrixD accErrv, TString producedBy)
{
   latexPrintoutOneValue2D(accv, accErrv, producedBy, 
                           "acceptance", "acceptance2D" , 
                           "Numerical values of the post-FSR acceptance for %4.0f-%4.0f GeV mass slice of \\DYee candidates");
}


void latexPrintoutAcceptance1D(TMatrixD accv, TMatrixD accErrv, TString producedBy)
{
   latexPrintoutOneValue1D(accv, accErrv, producedBy, 
                           "acceptance", "acceptance1D" , 
                           "Numerical values of the post-FSR acceptance for 1D measurement of \\DYee candidates" );
}

void latexPrintoutEfficiency2D(TMatrixD effv, TMatrixD effErrv, TString producedBy)
{
   latexPrintoutOneValue2D(effv, effErrv, producedBy, 
                           "efficiency", "efficiency2D" , 
                           "Reconstruction and selection efficiency $\\epsilon^{mc}$ for %4.0f-%4.0f GeV mass slice of \\DYee candidates");
}


void latexPrintoutEfficiency1D(TMatrixD effv, TMatrixD effErrv, TString producedBy)
{
   latexPrintoutOneValue1D(effv, effErrv, producedBy, 
                           "efficiency", "efficiency1D" , 
                           "Reconstruction and selection efficiency $\\epsilon^{mc}$ of \\DYee candidates" );
}

void latexPrintoutScaleFactors2D(TMatrixD scalev, TMatrixD scaleErrv, TString producedBy)
{
   latexPrintoutOneValue2D(scalev, scaleErrv, producedBy, 
                           "$\rho_{data/mc}$", "event-sf2D" , 
                           "Scale factors for correcting MC event efficiency for %4.0f-%4.0f GeV mass slice of \\DYee candidates");
}

void latexPrintoutScaleFactors1D(TMatrixD scalev, TMatrixD scaleErrv, TString producedBy)
{
   latexPrintoutOneValue1D(scalev, scaleErrv, producedBy, 
                           "$\rho_{data/mc}$", "event-sf1D" , 
                           "Scale factors for correcting MC event efficiency" );
}

void latexPrintoutFsr2D(TMatrixD corrv, TMatrixD corrErrv, TString producedBy)
{
   latexPrintoutOneValue2D(corrv, corrErrv, producedBy, 
                           "FSR", "fsr-binbybin-2D" , 
                           "Numerical values of the Fsr corrections in full phase space for %4.0f-%4.0f GeV mass slice of \\DYee candidates");
}

void latexPrintoutFsr1D(TMatrixD corrv, TMatrixD corrErrv, TString producedBy)
{
   latexPrintoutOneValue1D(corrv, corrErrv, producedBy, 
                           "FSR", "fsr-binbybin-1D" , 
                           "Numerical values of the Fsr corrections in full phase space of \\DYee candidates");
}

void latexPrintoutFsrInAcceptance2D(TMatrixD corrv, TMatrixD corrErrv, TString producedBy)
{
   latexPrintoutOneValue2D(corrv, corrErrv, producedBy, "FSR in acceptance", 
                           "fsrInAcc-binbybin-2D" , 
                           "Numerical values of the Fsr corrections in detector phase space  (i.e within acceptance) for %4.0f-%4.0f GeV mass slice of \\DYee candidates");
}

void latexPrintoutFsrInAcceptance1D(TMatrixD corrv, TMatrixD corrErrv, TString producedBy)
{
   latexPrintoutOneValue1D(corrv, corrErrv, producedBy, "FSR in acceptance", 
                            "fsrInAcc-binbybin-2D" , 
                            "Numerical values of the Fsr corrections in full phase space of \\DYee candidates");
}

void latexPrintoutBackgroundRates2D(TMatrixD observedYields, TMatrixD observedYieldsErr, 
                                    TMatrixD totalBackground, TMatrixD totalBackgroundError, 
                                    TMatrixD totalBackgroundErrorSyst, TMatrixD bkgRatesUsual, 
                                    TString producedBy)
{
   int nValues=3;
   TMatrixD* values[nValues];   values[0]=&observedYields; 
   values[1]=&totalBackground; values[2]=&bkgRatesUsual;
   int valuesType[nValues];   valuesType[0]=1;
   valuesType[1]=2; valuesType[2]=0;
   TMatrixD* valuesErr1[nValues];   valuesErr1[0]=&observedYieldsErr;
   valuesErr1[1]=&totalBackgroundError; valuesErr1[2]=0;
   TMatrixD* valuesErr2[nValues];   valuesErr2[0]=0;
   valuesErr2[1]=&totalBackgroundErrorSyst; valuesErr2[2]=0;
   TString valuesName[nValues];   valuesName[0]="observed yield";
   valuesName[1]="total background"; valuesName[2]="background fraction, \\%";
   TString floatFormats[nValues]; floatFormats[0]="$ %7.0f \\pm %6.0f $";
   floatFormats[1]="$ %7.1f \\pm %6.1f \\pm %6.1f$"; floatFormats[2]="$ %3.1f $";
   TString baseOfReferenceName="yields-signal-backgrounds-2D";
   TString tableName="Data yields vs total background levels predicted by Monte Carlo for %4.0f-%4.0f GeV mass slice.";
   latexPrintoutTwoColumns2D(nValues, valuesType, values, valuesErr1, valuesErr2, producedBy, valuesName, floatFormats, baseOfReferenceName, tableName);
}

void latexPrintoutBackgroundRates1D(TMatrixD observedYields, TMatrixD observedYieldsErr, 
                                    TMatrixD totalBackground, TMatrixD totalBackgroundError, 
                                    TMatrixD totalBackgroundErrorSyst, TMatrixD bkgRatesUsual, 
                                    TString producedBy)
{
   int nValues=3;
   TMatrixD* values[nValues];   values[0]=&observedYields; 
   values[1]=&totalBackground; values[2]=&bkgRatesUsual;
   int valuesType[nValues];   valuesType[0]=1;
   valuesType[1]=2; valuesType[2]=0;
   TMatrixD* valuesErr1[nValues];   valuesErr1[0]=&observedYieldsErr;
   valuesErr1[1]=&totalBackgroundError; valuesErr1[2]=0;
   TMatrixD* valuesErr2[nValues];   valuesErr2[0]=0;
   valuesErr2[1]=&totalBackgroundErrorSyst; valuesErr2[2]=0;
   TString valuesName[nValues];   valuesName[0]="observed yield";
   valuesName[1]="total background"; valuesName[2]="background fraction, \\%";
   TString floatFormats[nValues]; floatFormats[0]="$ %7.0f \\pm %6.0f $";
   floatFormats[1]="$ %7.1f \\pm %6.1f \\pm %6.1f$"; floatFormats[2]="$ %3.1f $";
   TString baseOfReferenceName="yields-signal-backgrounds-1D";
   TString tableName="Data yields vs total background levels predicted by Monte Carlo.";
   latexPrintoutTwoColumns1D(nValues, valuesType, values, valuesErr1, valuesErr2, producedBy, valuesName, floatFormats, baseOfReferenceName, tableName); 
}

void latexPrintoutCrossSection(TMatrixD signalYields      , TMatrixD signalYieldsStatErr, 
		               TMatrixD unfoldedYields    , TMatrixD unfoldedYieldsStatErr,
		               TMatrixD effCorrectedYields, TMatrixD effCorrectedYieldsStatErr,
		               TMatrixD accCorrectedYields, TMatrixD accCorrectedYieldsStatErr,
		               TMatrixD preFsrYields      , TMatrixD preFsrYieldsStatErr, 
                               TMatrixD relCrossSection,           TMatrixD relCrossSectionStatErr, 
                                                                   TMatrixD relCrossSectionSystErr,
                               TMatrixD relCrossSectionDET,        TMatrixD relCrossSectionStatErrDET, 
                                                                   TMatrixD relCrossSectionSystErrDET,
                               TMatrixD relPostFsrCrossSection,    TMatrixD relPostFsrCrossSectionStatErr, 
                                                                   TMatrixD relPostFsrCrossSectionSystErr,
                               TMatrixD relPostFsrCrossSectionDET, TMatrixD relPostFsrCrossSectionStatErrDET, 
                                                                   TMatrixD relPostFsrCrossSectionSystErrDET,
                               TMatrixD absCrossSection,           TMatrixD absCrossSectionStatErr, 
                                                                   TMatrixD absCrossSectionSystErr,
                               TMatrixD absCrossSectionDET,        TMatrixD absCrossSectionStatErrDET, 
                                                                   TMatrixD absCrossSectionSystErrDET,
                               TMatrixD absPostFsrCrossSection,    TMatrixD absPostFsrCrossSectionStatErr, 
                                                                   TMatrixD absPostFsrCrossSectionSystErr,
                               TMatrixD absPostFsrCrossSectionDET, TMatrixD absPostFsrCrossSectionStatErrDET, 
                                                                   TMatrixD absPostFsrCrossSectionSystErrDET,
                               TString  producedBy)
{
   int nValues=5;
   TMatrixD* values[nValues];       values[0]=&signalYields; 
   values[1]=&unfoldedYields;       values[2]=&effCorrectedYields;
   values[3]=&accCorrectedYields;   values[4]=&preFsrYields;

   int valuesType[nValues];         valuesType[0]=1;
   valuesType[1]=1;                 valuesType[2]=1;
   valuesType[3]=1;                 valuesType[4]=1;

   TMatrixD* valuesErr1[nValues];           
   valuesErr1[0]=&signalYieldsStatErr;
   valuesErr1[1]=&unfoldedYieldsStatErr;
   valuesErr1[2]=&effCorrectedYieldsStatErr;
   valuesErr1[3]=&accCorrectedYieldsStatErr;
   valuesErr1[4]=&preFsrYieldsStatErr;

   TMatrixD* valuesErr2[nValues];   valuesErr2[0]=0;
   valuesErr2[1]=0;                 valuesErr2[2]=0;
   valuesErr2[3]=0;                 valuesErr2[4]=0;   

   TString valuesName[nValues];     valuesName[0]="raw signal";
   valuesName[1]="unfolded";        valuesName[2]="eff corrected";
   valuesName[3]="acc corrected";   valuesName[4]="FSR corrected";

   TString floatFormats[nValues];          floatFormats[0]="$%8.1f \\pm %7.1f$";
   floatFormats[1]="$%8.1f \\pm %7.1f$";   floatFormats[2]="$%8.1f \\pm %7.1f$";
   floatFormats[3]="$%8.1f \\pm %7.1f$";   floatFormats[4]="$%8.1f \\pm %7.1f$";

   TString baseOfReferenceName;
   TString tableName;   

   if (DYTools::study2D==0)
     {
       baseOfReferenceName="yields-with-corrections-1D";
       tableName="The \\DYee candidate yields with successive corrections applied.";
       latexPrintoutOneColumn1D(nValues, valuesType, 
                                values, valuesErr1, 
                                valuesErr2, producedBy, 
                                valuesName, floatFormats, 
                                baseOfReferenceName, tableName);
     }  
   else if (DYTools::study2D==1)
     {
       baseOfReferenceName="yields-with-corrections-2D";
       tableName="The \\DYee candidate yields with successive corrections applied for %4.0f-%4.0f GeV mass slice";
       latexPrintoutOneColumn2D(nValues, valuesType, 
                                values, valuesErr1, 
                                valuesErr2, producedBy, 
                                valuesName, floatFormats, 
                                baseOfReferenceName, tableName);
     } 

   if (DYTools::study2D==1)
     {
       baseOfReferenceName="cross-sections-preFSR-Full-2D";
       tableName="Absolute and normalized to Z peak (60-120\\GeVcc) differential \\DYee pre-FSR cross-sections in Full phase space measured in this analysis for %4.0f-%4.0f GeV mass slice";
     }
   else if (DYTools::study2D==0)
     {
       baseOfReferenceName="cross-sections-preFSR-Full-1D";
       tableName="Absolute and normalized to Z peak (60-120\\GeVcc) differential \\DYee pre-FSR cross-sections in Full phase space measured in this analysis";
     }

   latexPrintoutCrossSectionItself(relCrossSection, relCrossSectionStatErr, 
                                   relCrossSectionSystErr, absCrossSection,           
                                   absCrossSectionStatErr, absCrossSectionSystErr,
                                   baseOfReferenceName, tableName,
                                   producedBy); 

   if (DYTools::study2D==1)
     {
       baseOfReferenceName="cross-sections-preFSR-Det-2D";
       tableName="Absolute and normalized to Z peak (60-120\\GeVcc) differential \\DYee pre-FSR cross-sections measured in Detector phase space this analysis for %4.0f-%4.0f GeV mass slice";
     }
   else if (DYTools::study2D==0)
     {
       baseOfReferenceName="cross-sections-preFSR-Det-1D";
       tableName="Absolute and normalized to Z peak (60-120\\GeVcc) differential \\DYee pre-FSR cross-sections measured in Detector phase space this analysis";
     }

   latexPrintoutCrossSectionItself(relCrossSectionDET, relCrossSectionStatErrDET, 
                                   relCrossSectionSystErrDET, absCrossSectionDET,           
                                   absCrossSectionStatErrDET, absCrossSectionSystErrDET,
                                   baseOfReferenceName, tableName,
                                   producedBy);

   if (DYTools::study2D==1)
     {
       baseOfReferenceName="cross-sections-postFSR-Full-2D";
       tableName="Absolute and normalized to Z peak (60-120\\GeVcc) differential \\DYee post-FSR cross-sections measured in Full phase space in this analysis for %4.0f-%4.0f GeV mass slice";
     }
   else if (DYTools::study2D==0)
     {
       baseOfReferenceName="cross-sections-postFSR-Full-1D";
       tableName="Absolute and normalized to Z peak (60-120\\GeVcc) differential \\DYee post-FSR cross-sections measured in Full phase space in this analysis";
     }

   latexPrintoutCrossSectionItself(relPostFsrCrossSection, relPostFsrCrossSectionStatErr, 
                                   relPostFsrCrossSectionSystErr, absPostFsrCrossSection,           
                                   absPostFsrCrossSectionStatErr, absPostFsrCrossSectionSystErr,
                                   baseOfReferenceName, tableName,
                                   producedBy); 

   if (DYTools::study2D==1)
     {
       baseOfReferenceName="cross-sections-postFSR-Det-2D";
       tableName="Absolute and normalized to Z peak (60-120\\GeVcc) differential \\DYee post-FSR cross-sections measured in Detector phase space in this analysis for %4.0f-%4.0f GeV mass slice";
     }
   else if (DYTools::study2D==0)
     {
       baseOfReferenceName="cross-sections-postFSR-Det-1D";
       tableName="Absolute and normalized to Z peak (60-120\\GeVcc) differential \\DYee post-FSR cross-sections measured in Detector phase space in this analysis";
     }

   latexPrintoutCrossSectionItself(relPostFsrCrossSectionDET, relPostFsrCrossSectionStatErrDET, 
                                   relPostFsrCrossSectionSystErrDET, absPostFsrCrossSectionDET,           
                                   absPostFsrCrossSectionStatErrDET, absPostFsrCrossSectionSystErrDET,
                                   baseOfReferenceName, tableName,
                                   producedBy); 

}

void latexPrintoutCrossSectionItself(TMatrixD relCrossSection,           
                                     TMatrixD relCrossSectionStatErr,
                                     TMatrixD relCrossSectionSystErr,
                                     TMatrixD absCrossSection,           
                                     TMatrixD absCrossSectionStatErr,
                                     TMatrixD absCrossSectionSystErr,
                                     TString  baseOfReferenceName,
                                     TString  tableName,
                                     TString  producedBy)
{
   int nValues=2;

   TMatrixD* values[nValues];       
   values[0]=&absCrossSection; values[1]=&relCrossSection;    

   int valuesType[nValues];         
   valuesType[0]=2; valuesType[1]=2;    

   TMatrixD* valuesErr1[nValues];           
   valuesErr1[0]=&absCrossSectionStatErr;
   valuesErr1[1]=&relCrossSectionStatErr;   

   TMatrixD* valuesErr2[nValues];           
   valuesErr2[0]=&absCrossSectionSystErr;
   valuesErr2[1]=&relCrossSectionSystErr; 

   TString valuesName[nValues];     
   valuesName[0]="absolute CS, pb";    valuesName[1]="normalized to Z peak";  

   TString floatFormats[nValues];          
   floatFormats[0]="$%4.4f \\pm %4.4f \\pm %4.4f$";    floatFormats[1]="$%2.8f \\pm %2.8f \\pm %2.8f$";

   if (DYTools::study2D==0)
     latexPrintoutOneColumn1D(nValues, valuesType, 
                                values, valuesErr1, 
                                valuesErr2, producedBy, 
                                valuesName, floatFormats, 
                                baseOfReferenceName, tableName);
   else if (DYTools::study2D==1)
     latexPrintoutOneColumn2D(nValues, valuesType, 
                                values, valuesErr1, 
                                valuesErr2, producedBy, 
                                valuesName, floatFormats, 
                                baseOfReferenceName, tableName);
}


void latexPrintoutOneValue2D(TMatrixD value, TMatrixD valueErr, TString producedBy, TString valueName, TString baseOfReferenceName, TString tableName)
{

   int nValues=1;
   TMatrixD* values[1];   values[0]=&value;
   int valuesType[1];   valuesType[0]=1;
   TMatrixD* valuesErr1[1];   valuesErr1[0]=&valueErr;
   TMatrixD* valuesErr2[1];   valuesErr2[0]=0;
   TString valuesName[1];   valuesName[0]=valueName;
   TString floatFormats[1]; floatFormats[0]=" $%7.4f \\pm %6.4f$ ";
   latexPrintoutTwoColumns2D(nValues, valuesType, values, valuesErr1, valuesErr2, producedBy, valuesName, floatFormats, baseOfReferenceName, tableName);

}


void latexPrintoutOneValue1D(TMatrixD value, TMatrixD valueErr, TString producedBy, TString valueName, TString baseOfReferenceName, TString tableName)
{
 
   int nValues=1;
   TMatrixD* values[1];   values[0]=&value;
   int valuesType[1];   valuesType[0]=1;
   TMatrixD* valuesErr1[1];   valuesErr1[0]=&valueErr;
   TMatrixD* valuesErr2[1];   valuesErr2[0]=0;
   TString valuesName[1];   valuesName[0]=valueName;
   TString floatFormats[1]; floatFormats[0]=" $%7.4f \\pm %6.4f$ ";
   latexPrintoutTwoColumns1D(nValues, valuesType, values, valuesErr1, valuesErr2, producedBy, valuesName,  floatFormats, baseOfReferenceName, tableName);
}

void latexPrintoutTwoColumns2D(const int nValues, int* valuesType, TMatrixD** values, TMatrixD** valuesErr1, TMatrixD** valuesErr2, TString producedBy, TString* valuesName, TString* floatFormats, TString baseOfReferenceName, TString tableName)
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
   TString path="tables2D/";
   gSystem->mkdir(path,kTRUE);
   TString txtFileName=path + "tables2D"+valueNameForSaving+".txt";
   FILE* file;
   if ((file = fopen(txtFileName, "r")))
     {
       for (int i=1; 1; i++)
         {
           TString probeFileName=path + "tables1D"+valueNameForSaving+"-";
           probeFileName+=i;
           probeFileName+=".txt";
           if ((file = fopen(probeFileName, "r")));
           else
             {
               txtFileName=probeFileName;
               break;
             }
         }
     }
   txtFile = fopen(txtFileName,"w");
   TString str=DayAndTimeTag();
   fprintf(txtFile,"%50s",str.Data());
   TString strComment="% Tables produced by "+ producedBy + " in 2D measurement\n";
   fprintf(txtFile,strComment);  
   fprintf(txtFile,"\n\n\n");
   for (int mslice=1; mslice<DYTools::nMassBins; mslice++) 
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
       strQuantities += "& $|Y|$  ";
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
                fprintf(txtFile," &");
                if (valuesType[i]==0) 
                  fprintf(txtFile,floatFormats[i], temp0(mslice,j));
                else if (valuesType[i]==1) 
                  fprintf(txtFile,floatFormats[i], temp0(mslice,j), temp1(mslice,j));
                else if (valuesType[i]==2) 
                  fprintf(txtFile,floatFormats[i], temp0(mslice,j), temp1(mslice,j), temp2(mslice,j));
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
                   fprintf(txtFile," &");
                   if (valuesType[i]==0) 
                     fprintf(txtFile,floatFormats[i], temp0(mslice,j));
                   else if (valuesType[i]==1) 
                     fprintf(txtFile,floatFormats[i], temp0(mslice,j), temp1(mslice,j));
                   else if (valuesType[i]==2) 
                     fprintf(txtFile,floatFormats[i], temp0(mslice,j), temp1(mslice,j), temp2(mslice,j));
                 }

               fprintf(txtFile," \\\\ \n ");
             }
           else fprintf(txtFile," \\\\ \n ");
        }
        fprintf(txtFile,"\\hline\n");
        fprintf(txtFile,"\\end{tabular}}\n");
        fprintf(txtFile,"\\end{center}\n");
        fprintf(txtFile,"\\end{table}\n");
        fprintf(txtFile,"\n\n");
     }
   fclose(txtFile);
}

void latexPrintoutTwoColumns1D(const int nValues, int* valuesType, TMatrixD** values, TMatrixD** valuesErr1, TMatrixD** valuesErr2, TString producedBy, TString* valuesName, TString* floatFormats, TString baseOfReferenceName, TString tableName)
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
   TString path="tables1D/";
   gSystem->mkdir(path,kTRUE);
   TString txtFileName=path + "tables1D"+valueNameForSaving+".txt";
   FILE *file;
   if ((file = fopen(txtFileName, "r")))
     {
       for (int i=1; 1; i++)
         {
           TString probeFileName=path + "tables1D"+valueNameForSaving+"-";
           probeFileName+=i;
           probeFileName+=".txt";
           if ((file = fopen(probeFileName, "r")));
           else
             {
               txtFileName=probeFileName;
               break;
             }
         }
     }
   txtFile = fopen(txtFileName,"w");
   TString str=DayAndTimeTag();
   fprintf(txtFile,"%50s",str.Data());
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
            fprintf(txtFile," &");
            if (valuesType[i]==0) 
              fprintf(txtFile,floatFormats[i], temp0(mslice,0));
            else if (valuesType[i]==1) 
              fprintf(txtFile,floatFormats[i], temp0(mslice,0), temp1(mslice,0));
            else if (valuesType[i]==2) 
              fprintf(txtFile,floatFormats[i], temp0(mslice,0), temp1(mslice,0), temp2(mslice,0));
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
               fprintf(txtFile," &");
               if (valuesType[i]==0) 
                 fprintf(txtFile,floatFormats[i], temp0(mslice+halfBins,0));
               else if (valuesType[i]==1) 
                 fprintf(txtFile,floatFormats[i], temp0(mslice+halfBins,0), temp1(mslice+halfBins,0));
               else if (valuesType[i]==2) 
                 fprintf(txtFile,floatFormats[i], temp0(mslice+halfBins,0), temp1(mslice+halfBins,0), temp2(mslice+halfBins,0));
             }
           fprintf(txtFile," \\\\ \n ");

           //fprintf(txtFile,"$ %7.4f \\pm %6.4f $ \\\\ \n ", value(mslice+halfBins,0), valueErr(mslice+halfBins,0));
         }
       else fprintf(txtFile," \\\\ \n ");
     }
   fprintf(txtFile,"\\hline\n");
   fprintf(txtFile,"\\end{tabular}}\n");
   fprintf(txtFile,"\\end{center}\n");
   fprintf(txtFile,"\\end{table}\n");
   fprintf(txtFile,"\n\n");

   fclose(txtFile);
}

void latexPrintoutOneColumn2D(const int nValues, int* valuesType, TMatrixD** values, TMatrixD** valuesErr1, TMatrixD** valuesErr2, TString producedBy, TString* valuesName, TString* floatFormats, TString baseOfReferenceName, TString tableName)
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
   TString path="tables2D/";
   gSystem->mkdir(path,kTRUE);
   TString txtFileName=path + "tables2D"+valueNameForSaving+".txt";
   FILE *file;
   if ((file = fopen(txtFileName, "r")))
     {
       for (int i=1; 1; i++)
         {
           TString probeFileName=path + "tables2D"+valueNameForSaving+"-";
           probeFileName+=i;
           probeFileName+=".txt";
           if ((file = fopen(probeFileName, "r")));
           else
             {
               txtFileName=probeFileName;
               break;
             }
         }
     }
   txtFile = fopen(txtFileName,"w");
   TString str=DayAndTimeTag();
   fprintf(txtFile,"%50s",str.Data());
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
       strBeginTabular+="}\n";   
       fprintf(txtFile,strBeginTabular);

       fprintf(txtFile,"\\hline\n");

       TString strQuantities = " $|Y|$ ";
       for (int i=0; i<nValues; i++)
         {
           strQuantities+=" & ";
           strQuantities+=valuesName[i];       
         }   
       strQuantities +=  " \\\\ \n";      
       fprintf(txtFile,strQuantities);
       fprintf(txtFile,"\\hline\n");

       for (int j=0; j<DYTools::nYBins[mslice]; j++)
         {
           fprintf(txtFile,"%1.1f-%1.1f", j*(DYTools::yRangeMax-DYTools::yRangeMin)/DYTools::nYBins[mslice], (j+1)*(DYTools::yRangeMax-DYTools::yRangeMin)/DYTools::nYBins[mslice]);
           for (int i=0; i<nValues; i++)
             {
                TMatrixD& temp0= *values[i];
                TMatrixD& temp1= *valuesErr1[i];
                TMatrixD& temp2= *valuesErr2[i];
                fprintf(txtFile," &");
                if (valuesType[i]==0) 
                  fprintf(txtFile,floatFormats[i], temp0(mslice,j));
                else if (valuesType[i]==1) 
                  fprintf(txtFile,floatFormats[i], temp0(mslice,j), temp1(mslice,j));
                else if (valuesType[i]==2) 
                  fprintf(txtFile,floatFormats[i], temp0(mslice,j), temp1(mslice,j), temp2(mslice,j));
             }
          fprintf(txtFile," \\\\ \n ");

         }
       fprintf(txtFile,"\\hline\n");
       fprintf(txtFile,"\\end{tabular}}\n");
       fprintf(txtFile,"\\end{center}\n");
       fprintf(txtFile,"\\end{table}\n");
       fprintf(txtFile,"\n\n");
     }
   fclose(txtFile);
}

void latexPrintoutOneColumn1D(const int nValues, int* valuesType, TMatrixD** values, TMatrixD** valuesErr1, TMatrixD** valuesErr2, TString producedBy, TString* valuesName, TString* floatFormats, TString baseOfReferenceName, TString tableName)
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
   TString path="tables1D/";
   gSystem->mkdir(path,kTRUE);
   TString txtFileName=path + "tables1D"+valueNameForSaving+".txt";
   FILE * file;
   if ((file = fopen(txtFileName, "r")))
     {
       for (int i=1; 1; i++)
         {
           TString probeFileName=path + "tables1D"+valueNameForSaving+"-";
           probeFileName+=i;
           probeFileName+=".txt";
           if ((file = fopen(probeFileName, "r")));
           else
             {
               txtFileName=probeFileName;
               break;
             }
         }
     }

   txtFile = fopen(txtFileName,"w");
   TString str=DayAndTimeTag();
   fprintf(txtFile,"%50s",str.Data());
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
   strBeginTabular+="}\n";   

   fprintf(txtFile,strBeginTabular);

   fprintf(txtFile,"\\hline\n");
   TString strQuantities = " mass, GeV  ";
   for (int i=0; i<nValues; i++)
     {
       strQuantities+=" & ";
       strQuantities+=valuesName[i];       
     }   
   strQuantities +=  " \\\\ \n";      
   fprintf(txtFile,strQuantities);

   fprintf(txtFile,"\\hline\n");

   for (int mslice=0; mslice<DYTools::nMassBins; mslice++)
     {
       fprintf(txtFile,"%4.0f-%4.0f", DYTools::massBinLimits[mslice], DYTools::massBinLimits[mslice+1] );
       for (int i=0; i<nValues; i++)
         {
            TMatrixD& temp0= *values[i];
            TMatrixD& temp1= *valuesErr1[i];
            TMatrixD& temp2= *valuesErr2[i];
            fprintf(txtFile," &");
            if (valuesType[i]==0) 
              fprintf(txtFile,floatFormats[i], temp0(mslice,0));
            else if (valuesType[i]==1) 
              fprintf(txtFile,floatFormats[i], temp0(mslice,0), temp1(mslice,0));
            else if (valuesType[i]==2) 
              fprintf(txtFile,floatFormats[i], temp0(mslice,0), temp1(mslice,0), temp2(mslice,0));
         }

       fprintf(txtFile," \\\\ \n ");
     }
   fprintf(txtFile,"\\hline\n");
   fprintf(txtFile,"\\end{tabular}}\n");
   fprintf(txtFile,"\\end{center}\n");
   fprintf(txtFile,"\\end{table}\n");
   fprintf(txtFile,"\n\n");

   fclose(txtFile);
}
