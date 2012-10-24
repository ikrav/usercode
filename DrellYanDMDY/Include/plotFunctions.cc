#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TTree.h>                  // class to access ntuples
#include <TCanvas.h>                // class for drawing
#include <TH1F.h>                   // 1D histograms
#include <TH2D.h>                   // 2D histograms
#include <TStyle.h> 
#include <THStack.h>                  
#include <TMatrixD.h>
#include <TString.h>
#include <TMultiGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TBenchmark.h> 
#include <TVector3.h>               // 3D vector class
#include <TArrayD.h>
#include <TVectorD.h>
#include <TLorentzVector.h>         // 4-vector class
#include <TRandom.h>
#include <TDatime.h>                // time stamp


#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings

#include "../Include/CSample.hh"        // helper class for organizing input ntuple files
#include "../Include/DYTools.hh"        // helper class for organizing input ntuple files
#include "../Include/MyTools.hh"        // helper class for organizing input ntuple files
#include "../Include/CPlot.hh"          // helper class for plots
#include "../Include/MitStyleRemix.hh"  // style settings for drawing
#include "../Include/ZeeData.hh"
#include "../Include/ElectronEnergyScale.hh"        // energy scale correction
#include "../Include/DYTools.hh"
#include "../Include/plotFunctions.hh"
#include "../Include/UnfoldingTools.hh"
     
#endif

// -----------------------------------------------------------------------------

void Plot1D(const TMatrixD &matrValues, const TMatrixD &matrErrors, TString name, TString title)
{

//acceptance, bkgRatesPercent, LEGO2,COLZ 
   TCanvas* canv=new TCanvas(name,name);
   canv->SetFillColor(0);
   canv->SetLeftMargin     (0.1);
   canv->SetRightMargin    (0.1);
   
   if (DYTools::study2D==1) 
     {
        std::cout<<"Plot1D is not for 2D studies"<<std::endl;
        return;
     }

   if (DYTools::nMassBins!=matrValues.GetNrows() || DYTools::nMassBins!=matrErrors.GetNrows()) 
     {
        std::cout<<"for 1D study, length of vector is not equal to the number of mass bins"<<std::endl;
        return;
     }
   double yVals[DYTools::nMassBins];
   double yErrs[DYTools::nMassBins];
   double xVals[DYTools::nMassBins];
   double xErrs[DYTools::nMassBins];

   for (int i=0; i<DYTools::nMassBins; i++)
     {
        yVals[i]=matrValues(i,0);
        yErrs[i]=matrErrors(i,0);
        xVals[i]=(DYTools::massBinLimits[i+1]+DYTools::massBinLimits[i])/2;
        xErrs[i]=(DYTools::massBinLimits[i+1]-DYTools::massBinLimits[i])/2;
     }

   TGraphErrors *gr = new TGraphErrors(DYTools::nMassBins,xVals,yVals,xErrs,yErrs);
   gr->SetTitle(title);
   TAxis* xAx=gr->GetXaxis();
   xAx->SetTitle("mass, GeV");
   canv->SetLogx();
   gr->SetMarkerStyle(1);
   gr->Draw("AP");
   SaveCanvas(canv, name, CPlot::sOutDir);

}

void PlotMatrixVariousBinning(const TMatrixD &matr, TString name, 
			      TString drawOption, TFile *histoFile, TString title, 
			      bool SetLogz)
{

//acceptance, bkgRatesPercent, LEGO2,COLZ 
   TCanvas* canv=new TCanvas(name,name);
   canv->SetFillColor(0);
   canv->SetLeftMargin     (0.1);
   canv->SetRightMargin    (0.1);

   TMatrixD matrDraw = AdjustMatrixBinning(matr);
   gStyle->SetPalette(1);

   int nM=matrDraw.GetNrows();
   int nY=matrDraw.GetNcols();

   TString histName = "Hist" + name;
   TH2D* Hist= new TH2D(histName,title, nM, DYTools::massBinLimits, nY, 
			DYTools::yRangeMin, DYTools::yRangeMax);
   Hist->SetDirectory(0);
   for (int i=0; i<nM; i++)
     for (int j=0; j<nY; j++)
       {
         Hist->SetBinContent(i+1,j+1,matrDraw(i,j));
       }
  
   Hist->SetStats(0);

   Hist->SetTitle(title);
   TAxis* xAx=Hist->GetXaxis();
   xAx->SetTitle("mass, GeV");
   //xAx->SetMoreLogLabels(kTRUE);
   TAxis* yAx=Hist->GetYaxis();
   yAx->SetTitle("|Y|");
   canv->SetLogx();
   if (SetLogz) canv->SetLogz(); //SetLogz=1 for calcCrossSection

   Hist->Draw(drawOption);
  //CPlot::sOutDir="plots" + DYTools::analysisTag;
  SaveCanvas(canv, name, CPlot::sOutDir);

  canv->SetPhi(120);
  canv->SetTheta(30);
  SaveCanvas(canv, name + TString("_Rotated"), CPlot::sOutDir);

  if (histoFile) canv->Write();

  std::cout << "DYTools::study2D=" << DYTools::study2D << "\n";
  if (DYTools::study2D==0) {
    const TString str1D="_1D";
    TMatrixD zeroErr=matr; zeroErr=0;
    TH1F* h1D= extractMassDependence(name + str1D, title + str1D,
				     matr, zeroErr, 0, 0,0);
    h1D->SetDirectory(0);
    TCanvas *c1=MakeCanvas(name+str1D, name+str1D, 800,600);
    c1->SetLogx(1);
    h1D->Draw();
    c1->Update();
    SaveCanvas(c1, name + str1D, CPlot::sOutDir);
    if (histoFile) c1->Write();
  }

}

// -----------------------------------------------------------------------------

void PlotMatrixVariousBinning(const TMatrixD &matr, TString name, 
			      TString drawOption, TFile *histoFile, bool SetLogz)
{
  // in case if we want names that we place on the plot and under 
  // what we save the plot to be the same
  PlotMatrixVariousBinning(matr, name, drawOption, histoFile, name, SetLogz);
}

// -----------------------------------------------------------------------------


TMatrixD AdjustMatrixBinning(const TMatrixD &matrUsualBinning)
{
  TMatrixD matrOut(DYTools::nMassBins,minMutualMultiple());
  for (int i=0; i<DYTools::nMassBins; i++)
    {
      int nTheSameCells=minMutualMultiple()/DYTools::nYBins[i];
      for (int j=0; j<DYTools::nYBins[i]; j++)
        for (int l=0; l<nTheSameCells; l++)          
           matrOut(i,j*nTheSameCells+l)=matrUsualBinning(i,j); 

     }
  return matrOut;
}

// -----------------------------------------------------------------------------

Int_t minMutualMultiple()
{
  int mult=DYTools::nYBins[0];
  for (int i=1; i<DYTools::nMassBins; i++)
  {
    mult=minMutualMultipleTwo(mult,DYTools::nYBins[i]);

  }
  return mult;
}

Int_t minMutualMultipleTwo(Int_t n1, Int_t n2)
{
  Int_t subMultiple=1;
  Int_t min; 
  if (n1<n2) min=n1;
  else if (n2<=n1) min=n2;
  for (int i=1; i<=min; i++)
    {
      if (((n1%i)==0) && ((n2%i)==0))
        subMultiple=i;
    }
  return subMultiple*(n1/subMultiple)*(n2/subMultiple);
}

  void RShapePlot (TMatrixD relCrossSection, TMatrixD relCrossSectionStatErr, 
TMatrixD relCrossSectionDET, TMatrixD relCrossSectionStatErrDET, 
TMatrixD relPostFsrCrossSection, TMatrixD relPostFsrCrossSectionStatErr, 
TMatrixD relPostFsrCrossSectionDET, TMatrixD relPostFsrCrossSectionStatErrDET)
{

   double* x;
   double* ex;
   double* y1;
   double* ey1;
   double* y2;
   double* ey2;
   double* y3;
   double* ey3;
   double* y4;
   double* ey4;

   Int_t n;
   if (DYTools::study2D==0)
     {
       x  = new double[DYTools::nMassBins];
       ex = new double[DYTools::nMassBins];
       for (int i=0; i<DYTools::nMassBins; i++)
         {
           x[i]=(DYTools::massBinLimits[i+1]+DYTools::massBinLimits[i])/2;
           ex[i]=(DYTools::massBinLimits[i+1]-DYTools::massBinLimits[i])/2;
         }
       y1=relCrossSection.GetMatrixArray();
       ey1=relCrossSectionStatErr.GetMatrixArray();
       y2=relCrossSectionDET.GetMatrixArray();
       ey2=relCrossSectionStatErrDET.GetMatrixArray();
       y3=relPostFsrCrossSection.GetMatrixArray();
       ey3=relPostFsrCrossSectionStatErr.GetMatrixArray();
       y4=relPostFsrCrossSectionDET.GetMatrixArray();
       ey4=relPostFsrCrossSectionStatErrDET.GetMatrixArray();

       n = DYTools::nMassBins;
       RShapeDrawAndSave(n,x,ex,y1,ey1,y2,ey2,y3,ey3,y4,ey4,"RShape1D",0,0);
       delete x;
       delete ex;
     }
   else if (DYTools::study2D==1)
     {
       x  = new double[DYTools::nYBinsMax];
       ex = new double[DYTools::nYBinsMax];       
       y1  = new double[DYTools::nYBinsMax];
       ey1 = new double[DYTools::nYBinsMax];
       y2  = new double[DYTools::nYBinsMax];
       ey2 = new double[DYTools::nYBinsMax];
       y3  = new double[DYTools::nYBinsMax];
       ey3 = new double[DYTools::nYBinsMax];
       y4  = new double[DYTools::nYBinsMax];
       ey4 = new double[DYTools::nYBinsMax];

        for (int i=0; i<DYTools::nMassBins; i++)
         {
           n = DYTools::nYBins[i];
           for (int j=0; j<DYTools::nYBins[i]; j++)
             {
               x[j]=DYTools::yRangeMin+
		 (j+0.5)*(DYTools::yRangeMax-DYTools::yRangeMin)/n;
               ex[j]=0.5*(DYTools::yRangeMax-DYTools::yRangeMin)/n;
               y1[j]=relCrossSection(i,j); 
               ey1[j]=relCrossSectionStatErr(i,j);
               y2[j]=relCrossSectionDET(i,j);
               ey2[j]=relCrossSectionStatErrDET(i,j);
               y3[j]=relPostFsrCrossSection(i,j);
               ey3[j]=relPostFsrCrossSectionStatErr(i,j);
               y4[j]=relPostFsrCrossSectionDET(i,j);
               ey4[j]=relPostFsrCrossSectionStatErrDET(i,j);
             }
           
           TString name2D="RShape2D_mass_";
           name2D+=DYTools::massBinLimits[i];
           name2D+="-";
           name2D+=DYTools::massBinLimits[i+1];

           RShapeDrawAndSave(n,x,ex,y1,ey1,y2,ey2,y3,ey3,y4,ey4,name2D,DYTools::massBinLimits[i],DYTools::massBinLimits[i+1]);

         }

       delete x;
       delete ex;
       delete y1;
       delete ey1;
       delete y2;
       delete ey2;
       delete y3;
       delete ey3;
       delete y4;
       delete ey4;

     }
}

void RShapeDrawAndSave(Int_t n, double* x,double* ex,double* y1,double* ey1,double* y2,double* ey2,double* y3,double* ey3,double* y4,double* ey4, TString name, double mass1, double mass2)
{

   TCanvas* canv=new TCanvas(name,name);
   canv->SetGrid();
   //for 1D analysis x-axis is mass, for 2D analysis x-axis is |Y|
   if (DYTools::study2D==0) canv->SetLogx(1);
   canv->SetLogy(1);
   canv->SetFillColor(0);
   // draw a frame to define the range
   TMultiGraph *mg = new TMultiGraph();
   TGraphErrors *gr1 = new TGraphErrors(n,x,y1,ex,ey1);
   gr1->SetName("PreFSR (d#sigma /dM)/ (d#sigma /dM)_{z}");
   gr1->SetTitle("Pre-FSR (d#sigma /dM)/ (d#sigma /dM)_{z}");
   gr1->SetFillColor(0);
   gr1->SetMarkerColor(kBlack);
   gr1->SetMarkerStyle(20);
   gr1->SetMarkerSize(1.0);

   TGraphErrors *gr2 = new TGraphErrors(n,x,y2,ex,ey2);
   gr2->SetName("PreFSR DET (d#sigma /dM)/ (d#sigma /dM)_{z}");
   gr2->SetTitle("Pre-FSR DET (d#sigma /dM)/ (d#sigma /dM)_{z}");
   gr2->SetFillColor(0);
   gr2->SetMarkerColor(kBlue);
   gr2->SetMarkerStyle(20);
   gr2->SetMarkerSize(1.0);
   gr2->SetLineColor(kBlue);

   TGraphErrors *gr3 = new TGraphErrors(n,x,y3,ex,ey3);
   gr3->SetName("PostFSR (d#sigma /dM)/ (d#sigma /dM)_{z}");
   gr3->SetTitle("Post-FSR (d#sigma /dM)/ (d#sigma /dM)_{z}");
   gr3->SetFillColor(0);
   gr3->SetMarkerColor(kRed);
   gr3->SetMarkerStyle(20);
   gr3->SetMarkerSize(1.0);
   gr3->SetLineColor(kRed);

   TGraphErrors *gr4 = new TGraphErrors(n,x,y4,ex,ey4);
   gr4->SetName("PostFSR DET(d#sigma /dM)/ (d#sigma /dM)_{z}");
   gr4->SetTitle("Post-FSR DEt(d#sigma /dM)/ (d#sigma /dM)_{z}");
   gr4->SetFillColor(0);
   gr4->SetMarkerColor(kGreen);
   gr4->SetMarkerStyle(20);
   gr4->SetMarkerSize(1.0);
   gr4->SetLineColor(kGreen);


   TString massRange;
   if (DYTools::study2D==0) massRange="";
   else if (DYTools::study2D==1) 
     {
       massRange="Mee: ";
       massRange+=int(mass1);
       massRange+=" - ";
       massRange+=int(mass2);
       massRange+=" GeV";
     }


   mg->SetTitle(massRange);

   mg->Add(gr4);
   mg->Add(gr3);
   mg->Add(gr2);
   mg->Add(gr1);

   mg->Draw("ap");


   TAxis* yax=mg->GetYaxis();

   if (DYTools::study2D==0) yax->SetRangeUser(5e-6,2);
   else if (DYTools::study2D==1) yax->SetRangeUser(3e-5,10);
   if (DYTools::study2D==0) mg->GetXaxis()->SetTitle("M_{ee}");
   else if (DYTools::study2D==1) mg->GetXaxis()->SetTitle("|Y|");
   
   if (DYTools::study2D==0) mg->GetYaxis()->SetTitle("1/#sigma_{z} d#sigma /dM");
   else if (DYTools::study2D==1) mg->GetYaxis()->SetTitle("1/#sigma_{z} d#sigma /dMd|Y|");
   mg->GetYaxis()->SetTitleOffset(1.20);

   TLegend *leg;
   if (DYTools::study2D==0 || 
       mass2>DYTools::massBinLimits[DYTools::nMassBins-2]) 
      leg = new TLegend(.50,.55,.95,.95);
   else leg = new TLegend(.60,.75,.95,.95);
   leg->AddEntry(gr1,"Pre FSR All Phase Space");
   leg->AddEntry(gr2,"Pre FSR Detector Phase space");
   leg->AddEntry(gr3,"Post FSR All Phase Space");
   leg->AddEntry(gr4,"Post FSR Detector Phase space");
   leg->SetFillColor(kWhite);
   leg->Draw();
   
   //CPlot::sOutDir="plots" + DYTools::analysisTag;
   SaveCanvas(canv, name, CPlot::sOutDir);


}

//for prepareYields
void DrawMassPeak(vector<TH1F*> hMassv, vector<CSample*> samplev, vector<TString> snamev, TH1F* hMassDibosons, bool hasData, 
		  bool mergeDibosons, TString labelDibosons, Int_t colorDibosons, Double_t lumi, char* lumitext, bool actualBinning, TFile *histoFile)
//for "actual binning": hMassv -> hMassBinsv, hMassDibosons -> hMassBinsDibosons
{

  TString canvName="";
  if (actualBinning) canvName="massBinning";
  else canvName="massPeak";
  TCanvas *c1 = MakeCanvas(canvName,canvName,800,600);
  CPlot plotMass("mass","","m(e^{+}e^{-}) [GeV/c^{2}]", "Events");
  plotMass.SetLogx();
  if(hasData) { plotMass.AddHist1D(hMassv[0],samplev[0]->label,"E"); }
  // Do not draw separately dibosons, but draw the merged histogram if needed
  if(mergeDibosons)
    plotMass.AddToStack(hMassDibosons, labelDibosons, colorDibosons);
  for(UInt_t isam=(hasData) ? 1:0; isam<samplev.size(); isam++){
    if( !(mergeDibosons && (snamev[isam]=="ww" || snamev[isam]=="wz" || snamev[isam]=="zz")))
      plotMass.AddToStack(hMassv[isam],samplev[isam]->label,samplev[isam]->color);
  }
  plotMass.SetLegend(0.75,0.55,0.98,0.9);
  if(lumi>0) plotMass.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
  if(hasData){
    hMassv[0]->GetXaxis()->SetMoreLogLabels();
    hMassv[0]->GetXaxis()->SetNoExponent();
  }
  SetSomeHistAttributes(hMassDibosons,"ewk");
  for(UInt_t isam=0; isam<samplev.size(); isam++) {
    SetSomeHistAttributes(hMassv[isam],snamev[isam]);
  }
  plotMass.SetLogy();
  plotMass.SetLogx();
  plotMass.Draw(c1);
  if (actualBinning)
    {
      plotMass.SetYRange(1.0,2000000);  
      plotMass.SetXRange(20,1500);
      double plotMassMin=(DYTools::study2D) ? 20. : 15.;
      plotMass.SetXRange(plotMassMin,1500);
    }
  else
    {
      plotMass.SetYRange(1.0,300000);
      double plotMassMin=(DYTools::study2D) ? 20. : 15.;
      plotMass.SetXRange(plotMassMin,1500);
    }  
  plotMass.Draw(c1,kFALSE);
  SaveCanvas(c1, canvName);
  if (histoFile) c1->Write();

  //THStack instead of plotMass

  hMassv[0]->SetLineColor(samplev[0]->color);
  SetSomeHistAttributes(hMassDibosons,"ewk");
  if (actualBinning)
    {
      hMassDibosons->GetXaxis()->SetRangeUser((DYTools::study2D) ? 20. : 15.,1500);
      hMassDibosons->GetYaxis()->SetRangeUser(1.0,2000000);
    }
  else
    {
      hMassDibosons->GetXaxis()->SetRangeUser((DYTools::study2D) ? 20. : 15.,1500);
      hMassDibosons->GetYaxis()->SetRangeUser(1.0,300000);
    }  
  THStack* mcHists=new THStack("mcHists","MC hists");
  if(mergeDibosons) mcHists->Add(hMassDibosons);
  for(UInt_t isam=(hasData) ? 1:0; isam<samplev.size(); isam++){
    SetSomeHistAttributes(hMassv[isam],snamev[isam]);
    if (actualBinning)
    {
      hMassv[isam]->GetXaxis()->SetRangeUser((DYTools::study2D) ? 20. : 15.,1500);
      hMassv[isam]->GetYaxis()->SetRangeUser(1.0,2000000);
    }
    else
    {
      hMassv[isam]->GetXaxis()->SetRangeUser((DYTools::study2D) ? 20. : 15.,1500);
      hMassv[isam]->GetYaxis()->SetRangeUser(1.0,300000);
    } 
    if( !(mergeDibosons && (snamev[isam]=="ww" || snamev[isam]=="wz" || snamev[isam]=="zz")))
      mcHists->Add(hMassv[isam]);
  }

  TH1F* ratioHist=(TH1F*)hMassv[0]->Clone("ratioHist");
  for(UInt_t isam=(hasData) ? 1:0; isam<samplev.size(); isam++)
    ratioHist->Add(hMassv[isam],-1);
  ratioHist->Divide(ratioHist,hMassv[0],1.0,1.0);



  if (actualBinning) canvName="massBinning-2";
  else canvName="massPeak-2";
  TCanvas *c2 = new TCanvas(canvName,canvName,800,800); 
  c2->Divide(1,2,0,0);
  TPad* pad1 = (TPad*)c2->GetPad(1);
  TPad* pad2 = (TPad*)c2->GetPad(2);
  pad1->SetPad(0,0.3,1.0,1.0);
  pad2->SetPad(0,0,  1.0,0.28);
  pad1->SetLeftMargin(0.18);
  pad1->SetTopMargin(0.08);
  pad1->SetRightMargin(0.07);
  pad1->SetBottomMargin(0.01); // All X axis labels and titles are thus cut off
  pad2->SetLeftMargin(0.18);
  pad2->SetTopMargin(0.01);
  pad2->SetRightMargin(0.07);
  pad2->SetBottomMargin(0.45);

  pad1->SetLogx();
  pad1->SetLogy();
  pad1->cd();

  mcHists->Draw("hist");
  hMassv[0]->Draw("pe,same");

  pad2->SetLogx();
  pad2->cd();
  ratioHist->Draw();

  SaveCanvas(c2, canvName);
  if (histoFile) c2->Write();

}

// -----------------------------------------------------------

//for prepareYields
void DrawFlattened(vector<TMatrixD*> yields, vector<TMatrixD*> yieldsSumw2, vector<CSample*> samplev, vector<TString> snamev, bool hasData, 
                   bool mergeDibosons, TString labelDibosons, Int_t colorDibosons, Double_t lumi, char* lumitext,
		   TFile *histoFile)
{ 
 //
  // Draw flattened distribution
  //
  // Create the histograms from the yields arrays
  int flatIndexMax = DYTools::getTotalNumberOfBins();
  TH1F *hFlattened[samplev.size()];
  for(UInt_t isam=0; isam < samplev.size(); isam++){
    //sprintf(hname, "hFlattanedMass_%i", isam);
    TString hname="flat-";
    hname+=isam;
    hFlattened[isam] = new TH1F(hname,"",flatIndexMax, 0.0, 1.0*flatIndexMax);
    TMatrixD *thisSampleYields = yields.at(isam);
    TMatrixD *thisSampleYieldsSumw2 = yieldsSumw2.at(isam);
    for(int im = 0; im < DYTools::nMassBins; im++){
      for(int iy = 0; iy < DYTools::nYBins[im]; iy++){
	int iflat = DYTools::findIndexFlat(im, iy);
	hFlattened[isam]->SetBinContent(iflat, (*thisSampleYields)(im,iy) );
	hFlattened[isam]->SetBinError(iflat, sqrt((*thisSampleYieldsSumw2)(im,iy)) );
      }
    }
  }
  // Merge dibosons
  int clone_idx=(hasData) ? 1:0;
  TH1F *hFlattenedDibosons = (TH1F*)hFlattened[clone_idx]->Clone("hFlattenedDibosons");
  hFlattenedDibosons->Reset();
  for(UInt_t isam=0; isam<samplev.size(); isam++) {
    if( snamev[isam] == "ww" || snamev[isam] == "wz" || snamev[isam] == "zz"){
      hFlattenedDibosons->Add(hFlattened[isam]);
    }
  }
  
  // Draw the flattened figure.
  TCanvas *c3 = MakeCanvas("flattened","flattened",800,600);
  CPlot plotFlattened("mass","","m(e^{+}e^{-}) [GeV/c^{2}]","Events");
  if(hasData) { plotFlattened.AddHist1D(hFlattened[0],samplev[0]->label,"E"); }
  // Do not draw separately dibosons, but draw the merged histogram if needed
  if(mergeDibosons)
    plotFlattened.AddToStack(hFlattenedDibosons, labelDibosons, colorDibosons);
  for(UInt_t isam=clone_idx; isam<samplev.size(); isam++){
    if( !(mergeDibosons && (snamev[isam]=="ww" || snamev[isam]=="wz" || snamev[isam]=="zz")))
      plotFlattened.AddToStack(hFlattened[isam],samplev[isam]->label,samplev[isam]->color);
  }
  SetSomeHistAttributes(hFlattenedDibosons,"ewk");
  for(UInt_t isam=0; isam<samplev.size(); isam++) {
    SetSomeHistAttributes(hFlattened[isam],snamev[isam]);
  }
  plotFlattened.SetLegend(0.75,0.55,0.98,0.9);
  if(lumi>0) plotFlattened.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
  if(hasData){
    hFlattened[0]->GetXaxis()->SetMoreLogLabels();
    hFlattened[0]->GetXaxis()->SetNoExponent();
  }
  plotFlattened.SetLogy();
  plotFlattened.Draw(c3);
  //plotFlattened.SetYRange(1.0,200000);  
  plotFlattened.SetYMin(1.0);
  plotFlattened.Draw(c3,kFALSE);
  SaveCanvas(c3, "flattened");
  if (histoFile) c3->Write();
}

// -----------------------------------------------------------

//for prepareYields
void Draw6Canvases(vector<TMatrixD*> yields, vector<TMatrixD*> yieldsSumw2,
                    vector<CSample*> samplev, vector<TString> snamev, 
		   bool hasData, double dataOverMc, double* dataOverMcEachBin, 
		   bool normEachBin, int singleCanvas,
		   TFile *histoFile) 
{

  int jStart;
  if (hasData) jStart=0; 
  else jStart=1;
  TH1F* allHists[DYTools::nMassBins][samplev.size()];

  TH1F *totalMc[DYTools::nMassBins];
  TH1F *ewkHist[DYTools::nMassBins];
  TH1F *ratioHist[DYTools::nMassBins];
  TH1F *ratioClone[DYTools::nMassBins];
  THStack *mcHists[DYTools::nMassBins]; 
  TString normStr =(normEachBin) ? "normEachBin-" : "normZpeak-";
  normStr.Append((singleCanvas!=0) ? "sngl-" : "multi-");
  if (singleCanvas==-1) normStr.Append("rot-");
  std::cout << "normStr=" << normStr << "\n";

  for (int i=1; i<DYTools::nMassBins; i++)
  // loop over mass bins starting from 1st (excluding underflow bin)
    {
      TString stackName="mcStack-" + normStr;
      stackName+=i;
      mcHists[i]=new THStack(stackName,"");
      TString totalMcName="totalMC-" + normStr;
      totalMcName+=i;
      totalMc[i]=new TH1F(totalMcName,"",DYTools::nYBins[i],DYTools::yRangeMin,DYTools::yRangeMax); 
      TString ewkName="ewk-" + normStr;
      ewkName+=i;
      ewkHist[i]=new TH1F(ewkName,"",DYTools::nYBins[i],DYTools::yRangeMin,DYTools::yRangeMax);


      // add entries taking into accout that the order will be reversed
      for (UInt_t j=samplev.size()-1; (j<samplev.size()); j--)
      //for (UInt_t j=0; j<samplev.size(); j++)
      // loop over data, signal MC and background MC samples
        {
          TString histName="hist-" + normStr;
          histName+=i; histName+="-"; histName+=j; 
          allHists[i][j]=new TH1F(histName,"",DYTools::nYBins[i],DYTools::yRangeMin,DYTools::yRangeMax);
          //nYBins, YRangeMin, yRangeMax - from DYTools.hh

          TMatrixD *thisSampleYields = yields.at(j);
          TMatrixD *thisSampleYieldsSumw2 = yieldsSumw2.at(j);

          for (int k=0; k<DYTools::nYBins[i]; k++)
          // loop over allHists bins
            {
              allHists[i][j]->SetBinContent(k+1,(*thisSampleYields)(i,k));
              allHists[i][j]->SetBinError(k+1,sqrt((*thisSampleYieldsSumw2)(i,k)));        
            }

        SetSomeHistAttributes(allHists[i][j],snamev[j]);         

          if (j!=0) 
            {
              if (normEachBin)  allHists[i][j]->Scale(dataOverMcEachBin[i]);
              else allHists[i][j]->Scale(dataOverMc);             
              totalMc[i]->Add(allHists[i][j]);
              if (snamev[j]=="wjets" || snamev[j]=="ww" || snamev[j]=="wz" || snamev[j]=="zz")
                 ewkHist[i]->Add(allHists[i][j]);
              else if (snamev[j]!="zee")
                 mcHists[i]->Add(allHists[i][j]);
            }
        }
      TString ratioHistName="ratioHist-" + normStr;
      ratioHistName+=i;
      ratioHist[i]=(TH1F*)allHists[i][0]->Clone(ratioHistName);
      ratioHist[i]->Add(totalMc[i],-1);
      ratioHist[i]->Divide(ratioHist[i],allHists[i][0],1.0,1.0);

      SetSomeHistAttributes(ewkHist[i],"ewk");

      mcHists[i]->Add(ewkHist[i]);
      for (UInt_t j=1; j<samplev.size(); j++)
        {
          if (snamev[j]=="zee") mcHists[i]->Add(allHists[i][j]);
        }
    }

  TLegend *legend = new TLegend(0.70,0.38, 0.90,0.78);
  legend->SetTextFont(42);
  legend->SetTextSize(0.06);
  legend->SetLineColor(0);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);
  if (hasData)
      legend->AddEntry(allHists[1][0],"data","LP");

  // Add signal MC first
  for (UInt_t j=1; j<samplev.size(); j++)
    {
      if (snamev[j]=="zee") 
	legend->AddEntry(allHists[1][j],"#gamma*/Z#rightarrow ee","PF");
    }

  legend->AddEntry(ewkHist[1],"EWK","PF");
  for (UInt_t j=1; j<samplev.size(); j++)
    {
      if (snamev[j]=="zee") ;
	// Do nothing, already handled signal MC before
      else if (snamev[j]=="wjets" || snamev[j]=="ww" || snamev[j]=="wz" || snamev[j]=="zz") ;
        // Do nothing, already listed
      else if (snamev[j] == "qcd" )
	legend->AddEntry(allHists[1][j],"QCD","PF");
      else if (snamev[j] == "ttbar" )
	legend->AddEntry(allHists[1][j],"t#bar{t}","PF");
      else if (snamev[j] == "wtop" )
	legend->AddEntry(allHists[1][j],"Wt + W#bar{t}","PF");
      else if (snamev[j] == "ztt" )
	legend->AddEntry(allHists[1][j],"#gamma*/Z#rightarrow#tau#tau","PF");
    }

  TLatex *cmsText = new TLatex();
  cmsText->SetTextFont(42);
  cmsText->SetTextSize(0.055);
  cmsText->SetTextAlign(31);
  cmsText->SetNDC();
  cmsText->SetText(0.93, 0.94, "CMS Preliminary");

  TLatex *lumiText = new TLatex();
  lumiText->SetTextFont(42);
  lumiText->SetTextSize(0.05);
  lumiText->SetTextAlign(33);
  lumiText->SetNDC();
  lumiText->SetText(0.91, 0.90, DYTools::strLumiAtECMS);

  TLatex *massLabels[DYTools::nMassBins];
  for (int i=1; i<DYTools::nMassBins; i++)
    {
      TString massStr="";
      massStr+=DYTools::massBinLimits[i];
      massStr+="<m<";
      massStr+=DYTools::massBinLimits[i+1];
      massStr+=" GeV";
      massLabels[i]=new TLatex();
      massLabels[i]->SetTextFont(42);
      massLabels[i]->SetTextSize(0.05);
      massLabels[i]->SetTextAlign(33);
      massLabels[i]->SetNDC();
      massLabels[i]->SetText(0.91, 0.82, massStr);
    }

  TCanvas* canv[DYTools::nMassBins];
  TString canvName;
  TString canvNames="y-Single-";
  if (normEachBin) canvNames+="norm-each-mass-bin";
  else canvNames+="norm-Z-peak";
  TCanvas* canvSingle=NULL;
  if (singleCanvas!=0) {
    if (singleCanvas==1) {
      canvSingle= MakeCanvas(canvNames,canvNames,1200,1800);
      canvSingle->Divide(2,6,0,0);
    }
    else {
      canvNames.Append("-rot");
      canvSingle= MakeCanvas(canvNames,canvNames,1800,1200);
      canvSingle->Divide(3,4,0,0);
    }
  }
  TPad* pad1[DYTools::nMassBins];
  TPad* pad2[DYTools::nMassBins];

  TLine *lineAtOne = new TLine(15.0, 1.0, 1500, 1.0);
  lineAtOne->SetLineStyle(kDashed);
  lineAtOne->SetLineWidth(1);
  lineAtOne->SetLineColor(kBlue);
    
  for (int i=1; i<DYTools::nMassBins; i++) {
    int idx=(i-1)%6+1;

    switch (singleCanvas) {
      case 1: { // portrait
          double xi=0,xf=0,yi=0,yf=0;
          if ((idx%2)==1) {xi=0.0; xf=0.5;}
          if ((idx%2)==0) {xi=0.5; xf=1.0;}
          if (idx==1 || idx==2) {yi=0.66; yf=1.00;}
          if (idx==3 || idx==4) {yi=0.33; yf=0.66;}
          if (idx==5 || idx==6) {yi=0.00; yf=0.33;}

          pad1[i] = (TPad*)canvSingle->GetPad(2*idx-1);
          pad1[i]->SetPad(xi,yi+0.3*(yf-yi),xf,yf);
          pad2[i] = (TPad*)canvSingle->GetPad(2*idx);
          pad2[i]->SetPad(xi,yi,xf,yi+0.28*(yf-yi));
      }
	break;
      case -1: { // landscape
          double xi=0,xf=0,yi=0,yf=0;
          if (((idx-1)/3)==1) {yi=0.0; yf=0.5;}
          if (((idx-1)/3)==0) {yi=0.5; yf=1.0;}
          if (idx==1 || idx==4) {xi=0.00; xf=0.33;}
          if (idx==2 || idx==5) {xi=0.33; xf=0.66;}
          if (idx==3 || idx==6) {xi=0.66; xf=1.00;}

          pad1[i] = (TPad*)canvSingle->GetPad(2*idx-1);
          pad1[i]->SetPad(xi,yi+0.3*(yf-yi),xf,yf);
          pad2[i] = (TPad*)canvSingle->GetPad(2*idx);
          pad2[i]->SetPad(xi,yi,xf,yi+0.28*(yf-yi));
      }
	break;
    case 0: {
          canvName="y-";
          canvName+=i;
          canvName+="-";
          if (normEachBin) canvName+="norm-each-mass-bin";
          else canvName+="norm-Z-peak";
          canv[i]=MakeCanvas(canvName,canvName,600,600);
          canv[i]->Divide(1,2,0,0);
          pad1[i] = (TPad*)canv[i]->GetPad(1);
          pad1[i]->SetPad(0,0.3,1.0,1.0);
          pad2[i] = (TPad*)canv[i]->GetPad(2);
          pad2[i]->SetPad(0,0,1.0,0.28);
        }
      break;
    default:
      std::cout << "\n\n\tDraw6Canvases: do not know what to do with pads for singleCanvas=" << singleCanvas << "\n\n";
      return;
    }

      pad1[i]->SetLeftMargin(0.18);
      pad1[i]->SetTopMargin(0.08);
      pad1[i]->SetRightMargin(0.07);
      pad1[i]->SetBottomMargin(0.01); // All X axis labels and titles are thus cut off

      pad2[i]->SetLeftMargin(0.18);
      pad2[i]->SetTopMargin(0.01);
      pad2[i]->SetRightMargin(0.07);
      pad2[i]->SetBottomMargin(0.45);

      pad1[i]->cd();
      pad1[i]->SetTicky();
      pad1[i]->SetFrameLineWidth(2);
      // Have to draw stack to have axes defined
      // but first ensure all data is visible
      if (allHists[i][0]->GetMaximum() > mcHists[i]->GetMaximum()) {
	mcHists[i]->SetMaximum(allHists[i][0]->GetMaximum());
      }
      mcHists[i]->Draw("hist");
      mcHists[i]->GetXaxis()->SetMoreLogLabels();
      mcHists[i]->GetXaxis()->SetNoExponent();
      mcHists[i]->GetXaxis()->SetTitle("|Y|");
      mcHists[i]->GetYaxis()->SetTitle("entries per bin");
      mcHists[i]->GetYaxis()->SetTitleOffset(1.2);
      mcHists[i]->GetXaxis()->SetLabelOffset(99);
      mcHists[i]->SetMinimum(1.5);
      mcHists[i]->Draw("hist");
      allHists[i][0]->Draw("pe,same");
      cmsText->Draw();
      lumiText->Draw();
      massLabels[i]->Draw();
      legend->Draw();
      
      pad2[i]->cd();
      pad2[i]->SetFrameLineWidth(2);
    
      // Set attributes
      ratioHist[i]->SetLineColor(kBlack);
      ratioHist[i]->SetLineWidth(1);
      ratioHist[i]->SetMarkerSize(1.0);
    
      double ymin = -0.7;
      double ymax = 0.7;
      if(i == 6){
	ymin = -1.5;
	ymax = 1.5;
      }
      ratioHist[i]->GetYaxis()->SetRangeUser(ymin,ymax);
      ratioHist[i]->GetYaxis()->SetNdivisions(3);
      ratioHist[i]->GetYaxis()->SetLabelSize(0.15);
      ratioHist[i]->GetYaxis()->SetTitle("(data-MC)/data");
      ratioHist[i]->GetYaxis()->SetTitleSize(0.1);
      ratioHist[i]->GetYaxis()->SetTitleOffset(0.7);
      //ratioHist[i]->GetXaxis()->SetLabelOffset(1.2);
      ratioHist[i]->GetXaxis()->SetLabelSize(0.15);
      ratioHist[i]->GetXaxis()->SetTickLength(0.10);
      ratioHist[i]->GetXaxis()->SetTitle(mcHists[i]->GetXaxis()->GetTitle());
      ratioHist[i]->GetXaxis()->SetTitleSize(3 * mcHists[i]->GetXaxis()->GetTitleSize());
      ratioHist[i]->SetMarkerColor(kOrange+7);
      ratioHist[i]->SetMarkerStyle(kFullCircle);
      ratioClone[i] = (TH1F*)ratioHist[i]->Clone("ratioClone");
      ratioClone[i]->SetMarkerColor(kBlack);
      ratioClone[i]->SetMarkerStyle(kOpenCircle);

      ratioHist[i]->Draw("PE");
      ratioClone[i]->Draw("PE,same");
    
      lineAtOne->Draw();
      
      if (!singleCanvas)
        {
	  if (histoFile) canv[i]->Write();
          SaveCanvas(canv[i], canvName);
        }

      if (singleCanvas && ((i==DYTools::nMassBins-1) || (idx==6)))
	{
	  TString name=canvNames;
	  char buf[20];
	  if (idx==6) sprintf(buf,"-%d",i/6);
	  else if (i!=DYTools::nMassBins-1) sprintf(buf,"-%d",i/6+1);
	  else buf[0]='\0';
	  name.Append(buf);
	  canvSingle->SetName(name);
	  canvSingle->SetTitle(name);
	  if (idx!=6) {
	    for (int pad_i=2*idx+1; pad_i<=12; ++pad_i) {
	      std::cout << "clearing pad_i=" << pad_i << std::endl;
	      canvSingle->GetPad(pad_i)->Clear();
	    }
	  }
	  if (histoFile) canvSingle->Write();
	  SaveCanvas(canvSingle, name);
	}

   }
   
  //if (singleCanvas)
  //   {
  //     if (histoFile) canvSingle->Write();
  //     SaveCanvas(canvSingle, canvNames);
  //   }

   HERE("leaving");
}

// -----------------------------------------------------------

//for prepareYields
void SetSomeHistAttributes (TH1F* hist, TString samplename)
{

  if (samplename=="data") return;
  else 
    {
       hist->SetFillStyle(1001);
       hist->SetMarkerStyle(21);
       hist->SetMarkerSize(1.8);
       if (samplename=="zee")
	{
	  hist->SetFillColor(kOrange-2);
          hist->SetMarkerColor(kOrange-2);
          hist->SetLineColor(kOrange+3);
	}        
       else if (samplename=="ttbar")
        {
	  hist->SetFillColor(kRed+2);
	  hist->SetMarkerColor(kRed+2);
          hist->SetLineColor(kRed+4);
	}
       else if (samplename=="wjets" || samplename=="ww" || samplename=="wz" || samplename=="zz" || samplename=="ewk")
	{
          hist->SetFillColor(kOrange+10);
	  hist->SetMarkerColor(kOrange+10);
          hist->SetLineColor(kOrange+3);
        }
      else if (samplename=="wtop") 
        {
	  hist->SetFillColor(46);
	  hist->SetMarkerColor(46);
	  hist->SetLineColor(kOrange+3);
        }
      else if (samplename=="ztt") 
        {
	  hist->SetFillColor(kOrange+7);
	  hist->SetMarkerColor(kOrange+7);
          hist->SetLineColor(kOrange+3);
        }
      else if (samplename=="qcd") 
        {
	  hist->SetFillColor(kViolet+5);
	  hist->SetMarkerColor(kViolet+5);
          hist->SetLineColor(kViolet+3);
        }
    }  
  return;  
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

void PlotMatrixMYSlices(const std::vector<int> &indices, int functionOfRapidity,
			const std::vector<TMatrixD> &matrV,
			const std::vector<TMatrixD> &matrErrV,
			const std::vector<TString> &labelV,
			TString name, TString drawOption, TFile *histoFile, 
			TString title, TString yAxisLabel,
			int ncolX, int ncolY,
			double yAxisMin, double yAxisMax)
{
  if (indices.size()==0) {
    std::cout << "PlotMatrixMYSlices: indices.size=0" << std::endl;
    throw 2;
  }
  if (matrV.size()==0) {
    std::cout << "PlotMatrixMYSlices: matrV.size=0" << std::endl;
    throw 2;
  }
  if (matrErrV.size()==0) {
    std::cout << "PlotMatrixMYSlices: matrErrV.size=0" << std::endl;
    throw 2;
  }
  if ((matrV.size()!=matrErrV.size()) || (matrV.size()!=labelV.size())) {
    std::cout << "PlotMatrixMYSlices: matrV.size=" << matrV.size() 
	      << ", matrErrV.size=" << matrErrV.size() 
	      << ", labelV.size=" << labelV.size()
	      << std::endl;
    throw 2;
  }

  gStyle->SetPalette(1);
 
  if ((ncolX<=0) && (ncolY<=0)) {
    ncolX=int(sqrt(indices.size()));
    ncolY=ncolX;
  }
  else {
    if (ncolX<=0) ncolX=indices.size()/ncolY;
    else if (ncolY<=0) ncolY=indices.size()/ncolX;
    if (ncolX==0) ncolX=1;
    if (ncolY==0) ncolY=1;
  }
  
  { 
    int loop=0;
    while (ncolX * ncolY < int(indices.size())) {
      if (loop%2==0) ncolX++;
      else ncolY++;
      loop++;
    }
  }

  int cX=800, cY=800;
  if ((ncolX==2) && (ncolY==3)) { cX=600; cY=900; }
  else if ((ncolX==3) && (ncolY==2)) { cX=900; cY=600; }
  TCanvas* canv=MakeCanvas(name,name, cX, cY);
  canv->Divide(ncolX,ncolY);
  //canv->SetFillColor(0);
  //canv->SetLeftMargin     (0.1);
  //canv->SetRightMargin    (0.1);

  int nM=matrV[0].GetNrows();
  if (nM!=DYTools::nMassBins) {
    std::cout << "PlotMatrixMYSlices: nM!=DYTools::nMassBins" << std::endl;
    throw 2;
  }
  int nY=matrV[0].GetNcols();
  if (nY!=DYTools::nYBinsMax) {
    std::cout << "PlotMatrixMYSlices: nN!=DYTools::nYBinsMax" << std::endl;
    throw 2;
  }
  for (unsigned int i=0; i<matrV.size(); ++i) {
    if (nM!=matrV[i].GetNrows()) {
      std::cout << "PlotMatrixMYSlices: matrV[" << i << "].GetNrows!=DYTools::nMassBins" << std::endl;
      throw 2;
    }
    if (nY!=matrV[i].GetNcols()) {
      std::cout << "PlotMatrixMYSlices: matrV[" << i << "].GetNcols!=DYTools::nYBinsMax" << std::endl;
      throw 2;
    }
  }
  for (unsigned int i=0; i<matrErrV.size(); ++i) {
    if (nM!=matrErrV[i].GetNrows()) {
      std::cout << "PlotMatrixMYSlices: matrErrV[" << i << "].GetNrows!=DYTools::nMassBins" << std::endl;
      throw 2;
    }
    if (nY!=matrErrV[i].GetNcols()) {
      std::cout << "PlotMatrixMYSlices: matrErrV[" << i << "].GetNcols!=DYTools::nYBinsMax" << std::endl;
      throw 2;
    }
  }

  const int colorCount=4;
  const int colors[colorCount] = { kBlack, kRed, kBlue, kGreen+2 };

  TString histName = "Hist" + name;
  const int hBufLen=100;
  char histTextBuf[hBufLen+1],histBuf[hBufLen+1];
  std::vector<CPlot*> cplots;
  cplots.reserve(indices.size());

  const int printData=0;
  if (!functionOfRapidity) {
    for (unsigned int ii=0; ii<indices.size(); ++ii) {
      int iY=indices[ii];
      double yMin=0, yMax=0;
      DYTools::findAbsYValueRange(0, iY, yMin, yMax);
      double yCenter=0.5*(yMin+yMax);

      CPlot* cp=new CPlot();
      cplots.push_back(cp);
      snprintf(histTextBuf, hBufLen, "|Y|_{center}= %5.3f",yCenter);
      cp->AddTextBox(TString(histTextBuf), 0.68, 0.4, 0.95, 0.46, 0);
      cp->SetXTitle("mass, GeV");
      cp->SetYTitle(yAxisLabel);
      if (yAxisMin<yAxisMax) cp->SetYRange(yAxisMin,yAxisMax);

      for (unsigned int iSet=0; iSet<matrV.size(); ++iSet) {
	if (printData) printf("iSet=%d: %s, %s\n",iSet,labelV[iSet].Data(),histTextBuf);
	TString drOpt=(iSet==0) ? drawOption : "LP";
	snprintf(histBuf, hBufLen,"h%d_%s_y%5.3f",iSet,histName.Data(),yCenter);
	TH1F* Hist= new TH1F(histName,title, DYTools::nMassBins, DYTools::massBinLimits);
	Hist->SetDirectory(0);
	//Hist->GetXaxis()->SetTitle("mass, GeV");
	Hist->SetMarkerColor(colors[iSet%colorCount]);
	for (int i=0; i<nM; i++) {
	  int iRap=DYTools::findAbsYBin(i, yCenter);
	  if (printData) printf(" %6.1f-%6.1f  %8.2lf  %8.2lf\n",DYTools::massBinLimits[i],DYTools::massBinLimits[i+1],matrV[iSet](i,iRap),matrErrV[iSet](i,iRap));
	  Hist->SetBinContent(i+1,matrV[iSet](i,iRap));
	  Hist->SetBinError(i+1,matrErrV[iSet](i,iRap));
	}
	cp->AddHist1D(Hist,labelV[iSet],drOpt,colors[iSet%colorCount]);
	cp->SetLogx();
      }
    }
  }
  else {
    for (unsigned int ii=0; ii<indices.size(); ++ii) {
      int iMass=indices[ii];
      CPlot* cp=new CPlot();
      cplots.push_back(cp);
      snprintf(histTextBuf, hBufLen, "mass= %1.0f - %1.0f",
	       DYTools::massBinLimits[iMass], DYTools::massBinLimits[iMass+1]);
      cp->AddTextBox(TString(histTextBuf), 0.58, 0.4, 0.95, 0.46, 0);
      cp->SetXTitle("|Y|");
      cp->SetYTitle(yAxisLabel);
      if (yAxisMin<yAxisMax) cp->SetYRange(yAxisMin,yAxisMax);

      for (unsigned int iSet=0; iSet<matrV.size(); ++iSet) {
	if (printData) printf("iSet=%d: %s, %s\n",iSet,labelV[iSet].Data(),histTextBuf);
	TString drOpt=(iSet==0) ? drawOption : "LP";
	snprintf(histBuf, hBufLen,"h%d_%s_mass%1.0f_%1.0f",iSet,histName.Data(),
		 DYTools::massBinLimits[iMass], DYTools::massBinLimits[iMass+1]);
	TH1F* Hist= new TH1F(histName,"", DYTools::nYBins[iMass], DYTools::yRangeMin, DYTools::yRangeMax);
	Hist->SetDirectory(0);
	//Hist->GetXaxis()->SetTitle("|Y|");
	Hist->SetMarkerColor(colors[iSet%colorCount]);
	for (int j=0; j<DYTools::nYBins[iMass]; j++) {
	  if (printData) {
	    double yMin,yMax;
	    DYTools::findAbsYValueRange(0, j, yMin,yMax);
	    printf(" %6.1f-%6.1f  %8.2lf  %8.2lf\n",yMin,yMax,matrV[iSet](iMass,j),matrErrV[iSet](iMass,j));
	  }
	  Hist->SetBinContent(j+1,matrV[iSet](iMass,j));
	  Hist->SetBinError(j+1,matrErrV[iSet](iMass,j));
	}
	cp->AddHist1D(Hist,labelV[iSet],drOpt,colors[iSet%colorCount]);
      }
    }
  }
  if (printData) fflush(stdout);

  canv->cd();
  for (unsigned int i=0; i<cplots.size(); ++i) {
    cplots[i]->Draw(canv, false, "png", i+1);
  }
  //CPlot::sOutDir="plots" + DYTools::analysisTag;
  SaveCanvas(canv, name, CPlot::sOutDir);
  if (histoFile) canv->Write();

}

// -----------------------------------------------------------------------------

void PlotMatrixMYSlices(const std::vector<int> &indices, int functionOfRapidity,
			const std::vector<TVectorD> &matrFIV,
			const std::vector<TVectorD> &matrErrFIV,
			const std::vector<TString> &labelV,
			TString name, TString drawOption, TFile *histoFile, 
			TString title, TString yAxisLabel,
			int ncolX, int ncolY,
			double yAxisMin, double yAxisMax) {
  if (matrFIV.size()!=matrErrFIV.size()) {
    std::cout << "PlotMatrixMYSlices(FIV): matrFIV.size=" << matrFIV.size() 
	      << ", matrErrFIV.size=" << matrErrFIV.size() << std::endl;
  }
  std::vector<TMatrixD> matrV, matrErrV;
  matrV.reserve(matrFIV.size()); matrErrV.reserve(matrErrFIV.size());
  int res=1;
  for (unsigned int i=0; res && (i<matrFIV.size()); ++i) {
    TMatrixD temp(DYTools::nMassBins,DYTools::nYBinsMax);
    res=unfolding::deflattenMatrix(matrFIV[i], temp);
    if (res) matrV.push_back(temp);
  }
  for (unsigned int i=0; res && (i<matrErrFIV.size()); ++i) {
    TMatrixD temp(DYTools::nMassBins,DYTools::nYBinsMax);
    res=unfolding::deflattenMatrix(matrErrFIV[i], temp);
    if (res) matrErrV.push_back(temp);
  }
  PlotMatrixMYSlices(indices,functionOfRapidity,matrV,matrErrV,labelV,
		   name, drawOption,histoFile,title, 
		   yAxisLabel, ncolX,ncolY,
		   yAxisMin, yAxisMax);
}

// -----------------------------------------------------------------------------

void PrintMatrixMYSlices(const std::vector<int> &indices, int functionOfRapidity,
			 const std::vector<TMatrixD> &matrV,
			 const std::vector<TMatrixD> &matrErrV,
			 const std::vector<TString> &labelV)
{
  if (indices.size()==0) {
    std::cout << "PrintMatrixMYSlices: indices.size=0" << std::endl;
    throw 2;
  }
  if (matrV.size()==0) {
    std::cout << "PrintMatrixMYSlices: matrV.size=0" << std::endl;
    throw 2;
  }
  if (matrErrV.size()==0) {
    std::cout << "PrintMatrixMYSlices: matrErrV.size=0" << std::endl;
    throw 2;
  }
  if ((matrV.size()!=matrErrV.size()) || (matrV.size()!=labelV.size())) {
    std::cout << "PrintMatrixMYSlices: matrV.size=" << matrV.size() 
	      << ", matrErrV.size=" << matrErrV.size() 
	      << ", labelV.size=" << labelV.size()
	      << std::endl;
    throw 2;
  }

  int nM=matrV[0].GetNrows();
  if (nM!=DYTools::nMassBins) {
    std::cout << "PrintMatrixMYSlices: nM!=DYTools::nMassBins" << std::endl;
    throw 2;
  }
  int nY=matrV[0].GetNcols();
  if (nY!=DYTools::nYBinsMax) {
    std::cout << "PrintMatrixMYSlices: nN!=DYTools::nYBinsMax" << std::endl;
    throw 2;
  }
  for (unsigned int i=0; i<matrV.size(); ++i) {
    if (nM!=matrV[i].GetNrows()) {
      std::cout << "PrintMatrixMYSlices: matrV[" << i << "].GetNrows!=DYTools::nMassBins" << std::endl;
      throw 2;
    }
    if (nY!=matrV[i].GetNcols()) {
      std::cout << "PrintMatrixMYSlices: matrV[" << i << "].GetNcols!=DYTools::nYBinsMax" << std::endl;
      throw 2;
    }
  }
  for (unsigned int i=0; i<matrErrV.size(); ++i) {
    if (nM!=matrErrV[i].GetNrows()) {
      std::cout << "PrintMatrixMYSlices: matrErrV[" << i << "].GetNrows!=DYTools::nMassBins" << std::endl;
      throw 2;
    }
    if (nY!=matrErrV[i].GetNcols()) {
      std::cout << "PrintMatrixMYSlices: matrErrV[" << i << "].GetNcols!=DYTools::nYBinsMax" << std::endl;
      throw 2;
    }
  }

  if (!functionOfRapidity) {
    for (unsigned int ii=0; ii<indices.size(); ++ii) {
      int iY=indices[ii];
      double yMin=0, yMax=0;
      DYTools::findAbsYValueRange(0, iY, yMin, yMax);
      double yCenter=0.5*(yMin+yMax);

      printf("\n\n|Y|_{center}= %5.3f\n",yCenter);
      printf("   mass range   |");
      for (unsigned int iSet=0; iSet<labelV.size(); ++iSet) {
	printf(" %s |",labelV[iSet].Data());
      }
      printf("\n");

      for (int i=0; i<nM; i++) {
	int iRap=DYTools::findAbsYBin(i, yCenter);
	printf(" %6.1f-%6.1f  ",DYTools::massBinLimits[i],DYTools::massBinLimits[i+1]);
	for (unsigned int iSet=0; iSet<matrV.size(); ++iSet) {
	  printf("   %8.2lf +/- %8.2lf",matrV[iSet](i,iRap),matrErrV[iSet](i,iRap));
	}
	printf("\n");
      }
    }
  }
  else {
    for (unsigned int ii=0; ii<indices.size(); ++ii) {
      int iMass=indices[ii];
      printf("\n\nmass= %1.0f - %1.0f\n",
	     DYTools::massBinLimits[iMass], DYTools::massBinLimits[iMass+1]);
      printf("   rapidity rng |");
      for (unsigned int iSet=0; iSet<labelV.size(); ++iSet) {
	printf(" %s |",labelV[iSet].Data());
      }
      printf("\n");

      for (int j=0; j<DYTools::nYBins[iMass]; j++) {
	double yMin,yMax;
	DYTools::findAbsYValueRange(0, j, yMin,yMax);
	printf(" %6.1f-%6.1f  ",yMin,yMax);
	for (unsigned int iSet=0; iSet<matrV.size(); ++iSet) {
	  printf("   %8.2lf +/- %8.2lf",matrV[iSet](iMass,j),matrErrV[iSet](iMass,j));
	}
	printf("\n");
      }
    }
  }
  fflush(stdout);
}

// -----------------------------------------------------------------------------

void PrintMatrixMYSlices(const std::vector<int> &indices, int functionOfRapidity,
			 const std::vector<TVectorD> &matrFIV,
			 const std::vector<TVectorD> &matrErrFIV,
			 const std::vector<TString> &labelV) {
  if (matrFIV.size()!=matrErrFIV.size()) {
    std::cout << "PrintMatrixMYSlices(FIV): matrFIV.size=" << matrFIV.size() 
	      << ", matrErrFIV.size=" << matrErrFIV.size() << std::endl;
  }
  std::vector<TMatrixD> matrV, matrErrV;
  matrV.reserve(matrFIV.size()); matrErrV.reserve(matrErrFIV.size());
  int res=1;
  for (unsigned int i=0; res && (i<matrFIV.size()); ++i) {
    TMatrixD temp(DYTools::nMassBins,DYTools::nYBinsMax);
    res=unfolding::deflattenMatrix(matrFIV[i], temp);
    if (res) matrV.push_back(temp);
  }
  for (unsigned int i=0; res && (i<matrErrFIV.size()); ++i) {
    TMatrixD temp(DYTools::nMassBins,DYTools::nYBinsMax);
    res=unfolding::deflattenMatrix(matrErrFIV[i], temp);
    if (res) matrErrV.push_back(temp);
  }
  PrintMatrixMYSlices(indices,functionOfRapidity,matrV,matrErrV,labelV);
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

//#include "../YieldsAndBackgrounds/plotFunctionsPrepareYields.C"
