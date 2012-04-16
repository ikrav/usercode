#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TCanvas.h>                // class for drawing
#include <TH1F.h>                   // 1D histograms
#include <THStack.h>
#include <TBenchmark.h>             // class to track macro running statistics
#include <TVector3.h>               // 3D vector class
#include <TArrayD.h>
#include <TVectorD.h>
#include <TMatrixD.h>
#include <TLorentzVector.h>         // 4-vector class
#include <TRandom.h>
#include <TDatime.h>                // time stamp

#include <TLatex.h> 
#include <TStyle.h> 

#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings

#include "../Include/CPlot.hh"          // helper class for plots
#include "../Include/MitStyleRemix.hh"  // style settings for drawing
#include "../Include/MyTools.hh"        // miscellaneous helper functions
#include "../Include/CSample.hh"        // helper class for organizing input ntuple files
     
// define structures to read in ntuple
#include "../Include/ZeeData.hh"

#include "../Include/ElectronEnergyScale.hh"        // energy scale correction

#include "../Include/DYTools.hh"
#include "../Include/plotFunctions.hh"

#endif

//for prepareYields
void DrawMassPeak(vector<TH1F*> hMassv, vector<CSample*> samplev, vector<TString> snamev, TH1F* hMassDibosons, bool hasData, 
                   bool mergeDibosons, TString labelDibosons, Int_t colorDibosons, Double_t lumi, char* lumitext, bool actualBinning)
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
  for(UInt_t isam=1; isam<samplev.size(); isam++){
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
    }
  else
    {
      plotMass.SetYRange(1.0,300000);
    }  
  plotMass.Draw(c1,kFALSE);
  SaveCanvas(c1, canvName);
}

//for prepareYields
void DrawFlattened(vector<TMatrixD*> yields, vector<TMatrixD*> yieldsSumw2, vector<CSample*> samplev, vector<TString> snamev, bool hasData, 
                   bool mergeDibosons, TString labelDibosons, Int_t colorDibosons, Double_t lumi, char* lumitext)
{ 
 //
  // Draw flattened distribution
  //
  // Create the histograms from the yields arrays
  int flatIndexMax = DYTools::getNumberOf2DBins();
  TH1F *hFlattened[samplev.size()];
  for(UInt_t isam=0; isam < samplev.size(); isam++){
    //sprintf(hname, "hFlattanedMass_%i", isam);
    TString hname="flat-";
    hname+=isam;
    hFlattened[isam] = new TH1F(hname,"",flatIndexMax, 0.0, 1.0*flatIndexMax);
    TMatrixD *thisSampleYields = yields.at(isam);
    TMatrixD *thisSampleYieldsSumw2 = yieldsSumw2.at(isam);
    for(int im = 0; im < DYTools::nMassBins2D; im++){
      for(int iy = 0; iy < DYTools::nYBins[im]; iy++){
	int iflat = findIndexFlat(im, iy);
	hFlattened[isam]->SetBinContent(iflat, (*thisSampleYields)(im,iy) );
	hFlattened[isam]->SetBinError(iflat, sqrt((*thisSampleYieldsSumw2)(im,iy)) );
      }
    }
  }
  // Merge dibosons
  TH1F *hFlattenedDibosons = (TH1F*)hFlattened[1]->Clone("hFlattenedDibosons");
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
  for(UInt_t isam=1; isam<samplev.size(); isam++){
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
  plotFlattened.SetYRange(1.0,200000);  
  plotFlattened.Draw(c3,kFALSE);
  SaveCanvas(c3, "flattened");
}

//for prepareYields
void Draw6Canvases(vector<TMatrixD*> yields, vector<TMatrixD*> yieldsSumw2,
                    vector<CSample*> samplev, vector<TString> snamev, 
                    bool hasData, double dataOverMc, double* dataOverMcEachBin, bool normEachBin, bool singleCanvas) 
{

  int jStart;
  if (hasData) jStart=0; 
  else jStart=1;
  TH1F* allHists[nMassBins2D][samplev.size()];

  TH1F *totalMc[nMassBins2D];
  TH1F *ewkHist[nMassBins2D];
  TH1F *ratioHist[nMassBins2D];
  TH1F *ratioClone[nMassBins2D];
  THStack *mcHists[nMassBins2D]; 

  for (int i=1; i<nMassBins2D; i++)
  // loop over mass bins starting from 1st (excluding underflow bin)
    {
      TString stackName="mcStack-";
      stackName+=i;
      mcHists[i]=new THStack(stackName,"");
      TString totalMcName="totalMC-";
      totalMcName+=i;
      totalMc[i]=new TH1F(totalMcName,"",nYBins[i],yRangeMin,yRangeMax); 
      TString ewkName="ewk-";
      ewkName+=i;
      ewkHist[i]=new TH1F(ewkName,"",nYBins[i],yRangeMin,yRangeMax);


 
      for (int j=0; j<samplev.size(); j++)
      // loop over data, signal MC and background MC samples
        {

          TString histName="hist-";
          histName+=i; histName+="-"; histName+=j; 
          allHists[i][j]=new TH1F(histName,"",nYBins[i],yRangeMin,yRangeMax);
          //nYBins, YRangeMin, yRangeMax - from DYTools.hh

          TMatrixD *thisSampleYields = yields.at(j);
          TMatrixD *thisSampleYieldsSumw2 = yieldsSumw2.at(j);

          for (int k=0; k<nYBins[i]; k++)
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
      TString ratioHistName="ratioHist-";
      ratioHistName+=i;
      ratioHist[i]=(TH1F*)allHists[i][0]->Clone(ratioHistName);
      ratioHist[i]->Add(totalMc[i],-1);
      ratioHist[i]->Divide(ratioHist[i],allHists[i][0],1.0,1.0);

      SetSomeHistAttributes(ewkHist[i],"ewk");

      mcHists[i]->Add(ewkHist[i]);
      for (int j=1; j<samplev.size(); j++)
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
  for (int j=1; j<samplev.size(); j++)
    {
      if (snamev[j]=="zee") 
	legend->AddEntry(allHists[1][j],"#gamma*/Z#rightarrow ee","PF");
    }

  for (int j=1; j<samplev.size(); j++)
    {
      if (snamev[j]=="zee") ;
	// Do nothing, already handled signal MC before
      else if (snamev[j]=="wjets" || snamev[j]=="ww" || snamev[j]=="wz" || snamev[j]=="zz"); 
      else if (snamev[j] == "qcd" )
	legend->AddEntry(allHists[1][j],"QCD","PF");
      else if (snamev[j] == "ttbar" )
	legend->AddEntry(allHists[1][j],"t#bar{t}","PF");
      else if (snamev[j] == "ztt" )
	legend->AddEntry(allHists[1][j],"#gamma*/Z#rightarrow#tau#tau","PF");
    }

  legend->AddEntry(ewkHist[1],"EWK","PF");

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
  lumiText->SetText(0.91, 0.90, "4.7 fb^{-1} at #sqrt{s} = 7 TeV");

  TLatex *massLabels[nMassBins2D];
  for (int i=1; i<nMassBins2D; i++)
    {
      TString massStr="";
      massStr+=massBinLimits2D[i];
      massStr+="<m<";
      massStr+=massBinLimits2D[i+1];
      massStr+=" GeV";
      massLabels[i]=new TLatex();
      massLabels[i]->SetTextFont(42);
      massLabels[i]->SetTextSize(0.05);
      massLabels[i]->SetTextAlign(33);
      massLabels[i]->SetNDC();
      massLabels[i]->SetText(0.91, 0.82, massStr);
    }

  TCanvas* canv[nMassBins2D];
  TString canvName;
  TString canvNames="y-Single-";
  if (normEachBin) canvNames+="norm-each-mass-bin";
  else canvNames+="norm-Z-peak";
  TCanvas* canvSingle=new TCanvas(canvNames,canvNames,1200,1800);
  canvSingle->Divide(2,6,0,0);
  TPad* pad1[nMassBins2D];
  TPad* pad2[nMassBins2D];

  TLine *lineAtOne = new TLine(15.0, 1.0, 1500, 1.0);
  lineAtOne->SetLineStyle(kDashed);
  lineAtOne->SetLineWidth(1);
  lineAtOne->SetLineColor(kBlue);
    
  for (int i=1; i<nMassBins2D; i++)
    {

      if (singleCanvas)
       {
          double xi,xf,yi,yf;
          if ((i%2)==1) {xi=0.0; xf=0.5;}
          if ((i%2)==0) {xi=0.5; xf=1.0;}
          if (i==1 || i==2) {yi=0.66; yf=1.00;}
          if (i==3 || i==4) {yi=0.33; yf=0.66;}
          if (i==5 || i==6) {yi=0.00; yf=0.33;}

          pad1[i] = (TPad*)canvSingle->GetPad(2*i-1);
          pad1[i]->SetPad(xi,yi+0.3*(yf-yi),xf,yf);
          pad2[i] = (TPad*)canvSingle->GetPad(2*i);
          pad2[i]->SetPad(xi,yi,xf,yi+0.28*(yf-yi));
       }
      else
       {
          canvName="y-";
          canvName+=i;
          canvName+="-";
          if (normEachBin) canvName+="norm-each-mass-bin";
          else canvName+="norm-Z-peak";
          canv[i]=new TCanvas(canvName,canvName,600,600);
          canv[i]->Divide(1,2,0,0);
          pad1[i] = (TPad*)canv[i]->GetPad(1);
          pad1[i]->SetPad(0,0.3,1.0,1.0);
          pad2[i] = (TPad*)canv[i]->GetPad(2);
          pad2[i]->SetPad(0,0,1.0,0.28);
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
          SaveCanvas(canv[i], canvName);
        }

   }
   
   if (singleCanvas)
     {
       SaveCanvas(canvSingle, canvNames);
     }

}

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
