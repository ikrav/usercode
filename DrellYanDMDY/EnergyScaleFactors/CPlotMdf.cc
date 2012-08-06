#include "CPlotMdf.hh"
#include <TLatex.h>

#if !defined(__CINT__) || defined(__MAKECINT__)
  int CPlot::sCount = 0;
  TString CPlot::sOutDir = ".";
#endif


CPlot::CPlot()
{
  TString name = "plot"; 
  name += sCount;
  CPlot(name,"","",""); 
}

CPlot::CPlot(TString name, TString title, TString xtitle, TString ytitle):
fStack(0),
fName(name),
fTitle(title),
fXTitle(xtitle),
fYTitle(ytitle),
fXmin(0),
fXmax(0),
fYmin(0),
fYmax(0),
fLogx(0),
fLogy(0),
fGridx(0),
fGridy(0),
fRebin(1),
fLeg(0),
fShowStats(0),
fStatsX(0.68),
fStatsY(0.90),
fRooPlot(0)
{
  sCount++;
}

CPlot::CPlot(TString name, RooPlot* frame, TString title, TString xtitle, TString ytitle):
fStack(0),
fName(name),
fTitle(title),
fXTitle(xtitle),
fYTitle(ytitle),
fXmin(0),
fXmax(0),
fYmin(0),
fYmax(0),
fLogx(0),
fLogy(0),
fGridx(0),
fGridy(0),
fRebin(1),
fLeg(0),
fShowStats(0),
fStatsX(0.68),
fStatsY(0.90),
fRooPlot(frame) 
{
  fRooPlot->SetTitle(title);
  fRooPlot->GetXaxis()->SetTitle(xtitle);
  fRooPlot->GetYaxis()->SetTitle(ytitle);
  sCount++;
}


//--------------------------------------------------------------------------------------------------
void CPlot::AddHist1D(TH1F *h, TString drawopt, int color, int linesty, int fillsty)
{
  if(!h)
    return;
   
  h->SetLineColor(color);
  h->SetLineStyle(linesty);
  h->SetFillColor(color);
  h->SetFillStyle(fillsty);
  
  //std::cout << "setLineStyle=" << linesty << "\n";
  CPlotItem item;
  item.hist1F = h;
  item.drawopt = drawopt;
  fItems.push_back(item);
  //h->DrawCopy();
}

void CPlot::AddHist1D(TH1F *h, TString label, TString drawopt, int color, int linesty, int fillsty)
{
  if(!h)
    return;

  if(!fLeg)
    fLeg = new TLegend(0.6,0.84,0.93,0.9);
  else
    fLeg->SetY1(fLeg->GetY1()-0.06);
 
  if(drawopt.CompareTo("E",TString::kIgnoreCase)==0) {
    //fLeg->AddEntry(h,label,"P");
    fLeg->AddEntry(h,label,"PL");
  } else {
    if(fillsty>0) fLeg->AddEntry(h,label,"F");
    else          fLeg->AddEntry(h,label,"L");
  } 
  
  AddHist1D(h,drawopt,color,linesty,fillsty);
}

void CPlot::AddHist1D(TH1D *h, TString drawopt, int color, int linesty, int fillsty)
{
  if(!h)
    return;
   
  h->SetLineColor(color);
  h->SetLineStyle(linesty);
  h->SetFillColor(color);
  h->SetFillStyle(fillsty);
  
  //std::cout << "setLineStyle=" << linesty << "\n";
  CPlotItem item;
  item.hist1D = h;
  item.drawopt = drawopt;
  fItems.push_back(item);
  //h->DrawCopy();
}

void CPlot::AddHist1D(TH1D *h, TString label, TString drawopt, int color, int linesty, int fillsty)
{
  if(!h)
    return;

  if(!fLeg)
    fLeg = new TLegend(0.6,0.84,0.93,0.9);
  else
    fLeg->SetY1(fLeg->GetY1()-0.06);
 
  if(drawopt.CompareTo("E",TString::kIgnoreCase)==0) {
    //fLeg->AddEntry(h,label,"P");
    fLeg->AddEntry(h,label,"PL");
  } else {
    if(fillsty>0) fLeg->AddEntry(h,label,"F");
    else          fLeg->AddEntry(h,label,"L");
  } 
  
  AddHist1D(h,drawopt,color,linesty,fillsty);
}

void CPlot::AddHist1D(TFile *f, TString histName, TString drawopt, int color, int linesty, int fillsty)
{
  if(!f)
    return;
  
  TH1F *h = (TH1F*)f->FindObjectAny(histName);
  AddHist1D(h,drawopt,color,linesty,fillsty);
}

void CPlot::AddHist1D(TFile *f, TString histName, TString label, TString drawopt, int color, int linesty, int fillsty)
{
  if(!f)
    return;
  
  TH1F *h = (TH1F*)f->FindObjectAny(histName);
  AddHist1D(h,label,drawopt,color,linesty,fillsty);
}

//--------------------------------------------------------------------------------------------------
void CPlot::AddToStack(TH1F *h, int color)
{
  if(!h)
    return;
    
  if(!fStack)
    fStack = new THStack(fName+TString("_stack"),"");
  
  fStack->Add(h);
  AddHist1D(h,"",color,1,1001);
}

void CPlot::AddToStack(TH1F *h, TString label, int color)
{
  if(!h)
    return;
  
  if(!fStack)
    fStack = new THStack(fName+TString("_stack"),"");
  
  if(!fLeg)
    fLeg = new TLegend(0.6,0.84,0.93,0.9);
  else
    fLeg->SetY1(fLeg->GetY1()-0.06);  
    
  // make legend entries appear in reverse of the order the histograms are added
  fStackEntries.push_back(fLeg->AddEntry(h,label,"F"));
  for(Int_t ientry=(fStackEntries.size()-2); ientry>=0; ientry--) {
    TObject* hh = fStackEntries[ientry]->GetObject();
    TString ll  = fStackEntries[ientry]->GetLabel();
    fStackEntries[ientry+1]->SetObject(hh);
    fStackEntries[ientry+1]->SetLabel(ll);
  }
  fStackEntries[0]->SetObject(h);
  fStackEntries[0]->SetLabel(label);
     
  fStack->Add(h);  
  AddHist1D(h,"",color,1,1001);  
}

void CPlot::AddToStack(TFile *f, TString histName, int color)
{
  if(!f)
    return;
  
  TH1F *h = (TH1F*)f->FindObjectAny(histName);
  AddToStack(h,color);
}

void CPlot::AddToStack(TFile *f, TString histName, TString label, int color)
{
  if(!f)
    return;
  
  TH1F *h = (TH1F*)f->FindObjectAny(histName);
  AddToStack(h,label,color);
}  
  
//--------------------------------------------------------------------------------------------------
void CPlot::AddHist2D(TH2D *h, TString drawopt, int fillcolor, int linecolor)
{
  if(!h)
    return;
  
  h->SetLineColor(linecolor);
  h->SetFillColor(fillcolor);
  h->SetMarkerStyle(kFullDotMedium); 
  
  CPlotItem item;
  item.hist2D = h;
  item.drawopt = drawopt;
  fItems.push_back(item);
}

void CPlot::AddHist2D(TFile *f, TString histName, TString drawopt, int fillcolor, int linecolor)
{
  if(!f)
    return;
  
  TH2D *h = (TH2D*)f->FindObjectAny(histName);
  AddHist2D(h,drawopt,linecolor,fillcolor);
}

//--------------------------------------------------------------------------------------------------
void CPlot::AddGraph(TGraph *gr, TString drawopt, int color, int marksty, int linesty)
{
  if(!gr) {
    std::cout << "CPlot::AddGraph(TGraph,drawopt): null graph\n";
    return;
  }
  
  gr->SetMarkerColor(color);
  gr->SetLineColor(color);
  gr->SetLineStyle(linesty);
  gr->SetLineWidth(2);
  gr->SetMarkerStyle(marksty);
  gr->SetMarkerSize(1.2);
  
  CPlotItem item;
  item.graph = gr;
  item.drawopt = drawopt;  
  fItems.push_back(item);
}

void CPlot::AddGraph(TGraph *gr, TString label, TString drawopt, int color, int marksty, int linesty)
{
  if(!gr) {
    std::cout << "CPlot::AddGraph(TGraph,label,drawopt): null graph\n";
    return;
  }

  if(!fLeg)
    fLeg = new TLegend(0.6,0.84,0.93,0.9);
  else
    fLeg->SetY1(fLeg->GetY1()-0.06);
  if( (drawopt.Contains("L",TString::kIgnoreCase)==0) ||
      (drawopt.Contains("C",TString::kIgnoreCase)==0) ) {    
    fLeg->AddEntry(gr,label,"LP");
  } else {  
    fLeg->AddEntry(gr,label,"P"); 
  }
  
  AddGraph(gr,drawopt,color,marksty,linesty);
}

void CPlot::AddGraph(TFile *f, TString grName, TString drawopt, int color, int marksty, int linesty)
{
  if(!f)
    return;
  
  TGraph *gr = (TGraph*)f->FindObjectAny(grName);
  AddGraph(gr,drawopt,color,marksty,linesty);
}

void CPlot::AddGraph(TFile *f, TString grName, TString label, TString drawopt, int color, int marksty, int linesty)
{
  if(!f)
    return;
  
  TGraph *gr = (TGraph*)f->FindObjectAny(grName);
  AddGraph(gr,label,drawopt,color,marksty,linesty);
}

//--------------------------------------------------------------------------------------------------
void CPlot::AddProfile(TProfile *pr, TString drawopt, int color, int marksty, int linesty)
{
  if(!pr)
    return;
  
  pr->SetMarkerColor(color);
  pr->SetLineColor(color);
  pr->SetLineStyle(linesty);
  pr->SetLineWidth(2);
  pr->SetMarkerStyle(marksty);
  pr->SetMarkerSize(1.2);
  
  CPlotItem item;
  item.prof = pr;
  item.drawopt = drawopt;  
  fItems.push_back(item);
}

void CPlot::AddProfile(TProfile *pr, TString label, TString drawopt, int color, int marksty, int linesty)
{
  if(!pr)
    return;

  if(!fLeg)
    fLeg = new TLegend(0.6,0.84,0.93,0.9);
  else
    fLeg->SetY1(fLeg->GetY1()-0.06);
    
  fLeg->AddEntry(pr,label,"LP");
  
  AddProfile(pr,drawopt,color,marksty,linesty);
}

void CPlot::AddProfile(TFile *f, TString prName, TString drawopt, int color, int marksty, int linesty)
{
  if(!f)
    return;
  
  TProfile *pr = (TProfile*)f->FindObjectAny(prName);
  AddProfile(pr,drawopt,color,marksty,linesty);
}

void CPlot::AddProfile(TFile *f, TString prName, TString label, TString drawopt, int color, int marksty, int linesty)
{
  if(!f)
    return;
  
  TProfile *pr = (TProfile*)f->FindObjectAny(prName);
  AddProfile(pr,label,drawopt,color,marksty,linesty);
}


//--------------------------------------------------------------------------------------------------
void CPlot::AddTextBox(TString text, double x1, double y1, double x2, double y2,
                       int bordersize, int textcolor, int fillcolor)
{
  TPaveText *tb = new TPaveText(x1,y1,x2,y2,"NDC");
  tb->SetTextColor(textcolor);
  if(fillcolor==-1)
    tb->SetFillStyle(0);
  else
    tb->SetFillColor(fillcolor);
  tb->SetBorderSize(bordersize);
  tb->AddText(text);
  fTextBoxes.push_back(tb);
}

void CPlot::AddTextBox(double x1, double y1, double x2, double y2, 
                       int bordersize, int textcolor, int fillcolor, int nlines,...)
{
  TPaveText *tb = new TPaveText(x1,y1,x2,y2,"NDC");
  tb->SetTextColor(textcolor);
  if(fillcolor==-1)
    tb->SetFillStyle(0);
  else
    tb->SetFillColor(fillcolor);
  tb->SetBorderSize(bordersize);
  tb->SetTextAlign(12);
  
  va_list ap;
  va_start(ap,nlines);
  for(int i=0; i<nlines; i++) {
    TString textline(va_arg(ap,char*));
    tb->AddText(textline);
  }
  va_end(ap);
  
  fTextBoxes.push_back(tb);  
}

//--------------------------------------------------------------------------------------------------
void CPlot::AddLine(double x1, double y1, double x2, double y2, int color, int style)
{
  TLine *line = new TLine(x1,y1,x2,y2);
  line->SetLineColor(color);
  line->SetLineStyle(style);
  line->SetLineWidth(2);
  fLines.push_back(line);
}

void CPlot::AddLine(double x1, double y1, double x2, double y2, 
                    int color, int style, TString label)
{
  TLine *line = new TLine(x1,y1,x2,y2);
  line->SetLineColor(color);
  line->SetLineStyle(style);
  line->SetLineWidth(2);
  fLines.push_back(line);

  if(!fLeg)
    fLeg = new TLegend(0.6,0.84,0.93,0.9);
  else
    fLeg->SetY1(fLeg->GetY1()-0.06);
  fLeg->AddEntry(line,label,"L");
}

//--------------------------------------------------------------------------------------------------
void CPlot::AddBox(double x1, double y1, double x2, double y2, 
                   int linecolor, int linesty, int fillcolor)
{
  TBox *box = new TBox(x1,y1,x2,y2);
  box->SetLineColor(linecolor);
  box->SetLineStyle(linesty);
  if(fillcolor==-1)
    box->SetFillStyle(0);
  else
    box->SetFillColor(fillcolor);
  box->SetLineWidth(2);
  fBoxes.push_back(box);
}

void CPlot::AddBox(double x1, double y1, double x2, double y2, 
                   int linecolor, int linesty, int fillcolor, TString label)
{
  TBox *box = new TBox(x1,y1,x2,y2);
  box->SetLineColor(linecolor);
  box->SetLineStyle(linesty);
  if(fillcolor==-1)
    box->SetFillStyle(0);
  else
    box->SetFillColor(fillcolor);
  box->SetLineWidth(2);
  fBoxes.push_back(box);

  if(!fLeg)
    fLeg = new TLegend(0.6,0.84,0.93,0.9);
  else
    fLeg->SetY1(fLeg->GetY1()-0.06);
  
  if(fillcolor<0) fLeg->AddEntry(box,label,"L");
  else            fLeg->AddEntry(box,label,"F");
}

//--------------------------------------------------------------------------------------------------
void CPlot::AddFcn(TF1* fcn, int color, int linesty)
{
  if(!fcn)
    return;
  
  fcn->SetLineColor(color);
  fcn->SetLineStyle(linesty);
  fFcns.push_back(fcn);
}

void CPlot::AddFcn(TF1* fcn, TString label, int color, int linesty)
{
  if(!fcn)
    return;
    
  if(!fLeg)
    fLeg = new TLegend(0.6,0.84,0.93,0.9);
  else
    fLeg->SetY1(fLeg->GetY1()-0.06);
  fLeg->AddEntry(fcn,label,"L");
  
  AddFcn(fcn,color,linesty);
} 

//--------------------------------------------------------------------------------------------------
void CPlot::Draw(TCanvas *c, bool doSave, TString format)
{ 
  //c->cd();
  
  c->SetLogy(fLogy);
  c->SetLogx(fLogx);
  
  if(!fItems.size() && !fRooPlot)
    return;   
  
  if(fRooPlot) {        
    fRooPlot->Draw();   
  }
      
  int nHist1F=0, nHist1D=0, nHist2D=0, nGraph=0, nProf=0;
  for(UInt_t i=0; i<fItems.size(); i++) {
    if(fItems[i].hist1F != 0) nHist1F++;
    if(fItems[i].hist1D != 0) nHist1D++;
    if(fItems[i].hist2D != 0) nHist2D++;
    if(fItems[i].graph  != 0) nGraph++;
    if(fItems[i].prof != 0) nProf++;
  }
  
  //
  // Draw 2D histogram, save if necessary, then exit
  //   Suggested options for:
  //     contour plot -> "CONT4Z"
  //     lego plot    -> "LEGO1 0"
  //     color plot   -> "COLZ"
  //   Default is scatter plot
  //  
  if(nHist2D>0) {
    for(UInt_t i=0; i<fItems.size(); i++) {
      if(fItems[i].hist2D==0) continue;
      
      fItems[i].hist2D->Draw(fItems[i].drawopt);
      fItems[i].hist2D->SetTitle(fTitle);
      fItems[i].hist2D->GetXaxis()->SetTitle(fXTitle);
      fItems[i].hist2D->GetYaxis()->SetTitle(fYTitle);
    
      //
      // Set log scale if necessary
      // 
      c->SetLogx(fLogx);
      c->SetLogy(fLogy);
      
      for(UInt_t ii=0; ii<fLines.size(); ii++)
        fLines[i]->Draw();

      for(UInt_t ii=0; ii<fBoxes.size(); ii++)
        fBoxes[i]->Draw();
      
      for(UInt_t j=0; j<fTextBoxes.size(); j++)
        fTextBoxes[j]->Draw();
      
      for(UInt_t j=0; j<fLatex.size(); j++)
        fLatex[j]->Draw();
            
      if(doSave) {
        gSystem->mkdir(sOutDir,true);
        TString outname = sOutDir+TString("/")+fName+TString(".");
	if(format.CompareTo("all",TString::kIgnoreCase)==0) {
	  c->SaveAs(outname+TString("png"));
	  c->SaveAs(outname+TString("eps"));
	  c->SaveAs(outname+TString("C"));
	} else {
	  //std::cout << "format=" << format << ", format[0]=" << format[0] << ", is '*' ?"  << ((format[0]=='*')?1:0) << "\n";
	  if (format[0]=='*') {
	    const char *tst=format.Data();
	    if (strstr(tst,"gif")) c->SaveAs(outname+"gif");
	    if (strstr(tst,"png")) c->SaveAs(outname+"png");
	    if (strstr(tst,"jpg")) c->SaveAs(outname+"jpg");
	    if (strstr(tst,"eps")) c->SaveAs(outname+"eps");
	    if (strstr(tst,"C")) c->SaveAs(outname+"C");
	  }
	  else c->SaveAs(outname+format);
	}
      }
      
      return;
    }
  }    
  
  // 
  // Draw 1D histograms
  //   Histograms are cloned so that content and properties 
  //   of the original histograms are not changed
  //
  std::cout << "draw " << nHist1F << "(float) and " << nHist1D << "(double precision) 1D histograms\n";
  std::vector<TH1F*> vHistsF;
  std::vector<TH1D*> vHistsD;
  std::vector<TString> vHistOptsF, vHistOptsD;
  if(nHist1F+nHist1D>0) {   
    
    double ymax=0;
    UInt_t ifirstF=0, fFirst=1;
    
    for(UInt_t i=0; i<fItems.size(); i++) {
      int is_hf=(fItems[i].hist1F!=0) ? 1:0;
      int is_hd=(fItems[i].hist1D!=0) ? 1:0;
      if (!is_hf && !is_hd) continue;
      if(is_hf && fStack && fStack->GetHists()->Contains(fItems[i].hist1F)) continue;
      
      //fItems[i].hist1D->DrawCopy();
      TString hname = fName;
      hname += "_h_";
      hname += i;
    
      //std::cout << "fRebin=" << fRebin << "\n";
      
      if (is_hf) {
	TH1F *h; 
	if(fRebin>1) {
	  std::cout << "rebin\n";
	  h = (TH1F*)fItems[i].hist1F->Rebin(fRebin,hname);
	}
	else {
	  h = (TH1F*)fItems[i].hist1F->Clone(hname);
	}
	//h->SetDirectory(0);
	
	if(fXmin < fXmax)
	  h->GetXaxis()->SetRangeUser(fXmin,fXmax);
	
	if(fYmin < fYmax) {
	  h->GetYaxis()->SetRangeUser(fYmin,fYmax);
	} else {
	  if(ymax < h->GetMaximum()) {
	    ymax = h->GetMaximum();
	    ifirstF = vHistsF.size();
	    fFirst=1;
	  }
	}
	vHistsF.push_back(h);
	vHistOptsF.push_back(fItems[i].drawopt);
      }
      if (is_hd) {
	TH1D *h; 
	if(fRebin>1) {
	  std::cout << "rebin\n";
	  h = (TH1D*)fItems[i].hist1D->Rebin(fRebin,hname);
	}
	else {
	  h = (TH1D*)fItems[i].hist1D->Clone(hname);
	}
	//h->SetDirectory(0);
	
	if(fXmin < fXmax)
	  h->GetXaxis()->SetRangeUser(fXmin,fXmax);
	
	if(fYmin < fYmax) {
	  h->GetYaxis()->SetRangeUser(fYmin,fYmax);
	} else {
	  if(ymax < h->GetMaximum()) {
	    ymax = h->GetMaximum();
	    ifirstF = vHistsF.size();
	    fFirst=0;
	  }
	}
	vHistsD.push_back(h);
	vHistOptsD.push_back(fItems[i].drawopt);
      }
    }

    //
    // Draw histogram stack
    //
    if(fStack) {
      if(vHistsF.size()>0) { 
        if(fYmin < fYmax) {
	  fStack->Draw("hist same");
	} else {	
	  if(fStack->GetMaximum() > ymax) {
	    fStack->SetTitle(fTitle);
	    fStack->Draw("hist");
	  } else {	    
	    fStack->Draw("hist same"); 
	  }
        }
	
      } else {
        // NOTE: Must draw first before accessing axes
	fStack->Draw("hist"); 

	if(fXmin < fXmax)
          fStack->GetXaxis()->SetRangeUser(fXmin,fXmax);
	 
	if(fYmin < fYmax) {
          fStack->SetMaximum(fYmax);  
	  fStack->SetMinimum(fYmin);
        }
	 
        fStack->SetTitle(fTitle);
	fStack->GetXaxis()->SetTitle(fXTitle);
	fStack->GetYaxis()->SetTitle(fYTitle);	
        fStack->Draw("hist");	 
      } 
    }
     
    int idraw=0;
    for(UInt_t i=0; i<vHistsF.size(); i++,idraw++) {
      TH1F *h = vHistsF[i];              
      h->SetLineWidth(2);
      char opt[100];
      if (idraw==0) sprintf(opt,"%s",vHistOptsF[i].Data());
      else sprintf(opt,"same %s",vHistOptsF[i].Data());    
      h->Draw(opt);
      //std::cout << "opt for TH1F i=" << i << " is " << opt << "\n";
      //h->DrawCopy(opt);
    }
    for(UInt_t i=0; i<vHistsD.size(); i++,idraw++) {
      TH1D *h = vHistsD[i];
      h->SetLineWidth(2);
      char opt[100];
      if (idraw==0) sprintf(opt,"%s",vHistOptsD[i].Data());
      else sprintf(opt,"same %s",vHistOptsD[i].Data());    
      h->Draw(opt);
      //std::cout << "opt for TH1D i=" << i << " is " << opt << "\n";
      //h->DrawCopy(opt);
    }
  }
  
  //
  // Draw graphs
  //
  std::vector<TGraph*> vGraphs;
  std::vector<TString> vGraphOpts;
  if(nGraph>0) {    
    for(UInt_t i=0; i<fItems.size(); i++) {
      if(fItems[i].graph==0) continue;
    
      TString grName = fName;
      grName += "_gr_";
      grName += i;
      
      TGraph *gr = (TGraph*)fItems[i].graph->Clone(grName);
      
      if(fXmin < fXmax)
        gr->GetXaxis()->SetLimits(fXmin,fXmax);
//        gr->GetXaxis()->SetRangeUser(fXmin,fXmax);
    
      if(fYmin < fYmax)
        gr->GetYaxis()->SetRangeUser(fYmin,fYmax);
	
      vGraphs.push_back(gr);
      vGraphOpts.push_back(fItems[i].drawopt);
    }
    
    if(vHistsF.size()==0) {
      vGraphs[0]->SetTitle(fTitle);
      vGraphs[0]->GetXaxis()->SetTitle(fXTitle);
      vGraphs[0]->GetYaxis()->SetTitle(fYTitle);
    }
    
    for(UInt_t i=0; i<vGraphs.size(); i++) {
      TGraph *gr = vGraphs[i];
      char opt[100];
      (i==0 && nHist1D==0) ? sprintf(opt,"AP%s",vGraphOpts[i].Data()) : sprintf(opt,"P%s",vGraphOpts[i].Data());
      gr->Draw(opt);
    }
  }

  //
  // Draw profile histograms
  //
  std::vector<TProfile*> vProfiles;
  std::vector<TString> vProfileOpts;
  if(nProf>0) {    
    for(UInt_t i=0; i<fItems.size(); i++) {
      if(fItems[i].prof==0) continue;
    
      TString prName = fName;
      prName += "_pr_";
      prName += i;
      
      TProfile *pr = (TProfile*)fItems[i].prof->Clone(prName);
      
      if(fXmin < fXmax)
        pr->GetXaxis()->SetLimits(fXmin,fXmax);
//        pr->GetXaxis()->SetRangeUser(fXmin,fXmax);
    
      if(fYmin < fYmax)
        pr->GetYaxis()->SetRangeUser(fYmin,fYmax);
	
      vProfiles.push_back(pr);
      vProfileOpts.push_back(fItems[i].drawopt);
    }
    
    if(vHistsF.size()==0) {
      vProfiles[0]->SetTitle(fTitle);
      vProfiles[0]->GetXaxis()->SetTitle(fXTitle);
      vProfiles[0]->GetYaxis()->SetTitle(fYTitle);
    }
    
    for(UInt_t i=0; i<vProfiles.size(); i++) {
      TProfile *pr = vProfiles[i];
      char opt[100];
      if(i>0 || nHist1F>0 || nGraph>0) 
        sprintf(opt,"same%s",vProfileOpts[i].Data());
      else
        sprintf(opt,"%s",vProfileOpts[i].Data());
      pr->Draw(opt);
    }
  }
      
  //
  // Draw legend
  //
  if(fLeg) {
    fLeg->SetFillStyle(0);
    fLeg->SetBorderSize(0);
    fLeg->Draw();
  }
  
  //
  // Draw statistics box
  //
  TLatex *stat=0, *sval=0;
  if(fShowStats) {
    char buffer[20];
    stat = new TLatex[3*(vHistsF.size()+vHistsD.size())];
    sval = new TLatex[3*(vHistsF.size()+vHistsF.size())];
    for(UInt_t i=0; i<vHistsF.size()+vHistsD.size(); i++) {
      UInt_t idx=3*i;
      UInt_t ii=i-vHistsF.size();
      int isHF=(i<vHistsF.size())?1:0;
      Color_t lineColor=(isHF) ? vHistsF[i]->GetLineColor() : vHistsD[ii]->GetLineColor();
      double integral=(isHF) ? vHistsF[i]->Integral() : vHistsD[ii]->Integral();
      double mean=(isHF) ? vHistsF[i]->GetMean() : vHistsD[ii]->GetMean();
      double rms=(isHF) ? vHistsF[i]->GetRMS() : vHistsD[ii]->GetRMS();
      int x = fShowStats;
      
      // number of entries
      if(x / 100) {
        stat[idx].SetNDC(); stat[idx].SetTextAlign(13); stat[idx].SetTextSize(0.03);
        stat[idx].SetText(fStatsX,fStatsY-0.04*(idx)-0.005*i,"Entries");
	stat[idx].SetTextColor(lineColor);
        stat[idx].Draw();   
        sprintf(buffer,"%i",int(integral));
        sval[idx].SetNDC(); sval[idx].SetTextAlign(33); sval[idx].SetTextSize(0.03);
        sval[idx].SetText(fStatsX+0.25,fStatsY-0.04*(idx)-0.005*i,buffer);
	sval[idx].SetTextColor(lineColor);
        sval[idx].Draw();
      }
      
      // mean
      x = x % 100;
      if(x / 10) {
        stat[idx+1].SetNDC(); stat[idx+1].SetTextAlign(13); stat[idx+1].SetTextSize(0.03);
        stat[idx+1].SetText(fStatsX,fStatsY-0.04*(idx+1)-0.005*i,"Mean");
	stat[idx+1].SetTextColor(lineColor);
        stat[idx+1].Draw();   
        sprintf(buffer,"%g",mean);
        sval[idx+1].SetNDC(); sval[idx+1].SetTextAlign(33); sval[idx+1].SetTextSize(0.03);
        sval[idx+1].SetText(fStatsX+0.25,fStatsY-0.04*(idx+1)-0.005*i,buffer);
	sval[idx+1].SetTextColor(lineColor);
        sval[idx+1].Draw();
      }
      
      // RMS
      x = x % 10;
      if(x) {
        stat[idx+2].SetNDC(); stat[idx+2].SetTextAlign(13); stat[idx+2].SetTextSize(0.03);
        stat[idx+2].SetText(fStatsX,fStatsY-0.04*(idx+2)-0.005*i,"RMS");
	stat[idx+2].SetTextColor(lineColor);
        stat[idx+2].Draw();   
        sprintf(buffer,"%g",rms);
        sval[idx+2].SetNDC(); sval[idx+2].SetTextAlign(33); sval[idx+2].SetTextSize(0.03);
        sval[idx+2].SetText(fStatsX+0.25,fStatsY-0.04*(idx+2)-0.005*i,buffer);
	sval[idx+2].SetTextColor(lineColor);
        sval[idx+2].Draw();
      }
    }
  }
 
  //
  // Draw functions
  //
  for(UInt_t i=0; i<fFcns.size(); i++)
    (i==0 && vHistsF.size()==0 && vGraphs.size()==0) ? fFcns[i]->Draw() : fFcns[i]->Draw("sameC");
  
  //
  // Draw lines
  //
  for(UInt_t i=0; i<fLines.size(); i++)
    fLines[i]->Draw();
  
  //
  // Draw Boxes
  //
  for(UInt_t i=0; i<fBoxes.size(); i++)
    fBoxes[i]->Draw();
  
  //
  // Draw textboxes
  //
  for(UInt_t i=0; i<fTextBoxes.size(); i++)
    fTextBoxes[i]->Draw();    
        
  for(UInt_t j=0; j<fLatex.size(); j++)
    fLatex[j]->Draw();

  //
  // Set log scale if necessary
  // 
  c->SetLogx(fLogx);
  c->SetLogy(fLogy);
  
  //
  // Set grid lines if necessary
  //
  c->SetGridx(fGridx);
  c->SetGridy(fGridy);
  
  //
  // Save plot if necessary
  //  
  if(doSave) {
    gSystem->mkdir(sOutDir,true);
    TString outname = sOutDir+TString("/")+fName+TString(".");
    if(format.CompareTo("all",TString::kIgnoreCase)==0) {
      c->SaveAs(outname+TString("png"));
      c->SaveAs(outname+TString("eps"));
      c->SaveAs(outname+TString("C"));
    } else {
      if (format[0]=='*') {
	const char *tst=format.Data();
	if (strstr(tst,"gif")) c->SaveAs(outname+"gif");
	if (strstr(tst,"png")) c->SaveAs(outname+"png");
	if (strstr(tst,"jpg")) c->SaveAs(outname+"jpg");
	if (strstr(tst,"eps")) c->SaveAs(outname+"eps");
	if (strstr(tst,"C")) c->SaveAs(outname+"C");
      }
      else c->SaveAs(outname+format);
    }
    
    delete [] stat;
    delete [] sval;
//    for(UInt_t i=0; i<vHists.size(); i++)
//      delete vHists[i];
  }
}

// -------------------------------------------------------------------
// -------------------------------------------------------------------

int TestLogPlot(CPlot *CP, TCanvas *c, const std::vector<std::string> &otherPlotOptions) {
  int logXPlot=0;
  assert(CP); assert(c);
  for (UInt_t t=0; !logXPlot && (t<otherPlotOptions.size()); ++t) {
    if (otherPlotOptions[t].find("-logx")!=std::string::npos) {
      logXPlot=1;
      c->SetLogx();
      CP->SetLogx(1);
      //TH1 *hframe = c->DrawFrame(6,0.4, 1000, 1.0);
    }
  } 
  return logXPlot;
}

// -------------------------------------------------------------------

void ProcessOtherPlotOptions(CPlot *CP, TCanvas *c, int plot_case, const std::vector<std::string> &otherPlotOptions, Bool_t &doSave, TString &format) {
  assert(CP); assert(c);
  for (unsigned int t=0; t<otherPlotOptions.size(); ++t) {
    std::cout << "otherPlotOptions[" << t << "]=" << otherPlotOptions[t] << "\n";
    if (otherPlotOptions[t].find("-yrange")!=std::string::npos) {
      double ymin=atof(otherPlotOptions[t+1].c_str());
      double ymax=atof(otherPlotOptions[t+2].c_str());
      CP->SetYRange(ymin,ymax);
      t+=2;
    }
    else if (otherPlotOptions[t].find("-xrange")!=std::string::npos) {
      double xmin=atof(otherPlotOptions[t+1].c_str());
      double xmax=atof(otherPlotOptions[t+2].c_str());
      CP->SetXRange(xmin,xmax);
      //TH1 *hframe = c->DrawFrame(xmin,CP->GetYMin(), xmax, CP->GetYMax());
      c->DrawFrame(xmin,CP->GetYMin(), xmax, CP->GetYMax());
      t+=2;
    }
    else if (otherPlotOptions[t].find("-logx")!=std::string::npos) {
      //CP->SetLogx();
      //TH1 *hframe = c->DrawFrame(6,0.4, 700, 1.0);
      CP->SetXRange(9,1000);
    }
    else if (otherPlotOptions[t].find("-logy")!=std::string::npos) {
      CP->SetLogy();
    }
    else if (otherPlotOptions[t].find("-tag")!=std::string::npos) {
      switch(plot_case) {
      case 0: break; 
      case 1: CP->SetYTitle("p_{T,tag}"); break;
      case 2: CP->SetYTitle("\\phi_{tag}"); break;
      case 3: CP->SetYTitle("\\eta_{tag}"); break;
      }
    }
    else if (otherPlotOptions[t].find("-save-format")!=std::string::npos) {
      format=otherPlotOptions[t+1];
      std::cout << "saving to fig." << format << "\n";
      t++;
    }
    else if (otherPlotOptions[t].find("-save")!=std::string::npos) {
      doSave=true;
    }
    else if (otherPlotOptions[t].find("-legend-lb")!=std::string::npos) {
	      std::cout << "special legend left-bottom\n";
	      CP->SetLegend(0.15,0.12,0.56,0.25);
    }
    else if (otherPlotOptions[t].find("-legend-rb")!=std::string::npos) {
      std::cout << "special legend right-bottom\n";
      CP->SetLegend(0.50,0.12,0.96,0.25);
    }
    else if (otherPlotOptions[t].find("-legend-rt")!=std::string::npos) {
      std::cout << "special legend right-top\n";
      CP->SetLegend(0.50,0.78,0.96,0.85);
    }
  }
  return;
}

// -------------------------------------------------------------------
