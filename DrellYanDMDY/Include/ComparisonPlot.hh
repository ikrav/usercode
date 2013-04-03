#ifndef ComparisonPlot_HH
#define ComparisonPlot_HH


#include <TLatex.h>
#include "../Include/CPlot.hh"
#include "../Include/DYTools.hh"
#include "../Include/MyTools.hh"
//#include <TGraphAsymmErrors.h>


// -------------------------------------------------------------

class ComparisonPlot_t : public CPlot {
public:
  typedef enum { _ratioPlain=1, _ratioRel=2 } TRatioType_t;
public:
  unsigned int fRefIdx;   // index of histogram to use as a reference
  unsigned int fErrorsOnRatios;
  TRatioType_t fRatioType;
  vector<int> fHRatioIndices;
  vector<unsigned int> fExcludeIndices; // indices of histograms to exclude from a comparison
  TCanvas *canvas;
  TPad *padMain,*padRatio;
  TString fRatioLabel;
  vector<TH1F*> fHRatioItems;

  double fRatioYMin,fRatioYMax;
  int fRatioNdivisions, fPrintValues, fPrintRatios, fPrintRatioNames;
  double fRatioYTitleSize,fRatioYLabelSize;
  double fRatioYTitleOffset;
  double fRatioXTitleSize;

  ComparisonPlot_t(TRatioType_t ratioType, const TString &name, const TString &title, const TString &xtitle, const TString &ytitle, const TString &ratioYLabel) :
    CPlot(name,title,xtitle,ytitle),
    fRefIdx(0), fErrorsOnRatios(1),
    fRatioType(ratioType),
    fHRatioIndices(),
    fExcludeIndices(),
    canvas(NULL),
    padMain(NULL), padRatio(NULL),
    fRatioLabel(ratioYLabel), fHRatioItems(),
    fRatioYMin(0.), fRatioYMax(0.),
    fRatioNdivisions(805),
    fPrintValues(0),
    fPrintRatios(0), fPrintRatioNames(0),
    fRatioYTitleSize(0.15), 
    fRatioYLabelSize(0.13),
    fRatioYTitleOffset(0.45),
    fRatioXTitleSize(0.15)
  {}

  ComparisonPlot_t(const ComparisonPlot_t &cp, 
		   const TString &newName, const TString &newTitle="") :
    CPlot(newName,newTitle,cp.fXTitle,cp.fYTitle),
    fRefIdx(0), fErrorsOnRatios(cp.fErrorsOnRatios),
    fRatioType(cp.fRatioType),
    fHRatioIndices(),
    fExcludeIndices(),
    canvas(NULL),
    padMain(NULL), padRatio(NULL),
    fRatioLabel(cp.fRatioLabel), fHRatioItems(),
    fRatioYMin(cp.fRatioYMin), fRatioYMax(cp.fRatioYMax),
    fRatioNdivisions(cp.fRatioNdivisions),
    fPrintValues(cp.fPrintValues),
    fPrintRatios(cp.fPrintRatios),
    fPrintRatioNames(cp.fPrintRatioNames),
    fRatioYTitleSize(cp.fRatioYTitleSize),
    fRatioYLabelSize(cp.fRatioYLabelSize),
    fRatioYTitleOffset(cp.fRatioYTitleOffset),
    fRatioXTitleSize(cp.fRatioXTitleSize)
  {}


  void ErrorsOnRatios(unsigned int on=1) { fErrorsOnRatios=on; }
  void SetRatioYRange(double ymin, double ymax) { fRatioYMin=ymin; fRatioYMax=ymax; }
  void SetRatioNdivisions(int cnt) { fRatioNdivisions=cnt; }
  void SetRatioLabelSize(double size) { fRatioYLabelSize=size; }
  void SetRatioYTitleSize(double size, double offset=-99) { fRatioYTitleSize=size; if (offset>=0) fRatioYTitleOffset=offset; }
  void SetRatioXTitleSize(double size) { fRatioXTitleSize=size; }
  void SetRatioTitleOffset(double size) { fRatioYTitleOffset=size; }
  void SetRefIdx(int refIdx) { fRefIdx=refIdx; }
  void SetPrintValues(int yes=1) { fPrintValues=yes; }
  void SetPrintRatio(int yes=1) { fPrintRatios=yes; }
  void SetPrintRatios(int yes=1) { fPrintRatios=yes; }
  void SetPrintRatioNames(int yes=1) { fPrintRatioNames=yes; }

  void SetRefIdx(const TH1F* h) {
    if (!h) { 
      std::cout << "error ComparisonPlot::SetRefIdx(ptr=NULL)\n";
      assert(0);
    }
    fRefIdx=-1;
    for (unsigned int i=0; i<fItems.size(); ++i) {
      if (fItems[i].hist1D==h) {
	fRefIdx=i;
      }
    }
    if (fRefIdx==(unsigned int)(-1)) {
      std::cout << "ComparisonPlot::SetRefIdx(ptr): failed to determine fRefIdx\n";
      fRefIdx=0;
    }
  }
  
  void SkipInRatioPlots(const TH1F* h) {
    if (!h) {
      std::cout << "error ComparisonPlot::SkipInRatioPlots(ptr=NULL)\n";
      assert(0);
    }
    unsigned int idx=-1;
    for (unsigned int i=0; i<fItems.size(); ++i) {
      if (fItems[i].hist1D==h) {
	idx=i;
      }
    }
    if (idx==(unsigned int)(-1)) {
      std::cout << "ComparisonPlot::SkipInRatioPlots(ptr): failed to determine idx\n";
    }
    else {
      fExcludeIndices.push_back(idx);
    }
  }

  void Skip(TH1F* h) { SkipInRatioPlots(h); }

  TLatex* AddTextCMSPreliminary(double x=0.93, double y=0.94) {
    TLatex *cmsText = new TLatex();
    cmsText->SetTextFont(42);
    cmsText->SetTextSize(0.055);
    cmsText->SetTextAlign(31);
    cmsText->SetNDC();
    cmsText->SetText(x, y, "#it{CMS Preliminary}");
    cmsText->Draw();
    return cmsText;
  }

  TLatex* AddTextCMSSimulation(double x=0.93, double y=0.94) {
    TLatex* txt=this->AddTextCMSPreliminary();
    txt->SetText(x, y, "#it{CMS Simulation}");
    return txt;
  }

  TLatex* AddTextLumi(double x=0.91, double y=0.90) {
    TLatex *lumiText = new TLatex();
    lumiText->SetTextFont(42);
    lumiText->SetTextSize(0.05);
    lumiText->SetTextAlign(33);
    lumiText->SetNDC();
    lumiText->SetText(x, y, DYTools::strLumiAtECMS);
    lumiText->Draw();
    return lumiText;
  }

  void PreparePads(TCanvas *c, int subpad1=1, int subpad2=2, TString name="comp", double padYDivPoint=0.29) {
    const double dy=0.005;

    TPad *pad1 = (TPad*)c->GetPad(subpad1);
    TPad *pad2 = (TPad*)c->GetPad(subpad2);
    if (fLogx) { pad1->SetLogx(); pad2->SetLogx(); }
    if (fLogy) { pad1->SetLogy(); }
    //pad1->SetFillColor(11);
    //pad2->SetFillColor(11);

    double xmin1= pad1->GetXlowNDC();
    double xmax1= xmin1 + pad1->GetWNDC();
    double ymin1= pad1->GetYlowNDC();
    double ymax1= ymin1 + pad1->GetHNDC();
    double xmin2= pad2->GetXlowNDC();
    double xmax2= xmin2 + pad2->GetWNDC();
    double ymin2= pad2->GetYlowNDC();
    double ymax2= ymin2 + pad2->GetHNDC();
    double xmin=(xmin1 < xmin2) ? xmin1 : xmin2;
    double xmax=(xmax1 > xmax2) ? xmax1 : xmax2;
    double ymin=(ymin1 < ymin2) ? ymin1 : ymin2;
    double ymax=(ymax1 > ymax2) ? ymax1 : ymax2;
    double ydiv=ymin + (ymax-ymin)*padYDivPoint;

    pad1->SetPad(name+ TString("_main"),"",xmin,ydiv+dy,xmax,ymax);
    pad2->SetPad(name + TString("_ratio"),"",xmin,ymin,xmax,ydiv-dy);

    pad1->SetFillColor(0);
    pad2->SetFillColor(0);

    const double leftM= c->GetLeftMargin();
    const double rightM= c->GetRightMargin();
    pad1->SetTopMargin(0.08);
    pad1->SetBottomMargin(0.015); // All X axis labels and titles are thus cut off
    pad1->SetRightMargin(rightM);
    pad1->SetLeftMargin(leftM);
    //p->SetFillStyle(0);

    pad2->SetTopMargin(0.045);
    pad2->SetTopMargin(0.14);
    pad2->SetBottomMargin(0.45);
    pad2->SetRightMargin(rightM);
    pad2->SetLeftMargin(leftM);
    //p->SetFillStyle(0);

    pad1->SetTicky();
    pad1->SetFrameLineWidth(2);
    pad2->SetTicky();
    pad2->SetFrameLineWidth(2);

    pad1->Draw();
    pad2->Draw();
  }

  void Prepare2Pads(TCanvas *c) {
    c->Divide(2,1);
    PreparePads(c); 
  }

  void Prepare6Pads(TCanvas *c, int landscape, TString name="comp", double padYDivPoint=0.29) {
    TString sa="_a";
    TString sb="_b";
    TString sc="_c";
    TString sd="_d";
    TString se="_e";
    TString sf="_f";

    if (landscape) {
      c->Divide(3,4);
      PreparePads(c, 1,4, name + sa, padYDivPoint);
      PreparePads(c, 2,5, name + sb, padYDivPoint);
      PreparePads(c, 3,6, name + sc, padYDivPoint);
      PreparePads(c, 7,10, name + sd, padYDivPoint);
      PreparePads(c, 8,11, name + se, padYDivPoint);
      PreparePads(c, 9,12, name + sf, padYDivPoint);
    }
    else {
      c->Divide(2,6);
      PreparePads(c, 1,3, name + sa, padYDivPoint);
      PreparePads(c, 2,4, name + sb, padYDivPoint);
      PreparePads(c, 5,6, name + sc, padYDivPoint);
      PreparePads(c, 7,8, name + sd, padYDivPoint);
      PreparePads(c,  9,11, name + se, padYDivPoint);
      PreparePads(c, 10,12, name + sf, padYDivPoint);      
    }

    for (int i=1; i<=12; i++) {
      TPad *pad=(TPad*)c->GetPad(i);
      pad->SetFillColor(0);
      pad->SetTickx(1);
      pad->SetTicky(1);
    }
  }

  void Draw6(TCanvas *c, int landscape, int idx, bool doSave=false, TString format="png", int *subpad1=NULL, int *subpad2=NULL) {
    if ((landscape!=0) && (landscape!=1)) {
      std::cout << "Draw6 assumes landscape=0 or 1\n";
      return;
    }
    int shift=2+landscape;
    int padIdx1= idx + shift*int((idx-1)/shift);
    int padIdx2= padIdx1 + shift;
    Draw(c,doSave,format,padIdx1,padIdx2);
    if (subpad1) *subpad1=padIdx1;
    if (subpad2) *subpad2=padIdx2;
  }


  void Draw(TCanvas *c, bool doSave=false, TString format="png", int subpad1=1, int subpad2=2) {
    if (canvas!=c) canvas=c;
    padMain = (TPad*)c->GetPad(subpad1);
    padRatio = (TPad*)c->GetPad(subpad2);

    //padMain->Draw();
    //padRatio->Draw();
    
    std::vector<TString> savedTitles;
    savedTitles.reserve(fItems.size());
    for (unsigned int i=0; i<fItems.size(); ++i) {
      if (fItems[i].hist1D!=0) {
	fHRatioIndices.push_back(i);
	fItems[i].hist1D->GetXaxis()->SetLabelOffset(99.05);
	savedTitles.push_back(fItems[i].hist1D->GetXaxis()->GetTitle());
	fItems[i].hist1D->GetXaxis()->SetTitle("");
      }
    }
    if ( fHRatioIndices.size() ==0 ) return;

    TH1F *hRef=fItems[fHRatioIndices[fRefIdx]].hist1D;

    CPlot::Draw(c,false,"png",subpad1);

    padRatio->cd();
    TVectorD hratioActive(fHRatioIndices.size());
    hratioActive=0;

    if (1 && fHRatioIndices.size()) {
      //unsigned int imin=(fHRatioIndices.size()==1) ? 0 : 1;
      unsigned int imin=0;
      for (unsigned int k=0; k<fExcludeIndices.size(); ++k) {
	if (imin==fExcludeIndices[k]) imin++; else break;
      }
      std::cout << "imin=" << imin << "\n";

      if (fPrintValues) printHisto(fItems[fRefIdx].hist1D);

      // prepare ratios
      double yrMin=1e9, yrMax=-1e9;
      for (unsigned int i=imin; i<fHRatioIndices.size(); ++i) {
	CPlotItem *item= &fItems[fHRatioIndices[i]];
	TH1F *histo=item->hist1D;
	TH1F *hratio=(TH1F*)histo->Clone(TString("ratio_of_") + histo->GetName());
	hratio->SetMarkerStyle(histo->GetMarkerStyle());
	hratio->SetDirectory(0);
	hratio->GetXaxis()->SetLabelOffset(0.05);
	fHRatioItems.push_back(hratio);

	int skip=(i==fRefIdx) ? 1:0;
	for (unsigned int k=0; !skip && (k<fExcludeIndices.size()); ++k) 
	  skip= (i==fExcludeIndices[k]) ? 1:0;
	hratioActive[fHRatioItems.size()-1]=1-skip;
	if (skip) { 
	  std::cout << "skipping entry #" << i << " (" << histo->GetName() << ")\n";
	  continue;
	}

	if (fPrintValues) printHisto(histo,1);

	switch(fRatioType) {
	case _ratioPlain: 
	  hratio->Divide(histo,hRef); 
	  break;
	case _ratioRel:  // (ref-h)/ref
	  hratio->Scale(-1);
	  hratio->Add(hRef);
	  hratio->Divide(hRef);
	  break;
	default:
	  std::cout << "cannot prepare ratio\n";
	  return;
	}
	hratio->SetTitle("");
	if (!fErrorsOnRatios) {
	  for (int ibin=1; ibin<=hratio->GetNbinsX(); ++ibin) {
	    hratio->SetBinError(ibin,0);
	  }
	}
	//printHisto(std::cout,hratio);
	if (yrMin > hratio->GetMinimum()) yrMin = hratio->GetMinimum();
	if (yrMax < hratio->GetMaximum()) yrMax = hratio->GetMaximum();
      }

      double yrC=fabs(yrMax);

      if (fRatioYMin==fRatioYMax) {
	std::cout << "yrMin=" << yrMin << ", yrMax=" << yrMax << "\n";
	if (fRatioType==_ratioRel) {
	  if (fabs(yrMin)<1) yrMin=trunc(yrMin*110)*0.01;
	  else yrMin=trunc(yrMin*11)*0.1;
	  if (fabs(yrMax)<1) yrMax=trunc(yrMax*110)*0.01;
	  else yrMax=trunc(yrMax*11)*0.1;
	  std::cout << "1st correction: yrMin=" << yrMin << ", yrMax=" << yrMax << "\n";
	}
	if (yrC<fabs(yrMin)) {
	  yrC=fabs(yrMin);
	  if (yrMax>-1e-7) yrMax=yrC; else yrMax=-yrC;
	}
	else {
	  if (yrMin>-1e-7) {
	    if (yrMin*yrMax<0) yrMin=yrC; 
	  }
	  else yrMin=-yrC;
	}
	std::cout << "updated: yrMin=" << yrMin << ", yrMax=" << yrMax << "\n";
      }
      else {
	yrMin=fRatioYMin; yrMax=fRatioYMax;
	std::cout << "user: yrMin=" << yrMin << ", yrMax=" << yrMax << "\n";
      }

      // plot ratios
      int first=1;
      for (unsigned int i=imin; i<fHRatioIndices.size(); ++i) {
	int skip=(i==fRefIdx) ? 1:0;
	for (unsigned int k=0; !skip && (k<fExcludeIndices.size()); ++k) 
	  skip= (i==fExcludeIndices[k]) ? 1:0;
	if (skip) { 
	  std::cout << "skipping entry #" << i << "\n";
	  continue;
	}

	CPlotItem *item= &fItems[fHRatioIndices[i]];
	TH1F *hratio=fHRatioItems[i-imin];

	hratio->GetYaxis()->SetRangeUser(yrMin,yrMax);
	//hratio->SetMarkerStyle(kFullCircle);

	if (fPrintRatios) printHisto(std::cout,hratio);
	TString opt=item->drawopt;
	opt="P";
	if (first) {
	  first=0;
	  hratio->GetYaxis()->SetTitle(fRatioLabel);
	  hratio->GetXaxis()->SetTitle(fXTitle);
	  hratio->GetXaxis()->SetLabelSize( 3*hRef->GetXaxis()->GetLabelSize() );
	  hratio->GetXaxis()->SetTickLength(0.09);
	  hratio->GetXaxis()->SetTitleSize( 3*hRef->GetXaxis()->GetTitleSize() );
	  if (fLogx) {
	    hratio->GetXaxis()->SetMoreLogLabels();
	    hratio->GetXaxis()->SetNoExponent();
	  }

	  hratio->GetYaxis()->SetNdivisions(fRatioNdivisions);
	  hratio->GetYaxis()->SetLabelSize( fRatioYLabelSize );
	  hratio->GetYaxis()->SetTitleSize(fRatioYTitleSize);
	  hratio->GetYaxis()->SetTitleOffset(fRatioYTitleOffset);
	  hratio->GetXaxis()->SetLabelSize( fRatioYLabelSize );
	  hratio->GetXaxis()->SetTitleSize(fRatioXTitleSize);

	//hratio->SetMarkerColor(kOrange+8);

	  opt.ReplaceAll("same","");
	  opt.ReplaceAll("SAME","");
	}
	else opt.Append("same");
	hratio->Draw(opt);
      }

      double nominal=0;
      switch(fRatioType) {
      case _ratioPlain: nominal=1; break;
      case _ratioRel: nominal=0; break;
      default:
	nominal=yrC;
      }
      TLine *lineAtOne = 
	new TLine(hRef->GetXaxis()->GetXmin(), nominal,
		  hRef->GetXaxis()->GetXmax(), nominal);
      if ((nominal>yrMin) && (nominal<yrMax)) {
	lineAtOne->SetLineStyle(kDashed);
	lineAtOne->SetLineWidth(1);
	lineAtOne->SetLineColor(kBlue);
	lineAtOne->Draw();
      }
      TLine *linePlus10 = 
	new TLine(hRef->GetXaxis()->GetXmin(), nominal+0.1,
		  hRef->GetXaxis()->GetXmax(), nominal+0.1);
      if ((nominal+0.1>yrMin) && (nominal+0.1<yrMax)) {
	linePlus10->SetLineStyle(8);
	linePlus10->SetLineWidth(1);
	linePlus10->SetLineColor(kGreen+1);
	linePlus10->Draw();
      }
      TLine *lineMinus10 = 
	new TLine(hRef->GetXaxis()->GetXmin(), nominal-0.1,
		  hRef->GetXaxis()->GetXmax(), nominal-0.1);
      if ((nominal-0.1>yrMin) && (nominal-0.1<yrMax)) {
	lineMinus10->SetLineStyle(8);
	lineMinus10->SetLineWidth(1);
	lineMinus10->SetLineColor(kGreen+1);
	lineMinus10->Draw();
      }
      
      c->cd();
    }

    if(doSave) {
      cout << "DEBUG: saving" << endl;
      gSystem->mkdir(sOutDir,true);
      TString outname = sOutDir+TString("/")+fName+TString(".");
      if(format.CompareTo("all",TString::kIgnoreCase)==0) {
	c->SaveAs(outname+TString("png"));
	c->SaveAs(outname+TString("eps"));
	c->SaveAs(outname+TString("C"));
      } else {
	c->SaveAs(outname+format);
      }
    }

    if (fPrintRatioNames) {
      std::cout << "\nComparisonPlot::Draw -- printing ratio names:\n";
      for (unsigned int i=0; i<fHRatioItems.size(); ++i) {
	if (!hratioActive[i]) {
	  std::cout << "skipping entry #" << i << "\n";
	  continue;
	}
	std::cout << " - " << fHRatioItems[i]->GetName() << "\n";
      }
      std::cout << "Ratio computed with reference to " << fItems[fHRatioIndices[fRefIdx]].hist1D->GetName() << "\n";
      std::cout << "\n";
    }

    return;
  }

};

// -------------------------------------------------------------
// -------------------------------------------------------------

inline
void prepareMassRanges(std::vector<TLatex*> &labels, double x=0.91, double y=0.82, int color=kBlack) {
  labels.reserve(labels.size() + DYTools::nMassBins);
  TLatex *massLabel=NULL;
  for (int i=0; i<DYTools::nMassBins; i++) {
    TString massStr=Form("%2.0lf<M_{ee}<%2.0lf GeV",DYTools::massBinLimits[i],DYTools::massBinLimits[i+1]);
    massLabel=new TLatex();
    massLabel->SetTextFont(42);
    massLabel->SetTextSize(0.05);
    massLabel->SetTextAlign(33);
    massLabel->SetNDC();
    massLabel->SetText(x, y, massStr);
    massLabel->SetTextColor(color);
    labels.push_back(massLabel);
  }
}

// -------------------------------------------------------------


#endif
