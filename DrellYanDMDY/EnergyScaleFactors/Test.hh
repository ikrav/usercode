#ifndef Test_H
#define Test_H

#include <TH1F.h>
#include <TH1D.h>
#include <vector>
#include <string>
#include "CPlotMdf.hh"

// -------------------------------------------------------



// -------------------------------------------------------
// histoSetNames -- names of the histosV[i]
// histo_names -- names of the histosV[i][j]
// histo_descriptions -- several lines for histosV[i][j]
template<class HistoClass_t>
inline int PrintHTMLHistoLines_OverlayedDuos(const char *output_dir, const char *file_tag, const std::vector<const std::vector<HistoClass_t*>*> &histosV, const std::vector<std::string> &histoSetNames, const std::vector<std::string> &plotTableHeaders, const std::vector<TString> &histo_names, const std::vector<std::vector<std::string>*> &tables, int print_descriptions_separately, const char *description, std::vector<std::string> &html_lines, int start_doc=0, int end_doc=0, const std::vector<std::string> *extra_html_lines=NULL, const char *output_fname=NULL, const std::vector<std::string> *legend_placement=NULL) {
  const char *fncname="PrintHTMLHistoLines_OverlayedDuos";
#ifndef __myLib__
  const int saveCFormat=1;
#endif
  const int createCFiles=saveCFormat;

  CPlot::sOutDir=output_dir;
  //if (CPlot::sOutDir[CPlot::sOutDir.Length()-1]!='/') CPlot::sOutDir.Append("/");

  html_lines.reserve(html_lines.size()+1000);

  const int buflen=300;
  char buf[buflen];
  char plot_fname[100];
  const char *title=(description) ? description : "histos(overlayed duos)";
  unsigned int save_pos=0;

  if (start_doc) {
    sprintf(buf,"<!DOCTYPE html\n    PUBLIC ""-//W3C//DTD HTML 3.2//EN"">"); html_lines.push_back(buf);
    sprintf(buf,"<html>\n<head><title>%s (%s)</title></head>\n",title,file_tag); html_lines.push_back(buf);
    sprintf(buf,"<body bgcolor=""EEEEEE"">\n"); html_lines.push_back(buf);
  }
  sprintf(buf,"<h3 style=""text-align:left; color:DD6600;"">%s</h3>\n",title); html_lines.push_back(buf);
  save_pos=html_lines.size();


  const int histCount=(histosV.size()) ? histosV[0]->size() : 0;
  std::cout << "histCount=" << histCount << "\n";
  std::cout << "distribution: "; for (unsigned int i=0; i<histosV.size(); ++i) std::cout << " " << histosV[i]->size(); std::cout << "\n";
  if (print_descriptions_separately) {
    sprintf(buf,"<table border=""0"" cellspacing=""2"" width=""100%%"">\n"); html_lines.push_back(buf);
    if (plotTableHeaders.size()) {
      html_lines.push_back("<tr><td width=""15%%"">Histogram names</td>\n");
      for (unsigned int i=0; i<plotTableHeaders.size(); ++i) {
	const char *p=plotTableHeaders[i].c_str();
	int shift= (*p=='|') ? 2:0;
	sprintf(buf,"<td width=""15%%"">%s</td>\n",p+shift);
	html_lines.push_back(buf);
      }
      html_lines.push_back("</tr>\n");
      if (2*plotTableHeaders.size()!=histosV.size()) {
	std::cout << " WARNING in " << fncname << ": plotTableHeaders.size=" << plotTableHeaders.size() << ", while histosV.size=" << histosV.size() << " (should be double)\n";
      }
    }
    html_lines.push_back("</table>\n");
  }

  sprintf(buf,"<table border=""0"" cellspacing=""2"" width=""100%%"">\n"); html_lines.push_back(buf);
  for (int i=0; i<histCount; ++i) {
    html_lines.push_back("<tr>\n");
    sprintf(buf,"<td width=""10%%""><center>%s</center><br></td>\n",histo_names[i].Data());
    html_lines.push_back(buf);
    for (unsigned int hvi=0; hvi<histosV.size(); hvi+=2) {
      HistoClass_t *h1=((*histosV[hvi])[i]!=NULL) ? (HistoClass_t*)(*histosV[hvi])[i]->Clone("h1") : NULL;
      HistoClass_t *h2=(HistoClass_t*)(*histosV[hvi+1])[i]->Clone("h2");
      if (h1) h1->SetTitle(histo_names[i]);
      h2->SetTitle(histo_names[i]);
      TLatex *latex=NULL;
      //std::cout << "hvi=" << hvi << ", h1->YTitle=" << h1->GetYaxis()->GetTitle() << "\n";
      if (hvi<4) latex=CreateChiSquareText(h1->Chi2Test(h2,"WW CHI2"),0);
      if (strstr(h2->GetYaxis()->GetTitle(),"asymmetry")) {
	double ym1=(h1) ? h1->GetMaximum() : 0;
	double ym2=h2->GetMaximum();
	if (ym2>ym1) ym1=ym2;
	ym1 *= 1.05;
	if (h1) h1->SetMaximum(ym1);
	h2->SetMaximum(ym2);
      }
      const char *figFileExt_base="gif";
      sprintf(buf,"fig%s_%s_%d",file_tag,histo_names[i].Data(),hvi);
      char *p=strchr(buf,' '); while (p) { *p='_'; p=strchr(buf,' '); }
      p=strchr(buf,'-'); while (p) { *p='_'; p=strchr(buf,'-'); }
      sprintf(plot_fname,buf);
      
      TCanvas *c1= new TCanvas("fig","fig",640,640);
      gPad->SetLeftMargin(0.15);
      if (h1) h1->GetYaxis()->SetTitleOffset(1.6);
      h2->GetYaxis()->SetTitleOffset(1.6);
      h2->SetMarkerStyle(kFullDotLarge); h2->SetMarkerColor(kRed+2);
      CPlot *CP=new CPlot(plot_fname,plot_fname,h2->GetXaxis()->GetTitle(),h2->GetYaxis()->GetTitle());
      if (h1) CP->AddHist1D(h1, histoSetNames[hvi].c_str(), "hist ", kBlue+2); 
      CP->AddHist1D(h2, histoSetNames[hvi+1].c_str(), "LPE", kRed+2);
      //std::cout << "chi2=" << h1->Chi2Test(h2,"UU CHI2") << ", chi2/ndof=" << h1->Chi2Test(h2,"UU CHI2/NDOF") << "\n";
      CP->SetLegend(0.65,0.60,0.9,0.70);
      const std::vector<std::string>* lp=legend_placement;
      if (lp) {
	if (PosOk((*lp)[hvi].find("lt"))) CP->SetLegend(0.20,0.78,0.45,0.88);
      }
      if (latex) CP->AddLatex(latex);
      TString figFileExt=figFileExt_base;
      if (createCFiles) figFileExt=TString("*") + figFileExt_base + TString("C");
      //std::cout << "\tfigFileExt=" << figFileExt << "\n";
      CP->Draw(c1,1,figFileExt);
      //Pause();
      delete CP;
      delete c1;
      sprintf(plot_fname,"%s.%s",buf,figFileExt_base);
      sprintf(buf,"<td width=""10%%""><a target=""_blank"" href=""%s""><img src=""%s"" alt=""%s"" width=""100%%""></a></td>\n",plot_fname,plot_fname,plot_fname); html_lines.push_back(buf);
    }
    sprintf(buf,"</tr>\n"); html_lines.push_back(buf);
  }
  sprintf(buf,"</table>\n"); html_lines.push_back(buf);

  if (tables.size()>0) {
    html_lines.push_back("<hr />\n<H3>Tables</H3>\n");
    for (unsigned int ti=0; ti<tables.size(); ++ti) {
      const std::vector<std::string> *tbl=tables[ti];
      html_lines.push_back("<table>\n");
      for (unsigned int i=0; i<tbl->size(); ++i) {
	html_lines.push_back((*tbl)[i]);
      }
      //sprintf(buf,"</td>\n</tr>\n"); html_lines.push_back(buf);
      html_lines.push_back("</table>\n");
    }
  }


  if (end_doc) {
    sprintf(buf,"</body>\n</html>\n"); html_lines.push_back(buf);
  }
  if (output_fname) {
    sprintf(plot_fname,"%s%s",CPlot::sOutDir.Data(),output_fname);
    FILE *f=fopen(plot_fname,"w");
    if (!f) {
      std::cout << "failed to create a file <" << plot_fname << ">\n";
      return 0;
    }
    for (unsigned int i=0; i<html_lines.size(); ++i) {
      if ((i==save_pos) && extra_html_lines) {
	fprintf(f,"<hr />\n");
	for (unsigned int j=0; j<extra_html_lines->size(); ++j) {
	  fprintf(f,"%s",(*extra_html_lines)[j].c_str());
	} 
	fprintf(f,"<hr />\n");
      }
      fprintf(f,"%s",html_lines[i].c_str());
    }
    fclose(f);
  }
  return 1;
}

// -------------------------------------------------------

template<class HistoClass_t>
inline
int PrepareComparisonTablesAndPlots(const char *file_tag, const std::vector<HistoClass_t*> &unsmearedMCHV, const std::vector<HistoClass_t*> &unscaledExpHV, const std::vector<HistoClass_t*> &smearedMCHV, const std::vector<HistoClass_t*> &scaledExpHV, const std::vector<std::vector<int>*> &combine, const std::vector<TString> &combination_names, int print_discrepancy_table, int print_descriptions_separately, const char *output_dir, std::vector<std::string> &html_lines) {
  const char *fncname="PrepareComparisonTablesAndPlots";

  std::cout << "entered " << fncname << " with " << unsmearedMCHV.size() << ", " << unscaledExpHV.size() << ", " << smearedMCHV.size() << ", and " << scaledExpHV.size() << " plots\n";

  int massBinCount=DYTools::NumberOfMassBins();
  double *massBins=DYTools::getMassBins();
  const char *massGridFormat=(DYTools::massRangeSet==13) ? " %2.1lf" : " %2.0lf";

  std::vector<HistoClass_t*> unProcAsymHV, procAsymHV;
  if (!GetAsymmetryV(unsmearedMCHV,unscaledExpHV,unProcAsymHV) ||
      !GetAsymmetryV(smearedMCHV,scaledExpHV,procAsymHV)) {
    std::cout << fncname << ": failed with basic data collections\n";
    return 0;
  }
    
  std::vector<HistoClass_t*> unsmearedCombMCHV, unscaledCombExpHV;
  std::vector<HistoClass_t*> smearedCombMCHV, scaledCombExpHV;
  std::vector<HistoClass_t*> unProcCombAsymHV, procCombAsymHV;
  if (!combineEtaEtaHistos(unsmearedMCHV,combine,combination_names,unsmearedCombMCHV) ||
      !combineEtaEtaHistos(unscaledExpHV,combine,combination_names,unscaledCombExpHV) ||
      !combineEtaEtaHistos(smearedMCHV,combine,combination_names,smearedCombMCHV) ||
      !combineEtaEtaHistos(scaledExpHV,combine,combination_names,scaledCombExpHV) ||
      !GetAsymmetryV(unsmearedCombMCHV,unscaledCombExpHV,unProcCombAsymHV) ||
      !GetAsymmetryV(smearedCombMCHV,scaledCombExpHV,procCombAsymHV))
    {
      std::cout << fncname << ": failed with combined data collections\n";
      return 0;
    }

  std::cout << "prepared combination histos (counts): " << unsmearedCombMCHV.size() << ", " << unscaledCombExpHV.size() << ", " << smearedCombMCHV.size() << ", " << scaledCombExpHV.size() << "\n";

  std::vector<HistoClass_t*> unsmearedCombNormMCHV, unscaledCombNormExpHV;
  std::vector<HistoClass_t*> smearedCombNormMCHV, scaledCombNormExpHV;
  std::vector<HistoClass_t*> unProcCombNormAsymHV, procCombNormAsymHV;
  std::vector<double> unsmearedCombNormMCV, unscaledCombNormExpV;
  std::vector<double> smearedCombNormMCV, scaledCombNormExpV;
  
  const int normalize_to_1=0;
  // Normalization to 1 
  if (normalize_to_1) {
    if (!Normalize(unsmearedCombMCHV, unsmearedCombNormMCHV, unsmearedCombNormMCV) ||
	!Normalize(unscaledCombExpHV, unscaledCombNormExpHV, unscaledCombNormExpV) ||
	!Normalize(smearedCombMCHV, smearedCombNormMCHV, smearedCombNormMCV) ||
	!Normalize(scaledCombExpHV, scaledCombNormExpHV, scaledCombNormExpV) ||
	!GetAsymmetryV(unsmearedCombNormMCHV,unscaledCombNormExpHV, unProcCombNormAsymHV) ||
	!GetAsymmetryV(smearedCombNormMCHV,scaledCombNormExpHV, procCombNormAsymHV)) {
      std::cout << fncname << ": failed with normalized data collections (normalize to unity)\n";
      return 0; 
    }
  }
  else {
    // normalization to Exp
    if (!NormalizeTo(unsmearedCombMCHV,unscaledCombExpHV, unsmearedCombNormMCHV, unsmearedCombNormMCV) ||
	!NormalizeTo(smearedCombMCHV,scaledCombExpHV, smearedCombNormMCHV, smearedCombNormMCV)) {
      std::cout << fncname << ": failed with normalized data collections (normalize to Exp)\n";
      return 0;
    }
    unscaledCombNormExpHV=unscaledCombExpHV;
    scaledCombNormExpHV=scaledCombExpHV;
  }
  
  //std::vector<std::string> lines2a,lines2b,lines2c,lines2d;
  //SPrint(unProcCombAsymHV,lines2a,1,1,2);
  //SPrint(unProcCombNormAsymHV,lines2b,1,1,2);
  //SPrint(procCombAsymHV,lines2c,1,1,2);
  //SPrint(procCombNormAsymHV,lines2d,1,1,2);
  ////PrintLines(lines2a);
  ////PrintLines(lines2b);
  ////std::cout << "\n";
  ////PrintLines(lines2c);
  ////PrintLines(lines2d);
  
  // Prepare the container of plots to print
  std::vector<const std::vector<HistoClass_t*>*> HistosV;
  std::vector<std::string> histoSetNames;
  std::vector<std::string> plotTableHeaders;
  std::vector<std::string> legendPlacement;

  //plotTableHeaders.push_back("| Combination name ");
  plotTableHeaders.push_back("| processed exp & MC");
  plotTableHeaders.push_back("| unprocessed exp & MC");
  HistosV.push_back(&smearedCombNormMCHV);
  HistosV.push_back(&scaledCombNormExpHV);
  HistosV.push_back(&unsmearedCombNormMCHV);
  HistosV.push_back(&unscaledCombNormExpHV);
  histoSetNames.push_back("smearedCombNormMC");
  histoSetNames.push_back((normalize_to_1) ? "scaledCombNormExp" : "scaledCombExp");
  histoSetNames.push_back("unsmearedCombNormMC");
  histoSetNames.push_back((normalize_to_1) ? "unscaledCombNormExp" : "unscaledCombExp");
  for (unsigned int i=0; i<4; ++i) legendPlacement.push_back("rt1");
  
  plotTableHeaders.push_back("| change in asymmetry");
  HistosV.push_back(&unProcCombAsymHV); 
  HistosV.push_back(&procCombAsymHV); 
  histoSetNames.push_back("unproc.asym.");
  histoSetNames.push_back("proc.asym.");
  for (unsigned int i=0; i<2; ++i) legendPlacement.push_back("lt");
  

  HistoClass_t emptyAsymH("emptyAsymH","emptyAsymH",massBinCount,massBins);
  std::vector<HistoClass_t*> emptyAsymHV;
  emptyAsymH.GetXaxis()->SetTitle("mass");
  emptyAsymH.GetYaxis()->SetTitle("asymmmetry");
  for (unsigned int i=0; i<procCombAsymHV.size(); ++i) {
    emptyAsymHV.push_back(NULL);
  }

  plotTableHeaders.push_back("| new asymmetry");
  HistosV.push_back(&emptyAsymHV);
  HistosV.push_back(&procCombAsymHV); 
  histoSetNames.push_back("");
  histoSetNames.push_back("proc.asym.");
  for (unsigned int i=0; i<2; ++i) legendPlacement.push_back("lt");
  
  plotTableHeaders.push_back("| unproc.vs.proc. MC");
  plotTableHeaders.push_back("| unproc.vs.proc. Exp");
  HistosV.push_back(&unsmearedCombMCHV);
  HistosV.push_back(&smearedCombMCHV);
  HistosV.push_back(&unscaledCombExpHV);
  HistosV.push_back(&scaledCombExpHV);
  histoSetNames.push_back("unsmeared MC");
  histoSetNames.push_back("smeared MC");
  histoSetNames.push_back("unscaled Exp");
  histoSetNames.push_back("scaled Exp");
  for (unsigned int i=0; i<4; ++i) legendPlacement.push_back("rt1");


  // goodness of fit
  char buf[100];
  std::vector<double> ymaxProcV,ymaxUnProcV;
  std::vector<double> sumY2ProcV,sumY2UnProcV;
  std::vector<double> chi2ProcV,chi2UnProcV;
  const unsigned int csize=combine.size();

  ymaxProcV.reserve(csize); ymaxUnProcV.reserve(csize);
  sumY2ProcV.reserve(csize); sumY2UnProcV.reserve(csize);
  chi2ProcV.reserve(csize); chi2UnProcV.reserve(csize);
  for (unsigned int ci=0; ci<combine.size(); ++ci) {
    double ymaxA,sumy2A;
    double ymaxB,sumy2B;
    CalculateMaxValues(*procCombAsymHV[ci],ymaxA,sumy2A);
    CalculateMaxValues(*unProcCombAsymHV[ci],ymaxB,sumy2B);
    ymaxProcV.push_back(ymaxA); sumY2ProcV.push_back(sumy2A);
    ymaxUnProcV.push_back(ymaxB); sumY2UnProcV.push_back(sumy2B);
    chi2ProcV.push_back(smearedCombNormMCHV[ci]->Chi2Test(scaledCombNormExpHV[ci],"WW CHI2"));
    chi2UnProcV.push_back(unsmearedCombNormMCHV[ci]->Chi2Test(unscaledCombNormExpHV[ci],"WW CHI2"));
  }
  html_lines.push_back("<hr />\n<h3>Goodness of fit info</h3><br>\n");
  sprintf(buf,"Mass grid [%d]: ",massBinCount);
  html_lines.push_back(buf);
  for (int i=0; i<=massBinCount; ) {
    int pos=0;
    for (int ti=0; (ti<5) && (i<=massBinCount); ++ti, ++i) {
      sprintf(buf+pos,massGridFormat,massBins[i]);
      pos=strlen(buf);
    }
    html_lines.push_back(buf);
  }
  html_lines.push_back("\n<br>");
  html_lines.push_back("<table border=""1"" cellspacing=""2"" width=""100%"">\n");
  html_lines.push_back("<tr><td width=""15\%""Collection</td><td>&chi;<sup>2</sup> of processed</td><td>&chi;<sup>2</sup> of unprocessed</td><td>max <i>dy</i><sub>i</sub> of processed</td><td>max <i>dy</i><sub>i</sub> of unprocessed</td><td>&Sigma;(<i>y</i><sub>i</sub>)<sup>2</sup> of processed</td><td>&Sigma;(<i>y</i><sub>i</sub>)<sup>2</sup> of unprocessed</td></tr>\n");
  for (unsigned int ci=0; ci<combine.size(); ++ci) {
    sprintf(buf,"<tr><td width=""15%%"">%s</td><td>%12.4g</td><td>%12.4g</td><td>%12.4g</td><td>%12.4g</td><td>%12.4g</td><td>%12.4g</td></tr>\n",
	    combination_names[ci].Data(), 
	    chi2ProcV[ci],chi2UnProcV[ci],
	    ymaxProcV[ci], ymaxUnProcV[ci], sumY2ProcV[ci], sumY2UnProcV[ci]);
    html_lines.push_back(buf);
  }
  html_lines.push_back("</table>\n<small><i>Here</i> &chi;<sup>2</sup> obtained from TH1::Chi2Test(TH1*,opt=""WW CHI2"")</small><br>\n<hr />\n");


  // prepare tables
  std::vector<std::vector<std::string>*> tables;
  if (print_discrepancy_table) {
    tables.reserve(combine.size());
    
    const int print_range=1; const int print_no_range=0;
    const int print_err=1; const int print_no_err=0;
    const int print_acc2=2; const int print_acc1=1;
    const int print_name=1;
    for (unsigned int ci=0; ci<combine.size(); ++ci) {
      std::vector<std::string> infoLines;
      std::vector<std::string> *htmlLines=new std::vector<std::string>();
      if (normalize_to_1) htmlLines->push_back("<br><small>MC normalized to 1</small></br>\n");
      else htmlLines->push_back("<br><small>MC normalized to experimental data area (hExp->Integral())</small></br>\n");

      SPrint(*smearedCombNormMCHV[ci],infoLines,print_range,print_no_err,print_acc1,print_name,"smeared MC");
      SPrint(*scaledCombNormExpHV[ci],infoLines,print_no_range,print_no_err,print_acc1,print_name,"scaled Data");
      SPrint(*procCombAsymHV[ci],infoLines,print_no_range,print_err,print_acc2,print_name,"modified asymmetry");
      SPrint(*unsmearedCombNormMCHV[ci],infoLines,print_no_range,print_no_err,print_acc1,print_name,"original MC");
      SPrint(*unscaledCombNormExpHV[ci],infoLines,print_no_range,print_no_err,print_acc1,print_name,"original Data");
      SPrint(*unProcCombAsymHV[ci],infoLines,print_no_range,print_err,print_acc2,print_name,"original asymmetry.");

      if (1) {
	sprintf(buf,"| max abs(y<sub>i</sub>) | | | %12.4g | | | %12.4g |\n",ymaxProcV[ci],ymaxUnProcV[ci]);
	infoLines.push_back(buf);
	sprintf(buf,"| sum (y<sub>i</sub>)<sup>2</sup> | | | %12.4g | | | %12.4g |\n",sumY2ProcV[ci],sumY2UnProcV[ci]);
	infoLines.push_back(buf);
      }

      PrintLinesHTML(*htmlLines,infoLines,NULL,1,1,0,0); // make a table, not a doc
      
      //std::cout << "\n\n"; PrintLines(std::cout,infoLines); std::cout << "\n\n";      
      tables.push_back(htmlLines);
    }
  }

  sprintf(buf,"comparison_data_%s.html",file_tag);
  std::string output_html_name=buf;
  std::vector<std::string> html_comparison_lines; // dummy
  PrintHTMLHistoLines_OverlayedDuos(output_dir,file_tag,HistosV,histoSetNames,plotTableHeaders,combination_names,tables,print_descriptions_separately,"Comparison tables and plots",html_comparison_lines,1,1,&html_lines,output_html_name.c_str(),&legendPlacement);


  std::string info="<hr /><a href="""; info+=output_html_name; info+=""">Comparison tables and plots</a><hr />\n";
  html_lines.push_back(info);
  return 1;
}

// -------------------------------------------------------

void fillCases(std::vector<std::string> &etaRanges, std::vector<std::string> &models) {
  etaRanges.clear(); models.clear();
  etaRanges.push_back("2bin"); models.push_back("fit model gauss");
  etaRanges.push_back("2binNegs"); models.push_back("fit model gauss");
  etaRanges.push_back("3EB3EE"); models.push_back("fit model gauss");
  etaRanges.push_back("3EB3EENegs"); models.push_back("fit model gauss");
  etaRanges.push_back("4bin"); models.push_back("fit model gauss");
  etaRanges.push_back("4binNegs"); models.push_back("fit model gauss");
  etaRanges.push_back("4EB3EE"); models.push_back("fit model gauss");
  etaRanges.push_back("4EB3EENegs"); models.push_back("fit model gauss");
  etaRanges.push_back("5bin"); models.push_back("fit model gauss");
  etaRanges.push_back("5binNegs"); models.push_back("fit model gauss");
  etaRanges.push_back("6bin"); models.push_back("fit model gauss");
  etaRanges.push_back("6binNegs"); models.push_back("fit model gauss");
  etaRanges.push_back("6bin"); models.push_back("fit model voigtian");
  etaRanges.push_back("6binNegs"); models.push_back("fit model voigtian");
  etaRanges.push_back("6bin"); models.push_back("fit model breitwigner");
  etaRanges.push_back("6binNegs"); models.push_back("fit model breitwigner");
  return;
}

// -------------------------------------------------------

#endif
