{
//=========Macro generated from canvas: c1/Fake rates v #{eta}
//=========  (Tue Jul  3 18:20:42 2012) by ROOT version5.27/06
   TCanvas *c1 = new TCanvas("c1", "Fake rates v #{eta}",1,58,700,500);
   c1->SetHighLightColor(2);
   c1->Range(-75.00001,-0.02416073,341.6667,0.1369108);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetLeftMargin(0.18);
   c1->SetBottomMargin(0.15);
   c1->SetFrameBorderMode(0);
   c1->SetFrameBorderMode(0);
   Double_t xAxis1[16] = {0, 5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 90, 150, 200, 300}; 
   
   TH1F *fRate_reg = new TH1F("fRate_reg","",15, xAxis1);
   fRate_reg->SetBinContent(3,0.05544148);
   fRate_reg->SetBinContent(4,0.07267442);
   fRate_reg->SetBinContent(5,0.05136986);
   fRate_reg->SetBinContent(6,0.09814323);
   fRate_reg->SetBinContent(7,0.07929717);
   fRate_reg->SetBinContent(8,0.08181045);
   fRate_reg->SetBinContent(9,0.07429676);
   fRate_reg->SetBinContent(10,0.07295229);
   fRate_reg->SetBinContent(11,0.07320368);
   fRate_reg->SetBinContent(12,0.06964575);
   fRate_reg->SetBinContent(13,0.07234851);
   fRate_reg->SetBinContent(14,0.07405092);
   fRate_reg->SetBinContent(15,0.07699784);
   fRate_reg->SetBinError(3,0.007750951);
   fRate_reg->SetBinError(4,0.01505378);
   fRate_reg->SetBinError(5,0.01360005);
   fRate_reg->SetBinError(6,0.01690787);
   fRate_reg->SetBinError(7,0.003616134);
   fRate_reg->SetBinError(8,0.004336647);
   fRate_reg->SetBinError(9,0.004139939);
   fRate_reg->SetBinError(10,0.00357018);
   fRate_reg->SetBinError(11,0.004606693);
   fRate_reg->SetBinError(12,0.001496737);
   fRate_reg->SetBinError(13,0.001108618);
   fRate_reg->SetBinError(14,0.00166207);
   fRate_reg->SetBinError(15,0.002992548);
   fRate_reg->SetEntries(1059.645);
   
   TPaveStats *ptstats = new TPaveStats(0.78,0.835,0.98,0.995,"brNDC");
   ptstats->SetName("stats");
   ptstats->SetBorderSize(1);
   ptstats->SetTextAlign(12);
   TText *text = ptstats->AddText("fRate_reg");
   text->SetTextSize(0.0368);
   text = ptstats->AddText("Entries = 1060   ");
   text = ptstats->AddText("Mean  =  73.55");
   text = ptstats->AddText("RMS   =  68.18");
   ptstats->SetOptStat(1111);
   ptstats->SetOptFit(0);
   ptstats->Draw();
   fRate_reg->GetListOfFunctions()->Add(ptstats);
   ptstats->SetParent(fRate_reg->GetListOfFunctions());

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#00ff00");
   fRate_reg->SetLineColor(ci);

   ci = TColor::GetColor("#00ff00");
   fRate_reg->SetMarkerColor(ci);
   fRate_reg->GetXaxis()->SetTitle("electron E_{T} GeV/c");
   fRate_reg->GetYaxis()->SetTitle("Fake Rate");
   fRate_reg->Draw("");
   Double_t xAxis2[16] = {0, 5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 90, 150, 200, 300}; 
   
   TH1F *fRate_reg2 = new TH1F("fRate_reg2","",15, xAxis2);
   fRate_reg2->SetBinContent(3,0.01030928);
   fRate_reg2->SetBinContent(4,0.03278688);
   fRate_reg2->SetBinContent(5,0.1458333);
   fRate_reg2->SetBinContent(6,0.1621622);
   fRate_reg2->SetBinContent(7,0.08378217);
   fRate_reg2->SetBinContent(8,0.08285862);
   fRate_reg2->SetBinContent(9,0.07792836);
   fRate_reg2->SetBinContent(10,0.07163323);
   fRate_reg2->SetBinContent(11,0.06077694);
   fRate_reg2->SetBinContent(12,0.06384304);
   fRate_reg2->SetBinContent(13,0.06695604);
   fRate_reg2->SetBinContent(14,0.07544936);
   fRate_reg2->SetBinContent(15,0.1092362);
   fRate_reg2->SetBinError(3,0.01036228);
   fRate_reg2->SetBinError(4,0.02356083);
   fRate_reg2->SetBinError(5,0.05900224);
   fRate_reg2->SetBinError(6,0.05046523);
   fRate_reg2->SetBinError(7,0.007371546);
   fRate_reg2->SetBinError(8,0.006816534);
   fRate_reg2->SetBinError(9,0.00637643);
   fRate_reg2->SetBinError(10,0.006054692);
   fRate_reg2->SetBinError(11,0.006355724);
   fRate_reg2->SetBinError(12,0.002655307);
   fRate_reg2->SetBinError(13,0.001816264);
   fRate_reg2->SetBinError(14,0.002955232);
   fRate_reg2->SetBinError(15,0.007335187);
   fRate_reg2->SetEntries(155.9759);

   ci = TColor::GetColor("#ff0000");
   fRate_reg2->SetLineColor(ci);

   ci = TColor::GetColor("#ff0000");
   fRate_reg2->SetMarkerColor(ci);
   fRate_reg2->Draw("same");
   Double_t xAxis3[16] = {0, 5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 90, 150, 200, 300}; 
   
   TH1F *fRate_reg3 = new TH1F("fRate_reg3","",15, xAxis3);
   fRate_reg3->SetBinContent(3,0.01932367);
   fRate_reg3->SetBinContent(4,0.00729927);
   fRate_reg3->SetBinContent(5,0.06603774);
   fRate_reg3->SetBinContent(6,0.02247191);
   fRate_reg3->SetBinContent(7,0.04685573);
   fRate_reg3->SetBinContent(8,0.07098121);
   fRate_reg3->SetBinContent(9,0.05681818);
   fRate_reg3->SetBinContent(10,0.06757253);
   fRate_reg3->SetBinContent(11,0.05798969);
   fRate_reg3->SetBinContent(12,0.0612772);
   fRate_reg3->SetBinContent(13,0.07365035);
   fRate_reg3->SetBinContent(14,0.08447327);
   fRate_reg3->SetBinContent(15,0.09912927);
   fRate_reg3->SetBinError(3,0.00975474);
   fRate_reg3->SetBinError(4,0.007325861);
   fRate_reg3->SetBinError(5,0.02577089);
   fRate_reg3->SetBinError(6,0.01606759);
   fRate_reg3->SetBinError(7,0.004490076);
   fRate_reg3->SetBinError(8,0.004761523);
   fRate_reg3->SetBinError(9,0.00418283);
   fRate_reg3->SetBinError(10,0.005147069);
   fRate_reg3->SetBinError(11,0.005133631);
   fRate_reg3->SetBinError(12,0.002255965);
   fRate_reg3->SetBinError(13,0.001794267);
   fRate_reg3->SetBinError(14,0.003594309);
   fRate_reg3->SetBinError(15,0.008542706);
   fRate_reg3->SetEntries(421.2557);

   ci = TColor::GetColor("#0000ff");
   fRate_reg3->SetLineColor(ci);

   ci = TColor::GetColor("#0000ff");
   fRate_reg3->SetMarkerColor(ci);
   fRate_reg3->Draw("same");
   
   TLegend *leg = new TLegend(0.2,0.68,0.45,0.85,NULL,"brNDC");
   leg->SetBorderSize(1);
   leg->SetTextFont(62);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(19);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("fRate_reg","#eta < 1.4442","lp");

   ci = TColor::GetColor("#00ff00");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#00ff00");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(1);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("fRate_reg2","1.566< #eta < 2.0","lp");

   ci = TColor::GetColor("#ff0000");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#ff0000");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(1);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("fRate_reg3","2.0 < #eta < 2.5","lp");

   ci = TColor::GetColor("#0000ff");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#0000ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(1);
   entry->SetMarkerSize(1);
   leg->Draw();
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
