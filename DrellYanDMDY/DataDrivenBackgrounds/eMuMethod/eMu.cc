#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TText.h>
#include <THStack.h>
#include <math.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <algorithm>

#include "eMu.hh"
#include "DYTools.hh"

using std::cout;
using std::ostringstream;
using std::setprecision;
using std::make_pair;
using std::for_each;

/*
using DYTools::_massBinLimits2011;
using DYTools::_massBinLimits2D;
using DYTools::_nMassBins2D;
using DYTools::yRangeMin;
using DYTools::yRangeMax;
using DYTools::_nYBinsMax2D;
*/

eMu::eMu():doDMDY(false),doPUreWeight(false),saveToRootFile(false), verbose(false), reWeight(1.00){
  //set default parameters
  //These params can be overidden via python interface
  emuNtupleDir = "../../root_files/selected_events/DY_m10+pr+a05+o03+pr_4680pb/ntuples_emu";
  eeNtupleDir = "../../root_files/selected_events/DY_m10+pr+a05+o03+pr_4680pb/ntuples";
  //emuNtupleDir = "EMU_TEST";
  //eeNtupleDir = "EE_TEST";
  subDir = "/";
  filePostfix = "_select.root";
  outFileName = "true2eBkgDataPoints_tmp.root";
  xmax = 600;
  xmin = 0;
  nBins = 30;
  rap_Bins = _nYBinsMax2D;
  rap_Min = yRangeMin;
  rap_Max = yRangeMax;
  lumiVal = "4.7";
  dataMap.insert(make_pair("data","Data"));
  eebkgMap.insert(make_pair("ttbar","t#bar{t}"));
  eebkgMap.insert(make_pair("ztt","Z#rightarrow#tau#tau"));
  eebkgMap.insert(make_pair("ww","WW"));  
  eebkgMap.insert(make_pair("wtop","tW"));
  eebkgMap.insert(make_pair("wtopbar","#bar{t}W"));
  eebkgMap.insert(make_pair("wz","WZ"));
  eebkgMap.insert(make_pair("zz","ZZ"));
  emubkgMap.insert(make_pair("wjets","WX#nu (X=e/#mu/#tau)"));

  d_array_size = dataMap.size();
  emu_array_size = eebkgMap.size();
  emu_bkgarray_size = emubkgMap.size();
  
  //true2eBackgroundFromData(6,25); //This is hardwired, need to change this
  //true2eBackgroundFromDataError(6,25);
  //true2eBackgroundFromDataErrorSyst(6,25); 
}

//Need to be able to adapt this to constant binning. CONFIG FILE!
//const double _massBinLimits2011[] = 
  //  {15,20,30,40,50,60,76,86,96,106,120,150,200,600};
  //{15,20,25,30,35,40,45,50,55,60,64,68,72,76,81,86,91,96,101,106,110,115,120,126,133,141,150,160,171,185,200,220,243,273,320,380,440,510,600,1000,1500}; 
  //{20,40,60,80,100,120,140,160,180,200,220,240,260,280,300,320,340,360,380,400,420,440,460,480,500,520,540,560,580,600,620,640,660,680,700,720,740,760,780,800,820,840,860,880,900,920,940,960,980,1000,1020,1040,1060,1080,1100,1120,1140,1160,1180,1200};
//{10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,360,370,380,390,400,410,420,430,440,450,460,470,480,490,500,510,520,530,540,550,560,570,580,590,600,610,620,630,640,650,660,670,680,690,700,710,720,730,740,750,760,770,780,790,800,810,820,830,840,850,860,870,880,890,900,910,920,930,940,950,960,970,980,990,1000,1010,1020,1030,1040,1050,1060,1070,1080,1090,1100,1110,1120,1130,1140,1150,1160,1170,1180,1190,1200,1210,1220,1230,1240,1250,1260,1270,1280,1290,1300,1310,1320,1330,1340,1350,1360,1370,1380,1390,1400,1410,1420,1430,1440,1450,1460,1470,1480,1490,1500};

const int nBinsV = sizeof(_massBinLimits2011)/sizeof(double);
const int numMassBins = nBinsV - 1;

//function definitions
void sethistoStyle(TH1F *inputHisto);
void setEEhistoStyle(TH1F *inputHisto);
void setEMUhistoStyle(TH1F *inputHisto);
void printBinContents(vector<vector<shared_ptr<TH1F> >* >& vvHist);
void printBinContents(TH1F* myhist);
void printBinContents2(TH1F* inHisto, vector<vector<shared_ptr<TH1F> >* >& vvHist);


//unsigned int d_array_size = sizeof(d_array)/sizeof(string);
//unsigned int emu_array_size = sizeof(emu_array)/sizeof(string);
//unsigned int emu_bkgarray_size = sizeof(emu_bkgarray)/sizeof(string);

int eMu::run(){
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetHistMinimumZero(kTRUE) ;

  string d_array[dataMap.size()];// = {"data"};
  string emu_array[eebkgMap.size()];// = {"ttbar","ww","ztt"}; //need to add wtop
  string emu_bkgarray[emubkgMap.size()];// = {"wjets","wz"};
  //fill emu_array
  for_each(dataMap.begin(), dataMap.end(), Fill_array_from_map_functor(d_array));
  for_each(eebkgMap.begin(), eebkgMap.end(), Fill_array_from_map_functor(emu_array));
  for_each(emubkgMap.begin(), emubkgMap.end(), Fill_array_from_map_functor(emu_bkgarray));

  if (verbose) cout << numMassBins << "\t" << nBinsV << " bin numbers \n";

  vector<string>* emu_modesv = new vector<string>;
  //add in the names from the array,data,signa; MC and background MC
  //Why not add them as one array? Because I want to keep track of the individual components
  emu_modesv->insert(emu_modesv->end(), d_array, d_array + d_array_size);
  emu_modesv->insert(emu_modesv->end(), emu_array, emu_array + emu_array_size );
  emu_modesv->insert(emu_modesv->end(), emu_bkgarray, emu_bkgarray + emu_bkgarray_size);

  vector<string>* ee_modesv = new  vector<string>; //add data and background 
  ee_modesv->insert(ee_modesv->end(), d_array, d_array + d_array_size);
  ee_modesv->insert(ee_modesv->end(), emu_array, emu_array + emu_array_size );
  string s_path_array[] = {emuNtupleDir,eeNtupleDir};//directory paths. 
  vector<string> s_pathv(s_path_array,s_path_array+2);// initialise vector with string array

  THStack eMuMassTot("eMuMassTot","");
  THStack eMuVMassTot("eMuVMassTot","");
  THStack eeMassTot("eeMassTot","");
  THStack eeVMassTot("eeVMassTot","");

  map<double, shared_ptr<THStack> > eMuDMDYTot;
  map<double, shared_ptr<THStack> > eeDMDYTot;

  //data vector histos
  vector<shared_ptr<TH1F> > datahistv; //variable bin histo
  vector<shared_ptr<TH1F> > datahMassv; //fixed bin histo
  vector<shared_ptr<TH1F> > sumMassv;
  vector<shared_ptr<TH1F> > sumVMassv;
  vector<vector<shared_ptr<TH1F> >* > hist2Dv;
  vector<vector<shared_ptr<TH1F> >* > hMass2Dv;
  vector<vector<shared_ptr<TH1F> >* > stathist2Dv;
  vector<vector<shared_ptr<TH1F> >* > stathMass2Dv;

  //For 2D analysis

  //Fill massBins2D from _massBinLimits2D but ignore overflow and underflow

  const int numBinsWithOverFlowAndUnderFlow = _nMassBins2D+1;
  const int numBinsOfInterest =  _nMassBins2D-1;
  double massBins2D[numBinsWithOverFlowAndUnderFlow]; //upper edge values of mass bins

  if (doDMDY){
    for (int i  = 1; i < numBinsWithOverFlowAndUnderFlow; ++i){//start after underflow
      massBins2D[i-1] = _massBinLimits2D[i];
    }
    massBins2D[_nMassBins2D] = 8000; //Add overflow bin at the end to catch all
  }
  
  const int num_massBins2D = sizeof(massBins2D)/sizeof(double); //This should be equivalent to (_nMassBins2D + 1)

  TMatrixT<double> true2eBackgroundFromData(numBinsOfInterest,rap_Bins);
  TMatrixT<double> true2eBackgroundFromDataError(numBinsOfInterest,rap_Bins);
  TMatrixT<double> true2eBackgroundFromDataErrorSyst(numBinsOfInterest,rap_Bins);

  map<double, vector<shared_ptr<TH1F> > > data_dmdyHistv;
  map<double, vector<shared_ptr<TH1F> > > sum_dmdyHistv;
  map<double, vector<vector<shared_ptr<TH1F> >* > > mc_dmdyHistv;
  map<double, vector<vector<shared_ptr<TH1F> >* > > stats_mc_dmdyHistv;

  if (doDMDY){
    //loop over over dmdy mass bins
    for (int i = 0; i < num_massBins2D; ++i){
      double massB = massBins2D[i];
      ostringstream upperMass;
      upperMass << massB;
      data_dmdyHistv.insert(make_pair(massB,vector<shared_ptr<TH1F> >()));
      sum_dmdyHistv.insert(make_pair(massB,vector<shared_ptr<TH1F> >()));
      mc_dmdyHistv.insert(make_pair(massB,vector<vector<shared_ptr<TH1F> >*>()));
      stats_mc_dmdyHistv.insert(make_pair(massB,vector<vector<shared_ptr<TH1F> >*>()));
      TString hist_name("eMuDMDYTot"+upperMass.str());
      TString hist_name2("eMuDMDYTot"+upperMass.str());
      eMuDMDYTot.insert(make_pair(massB,shared_ptr<THStack>(new THStack(hist_name,"" ))) );
      eeDMDYTot.insert(make_pair(massB,shared_ptr<THStack>(new THStack(hist_name2,""))));
    }
  }

  vector<TCanvas*> canv;  

  typedef vector<string>::iterator sIter;

  for(sIter iter_path = s_pathv.begin(); iter_path !=s_pathv.end(); ++iter_path){

    int counter(0);
    int colour(2);//start with red  

    if (doDMDY){
      //sort out binning in case of dmdy analysis. should but in dmdy condition
      for (int i = 0; i < num_massBins2D; ++i){
	double massB = massBins2D[i];
	//for each map index which corresponds to a vector of vectors containing TH1F pointers, instantiate!
	mc_dmdyHistv[massB].push_back(new vector<shared_ptr<TH1F> >());
	stats_mc_dmdyHistv[massB].push_back(new vector<shared_ptr<TH1F> >());
	//for each map index which corresponds to a vectors containing TH1F pointers, instantiate!
	ostringstream upperMass;
	upperMass << massB;
	TString histoName(*iter_path+"upperMass"+upperMass.str());
	//need to change number of bins for last mass range
	if (i == 6) rap_Bins = 10;
	data_dmdyHistv[massB].push_back(shared_ptr<TH1F>(new TH1F(histoName,"",rap_Bins, rap_Min, rap_Max)));
	sum_dmdyHistv[massB].push_back(shared_ptr<TH1F>(new TH1F(histoName+"sum","",rap_Bins, rap_Min, rap_Max)));
	if (i == 6) rap_Bins = _nYBinsMax2D;
	// need to use correct binning
	//need to change binning in last vector
      }
    }

    vector<shared_ptr<TH1F> >* histv = new vector<shared_ptr<TH1F> >();
    vector<shared_ptr<TH1F> >* hMassv = new vector<shared_ptr<TH1F> >();
    vector<shared_ptr<TH1F> >* stathistv = new vector<shared_ptr<TH1F> >();
    vector<shared_ptr<TH1F> >* stathMassv = new vector<shared_ptr<TH1F> >();
    hist2Dv.push_back(histv);
    hMass2Dv.push_back(hMassv);
    stathist2Dv.push_back(stathistv);
    stathMass2Dv.push_back(stathMassv);

    datahistv.push_back(shared_ptr<TH1F>(new TH1F("datahistv","",numMassBins,_massBinLimits2011)));
    datahMassv.push_back(shared_ptr<TH1F>(new TH1F("datahMassv","",nBins,xmin,xmax)));

    sumVMassv.push_back(shared_ptr<TH1F>(new TH1F("sumVMassv","",numMassBins,_massBinLimits2011)));
    sumMassv.push_back(shared_ptr<TH1F>(new TH1F("sumMassv","",nBins,xmin,xmax)));

    //select the list of strings to loop over;  emu_modesv or ee_modesv?
    vector<string>* sv;
    sv = (s_pathv.at(0) == *iter_path) ? emu_modesv : ee_modesv; 

    for(sIter iter = sv->begin(); iter != sv->end(); ++iter){
       
      bool isData(false);
      //check if we are looking at data
      if (*iter == "data"){//need to remove this hardwiring and specify the map name instead
	isData = true;
      } else {
	isData = false; 
      }
 
      //should use smart pointers everywhere
      // memory leaked for 2nd iteration and above
      //not a big deal right now. Fix later! 
      vector<TFile*> filev;
      vector<TTree*> treev;

      string fname = *iter_path + subDir + *iter + filePostfix;//make the path configurable. ie. be able to choose "/ntuples/" and "_select.root"
      TString TS_fname(fname);
      //should use smart pointers
      filev.push_back(new TFile(TS_fname));
      treev.push_back((TTree*) filev.back()->Get("Events")); //make "Events" configurable rather than hardwiring it   
             
      if (!isData) {
        ostringstream scounter;//hold number as a string
        scounter << counter;
        ++counter;
        
        string hname = string("histo") + scounter.str();
        TString TS_hname(hname);
        hist2Dv.back()->push_back(shared_ptr<TH1F>(new TH1F(TS_hname,"",numMassBins,_massBinLimits2011)));
        hMass2Dv.back()->push_back(shared_ptr<TH1F>(new TH1F(TS_hname+TString("a"),"",nBins,xmin,xmax)));
        stathist2Dv.back()->push_back(shared_ptr<TH1F>(new TH1F(TS_hname+TString("stat"),"",numMassBins,_massBinLimits2011)));
        stathMass2Dv.back()->push_back(shared_ptr<TH1F>(new TH1F(TS_hname+TString("stat_a"),"",nBins,xmin,xmax)));

	if (doDMDY){
          for (int i = 0; i < num_massBins2D; ++i){
	    double massB = massBins2D[i];
	    ostringstream upperMass;
	    upperMass << massB;
	    //need to change number of bins for last mass range
	    if (i == 6) rap_Bins = 10;
	    TString h_name(TS_hname+upperMass.str());
	    mc_dmdyHistv[massB].back()->push_back(shared_ptr<TH1F>(new TH1F(h_name,"",rap_Bins, rap_Min, rap_Max)));
	    stats_mc_dmdyHistv[massB].back()->push_back(shared_ptr<TH1F>(new TH1F(h_name+"_stats","",rap_Bins, rap_Min, rap_Max)));
	    if (i == 6) rap_Bins = _nYBinsMax2D;
	  }
	}

      } 
   
      float mass, weight, pu_weight(1.00),rapidity;
      //double reWeight(1.00);//rather than running the MC everytime a new weight is required, just use the reweight variable
      
      treev.back()->SetBranchAddress("mass", &mass);
      treev.back()->SetBranchAddress("weight", &weight);
      if (doPUreWeight){
        treev.back()->SetBranchAddress("pu_weight", &(pu_weight));
      }
      if (doDMDY){
	treev.back()->SetBranchAddress("rapidity", &(rapidity));
      }
      //number of entries in ntuple
      const unsigned int numEntries = treev.back()->GetEntries(); 
    
      //loop over entries in ntuple
      for(unsigned int i =0; i < numEntries; ++i){//loop     
        treev.back()->GetEntry(i);
        //fill histo with events
	if (isData) {
	  datahMassv.back()->Fill(mass,weight);
	  datahistv.back()->Fill(mass,weight);
          if (doDMDY){
	    (data_dmdyHistv.upper_bound(mass))->second.back()->Fill(fabs(rapidity),weight);
	  }
	} else {
	  
	  hist2Dv.back()->back()->Fill(mass,weight*reWeight*pu_weight);//only rewight MC
	  hMass2Dv.back()->back()->Fill(mass,weight*reWeight*pu_weight);//only rewight MC
          stathist2Dv.back()->back()->Fill(mass);
	  stathMass2Dv.back()->back()->Fill(mass);

          if (doDMDY){

            mc_dmdyHistv.upper_bound(mass)->second.back()->back()->Fill(fabs(rapidity),weight*reWeight*pu_weight);
	    stats_mc_dmdyHistv.upper_bound(mass)->second.back()->back()->Fill(fabs(rapidity),weight*reWeight*pu_weight);
	  }

	}
      }//end of loop over entries

      //if we are looking at a datafile we want to jump to the next file 
      // we only want to fill datahMassv and datahistv in the case of data
      if (isData) continue;

      //sum histo
      sumVMassv.back()->Add(hist2Dv.back()->back().get());
      sumMassv.back()->Add(hMass2Dv.back()->back().get());
    
      //sum histos
      if (doDMDY){
	//loop over over dmdy mass bins
	for (int i = 0; i < num_massBins2D; ++i){
	  double massB = massBins2D[i];
	  sum_dmdyHistv[massB].back()->Add(mc_dmdyHistv[massB].back()->back().get());
	  //set colours of histos
	  mc_dmdyHistv[massB].back()->back()->SetFillColor(colour);

	  //stack histos
	  if (*iter_path == s_pathv.at(0)){ //Emu case only
	    eMuDMDYTot[massB]->Add(mc_dmdyHistv[massB].back()->back().get());
	  }
	  if (*iter_path == s_pathv.at(1)){ //EE case only
	    eeDMDYTot[massB]->Add(mc_dmdyHistv[massB].back()->back().get());
	  }
	}
      }
    
      //Set colours of histos    
      hist2Dv.back()->back()->SetFillColor(colour);
      hMass2Dv.back()->back()->SetFillColor(colour);
      ++colour;//move to next colour         
      
      //Stack histos
      if (*iter_path == s_pathv.at(0)){ //Emu case only
        eMuVMassTot.Add(hist2Dv.back()->back().get());
        eMuMassTot.Add(hMass2Dv.back()->back().get());
      }

      if (*iter_path == s_pathv.at(1)){ //EE case only
        eeVMassTot.Add(hist2Dv.back()->back().get());
        eeMassTot.Add(hMass2Dv.back()->back().get());
      }
     
    }//end of loop over channel names
  }//end of loop over directories 

  TCanvas* emuCan = new TCanvas("emuCan","emu distro");
  emuCan->cd();
  eMuVMassTot.Draw("");
  
  TCanvas* eeCan = new TCanvas("eeCan","ee distro");
  eeCan->cd();
  eeVMassTot.Draw("");

  //fill vector of labels to use in plot
  vector<const char*> entryLabels;
  for_each(eebkgMap.begin(), eebkgMap.end(), Fill_vector_from_map_functor(entryLabels));
  for_each(emubkgMap.begin(), emubkgMap.end(), Fill_vector_from_map_functor(entryLabels));
  
  //loop over each mass bin and plot
  if (doDMDY){

    //plot 2D histos
    TCanvas* rapCan = new TCanvas("rapCan","emu Rapidity distributions");
    rapCan->cd();

    TPad *rapPad1 = new TPad("rapPad1","rapPad1",0.1,0.1,0.95,0.995);
    rapPad1->Divide(2,3);
    rapPad1->Draw();
    rapPad1->SetLogy();
    //rapPpad1->SetLogx();
    rapPad1->cd();

    TString massBinLabels[] = {"0 to 20 GeV", "20 to 30 GeV", "30 to 45 GeV","45 to 60 GeV",
			       "60 to 120 GeV", "120 to 200 GeV", "200 to 1500 GeV", "1500 to 8000 GeV"};

    shared_ptr<TLegend> rapLeg[num_massBins2D];
    for (int i = 1; i < num_massBins2D-1; ++i){//skip plotting the 0->20 and 1500->8000 region
      double massB = massBins2D[i];
      //change section in pad
      rapPad1->cd(i);
      data_dmdyHistv[massB].at(0)->GetXaxis()->SetTitleSize(0.08);
      data_dmdyHistv[massB].at(0)->GetXaxis()->SetLabelSize(0.07);
      data_dmdyHistv[massB].at(0)->SetXTitle("rapidity"); 
      data_dmdyHistv[massB].at(0)->SetMaximum(data_dmdyHistv[massB].at(0)->GetMaximum()*1.3);
      data_dmdyHistv[massB].at(0)->Draw("e ");
      eMuDMDYTot[massB]->Draw("same");
      data_dmdyHistv[massB].at(0)->Draw("e same");

      rapLeg[i] = shared_ptr<TLegend> (new TLegend(0.6 ,0.72,0.85,0.85));
      rapLeg[i]->SetFillColor(3);
      rapLeg[i]->SetTextSize(0.07);
      rapLeg[i]->SetHeader(massBinLabels[i]);
      if (i == 6){ //add key to the last plot
	rapLeg[i] =  shared_ptr<TLegend> (new TLegend(0.55 ,0.35,0.85,0.85));//previous memory allocation should have been freed by smart pointer
        rapLeg[i]->SetFillColor(3);
	rapLeg[i]->SetTextSize(0.07);
	rapLeg[i]->SetHeader(massBinLabels[i]);       
	for (unsigned int j = 0; j < hist2Dv.at(0)->size(); ++j){
	  if (j < entryLabels.size()){//make sure there is an label entry for all histos
	    rapLeg[i]->AddEntry(hist2Dv.at(0)->at(j).get(), entryLabels.at(j),"f");;
	  }
	}
      }
      rapLeg[i]->Draw();
    }
  
    //rapCan->Write();
    rapCan->SaveAs("rapidity.png");

    //apply emu method to get ee rapidity distributions

    TCanvas* rapCan2 = new TCanvas("rapCan2","ee Rapidity distributions");
    rapCan2->cd();

    TPad *rapPad2 = new TPad("rapPad2","rapPad2",0.1,0.1,0.95,0.995);
    rapPad2->Divide(2,3);
    rapPad2->Draw();
    rapPad2->SetLogy();
    //rapPpad1->SetLogx();
    rapPad2->cd();

    TH1F *eeRapidity[num_massBins2D];
    for (int i = 1; i < num_massBins2D-1; ++i){
      double massB = massBins2D[i];
      eeRapidity[i] = subtractEMubackground3(data_dmdyHistv[massB].at(0).get(), mc_dmdyHistv[massB],stats_mc_dmdyHistv[massB]);
      //eeRapidity[i] = subtractEMubackground3(sum_dmdyHistv[massB].at(0).get(), mc_dmdyHistv[massB],stats_mc_dmdyHistv[massB]); //closure test
      //change section in pad
      rapPad2->cd(i);
      eeRapidity[i]->SetXTitle("rapidity");
      eeRapidity[i]->GetXaxis()->SetTitleSize(0.08);
      eeRapidity[i]->GetXaxis()->SetLabelSize(0.07);
      eeRapidity[i]->SetMaximum(eeRapidity[i]->GetMaximum()*1.3);
      eeRapidity[i]->Draw("e");
      eeDMDYTot[massB]->Draw("same");
      eeRapidity[i]->Draw("e same");
      rapLeg[i]->Draw();

      //Fill matrix with 2D numbers
      unsigned int numXbins = eeRapidity[i]->GetNbinsX();
      for (unsigned int j = 0; j < numXbins; ++j){ //loop over rapidity bins
	true2eBackgroundFromData[i-1][j] = eeRapidity[i]->GetBinContent(j+1);
	true2eBackgroundFromDataError[i-1][j] = eeRapidity[i]->GetBinError(j+1);
	true2eBackgroundFromDataErrorSyst[i-1][j] = 0.5 * eeRapidity[i]->GetBinContent(j+1);
      }
    }

    rapCan2->SaveAs("eeRapidity.png");
    //should delete all the eeRapidity objects I have newed (memory leaks)
  }

  //==============================================
  //Calculate overall acceptance
  //==============================================

  double toteeEvts = sumVMassv.at(1)->Integral(1,nBinsV);
  double toteMuEvts = sumVMassv.at(0)->Integral(1,nBinsV);

  double acceptance = (2*toteeEvts)/toteMuEvts;

  if (verbose) cout << "Number of MC e-e events is: " << toteeEvts  << "\n";
  if (verbose) cout << "Number of MC e-mu events is: " << toteMuEvts  << "\n";
  if (verbose) cout << "Acceptance is: " << acceptance << "\n"; 

  //=======================================================================
  //Over lay ee and emu events with emu events rescaled by A/2
  //===================================================

  if (verbose) cout << "Number of MC e-mu events is: " << toteMuEvts  << "\n";
  if (verbose) cout << "Number of data e-mu events is: " << datahistv.at(0)->Integral(1,nBinsV)   << "\n";

  if (verbose) cout << "Bin by bin values for ee MC is:\n";
  if (verbose) printBinContents(sumVMassv.at(1).get());

  vector<shared_ptr<TH1F> > emuDistov;
  emuDistov.push_back(sumMassv.at(0));
  emuDistov.push_back(sumVMassv.at(0));
  emuDistov.push_back(datahMassv.at(0));
  emuDistov.push_back(datahistv.at(0));

  vector<shared_ptr<TH1F> > eMuMassScaledv;// this becomes the calculated ee distro
  eMuMassScaledv.push_back(shared_ptr<TH1F>(new TH1F("histo1","",nBins,xmin,xmax)));
  eMuMassScaledv.push_back(shared_ptr<TH1F>(new TH1F("histo2","",numMassBins,_massBinLimits2011)));
  eMuMassScaledv.push_back(shared_ptr<TH1F>(new TH1F("histo3","",nBins,xmin,xmax)));
  eMuMassScaledv.push_back(shared_ptr<TH1F>(new TH1F("histo4","",numMassBins,_massBinLimits2011)));
  //Fill above using a loop

  if ( emuDistov.size() != eMuMassScaledv.size()) 
    throw("vector lengths don't match");

  for (vector<shared_ptr<TH1F> >::size_type i = 0; i< emuDistov.size(); ++i){
    emuDistov.at(i)->Copy(*(eMuMassScaledv.at(i)));
    eMuMassScaledv.at(i)->Sumw2();// errors are correct. ie error/bincontent stays the same when scaled
    eMuMassScaledv.at(i)->Scale(acceptance*0.5);
    eMuMassScaledv.at(i)->SetLineColor(kRed);
  }

  TCanvas *eePad = new TCanvas("eePad","ee distributions",0,0,900,700);
  eePad->cd();
  
  eeVMassTot.Draw("");
  eMuMassScaledv.at(1)->Draw("same e");

  //===============================================================
  // compare data and emu MC
  //================================================================
  
  TCanvas *emuDataVMC = new TCanvas("emuDataVMC","Data v MC for emu",0,0,900,700);
  emuDataVMC->cd();

  TPad *emupad1 = new TPad("emupad1","emupad1",0.1,0.1,0.95,0.995);
  emupad1->Draw();
  emupad1->SetLogy();
  //emupad1->SetLogx();
  emupad1->cd();
 
  setEMUhistoStyle(datahistv.at(0).get());
  eMuVMassTot.Draw("");
  //datahistv.at(0)->SetXTitle("e^{+}e^{-} mass GeV/c^{2}");
  eMuVMassTot.GetXaxis()->SetTitle("e#mu mass GeV/c^{2}");
  datahistv.at(0)->Draw("e same");

  //add legend
  //Set up legend
  TLegend *stLeg = new TLegend(0.4 ,0.62,0.85,0.85);
  stLeg->SetFillColor(0);
  stLeg->SetTextSize(0.035);
  stLeg->SetHeader("e#mu data v MC distributions");
  stLeg->AddEntry(datahistv.at(0).get(), "data","lp");
  for (unsigned int ihist = 0; ihist < hist2Dv.at(0)->size(); ++ihist){
    if (ihist < entryLabels.size()){//make sure there is an label entry for all histos
      stLeg->AddEntry(hist2Dv.at(0)->at(ihist).get(), entryLabels.at(ihist),"f");;
    }
  }
  stLeg->Draw();

  //add text stating lumi
  TPaveText *txtptr = new TPaveText(0.35,0.92,0.55,0.97, "NDC"); // NDC sets coords
  txtptr->SetFillColor(0); // text is black on white
  txtptr->SetTextSize(0.03); 
  txtptr->SetBorderSize(0);
  txtptr->SetTextAlign(12);

  string lumiLabel = string("#int L dt = ") + lumiVal + string("fb^{-1}");

  txtptr->AddText(lumiLabel.c_str());
  //txtptr->AddText("#int L dt = 36.1 pb^{-1}");
  txtptr->Draw();
  //gPad->RedrawAxis();

  emuDataVMC->SaveAs("emu_emumass_MC_v_Data.png");
  //emuDataVMC->SaveAs("emu_emumass_MC_v_Data.C");
  //===============================================================
  //Estimate ee distro using emu MC
  //================================================================

  TH1F *eeEstimate = subtractEMubackground3(sumVMassv.at(0).get(),hist2Dv,stathist2Dv);

  if (verbose) cout << setprecision (6) << "Number of estimated ee events from e-mu MC is: " << eeEstimate->Integral(1,nBinsV)  << "\n";

  TCanvas *eeEstimatePad2 = new TCanvas("eeEstimatePad2","ee estimate using emu MC",0,0,900,700);
  eeEstimatePad2->cd();

  TPad *emupad2a = new TPad("emupad2a","",0.1,0.1,0.95,0.995);
  emupad2a->Draw();
  emupad2a->SetLogy();
  emupad2a->cd();

  setEEhistoStyle(eeEstimate);

  eeEstimate->Draw("e");
  eeVMassTot.Draw("same");
  eeEstimate->Draw("same e");

  TLegend *stLeg1a = new TLegend(0.4 ,0.62,0.85,0.85);
  stLeg1a->SetFillColor(0);
  stLeg1a->SetTextSize(0.035);
  stLeg1a->SetHeader("e^{+}e^{-} MC estimate v MC distributions");
  stLeg1a->AddEntry(datahistv.at(0).get(), "MC (from e#mu)","lp");
  for (unsigned int ihist = 0; ihist < hist2Dv.at(1)->size(); ++ihist){
    if (ihist < entryLabels.size()){//make sure there is an label entry for all histos
      stLeg1a->AddEntry(hist2Dv.at(0)->at(ihist).get(), entryLabels.at(ihist),"f");;
    }
  }
  stLeg1a->Draw(); 

  txtptr->Draw();

  eeEstimatePad2->SaveAs("emu_eemass_MC.png");
  //eeEstimatePad2->SaveAs("emu_eemass_MC.C");

  //==============================================
  //Plot emu MC
  //===============================================

  TCanvas *emuEst = new TCanvas("emuEst","emu stacked MC",0,0,900,700);
  emuEst->cd();

  TPad *emuEstpad = new TPad("emuEstpad","",0.1,0.1,0.95,0.995);
  emuEstpad->Draw();
  emuEstpad->SetLogy();
  emuEstpad->cd();
  
  //eMuVMassTot.SetYTitle("No. of Events");
  //eMuVMassTot.SetXTitle("e#mu mass GeV/c^{2}");

  eMuVMassTot.Draw("");
  emuEst->SaveAs("emuMassStackedMC.png");

  //===============================================================
  //Estimate ee distro using emu Data
  //================================================================

  TH1F *eeEstimate5 = subtractEMubackground3(datahistv.at(0).get(),hist2Dv,stathist2Dv);

  if (verbose) cout << setprecision (5) << "Number of estimated ee events from e-mu data is: " << eeEstimate5->Integral(1,nBinsV)  << "\n";

  TCanvas *eeEstimatePad6 = new TCanvas("eeEstimatePad6","ee estimate using emu Data",0,0,900,700);
  eeEstimatePad6->cd();

  TPad *emupad4 = new TPad("emupad4","",0.1,0.1,0.95,0.995);
  emupad4->Draw();
  emupad4->SetLogy();
  emupad4->cd();

  setEEhistoStyle(eeEstimate5);
  eeEstimate5->Draw("e");
  eeVMassTot.Draw("same");
  eeEstimate5->Draw("same e");

  //add legend
  //Set up legend
  TLegend *stLeg3 = new TLegend(0.4 ,0.62,0.85,0.85);
  stLeg3->SetFillColor(0);
  stLeg3->SetTextSize(0.035);
  stLeg3->SetHeader("e^{+}e^{-} data estimate v MC distributions");
  stLeg3->AddEntry(datahistv.at(0).get(), "data (from e#mu)","lp");
  for (unsigned int ihist = 0; ihist < hist2Dv.at(1)->size(); ++ihist){
    if (ihist < entryLabels.size()){//make sure there is an label entry for all histos
      stLeg3->AddEntry(hist2Dv.at(0)->at(ihist).get(), entryLabels.at(ihist),"f");;
    }
  }
  stLeg3->Draw(); 

  txtptr->Draw();
  eeEstimatePad6->SaveAs("emu_eemass_MC_v_Data.png");
  //Write out ee data info to a file
  if (saveToRootFile){
    TFile *outFile = new TFile(outFileName,"RECREATE");
    outFile->cd();    
    //eMuVMassTot.Write();
    //datahistv.at(0)->Write();
    true2eBackgroundFromData.Write("true2eBackgroundFromData");
    true2eBackgroundFromDataError.Write("true2eBackgroundFromDataError");
    true2eBackgroundFromDataErrorSyst.Write("true2eBackgroundFromDataErrorSyst");

  }

  cout << "press \"Ctrl + c\" to return prompt\n";
  //  graphicsPlease.Run();
  
  //=====================================================================
  return 0;
}

TH1F* eMu::subtractEMubackground3(TH1F *inputHisto, vector<vector<shared_ptr<TH1F> >* >& vvHist, vector<vector<shared_ptr<TH1F> >* >& statHist)
{
  //find the emu index
  //int eMuIndex = vHist.size()-1;
  vector<shared_ptr<TH1F> >* vHistMuon;
  vector<shared_ptr<TH1F> >* vHistElec;

  if (vvHist.size()< 2) 
    throw ("Need vector of size greater than 2");

  vHistMuon = vvHist.at(0);//point vvHist to the vector of pointers to TH1F muons
  vHistElec = vvHist.at(1);//point vvHist to the vector of pointers to TH1F electron

  //TH1F *acceptance;

  // int numXbins = vHistMuon[0]->GetNbinsX();
  unsigned int numXbins = vvHist.at(0)->at(0)->GetNbinsX();
  double xMin = vHistMuon->at(0)->GetXaxis()->GetXmin();
  double xMax = vHistMuon->at(0)->GetXaxis()->GetXmax();

  static int callCounter(0);// have a different index for each call of the functions
  ++callCounter;

  ostringstream sstreamCounter;
  sstreamCounter << callCounter;

  // Need to change below, very ugly!
  TH1F *outHisto;
  string histName("emuHisto");
  outHisto = new TH1F(TString(histName+sstreamCounter.str()),"",numMassBins,_massBinLimits2011);

  //need to check if this is the variable bin histo
  //ugly!

  {
    unsigned int chkNbins = outHisto->GetNbinsX();
    double chkxMin = outHisto->GetXaxis()->GetXmin();
    double chkxMax = outHisto->GetXaxis()->GetXmax();

    if ((chkNbins != numXbins) || (chkxMin != xMin) || (chkxMax != xMax)) {
      delete outHisto;
      outHisto = new TH1F(TString("outHisto")+sstreamCounter.str(),"",numXbins,xMin,xMax);
      if (verbose) cout << "Changing Binning\n";
    }

  }

   double sumOfHistBinIndex;
   //loop over all bins and calc R for each bin
   for (unsigned int i = 1; i < (numXbins+1); ++i){
          
     double dN_sq(0);
     double binContent(0);
         
     for (unsigned int k = 0; k <  emu_array_size; ++k) {//loop over each channel that exists in ee and emu
        //add emu events for each channel. Do same for ee.
       double emuContents(0);
       double eeContents(0);
       emuContents = vHistMuon->at(k)->GetBinContent(i);
       eeContents = vHistElec->at(k)->GetBinContent(i);
       
     
       sumOfHistBinIndex = 0; //initialise
       for (unsigned int j = 0; j < vHistMuon->size(); ++j){
         if (j != k){//add contents of all the other bins
       	 sumOfHistBinIndex += vHistMuon->at(j)->GetBinContent(i);
         }
       }//end of loop over each emu channel
       
       //calulate acceptance and r for each bin
       double r;
       double accept(0);
       if ((emuContents > 0) && (eeContents > 0)){
         accept = (2*eeContents)/emuContents;
         double numerator = sumOfHistBinIndex/(sumOfHistBinIndex + emuContents);
         double denominator = emuContents/(sumOfHistBinIndex + emuContents);
         r = numerator/denominator;
       } else {
         r = 0;
         accept = 0;
       }

       //calculate the error on ee
       
       //first calc error on A
       double statEmuContents(0);
       double statEeContents(0);
       statEmuContents = statHist.at(0)->at(k)->GetBinContent(i);
       statEeContents = statHist.at(1)->at(k)->GetBinContent(i);
       
       
       //Error on A explaination
       //A=2Nee/Nemu
       // The error on A squared (dA^2) is (dA/dNee)^2*(dNee)^2 + (dA/dNemu)^2*(dNemu)^2
       // This equates to (2/Nemu)^2*(dNee)^2 + (-2Nee/Nemu^2)^2*dNemu^2
       // Remember Nee = dNee^2 and Nemu = dNemu^2
       // Therefore we can write error on A squared as
       //     double errorA_sq = (2/statEmuContents) * (2/statEmuContents) * statEeContents
       // + ((2*statEeContents)/(statEmuContents* statEmuContents))*((2*statEeContents)/(statEmuContents* statEmuContents))*statEmuContents;
       
       //more efficient (simplifies) as
       double errorA_sq; 
       if ((statEeContents == 0) || (statEmuContents == 0)){	   
	 errorA_sq = 0;
       } else {
	 errorA_sq =  4 * statEeContents * (statEmuContents + statEeContents) / (statEmuContents* statEmuContents *statEmuContents);
       }
       
       double eMuObs = inputHisto->GetBinContent(i);
       double dNdA = eMuObs*(1/(1+r))*0.5;
       double dNdNemu = accept*(1/(1+r))*0.5;
       
       dN_sq += dNdA*dNdA*errorA_sq + dNdNemu*dNdNemu*eMuObs;
       binContent += dNdA*accept;   

       //Need to use limit on 0 events as the error for 0 events
       // at 90% confidence the limit for 0 events is 2.3
       //90% is 1.64 sigma. 
       //68% confidence or 1 sigma is 1.145. ie. What mean value will have a Poisson probability of 0.318 at 0
       if (binContent==0.0){
	 if (dNdNemu == 0.0){
	   dN_sq += (0.93/2.0)*1.145*1.145; //need to update this hardwired in the average exceptance
	 } else {
           dN_sq += dNdA*dNdA*errorA_sq + dNdNemu*dNdNemu*1.145*1.145;
	 }
       }
     } //end of loop over each channel that exists in ee and emu    
     outHisto->SetBinContent(i,binContent);
     
     double err = sqrt(dN_sq+(binContent));
     outHisto->SetBinError(i,err);
     
     //true2eBackgroundFromData[i-1] = binContent;
     //true2eBackgroundFromDataError[i-1] = err;
     //true2eBackgroundFromDataErrorSyst[i-1] = binContent * 0.05;// this is 1.3% + 4.8% error for ee and emu MC agreement respectively

     if (verbose) cout << setprecision (4) << binContent << "$\\pm$" << err << "\n";
   }
   return outHisto;
} 

double eMu::calcError(TH1F *inputHisto, vector<vector<shared_ptr<TH1F> >* >& vvHist)
{
  unsigned int numXbins = vvHist.at(0)->at(0)->GetNbinsX();
  double eMuObs = inputHisto->Integral(1,numXbins);
  return calcError(eMuObs, vvHist);
}

double eMu::calcError(const double& eMuObs, vector<vector<shared_ptr<TH1F> >* >& vvHist)
{
  //find the emu index
  //int eMuIndex = vHist.size()-1;
  vector<shared_ptr<TH1F> >* vHistMuon;
  vector<shared_ptr<TH1F> >* vHistElec;

  if (vvHist.size()< 2) 
    throw ("Need vector of size greater than 2");

  vHistMuon = vvHist.at(0);//point vvHist to the vector of pointers to TH1F muons
  vHistElec = vvHist.at(1);//point vvHist to the vector of pointers to TH1F electron

  // int numXbins = vHistMuon[0]->GetNbinsX();
  unsigned int numXbins = vvHist.at(0)->at(0)->GetNbinsX();

  //add emu events for each channel. Do same for ee.
  double emuContents(0);
  double eeContents(0);
  double accept(0);
         
  for (unsigned int k = 0; k < emu_array_size; ++k) {       
    emuContents += vHistMuon->at(k)->Integral(1,numXbins);
    eeContents += vHistElec->at(k)->Integral(1, numXbins);
  }  
  accept = (2*eeContents)/emuContents;
    
  double emuBkg(0);
  for (unsigned int k =  emu_array_size; k < vHistMuon->size(); ++k) {
    emuBkg += vHistMuon->at(k)->Integral(1,numXbins);
  }

  double r = emuBkg/emuContents;

  //calculate the error on ee
    
  //first calc error on A
  //Error on A explaination
  //A=2Nee/Nemu
  // The error on A squared (dA^2) is (dA/dNee)^2*(dNee)^2 + (dA/dNemu)^2*(dNemu)^2
  // This equates to (2/Nemu)^2*(dNee)^2 + (-2Nee/Nemu^2)^2*dNemu^2
  // Remember Nee = dNee^2 and Nemu = dNemu^2
  // Therefore we can write error on A squared as
  //     double errorA_sq = (2/statEmuContents) * (2/statEmuContents) * statEeContents
  // + ((2*statEeContents)/(statEmuContents* statEmuContents))*((2*statEeContents)/(statEmuContents* statEmuContents))*statEmuContents;
      
      //more efficient (simplifies) as
  double errorA_sq =  4 * eeContents * (emuContents + eeContents) / (emuContents* emuContents *emuContents);
      
  double dNdA = eMuObs*(1/(1+r))*0.5;
  double dNdNemu = accept*(1/(1+r))*0.5;       
  double dN_sq = dNdA*dNdA*errorA_sq + dNdNemu*dNdNemu*eMuObs;
  
  return sqrt(dN_sq);
} 

void printBinContents(TH1F* myhist)
{
  unsigned int numXbins = myhist->GetNbinsX();
  for (unsigned int i = 1; i < (numXbins+1); ++i){
      cout << setprecision (3) << "&" << myhist->GetBinContent(i) << "";
      cout << setprecision (2) << "$\\pm$" << myhist->GetBinError(i) << "\n";
  }
  cout << setprecision (8) << "\\\\ \\hline";
}

void printBinContents(vector<vector<shared_ptr<TH1F> >* >& vvHist)
{
  unsigned int numXbins = vvHist.at(1)->at(0)->GetNbinsX();

  unsigned int numChans = vvHist.at(1)->size();

  for (unsigned int k = 0; k < numChans; ++k) {
    for (unsigned int i = 1; i < (numXbins+1); ++i){
      cout << setprecision (2) << "&" << vvHist.at(1)->at(k)->GetBinContent(i) << "";
    }
    cout << "\\\\ \\hline \n";
  }
}

void printBinContents2(TH1F* inHisto, vector<vector<shared_ptr<TH1F> >* >& vvHist)
{
  unsigned int numXbins = vvHist.at(1)->at(0)->GetNbinsX();

  unsigned int numChans = vvHist.at(1)->size();

  for (unsigned int k = 0; k < numChans; ++k) {
    for (unsigned int i = 1; i < (numXbins+1); ++i){
      //cout << setprecision (2) << "&" << vvHist.at(1)->at(k)->GetBinContent(i) << "";
      double numConts(0); 
      for (unsigned int j = 0; j < numChans; ++j) {
	numConts +=  vvHist.at(1)->at(j)->GetBinContent(i);
      }      
      cout << setprecision (2) << "&" << vvHist.at(1)->at(k)->GetBinContent(i)/numConts *inHisto->GetBinContent(i) << "";
    }
    cout << setprecision (8) << "\\\\ \\hline";
  }
}

void sethistoStyle(TH1F *inputHisto)
{
  inputHisto->SetMarkerStyle(20);
  inputHisto->SetLineWidth(3);
  inputHisto->SetYTitle("No. of Events");
}

void setEEhistoStyle(TH1F *inputHisto)
{
  inputHisto->SetXTitle("e^{+}e^{-} mass GeV/c^{2}");
  sethistoStyle(inputHisto);
}

void setEMUhistoStyle(TH1F *inputHisto)
{
  inputHisto->SetXTitle("e#mu mass GeV/c^{2}");
  sethistoStyle(inputHisto);
}


