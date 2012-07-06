#include "../Include/InputFileMgr.hh"
#include <assert.h>
#include <fstream>
#include <sstream>

// -----------------------------------------------------------
// -----------------------------------------------------------

void InitialInputMgr_t::Clear() {
  FLoadedFileName.Clear(); FOutputDir.Clear();
  FSavePlotFormat.Clear();
  FTotLumi=0.;
  FWeightEvents=1;
  FEnergyScaleTag.Clear();
  FGenerateEEMFile.Clear();
  FSampleNames.clear();
  ClearVec(FSampleInfos);
}

// -----------------------------------------------------------

int InitialInputMgr_t::Load(const TString &inputfname) {
  std::ifstream ifs;
  ifs.open(inputfname.Data());
  if (!ifs.is_open()) {
    std::cout << "failed to load input file <" << inputfname << ">\n";
    throw 2;
  }
  std::string line;
  Int_t state=0;
  CSample *sample=NULL;
  while(getline(ifs,line)) {
    if ((line[0]=='#') && (line[1]=='$') && (line[2]=='$')) {
      if (line.find("generate_EEM_files=") != string::npos) {
	FGenerateEEMFile=line.substr(line.find('=')+1);
	std::cout << "\n\tEEM files will be generated, tag=<" << FGenerateEEMFile << ">\n\n";
	continue;
      }
    }
    if(line[0]=='#') continue;
    if(line[0]=='%') { 
      state++; 
      continue; 
    }
    if(line[0]=='$') {
      sample=new CSample();
      FSampleInfos.push_back(sample);
      stringstream ss(line);
      string chr;
      string sname;
      Int_t color;
      ss >> chr >> sname >> color;
      string label = line.substr(line.find('@')+1);
      sample->label = label;
      sample->color = color;
      FSampleNames.push_back(sname);
      continue;
    }
    
    if(state==0) {  // general settings
      stringstream ss1(line); ss1 >> FTotLumi;
      getline(ifs,line);
      stringstream ss2(line); ss2 >> FWeightEvents;
      getline(ifs,line);
      FOutputDir = TString(line);
      getline(ifs,line);
      // backwards compatibility for the input file
      if (line.size()>3) {  // escale is defined
	FEnergyScaleTag=TString(line);
	getline(ifs,line);
	// check that it was correct to use this work-around
	if (line.find('%')!=std::string::npos) {
	  std::cout << "backwards-compatibility code failure\n";
	  return 0;
	}
      }
      FSavePlotFormat = TString(line);      

    } else if(state==1) {  // define data sample
      string samplefname;
      Double_t xsec;
      string json;
      stringstream ss(line);
      ss >> samplefname >> xsec >> json;
      sample->fnamev.push_back(samplefname);
      sample->xsecv.push_back(xsec);
      sample->jsonv.push_back(json);
    
    } else if(state==2) {  // define MC samples
      string samplefname;
      Double_t xsec;
      stringstream ss(line);
      ss >> samplefname >> xsec;
      sample->fnamev.push_back(samplefname);
      sample->xsecv.push_back(xsec);
    }
  }
  ifs.close();
  
  FLoadedFileName=inputfname;
  return 1;
}

// -----------------------------------------------------------

// -----------------------------------------------------------
// -----------------------------------------------------------

int TnPInputFileMgr2011_t::Load(const TString &configFile) {
  this->clear();
  std::ifstream ifs;
  ifs.open(configFile.Data());
  if (!ifs.is_open()) {
    std::cout << "TnPInputFileMgr2011::Load  tried to open the configuration file <" << configFile << ">\n";
    return 0;
  }
  std::string line;
  int state=0;
  while(getline(ifs,line)) {
    if(line[0]=='#') continue;
    if(line[0]=='%') break;
    if(state==0){
      // Read 1st line of content: data or MC?
      FSampleTypeStr = TString(line);
      state++;
    }else if(state==1) {
      // Read 2d content line: efficiency type string
      FEffTypeStr = TString(line);
      state++;
    }else if(state==2) {
      // Read 3d content line: fitting mode
      FCalcMethodStr = TString(line);
      state++;
    }else if(state==3) {
      // Read 4th content line: SC ET binning
      FEtBinsKindStr = TString(line);
      state++;
    }else if(state==4) {
      // Read 5th content line: SC eta binning
      FEtaBinsKindStr = TString(line);
      state++;
    }else if(state==5) {
      // Read 5th content line: SC eta binning
      FDirTag = TString(line);
      state++;
    }else if(state==6) {
      FFileNames.push_back(TString(line));
    }
  }
  ifs.close();
  int res=((state==6) && FFileNames.size()) ? 1:0;
  if (!res) {
    std::cout << "Failed to load file <" << configFile << ">\n";
  }
  return res;
}

// -----------------------------------------------------------
// -----------------------------------------------------------

int TnPInputFileMgr_t::Load(const TString &configFile) {
  this->clear();
  std::ifstream ifs;
  ifs.open(configFile.Data());
  if (!ifs.is_open()) {
    std::cout << "TnPInputFileMgr::Load tried to open configFile=<" << configFile << ">\n";
    return 0;
  }
  string line;
  Int_t state=0;
  Int_t subState=0;
  while(getline(ifs,line)) {
    if(line[0]=='#') continue;
    if(line[0]=='%') break;
    if(state==0){
      // Read 1st line of content: data or MC?
      FSampleTypeStr = TString(line);
      state++;
    }else if(state==1) {
      // Read 2d content line: efficiency fitting mode
      size_t pos=line.find(':');
      if (pos==string::npos) {
        std::cout << "expected format is EFFICIENCY:fitting_mode\n";
        return 0;
      }
      subState++;
      int orderErr=0;
      switch(subState) {
      case 1: if (line.find("RECO")==string::npos) orderErr=1; break;
      case 2: if (line.find("ID"  )==string::npos) orderErr=1; break;
      case 3: if (line.find("HLT" )==string::npos) orderErr=1; break;
      default:
	std::cout << "incorrect value of subState=" << subState << "\n";
	return 0;
      }
      if (orderErr) {
	std::cout << "EfficiencyKind:CalculationMethod should be ordered RECO,ID,HLT\n";
	return 0;
      }
      FEffTypeStrV.push_back(TString(line.substr(0,pos)));
      FCalcMethodStrV.push_back(TString(line.c_str()+pos+1));
      if (subState==3) state++;
    }else if(state==2) {
      // Read 3rd content line: SC ET binning
      FEtBinsKindStr = TString(line);
      state++;
    }else if(state==3) {
      // Read 4th content line: SC eta binning
      FEtaBinsKindStr = TString(line);
      state++;
    }else if(state==4) {
      // Read 5th content line: directory tag
      FDirTag = TString(line);
      state++;
    }else if(state==5) {
      FFileNames.push_back(TString(line));
    }
  }
  ifs.close();
  if (state!=5) {
    std::cout << "Load: did not reach state=5\n";
    return 0;
  }

  //std::cout << "Loaded:\n" << *this << "\n";
  return 1;
}

// -----------------------------------------------------------

bool TnPInputFileMgr_t::hasSameBinCounts(const TnPInputFileMgr_t &mgr) const {
  return ((FEtBinsKindStr==mgr.FEtBinsKindStr) &&
	  (FEtaBinsKindStr==mgr.FEtaBinsKindStr) &&
	  (FEffTypeStrV.size() == mgr.FEffTypeStrV.size()) &&
	  (FCalcMethodStrV.size() == mgr.FCalcMethodStrV.size())) ? kTRUE : kFALSE;
}

// -----------------------------------------------------------
// -----------------------------------------------------------

int MCInputFileMgr_t::Load(const TString& inputFileName) {
  std::ifstream ifs;
  ifs.open(inputFileName.Data());
  if (!ifs.is_open()) {
    std::cout << "MCInputFileMgr: failed to open file <" << inputFileName << ">\n";
    return 0;
  }
  std::string line;
  Int_t state=0;
  while(getline(ifs,line)) {
    if(line[0]=='#') continue;
    if(state == 0){
      FDirTag = TString(line);
      state++;
      continue;
    }else{
      std::string fname;
      Int_t color1, linesty;
      std::stringstream ss(line);
      Double_t xsec1;
      ss >> fname >> xsec1 >> color1 >> linesty;
      string label1 = line.substr(line.find('@')+1);
      if (label1.size()) {
	FFileNames.push_back(fname);
	FLabels.push_back(label1);
	FColors.push_back(color1);
	FLineStyles.push_back(linesty);
	FXSecs.push_back(xsec1);
	FLumis.push_back(0);
      }
    }
  }
  ifs.close();
  std::cout << "Loaded:\n" << *this << "\n";
  return FFileNames.size();
}

// -----------------------------------------------------------

// -----------------------------------------------------------
// -----------------------------------------------------------

