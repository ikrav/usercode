#ifndef InputFileMgr_HH
#define InputFileMgr_HH

#include <TString.h>
#include <vector>
#include <iostream>
#include "../Include/DYToolsUI.hh"
#include "../Include/CSample.hh"
#include "../Include/MyTools.hh"

// --------------------------------------------------------

std::ostream& operator<<(std::ostream& out, const CSample &s) {
  out << "   CSample{ label=" << s.label << ", color=" << s.color << ", " << s.fnamev.size() << " files }:\n";
  for (unsigned int i=0; i<s.fnamev.size(); i++) {
    out << " " << (i+1) << ". xsec=" << s.xsecv[i] << "  , weight=" << s.weightv[i] << ", fname=<" << s.fnamev[i] << ">, json=<" << s.jsonv[i] << "\n";
  }
  return out;
}

// --------------------------------------------------------

class InitialInputMgr_t {
protected:
  TString FLoadedFileName;
  TString FOutputDir, FSavePlotFormat;
  double FTotLumi;
  int FWeightEvents;
  TString FEnergyScaleTag;
  TString FGenerateEEMFile;
  std::vector<TString> FSampleNames;
  std::vector<CSample*> FSampleInfos;
public:
  InitialInputMgr_t() : 
    FLoadedFileName(),
    FOutputDir(), FSavePlotFormat(), FTotLumi(0.),
    FWeightEvents(1), FEnergyScaleTag(), FGenerateEEMFile(),
    FSampleNames(), FSampleInfos() 
  {}

  void Clear();

  // Access
  const TString& loadedFileName() const { return FLoadedFileName; }
  const TString& outputDir() const { return FOutputDir; }
  const TString& savePlotFormat() const { return FSavePlotFormat; }
  double totalLumi() const { return FTotLumi; }
  const TString& energyScaleTag() const { return FEnergyScaleTag; }
  const TString& generateEEMFile() const { return FGenerateEEMFile; }
  unsigned int sampleCount() const { return FSampleNames.size(); }
  const std::vector<TString>& sampleNames() const { return FSampleNames; }
  const std::vector<CSample*>& sampleInfos() const { return FSampleInfos; }

  const TString& sampleName(unsigned int i) const { return FSampleNames[i]; }
  const CSample* sampleInfo(unsigned int i) const { return FSampleInfos[i]; }

  const TString dirTag() const {
    int idx=FOutputDir.Index("/DY");
    TString tag= (idx==-1) ? FOutputDir : FOutputDir(idx,FOutputDir.Length());
    return tag;
  }

  // Load
  int Load(const TString &inputFile);

  // output
  friend std::ostream& operator<<(std::ostream &out, const InitialInputMgr_t &m) {
    out << "InitialInputMgr(loadedFile=<" << m.FLoadedFileName << ">, ";
    out << "outputDir=<" << m.FOutputDir << ">, savePlotFormat=<" << m.FSavePlotFormat << ">, totalLumi=" << m.FTotLumi << ", weightEvents=" << m.FWeightEvents << ", energyScaleTag=<" << m.FEnergyScaleTag << ">, generateEEMFile=<" << m.FGenerateEEMFile << ">, " << m.FSampleNames.size() << " sample names and " << m.FSampleInfos.size() << " sample infos):\n";
    for (unsigned int i=0; i<m.FSampleNames.size(); ++i) {
      out << " " << (i+1) << " sampleName=" << m.FSampleNames[i] << ":\n";
      out << ( *(m.FSampleInfos[i]) ) << "\n";
    }
    return out;
  }
};

// --------------------------------------------------------

// --------------------------------------------------------
// 2011-style input file

class TnPInputFileMgr2011_t { 
protected:
  TString FSampleTypeStr;
  TString FEffTypeStr, FCalcMethodStr;
  TString FEtBinsKindStr, FEtaBinsKindStr;
  TString FDirTag;
  std::vector<TString> FFileNames;
public:
  TnPInputFileMgr2011_t() : FSampleTypeStr(), FEffTypeStr(), FCalcMethodStr(),
			    FEtBinsKindStr(), FEtaBinsKindStr(), FDirTag(),
			    FFileNames() {}

  // cleanup
  void clear() {
    FSampleTypeStr.Clear(); FEffTypeStr.Clear(); FCalcMethodStr.Clear();
    FEtBinsKindStr.Clear(); FEtaBinsKindStr.Clear(); FDirTag.Clear();
    FFileNames.clear();
  }

  // access
  const TString& sampleTypeString() const { return FSampleTypeStr; }
  const TString& effTypeString() const { return FEffTypeStr; }
  const TString& calcMethodString() const { return FCalcMethodStr; }
  const TString& etBinsKindString() const { return FEtBinsKindStr; }
  const TString& etaBinsKindString() const { return FEtaBinsKindStr; }
  const TString& dirTag() const { return FDirTag; }
  const std::vector<TString> fileNames() const { return FFileNames; }
  unsigned int fileCount() const { return FFileNames.size(); }

  const TString& fileName(unsigned int i) const { return FFileNames[i]; }
  const TString& operator[](unsigned int i) const { return FFileNames[i]; }

  // Access with conversion
#ifdef DYToolsUI_HH
  DYTools::TDataKind_t sampleType() const { return DetermineDataKind(FSampleTypeStr); }
  DYTools::TEfficiencyKind_t effType() const { return DetermineEfficiencyKind(FEffTypeStr); }
  DYTools::TTnPMethod_t calcMethod() const { return DetermineTnPMethod(FCalcMethodStr); }
  DYTools::TEtBinSet_t etBinsKind() const { return DetermineEtBinSet(FEtBinsKindStr); }
  DYTools::TEtaBinSet_t etaBinsKind() const { return DetermineEtaBinSet(FEtaBinsKindStr); }
#endif

  // Load 

  int Load(const TString &inputFile);

  friend std::ostream& operator<< (std::ostream& out, const TnPInputFileMgr2011_t &mgr) {
    out << "EffStudyInputMgr2011 (sampleType=" << mgr.FSampleTypeStr << ", effType=" << mgr.FEffTypeStr << ", calcMethod=" << mgr.FCalcMethodStr << ", etBinning=" << mgr.FEtBinsKindStr << ", etaBinning=" << mgr.FEtaBinsKindStr << ", dirTag=<" << mgr.FDirTag << ">, " << mgr.FFileNames.size() << " ntuple files):\n";
    for (unsigned int i=0; i<mgr.FFileNames.size(); ++i) {
      out << (i+1) << " -  " << mgr.FFileNames[i] << "\n";
    }
    return out;
  }
};

// --------------------------------------------------------

class TnPInputFileMgr_t {
protected:
  TString FSampleTypeStr;
  std::vector<TString> FEffTypeStrV, FCalcMethodStrV;
  TString FEtBinsKindStr, FEtaBinsKindStr;
  TString FDirTag;
  std::vector<TString> FFileNames;
public:
  TnPInputFileMgr_t() : FSampleTypeStr(), FEffTypeStrV(), FCalcMethodStrV(),
			 FEtBinsKindStr(), FEtaBinsKindStr(), FDirTag(),
			 FFileNames() {}

  // cleanup
  void clear() {
    FSampleTypeStr.Clear(); FEffTypeStrV.clear(); FCalcMethodStrV.clear();
    FEtBinsKindStr.Clear(); FEtaBinsKindStr.Clear(); FDirTag.Clear();
    FFileNames.clear();
  }

  // access
  const TString& sampleTypeString() const { return FSampleTypeStr; }
  template<class idx_t>
  const TString& effTypeString(idx_t idx) const { return FEffTypeStrV[idx]; }
  template<class idx_t>
  const TString& calcMethodString(idx_t idx) const { return FCalcMethodStrV[idx]; }
  const TString& etBinsKindString() const { return FEtBinsKindStr; }
  const TString& etaBinsKindString() const { return FEtaBinsKindStr; }
  const TString& dirTag() const { return FDirTag; }
  const std::vector<TString> fileNames() const { return FFileNames; }
  unsigned int fileCount() const { return FFileNames.size(); }

  const TString& fileName(unsigned int i) const { return FFileNames[i]; }
  const TString& operator[](unsigned int i) const { return FFileNames[i]; }

  bool hasSameBinCounts(const TnPInputFileMgr_t &mgr) const;

  // Access with conversion
#ifdef DYToolsUI_HH
  DYTools::TDataKind_t sampleType() const { return DetermineDataKind(FSampleTypeStr); }
  template<class idx_t>
  DYTools::TEfficiencyKind_t effType(idx_t idx) const { return DetermineEfficiencyKind(FEffTypeStrV[idx]); }
  template<class idx_t>
  DYTools::TTnPMethod_t effCalcMethod(idx_t idx) const { return DetermineTnPMethod(FCalcMethodStrV[idx]); }
  int etBinsCount() const { return DYTools::getNEtBins(int(this->etBinsKind())); }
  DYTools::TEtBinSet_t etBinsKind() const { return DetermineEtBinSet(FEtBinsKindStr); }
  DYTools::TEtaBinSet_t etaBinsKind() const { return DetermineEtaBinSet(FEtaBinsKindStr); }
#endif

  // Load 

  int Load(const TString &inputFile);

  friend std::ostream& operator<< (std::ostream& out, const TnPInputFileMgr_t &mgr) {
    out << "EffStudyInputMgr: sampleType=" << mgr.FSampleTypeStr;
    out << ", etBinning=" << mgr.FEtBinsKindStr 
	<< ", etaBinning=" << mgr.FEtaBinsKindStr 
	<< ", dirTag=<" << mgr.FDirTag << ">"; 
    for (unsigned int i=0; i<mgr.FEffTypeStrV.size(); ++i) {
      out << "; effType=" << mgr.FEffTypeStrV[i] 
	  << ", calcMethod=" << mgr.FCalcMethodStrV[i];
    }
    out << "; " << mgr.FFileNames.size() << " ntuple files:\n";
    for (unsigned int i=0; i<mgr.FFileNames.size(); ++i) {
      out << "   " << mgr.FFileNames[i] << "\n";
    }
    return out;
  }
};

// --------------------------------------------------------
// --------------------------------------------------------

class MCInputFileMgr_t {
protected:
  TString FDirTag, FEScaleTag;
  std::vector<TString> FFileNames,FLabels;
  std::vector<Int_t> FColors,FLineStyles;
  std::vector<Double_t> FXSecs,FLumis;
public:
  MCInputFileMgr_t() : FDirTag(),FEScaleTag("20120101_default"),
		       FFileNames(),FLabels(),FColors(),
		       FLineStyles(),FXSecs(),FLumis() {}

  unsigned int size() const { return FFileNames.size(); }
  const TString& dirTag() const { return FDirTag; }
  const TString& escaleTag() const { return FEScaleTag; }
  const std::vector<TString>& fileNames() const { return FFileNames; }
  const std::vector<TString>& labels() const { return FLabels; }
  const std::vector<Int_t>& colors() const { return FColors; }
  const std::vector<Int_t>& lineStyles() const { return FLineStyles; }
  const std::vector<Double_t>& xsecs() const { return FXSecs; }
  const std::vector<Double_t>& lumis() const { return FLumis; }
  template<class idx_t>
  const TString& fileName(idx_t idx) const { return FFileNames[idx]; }
  template<class idx_t>
  const TString& label(idx_t idx) const { return FLabels[idx]; }
  template<class idx_t>
  Int_t color(idx_t idx) const { return FColors[idx]; }
  template<class idx_t>
  Int_t lineStyle(idx_t idx) const { return FLineStyles[idx]; }
  template<class idx_t>
  Double_t xsec(idx_t idx) const { return FXSecs[idx]; }
  template<class idx_t>
  Double_t lumi(idx_t idx) const { return FLumis[idx]; }
  
  int Load(const TString &inputFileName);

  friend std::ostream& operator<<(std::ostream& out, const MCInputFileMgr_t &m) {
    out << "mcInputFileMgr> (" << m.size() << " items):\n";
    for (unsigned int i=0; i<m.size(); ++i) {
      out << "  #" << i << "  <" << m.FFileNames[i] 
	  << "> label=<" << m.FLabels[i]
	  << "> color=" << m.FColors[i]
	  << " lineStyle=" << m.FLineStyles[i]
	  << " xsec=" << m.FXSecs[i] 
	  << " lumi=" << m.FLumis[i]
	  << "\n";
    }
    std::cout << " dirTag=<" << m.FDirTag << ">"
	      << ", escaleTag=<" << m.FEScaleTag 
	      << ">\n";
    return out;
  }
};

// --------------------------------------------------------
// --------------------------------------------------------

class XSecInputFileMgr_t {
protected:
  double FTotLumi;
  TString FName;
  TString FYieldsTag, FConstTag, FEvtEffScaleTag;
  TString FTrigSet;
public:
  XSecInputFileMgr_t() : FTotLumi(0.), FName(),
			 FYieldsTag(), FConstTag(), 
			 FEvtEffScaleTag(),
			 FTrigSet() {}
  void Clear();

  double totLumi() const { return FTotLumi; }
  const TString& name() const { return FName; }
  const TString& yieldsTag() const { return FYieldsTag; }
  const TString& constTag() const { return FConstTag; }
  const TString& evtEffScaleTag() const { return FEvtEffScaleTag; }
  const TString& triggerSetName() const { return FTrigSet; }
  
  int Load(const TString &fname);

  friend std::ostream& operator<<(std::ostream& out, const XSecInputFileMgr_t &mgr) {
    out << "XSecInputFileMgr:\n";
    out << "  -- name=<" << mgr.name() << ">\n";
    out << "  -- yieldsTag=<" << mgr.yieldsTag() << ">\n";
    out << "  -- constTag=<" << mgr.constTag() << ">\n";
    out << "  -- evtEffScaleTag=<" << mgr.evtEffScaleTag() << ">\n";
    out << "  -- triggerSet=<" << mgr.triggerSetName() << ">";
    out << std::endl;
    return out;
  }
};

// --------------------------------------------------------
// --------------------------------------------------------

class EScaleTagFileMgr_t {
protected:
  std::vector<TString> FEScaleTags;
public:
  EScaleTagFileMgr_t() : FEScaleTags() {}
  unsigned int size() const { return FEScaleTags.size(); }
  const TString& escaleTag(int i) const { return FEScaleTags[i]; }
  const std::vector<TString>& getEScaleTags() const { return FEScaleTags; }

  int Load(const TString &inputFileName);
  
  friend std::ostream& operator<<(std::ostream& out, const EScaleTagFileMgr_t &m) {
    out << "eScaleTagFileMgr> (" << m.size() << " items):\n";
    for (unsigned int i=0; i<m.size(); i++) {
      out << "  #" << i << " <" << m.FEScaleTags[i] << ">\n";
    }
    return out;
  }
};


// --------------------------------------------------------
// --------------------------------------------------------


class TDescriptiveInfo_t : public TObject {
public:
  TDescriptiveInfo_t() : TObject(), _info() {}
  ~TDescriptiveInfo_t() {}
  
  std::vector<std::string> _info;

  void print() const {
    std::cout << std::string(80,'-') << "\n";
    std::cout << "File info " << _info.size() << " lines:\n";
    for (unsigned int i=0; i<_info.size(); ++i) {
      std::cout << _info[i];
      if (_info[i].find('\n')==std::string::npos) std::cout << "\n";
    }
    std::cout << std::string(80,'-') << "\n";
  }

  ClassDef(TDescriptiveInfo_t,1)
};

// --------------------------------------------------------

#endif
