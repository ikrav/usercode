#ifndef SAMPLEITALY_H
#define SAMPLEITALY_H

#include "SampleBase.h"

#include <TString.h>
#include <TROOT.h>
#include <TChain.h>
#include <TSystem.h>

class SampleItaly : public SampleBase {

 public:

  SampleItaly(TString dir = "", TString base = "", 
	    TString end = "", int nfiles = 0);
  ~SampleItaly();

  void initialize();

  //inline void setCategory(int cat){_category = cat;};
  //inline void setInitialized(bool b){_isInitialized = b;};

  //inline int getCategory(){return _category;};
  //inline bool isInitialized(){return _isInitialized;};

 protected:

  // Input files info
  TString   _dataDir;
  TString   _fileNameBase;
  TString   _fileNameEnding;
  int       _nFiles;

  // Access
  TChain   *_chain;


  // Variable naming
  int       _category;
  TString   _prefix;
  TString   _suffix;

  // Other
  //bool      _isInitialized;

  ClassDef(SampleItaly,1)

};
#endif
