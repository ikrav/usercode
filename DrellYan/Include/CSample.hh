#ifndef CSAMPLE_HH
#define CSAMPLE_HH

//
// helper class to handle sample inputs
//
class CSample
{
public:
  CSample(){}
  ~CSample(){}
  
  TString          label;    // plot item label
  Int_t            color;    // plot item color
  vector<TString>  fnamev;   // ntuple files
  vector<Double_t> xsecv;    // per file cross section
  vector<TString>  jsonv;    // per file JSON file
  vector<Double_t> weightv;  // per file event weight
};

#endif
