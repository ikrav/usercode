#ifndef JSONPARSER
#define JSONPARSER

#include <TString.h>
#include <vector>

class JsonParser {

public:
  
  JsonParser();
  virtual ~JsonParser(){};
  
  void Initialize(TString filename);
  bool HasRunLumi(int run, int lumi);
  void Reset();
  void Print();
  int AssembleNumber(vector<int> data);

  enum Depth {
    d_outside           = 0,
    d_insideOuterBracket = 1,
    d_insideInnerBracketNumberOne = 2,
    d_insideInnerBracketNumberTwo = 3
  };


private:

  vector <int>            _runList;
  vector < vector<int> >  _lumiStatusList;

  bool                    _isInitialized;

  ClassDef(JsonParser,1)

};
#endif
