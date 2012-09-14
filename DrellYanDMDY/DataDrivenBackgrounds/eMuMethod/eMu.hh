#ifndef EMU_HH
#define EMU_HH

//STL Headers
#include <string>
#include <vector>
#include <map>
#include <utility>

//BOOST Headers
#include <boost/shared_ptr.hpp>

//ROOT Headers
#include <TH1.h>
#include <TString.h>
#include <TVectorT.h>
#include <TMatrixT.h>

using std::string;
using std::vector;
using std::map;
using std::pair;
using boost::shared_ptr;

struct eMu {

  eMu(const string &directoryTag="DY_m10+pr+a05+o03+pr_4680pb");
  virtual ~eMu(){};
  int run();
  bool doDMDY; ///< flag to look at rapidity in mass bins
  bool doPUreWeight;///< flag to run pu reweighting reading values from histogram 
  bool saveToRootFile;///< flag to save output to ROOT file
  bool verbose;///< flag to set verbose option
  string emuNtupleDir; ///< eMu Ntuple directory
  string eeNtupleDir; ///< ee Ntuple directory
  string subDir; ///< sub directory for both eMu and ee
  string dirTag; ///< directory tag where selected events are placed
  string filePostfix; ///<
  string lumiVal;
  TString outFileName;
  double reWeight; ///< reWeight MC rather than having to rerun
  map<string, const char*> dataMap;
  map<string, const char*> eebkgMap;
  map<string, const char*> emubkgMap;
  float xmax; ///< maximum emu and ee mass bin value
  float xmin; ///< minimum emu and ee mass bin value
  unsigned int nBins; ///< number of mass histogram bins
  float rap_Max; ///< maximum emu and ee rapidity bin value
  float rap_Min; ///< minimum emu and ee rapidity bin value
  unsigned int rap_Bins; ///< number of rapidity histogram bins

  unsigned int d_array_size;
  unsigned int emu_array_size;
  unsigned int emu_bkgarray_size;

  //TVectorT<double> true2eBackgroundFromData;
  //TVectorT<double> true2eBackgroundFromDataError;
  //TVectorT<double> true2eBackgroundFromDataErrorSyst;

  //double massBins2D[];

  TMatrixT<double> true2eBackgroundFromData;
  TMatrixT<double> true2eBackgroundFromDataError;
  TMatrixT<double> true2eBackgroundFromDataErrorSyst;

  void doDMDYanal(const bool& runAnal) {doDMDY = runAnal;}
  void runPUreWeighting(const bool& runReWeighting) {doPUreWeight  = runReWeighting;}
  void setSaveToRootFile(const bool& saveFlag) {saveToRootFile = saveFlag;}
  void setMCreWeight(const double& reWeightVal) { reWeight = reWeightVal;}
  void setVerbose(const bool& verboseFlag) {verbose = verboseFlag;}
  TH1F* subtractEMubackground3(TH1F *inputHisto, vector<vector<shared_ptr<TH1F> >* >& sHist, vector<vector<shared_ptr<TH1F> >* >& statHist);
  double calcError(const double& eMuObs, vector<vector<shared_ptr<TH1F> >* >& vvHist);
  double calcError(TH1F *inputHisto, vector<vector<shared_ptr<TH1F> >* >& vvHist);

};

class Fill_array_from_map_functor {
public:
  Fill_array_from_map_functor(string stringArray[]): particle_name_array(stringArray){}
  //need to implement something that checks the size of the array and doesn't overrun
  void operator()(pair<string, const char*> inputPair){
    //    particle_name_array[index] = inputMap->first;
    *(particle_name_array) = inputPair.first;
    ++particle_name_array;
    //if (index < (aSize -1)) ++index; //make sure index doesn't go above array size
  }  
  string* particle_name_array;
};


class Fill_vector_from_map_functor {
public:
  Fill_vector_from_map_functor(vector<const char *> & vTemp): vMod(vTemp){}
  void operator()(pair<string, const char*> inputPair){
    vMod.push_back(inputPair.second);
  }  
  vector<const char *>& vMod;
};

#endif
