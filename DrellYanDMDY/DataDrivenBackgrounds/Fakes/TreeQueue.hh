#ifndef TREEQUEUE_HH
#define TREEQUEUE_HH

#include <string>
#include <vector>
#include <fstream>                  // file stream

#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TString.h>                // ROOT string class

using std::vector;
using std::string;
using std::fstream;


//===============================================================================
//Class to read a txt file containing a list of root files
//and load in the ntuples
//==============================================================================

class TreeQueue{
  public:
  TreeQueue(){};
  TreeQueue(string listOfRootFiles);
  virtual ~TreeQueue();
  TTree* getTree(const TString& treeName);
  virtual bool nextFile(); //method to jump to next file
  virtual void reset();
  void printOut() const;
  void status() const;

  protected:
  TreeQueue(string listOfRootFiles,  bool isInherited);
  vector<string> filelist;
  vector<string>::iterator s_iter;
  vector<string>::iterator current;
  ifstream rootFileStream;

  private:
  TTree* evtTree;
  TFile* infile;//should definitely use smart pointers here, currently have a memory leak  
  virtual bool storeFilenames();
  
};

#endif

