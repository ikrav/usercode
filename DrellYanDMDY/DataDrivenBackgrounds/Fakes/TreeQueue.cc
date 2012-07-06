
#include "TreeQueue.hh"
#include <iostream>                 // stl iostream
#include <assert.h>                 // assert

using std::cout;

//=================
//Class definitions
//=================

TreeQueue::TreeQueue(string listOfRootFiles){ 
  rootFileStream.open(listOfRootFiles.c_str(), ifstream::in);
  //Need to test if file exists. If it doesn't say so in an error
  if (storeFilenames()){//anyfilenames stored?
    reset();//Then point to the first file
  } else {
    //should I throw an exception
    throw "empty list of files";
  }
}

TreeQueue::TreeQueue(string listOfRootFiles, bool isInherited){
  //partially construct members as Derived class will do the rest 
  rootFileStream.open(listOfRootFiles.c_str(), ifstream::in);
  isInherited = true;
}


TreeQueue::~TreeQueue(){
  //delete infile; memory leak, need to implement this properly!
}

bool TreeQueue::nextFile(){
  if (s_iter != filelist.end()){
    infile = new TFile(TString(*s_iter));
    assert(infile);
    current = s_iter;
    ++s_iter;
    return true;
  } else {
    return false;
  }
}

void TreeQueue::reset(){
  s_iter = filelist.begin();
  infile = new TFile(TString(*s_iter));  
  assert(infile);
  current = s_iter;
}

TTree* TreeQueue::getTree(const TString& treeName){
  // Get the TTree
  evtTree = (TTree*)infile->Get(treeName); 
  assert(evtTree);
  return evtTree;
}

bool TreeQueue::storeFilenames(){  
  string line;
  filelist.clear();//clear vector
  while (rootFileStream.good()){
    getline(rootFileStream,line);
    if (line[0] != '#'){//ignore comment lines
      if (line.size() != 0) {//ignore empty lines
        filelist.push_back(line);
      }
    }
  }
  //need to clear file in order to reuse
  rootFileStream.clear();

  if (!filelist.empty()) {
    return true;
  } else {
    return false;
  }
}

void TreeQueue::printOut() const{
  for (vector<string>::const_iterator iter = filelist.begin(); iter != filelist.end(); ++iter){
    cout << *iter << "\n";
  }

  throw "End of printout \n";
}

void TreeQueue::status() const{
  cout << "Currently pointing to: " << *current << "\n";
}
