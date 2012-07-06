
#include "ConfList.hh"
#include <sstream>

using std::istringstream;

//================================
//Conflist definitions
//================================


ConfList::ConfList(string listOfRootFiles): TreeQueue(listOfRootFiles, true){
  if (storeFilenames()){//anyfilenames stored?
    reset();//Then point to the first file
  } else {
    //should I throw an exception
    throw "empty list of files";
  }
}

bool ConfList::storeFilenames(){  
  string line;
  xseclist.clear();//clear vector
  filelist.clear();//clear vector
  while (rootFileStream.good()){
    getline(rootFileStream,line);
    if (line[0] != '#'){//ignore comment lines
      if (line.size() != 0) {//ignore empty lines
	istringstream concat(line);
	string rootfile;
	double xsec;
	concat >> rootfile >> xsec;
	filelist.push_back(rootfile);
        xseclist.push_back(xsec);
      }
    }
  }
  //need to clear file for reuse
  rootFileStream.clear();

  if (filelist.size() != xseclist.size()){
    throw "size mismatch";
  }

  if (!filelist.empty()) {//same as checking if xseclist is empty
    return true;
  } else {
    return false;
  }
}

bool ConfList::nextFile(){
  bool outcome;
  if (TreeQueue::nextFile()){
    if (xsec_iter != xseclist.end()){
      current_xsec = xsec_iter;
      ++xsec_iter;
      outcome = true;
    } else {
      outcome = false;
    }
  } else {
    outcome = false;
  }

  return outcome;
}

void ConfList::reset(){
  TreeQueue::reset();
  xsec_iter = xseclist.begin();
  current_xsec =  xsec_iter;
}
