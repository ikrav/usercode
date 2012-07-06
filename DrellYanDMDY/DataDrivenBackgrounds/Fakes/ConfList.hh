#ifndef CONFLIST_HH
#define CONFLIST_HH

#include "TreeQueue.hh"

//===============================================
//conflist class
//=================================================


class ConfList : public TreeQueue{
  public:
  ConfList(string listOfRootFiles);
  bool nextFile(); //method to jump to next file
  void reset();
  double getVal() const {return *current_xsec;}

  private:
  vector<double> xseclist;
  vector<double>::iterator xsec_iter;
  bool storeFilenames();
  vector<double>::iterator current_xsec;
};

#endif

