#ifndef CMDLINEOPTS_HH
#define CMDLINEOPTS_HH

#include <vector>
#include <string>
#include <map>

using std::string;
using std::vector;
using std::map;


/*!
 * Class for accepting and processing command line options 
*/


class CmdLineOpts{

public:

  CmdLineOpts(int, char**);/*!< Class constructor*/

  void addOption(const string& optionName, bool &flag, string message="Option implemented");
  void addOption(const string& optionName, int &intValue, string message="Option implemented");
  void addOption(const string& optionName, float &floatValue, string message="Option implemented");
  void addOption(const string& optionName, double &doubleValue, string message="Option implemented");
  void addOption(const string& optionName, string &stringValue, string message="Option implemented");

  void readCmdLine();
  //~CmdLineOpts();/*!< Class destructor*/

private:

  struct CmdVar{

  public:

    CmdVar(const string& optionName, bool& flag, string type, string message);
    CmdVar(const string& optionName, int& intValue, string type, string message);
    CmdVar(const string& optionName, float& floatValue, string type, string message);
    CmdVar(const string& optionName, double& doubleValue, string type, string message);
    CmdVar(const string& optionName, string& stringValue, string type, string message);

    //bool* getBool(){return flagPtr_;}
    //int* getInt(){return intPtr_;}
    //float* getFloat(){return floatPtr_;}
    //double* getDouble(){return doublePtr_;}

    string optionName_;
    bool *flagPtr_;
    int *intPtr_;
    float *floatPtr_;
    double *doublePtr_;    
    string *stringPtr_;    
    string type_;
    string message_;


  };
  
  map<string,CmdVar> cmdMap_;
  vector<string> message_list;
  vector<string> unrequestedOption_list;
  string progName_; 
  int argc_;
  char **argv_;

  void help() const;




}; 

#endif
