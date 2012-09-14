#include "CmdLineOpts.hh"
#include <iostream>
#include <sstream>
#include <utility>

using std::cout;
using std::cerr;
using std::istringstream;
using std::make_pair;

CmdLineOpts::CmdLineOpts(int argc, char** argv):argc_(argc), argv_(argv){
  progName_ = argv[0];
}


void CmdLineOpts::addOption(const string& optionName, bool& flag, string message){  
  string type = "bool";
  CmdVar inVar(optionName, flag, type, message);
  cmdMap_.insert(make_pair(optionName,inVar));
}

void CmdLineOpts::addOption(const string& optionName, int& intValue, string message){
  string type = "int";
  CmdVar inVar(optionName, intValue, type ,message);
  cmdMap_.insert(make_pair(optionName,inVar));
}

void CmdLineOpts::addOption(const string& optionName, float& floatValue, string message){
  string type = "float";
  CmdVar inVar(optionName, floatValue, type ,message);
  cmdMap_.insert(make_pair(optionName,inVar));
}

void CmdLineOpts::addOption(const string& optionName, double& doubleValue, string message){
  string type = "double";
  CmdVar inVar(optionName, doubleValue, type ,message);
  cmdMap_.insert(make_pair(optionName,inVar));
}

void CmdLineOpts::addOption(const string& optionName, string& stringValue, string message){
  string type = "string";
  CmdVar inVar(optionName, stringValue, type ,message);
  cmdMap_.insert(make_pair(optionName,inVar));
}

void CmdLineOpts::readCmdLine(){
  for (int i=1; i<argc_; ++i){
    string inputArg = argv_[i];
    if ((inputArg == "--help") || (inputArg == "--Help")) help();
    if (inputArg.find("--") == 0) {// check that "--" is at the start of the option
      if (cmdMap_.count(inputArg)){ // count will return 1 if the key exists and 0 if it doesn't
	map<string, CmdVar>::iterator map_iter = cmdMap_.find(inputArg);
	if (map_iter->second.type_ == "bool"){
	  *(map_iter->second.flagPtr_) = true;
          message_list.push_back(map_iter->second.message_);
	}
        if (map_iter->second.type_ == "int"){
	  string inputVal = argv_[++i];
	  size_t pos = inputVal.find("--");// check that option isn't followed by another option. Returns npos if there is no match
	  if (pos != string::npos){ 
            std::cerr << "Need to specify number after argument and not " << inputVal  << "\n";
	    //throw "Need to specify number after argument";
	  }
          istringstream issValue(inputVal);	
	  issValue >> *(map_iter->second.intPtr_);
          message_list.push_back(map_iter->second.message_);
	}

        if (map_iter->second.type_ == "float"){
	  string inputVal = argv_[++i];
          size_t pos = inputVal.find("--");// check that option isn't followed by another option. Returns npos if there is no match
	  if (pos != string::npos){ 
	    std::cerr << "Need to specify number after argument and not " << inputVal  << "\n";
	    //throw "Need to specify number after argument";
	  }
          istringstream issValue(inputVal);	
	  issValue >> *(map_iter->second.floatPtr_);
          message_list.push_back(map_iter->second.message_);
	}
        
        if (map_iter->second.type_ == "double"){
	  string inputVal = argv_[++i];
          size_t pos = inputVal.find("--");// check that option isn't followed by another option. Returns npos if there is no match
	  if (pos != string::npos){ 
            std::cerr << "Need to specify number after argument and not " <<  inputVal  << "\n";
	    //throw "Need to specify number after argument";
	  }
	  istringstream issValue(inputVal);	
	  issValue >> *(map_iter->second.doublePtr_);
          message_list.push_back(map_iter->second.message_);
	}

        if (map_iter->second.type_ == "string"){
	  string inputVal = argv_[++i];
          size_t pos = inputVal.find("--");// check that option isn't followed by another option. Returns npos if there is no match
	  if (pos != string::npos){ 
            std::cerr << "Need to specify number after argument and not " <<  inputVal  << "\n";
	    //throw "Need to specify number after argument";
	  }
	  istringstream issValue(inputVal);	
	  issValue >> *(map_iter->second.stringPtr_);
          message_list.push_back(map_iter->second.message_);
	}
      } else {
        unrequestedOption_list.push_back(inputArg);
      }
    }
  }

  //Print out messages
  cout << "\nRunning program " << progName_ << " \n";
  cout  <<"Detected " << argc_-1 << " arguments\n\n";
  for (unsigned int i=0; i<message_list.size(); i++){
    cout << message_list[i] << "\n";
  }
  cout << "\n";
 
  for (unsigned int i=0; i<unrequestedOption_list.size();++i){
    cerr << "Not requested option " << unrequestedOption_list[i]  << "; option ignored\n";
  }
  cout << "\n";

}

void CmdLineOpts::help() const{
  //Need to be able to add to this from one's own main program
  cout <<"\nUsage: " << progName_<< " [OPTIONS]\n\n"
       <<"OPTIONS have to be specified with -- for example --totEvts.\n"
       <<"The following value is taken and assigned to a variable.\n" 
       <<" The variable to be assigned has to be specified in your main\n"
       << "code with the CmdLineOpts::addOption method. Values can be assigned\n"
       <<"to  the types int, float and double. bools can be set to true where\n"
       <<" no value is needed after the -- option\n\n"
       <<"IMPORTANT: must run CmdLineOpts.readCmdLine() in main prog\n"
       <<"to process arguments\n"
       <<"IMPORTANT: options are not checked for duplications\n"
       <<"which options is read in this case is not specified\n";

  throw;
}

CmdLineOpts::CmdVar::CmdVar(const string& optionName, bool& flag, string type, string message):
  optionName_(optionName), flagPtr_(&flag),type_(type),message_(message){
}
CmdLineOpts::CmdVar::CmdVar(const string& optionName, int& intValue, string type, string message):
 optionName_(optionName), intPtr_(&intValue),type_(type),message_(message){
}
CmdLineOpts::CmdVar::CmdVar(const string& optionName, float& floatValue, string type, string message):
  optionName_(optionName), floatPtr_(&floatValue),type_(type),message_(message){
}
CmdLineOpts::CmdVar::CmdVar(const string& optionName, double& doubleValue, string type, string message):
  optionName_(optionName), doublePtr_(&doubleValue),type_(type),message_(message){
}
CmdLineOpts::CmdVar::CmdVar(const string& optionName, string& stringValue, string type, string message):
  optionName_(optionName), stringPtr_(&stringValue),type_(type),message_(message){
}




