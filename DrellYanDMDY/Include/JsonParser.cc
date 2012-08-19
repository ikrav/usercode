#include "JsonParser.hh"
#include <TClass.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctype.h>

using namespace std;

ClassImp(JsonParser)

JsonParser::JsonParser(){

  Reset();

}

void JsonParser::Initialize(TString filename){

  Reset();

  ifstream jsonFile(filename.Data());
  if (!jsonFile) {
    cout << "JsonParser::Error: Unable to open JSON file";
    assert(0);
  }
  
  char x;

  // First char should be opening curly bracket
  jsonFile >> x;
  if( x != '{') {
    cout << "JsonParser: ERROR: parse error 1" << endl;
    return;
  }

  bool readingRunNumber = false;
  vector<int> runNumberDigits;
  runNumberDigits.clear();
  int runNumber = 0;
  int depth = 0;
  vector<int> firstLumiDigits;
  firstLumiDigits.clear();
  int firstLumi = 0;
  vector<int> secondLumiDigits;
  secondLumiDigits.clear();
  int secondLumi = 0;
  vector<int> allRunLumis;
  allRunLumis.clear();
  while ( ! jsonFile.eof() ) {
    
    jsonFile >> x;
    if( x == '"' ) {
      // Double quotes mean beginning or ending of the sequnce of digits
      // for run number
      if( readingRunNumber ) {
	readingRunNumber = false;
	// Assemble the right number. (Remember, we have read digits in the reverse order)
	runNumber = AssembleNumber(runNumberDigits);
	_runList.push_back(runNumber);
	allRunLumis.clear();
// 	cout << "Run number " << runNumber << endl;
      }else{
	readingRunNumber = true;
	runNumberDigits.clear();
      }
    } else if( readingRunNumber && isdigit(x) ){
      // We should be now reading digits of the run number
      int digit = x - '0'; // convert into int
      runNumberDigits.push_back(digit);
    } else if( readingRunNumber && ! isdigit(x) ){
      // only digits should be found between the double quotes
      cout << "JsonParser: ERROR: parse error 2" << endl;
      return;
    } else if( x == ':' ){
      // We are about to start reading lumi sections for this run
      if (depth != d_outside ){
	cout << "JsonParser: ERROR: parse error 3" << endl;
	return;
      }
    } else if( x == '[' && depth == d_outside) {
      depth = d_insideOuterBracket;
    } else if( x == '[' && depth == d_insideOuterBracket) {
      depth = d_insideInnerBracketNumberOne;
      firstLumiDigits.clear();
      firstLumi = 0;
    } else if( depth == d_insideInnerBracketNumberOne && isdigit(x) ) {
     // read digits of the starting lumi section number
      int digit = x - '0';
      firstLumiDigits.push_back(digit);
    } else if( depth == d_insideInnerBracketNumberOne && !isdigit(x) && !(x==',') ){
      cout << "JsonParser: ERROR: parse error 4" << endl;
      return;
    } else if( depth == d_insideInnerBracketNumberOne && x == ',' ){
      // The starting lumi section is completed
      depth = d_insideInnerBracketNumberTwo;
      firstLumi = AssembleNumber(firstLumiDigits);
      firstLumiDigits.clear();
      secondLumiDigits.clear();
      secondLumi = 0;
    } else if( depth == d_insideInnerBracketNumberTwo && isdigit(x) ){
      // read digits of the ending lumi section
      int digit = x - '0';
      secondLumiDigits.push_back(digit);
    } else if( depth == d_insideInnerBracketNumberTwo && !isdigit(x) && !(x==']') ){
      cout << "JsonParser: ERROR: parse error 5" << endl;
      return;
    } else if( depth == d_insideInnerBracketNumberTwo && x==']' ){
      // The starting lumi section is completed
      secondLumi = AssembleNumber(secondLumiDigits);
      firstLumiDigits.clear();
      secondLumiDigits.clear();
      depth = d_insideOuterBracket;
      // Now we have starting and ending lumi for this lumi block. Put all relevant
      // lumis into the list of lumis for this run
      for(int ilumi = firstLumi; ilumi <= secondLumi; ilumi++)
	allRunLumis.push_back(ilumi);
//       cout << "       lumi " << firstLumi << "    " << secondLumi << endl;
    } else if( depth == d_insideOuterBracket && x==','){
      // the next range of lumis is about to start, do nothing
    } else if( depth == d_insideOuterBracket && x==']' ){
      // Completed all blocks of lumis for this run
      _lumiStatusList.push_back(allRunLumis);
      depth = d_outside;
    } else if( depth == d_outside && x==',' ){
      // about to move to the next run block
    } else if( x == '}') {
      // End of JSON data is a closing curly bracket
      if ( x == '}' ) break;
    } else {
      cout << "JsonParser: ERROR: parse error 6" << endl;
      return;
    }
//     cout << ">> "<< x << " <<" << endl;
    
  }
  jsonFile.close();
  
  _isInitialized = true;
  
}

bool JsonParser::HasRunLumi(int run, int lumi){

  bool result = false;

  if( ! _isInitialized ) {
    cout << "JsonParser::ERROR: attempt to use without initialization" << endl;
    return result;
  }

  if( _runList.size() != _lumiStatusList.size() ){
    cout << "JsonParser::ERROR: internal bookkeeping inconsistency" << endl;
    return result;
  }
  
  vector < vector < int > >::iterator iter_irun;
  vector<int>        ::iterator iter_ilumi;
  
  int irun = -1;
  for(iter_irun = _lumiStatusList.begin(); iter_irun != _lumiStatusList.end(); iter_irun++){
    irun++;
    if( result ) break;
    if ( _runList[irun] != run ) continue;
    for(iter_ilumi = (*iter_irun).begin(); iter_ilumi != (*iter_irun).end(); iter_ilumi++){
      if( (*iter_ilumi) == lumi ) {
	result = true;
	break;
      }
    }
  }

  return result;
}

void JsonParser::Reset(){

  _runList.clear();
  _lumiStatusList.clear();
  _isInitialized = false;

}

int JsonParser::AssembleNumber(vector<int> data){

  int result = 0;
  for(UInt_t i = 0; i < data.size(); i++){
    result += data[i] * (int)pow(10,data.size()-i-1);
  }
  return result;
  
 }

void JsonParser::Print(){

  if( ! _isInitialized ) {
    cout << "JsonParser::ERROR: attempt to use without initialization" << endl;
    return;
  }

  if( _runList.size() != _lumiStatusList.size() ){
    cout << "JsonParser::ERROR: internal bookkeeping inconsistency" << endl;
    return;
  }
  
  vector < vector < int > >::iterator iter_irun;
  vector<int>        ::iterator iter_ilumi;
  
  int irun = -1;
  for(iter_irun = _lumiStatusList.begin(); iter_irun != _lumiStatusList.end(); iter_irun++){
    irun++;
    printf("Run %d\n", _runList[irun]);
    printf("  lumis: ");
    bool printedDash = false;
    bool previousNotPrinted = false;
    for(iter_ilumi = (*iter_irun).begin(); iter_ilumi != (*iter_irun).end(); iter_ilumi++){
      bool first = ( iter_ilumi == (*iter_irun).begin() );
      bool last  = ( iter_ilumi == (*iter_irun).end() - 1 );
      if( ! ( (*(iter_ilumi-1)) == (*iter_ilumi) - 1) || first || last ){
	if(previousNotPrinted && !last) {
	  printf(" %d", (*(iter_ilumi-1)));
	}
	printf("  %d", (*iter_ilumi));
	printedDash = false;
	previousNotPrinted = false;
      } else {
	previousNotPrinted = true;
	if( !printedDash){
	  printf(" -");
	  printedDash = true;
	}
      }
//       prev = (*iter_ilumi);
    }
    printf("\n");
  }

}

