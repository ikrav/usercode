#include "eMu_pyConf.hh"

#include <iostream>

using std::cout;

///Constructor
eMu_If::eMu_If() : eMu(){}

///Help information for members
int eMu_If::help(){
  cout << "Class members:\n"
       << "\t\temuNtupleDir (sets the emu directory)\n"
       << "\t\teeNtupleDir (sets the emu directory)\n"
       << "\t\tfilePostfix (NEED TO CHANGE IMPLEMENTATION)\n"
       << "\t\thelp (prints this help output for class members)\n"
       << "\t\tlumiVal (luminosity value to write on plots)\n"
       << "\t\tnBins (number of xbins in plots)\n"
       << "\t\tsubDir (sets the sub directory under emu and ee dirs)\n"
       << "\t\txmax (maximum xvalue on plots)\n"
       << "\t\txmin (minimum xvalue on plots)\n"
       << "\n\nClass member functions:\n"
       << "\t\tdoDMDY (bool) (run 2D plots)\n"
       << "\t\trunPUreWeighting(bool) (apply pile up reweighting)\n"
       << "\t\tsetMCreWeighting(double) (apply specified weight to MC events)\n"
       << "\t\tsaveRootFile(bool) (save output to ROOT file)\n";
  return 0;
}
