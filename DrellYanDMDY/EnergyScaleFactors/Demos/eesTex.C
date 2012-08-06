
#include "../Include/ElectronEnergyScale.hh"



void eesTex(const TString &electronEnergyScaleString) {
  std::cout << "\nPrinting ElectronEnergyScale correction parameters as a tex table\n";
  TString eesString=electronEnergyScaleString;

  ElectronEnergyScale esc(electronEnergyScaleString);
  TString fname=electronEnergyScaleString;
  fname.Append(".tex");
  esc.printAsTexTable(fname);
  return;
}
