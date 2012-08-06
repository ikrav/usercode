#include "example_run_fit_20120802.C"
#include "example_run_rec_20120802.C"

void loop_fit() {
  std::vector<std::string> etaRanges,models;
  fillCases(etaRanges,models);

  if (etaRanges.size()!=models.size()) {
    std::cout << "etaRanges.size=" << etaRanges.size() << ", models.size=" << models.size() << "\n";
    return;
  }

  for (unsigned int i=0; i<etaRanges.size(); ++i) {
    example_fit_run(etaRanges[i].c_str(),models[i].c_str());
    example_rec_run(etaRanges[i].c_str(),models[i].c_str());
  }

  return;

}

