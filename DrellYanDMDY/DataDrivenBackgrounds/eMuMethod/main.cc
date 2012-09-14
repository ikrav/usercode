#include "eMu.hh"
#include "CmdLineOpts.hh"
#include <iostream>

int main(int argc, char** argv){

  CmdLineOpts myOpts(argc, argv); // read in command line options

  double reWeight(1.0);//rather than running the MC everytime a new weight is required, just use the reweight variable
  bool doPUreWeight(false);//Do pileup reweighting
  bool doDMDY(false); //run 2D rapidity analysis
  bool saveRootFile(false);// save output to a root file
  bool activateVerbose(false);
  string dirTag="DY_m10+pr+a05+o03+pr_4680pb";

  myOpts.addOption("--doDMDY", doDMDY,
                   "--doDMDY: Will run 2D rapidity analysis, need dilepton rapidity value in the ntuple");
  myOpts.addOption("--doPUreWeight", doPUreWeight,
		   "--doPUreWeight: Will run pileUpreWeighting, need pile-up weights in the ntuple");
  myOpts.addOption("--reWeight", reWeight,
		   "--reWeight: Will reweight events by specified factor");
  myOpts.addOption("--saveRootFile", saveRootFile,
		   "--saveRootFile: Will save histos etc to root file");
  myOpts.addOption("--dirTag",dirTag,
		   "--dirTag: directory tag of a form DY_....pb");
  myOpts.addOption("--verbose", activateVerbose,
                   "--verbose: Will turn on verbosity");

  myOpts.readCmdLine();// process the command line options

  eMu emuMethod(dirTag);

  //set options
  emuMethod.doDMDYanal(doDMDY);
  emuMethod.runPUreWeighting(doPUreWeight);
  emuMethod.setSaveToRootFile(saveRootFile);
  emuMethod.setMCreWeight(reWeight);
  emuMethod.setVerbose(activateVerbose);

  //run the program
  int res=emuMethod.run();
  if (res!=1) std::cout << "error from emuMethod.run. code=" << res << "\n";
  return res;
}
