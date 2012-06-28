#include "eMu.hh"
#include "CmdLineOpts.hh"

int main(int argc, char** argv){

  CmdLineOpts myOpts(argc, argv); // read in command line options

  double reWeight(1.0);//rather than running the MC everytime a new weight is required, just use the reweight variable
  bool doPUreWeight(false);//Do pileup reweighting
  bool doDMDY(false); //run 2D rapidity analysis
  bool saveRootFile(false);// save output to a root file
  bool activateVerbose(false);

  myOpts.addOption("--doDMDY", doDMDY,
                   "--doDMDY: Will run 2D rapidity analysis, need dilepton rapidity value in the ntuple");
  myOpts.addOption("--doPUreWeight", doPUreWeight,
		   "--doPUreWeight: Will run pileUpreWeighting, need pile-up weights in the ntuple");
  myOpts.addOption("--reWeight", reWeight,
		   "--reWeight: Will reweight events by specified factor");
  myOpts.addOption("--saveRootFile", saveRootFile,
		   "--saveRootFile: Will save histos etc to root file");
  myOpts.addOption("--verbose", activateVerbose,
                   "--verbose: Will turn on verbosity");

  myOpts.readCmdLine();// process the command line options

  eMu emuMethod;

  //set options
  emuMethod.doDMDYanal(doDMDY);
  emuMethod.runPUreWeighting(doPUreWeight);
  emuMethod.setSaveToRootFile(saveRootFile);
  emuMethod.setMCreWeight(reWeight);
  emuMethod.setVerbose(activateVerbose);

  //run the program
  emuMethod.run();
  return 0;
}
