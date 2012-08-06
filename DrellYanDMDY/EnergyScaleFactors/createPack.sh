date="`date +%Y%m%d`"
dir="dir-Pack${date}/"

files="ElectronEnergyScaleAdv.cc ElectronEnergyScaleAdv.hh ElectronEnergyScaleAdv.C"
files="${files} HelpingTools.hh HelpingTools.cc HelpingTools.C"
files="${files} MyFitModels.hh MyFitModels.cc"
files="${files} example_run*.C"
files="${files} fitMassScaleFactorsE.cc fitMassScaleFactorsE.hh fitMassScaleFactorsE.C"
files="${files} createPack.sh Test.hh CPlotMdf.hh CPlotMdf.cc"
#files="${files} dir-TestPackage-*Gauss-*/dummy"
files="${files} Demos/eesCorrectionDemo.C Demos/eesSmearEtaEtaBinDemo.C Demos/eesSmearEventDemo.C Demos/eesTex.C Demos/plotTwoESFSetsDemo.C Demos/rootlogon.C"

tar cvfz ~/pack${date}.tz ${files}
