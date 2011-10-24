#!/bin/sh

cd ../Include
rm *.so
cd ../FullChain

filename_data='"../config_files/data.conf"'
filename_mc='"../config_files/summer11mc.input"'
filename_cs='"../config_files/xsecCalc.conf"'

#Selection
statusSelection=OK
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "WILL DO: plotDY(\"../config_files/data.conf\")"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
cd ../Selection
rm *.so
echo
root -b -q -l plotDY.C+\($filename_data\)
if [ $? != 0 ]; then 
   statusSelection=FAILED
fi
cd ../FullChain
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "DONE: plotDY(\"../config_files/data.conf\")"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"

#PlotSelectDY
statusPlotSelectDY=OK
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "WILL DO: plotSelectDY.C(\"../config_files/data.conf\")"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
cd ../MassSpectrum
rm *.so
echo
root -b -q -l plotSelectDY.C+\($filename_data\)
if [ $? != 0 ]; then 
   statusPlotSelectDY=FAILED
fi
cd ../FullChain
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "DONE: plotSelectDY.C(\"../config_files/data.conf\")"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"

#SubtractBackground
statusSubtractBackground=OK
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "WILL DO: subtractBackground.C(\"../config_files/data.conf\")"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
cd ../MassSpectrum
rm *.so
echo
root -b -q -l subtractBackground.C+\($filename_data\)
if [ $? != 0 ]; then 
   statusSubtractBackground=FAILED
fi
cd ../FullChain
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "DONE: subtractBackground.C(\"../config_files/data.conf\")"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"

#Unfolding
statusUnfolding=OK
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "WILL DO: plotDYUnfoldingMatrix(\"../config_files/summer11mc.input\")"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
cd ../Unfolding
rm *.so
echo
root -b -q -l plotDYUnfoldingMatrix.C+\($filename_mc\)
if [ $? != 0 ]; then 
   statusUnfolding=FAILED
fi
cd ../FullChain
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "DONE: plotDYUnfoldingMatrix(\"../config_files/summer11mc.input\")"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"

#Unfolding Systematics
statusUnfoldingSyst=OK
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "WILL DO: plotDYUnfoldingMatrix(\"../config_files/summer11mc.input\")"
echo "in systematics mode"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
cd ../Unfolding
rm *.so
echo
root -b -q -l auxScriptSyst.C
if [ $? != 0 ]; then 
   statusUnfoldingSyst=FAILED
fi 
cd ../FullChain
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "DONE: plotDYUnfoldingMatrix(\"../config_files/summer11mc.input\")"
echo "in systematics mode"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"

#Acceptance
statusAcceptance=OK
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "WILL DO: plotDYAcceptance(\"../config_files/summer11mc.input\")"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
cd ../Acceptance
rm *.so
echo
root -b -q -l plotDYAcceptance.C+\($filename_mc\)
if [ $? != 0 ]; then 
   statusAcceptance=FAILED
fi
cd ../FullChain
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "DONE: plotDYAcceptance(\"../config_files/summer11mc.input\")"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"

#Acceptance Systematics
statusAcceptanceSyst=OK
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "WILL DO: plotDYAcceptance(\"../config_files/summer11mc.input\")"
echo "in systematics mode"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
cd ../Acceptance
rm *.so
echo
root -b -q -l auxScriptSyst.C
if [ $? != 0 ]; then 
   statusUnfoldingSyst=FAILED
fi 
cd ../FullChain
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "DONE: plotDYAcceptance(\"../config_files/summer11mc.input\")"
echo "in systematics mode"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"

#Efficiency
statusEfficiency=OK
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "WILL DO: plotDYEfficiency(\"../config_files/summer11mc.input\")"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
cd ../Efficiency
rm *.so
echo
root -b -q -l plotDYEfficiency.C+\($filename_mc\)
if [ $? != 0 ]; then 
   statusEfficiency=FAILED
fi
cd ../FullChain
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "DONE: plotDYEfficiency(\"../config_files/summer11mc.input\")"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"

#Efficiency Scale Factors
statusEfficiencyScaleFactors=OK
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "WILL DO: seven scripts from EfficiencyScaleFactors"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
cd ../EfficiencyScaleFactors
rm *.so
echo
root -b -q -l auxScript.C+
if [ $? != 0 ]; then 
   statusEfficiencyScaleFactors=FAILED
fi
cd ../FullChain
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "DONE: seven scripts from EfficiencyScaleFactors"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"

#PlotDYFSRCorrections
statusPlotDYFSRCorrections=OK
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "WILL DO: plotDYFSRCorrections(\"../config_files/summer11mc.input\")"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
cd ../Fsr
rm *.so
echo
root -b -q -l plotDYFSRCorrections.C+\($filename_mc\)
if [ $? != 0 ]; then 
   statusPlotDYFSRCorrections=FAILED
fi
cd ../FullChain
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "DONE: plotDYFSRCorrections(\"../config_files/summer11mc.input\")"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"

#PlotDYFSRCorrectionsSansAcc
statusPlotDYFSRCorrectionsSansAcc=OK
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "WILL DO: plotDYFSRCorrectionsSansAcc(\"../config_files/summer11mc.input\")"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
cd ../Fsr
rm *.so
echo
root -b -q -l plotDYFSRCorrectionsSansAcc.C+\($filename_mc\)
if [ $? != 0 ]; then 
   statusPlotDYFSRCorrectionsSansAcc=FAILED
fi
cd ../FullChain
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "DONE: plotDYFSRCorrectionsSansAcc(\"../config_files/summer11mc.input\")"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"

#TheoryErrors
statusTheoryErrors=OK
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "WILL DO: TheoryErrors(\"../config_files/summer11mc.input\")"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
cd ../Theory
rm *.so
echo
root -b -q -l TheoryErrors.C+\($filename_mc\)
if [ $? != 0 ]; then 
   statusTheoryErrors=FAILED
fi
cd ../FullChain
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "DONE: TheoryErrors(\"../config_files/summer11mc.input\")"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"

#CrossSection
statusCrossSection=OK
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "WILL DO: CrossSection(\"../config_files/xsecCalc.conf\")"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
cd ../CrossSection
rm *.so
echo
root -b -q -l CrossSection.C+\($filename_cs\)
if [ $? != 0 ]; then 
   statusCrossSection=FAILED
fi
cd ../FullChain
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "DONE: CrossSection(\"../config_files/xsecCalc.conf\")"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"

echo "Full chain summary:"
echo "              Selection:    " $statusSelection
echo "           PlotSelectDY:    " $statusPlotSelectDY
echo "     SubtractBackground:    " $statusSubtractBackground
echo "              Unfolding:    " $statusUnfolding
echo "         Syst Unfolding:    " $statusUnfoldingSyst
echo "             Acceptance:    " $statusAcceptance
echo "        Syst Acceptance:    " $statusAcceptanceSyst
echo "             Efficiency:    " $statusEfficiency
echo " EfficiencyScaleFactors:    " $statusEfficiencyScaleFactors
echo "         FSRCorrections:    " $statusPlotDYFSRCorrections
echo "  FSRCorrectionsSansAcc:    " $statusPlotDYFSRCorrectionsSansAcc
echo "           TheoryErrors:    " $statusTheoryErrors
echo "           CrossSection:    " $statusCrossSection

