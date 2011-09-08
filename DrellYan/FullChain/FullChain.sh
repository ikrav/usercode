#!/bin/sh

cd ../Include
rm *.so
cd ../FullChain

echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "WILL DO: plotDY(\"../config_files/data.conf\")"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
cd ../Selection
rm *.so
#root -b -q -l auxScript.C+
#still crashes for some reason
cd ../FullChain
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "DONE: plotDY(\"../config_files/data.conf\")"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"

echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "WILL DO: plotSelectDY.C(\"../config_files/data.conf\")"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
cd ../MassSpectrum
rm *.so
root -b -q -l auxScript1.C+
cd ../FullChain
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "DONE: plotSelectDY.C(\"../config_files/data.conf\")"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"

echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "WILL DO: subtractBackground(\"../config_files/data.conf\")"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
cd ../MassSpectrum
rm *.so
root -b -q -l auxScript2.C+
cd ../FullChain
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "DONE: subtractBackground(\"../config_files/data.conf\")"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"

echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "WILL DO: plotDYUnfoldingMatrix(\"../config_files/summer11mc.input\")"
echo "20 times; and calcUnfoldingSystematics"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
cd ../Unfolding
rm *.so
root -b -q -l auxScript.C+
cd ../FullChain
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "DONE: plotDYUnfoldingMatrix(\"../config_files/summer11mc.input\")"
echo "20 times; and calcUnfoldingSystematics"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"

echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "WILL DO: plotDYAcceptance(\"../config_files/summer11mc.input\")"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
cd ../Acceptance
root -b -q -l auxScript.C+
cd ../FullChain
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "DONE: plotDYAcceptance(\"../config_files/summer11mc.input\")"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"

echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "WILL DO: plotDYEfficiency(\"../config_files/summer11mc.input\")"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
cd ../Efficiency
rm *.so
root -b -q -l auxScript.C+
cd ../FullChain
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "DONE: plotDYEfficiency(\"../config_files/summer11mc.input\")"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"

echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "WILL DO: seven scripts from EfficiencyScaleFactors"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
cd ../EfficiencyScaleFactors
rm *.so
root -b -q -l auxScript.C+
cd ../FullChain
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "DONE: seven scripts from EfficiencyScaleFactors"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"

echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "WILL DO: plotDYFSRCorrections(\"../config_files/summer11mc.input\")"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
cd ../Fsr
rm *.so
root -b -q -l auxScript1.C+
cd ../FullChain
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "DONE: plotDYFSRCorrections(\"../config_files/summer11mc.input\")"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"

echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "WILL DO: plotDYFSRCorrectionsSansAcc(\"../config_files/summer11mc.input\")"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
cd ../Fsr
rm *.so
root -b -q -l auxScript2.C+
cd ../FullChain
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "DONE: plotDYFSRCorrectionsSansAcc(\"../config_files/summer11mc.input\")"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"

echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "WILL DO: calcCrossSection(\"../config_files/xsecCalc.conf\")"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
cd ../CrossSection
rm *.so
root -b -q -l auxScript.C+
cd ../FullChain
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "DONE: calcCrossSection(\"../config_files/xsecCalc.conf\")"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
