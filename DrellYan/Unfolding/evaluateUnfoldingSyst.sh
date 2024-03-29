#!/bin/bash

#
# Check if the environment variables are set. Assign values if they are empty
#

if [ -s ${triggerSet} ] ; then  
    triggerSet="Full2011_hltEffNew"  # not used in this script yet
fi
if [ -s ${mcConfInputFile} ] ; then
    mcConfInputFile="../config_files/fall11mc.input" 
fi
if [ -s ${xsecConfInputFile} ] ; then
    xsecConfInputFile="../config_files/xsecCalc.conf"
fi

# check whether the full run was requested, overriding internal settings
if [ -s ${fullRun} ] ; then
  fullRun=0
fi

echo
echo
echo "evaluateUnfoldingSyst.sh:"
echo "    triggerSet=${triggerSet} (not used)"
echo "    mcConfInputFile=${mcConfInputFile}"
echo "    xsecConfInputFile=${xsecConfInputFile}"
echo 
echo

# 
#  Individual flags to control the calculation
#

doFsrStudy=1
doResolutionStudy=1
doCalcUnfoldingSyst=1

#
#  Modify flags if fullRun=1
#

if [ ${fullRun} -eq 1 ] ; then
    doFsrStudy=1; doResolutionStudy=1; doCalcUnfoldingSyst=1
fi

#
#  Flag of an error
#
noError=1

# --------------------------------
#    Define functions to run
# --------------------------------

runPlotDYUnfoldingMatrix() {
  loc_massLimit=-1
  root -b -q -l ${LXPLUS_CORRECTION} plotDYUnfoldingMatrix.C+\(\"${mcConfInputFile}\",${StudyFlag},${RandomSeed},${ReweightFsr},${loc_massLimit}\)
  if [ $? != 0 ] ; then noError=0;
  else 
     echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
     echo 
     echo "DONE: plotDYUnfoldingMatrix.C(\"${mcConfInputFile}\",${StudyFlag},${RandomSeed},${ReweightFsr},${loc_massLimit})"
     echo 
     echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
  fi
}

#info:
# // void plotDYUnfoldingMatrix(const TString input, int systematicsMode = 0, int randomSeed = 1, double reweightFsr = 1.00, double massLimit = -1)
# // systematicsMode 0 - no systematic calc, no reweighting
# // 1 - systematic mode, 2 - (reweighting of mass diff < -1 GeV) mode
#  //check mass spectra with reweight = 95%; 100%; 105%  
#  //mass value until which do reweighting


runCalcUnfoldingSystematics() {
  root -b -q -l ${LXPLUS_CORRECTION} calcUnfoldingSystematics.C+\(\"${xsecConfInputFile}\"\)
  if [ $? != 0 ] ; then noError=0;
  else 
     echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
     echo 
     echo "DONE: calcUnfoldingSystematics.C+(\"${xsecConfInputFile}\")"
     echo 
     echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
  fi
}


# --------------------------------
#    Main sequence
# --------------------------------

#
#  Compile header files
#
root -b -q -l rootlogon.C+
if [ $? != 0 ] ; then noError=0; fi

# 
#   Check that the codes compile
#

storeMCConfInputFile=${mcConfInputFile}
mcConfInputFile="_DebugRun_"
StudyFlag=0; RandomSeed=0; ReweightFsr=0
if [ ${noError} -eq 1 ] ; then runPlotDYUnfoldingMatrix; fi
mcConfInputFile=${storeMCConfInputFile}

storeXSecConfInputFile=${xsecConfInputFile}
xsecConfInputFile="_DebugRun_"
if [ ${noError} -eq 1 ] ; then runCalcUnfoldingSystematics; fi
xsecConfInputFile=${storeXSecConfInputFile}


#
#   Calculations
#

if [ ${doFsrStudy} -eq 1 ] && [ ${noError} -eq 1 ] ; then
  StudyFlag="DYTools::FSR_STUDY"; RandomSeed=1
  loopReweightFsr="1.05 0.95"
  for ReweightFsr in ${loopReweightFsr} ; do
    if [ ${noError} -eq 1 ] ; then
      runPlotDYUnfoldingMatrix
    fi
  done
fi


if [ ${doResolutionStudy} -eq 1 ] && [ ${noError} -eq 1 ] ; then
  StudyFlag="DYTools::RESOLUTION_STUDY"
  ReweightFsr="1.0"
  i=0
  while [ ${i} -le 20 ] ; do
      RandomSeed=$(( 1000 + $i ))
      runPlotDYUnfoldingMatrix
      i=$(( $i + 1 ))
  done
fi

if [ ${doCalcUnfoldingSyst} -eq 1 ] && [ ${noError} -eq 1 ] ; then
  runCalcUnfoldingSystematics
fi


# return the error code
exit ${noError}
