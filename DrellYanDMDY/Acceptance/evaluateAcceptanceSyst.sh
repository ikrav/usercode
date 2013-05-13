#!/bin/bash

mcConfInputFile="../config_files/summer12mc.input" 
debugMode=1

if [ ${#1} -gt 0 ] ; then mcConfInputFile=$1; fi
if [ ${#2} -gt 0 ] ; then debugMode=$2; fi

#
# Check if the environment variables are set. Assign values if they are empty
#

#if [ -s ${triggerSet} ] ; then  
#    triggerSet="Full2011_hltEffNew"  # not used in this script yet
#fi
if [ -s {mcConfInputFile} ] || [ ${#mcConfInputFile} -eq 0 ] ; then
    mcConfInputFile="../config_files/summer12mc.input" 
fi
if [ -s ${xsecConfInputFile} ] || [ ${#xsecConfInputFile} -eq 0 ] ; then
    xsecConfInputFile="../config_files/xsecCalc8TeV.conf"
fi
if [ -s ${debugMode} ] || [ ${#debugMode} -eq 0 ] ; then
    debugMode=0
fi

# check whether the full run was requested, overriding internal settings
if [ -s ${fullRun} ] ; then
  fullRun=0
fi

echo
echo
echo "evaluateAcceptanceSyst.sh:"
#echo "    triggerSet=${triggerSet} (not used)"
echo "    mcConfInputFile=${mcConfInputFile}"
echo "    xsecConfInputFile=${xsecConfInputFile}"
echo 
echo

# 
#  Individual flags to control the calculation
#

doFsrStudy=1
doCalcAcceptanceSyst=1

#
#  Modify flags if fullRun=1
#

if [ ${fullRun} -eq 1 ] ; then
    doFsrStudy=1; doCalcAcceptanceSyst=1
fi

#
#  Flag of an error
#
noError=1

# --------------------------------
#    Define functions to run
# --------------------------------

runPlotDYAcceptance() {
  loc_massLimit=-1
  root -b -q -l ${LXPLUS_CORRECTION} plotDYAcceptance.C+\(\"${mcConfInputFile}\",${StudyFlag},${ReweightFactor},${loc_massLimit},${debugMode}\)
  if [ $? != 0 ] ; then noError=0;
  else 
     echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
     echo 
     echo "DONE: plotDYAcceptance.C(\"${mcConfInputFile}\",${StudyFlag},${ReweightFactor},${loc_massLimit},debug=${debugMode})"
     echo 
     echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
  fi
}

#void plotDYAcceptance(const TString input, int systematicsMode = DYTools::NORMAL, double reweightFsr = 1.0, double massLimit=-1)
#//systematicsMode 0 (NORMAL) - no systematic calc
#//2 (FSR_STUDY) - systematics due to FSR, reweighting
#//check mass spectra with reweight = 95%; 100%; 105%  
#//mass value until which do reweighting 

runCalcAcceptanceSystematics() {
  root -b -q -l ${LXPLUS_CORRECTION} calcAcceptanceSystematics.C+\(\"${xsecConfInputFile}\"\)
  if [ $? != 0 ] ; then noError=0;
  else 
     echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
     echo 
     echo "DONE: calcAcceptanceSystematics.C+(\"${xsecConfInputFile}\")"
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
StudyFlag=0; ReweightFactor=0
if [ ${noError} -eq 1 ] ; then runPlotDYAcceptance; fi
mcConfInputFile=${storeMCConfInputFile}

storeXSecConfInputFile=${xsecConfInputFile}
xsecConfInputFile="_DebugRun_"
if [ ${noError} -eq 1 ] ; then runCalcAcceptanceSystematics; fi
xsecConfInputFile=${storeXSecConfInputFile}


#
#   Calculations
#

if [ ${doFsrStudy} -eq 1 ] && [ ${noError} -eq 1 ] ; then
  StudyFlag="DYTools::FSR_STUDY"
  loopReweightFactor="1.05 0.95"
  for ReweightFactor in ${loopReweightFactor} ; do
    if [ ${noError} -eq 1 ] ; then
      runPlotDYAcceptance
    fi
  done
fi


if [ ${doCalcAcceptanceSyst} -eq 1 ] && [ ${noError} -eq 1 ] ; then
  runCalcAcceptanceSystematics
fi


# return the error code
exit ${noError}
