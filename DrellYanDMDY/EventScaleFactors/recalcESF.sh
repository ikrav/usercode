#!/bin/bash

#
# use for debug, tuning and recalculation of the efficiencies or scale factors
# evaluateESF.sh is the main script
#

debugMode=0
fullRun=1
puReweight=1

if [ ${#1} -gt 0 ] ; then mcConfInputFile=$1; fi
if [ ${#2} -gt 0 ] ; then triggerSet=$2; fi
if [ ${#3} -gt 0 ] ; then debugMode=$3; fi

tnpMCFile="../config_files/sf_mc_eta2.conf"
#tnpMCFile="../config_files/sf_mc_evtTrig.conf"
#tnpMCFile="../config_files/sf_mc_evtTrig_eta2.conf"
#tnpMCFile="../config_files/sf_mc_spektras_evtTrig_eta2.conf"
#tnpMCFile="../config_files/sf_mc_evtTrig_eta2_test11.conf"

tnpDataFile="../config_files/sf_data_eta2.conf"
#tnpDataFile="../config_files/sf_data_evtTrig.conf"
#tnpDataFile="../config_files/sf_data_evtTrig_eta2.conf"
#tnpDataFile="../config_files/sf_data_spektras_evtTrig_eta2.conf"
#tnpDataFile="../config_files/sf_data_evtTrig_eta2_test11.conf"

collectEvents=0  # it is recommended to have collectEvents=1 in evaluateESF!

# if you do not want to have the time stamp, comment the line away 
# or set timeStamp=
timeStamp="-`date +%Y%m%d-%H%M`"
#timeStamp=

#
# Check if the environment variables are set. Assign values if they are empty
#
if [ ${#triggerSet} -eq 0 ] ; then  
    triggerSet="Full2011_hltEffOld"
fi
if [ ${#mcConfInputFile} -eq 0 ] ; then
    mcConfInputFile="../config_files/fall11mc.input" # used in CalcEventEff.C
fi
if [ ${#tnpMCFile} -eq 0 ] ; then
  tnpMCFile="../config_files/sf_mc.conf"   # file for eff_*.C
fi
if [ ${#tnpDataFile} -eq 0 ] ; then
  tnpDataFile="../config_files/sf_data.conf"   # file for eff_*.C
fi
if [ ${#debugMode} -eq 0 ] ; then
    debugMode=1
fi

# check whether the full run was requested, overriding internal settings
if [ -s ${fullRun} ] ; then
  fullRun=0
fi

echo
echo
echo "calcEffScaleFactors.sh:"
echo "    triggerSet=${triggerSet}"
echo "    mcConfInputFile=${mcConfInputFile}"
echo "    tnpMCFile=${tnpMCFile}"
echo "    tnpDataFile=${tnpDataFile}"
echo "    timeStamp=${timeStamp}"
echo "    debugMode=${debugMode}"
echo 
echo

# 
#  Individual flags to control the calculation
#

puDependence=0  # as a function of nPV's


runMC_Reco=1
runMC_Id=0
runMC_Hlt=0
runData_Reco=0
runData_Id=0
runData_Hlt=0
runCalcEventEff=0

#
#  Modify flags if fullRun=1
#

if [ ${fullRun} -eq 1 ] ; then
  runMC_Reco=1; runMC_Id=1; runMC_Hlt=1
  runData_Reco=1; runData_Id=1; runData_Hlt=1
  runCalcEventEff=1
fi


#
#  Flag of an error
#
noError=1


# determine whether a triple run on data is required
lumiWeighting=0
# Lumi weighting is disabled from April 01, 2012
#tmp1=${triggerSet/hltEffNew/}  # replace hltEffNew with nothing
#tmp2=${triggerSet/Full2011/}   # replace Full2011 with nothing
##echo "lengths=${#triggerSet}, ${#tmp1}, ${#tmp2}"
## compare the lengths
#if [ ${#triggerSet} -ne ${#tmp1} ] && 
#   [ ${#triggerSet} -ne ${#tmp2} ] ; then
#  lumiWeighting=1
#else
#  lumiWeighting=0
#fi


# --------------------------------
#    Define functions to run
# --------------------------------

checkFile() { 
  if [ ${noError} -eq 0 ] ; then 
      noError=0
      echo "error present before checking for the file(s) $@"
      return
  fi

  # double underscore not to mess a possible variable with the same name in the script
  for __fname in $@ ; do
# if "-f" does not work (check plain file), 
# one can use -e (name exists), but then directory will return 'true' as well
# 2. else... echo ".. ok" can be removed
  if [ ${#__fname} -gt 0 ] ; then 
      if [ ! -f ${__fname} ] ; then 
	  echo "file ${__fname} is missing"
	  noError=0
      else
	  echo "file <${__fname}> checked ok"
      fi
  fi
  done
}

runCalcEff() {
 effKind=$1
 dataKind=${inpFile/data/}
 if [ ${#dataKind} -eq ${#inpFile} ] ; then dataKind="mc"; else dataKind="data"; fi
# calculate
 root -l -q -b  ${LXPLUS_CORRECTION} calcEff.C+\(\"${inpFile}\",\"${effKind}\",\"${triggerSet}\",${puReweight},${puDependence}\) \
     | tee log${timeStamp}-calcEff-${dataKind}-${effKind}-puW${puReweight}.out
  if [ $? != 0 ] ; then noError=0;
  else 
     checkFile calcEff_C.so
     echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
     echo 
     echo "DONE: calcEff(\"$inpFile\",\"${effKind}\",\"${triggerSet}\",puReweight=${puReweight},puDep=${puDependence})" #,debug=${debugMode})"
     echo 
     echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
  fi
}

runCalcEventEff() {
 _collectEvents=$1
 echo "_collectEvents=${_collectEvents}"
 if [ ${#_collectEvents} -eq 0 ] ; then _collectEvents=1; fi
 root -l -b -q ${LXPLUS_CORRECTION} calcEventEff.C+\(\"${mcConfInputFile}\",\"${tnpDataFile}\",\"${tnpMCFile}\",\"${triggerSet}\",${_collectEvents},${puReweight},${debugMode}\) \
     | tee log${timeStamp}-calcEventEff-puW${puReweight}.out
  if [ $? != 0 ] ; then noError=0;
  else 
      checkFile calcEventEff_C.so
     echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
     echo 
     echo "DONE: calcEventEff(\"${mcConfInputFile}\",\"${tnpDataFile}\",\"${tnpMCFile}\",\"${triggerSet}\",collectEvents=${_collectEvents},puReweight=${puReweight},debug=${debugMode})"
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
checkFile tnpSelectEvents_hh.so

# 
#   Check that the codes compile
#

storeTriggerSet=${triggerSet}
triggerSet="_DebugRun_"
runAny=$(( ${runMC_Reco} + ${runMC_Id} + ${runMC_Hlt} + ${runData_Reco} + ${runData_Id} + ${runData_Hlt} ))
if [ ${runAny} -gt 0 ] && [ ${noError} -eq 1 ] ; then runCalcEff; fi
if [ ${runCalcEventEff} -eq 1 ] && [ ${noError} -eq 1 ] ; then runCalcEventEff; fi
if [ ${noError} -eq 1 ] ; then echo; echo "  -=- Resuming normal calculation -=-"; echo; fi
triggerSet=${storeTriggerSet}


# Process MC

if [ ${runMC_Reco} -eq 1 ] && [ ${noError} -eq 1 ] ; then
  inpFile="${tnpMCFile}"
  runCalcEff "RECO"
fi

if [ ${runMC_Id} -eq 1 ] && [ ${noError} -eq 1 ] ; then
  inpFile="${tnpMCFile}"
  runCalcEff "ID"
fi

if [ ${runMC_Hlt} -eq 1 ] && [ ${noError} -eq 1 ] ; then
  inpFile="${tnpMCFile}"
  runCalcEff "HLT"
fi


# Process data

storeTriggerSet=${triggerSet}

if [ ${lumiWeighting} -eq 0 ] ; then
  loopTriggers="${triggerSet}"
else
  loopTriggers="2011A_SingleEG_hltEffNew 2011A_DoubleEG_hltEffNew 2011B_DoubleEG_hltEffNew"
fi

for triggerSet in ${loopTriggers} ; do
  if [ ${runData_Reco} -eq 1 ] && [ ${noError} -eq 1 ] ; then
    inpFile="${tnpDataFile}"
    runCalcEff "RECO"
  fi
done

for triggerSet in ${loopTriggers} ; do
  if [ ${runData_Id} -eq 1 ] && [ ${noError} -eq 1 ] ; then
    inpFile="${tnpDataFile}"
    runCalcEff "ID"
  fi
done

for triggerSet in ${loopTriggers} ; do
  if [ ${runData_Hlt} -eq 1 ] && [ ${noError} -eq 1 ] ; then
    inpFile="${tnpDataFile}"
    runCalcEff "HLT"
  fi
done

triggerSet=${storeTriggerSet}


# Calculate efficiency scale factors

if [ ${runCalcEventEff} -eq 1 ] && [ ${noError} -eq 1 ] ; then
    runCalcEventEff ${collectEvents}
fi


# return the error code
exit ${noError}

