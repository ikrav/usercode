#!/bin/bash

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

tnpDataFile="../config_files/sf_data_eta2.conf"
#tnpDataFile="../config_files/sf_data_evtTrig.conf"
#tnpDataFile="../config_files/sf_data_evtTrig_eta2.conf"
#tnpDataFile="../config_files/sf_data_spektras_evtTrig_eta2.conf"


collectEvents=1 # recommended to have it set to 1. calcEventEff prepares skim fil

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

runMC_Reco=0
runMC_Id=0
runMC_Hlt=0
runData_Reco=0
runData_Id=1
runData_Hlt=1
runCalcEventEff=0

#
#  Modify flags if fullRun=1
#

if [ ${fullRun} -eq 1 ] ; then
  runMC_Reco=1; runMC_Id=1; runMC_Hlt=1
  runData_Reco=1; runData_Id=1; runData_Hlt=1
  runCalcEventEff=1   # it prepares the skim for event efficiencies
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


runEffReco() {
 dataKind=${inpFile/data/}
 if [ ${#dataKind} -eq ${#inpFile} ] ; then dataKind="mc"; else dataKind="data"; fi
# calculate
 root -b -q -l  eff_Reco.C+\(\"${inpFile}\",\"RECO\",\"${triggerSet}\",${puReweight},${debugMode}\) \
     | tee log${timeStamp}-${dataKind}-RECO-puW${puReweight}.out
  if [ $? != 0 ] ; then noError=0;
  else
     checkFile eff_Reco_C.so
     echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
     echo 
     echo "DONE: eff_Reco(\"$inpFile\",\"RECO\",\"${triggerSet}\",puReweight=${puReweight},debug=${debugMode})"
     echo 
     echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
  fi
}


runEffIdHlt() {
 effKind=$1
 if [ ${#dataKind} -eq ${#inpFile} ] ; then dataKind="mc"; else dataKind="data"; fi
# calculate
 root -b -q -l  eff_IdHlt.C+\(\"${inpFile}\",\"${effKind}\",\"${triggerSet}\",${puReweight},${debugMode}\) \
     | tee log${timeStamp}-${dataKind}-${effKind}-puW${puReweight}.out
  if [ $? != 0 ] ; then noError=0;
  else 
     checkFile eff_IdHlt_C.so
     echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
     echo 
     echo "DONE: eff_IdHlt(\"$inpFile\",\"${effKind}\",\"${triggerSet}\",puReweight=${puReweight},debug=${debugMode})"
     echo 
     echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
  fi
}

runCalcEventEff() {
 _collectEvents=$1
 if [ ${#_collectEvents} -eq 0 ] ; then _collectEvents=1; fi
 root -b -q -l  calcEventEff.C+\(\"${mcConfInputFile}\",\"${tnpDataFile}\",\"${tnpMCFile}\",\"${triggerSet}\",${_collectEvents},${puReweight},${debugMode}\) \
     | tee log${timeStamp}-calcEventEff-puW${puReweight}.out
  if [ $? != 0 ] ; then noError=0;
  else 
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
if [ $(( ${runMC_Reco} + ${runData_Reco} )) -gt 0 ] && [ ${noError} -eq 1 ] ; then runEffReco; fi
doIdHlt=$(( ${runMC_Id} + ${runMC_Hlt} + ${runData_Id} + ${runData_Hlt} ))
if [ ${doIdHlt} -gt 0 ] && [ ${noError} -eq 1 ] ; then runEffIdHlt "ID"; fi
if [ ${runCalcEventEff} -eq 1 ] && [ ${noError} -eq 1 ] ; then runCalcEventEff; fi
if [ ${noError} -eq 1 ] ; then echo; echo "  -=- Resuming normal calculation -=-"; echo; fi
triggerSet=${storeTriggerSet}



# Process MC

if [ ${runMC_Reco} -eq 1 ] && [ ${noError} -eq 1 ] ; then
  inpFile="${tnpMCFile}"
  runEffReco
fi

if [ ${runMC_Id} -eq 1 ] && [ ${noError} -eq 1 ] ; then
  inpFile="${tnpMCFile}"
  runEffIdHlt "ID"
fi

if [ ${runMC_Hlt} -eq 1 ] && [ ${noError} -eq 1 ] ; then
  inpFile="${tnpMCFile}"
  runEffIdHlt "HLT"
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
    runEffReco
  fi
done

for triggerSet in ${loopTriggers} ; do
  if [ ${runData_Id} -eq 1 ] && [ ${noError} -eq 1 ] ; then
    inpFile="${tnpDataFile}"
    runEffIdHlt "ID"
  fi
done

for triggerSet in ${loopTriggers} ; do
  if [ ${runData_Hlt} -eq 1 ] && [ ${noError} -eq 1 ] ; then
    inpFile="${tnpDataFile}"
    runEffIdHlt "HLT"
  fi
done

triggerSet=${storeTriggerSet}


# Calculate efficiency scale factors

if [ ${runCalcEventEff} -eq 1 ] && [ ${noError} -eq 1 ] ; then
    runCalcEventEff ${collectEvents}
fi


# return the error code
exit ${noError}

