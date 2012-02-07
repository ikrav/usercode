#!/bin/bash

#
# Check if the environment variables are set. Assign values if they are empty
#
if [ -s ${triggerSet} ] ; then  
    triggerSet="Full2011_hltEffNew"
fi
if [ -s ${mcConfInputFile} ] ; then
    mcConfInputFile="../config_files/fall11mc.input" # used in CalcEventEff.C
fi
if [ -s ${tnpFileStart} ] ; then
  tnpFileStart="../config_files/sf"   # beginning of the file used in eff_*.C 
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
echo "    tnpFileStart=${tnpFileStart}"
echo 
echo

# 
#  Individual flags to control the calculation
#

runMC_Reco=1
runMC_Id=1
runMC_Hlt=1
runData_Reco=1
runData_Id=1
runData_Hlt=1
runCalcEventEff=1

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
tmp1=${triggerSet/hltEffNew/}  # replace hltEffNew with nothing
tmp2=${triggerSet/Full2011/}   # replace Full2011 with nothing
#echo "lengths=${#triggerSet}, ${#tmp1}, ${#tmp2}"
# compare the lengths
if [ ${#triggerSet} -ne ${#tmp1} ] && 
   [ ${#triggerSet} -ne ${#tmp2} ] ; then
  lumiWeighting=1
else
  lumiWeighting=0
fi


# --------------------------------
#    Define functions to run
# --------------------------------

runEffReco() {
 root -b -q -l rootlogon.C+ eff_Reco.C+\(\"${inpFile}\",\"${triggerSet}\"\)
  if [ $? != 0 ] ; then noError=0;
  else 
     echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
     echo 
     echo "DONE: eff_Reco(\"$inpFile\",\"${triggerSet}\")"
     echo 
     echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
  fi
}


runEffIdHlt() {
 root -b -q -l rootlogon.C+ eff_IdHlt.C+\(\"${inpFile}\",\"${triggerSet}\"\)
  if [ $? != 0 ] ; then noError=0;
  else 
     echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
     echo 
     echo "DONE: eff_IdHlt(\"$inpFile\",\"${triggerSet}\")"
     echo 
     echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
  fi
}

runCalcEventEff() {
 root -b -q -l rootlogon.C+ calcEventEff.C+\(\"${mcConfInputFile}\",\"${triggerSet}\"\)
  if [ $? != 0 ] ; then noError=0;
  else 
     echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
     echo 
     echo "DONE: eff_calcEventEff(\"$mcConfInputFile\",\"${triggerSet}\")"
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

storeTriggerSet=${triggerSet}
triggerSet="_DebugRun_"
if [ ${noError} -eq 1 ] ; then runEffReco; fi
if [ ${noError} -eq 1 ] ; then runEffIdHlt; fi
if [ ${noError} -eq 1 ] ; then runCalcEventEff; fi
if [ ${noError} -eq 1 ] ; then echo; echo "  -=- Resuming normal calculation -=-"; echo; fi
triggerSet=${storeTriggerSet}


# Process MC

if [ ${runMC_Reco} -eq 1 ] && [ ${noError} -eq 1 ] ; then
  inpFile="${tnpFileStart}_mc_RECO.conf"
  runEffReco
fi

if [ ${runMC_Id} -eq 1 ] && [ ${noError} -eq 1 ] ; then
  inpFile="${tnpFileStart}_mc_ID.conf"
  runEffIdHlt
fi

if [ ${runMC_Hlt} -eq 1 ] && [ ${noError} -eq 1 ] ; then
  inpFile="${tnpFileStart}_mc_HLT.conf"
  runEffIdHlt
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
    inpFile="${tnpFileStart}_data_RECO.conf"
    runEffReco
  fi
done

for triggerSet in ${loopTriggers} ; do
  if [ ${runData_Id} -eq 1 ] && [ ${noError} -eq 1 ] ; then
    inpFile="${tnpFileStart}_data_ID.conf"
    runEffIdHlt
  fi
done

for triggerSet in ${loopTriggers} ; do
  if [ ${runData_Hlt} -eq 1 ] && [ ${noError} -eq 1 ] ; then
    inpFile="${tnpFileStart}_data_HLT.conf"
    runEffIdHlt
  fi
done

triggerSet=${storeTriggerSet}


# Calculate efficiency scale factors

if [ ${runCalcEventEff} -eq 1 ] && [ ${noError} -eq 1 ] ; then
    runCalcEventEff
fi


# return the error code
exit ${noError}

