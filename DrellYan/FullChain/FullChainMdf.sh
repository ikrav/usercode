#!/bin/bash

# redirect stderr to stdout. 
# Note: This redirection makes it impossible to kill the script easily!
#exec 2>&1

# ------------------  Define some variables

filename_data="../config_files/data.conf"
filename_mc="../config_files/summer11mc.input"
filename_cs="../config_files/xsecCalc.conf"
expectEfficiencyScaleFactorsFile="../root_files/constants/DY_m10+pr+a05+o03+pr_4680pb/scale_factors_${triggerSet}.root"
expectUnfoldingSystematicsFile="../root_files/systematics/DY_m10+pr+a05+o03+pr_4680pb/unfolding_systematics.root"

# export some variables

export mainConfInputFile=${filename_data}
export mcConfInputFile=${filename_mc}
export xsecConfInputFile=${filename_cs}

export triggerSet="Full2011_hltEffNew"
export tnpFileStart="../config_files/sf"



# specify whether you want to clean the old logs
clear_old_logs=0

# specify whether the support files need to be rebuilt
force_rebuild_include_files=0

## controlling your work
# catch-all flag
do_all_steps=0

# main procedures: calculate and make plots
do_calcs=0
do_plots=0

# individual flags. 
# Note: all the above flags have to be 0 for these individual flags 
# to be effective
do_selection=1
do_plotSelectDY=0
do_subtractBackground=0
do_unfolding=0
do_unfoldingSyst=0
do_acceptance=0
do_acceptanceSyst=0
do_efficiency=0
do_efficiencyScaleFactors=0
do_plotFSRCorrections=0
do_plotFSRCorrectionsSansAcc=0
do_theoryErrors=0
do_crossSection=0

# use logDir="./" if you want that the log files are placed in the directory
# where the producing script resides
logDir="./"
logDir="../logs"
if [ ! -d ${logDir} ] ; then mkdir ${logDir}; fi
if [ ${clear_old_logs} -eq 1 ] && [ "${logDir}"="../logs" ] ; then
    echo -e "\tclean old logs\n"
    rm -f ${logDir}/*log
fi

# if you do not want to have the time stamp, comment the line away 
# or set timeStamp=
timeStamp="-`date +%Y%m%d-%H%M`"
#timeStamp=

#
# no error flag
#

noError=1


# -------------  Change individual flags in agreement with superior flags

if [ ${do_all_steps} -eq 1 ] || [ ${do_calcs} -eq 1 ] ; then
    do_selection=1
    do_subtractBackground=1
    do_unfolding=1
    do_unfoldingSyst=1
    do_acceptance=1
    do_acceptanceSyst=1
    do_efficiency=1
    do_efficiencyScaleFactors=1
    do_theoryErrors=1
    do_crossSection=1
fi
if [ ${do_all_steps} -eq 1 ] || [ ${do_plots} -eq 1 ] ; then
    do_plotSelectDY=1
    do_plotFSRCorrections=1
    do_plotFSRCorrectionsSansAcc=1
fi


# ------------------------- Define a few functions

checkFile() { 
  __fname="$1"    # double underscore not to mess a possible variable with the same name in the script
# if "-f" does not work (check plain file), 
# one can use -e (name exists), but then directory will return 'true' as well
# 2. else... echo ".. ok" can be removed
  if [ ! -f ${__fname} ] ; then 
      echo "file ${__fname} is missing"
      noError=0
  else
      echo "file ${__fname} checked ok"
  fi
}


# -------------------- Preparatory checks

checkFile ${filename_data}
checkFile ${filename_mc}
checkFile ${filename_cs}


if [ ${force_rebuild_include_files} -eq 1 ] ; then
    echo -e " All libraries in Include will be rebuilt\n"
    cd ../Include
    rm -f *.so
    cd ../FullChain
fi


# -------------------- Main work

# prepare support libraries
cd ../Include
root -b -q -l rootlogon.C+              | tee ${logDir}/out${timeStamp}-00-include.log


#Selection
if [ ${do_selection} -eq 1 ] && [ ${noError} -eq 1 ] ; then
statusSelection=OK
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "WILL DO: selectEvents1D(\"${filename_data}\")"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
cd ../Selection
rm -f *.so
echo
checkFile selectEvents1D.C
root -b -q -l rootlogon.C+ selectEvents1D.C+\(\"$filename_data\"\)           | tee ${logDir}/out${timeStamp}-01-selectEvents1D.log
#echo "exit code {$?}"
if [ $? != 0 ]; then 
   statusSelection=FAILED
   noError=0
fi
cd ../FullChain
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "DONE: selectEvents1D(\"${filename_data}\")"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
else
  statusSelection=skipped
fi

#PlotSelectDY
if [ ${do_plotSelectDY} -eq 1 ] && [ ${noError} -eq 1 ]  ; then
statusPlotSelectDY=OK
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "WILL DO: prepareYields1D.C(\"${filename_data}\")"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
cd ../MassSpectrum
rm -f *.so
echo
checkFile prepareYields1D.C
root -b -q -l rootlogon.C+ prepareYields1D.C+\(\"$filename_data\"\)       | tee ${logDir}/out${timeStamp}-02-prepareYields1D.log
if [ $? != 0 ]; then 
   statusPlotSelectDY=FAILED
   noError=0
fi
cd ../FullChain
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "DONE: prepareYields1D.C(\"${filename_data}\")"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
else 
    statusPlotSelectDY=skipped
fi

#SubtractBackground
if [ ${do_subtractBackground} -eq 1 ] && [ ${noError} -eq 1 ] ; then
statusSubtractBackground=OK
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "WILL DO: subtractBackground.C(\"${filename_data}\")"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
cd ../MassSpectrum
rm -f *.so
echo
checkFile subtractBackground.C
root -b -q -l rootlogon.C+ subtractBackground.C+\(\"$filename_data\"\)      | tee ${logDir}/out${timeStamp}-03-subtractBackground.log
if [ $? != 0 ]; then 
   statusSubtractBackground=FAILED
   noError=0
fi
cd ../FullChain
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "DONE: subtractBackground.C(\"${filename_data}\")"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
else
    statusSubtractBackground=skipped
fi


#Unfolding
if [ ${do_unfolding} -eq 1 ] && [ ${noError} -eq 1 ] ; then
statusUnfolding=OK
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "WILL DO: plotDYUnfoldingMatrix(\"${filename_mc}\")"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
cd ../Unfolding
rm -f *.so
echo
checkFile plotDYUnfoldingMatrix.C
root -b -q -l rootlogon.C+ plotDYUnfoldingMatrix.C+\(\"$filename_mc\"\)    | tee ${logDir}/out${timeStamp}-04-plotDYUnfoldingMatrix.log
if [ $? != 0 ]; then 
   statusUnfolding=FAILED
   noError=0
fi
cd ../FullChain
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "DONE: plotDYUnfoldingMatrix(\"${filename_mc}\")"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
else
    statusUnfolding=skipped
fi

#Unfolding Systematics
if [ ${do_unfoldingSyst} -eq 1 ] && [ ${noError} -eq 1 ] ; then
statusUnfoldingSyst=OK
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "WILL DO: evaluateUnfoldingSyst.sh"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
cd ../Unfolding
rm -f *.so
echo
checkFile evaluateUnfoldingSyst.sh
source evaluateUnfoldingSyst.sh   | tee ${logDir}/out${timeStamp}-05-evaluateUnfoldingSyst.log
if [ $? != 0 ]; then 
   statusUnfoldingSyst=FAILED
   noError=0
else
# check that there was no crash and the expected file was created
if [ ! -f ${expectUnfoldingSystematicsFile} ] ; then
  echo "  ERROR:"
  echo "  !! expectUnfoldingSystematicsFile=${expectUnfoldingSystematicsFile} not created"
  echo
  statusUnfoldingSyst=FAILED
  noError=0
fi
fi 
cd ../FullChain
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "DONE: evaluateUnfoldingSyst.sh"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
else
    statusUnfoldingSyst=skipped
fi

#Acceptance
if [ ${do_acceptance} -eq 1 ] && [ ${noError} -eq 1 ] ; then
statusAcceptance=OK
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "WILL DO: plotDYAcceptance(\"${filename_mc}\")"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
cd ../Acceptance
rm -f *.so
echo
checkFile plotDYAcceptance.C
root -b -q -l rootlogon.C+ plotDYAcceptance.C+\(\"$filename_mc\"\)       | tee ${logDir}/out${timeStamp}-06-plotDYAcceptance.log
if [ $? != 0 ]; then 
   statusAcceptance=FAILED
   noError=0
fi
cd ../FullChain
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "DONE: plotDYAcceptance(\"${filename_mc}\")"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
else
    statusAcceptance=skipped
fi

#Acceptance Systematics
if [ ${do_acceptanceSyst} -eq 1 ] && [ ${noError} -eq 1 ] ; then
statusAcceptanceSyst=OK
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "WILL DO: plotDYAcceptance(\"${filename_mc}\")"
echo "in systematics mode"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
cd ../Acceptance
rm -f *.so
echo
checkFile auxScriptSyst.C
root -b -q -l rootlogon.C+ auxScriptSyst.C       | tee ${logDir}/out${timeStamp}-07-auxScriptSyst-acceptance.log
if [ $? != 0 ]; then 
#   statusUnfoldingSyst=FAILED         ---- BUG!!!
    statusAcceptanceSyst=FAILED
    noError=0
fi 
cd ../FullChain
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "DONE: plotDYAcceptance(\"${filename_mc}\")"
echo "in systematics mode"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
else
    statusAcceptanceSyst=skipped
fi

#Efficiency
if [ ${do_efficiency} -eq 1 ] && [ ${noError} -eq 1 ] ; then
statusEfficiency=OK
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "WILL DO: plotDYEfficiency(\"${filename_mc}\")"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
cd ../Efficiency
rm -f *.so
echo
checkFile plotDYEfficiency.C
root -b -q -l rootlogon.C+ plotDYEfficiency.C+\(\"$filename_mc\"\)       | tee ${logDir}/out${timeStamp}-08-plotDYEfficiency.log
if [ $? != 0 ]; then 
   statusEfficiency=FAILED
   noError=0
fi
cd ../FullChain
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "DONE: plotDYEfficiency(\"${filename_mc}\")"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
else
    statusEfficiency=skipped
fi

#Efficiency Scale Factors
if [ ${do_efficiencyScaleFactors} -eq 1 ] && [ ${noError} -eq 1 ] ; then
statusEfficiencyScaleFactors=OK
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "WILL DO: evaluateESF.sh in EfficiencyScaleFactors"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
cd ../EfficiencyScaleFactors
rm -f *.so
echo
checkFile evaluateESF.sh
#source evaluateESF.sh  | tee ${logDir}/out${timeStamp}-09-evaluateESF-efficiencyScaleFactors.log
if [ $? != 0 ]; then 
   statusEfficiencyScaleFactors=FAILED
   noError=0
else 
# check that there was no crash and the expected file was created
if [ ! -f ${expectEfficiencyScaleFactorsFile} ] ; then
  echo "  ERROR:"
  echo "  !! expectEfficiencyScaleFactorsFile=${expectEfficiencyScaleFactorsFile} not created"
  echo
  statusEfficiencyScaleFactors=FAILED
  noError=0
fi
fi
cd ../FullChain
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "DONE: evaluateESF.sh from EfficiencyScaleFactors"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
else 
    statusEfficiencyScaleFactors=skipped
fi


#PlotDYFSRCorrections
if [ ${do_plotFSRCorrections} -eq 1 ] && [ ${noError} -eq 1 ] ; then
statusPlotDYFSRCorrections=OK
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "WILL DO: plotDYFSRCorrections(\"${filename_mc}\")"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
cd ../Fsr
rm -f *.so
echo
checkFile plotDYFSRCorrections.C
root -b -q -l rootlogon.C+ plotDYFSRCorrections.C+\(\"$filename_mc\"\)     | tee ${logDir}/out${timeStamp}-10-plotDYFSRCorrections.log
if [ $? != 0 ]; then 
   statusPlotDYFSRCorrections=FAILED
   noError=0
fi
cd ../FullChain
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "DONE: plotDYFSRCorrections(\"${filename_mc}\")"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
else
    statusPlotDYFSRCorrections=skipped
fi

#PlotDYFSRCorrectionsSansAcc
if [ ${do_plotFSRCorrectionsSansAcc} -eq 1 ] && [ ${noError} -eq 1 ] ; then
statusPlotDYFSRCorrectionsSansAcc=OK
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "WILL DO: plotDYFSRCorrectionsSansAcc(\"${filename_mc}\")"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
cd ../Fsr
rm -f *.so
echo
checkFile plotDYFSRCorrectionsSansAcc.C
root -b -q -l rootlogon.C+ plotDYFSRCorrectionsSansAcc.C+\(\"$filename_mc\"\)     | tee ${logDir}/out${timeStamp}-11-plotDYFSRCorrectionsSansAcc${timeStamp}.out
if [ $? != 0 ]; then 
   statusPlotDYFSRCorrectionsSansAcc=FAILED
   noError=0
fi
cd ../FullChain
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "DONE: plotDYFSRCorrectionsSansAcc(\"${filename_mc}\")"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
else
    statusPlotDYFSRCorrectionsSansAcc=skipped
fi

#TheoryErrors
if [ ${do_theoryErrors} -eq 1 ] && [ ${noError} -eq 1 ] ; then
statusTheoryErrors=OK
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "WILL DO: TheoryErrors(\"${filename_mc}\")"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
cd ../Theory
rm -f *.so
echo
checkFile TheoryErrors.C
root -b -q -l rootlogon.C+ TheoryErrors.C+\(\"$filename_mc\"\)     | tee ${logDir}/out${timeStamp}-12-TheoryErrors${timeStamp}.out
if [ $? != 0 ]; then 
   statusTheoryErrors=FAILED
   noError=0
fi
cd ../FullChain
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "DONE: TheoryErrors(\"${filename_mc}\")"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
else
    statusTheoryErrors=skipped
fi

#CrossSection
if [ ${do_crossSection} -eq 1 ] && [ ${noError} -eq 1 ] ; then
statusCrossSection=OK
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "WILL DO: CrossSection(\"${filename_cs}\")"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
cd ../CrossSection
rm -f *.so
echo
checkFile calcCrossSection.C
root -b -q -l rootlogon.C+ calcCrossSection.C+\(\"$filename_cs\"\)     | tee ${logDir}/out${timeStamp}-13-CrossSection${timeStamp}.out
if [ $? != 0 ]; then 
   statusCrossSection=FAILED
   noError=0
fi
cd ../FullChain
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
echo "DONE: CrossSection(\"${filename_cs}\")"
echo "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD"
else
    statusCrossSection=skipped
fi

# ------------------------------ final summary

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
if [ ${noError} -eq 0 ] ; then 
  echo
  echo " !! ERROR was detected !!"
  echo
fi

