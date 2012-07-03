#!/bin/bash 

#
# some variables
#

debugMode=1
fullRun=0
study2D=1

# 1) user-defined

reselectEvents=1   # study2D does not require reselection, in general
#extraPath="../"
extraPath="./"
dirTag="DY_m10+pr+a05+o03+pr_4680pb";
workConfFile="../config_files/data_emu.conf"    # needed to select events

anTagUser="" #"ymax9"
#triggerSet="Full2011_hltEffOld"

copySelectedNTuples=1   # the main sequence has to be run before this script
   # but the ntuples needs to be copied only once

# two individual parts 
doTrue2eBkg=1    # part 1
doFake2eBkg=0    # part 2

# 2) script-internal

plotsDirExtraTag=""

if [ ${study2D} -eq 1 ] ; then
  anTag="2D${anTagUser}"
else
  anTag="1D${anTagUser}"
fi

selectionRunMode=DYTOOLS::NORMAL

eeSelEventsDirT="../root_files/selected_events/TMPDIR/ntuples"
emuSelEventsDirT="../root_files/selected_events/TMPDIR/ntuples_emu"
resultDirT="../root_files/resultants/TMPDIR/"

eeNtuplesDirMain=${eeSelEventsDirT/TMPDIR/${dirTag}}
emuNtuplesDirMain=${emuSelEventsDirT/TMPDIR/${dirTag}}
resultDirMain=${resultDirT/TMPDIR/${dirTag}}

runPath=${PWD}
logPath=${PWD}/dir-dataDrivenBkg-logs
if [ ! -e ${logPath} ] ; then  mkdir ${logPath}; fi

if [ ${fullRun} -eq 1 ] ; then
  doTrue2eBkg=1
  doFake2eBkg=1
fi

err=0


#
#  some functions
#

# -------------------

checkDirs() {
  if [ ${err} -eq 1 ] ; then return; fi
  if [ ${#ntuplesDir} -eq 0 ] || [ ${#yieldsDir} -eq 0 ] || [ ${#resultDir} -eq 0 ] ; then
    echo " one of the directories' variables is empty:"
    echo "   ntuplesDir=${ntuplesDir}"
    echo "   resultDir=${resultDir}"
    exit 1
  fi
}


# -------------------- main cycle
#
# the calculation is performed in temporary directory, then the files
# may be moved to the destination directory
# In some cases the calculation is immediately performed in destination dirs
#


cleanFiles() {
#
#  clean up the files to make sure we are producing them now
#
echo -e "\n\tcleanFiles inactive\n"
return
  rm -f ${resultDirMain}/
  rm -f ${yieldsDirTmp}yields${anTag}.root
  rm -f ${yieldsDirTmp}yields_bg-subtracted${anTag}.root
  rm -f ${resultDirTmp}unfolding_resultants${anTag}.root
  #rm -f ${ntuplesDir}*root 
  rm -f ${yieldsDir}yields${anTag}.root
  rm -f ${yieldsDir}yields_bg-subtracted${anTag}.root
  rm -f ${resultDir}unfolding_resultants${anTag}.root
}


# -----------------  selectEvents ----------------
evaluateTrue2eBkg() {

  cleanFiles
#
#  first try to compile the needed code
#
  cd ../DataDrivenBackgrounds
  if [ ${err} -eq 0 ] ; then
      cd eMuMethod
      rm -f eMuBkgExe
      gmake eMuBkgExe 2>&1 | tee ${logPath}/log-gMake-eMu.out
      testFileExists eMuBkgExe
      cd ..
  fi
#
#   select events
# 
  if [ ${err} -eq 0 ] && [ ${reselectEvents} -eq 1 ] ; then 
    root -l -b -q selectEmuEvents.C+\(\"${workConfFile}\",${debugMode}\) \
        | tee ${logPath}/log-selectEmuEvents.out
    testFileExists ${emuNtuplesDirMain}data${anTagUser}_select.root \
	${emuNtuplesDirMain}zee${anTagUser}_select.root
  fi
#
#   evaluate true2e backgrounds
#
  if [ ${err} -eq 0 ] ; then
      cd eMuMethod
      flags=" --saveRootFile"
      if [ ${study2D} -eq 1 ] ; then
	  flags="--doDMDY ${flags}"
      fi
      ./eMuBkgExe ${flags} 2>&1 | tee ${logPath}/log-emu.out
      testFileExists ${resultDir}/true2eBkgDataPoints_tmp.root
      cd ..
  fi
  
  cd ${runPath}
}


# -------------------

testFileExists() {
  if [ $? != 0 ] ; then err=1; fi
  file=$1
  echo "test file=${file}"
  if [ ! -f ${file} ] ; then echo "${file} is missing"; err=1; fi
  if [ ! -z ${2} ] && [ "${2}" = "terminate" ] && [ ${err} -eq 1 ] ; then echo "stopping"; exit 1; fi
  echo "....ok"
}

# -------------------

#
#  some checks and preliminaries
#

cd ${extraPath}


#
# Main code
#

if [ ${err} -eq 0 ] && [ ${doTrue2eBkg} -eq 1 ] ; then
  evaluateTrue2eBkg
fi

if [ ${err} -eq 0 ] && [ ${doFake2eBkg} -eq 1 ] ; then
  evaluateFake2eBkg
fi


if [ ${err} -eq 1 ] ; then 
    echo 
    echo "   Error was encountered in evaluateDataDrivenBkg.sh"
    echo
    exit 1
    break; 
fi


