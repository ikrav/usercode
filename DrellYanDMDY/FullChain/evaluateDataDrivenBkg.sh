#!/bin/bash 

#
# some variables
#

debugMode=0    # affects only event (re)selectionStep 
fullRun=0
study2D=1        # study2D does not match DYTools.hh!
reselectEvents=1   # change of study2D does not require reselection, in general
                   # note that only emu events are selected

# 1) user-defined

#extraPath="../"
extraPath="./"
dirTag="DY_j22_19789pb";
workConfFile="../config_files/data_emu.conf"    # needed to select events

anTagUser="" #"ymax9"
#triggerSet="Full2012_hltEffOld"

# if you do not want to have the time stamp, comment the line away 
# or set timeStamp=
timeStamp="`date +%Y%m%d-%H%M`-"
#timeStamp=

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
#resultDirT="../root_files/resultants/TMPDIR/"

eeNtuplesDirMain=${eeSelEventsDirT/TMPDIR/${dirTag}}
emuNtuplesDirMain=${emuSelEventsDirT/TMPDIR/${dirTag}}
#resultDirMain=${resultDirT/TMPDIR/${dirTag}}
resultDirMain="./"

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
  if [ ${#eeNtuplesDirMain} -eq 0 ] || [ ${#emuNtuplesDirMain} -eq 0 ] || [ ${#yieldsDir} -eq 0 ] || [ ${#resultDirMain} -eq 0 ] ; then
    echo " one of the directories' variables is empty:"
    echo "   ntuplesDirs: <${eeNtuplesDirMain}>, <${emuNtuplesDirMain}>"
    echo "   resultDir: <${resultDirMain}>"
    exit 1
  fi
}


# -------------------- main cycle

cleanFiles() {
#
#  clean up the files to make sure we are producing them now
#
    runCase=$1
    if [ ${runCase} == "true2e" ] ; then
	echo "remove true2e file"
	rm -f ${resultDirMain}/true2eBkgDataPoints_${anTag}.root
    fi
}


# -----------------  selectEvents ----------------
evaluateTrue2eBkg() {

  cleanFiles "true2e"
#
#  first try to compile the needed code
#
  cd ../DataDrivenBackgrounds
  if [ ${err} -eq 0 ] ; then
      cd eMuMethod
      rm -f eMuBkgExe
      gmake eMuBkgExe 2>&1 | tee ${logPath}/log${timeStamp}-gMake-eMu.out
      testFileExists eMuBkgExe
      cd ..
  fi
#
#   select events
# 
  if [ ${err} -eq 0 ] && [ ${reselectEvents} -eq 1 ] ; then 
    rm selectEmuEvents_C.so
    root -l -b -q selectEmuEvents.C+\(\"${workConfFile}\",${debugMode}\) \
        | tee ${logPath}/log${timeStamp}-selectEmuEvents.out
    testFileExists selectEmuEvents_C.so
  fi
  testFileExists ${emuNtuplesDirMain}/data${anTagUser}_select.root \
      ${emuNtuplesDirMain}/zee${anTagUser}_select.root
  testFileExists ${eeNtuplesDirMain}/data${anTagUser}_select.root \
      ${eeNtuplesDirMain}/zee${anTagUser}_select.root

#
#   evaluate true2e backgrounds
#
  if [ ${err} -eq 0 ] ; then
      cd eMuMethod
      flags=" --saveRootFile"
      flags="${flags} --doPUreWeight"
      flags="${flags} --verbose" # debug
      flags="${flags} --dirTag ${dirTag}"
      if [ ${study2D} -eq 1 ] ; then
	  flags="--doDMDY ${flags}"
      fi
      ./eMuBkgExe ${flags} 2>&1 | tee ${logPath}/log${timeStamp}-emu.out
      testFileExists ${resultDirMain}/true2eBkgDataPoints_${anTag}.root
      cd ..
  fi
  
  cd ${runPath}
}


# -------------------

testFileExists() {
  if [ $? != 0 ] ; then err=1; fi
  file=$1
  echo "evaluateDataDrivenBkg.sh: test file <${file}>"
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


