#!/bin/bash 

#
# some variables
#

# 1) user-defined

#extraPath="../"
extraPath="./"
debugMode=0
dirTag="DY_m10+pr+a05+o03+pr_4680pb";
sourceConfigFile="../config_files/data.conf"       # needed to prepare yields
#sourceConfigFile="../config_files/data_evtTrig.conf"       # needed to prepare yields
mcConfInputFile="../config_files/fall11mc.input"   # needed for unfolding
#mcConfInputFile="../config_files/fall11mc_evtTrig.input"   # needed for unfolding
tmpConfFilePath="../config_files/escale"

anTagUser="ymax9"
anTag="1D${anTagUser}"
triggerSet="Full2011_hltEffOld"

copySelectedNTuples=1   # the main sequence has to be run before this script
   # but the ntuples needs to be copied only once

evaluateSystematics=1   # final calculation
# four individual parts 
doStatisticalStudy=1    # part 1
doShapeSystematics=1    # part 2
doEtaSystematics=1      # part 3
doResidualShapeSystStudy=1   # part 4

# at present (June 20, 2012) doCalculateReferenceShape is untested
doCalculateReferenceShape=0   # if the main sequence was already run
                              #   doCalculateReferenceShape can be set to 0

seedMin=1001
seedMax=1020

defaultEtaDistribution="6binNegs_"
shapeDependenceStudy="Voigtian BreitWigner"
#shapeDependenceStudy="Voigtian"

etaDistrArr="${etaDistrArr} 6bins_Gauss_20120119"
etaDistrArr="${etaDistrArr} 6binNegs_Gauss_20120119"
etaDistrArr="${etaDistrArr} 2binNegs_Gauss 4binNegs_Gauss"
etaDistrArr="${etaDistrArr} 3EB3EENegs_Gauss 4EB3EENegs_Gauss" 
etaDistrArr="${etaDistrArr} 5binNegs_Gauss"


# 2) script-internal


workConfFileT="${tmpConfFilePath}/data_escale_MODEL.conf"
ntuplesDirTag="${dirTag}_escale"
tmpDir="${dirTag}_escale_tmp"
destDirT="${dirTag}_escale_STUDY"
workMCConfInputFileT="../config_files/fall11mc_escale_MODEL.input"
plotsDirExtraTag=""
unfoldingStudy=DYTOOLS::NORMAL
selectionRunMode=DYTOOLS::NORMAL

selEventsDirT="../root_files/selected_events/TMPDIR/"
ntuplesDirT="../root_files/selected_events/TMPDIR/ntuples/"
yieldsDirT="../root_files/yields/TMPDIR/"
#yieldsPlotsDirT="../YieldsAndBackgrounds/TMPDIR"
constDirT="../root_files/constants/TMPDIR/"

ntuplesDirMain=${ntuplesDirT/TMPDIR/${dirTag}_escale}
selEventsDirMain=${selEventsDirT/TMPDIR/${dirTag}_escale}
yieldsDirTmp=${yieldsDirT/TMPDIR/${dirTag}_escale}
constDirTmp=${constDirT/TMPDIR/${dirTag}_escale}
#yieldsPlotsTmp=${yieldsPlotsDirT/TMPDIR/${dirTag}_escale}

runPath=${PWD}
logPath=${PWD}/dir-escale-logs
if [ ! -e ${logPath} ] ; then  mkdir ${logPath}; fi

err=0
calculateUnfolding=1
model=

#
#  some functions
#

# -------------------

checkDirs() {
  if [ ${err} -eq 1 ] ; then return; fi
  if [ ${#ntuplesDir} -eq 0 ] || [ ${#yieldsDir} -eq 0 ] || [ ${#constDir} -eq 0 ] ; then
    echo " one of the directories' variables is empty:"
    echo "   ntuplesDir=${ntuplesDir}"
    echo "   yieldsDir=${yieldsDir}"
    echo "   constDir=${constDir}"
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
  #rm -f ${ntuplesDirTmp}*root 
  rm -f ${yieldsDirTmp}yields${anTag}.root
  rm -f ${yieldsDirTmp}yields_bg-subtracted${anTag}.root
  rm -f ${constDirTmp}unfolding_constants${anTag}.root
  #rm -f ${ntuplesDir}*root 
  rm -f ${yieldsDir}yields${anTag}.root
  rm -f ${yieldsDir}yields_bg-subtracted${anTag}.root
  rm -f ${constDir}unfolding_constants${anTag}.root
}


# Prepare unfolded spectrum

deriveUnfoldedSpectrum() {
  if [ ${err} -eq 1 ] ; then return; fi
  cd ${extraPath}
  checkDirs

#
#  prepareYields
#
  if [ ${calculateUnfolding} -ne 2 ] ; then

  cleanFiles
#
  cd ../Selection
  if [ ${err} -eq 0 ] ; then 
    root -l -b -q selectEvents.C+\(\"${workConfFile}\",\"${triggerSet}\",${selectionRunMode},${debugMode}\) \
        | tee ${logPath}/log-${model}-selectEvents.out
    testFileExists ${ntuplesDir}data${anTagUser}_select${model}.root
  fi
  cd ../YieldsAndBackgrounds
  if [ ${err} -eq 0 ] ; then
    root -l -b -q prepareYields.C+\(\"${workConfFile}\",${selectionRunMode},\"${plotsDirExtraTag}\"\) \
      | tee ${logPath}/log-${model}-prepareYields.out
    testFileExists ${yieldsDir}yields${anTag}.root
  fi
  if [ ${err} -eq 0 ] ; then
    root -l -b -q subtractBackground.C+\(\"${workConfFile}\",${selectionRunMode},\"${plotsDirExtraTag}\"\) \
      | tee ${logPath}/log-${model}-subtractBackground.out
    testFileExists ${yieldsDir}yields_bg-subtracted${anTag}.root
  fi
  fi  # skip derivation
#
#  unfolding
#
  if [ ${err} -eq 0 ] && [ ${calculateUnfolding} -ne 0 ] ; then
    cd ../Unfolding
    root -b -q -l makeUnfoldingMatrix.C+\(\"${workMCConfInputFile}\",\"${triggerSet}\",${unfoldingStudy},1,1.0,-1.0,${debugMode}\) \
	| tee ${logPath}/log-${model}-unfolding.out
    testFileExists ${constDir}/unfolding_constants${anTag}.root
  fi
  
  cd ${runPath}
}


# ----------------------- ntuple copier

cloneNTuples() {
  selEventsDestDir=$1
  if [ ${#selEventsDestDir} -eq 0 ] ; then
    echo "cloneNTuples: argument selEventsDestDir is empty"
    exit
  fi
  srcDir=${selEventsDirT/TMPDIR/${dirTag}}
  echo "copySelectedNTuples: srcDir=\"${srcDir}\"" 
  echo "copySelectedNTuples: destDir=\"${selEventsDestDir}\"" 
  if [ ! -e ${selEventsDestDir} ] ; then mkdir ${selEventsDestDir}; fi
  if [ ! -e "${selEventsDestDir}/ntuples" ] ; 
      then mkdir ${selEventsDestDir}/ntuples; fi

  sets="data qcd ttbar wjets ww wz ztt zz zee"
  for s in ${sets} ; do
    fname="${s}${anTagUser}_select.root"
    set -x  # debug on
    cp ${srcDir}ntuples/${fname} ${selEventsDestDir}ntuples/${fname}
    set +x  # debug off
  done
  set -x  # debug on
  cp ${srcDir}/npv${anTagUser}.root ${selEventsDestDir}/npv${anTagUser}.root
  cp ${selEventsDestDir}npv${anTagUser}.root \
      ${selEventsDestDir}npv${anTagUser}-copy.root
  set +x  # debug off
}

# --------------------  file renamer

renameYields() {
  if [ ${err} -eq 1 ] ; then return; fi
  cd ${extraPath}
  checkDirs
  #if [ ! -e ${yieldsDir} ] ; then mkdir ${yieldsDir} ; fi
  files="yields${anTag}.root  yields_bg-subtracted${anTag}.root"
  files="${files} yield_plots${anTag}.root"
  files="${files} yields_bg-subtracted${anTag}-plots.root"
  for f in ${files} ; do
      testFileExists  ${yieldsDir}/${f}
  done 
  if [ ${err} -eq 0 ] ; then
    for f in ${files} ; do
      destFile=${f/.root/_${model}.root}
      mv  ${yieldsDir}/${f}  ${yieldsDir}/${destFile}
    done
  fi
  cd ${runPath}
}

# -------------------

renameUnfoldedConstants() {
  if [ ${err} -eq 1 ] ; then return; fi
  cd ${extraPath}
  checkDirs
  testFileExists  ${constDir}/unfolding_constants${anTag}.root
  testFileExists  ${constDir}/yields_MC_unfolding_reference_${anTag}.root
  if [ ${err} -eq 0 ] ; then
      if [ ! -d ${constDir} ] ; then mkdir ${constDir}; fi
      echo "renaming unfolded constants"
      mv ${constDir}/unfolding_constants${anTag}.root ${constDir}/unfolding_constants${anTag}_${model}.root
      mv ${constDir}/yields_MC_unfolding_reference_${anTag}.root ${constDir}/yields_MC_unfolding_reference_${anTag}_${model}.root
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

prepareEScaleString() {
    tmp=$escale
    escale=${tmp//\//\\\/}   # replace '/' with '\/' for sed
}

# -------------------

prepareConfFile() {
  escale=$1
  echo "prepareConfFile: escale=<${escale}>"
  if [ ${#escale} -eq 0 ] ; then
    echo "cannot prepareConfFile. Supplied escale is empty"
    exit 2
  fi
  cd ${extraPath}
  echo "pwd=${PWD}"
  if [ ! -z ${workConfFile} ] ; then
    sed "s/#Date20120101_default/${escale}/" ${sourceConfigFile} | \
	sed "s/Date20120101_default/${escale}/" | \
	sed "s/${dirTag}/${tag}/" \
	> ${workConfFile}
  fi
  if [ ! -z ${workMCConfInputFile} ] ; then
     sed "s/#Date20120101_default/${escale}/" ${mcConfInputFile} | \
      sed "s/Date20120101_default/${escale}/" | \
      sed "s/${dirTag}/${tag}/" \
        > ${workMCConfInputFile}
  fi
  cd ${runPath}
}

# -------------------

#
#  some checks and preliminaries
#

cd ${extraPath}

if [ ! -e ${yieldsDirTmp} ] ; then 
    mkdir ${yieldsDirTmp} ; 
fi
if [ ! -e ${tmpConfFilePath} ] ; then
    mkdir ${tmpConfFilePath} ;
fi

# ---
# Drell-Yan DMDY code uses MC evaluation of the backgrounds
# They are available programically
# ---
#bkgrounds="true2eBkgDataPoints.root fakeBkgDataPoints.root"
#for f in  ${bkgrounds} ; do
#  testFileExists ../root_files/yields/${dirTag}/${f}  terminate
#  cp  ../root_files/yields/${dirTag}/${f} \
#      ${yieldsDirTmp}/${f}
#done

#
# Main code
#

#
# Copy n-tuples if needed
#
if [ ${copySelectedNTuples} -gt 0 ] ; then
  srcDir=${ntuplesDirT/TMPDIR/${dirTag}}
  cloneNTuples  ${selEventsDirMain}
fi


if [ ${err} -eq 0 ] && [ ${doStatisticalStudy} -eq 1 ] ; then
  tag="${dirTag}_escale_randomized"
  selectionRunMode=DYTools::ESCALE_STUDY_RND
  selEventsDir=${selEventsDirT/TMPDIR/${tag}}
  ntuplesDir=${ntuplesDirT/TMPDIR/${tag}}
  yieldsDir=${yieldsDirT/TMPDIR/${tag}}
  constDir=${constDirT/TMPDIR/${tag}}
  workConfFile=${workConfFileT/MODEL/randomized}
  workMCConfInputFile=
  seed=${seedMin}
  cloneNTuples  ${selEventsDir}
  while [ ${seed} -le ${seedMax} ] && [ ${err} -eq 0 ] ; do
    model="seed${seed}"
    plotsDirExtraTag="seed${seed}"
    prepareConfFile "Date20120101_default_RANDOMIZED${seed}"
    model="_default20120101_DtRND${seed}"
    calculateUnfolding=0
    deriveUnfoldedSpectrum
    renameYields
    seed=$(( seed + 1 ))
  done
fi

if [ ${err} -eq 0 ] && [ ${doShapeSystematics} -eq 1 ] ; then
  tag="${dirTag}_escale_shape"
  selectionRunMode=DYTools::ESCALE_STUDY
  selEventsDir=${selEventsDirT/TMPDIR/${tag}}
  ntuplesDir=${ntuplesDirT/TMPDIR/${tag}}
  yieldsDir=${yieldsDirT/TMPDIR/${tag}}
  constDir=${constDirT/TMPDIR/${tag}}
  workConfFile=${workConfFileT/MODEL/shape}
  workMCConfInputFile=${workMCConfInputFileT/MODEL/shape}
  cloneNTuples  ${selEventsDir}
  for shape in ${shapeDependenceStudy} ; do
    if [ ${err} -eq 0 ] ; then
      model=${defaultEtaDistribution}${shape}
      plotsDirExtraTag="${shape}"
      prepareConfFile "File${shape}_..\/root_files\/constants\/EScale\/testESF_${model}.inp"
      model="_${model}"
      calculateUnfolding=1
      deriveUnfoldedSpectrum
      renameYields
      renameUnfoldedConstants
    fi
  done
fi


 if [ ${err} -eq 0 ] && [ ${doEtaSystematics} -eq 1 ] ; then
  tag="${dirTag}_escale_eta"
  selectionRunMode=DYTools::ESCALE_STUDY
  selEventsDir=${selEventsDirT/TMPDIR/${tag}}
  ntuplesDir=${ntuplesDirT/TMPDIR/${tag}}
  yieldsDir=${yieldsDirT/TMPDIR/${tag}}
  constDir=${constDirT/TMPDIR/${tag}}
  workConfFile=${workConfFileT/MODEL/eta}
  workMCConfInputFile=${workMCConfInputFileT/MODEL/eta}
  cloneNTuples  ${selEventsDir}
  for aModel in ${etaDistrArr} ; do
    model=${aModel}
    prepareConfFile "FileGauss_..\/root_files\/constants\/EScale\/testESF_${model}.inp"
    plotsDirExtraTag="${model}"
    model="_${model}"
    calculateUnfolding=1
    deriveUnfoldedSpectrum
    renameYields
    renameUnfoldedConstants
  done
fi


if [ ${err} -eq 0 ] && [ ${doResidualShapeSystStudy} -eq 1 ] ; then
  tag="${dirTag}"
  selectionRunMode=DYTools::ESCALE_STUDY
  ntuplesDir=${ntuplesDirT/TMPDIR/${tag}}
  yieldsDir=${yieldsDirT/TMPDIR/${tag}}
  constDir=${constDirT/TMPDIR/${tag}}
  ntuplesDirTmp=${ntuplesDir}
  yieldsDirTmp=${yieldsDir}
  constDirTmp=${constDir}
  workConfFile=${sourceConfigFile}
  workMCConfInputFile=${mcConfInputFile}
  if [ ${doCalculateReferenceShape} -eq 1 ] ; then
    # no need to prepare input files
    calculateUnfolding=1
    selectionRunMode=DYTools::NORMAL
    plotsDirExtraTag="_ref"
    deriveUnfoldedSpectrum
    # no need to rename final files
  fi
  if [ ${err} -eq 0 ] ; then
    calculateUnfolding=2
    unfoldingStudy=DYTOOLS::ESCALE_RESIDUAL
    constDir=${constDirT/TMPDIR/${dirTag}_escale_residual}
    plotsDirExtraTag="_residual"
    saveAnTag=${anTag}
    anTag="${anTag}_escaleResidual"
    deriveUnfoldedSpectrum
    anTag=${saveAnTag}
  fi
  cd ${runPath}
fi


if [ ${err} -eq 0 ] && [ ${evaluateSystematics} -eq 1 ] ; then
  cd ${extraPath}
  cd ../Unfolding
  root -l -q -b calcEscaleSystematics.C+\(\"${dirTag}\",1\)   \
    | tee ${logPath}/log-calcEscaleSystematics${anTag}.log
  cd ${runPath}
fi


 if [ ${err} -eq 1 ] ; then 
     echo 
     echo "   Error was encountered in evaluateEScaleSyst.sh"
     echo
     exit 1
     break; 
  fi


