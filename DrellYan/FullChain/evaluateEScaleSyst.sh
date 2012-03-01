#!/bin/bash

#
# some variables
#

# 1) user-defined

extraPath="../"
extraPath="./"
dirTag="DY_m10+pr+a05+o03+pr_4680pb";
sourceConfigFile="../config_files/data.conf"       # needed to prepare yields
mcConfInputFile="../config_files/fall11mc.input"   # needed for unfolding

evaluateSystematics=1   # final calculation
# four individual parts 
doStatisticalStudy=1    # part 1
doShapeSystematics=1    # part 2
doEtaSystematics=1      # part 3
doResidualShapeSystStudy=1   # part 4

doCalculateReferenceShape=1   # if the main sequence was already run
                              # doCalculateReferenceShape can be set to 0

seedMin=1000
seedMax=1020

defaultEtaDistribution="6binNegs_"
shapeDependenceStudy="Voigtian BreitWigner"

etaDistrArr="${etaDistrArr} 6bins_Gauss_20120119"
etaDistrArr="${etaDistrArr} 6binNegs_Gauss_20120119"
#etaDistrArr="${etaDistrArr} 2binNegs_Gauss 4binNegs_Gauss"
#etaDistrArr="${etaDistrArr} 3EB3EENegs_Gauss 4EB3EENegs_Gauss" 
#etaDistrArr="${etaDistrArr} 5binNegs_Gauss"


# 2) script-internal


workConfFileT="../config_files/data_escale_MODEL.conf"
tmpDir="${dirTag}_escale_tmp"
destDirT="${dirTag}_escale_STUDY"
workMCConfInputFileT="../config_files/fall11mc_escale_MODEL.input"
unfoldingStudy=DYTOOLS::NORMAL

ntuplesDirTmp="../root_files/selected_events/${tmpDir}/ntuples/"
yieldsDirTmp="../root_files/yields/${tmpDir}/"
constDirTmp="../root_files/constants/${tmpDir}/"

ntuplesDirT="../root_files/selected_events/TMPDIR/ntuples/"
yieldsDirT="../root_files/yields/TMPDIR/"
constDirT="../root_files/constants/TMPDIR/"

runPath=${PWD}
logPath=${PWD}/dir-Logs
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


deriveUnfoldedSpectrum() {
  if [ ${err} -eq 1 ] ; then return; fi
  cd ${extraPath}
  checkDirs

  if [ ${calculateUnfolding} -ne 2 ] ; then
#
#  clean up the files to make sure we are producing them now
#
  rm -f ${ntuplesDirTmp}*root 
  rm -f ${yieldsDirTmp}yields.root
  rm -f ${yieldsDirTmp}massHist.root
  rm -f ${yieldsDirTmp}yields_bg-subtracted.root
  rm -f ${constDirTmp}unfolding_constants.root
  rm -f ${ntuplesDir}*root 
  rm -f ${yieldsDir}yields.root
  rm -f ${yieldsDir}massHist.root
  rm -f ${yieldsDir}yields_bg-subtracted.root
  rm -f ${constDir}unfolding_constants.root

#
  cd ../Selection
  if [ ${err} -eq 0 ] ; then
    root -l -b -q selectEvents1D.C+\(\"${workConfFile}\"\) \
        | tee ${logPath}/log-${model}-selectEvents.out
    testFileExists ${ntuplesDirTmp}zee_select.root
  fi
  cd ../MassSpectrum
  if [ ${err} -eq 0 ] ; then
    root -l -b -q prepareYields1D.C+\(\"${workConfFile}\"\) \
      | tee ${logPath}/log-${model}-prepareYields.out
    testFileExists ${yieldsDirTmp}yields.root
    testFileExists ${yieldsDirTmp}massHist.root
  fi
  if [ ${err} -eq 0 ] ; then
    root -l -b -q subtractBackground.C+\(\"${workConfFile}\"\) \
      | tee ${logPath}/log-${model}-subtractBackground.out
    testFileExists ${yieldsDirTmp}yields_bg-subtracted.root
  fi
  fi  # skip derivation
#
  if [ ${err} -eq 0 ] && [ ${calculateUnfolding} -ne 0 ] ; then
    cd ../Unfolding
    root -b -q -l plotDYUnfoldingMatrix.C+\(\"${workMCConfInputFile}\",${unfoldingStudy}\)
    testFileExists ${constDirTmp}/unfolding_constants.root
  fi
  
  cd ${runPath}
}

# --------------------  file renamer

renameYields() {
  if [ ${err} -eq 1 ] ; then return; fi
  cd ${extraPath}
  checkDirs
  if [ ! -e ${yieldsDir} ] ; then mkdir ${yieldsDir} ; fi
  testFileExists ${yieldsDirTmp}/massHist.root
  testFileExists ${yieldsDirTmp}/yields.root
  testFileExists ${yieldsDirTmp}/yields_bg-subtracted.root
  mv ${yieldsDirTmp}/massHist.root ${yieldsDir}/massHist_${model}.root
  mv ${yieldsDirTmp}/yields.root ${yieldsDir}/yields_${model}.root
  mv ${yieldsDirTmp}/yields_bg-subtracted.root ${yieldsDir}/yields_bg-subtracted_${model}.root
  cd ${runPath}
}

# -------------------

renameUnfoldedConstants() {
  if [ ${err} -eq 1 ] ; then return; fi
  cd ${extraPath}
  checkDirs
  testFileExists  ${constDirTmp}/unfolding_constants.root
  testFileExists  ${constDirTmp}/yields_MC_unfolding_reference.root
  if [ ${err} -eq 0 ] ; then
      if [ ! -d ${constDir} ] ; then mkdir ${constDir}; fi
      echo "renaming unfolded constants"
      mv ${constDirTmp}/unfolding_constants.root ${constDir}/unfolding_constants_${model}.root
      mv ${constDirTmp}/yields_MC_unfolding_reference.root ${constDir}/yields_MC_unfolding_reference_${model}.root
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
	sed "s/${dirTag}/${tmpDir}/" \
	> ${workConfFile}
  fi
  if [ ! -z ${workMCConfInputFile} ] ; then
     sed "s/#Date20120101_default/${escale}/" ${mcConfInputFile} | \
      sed "s/Date20120101_default/${escale}/" | \
      sed "s/${dirTag}/${tmpDir}/" \
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

#tmpFullPath="../root_files/
#if [ ! -e ${

bkgrounds="true2eBkgDataPoints.root fakeBkgDataPoints.root"
for f in  ${bkgrounds} ; do
  testFileExists ../root_files/yields/${dirTag}/${f}  terminate
  cp  ../root_files/yields/${dirTag}/${f} \
      ${yieldsDirTmp}/${f}
done

#
# Main code
#

if [ ${err} -eq 0 ] && [ ${doStatisticalStudy} -eq 1 ] ; then
  tag="${dirTag}_escale_randomized"
  ntuplesDir=${ntuplesDirT/TMPDIR/${tag}}
  yieldsDir=${yieldsDirT/TMPDIR/${tag}}
  constDir=${constDirT/TMPDIR/${tag}}
  workConfFile=${workConfFileT/MODEL/randomized}
  workMCConfInputFile=
  seed=${seedMin}
  while [ ${seed} -le ${seedMax} ] && [ ${err} -eq 0 ] ; do
    model="seed${seed}"
    prepareConfFile "Date20120101_default RANDOMIZED${seed}"
    calculateUnfolding=0
    deriveUnfoldedSpectrum
    renameYields
    seed=$(( seed + 1 ))
  done
fi

if [ ${err} -eq 0 ] && [ ${doShapeSystematics} -eq 1 ] ; then
  tag="${dirTag}_escale_shape"
  ntuplesDir=${ntuplesDirT/TMPDIR/${tag}}
  yieldsDir=${yieldsDirT/TMPDIR/${tag}}
  constDir=${constDirT/TMPDIR/${tag}}
  workConfFile=${workConfFileT/MODEL/shape}
  workMCConfInputFile=${workMCConfInputFileT/MODEL/shape}
  for shape in ${shapeDependenceStudy} ; do
    if [ ${err} -eq 0 ] ; then
      model=${defaultEtaDistribution}${shape}
      prepareConfFile "File${shape} ..\/root_files\/constants\/testESF_${model}.inp"
      calculateUnfolding=1
      deriveUnfoldedSpectrum
      renameYields
      renameUnfoldedConstants
    fi
  done
fi


 if [ ${err} -eq 0 ] && [ ${doEtaSystematics} -eq 1 ] ; then
  tag="${dirTag}_escale_eta"
  ntuplesDir=${ntuplesDirT/TMPDIR/${tag}}
  yieldsDir=${yieldsDirT/TMPDIR/${tag}}
  constDir=${constDirT/TMPDIR/${tag}}
  workConfFile=${workConfFileT/MODEL/eta}
  workMCConfInputFile=${workMCConfInputFileT/MODEL/eta}
  for model in ${etaDistrArr} ; do
    prepareConfFile "FileGauss ..\/root_files\/constants\/testESF_${model}.inp"
    calculateUnfolding=1
    deriveUnfoldedSpectrum
    renameYields
    renameUnfoldedConstants
  done
fi


if [ ${err} -eq 0 ] && [ ${doResidualShapeSystStudy} -eq 1 ] ; then
  tag="${dirTag}"
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
    deriveUnfoldedSpectrum
    # no need to rename final files
  fi
  if [ ${err} -eq 0 ] ; then
    cd ${extraPath}../Unfolding
    root -l -q -b create_mass_shape_weights.C+\(\"${dirTag}\"\)
    testFileExists ${yieldsDir}shape_weights.root
  fi
  if [ ${err} -eq 0 ] ; then
    calculateUnfolding=2
    unfoldingStudy=DYTOOLS::ESCALE_RESIDUAL
    constDirTmp=${constDirT/TMPDIR/${dirTag}_escale_residual}
    deriveUnfoldedSpectrum
  fi
  cd ${runPath}
fi


if [ ${err} -eq 0 ] && [ ${evaluateSystematics} -eq 1 ] ; then
  cd ${extraPath}
  cd ../Unfolding
  root -l -q -b calcEscaleSystematics.C+\(\"${dirTag}\",1\)
  cd ${runPath}
fi


 if [ ${err} -eq 1 ] ; then 
     echo 
     echo "   Error was encountered in evaluateEScaleSyst.sh"
     echo
     exit 1
     break; 
  fi


