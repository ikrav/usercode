#!/bin/bash

mcInput=$1
debugMode=0

if [ "$2" == "-debug" ] ; then
  debugMode=$3
  echo "setupTheory: setting debugMode=${debugMode}"
fi

if [ ${#mcInput} -eq 0 ] ; then
    mcInput="../config_files/fall11mc.input"
fi

# acceptance correction theoretical uncertainty used in 2011 summer
root -l -q -b TheoryErrors.C+\(\"${mcInput}\"\)

# theoretical cross section used in 2011
root -l -q -b createThXSec1Dsummer2011.C+


# theoretical cross section from Zee MC signal sample
useFEWZarr="true false"
#useFEWZarr="true"
#fineGridArr="0 1 2"
#fineGridArr="1"
fineGridArr="0 1 2 3"
for fineGrid in ${fineGridArr} ; do
  for useFEWZ in ${useFEWZarr} ; do
      echo "useFEWZ=${useFEWZ}, fineGrid=${fineGrid}"
      root -l -q -b getXsecExtended.C+\(\"${mcInput}\",${debugMode},${useFEWZ},${fineGrid}\) | tee log-xSec-FEWZ${useFEWZ}-fineGrid${fineGrid}-debug${debugMode}.out
  done
done

