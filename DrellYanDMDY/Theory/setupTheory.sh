#!/bin/bash

mcInput=$1
debugMode=0

if [ ${#mcInput} -eq 0 ] ; then
    mcInput="../config_files/fall11mc.input"
fi

# acceptance correction theoretical uncertainty used in 2011 summer
root -l -q -b TheoryErrors.C+

# theoretical cross section used in 2011
root -l -q -b createThXSec1Dsummer2011.C+


# theoretical cross section from Zee MC signal sample
useFEWZarr="true false"
fineGridArr="0 1"
for useFEWZ in ${useFEWZarr} ; do
  for fineGrid in ${fineGridArr} ; do
      echo "useFEWZ=${useFEWZ}, fineGrid=${fineGrid}"
      root -l -q -b getXsecExtended.C+\(\"${mcInput}\",${debugMode},${useFEWZ},${fineGrid}\)
  done
done

