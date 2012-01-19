#!/bin/bash

dir1=$PWD
dir2=$1

if [ ${#dir2} -eq 0 ] ; then
  echo "remote file name to compare is not provided"
  exit
fi

fname=${dir2##*/}
#echo fname=${fname}

if [ "$2" == "-r" ] ; then
  diff ${dir1}/${fname} ${dir2}
else
  diff ${dir2} ${dir1}/${fname}
fi

