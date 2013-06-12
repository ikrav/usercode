#!/bin/bash

files=`find $@ -maxdepth 1 -type f`
if [ ${#files} -eq 0 ] ; then
 echo "no files provided"
 exit
fi

for f in ${files} ; do 
  md5sum ${f}
done
