#!/bin/bash

# Check if grid proxy is there, it will be needed to run xrdcp
grid-proxy-info -e >& /dev/null
proxyError=$?
if [ $proxyError == 0 ]; then
    echo
    echo "Valid proxy is found, proceeding"
    echo
else
    echo 
    echo "No proxy found, can not run xrdcp"
    echo
    exit
fi

#flist="./file_list_test.txt"
flist="./files_to_copy_toEOS.txt"
echo "Process file list from $flist"
echo

targetHostEOS="eoscms.cern.ch"
targetDirEOS="/eos/cms/store/group/phys_smp/DYToEE_diff_xsec/ntuples/53X/"
sourceDirLocal="/home/hep/ikrav/work/ntuples/DrellYan_8TeV_53X_local/"

#
# Check that source files all exist
#

sourceOK=1
list=`cat $flist`
while read file
do
    fileLocal="${sourceDirLocal}/${file}"
    if [ -e $fileLocal ]; then
        echo "File source checked ok: $fileLocal"
    else
        echo "File NOT FOUND: $fileLocal"
	sourceOK=0
    fi
done <<< "$list"

if [ $sourceOK == 1 ]; then
   echo 
   echo "All source files are present, proceeding."
   echo
else
   echo
   echo "One or more source files are missing, exiting, no copy is done"
   echo
   exit
fi

#
# Check that target files DO NOT exist yet
#

targetOK=1
list=`cat $flist`
while read file
do
    fileEOS="${targetDirEOS}/${file}"
    xrd $targetHostEOS existfile $fileEOS >& /dev/null
    errorCode=$?
    if [ $errorCode == 0 ]; then
        echo "File target ALREADY EXISTS: root://${targetHostEOS}/${fileEOS}"
	targetOK=0
    else
        echo "File not yet copied, target location ok:  root://${targetHostEOS}/${fileEOS}"
    fi
done <<< "$list"

if [ $targetOK == 1 ]; then
   echo 
   echo "No target files exist yet, we are free to copy."
   echo
else
   echo
   echo "One or more target files is already there, exiting, no copy is done"
   echo
   exit
fi

#
# Execute the copying
#

list=`cat $flist`
copyOK=1
while read file
do
    fileEOS="${targetDirEOS}/${file}"
    urlEOS="root://${targetHostEOS}/${fileEOS}"
    fileLocal="${sourceDirLocal}/${file}"
    command="xrdcp ${fileLocal} $urlEOS"
    echo "Executing $command"
    $command
    # check that source and target have the same size
    # but before that, check that the source and the target both exist.
    # if either doesn't, the copy has obviously failed.
    # handle the source
    sourceSize=-1 # note, different from non-existing sourceSize=-1
    if [ -e $fileLocal ]; then
	sourceSize=`stat --format="%s" ${fileLocal}`
    fi
    # handle the target
    xrd $targetHostEOS existfile $fileEOS >& /dev/null
    errorCode=$?
    eosFileCheck=$?
    if [ $eosFileCheck == 0 ]; then
	xrdsize="xrd $targetHostEOS stat ${fileEOS}"
	targetSize=`$xrdsize | awk '{if($3=="Size:"){printf("%s",$4);}}'`
    fi
    # compare the sizes
    if [ $sourceSize == $targetSize ]; then
        echo "     copy is successful"
    else
	echo "     copy has FAILED"
	echo "        source size: $sourceSize"
	echo "        target size: $targetSize"
        copyOK=0
    fi
     
done <<< "$list"

if [ $copyOK == 1 ]; then
    echo
    echo "All files have been successfully copied"
    echo
else
    echo
    echo "One or more files failed to copy"
    echo
fi
