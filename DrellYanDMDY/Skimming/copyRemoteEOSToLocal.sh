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
flist="./files_to_copy.txt"
echo "Process file list from $flist"
echo

sourceHostEOS="eoscms.cern.ch"
sourceDirEOS="/eos/cms/store/group/phys_smp/DYToEE_diff_xsec/ntuples/53X/"
targetDirLocal="/mnt/hadoop/user/uscms01/pnfs/unl.edu/data4/cms/store/user/ikrav/DrellYan_8TeV/53X/"

#
# Check that source files all exist
#
sourceOK=1
list=`cat $flist`
while read file
do
    fileEOS="${sourceDirEOS}/${file}"
    xrd $sourceHostEOS existfile $fileEOS >& /dev/null
    errorCode=$?
    if [ $errorCode == 0 ]; then
        echo "File source checked ok: root://${sourceHostEOS}/${fileEOS}"
    else
        echo "File NOT FOUND:  root://${sourceHostEOS}/${fileEOS}"
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
    fileLocal="${targetDirLocal}/${file}"
    if [ -e $fileLocal ]; then
        echo "File ALREADY EXISTS: $fileLocal"
	targetOK=0
    else
        echo "File not yet copied, target location ok: $fileLocal"
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
    fileEOS="${sourceDirEOS}/${file}"
    urlEOS="root://${sourceHostEOS}/${fileEOS}"
    fileLocal="${targetDirLocal}/${file}"
    command="xrdcp $urlEOS ${fileLocal}"
    echo "Executing $command"
    $command
    # check that source and target have the same size
    # but before that, check that the source and the target both exist.
    # if either doesn't, the copy has obviously failed.
    # handle the source
    xrd $sourceHostEOS existfile $fileEOS >& /dev/null
    errorCode=$?
    eosFileCheck=$?
    if [ $eosFileCheck == 0 ]; then
	xrdsize="xrd $sourceHostEOS stat ${fileEOS}"
	sourceSize=`$xrdsize | awk '{if($3=="Size:"){printf("%s",$4);}}'`
    fi
    # handle the target
    targetSize=-2 # note, different from non-existing sourceSize=-1
    if [ -e $fileLocal ]; then
	targetSize=`stat --format="%s" ${fileLocal}`
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
