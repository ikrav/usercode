#!/bin/bash

#flist="./file_list_test.txt"
flist="./files_to_copy.txt"
echo "Process file list from $flist"
echo

sourceDirLocal="/data/blue/ksung/DYAna53X/"
#sourceDirLocal="/home/ikrav/"
targetDirEOS="/eos/cms/store/group/phys_smp/DYToEE_diff_xsec/ntuples/53X/"

# For some reason, in the shell of this script the EOS command is not defined,
# so explicitly set it to what the interactive shell shows under "which eos"
eoscommand="/afs/cern.ch/project/eos/installation/0.2.30/bin/eos.select"

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
    $eoscommand stat $fileEOS >& /dev/null
    result=$?
    if [ $result == 0 ]; then
        echo "File ALREADY EXISTS: $fileEOS"
	targetOK=0
    else
        echo "File not yet copied, target location ok: $fileEOS"
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
    fileLocal="${sourceDirLocal}/${file}"
    fileEOS="${targetDirEOS}/${file}"
    command="${eoscommand} cp ${fileLocal} $fileEOS"
    echo "Executing $command"
    $command
    # check that source and target have the same size
    # but before that, check that the source and the target both exist.
    # if either doesn't, the copy has obviously failed.
    # handle the source
    sourceSize=-1
    if [ -e $fileLocal ]; then
	sourceSize=`stat --format="%s" ${fileLocal}`
    fi
    # handle the target
    targetSize=-2 # note, different from non-existing sourceSize=-1
    $eoscommand stat $fileEOS >& /dev/null
    eosFileCheck=$?
    if [ $eosFileCheck == 0 ]; then
	targetSize=`${eoscommand} ls -l ${fileEOS} | awk '{print $5}'`
    fi
    # compare the sizes
    if [ $sourceSize == $targetSize ]; then
        echo "     copy is successful"
    else
	echo "     copy has FAILED"
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
