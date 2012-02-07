#!/bin/bash

# !! This script will be replaced by evaluateESF.sh

#timeStamp="-`date +%Y%m%d-%H%M`"
root -b -q -l auxScriptAdv.C+         | tee out${timeStamp}-09-auxScript-efficiencyScaleFactors.log
